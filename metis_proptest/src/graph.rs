use std::fs::OpenOptions;
use std::path::Path;
use std::sync::Arc;
use std::io::{self, prelude::*, BufReader};

use crate::utils::zip_first;

#[derive(Clone)]
pub struct Graph {
    // we use Arc<Vec<usize>> (with an added indirection) because many of our operations will
    // change the lengths of the arrays. If we used Arc<[usize]>, then converting to a Vec would
    // always entail a full clone, even if there was only one reference. I haven't tested it, but I
    // expect this to almost always be worth the added cost of indirection.

    pub xadj: Arc<Vec<usize>>,
    pub adjncy: Arc<Vec<usize>>,
    pub ncon: usize,
    pub adjwgt: Option<Arc<Vec<usize>>>,
    pub vwgt: Option<Arc<Vec<usize>>>,
    pub vsize: Option<Arc<Vec<usize>>>,
}

pub struct VertexFields<'a> {
    pub vtx: usize,
    pub edges: &'a [usize],
    pub adjwgt: Option<&'a [usize]>,
    pub vwgt: Option<&'a [usize]>,
    pub vsize: Option<usize>,
}

impl VertexFields<'_> {
    /// iterator of items in order they appear in the graph file
    fn file_order_iter(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.vsize.as_ref().into_iter().copied()
            .chain(self.vwgt.into_iter().flatten().copied())
            .chain(zip_first(self.edges.iter().map(|e| e + 1), self.adjwgt.into_iter().flatten().copied()).flat_map(|(v, e)| std::iter::once(v).chain(e)))
    }

    /// number of fields on the line in the graph file
    fn file_field_num(&self) -> usize {
        self.vsize.is_some() as usize
        + self.vwgt.map_or(0, |vwgt| vwgt.len())
        + self.edges.len() * (1 + self.adjwgt.is_some() as usize)
    }
}

impl Graph {
    pub fn nvtxs(&self) -> usize {
        self.xadj.len().checked_sub(1).unwrap()
    }

    pub fn nedges(&self) -> usize {
        self.adjncy.len() / 2
    }

    pub fn edges_iter(&self) -> impl Iterator<Item = &[usize]> {
        self.xadj.windows(2).map(|w| &self.adjncy[w[0]..w[1]])
    }

    pub fn vtx_field_iter(&self) -> impl Iterator<Item = VertexFields> {
        self.xadj.windows(2).enumerate().map(move |(vtx, w)| {
            VertexFields { 
                vtx,
                edges: &self.adjncy[w[0]..w[1]],
                adjwgt: self.adjwgt.as_deref().map(|a| &a[w[0]..w[1]]),
                vwgt: self.vwgt.as_deref().map(|v| &v[(vtx * self.ncon)..((vtx + 1) * self.ncon)]),
                vsize: self.vsize.as_deref().map(|v| v[vtx])
            }
        })
    }

    pub fn write(&self, mut w: impl Write) -> io::Result<()> {

        // write header
        let has_fmt = self.vwgt.is_some() || self.vsize.is_some() || self.adjwgt.is_some();
        let fmt_num = (self.adjwgt.is_some() as i32)
        + (10 * self.vwgt.is_some() as i32)
        + (100 * self.vsize.is_some() as i32);
        let nvtxs = self.nvtxs();
        let nedges = self.nedges();
        let ncon = self.ncon;
        if has_fmt && ncon > 0 {
            writeln!(w, "{nvtxs} {nedges} {fmt_num:03} {ncon}")?;
        } else if has_fmt {
            writeln!(w, "{nvtxs} {nedges} {fmt_num:03}")?;
        } else {
            assert!(fmt_num == 0);
            writeln!(w, "{nvtxs} {nedges}")?;
        }

        for fields in self.vtx_field_iter() {
            let n = fields.file_field_num();
            if n == 0 {
                writeln!(w)?;
                continue
            }
            let mut it = fields.file_order_iter().peekable();
            for field in (&mut it).take(n - 1) {
                write!(w, "{field} ")?;
            }
            writeln!(w, "{}", it.next().unwrap())?;
        }
        w.flush()
    }

    pub fn read_from_path(p: impl AsRef<Path>) -> io::Result<Self> {
        let f = OpenOptions::new().read(true).open(p)?;
        let r = BufReader::with_capacity(32 * 1024, f);
        Self::read(r)
    }

    pub fn read(mut r: impl BufRead) -> io::Result<Self> {
        let mut buf = String::with_capacity(256);

        // read header
        if r.read_line(&mut buf)? == 0 {
            panic!("premature EOF")
        }
        while buf.starts_with('%') {
            buf.clear();
            if r.read_line(&mut buf)? == 0 {
                panic!("premature EOF")
            };
        }
        let mut hfields = buf.split_ascii_whitespace().map(|f| f.parse::<usize>().expect("non-numeric header field"));
        let nvtxs = hfields.next().expect("no nvtxs header field");
        let nedges = hfields.next().expect("no nedges header field");
        let fmt = hfields.next().unwrap_or(0);
        let ncon = hfields.next().unwrap_or(0);

        assert!(fmt <= 111, "invalid fmt");
        assert!((fmt / 1  ) % 10 | 1 == 1, "invalid fmt");
        assert!((fmt / 10 ) % 10 | 1 == 1, "invalid fmt");
        assert!((fmt / 100) % 10 | 1 == 1, "invalid fmt");

        let has_adjwgt = (fmt / 1  ) % 10 == 1;
        let has_vwgt   = (fmt / 10 ) % 10 == 1;
        let has_vsize  = (fmt / 100) % 10 == 1;

        let mut xadj:   Vec<usize> = Vec::with_capacity(nvtxs + 1);
        let mut adjncy: Vec<usize> = Vec::with_capacity(nedges * 2);
        let mut adjwgt: Vec<usize> = Vec::with_capacity(if has_adjwgt { nedges * 2 } else { 0 });
        let mut vwgt:   Vec<usize> = Vec::with_capacity(if has_vwgt { nvtxs * ncon } else { 0 });
        let mut vsize:  Vec<usize> = Vec::with_capacity(if has_vsize { nvtxs } else { 0 });

        xadj.push(0);

        for _ in 0..nvtxs {
            buf.clear();
            if r.read_line(&mut buf)? == 0 {
                panic!("premature EOF")
            };
            while buf.starts_with('%') {
                buf.clear();
                if r.read_line(&mut buf)? == 0 {
                    panic!("premature EOF")
                };
            }

            let mut fields = buf.split_ascii_whitespace().map(|f| f.parse::<usize>().expect("non-numeric entry"));

            if has_vsize {
                vsize.push(fields.next().expect("missing vsize"));
            }

            if has_vwgt {
                for _ in 0..ncon {
                    vwgt.push(fields.next().expect("missing vwgt"))
                }
            }

            if has_adjwgt {
                while let Some(v) = fields.next() {
                    adjncy.push(v.checked_sub(1).expect("zero-indexed adjncy"));
                    adjwgt.push(fields.next().expect("missing corresponding adjwgt"));
                }
            } else {
                adjncy.extend(fields.map(|v| v.checked_sub(1).expect("zero-indexed adjncy")));
            }

            xadj.push(adjncy.len());
        }

        'nojunk: {
            buf.clear();
            if r.read_line(&mut buf)? == 0 {
                break 'nojunk
            };
            while buf.starts_with('%') {
                buf.clear();
                if r.read_line(&mut buf)? == 0 {
                    break 'nojunk
                };
            }

            panic!("junk line after graph: {buf:?}")
        }

        assert_eq!(xadj.len(), nvtxs + 1);
        assert_eq!(adjncy.len(), nedges * 2);
        if has_adjwgt {
            assert_eq!(adjwgt.len(), nedges * 2);
        }
        if has_vwgt {
            assert_eq!(vwgt.len(), nvtxs * ncon);
        }
        if has_vsize {
            assert_eq!(vsize.len(), nvtxs);
        }


        let ret = Graph {
            xadj: xadj.into(),
            adjncy: adjncy.into(),
            ncon,
            adjwgt: has_adjwgt.then(|| adjwgt.into()),
            vwgt: has_vwgt.then(|| vwgt.into()),
            vsize: has_vsize.then(|| vsize.into()),
        };

        Ok(ret)
    }
}
