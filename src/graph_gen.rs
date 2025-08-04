use std::{error::Error, io::{self, BufRead}, ops::Range};

use fastrand::Rng;

use crate::{
    dal::{DirectAccessList, DirectAccessMap},
    util::make_csr,
    *,
};

/// base representation of a graph - we can use it to avoid reparsing graphs
#[derive(Clone)]
pub struct Csr {
    xadj: Box<[idx_t]>,
    adjncy: Box<[idx_t]>,
}

impl Csr {

    /// returns (xadj, adjncy)
    #[cfg(test)]
    pub fn into_parts(self) -> (Box<[idx_t]>, Box<[idx_t]>) {
        let Csr { xadj, adjncy } = self;
        (xadj, adjncy)
    }

    #[cfg(test)]
    pub fn as_parts(&self) -> (&[idx_t], &[idx_t]) {
        let Csr { xadj, adjncy } = self;
        (&xadj, &adjncy)
    }

    #[cfg(test)]
    pub fn nvtxs(&self) -> usize {
        self.xadj.len() - 1
    }

    #[cfg(test)]
    pub(crate) fn test_graph(tg: tests::TestGraph) -> Self {
        GraphBuilder::test_graph(tg, Optype::Ometis, 3, 1).to_csr()
    }
}

#[derive(Clone)]
pub struct GraphBuilder {
    dbg_lvl: u32,
    op: GraphOpSettings,
    seed: idx_t,
    xadj: Vec<idx_t>,
    adjncy: Vec<idx_t>,
    vwgt: Option<Vec<idx_t>>,
    edge_match: Ctype,
    initial_part: Iptype,
}

#[derive(Clone)]
enum KwayObjective {
    Vol {
        vsize: Option<Vec<idx_t>>,
    },
    Cut {
        adjwgt: Option<Vec<idx_t>>,
    },
}

impl KwayObjective {
    fn objtype(&self) -> Objtype {
        match self {
            KwayObjective::Vol { .. } => Objtype::Vol,
            KwayObjective::Cut { .. } => Objtype::Cut,
        }
    }
}

/// settings common between partition routines (pmetis and kmetis)
#[derive(Clone)]
struct GraphPartSettings {
    ncuts: idx_t,
    nparts: usize,
    ncon: usize,
    ubvec: Option<Vec<real_t>>,
    tpwgts: Option<Vec<real_t>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Refinement algorithm for ometis
pub enum OmetisRtype {
    Sep1Sided,
    Sep2Sided,
}

/// refer to options.c for what's actually used
#[derive(Clone)]
enum GraphOpSettings {
    Pmetis {
        common: GraphPartSettings,
        adjwgt: Option<Vec<idx_t>>,
    },
    Kmetis {
        common: GraphPartSettings,
        objtype: KwayObjective,
        minconn: bool,
        contig: bool,
        niparts: idx_t,
    },
    Ometis {
        nseps: idx_t,
        pfactor: idx_t,
        ccorder: bool,
        compress: bool,
        compute_vertex_separator: bool,
        rtype: OmetisRtype,
    },
}

impl GraphOpSettings {
    fn op(&self) -> Optype {
        match self {
            GraphOpSettings::Pmetis { .. } => Optype::Pmetis,
            GraphOpSettings::Kmetis { .. } => Optype::Kmetis,
            GraphOpSettings::Ometis { .. } => Optype::Ometis,
        }
    }

    fn common(&mut self) -> Option<&mut GraphPartSettings> {
        match self {
            GraphOpSettings::Pmetis { common, .. } => Some(common),
            GraphOpSettings::Kmetis { common, .. } => Some(common),
            GraphOpSettings::Ometis { .. } => None,
        }
    }

    fn common_ref(&self) -> Option<&GraphPartSettings> {
        match self {
            GraphOpSettings::Pmetis { common, .. } => Some(common),
            GraphOpSettings::Kmetis { common, .. } => Some(common),
            GraphOpSettings::Ometis { .. } => None,
        }
    }
}

fn vec_ptr<T>(v: &mut Option<Vec<T>>) -> *mut T {
    match v.as_mut() {
        Some(v) => v.as_mut_ptr(),
        None => std::ptr::null_mut(),
    }
}

#[allow(unused)]
impl GraphBuilder {
    pub fn new(op: Optype, nparts: usize, ncon: usize) -> Self {
        if op.is_ometis() {
            assert_eq!(ncon, 1, "ometis requires ncon=1");
        }
        if op.is_ometis() {
            assert_eq!(nparts, 3, "ometis requires nparts=3");
        }
        Self {
            op: match op {
                Optype::Pmetis => GraphOpSettings::Pmetis { 
                    common: GraphPartSettings { 
                        ncuts: 1,
                        nparts,
                        ncon,
                        ubvec: None,
                        tpwgts: None,
                    },
                    adjwgt: None,
                },
                Optype::Kmetis => GraphOpSettings::Kmetis {
                    common: GraphPartSettings { 
                        ncuts: 1,
                        nparts,
                        ncon,
                        ubvec: None,
                        tpwgts: None,
                    },
                    objtype: KwayObjective::Cut { 
                        adjwgt: None
                    },
                    minconn: false,
                    contig: false,
                    niparts: -1,
                },
                Optype::Ometis => GraphOpSettings::Ometis { 
                    compute_vertex_separator: false,
                    compress: true,
                    nseps: 1,
                    pfactor: 0,
                    ccorder: false,
                    rtype: OmetisRtype::Sep1Sided
                },
            },
            seed: 4321,
            xadj: Vec::new(),
            adjncy: Vec::new(),
            vwgt: None,
            edge_match: Ctype::Shem,
            initial_part: match op {
                Optype::Pmetis => Iptype::Grow,
                Optype::Kmetis => Iptype::Rb,
                Optype::Ometis => Iptype::Edge,
            },
            dbg_lvl: 0,
        }
    }

    #[cfg(test)]
    pub(crate) fn test_graph(tg: tests::TestGraph, op: Optype, nparts: usize, ncon: usize) -> Self {
        tests::read_graph(tg, op, nparts, ncon)
    }

    pub fn from_csr(csr: Csr, op: Optype, nparts: usize, ncon: usize) -> Self {
        let mut ret = Self::new(op, nparts, ncon);
        ret.xadj = csr.xadj.into_vec();
        ret.adjncy = csr.adjncy.into_vec();
        ret
    }

    pub fn to_csr(&self) -> Csr {
        Csr {
            xadj: self.xadj.clone().into(),
            adjncy: self.adjncy.clone().into(),
        }
    }

    pub fn into_csr(self) -> Csr {
        Csr {
            xadj: self.xadj.into(),
            adjncy: self.adjncy.into(),
        }
    }

    pub fn new_basic(op: Optype, nparts: usize) -> Self {
        Self::new(op, nparts, 1)
    }

    /// validates degrees, returning a list of (Something about necessary edges)
    fn validate_degrees(mut deg: &[idx_t]) -> bool {
        let mut deg = Vec::from(deg);
        deg.sort_unstable();
        while let Some(d) = deg.pop() {
            if d < 0 || d as usize > deg.len() {
                return false;
            }
            if d == 0 {
                return true;
            }
            for item in deg.iter_mut().rev().take(d as usize) {
                *item -= 1;
                if *item < 0 {
                    return false;
                }
            }
        }
        true
    }

    /// approximate. uses <https://web.stanford.edu/~saberi/randgraphext.pdf>
    pub fn random_from_degrees(&mut self, mut deg: Vec<idx_t>) {
        if deg.len() <= 1 || !Self::validate_degrees(&deg) {
            return;
        }
        let mut xadj = {
            let mut a = deg.clone();
            util::make_csr(a.len() - 1, &mut a);
            a
        };
    }

    /// Allocate every vertex to have corresponding vertex degree
    pub fn with_vtx_degrees(&mut self, mut deg: Vec<idx_t>) {
        debug_assert!(deg.iter().all(|&x| x >= 0), "cannot have negative degree");
        {
            let max = deg.iter().copied().max().unwrap_or(0);
            let nvtxs = deg.len();
            debug_assert!(
                (max as usize) < nvtxs,
                "vertex has a degree of {max} but there are only {nvtxs} vertices"
            );
        }
        debug_assert!(
            deg.iter()
                .fold(0u32, |acc, &deg| acc.wrapping_add(deg as u32))
                % 2
                == 0,
            "total degree must be even"
        );
        deg.push(0);
        util::make_csr(deg.len() - 1, &mut deg);
        self.xadj = deg;
    }

    pub fn nvtxs(&self) -> usize {
        self.xadj.len() - 1
    }

    pub fn ncon(&self) -> usize {
        self.op.common_ref().map_or(1, |c| c.ncon)
    }

    fn adj(&self, vtx: idx_t) -> &[idx_t] {
        let vtx = vtx as usize;
        let istart = self.xadj[vtx] as usize;
        let iend = self.xadj[vtx + 1] as usize;
        debug_assert!(istart <= iend, "start > end ({istart} > {iend})");
        &self.adjncy[istart..iend]
    }

    fn adj_mut(&mut self, vtx: idx_t) -> &mut [idx_t] {
        let vtx = vtx as usize;
        let istart = self.xadj[vtx] as usize;
        let iend = self.xadj[vtx + 1] as usize;
        &mut self.adjncy[istart..iend]
    }

    fn insert_edge(&mut self, from_free: idx_t, free: &mut pqueue::rs::IPQueue, vtx_pair: [idx_t; 2]) {
        let [from, to] = vtx_pair;
        debug_assert_ne!(from, to);
        // eprint!("from {from} to {to} ");
        debug_assert!(from_free > 0);
        let from_off = self.adj(from).len() as idx_t - from_free;
        debug_assert!(from_off >= 0, "{from_off} < 0");
        let from_base = self.xadj[from as usize];
        let to_free = free.get(to);
        // eprintln!("({from_free} free to {to_free} free)");
        debug_assert!(to_free > 0);
        let to_off = self.adj(to).len() as idx_t - to_free;
        debug_assert!(to_off >= 0, "{to_off} < 0");
        let to_base = self.xadj[to as usize];
        self.adjncy[(from_base + from_off) as usize] = to;
        self.adjncy[(to_base + to_off) as usize] = from;
        if to_free <= 1 {
            free.delete(to)
        } else {
            free.update(to, to_free - 1);
        }
    }

    pub fn edge_list(&mut self, edges: impl IntoIterator<Item = (idx_t, idx_t)>) {
        fn inner(this: &mut GraphBuilder, mut edges: Vec<(idx_t, idx_t)>) {
            assert!(this.xadj.is_empty(), "already set edges");
            edges.retain(|&(l, r)| l != r);
            let part = edges.len();
            edges.extend_from_within(..);
            edges[part..].iter_mut().for_each(|e| *e = (e.1, e.0));
            edges.sort_unstable();
            edges.dedup();
            edges.sort_unstable_by_key(|&(l, r)| (r, l));
            edges.dedup();

            let max = *edges
                .iter()
                .flat_map(|(l, r)| [l, r])
                .max()
                .expect("empty edge list") as usize;
            let mut degress = vec![0; max + 1];
            for &(l, r) in &edges {
                degress[l as usize] += 1;
                // degress[r as usize] += 1;
            }
            // TODO: this clone is unnecessary, and can be done by modifying xadj and reverting
            let mut xadj = degress.clone();
            xadj.push(idx_t::MIN);
            make_csr(max as usize + 1, &mut xadj);
            let mut adjncy = vec![0; xadj[max as usize + 1] as usize];
            for &(l, r) in &edges {
                // only need one side since we duplicate
                let l_idx = (degress[l as usize] + xadj[l as usize]) as usize - 1;
                // let r_idx = (degress[r as usize] + xadj[r as usize]) as usize - 1;
                adjncy[l_idx] = r;
                // adjncy[r_idx] = l;
                degress[l as usize] -= 1;
                // degress[r as usize] -= 1;
            }
            this.adjncy = adjncy;
            this.xadj = xadj;
        }
        inner(self, edges.into_iter().collect())
    }

    /// using vertex degrees specified by `with_vtx_degrees`, make random edges
    pub fn random_edges(&mut self) {
        let mut rng = fastrand::Rng::new();
        let mut remaining_vtx = pqueue::rs::IPQueue::new(self.nvtxs());

        for (i, deg) in util::from_csr_iter(&self.xadj).enumerate() {
            if deg > 0 {
                remaining_vtx.insert(i as idx_t, deg);
            }
        }
        self.adjncy
            .resize(self.xadj.last().copied().unwrap_or(0) as usize, -1);
        while let Some((cnt, vtx)) = remaining_vtx.pop_node() {
            let mut list = remaining_vtx.clone().to_dal();
            eprintln!("assembling vertex {vtx}");
            for i in 0..cnt {
                debug_assert!(remaining_vtx.length() >= 1);
                eprintln!("list: {list:?}");
                let add = rng.usize(0..list.len());
                eprint!("adding idx {add} ");
                let add = list.as_slices().1[add];
                eprintln!("with val {add}");
                debug_assert_ne!(add, -1);
                self.insert_edge(cnt - i, &mut remaining_vtx, [vtx as idx_t, add]);
                list.remove(add);
            }
            debug_assert_eq!(remaining_vtx.try_get(vtx as idx_t), None);
        }
        debug_assert_eq!(self.adjncy.iter().position(|&a| a == -1), None);
    }

    pub fn random_tpwgts(&mut self) {
        let comm = self.op.common().unwrap();
        let mut rng = Rng::new();
        let mut one_con = || {
            let base = Vec::from_iter((0..comm.nparts).map(|_| rng.f32()));
            let total: f32 = base.iter().sum();
            base.into_iter().map(move |b| b / total)
        };
        let mut tpwgts = comm
            .tpwgts
            .take()
            .unwrap_or(Vec::with_capacity(comm.ncon * comm.nparts));
        tpwgts.resize(comm.ncon * comm.nparts, 0.0);
        for con in 0..comm.ncon {
            tpwgts
                .iter_mut()
                .skip(con)
                .step_by(comm.ncon)
                .zip(one_con())
                .for_each(|(tpwgt, wgt)| *tpwgt = wgt);
        }
        comm.tpwgts = Some(tpwgts);
    }

    pub fn random_ubvec(&mut self) {
        let mut rng = Rng::new();
        let comm = self.op.common().unwrap();
        let mut ubvec = comm.ubvec.take().unwrap_or(Vec::with_capacity(comm.ncon));
        ubvec.clear();
        ubvec.extend((0..comm.ncon).map(|_| rng.f32() / 15.0 + 1.001));
        comm.ubvec = Some(ubvec);
    }

    pub fn random_vwgt(&mut self) {
        let mut rng = Rng::new();
        let mut vwgt = self
            .vwgt
            .take()
            .unwrap_or(Vec::with_capacity(self.ncon() * self.nvtxs()));
        vwgt.clear();
        vwgt.extend((0..(self.ncon() * self.nvtxs())).map(|_| rng.u32(1..20) as idx_t));
        self.vwgt = Some(vwgt);
    }

    pub fn random_vsize(&mut self) {
        let mut rng = Rng::new();
        let nvtxs = self.nvtxs();
        match &mut self.op {
            GraphOpSettings::Kmetis { objtype: KwayObjective::Cut { adjwgt: Some(_) }, .. } => panic!("cannot set vsize -- already set adjwgt"),
            GraphOpSettings::Ometis { .. } | GraphOpSettings::Pmetis { .. } => panic!("cannot set vsize unless kway"),
            _ => ()
        }
        let mut vsize = Vec::with_capacity(nvtxs);
        vsize.extend((0..nvtxs).map(|_| rng.u32(1..20) as idx_t));
        match &mut self.op {
            GraphOpSettings::Kmetis { objtype, .. } => *objtype = KwayObjective::Vol { vsize: Some(vsize) },
            GraphOpSettings::Ometis { .. } | GraphOpSettings::Pmetis { .. } => panic!("cannot set vsize unless kway"),
            _ => ()
        }
    }

    fn vtxs(&self) -> impl Iterator<Item = usize> {
        0..self.nvtxs()
    }

    pub fn random_adjwgt(&mut self) {
        let mut rng = Rng::new();
        let nedges = self.adjncy.len();
        let orig = match &mut self.op {
            GraphOpSettings::Pmetis { adjwgt, .. } => adjwgt,
            GraphOpSettings::Kmetis { objtype: KwayObjective::Cut { adjwgt }, .. } => adjwgt,
            GraphOpSettings::Kmetis { objtype: KwayObjective::Vol { vsize: Some(_) }, .. } => panic!("already set vsize - can't also set adjwgt"),
            GraphOpSettings::Kmetis { .. } => &mut None,
            GraphOpSettings::Ometis { .. } => panic!("can't set adjwgt in ometis"),
        };
        let mut adjwgt = orig.as_mut().take().map_or_else(|| Vec::with_capacity(nedges), |a| std::mem::take(a));
        adjwgt.clear();
        adjwgt.extend((0..self.adjncy.len()).map(|_| rng.u32(1..50) as idx_t));
        for vtx in self.vtxs() {
            for (pos, &other) in self.adj(vtx as idx_t).iter().enumerate() {
                self.adj(other)
                    .iter()
                    .position(|&v| v == vtx as idx_t)
                    .map(|iadj| {
                        adjwgt[pos + self.xadj[vtx] as usize] =
                            adjwgt[self.xadj[other as usize] as usize + iadj]
                    });
            }
        }
        match &mut self.op {
            GraphOpSettings::Pmetis { adjwgt: dst, .. } => *dst = Some(adjwgt),
            GraphOpSettings::Kmetis { objtype, .. } => *objtype = KwayObjective::Cut { adjwgt: Some(adjwgt) },
            _ => unreachable!()
        }
    }

    pub fn set_contig(&mut self, contig: bool) {
        match &mut self.op {
            GraphOpSettings::Kmetis { contig: x, .. } => *x = contig,
            _ => panic!("Can only specify contig for kmetis")
        }
    }
    
    pub fn set_minconn(&mut self, minconn: bool) {
        match &mut self.op {
            GraphOpSettings::Kmetis { minconn: x, .. } => *x = minconn,
            _ => panic!("Can only specify minconn for kmetis")
        }
    }

    pub fn set_objective(&mut self, objtype: Objtype) {
        match (&mut self.op, objtype) {
            (GraphOpSettings::Kmetis {  .. }, Objtype::Node) => panic!("cannot set node objtype"),
            (GraphOpSettings::Kmetis { objtype: KwayObjective::Vol { .. }, .. }, Objtype::Vol) => (),
            (GraphOpSettings::Kmetis { objtype: KwayObjective::Cut { .. }, .. }, Objtype::Cut) => (),
            (GraphOpSettings::Kmetis { objtype: objtype @ KwayObjective::Vol { vsize: None }, .. }, Objtype::Cut) => *objtype = KwayObjective::Cut { adjwgt: None },
            (GraphOpSettings::Kmetis { objtype: objtype @ KwayObjective::Cut { adjwgt: None }, .. }, Objtype::Vol) => *objtype = KwayObjective::Vol { vsize: None },
            (GraphOpSettings::Kmetis {  .. }, _) => panic!("trying to set objtype to {objtype:?}, but we already populated in another objtype"),
            (GraphOpSettings::Pmetis { .. } 
                | GraphOpSettings::Ometis { .. }, _) => panic!("setting objective does nothing on pmetis or ometis"),
        }
    }

    pub fn set_ccorder(&mut self, cc: bool) {
        match &mut self.op {
            GraphOpSettings::Pmetis { ..} |
            GraphOpSettings::Kmetis { .. } => panic!("cannot set ccorder on pmetis or kmetis"),
            GraphOpSettings::Ometis { ccorder, .. } => *ccorder = cc,
        }
    }

    pub fn set_nseps(&mut self, n: idx_t) {
        match &mut self.op {
            GraphOpSettings::Pmetis { ..} |
            GraphOpSettings::Kmetis { .. } => panic!("cannot set nseps on pmetis or kmetis"),
            GraphOpSettings::Ometis { nseps, .. } => *nseps = n,
        }
    }

    /// sets `compress` for ometis (default `true`)
    pub fn set_compress(&mut self, do_compress: bool) {
        // TODO: make ometis tests without compress
        // TODO: write tests for compression routines
        match &mut self.op {
            GraphOpSettings::Pmetis { ..} |
            GraphOpSettings::Kmetis { .. } => panic!("cannot set nseps on pmetis or kmetis"),
            GraphOpSettings::Ometis { compress, .. } => *compress = do_compress,
        }
    }

    /// call [`parmetis::METIS_ComputeVertexSeparator`] instead of [`ometis::METIS_NodeND`]
    pub fn compute_vertex_separator(&mut self, do_vtx_sep: bool) {
        // TODO: make ometis tests without compress
        // TODO: write tests for compression routines
        match &mut self.op {
            GraphOpSettings::Pmetis { ..} |
            GraphOpSettings::Kmetis { .. } => panic!("cannot set nseps on pmetis or kmetis"),
            GraphOpSettings::Ometis { compute_vertex_separator, .. } => *compute_vertex_separator = do_vtx_sep,
        }
    }

    pub fn set_pfactor(&mut self, p: idx_t) {
        match &mut self.op {
            GraphOpSettings::Pmetis { ..} |
            GraphOpSettings::Kmetis { .. } => panic!("cannot set nseps on pmetis or kmetis"),
            GraphOpSettings::Ometis { pfactor, .. } => *pfactor = p,
        }
    }

    // pub fn set_refine_type(&mut self, rtype: Rtype) {
    //     assert_eq!(
    //         self.op,
    //         Optype::Ometis,
    //         "setting the refine type is only availiable on ometis (despite what the manual says)"
    //     );
    //
    //     self.refine_type = Some(rtype);
    // }

    pub fn set_edge_match(&mut self, ctype: Ctype) {
        self.edge_match = ctype;
    }

    pub fn set_initial_part_strategy(&mut self, iptype: Iptype) {
        self.initial_part = iptype;
    }

    pub fn set_seed(&mut self, seed: idx_t) {
        self.seed = seed;
    }

    pub fn enable_dbg(&mut self, lvl: DbgLvl) {
        self.dbg_lvl |= lvl as u32
    }

    pub fn set_from_options_arr(&mut self, options: &[idx_t; METIS_NOPTIONS as usize]) {
        macro_rules! match_setting {
            ($options:ident[$idx:ident] == $base:ident::{$($var:ident),* $(,)?}) => {
                {
                    let opt = $options[$idx as usize];
                    if opt == -1 {
                        None
                    } $(
                    else if opt == $base::$var as idx_t {
                        Some($base::$var)
                    }
                    )*
                    else {
                        panic!("unknown option variant {} for option {} ", opt, stringify!($idx))
                    }
                }
            };
        }
        self.seed = options[METIS_OPTION_SEED as usize];
        if let Some(iptype) =
            match_setting!(options[METIS_OPTION_IPTYPE] == Iptype::{Rb, Random, Grow, Edge, Node})
        {
            self.initial_part = iptype;
        }
        if let Some(ctype) = match_setting!(options[METIS_OPTION_CTYPE] == Ctype::{Shem, Rm}) {
            self.edge_match = ctype;
        }
        if let Some(rtype) = match_setting!(options[METIS_OPTION_RTYPE] == Rtype::{Fm, Greedy, Sep1Sided, Sep2Sided}) {
            match &mut self.op {
                GraphOpSettings::Ometis { rtype: dst, .. } => *dst = match rtype {
                    Rtype::Greedy |
                    Rtype::Fm => panic!("can't set rtype of ometis to {rtype:?}"),
                    Rtype::Sep1Sided => OmetisRtype::Sep1Sided,
                    Rtype::Sep2Sided => OmetisRtype::Sep2Sided,
                },
                _ => panic!("can't set rtype outside of ometis")
            }
        }
    }

    pub fn call(&mut self) -> Result<(idx_t, Vec<i32>), ()> {
        assert_eq!(self.adjncy.len(), *self.xadj.last().unwrap_or(&0) as usize);
        let mut nvtxs = self.nvtxs() as idx_t;
        let mut ncon = self.ncon() as idx_t;
        let mut part = vec![0; self.nvtxs()];
        let mut objval = 0;
        let mut options = [-1; METIS_NOPTIONS as usize];
        options[METIS_OPTION_CTYPE as usize] = self.edge_match as idx_t;
        options[METIS_OPTION_IPTYPE as usize] = self.initial_part as idx_t;
        options[METIS_OPTION_ONDISK as usize] = 1;
        options[METIS_OPTION_SEED as usize] = self.seed;
        options[METIS_OPTION_DBGLVL as usize] = self.dbg_lvl as idx_t;

        let res = match &mut self.op {
            GraphOpSettings::Pmetis { common: GraphPartSettings { ncuts, nparts, ncon, ubvec, tpwgts }, adjwgt } => unsafe {
                let mut nparts = *nparts as idx_t;
                let mut ncon = *ncon as idx_t;
                pmetis::METIS_PartGraphRecursive(
                    &mut nvtxs,
                    &mut ncon,
                    self.xadj.as_mut_ptr(),
                    self.adjncy.as_mut_ptr(),
                    vec_ptr(&mut self.vwgt),
                    std::ptr::null_mut(),
                    vec_ptr(adjwgt),
                    &mut nparts,
                    vec_ptr(tpwgts),
                    vec_ptr(ubvec),
                    options.as_mut_ptr(),
                    &mut objval,
                    part.as_mut_ptr(),
                )
            },
            GraphOpSettings::Kmetis { common: GraphPartSettings { ncuts, nparts, ncon, ubvec, tpwgts }, objtype, minconn, contig, niparts } => unsafe {
                let mut nparts = *nparts as idx_t;
                let mut ncon = *ncon as idx_t;
                let vsize;
                let adjwgt;
                match objtype {
                    KwayObjective::Vol { vsize: mv } => {
                        options[METIS_OPTION_OBJTYPE as usize] = Objtype::Vol as idx_t;
                        vsize = vec_ptr(mv);
                        adjwgt = std::ptr::null_mut();
                    },
                    KwayObjective::Cut { adjwgt: aw } => {
                        options[METIS_OPTION_OBJTYPE as usize] = Objtype::Cut as idx_t;
                        vsize = std::ptr::null_mut();
                        adjwgt = vec_ptr(aw);
                    },
                };
                options[METIS_OPTION_CONTIG as usize] = *contig as idx_t;
                options[METIS_OPTION_MINCONN as usize] = *minconn as idx_t;
                kmetis::METIS_PartGraphKway(
                    &mut nvtxs,
                    &mut ncon,
                    self.xadj.as_mut_ptr(),
                    self.adjncy.as_mut_ptr(),
                    vec_ptr(&mut self.vwgt),
                    vsize,
                    adjwgt,
                    &mut nparts,
                    vec_ptr(tpwgts),
                    vec_ptr(ubvec),
                    options.as_mut_ptr(),
                    &mut objval,
                    part.as_mut_ptr(),
                )
            },
            GraphOpSettings::Ometis { nseps, pfactor, ccorder, rtype, compress, compute_vertex_separator } => unsafe {
                options[METIS_OPTION_COMPRESS as usize] = *compress as idx_t;
                options[METIS_OPTION_CCORDER as usize] = *ccorder as idx_t;
                options[METIS_OPTION_PFACTOR as usize] = *pfactor as idx_t;
                options[METIS_OPTION_RTYPE as usize] = match rtype {
                    OmetisRtype::Sep1Sided => Rtype::Sep1Sided,
                    OmetisRtype::Sep2Sided => Rtype::Sep2Sided,
                } as idx_t;
                options[METIS_OPTION_NSEPS as usize] = *nseps as idx_t;
                let mut iperm = vec![0; nvtxs as usize];
                let ret;
                if *compute_vertex_separator {
                    ret = parmetis::METIS_ComputeVertexSeparator(
                        &mut nvtxs,
                        self.xadj.as_mut_ptr(),
                        self.adjncy.as_mut_ptr(),
                        vec_ptr(&mut self.vwgt),
                        options.as_mut_ptr(),
                        &mut objval,
                        part.as_mut_ptr()
                    );
                } else {
                    ret = ometis::METIS_NodeND(
                        &mut nvtxs,
                        self.xadj.as_mut_ptr(),
                        self.adjncy.as_mut_ptr(),
                        vec_ptr(&mut self.vwgt),
                        options.as_mut_ptr(),
                        part.as_mut_ptr(), // use part as permutation vector
                        iperm.as_mut_ptr()
                    );
                    part.extend_from_slice(&iperm);
                }
                ret
            },
        };
        if res != METIS_OK {
            Err(())
        } else {
            Ok((objval, part))
        }
    }

    /// Calls [`parmetis::METIS_NodeNDP`]. Must be set to [`Optype::Ometis`]. Returns `(perm,
    /// iperm, sizes)`
    pub fn call_ndp(&mut self, npes: idx_t) -> Result<(Vec<idx_t>, Vec<idx_t>, Vec<idx_t>), ()> {
        assert_eq!(self.adjncy.len(), *self.xadj.last().unwrap_or(&0) as usize);
        let mut nvtxs = self.nvtxs() as idx_t;
        let mut ncon = self.ncon() as idx_t;
        let mut part = vec![0; self.nvtxs()];
        let mut objval = 0;
        let mut options = [-1; METIS_NOPTIONS as usize];
        options[METIS_OPTION_CTYPE as usize] = self.edge_match as idx_t;
        options[METIS_OPTION_IPTYPE as usize] = self.initial_part as idx_t;
        options[METIS_OPTION_ONDISK as usize] = 1;
        options[METIS_OPTION_SEED as usize] = self.seed;
        options[METIS_OPTION_DBGLVL as usize] = self.dbg_lvl as idx_t;

        let Self { dbg_lvl, op: GraphOpSettings::Ometis {
            nseps,
            pfactor,
            ccorder,
            rtype,
            compress,
            compute_vertex_separator
        }, seed, xadj, adjncy, vwgt, edge_match, initial_part } = self else { panic!("not ometis") };
        assert!(!*compute_vertex_separator, "cannot use compute vertex separator when calling NodeNDP");
        options[METIS_OPTION_COMPRESS as usize] = *compress as idx_t;
        options[METIS_OPTION_CCORDER as usize] = *ccorder as idx_t;
        options[METIS_OPTION_PFACTOR as usize] = *pfactor as idx_t;
        options[METIS_OPTION_RTYPE as usize] = match rtype {
            OmetisRtype::Sep1Sided => Rtype::Sep1Sided,
            OmetisRtype::Sep2Sided => Rtype::Sep2Sided,
        } as idx_t;
        options[METIS_OPTION_NSEPS as usize] = *nseps as idx_t;
        let mut perm = vec![0; nvtxs as usize];
        let mut iperm = vec![0; nvtxs as usize];
        let mut sizes = vec![0; 2 * npes as usize - 1];

        let res = unsafe {
            parmetis::METIS_NodeNDP(nvtxs, xadj.as_mut_ptr(), adjncy.as_mut_ptr(), vec_ptr(vwgt),
                npes, options.as_ptr(), perm.as_mut_ptr(), iperm.as_mut_ptr(), sizes.as_mut_ptr())
        };
        if res != METIS_OK {
            Err(())
        } else {
            Ok((perm, iperm, sizes))
        }
    }

    pub fn write_graph<W: io::Write>(&self, mut w: W) -> io::Result<()> {
        // command line invocation
        {
            writeln!(w, "% partition with this command:")?;
            let command = match self.op.op() {
                Optype::Kmetis => "gpmetis -ptype=kway",
                Optype::Ometis => "ndmetis",
                Optype::Pmetis => "gpmetis -ptype=rb",
            };
            write!(w, "% {command}")?;

            let ctype = match self.edge_match {
                Ctype::Rm => "-ctype=rm",
                Ctype::Shem => "-ctype=shem",
            };
            write!(w, " {ctype}")?;

            match &self.op {
                GraphOpSettings::Pmetis { common: GraphPartSettings { ncuts, nparts, ncon, ..}, .. } => todo!(),
                GraphOpSettings::Kmetis { common, objtype, minconn, contig, niparts } => todo!(),
                GraphOpSettings::Ometis { nseps, pfactor, ccorder, rtype, compress, compute_vertex_separator } => todo!(),
            }

            {
                let iptype = match self.initial_part {
                    Iptype::Grow => "-iptype=grow",
                    Iptype::Random => "-iptype=random",
                    Iptype::Edge => "-iptype=edge",
                    Iptype::Node => "-iptype=node",
                    Iptype::Rb => unreachable!(),
                };
                write!(w, " {iptype}")?;
            }
            if let GraphOpSettings::Kmetis { objtype, contig, .. } = &self.op {
                let objtype = match objtype.objtype() {
                    Objtype::Cut => "-objtype=cut",
                    Objtype::Vol => "-objtype=vol",
                    Objtype::Node => unreachable!(),
                };
                write!(w, " {objtype}")?;

                if *contig {
                    write!(w, " -contig")?;
                }
            }

            write!(w, " -seed={}", self.seed)?;

            if let Some(comm) = self.op.common_ref() {
                if let Some(ubvec) = comm.ubvec.as_deref() {
                    write!(w, r#" -ubvec="{:.4}""#, space_sep(ubvec))?;
                }
                writeln!(w, " <GRAPH_FILE> {}", comm.nparts)?;
            } else {
                debug_assert!(self.op.op().is_ometis());
                writeln!(w, " <GRAPH_FILE>")?;
            }
        }

        // header line
        let has_vsize = self.vsize().is_some();
        let has_vwgt = self.vwgt.is_some() || self.ncon() > 1;
        let has_adjwgt = self.adjwgt().is_some();
        let fmt = [has_vsize, has_vwgt, has_adjwgt].map(|b| b as u8 + b'0');
        writeln!(w, "{nvtxs} {nedges} {fmt} {ncon}", 
            nvtxs = self.nvtxs(),
            nedges = self.adjncy.len() / 2,
            fmt = std::str::from_utf8(&fmt).unwrap(),
            ncon = self.ncon()
        )?;

        // main graph structure
        for vtx in self.vtxs() {
            // FORMAT: size? wgt{ncon} (edge, edgewgt?)* 
            let params = (self.vsize().map(|vs| vs[vtx]).into_iter())
            .chain({
                    self.vwgt.iter().flat_map(|vwgt| &vwgt[(vtx * self.ncon())..((vtx + 1) * self.ncon())]).copied()
                        .chain(std::iter::repeat(1))
                        .take(self.ncon())
                })
                .chain({
                    let range = (self.xadj[vtx] as usize)..(self.xadj[vtx + 1] as usize);
                    let edges = self.adjncy[range.clone()].iter().copied().map(|x| Some(x + 1));
                    let adjwgts = self.adjwgt().map(|a| &a[range]).into_iter().flatten().copied().map(Some).chain(std::iter::repeat(None));
                    edges.zip(adjwgts).flat_map(|(e, w)| [e, w]).flatten()
                })
            ;
            writeln!(w, "{}", space_sep(params))?;
        }

        Ok(())
    }

    pub fn read_graph(
        mut f: impl BufRead,
        op: Optype,
        nparts: usize,
        ncon: usize,
    ) -> Result<GraphBuilder, Box<dyn Error>> {
        Self::read_graph_cfg_index(f, op, nparts, ncon, false)
    }

    /// read a graph from a file in a simplified version of the format specified in the manual
    pub fn read_graph_cfg_index(
        mut f: impl BufRead,
        op: Optype,
        nparts: usize,
        ncon: usize,
        one_indexed: bool,
    ) -> Result<GraphBuilder, Box<dyn Error>> {
        let mut buf = String::with_capacity(80);
        f.read_line(&mut buf)?;

        let mut xadj = Vec::<idx_t>::new();
        let mut adjncy = Vec::<idx_t>::new();

        let mut x = 0;
        buf.clear();
        while 0 != f.read_line(&mut buf)? {
            if buf.trim_start().starts_with('#') || buf.trim().is_empty() {
                continue;
            }

            xadj.push(x);

            let split = buf.split_whitespace();
            if one_indexed {
                for a in split {
                    adjncy.push(a.parse::<idx_t>()?.checked_sub(1)
                        .expect("0 indexed but one indexed is set"));
                    x += 1;
                }
            } else {
                for a in split {
                    adjncy.push(a.parse::<idx_t>()?);
                    x += 1;
                }
            }

            buf.clear();
        }
        xadj.push(x);

        for (i, (start, end)) in xadj.windows(2).map(|w| (w[0] as usize, w[1] as usize)).enumerate() {
            assert!(start < adjncy.len());
            assert!(end <= adjncy.len());
            assert!(start < end);
            // dbg!(&adjncy[start..end]);
            for &j in &adjncy[start..end] {
                assert!(j >= 0, "no negatives");
                assert_ne!(i, j as usize, "no self loops");
                assert!(j < xadj.len() as idx_t - 1, "adj in bounds");
            }
        }

        Ok(Self {
            xadj,
            adjncy,
            ..Self::new(op, nparts, ncon)
        })
    }

    pub fn vsize(&self) -> Option<&[idx_t]> {
        if let GraphOpSettings::Kmetis { common, objtype: KwayObjective::Vol { vsize }, minconn, contig, niparts } = &self.op {
            vsize.as_deref()
        } else {
            None
        }
    }
    pub fn vwgt(&self) -> Option<&[idx_t]> {
        self.vwgt.as_deref()
    }
    pub fn adjwgt(&self) -> Option<&[idx_t]> {
        match &self.op {
            GraphOpSettings::Kmetis { objtype: KwayObjective::Cut { adjwgt }, ..} |
            GraphOpSettings::Pmetis { adjwgt, .. } => adjwgt.as_deref(),
            _ => None
        }
    }

    /// adapted from mtest.c: VerifyPart
    pub fn verify_part(&self, _objval: idx_t, part: &[idx_t]) {
        assert_eq!(
            self.xadj.len() - 1,
            part.len(),
            "part is the partition that each vertex goes to"
        );

        let Some(comm) = self.op.common_ref() else {
            // TODO: verify sepnd
            return
        };

        let mut pwgts = vec![0; comm.nparts];

        assert_eq!(
            *part.iter().max().unwrap_or(&0),
            comm.nparts as idx_t - 1,
            "total number of partitions eq to nparts"
        );

        let mut _cut = 0;
        for i in 0..self.nvtxs() {
            pwgts[part[i] as usize] += self.vwgt.as_ref().map(|v| v[i]).unwrap_or(1);
            for j in self.xadj[i]..self.xadj[i + 1] {
                if part[i] != part[self.adjncy[j as usize] as usize] {
                    _cut += self.adjwgt().map(|v| v[j as usize]).unwrap_or(1);
                }
            }
        }

        // eprintln!("todo: make this work always. This assumes edgecut but we call this for vol too");
        // assert_eq!(
        //     cut,
        //     2 * objval,
        //     "objval should be edgecut, and the calculated cut should be double it"
        // );
        let actual = (comm.nparts as idx_t * pwgts.iter().max().unwrap()) as f64;
        // should be 1.10, but it's annoying
        let expected = 1.13 * pwgts.iter().sum::<idx_t>() as f64;
        assert!(
            // (nparts * pwgts[iargmax(nparts, pwgts)]) as f64
            actual <= expected,
            "actual: {actual:.1}, expected: {expected:.1}.\n\tThis assert spuriously fails sometimes - rerun tests."
        );
    }
}

/// utility for writing space-separated lists
fn space_sep<I>(iter: I) -> impl std::fmt::Display + use<I>
where 
    I: IntoIterator + Clone,
    I::Item: std::fmt::Display
{
    use std::fmt::{Display, Write};

    struct D<I>(I);
    impl<I> Display for D<I>
    where 
        I: IntoIterator + Clone,
        I::Item: std::fmt::Display
    {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let mut it = self.0.clone().into_iter();
            {
                let Some(first) = it.next() else { return Ok(()) };
                first.fmt(f)?;
            }
            for item in it {
                f.write_char(' ')?;
                item.fmt(f)?;
            }
            Ok(())
        }
    }
    D(iter)
}
