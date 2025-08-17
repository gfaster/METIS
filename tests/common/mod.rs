use std::{fmt::{self, Debug}, io::{self, Read, Write}, os::fd::{AsFd, AsRawFd, BorrowedFd}, path::PathBuf, process::{Child, Command, Output, Stdio}, sync::LazyLock};

pub use metis::bindings::{
    Optype as Ptype,
    Objtype,
    Ctype,
    Iptype,
    Rtype,
    idx_t,
    real_t,
};


const PRINT_CMDS: bool = true;


#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ver {
    Rust,
    C
}

impl Ver {
    #![expect(dead_code)]

    /// Returns `true` if the ver is [`C`].
    ///
    /// [`C`]: Ver::C
    #[must_use]
    pub fn is_c(&self) -> bool {
        matches!(self, Self::C)
    }

    /// Returns `true` if the ver is [`Rust`].
    ///
    /// [`Rust`]: Ver::Rust
    #[must_use]
    pub fn is_rust(&self) -> bool {
        matches!(self, Self::Rust)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TestGraph {
    Elt4,
    Youtube,
    Webbase2004,
    WebSpam,
    Luxembourg,
    Orani,
    CitDblp,
}

impl TestGraph {
    #![allow(dead_code)]
    #![deny(clippy::wildcard_enum_match_arm)]

    pub const fn is_smallish(self) -> bool {
        match self {
            TestGraph::Elt4 |
            TestGraph::Webbase2004 |
            TestGraph::WebSpam |
            TestGraph::CitDblp => true,
            TestGraph::Orani |
            TestGraph::Youtube |
            TestGraph::Luxembourg => false,
        }
    }

    pub fn is_contiguous(self) -> bool {
        match self {
            TestGraph::Elt4
            | TestGraph::Youtube
            | TestGraph::Webbase2004
            | TestGraph::WebSpam
            | TestGraph::Luxembourg
            | TestGraph::Orani => true,
            TestGraph::CitDblp => false,
        }
    }

    const fn file(self) -> &'static str {
        match self {
            TestGraph::Elt4 => "graphs/4elt.graph",
            TestGraph::Youtube => "graphs/soc-youtube.graph",
            TestGraph::Webbase2004 => "graphs/web-webbase-2001.graph",
            TestGraph::WebSpam => "graphs/web-spam.graph",
            TestGraph::Luxembourg => "graphs/road-luxembourg-osm.graph",
            TestGraph::Orani => "graphs/econ-orani678.graph",
            TestGraph::CitDblp => "graphs/cit-DBLP.graph",
        }
    }

    pub fn test_suite() -> impl Iterator<Item = Self> {
        static DO_BIG: LazyLock<bool> =
            LazyLock::new(|| std::env::var_os("DO_BIG") == Some("1".into()));
        let do_big = *DO_BIG;
        Self::ALL
            .into_iter()
            .filter(move |g| g.is_smallish() || do_big)
    }

    const COUNT: usize = 7;

    const ALL: [Self; Self::COUNT] = [
        Self::Elt4,
        Self::Youtube,
        Self::Webbase2004,
        Self::WebSpam,
        Self::Luxembourg,
        Self::Orani,
        Self::CitDblp,
    ];
}

impl From<TestGraph> for PathBuf {
    fn from(value: TestGraph) -> Self {
        value.file().into()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bin {
    GpMetis,
}
pub use Bin::*;

impl Bin {
    fn exe(self, ver: Ver) -> &'static str {
        match (self, ver) {
            (GpMetis, Ver::Rust) => env!("CARGO_BIN_EXE_gpmetis"),
            (GpMetis, Ver::C) => concat!(env!("METIS_NORM"), "/bin/gpmetis"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Params {
    pub bin: Bin,
    pub filename: PathBuf,
    pub nparts: idx_t,

    pub ptype: Option<Ptype>,
    pub objtype: Option<Objtype>,
    pub ctype: Option<Ctype>,
    pub iptype: Option<Iptype>,
    pub rtype: Option<Rtype>,

    pub no2hop: bool,
    pub minconn: bool,
    pub contig: bool,

    pub ondisk: bool,

    pub dropedges: bool,

    pub nooutput: bool,

    pub balance: bool,
    pub ncuts: Option<idx_t>,
    pub niter: Option<idx_t>,
    pub niparts: Option<idx_t>,

    // mesh params
    // pub gtype: idx_t,
    // pub ncommon: idx_t,

    pub seed: Option<idx_t>,
    pub dbglvl: Option<metis::moptions_et>,


    // nd params
    pub nseps: Option<idx_t>,
    pub ufactor: Option<idx_t>,
    pub pfactor: Option<idx_t>,
    pub compress: bool,
    pub ccorder: bool,

    pub ubvec: Vec<real_t>,

    // pub numflag: idx_t,
    pub tpwgts: Vec<real_t>,

    pub overrides: &'static str,

    pub help: bool,
}

impl Params {
    pub fn new(bin: Bin, file: impl Into<PathBuf>, nparts: idx_t) -> Self {
        Self::new_inner(bin.into(), file.into(), nparts)
    }

    fn new_inner(bin: Bin, filename: PathBuf, nparts: idx_t) -> Self {
        Self {
            bin,
            filename,
            nparts,
            ptype: Default::default(),
            objtype: Default::default(),
            ctype: Default::default(),
            iptype: Default::default(),
            rtype: Default::default(),
            no2hop: Default::default(),
            minconn: Default::default(),
            contig: Default::default(),
            ondisk: Default::default(),
            dropedges: Default::default(),
            nooutput: Default::default(),
            balance: Default::default(),
            ncuts: Default::default(),
            niter: Default::default(),
            niparts: Default::default(),
            seed: Default::default(),
            dbglvl: Default::default(),
            nseps: Default::default(),
            ufactor: Default::default(),
            pfactor: Default::default(),
            compress: Default::default(),
            ccorder: Default::default(),
            ubvec: Default::default(),
            tpwgts: Default::default(),
            help: Default::default(),
            overrides: Default::default(),
        }
    }

    pub fn call(&self, ver: Ver) -> io::Result<(String, String, Vec<idx_t>)> {
        let (output, part) = self.call_no_success(ver)?;

        let stdout = String::from_utf8(output.stdout).unwrap();
        let stderr = String::from_utf8(output.stderr).unwrap();

        if !output.status.success() {
            println!("{ver:?} {}", output.status);
            println!("\n====== {ver:?} stdout =======\n");
            print!("{stdout}");
            println!("\n====== {ver:?} stderr =======\n");
            print!("{stderr}");
            panic!("{}", output.status);
        }

        Ok((stdout, stderr, part))
    }

    pub fn call_no_success(&self, ver: Ver) -> io::Result<(Output, Vec<idx_t>)> {
        let Self {
            bin,
            filename,
            nparts,
            ptype,
            objtype,
            ctype,
            iptype,
            rtype,
            no2hop,
            minconn,
            contig,
            ondisk,
            dropedges,
            nooutput,
            balance,
            ncuts,
            niter,
            niparts,
            seed,
            dbglvl,
            nseps,
            ufactor,
            pfactor,
            compress,
            ccorder,
            ubvec,
            tpwgts,
            help,
            overrides,
        } = self;

        println!("Calling {filename:?} with {ver:?}");

        let (pr, pw) = io::pipe().unwrap();
        
        let mut cmd = Command::new(bin.exe(ver));
        cmd.arg(filename);
        cmd.arg(format!("{nparts}"));

        cmd.arg(format!("--outfile=/proc/self/fd/{}", pw.as_raw_fd()));

        if !overrides.is_empty() && ver.is_rust() {
            cmd.env("METIS_OVERRIDE_SYMS", *overrides);
        }

        let tpr;
        let tpw;
        let tp_buf;
        if !tpwgts.is_empty() {
            use std::fmt::Write;
            let (tpr_o, tpw_o) = io::pipe().unwrap();
            unset_cloexec(&tpr_o).unwrap();
            let ncon = tpwgts.len() / *nparts as usize;
            let mut buf = String::new();
            'done: for p in 0..*nparts as usize {
                for c in 0..ncon {
                    let Some(tpwgt) = tpwgts.get(p * ncon + c) else {
                        break 'done
                    };
                    writeln!(buf, "{p} : {c} = {tpwgt}").unwrap();
                }
            }
            tp_buf = buf;

            cmd.arg(format!("--tpwgts=/proc/self/fd/{}", tpr_o.as_raw_fd()));
            tpr = Some(tpr_o);
            tpw = Some(tpw_o);
        } else {
            tp_buf = String::new();
            tpr = None;
            tpw = None;
        }

        if let Some(ptype) = ptype {
            cmd.arg(match ptype {
                Ptype::Kmetis => "-ptype=kway",
                Ptype::Pmetis => "-ptype=rb",
                Ptype::Ometis => unimplemented!(),
            });
        }

        if let Some(ctype) = ctype {
            cmd.arg(match ctype {
                Ctype::Rm => "-ctype=rm",
                Ctype::Shem => "-ctype=shem",
            });
        }

        cmd.args([
            contig.then_some("-contig"),
            minconn.then_some("-minconn"),
            balance.then_some("-balance"),
            ccorder.then_some("-ccorder"),
            compress.then_some("-compress"),
            ondisk.then_some("-ondisk"),
            nooutput.then_some("-nooutput"),
            no2hop.then_some("-no2hop"),
            dropedges.then_some("-dropedges"),
            help.then_some("-help"),
        ].iter().flatten());


        if let Some(ufactor) = ufactor {
            cmd.arg(format!("-ufactor={ufactor}"));
        }

        if let Some(objtype) = objtype {
            cmd.arg(match objtype {
                Objtype::Cut => "-objtype=cut",
                Objtype::Vol => "-objtype=vol",
                Objtype::Node => unimplemented!(),
            });
        }

        if let Some(rtype) = rtype {
            cmd.arg(match rtype {
                Rtype::Greedy => "-rtype=greedy",
                Rtype::Fm => "-rtype=fm",
                Rtype::Sep2Sided => "-rtype=1sided",
                Rtype::Sep1Sided => "-rtype=2sided",
            });
        }

        if let Some(iptype) = iptype {
            cmd.arg(match iptype {
                Iptype::Grow => "-iptype=grow",
                Iptype::Random => "-iptype=random",
                Iptype::Edge => "-iptype=edge",
                Iptype::Node => "-iptype=node",
                Iptype::Rb => "-iptype=rb",
            });
        }

        if let Some(niparts) = niparts
        {
            cmd.arg(format!("-niparts={}", niparts));
        }
        if let Some(niter) = niter
        {
            cmd.arg(format!("-niter={}", niter));
        }
        if let Some(nseps) = nseps
        {
            cmd.arg(format!("-nseps={}", nseps));
        }
        if let Some(ncuts) = ncuts
        {
            cmd.arg(format!("-ncuts={}", ncuts));
        }
        if let Some(seed) = seed {
            cmd.arg(format!("-seed={}", seed));
        }
        if let Some(dbglvl) = dbglvl
        {
            cmd.arg(format!("-dbglvl={}", dbglvl));
        }
        if let Some(pfactor) = pfactor
        {
            cmd.arg(format!("-pfactor={}", pfactor));
        }

        cmd.stdin(Stdio::null());
        cmd.stdout(Stdio::piped());
        cmd.stderr(Stdio::piped());

        // println!("{cmd:?}");
        unset_cloexec(&pw).unwrap();
        if PRINT_CMDS {
            eprintln!("{cmd:?}")
        }
        let child = cmd.spawn().expect("could not spawn process");
        drop(pw);
        drop(tpr);

        let mut partfile = String::new();
        let mut pipe = PipedChild::new();
        if let Some(tpw) = tpw {
            pipe.writer(tpw, tp_buf.as_bytes());
        }
        let output = pipe.reader(pr, &mut partfile).wait_with_output(child).expect("child failed");
        drop(pipe);

        let parts = partfile.split_whitespace().map(|i| i.parse().unwrap()).collect();


        Ok((output, parts))
    }

    pub fn call_suite_assert_eq(&mut self) {
        for graph in TestGraph::test_suite() {
            self.filename = graph.into();
            self.call_assert_eq();
        }
    }

    pub fn call_assert_eq(&self) {
        assert!(cfg!(feature = "normalized"), "normalized is not set, output is guaranteed to differ");
        let (co, ce, cp) = self.call(Ver::C).unwrap();
        let (ro, re, rp) = self.call(Ver::Rust).unwrap();
        assert_eq_diff_str(co, ro, "stdout");
        assert_eq_diff_str(ce, re, "stderr");
        assert_eq_diff(cp, rp, "partition", true);
    }
}

#[derive(Default)]
struct PipedChild<'a> {
    // need Box because we want to drop when done
    writers: Vec<(Box<dyn Write + Send + 'a>, &'a [u8])>,
    readers: Vec<(Box<dyn Read + Send + 'a>, &'a mut String)>,
}

impl<'a> PipedChild<'a> {
    fn new() -> Self {
        Default::default()
    }

    fn writer(&mut self, w: impl Write + Send + 'a, buf: &'a [u8]) -> &mut Self {
        self.writers.push((Box::new(w), buf.as_ref()));
        self
    }

    fn reader(&mut self, w: impl Read + Send + 'a, buf: &'a mut String) -> &mut Self {
        self.readers.push((Box::new(w), buf));
        self
    }

    /// use threads to read and write to pipes, and wait for child with output. Right now this spawns a
    /// bunch of threads, but I could make it faster in the future if it actually seems to be a problem
    fn wait_with_output(&mut self, child: Child) -> io::Result<Output> {
        // don't want to block and async is a pain
        // I could change this in the future to use non-blocking 
        let Self { writers, readers } = std::mem::take(self);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (mut w, buf) in writers {
                handles.push(scope.spawn(move || w.write_all(buf)))
            }
            for (mut r, buf) in readers {
                handles.push(scope.spawn(move || r.read_to_string(buf).map(|_| ())));
            }
            let res = child.wait_with_output().unwrap();
            for handle in handles {
                handle.join().expect("thread panicked (should not be possible)")?
            }
            Ok(res)
        })
    }
}


fn unset_cloexec(s: impl AsFd) -> io::Result<()> {
    fn inner(s: BorrowedFd) -> io::Result<()> {
        unsafe {
            let flags = libc::fcntl(s.as_raw_fd(), libc::F_GETFD);
            if flags < 0 {
                return Err(io::Error::last_os_error())
            }
            if libc::fcntl(s.as_raw_fd(), libc::F_SETFD, flags & !libc::FD_CLOEXEC) == -1 {
                return Err(io::Error::last_os_error())
            }
        }
        Ok(())
    }
    inner(s.as_fd())
}

#[inline(never)]
fn assert_eq_diff_failed(c_out: &str, rs_out: &str, name: &str, two_col: bool) -> ! {

    let (cpr, cpw) = io::pipe().expect("failed to make pipe");
    let (rpr, rpw) = io::pipe().expect("failed to make pipe");
    unset_cloexec(&cpr).unwrap();
    unset_cloexec(&rpr).unwrap();

    println!("c {name} and rs {name} differ");
    let c_fd = format!("/proc/self/fd/{}", cpr.as_raw_fd());
    let rs_fd = format!("/proc/self/fd/{}", rpr.as_raw_fd());
    println!("c {name} vs rs {name} ({c_fd} vs {rs_fd})");
    unset_cloexec(&cpr).unwrap();
    unset_cloexec(&rpr).unwrap();
    let child = Command::new("diff")
        .arg(c_fd)
        .arg(rs_fd)
        .arg("--minimal")
        .arg("--color=always")
        .args(two_col.then_some("-y"))
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn diff");
    drop(cpr);
    drop(rpr);

    let res = PipedChild::new()
        .writer(cpw, c_out.as_bytes())
        .writer(rpw, rs_out.as_bytes())
        .wait_with_output(child)
        .expect("diff failed");
    print!("{}", String::from_utf8_lossy(&res.stdout));
    eprint!("{}", String::from_utf8_lossy(&res.stderr));
    match res.status.code().expect("terminated by signal") {
        0 => panic!("c and rs {name} are actually the same!"),
        1 => (),
        2 => panic!("diff(1) failed in odd way ({})", res.status),
        _ => panic!("diff(1) failed in odd (and supposedly impossible) way: ({})", res.status),
    }
    panic!("c {name} and rs {name} are different!")
}

pub fn assert_eq_diff_str(c_out: impl AsRef<str>, rs_out: impl AsRef<str>, name: impl AsRef<str>) {
    fn inner(c_out: &str, rs_out: &str, name: &str) {
        if c_out == rs_out {
            return
        }
        assert_eq_diff_failed(c_out, rs_out, name, false);
    }
    inner(c_out.as_ref(), rs_out.as_ref(), name.as_ref());
}

pub fn assert_eq_diff<T: PartialEq + Debug>(c_out: T, rs_out: T, name: impl AsRef<str>, two_col: bool) {
    #[inline(never)]
    fn inner_failed(c_out: &dyn fmt::Debug, rs_out: &dyn fmt::Debug, name: &str, two_col: bool) {
        let c_out = format!("{c_out:#?}");
        let rs_out = format!("{rs_out:#?}");
        assert_eq_diff_failed(&c_out, &rs_out, name, two_col);
    }

    if c_out != rs_out {
        inner_failed(&c_out, &rs_out, name.as_ref(), two_col);
    }
}
