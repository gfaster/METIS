use std::path::PathBuf;

use crate::*;

#[allow(non_camel_case_types)]
#[derive(Debug)]
pub struct params_t {
    pub ptype: idx_t,
    pub objtype: idx_t,
    pub ctype: idx_t,
    pub iptype: idx_t,
    pub rtype: idx_t,

    pub no2hop: bool,
    pub minconn: bool,
    pub contig: bool,

    pub ondisk: bool,

    pub dropedges: idx_t,

    pub nooutput: bool,

    pub balance: idx_t,
    pub ncuts: idx_t,
    pub niter: idx_t,
    pub niparts: idx_t,

    pub gtype: idx_t,
    pub ncommon: idx_t,

    pub seed: idx_t,
    pub dbglvl: idx_t,

    pub nparts: usize,

    pub nseps: idx_t,
    pub ufactor: idx_t,
    pub pfactor: idx_t,
    pub compress: idx_t,
    pub ccorder: idx_t,

    pub filename: PathBuf,
    pub outfile: Option<PathBuf>,
    pub xyzfile: Option<PathBuf>,
    pub tpwgtsfile: Option<PathBuf>,
    pub ubvecstr: Option<String>,

    pub wgtflag: idx_t,
    pub numflag: idx_t,
    pub tpwgts: Vec<real_t>,
    pub ubvec: Option<Vec<real_t>>,

    pub iotimer: f64,
    pub parttimer: f64,
    pub reporttimer: f64,

    pub maxmemory: usize,
}

impl params_t {
    pub fn new() -> Self {
        Self {
            ptype: -1,
            objtype: -1,
            ctype: -1,
            iptype: -1,
            rtype: -1,
            no2hop: false,
            minconn: false,
            contig: false,
            ondisk: false,
            dropedges: -1,
            nooutput: false,
            balance: -1,
            ncuts: -1,
            niter: -1,
            niparts: -1,
            gtype: -1,
            ncommon: -1,
            seed: -1,
            dbglvl: 0,
            nparts: usize::MAX,
            nseps: -1,
            ufactor: -1,
            pfactor: -1,
            compress: -1,
            ccorder: -1,
            filename: PathBuf::new(),
            outfile: None,
            xyzfile: None,
            tpwgtsfile: None,
            ubvecstr: None,
            wgtflag: 3,
            numflag: 0,
            tpwgts: Vec::new(),
            ubvec: None,
            iotimer: 0.0,
            parttimer: 0.0,
            reporttimer: 0.0,
            maxmemory: 0,
        }
    }
}
