use std::path::PathBuf;

use crate::*;

#[allow(non_camel_case_types)]
pub struct params_t {
    pub ptype: idx_t,
    pub objtype: idx_t,
    pub ctype: idx_t,
    pub iptype: idx_t,
    pub rtype: idx_t,

    pub no2hop: idx_t,
    pub minconn: idx_t,
    pub contig: idx_t,

    pub ondisk: idx_t,

    pub dropedges: idx_t,

    pub nooutput: idx_t,

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
    pub outfile: PathBuf,
    pub xyzfile: PathBuf,
    pub tpwgtsfile: Option<PathBuf>,
    /// this was a char*, but I haven't found its use
    pub ubvecstr: String,

    pub wgtflag: idx_t,
    pub numflag: idx_t,
    pub tpwgts: *mut real_t,
    pub ubvec: *mut real_t,

    pub iotimer: real_t,
    pub parttimer: real_t,
    pub reporttimer: real_t,

    pub maxmemory: usize,
}
