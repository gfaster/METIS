use metis::{bindings::*, idx_t};
use std::process::ExitCode;

use crate::getopt;

// Alos have some misc utils thrown in here because why not

/// Yes, I really don't need this, but I'm keeping it for now to keep the printout semantically the
/// same
pub const SVNINFO: &str = "unknown";

pub const __DATE__: &str = "Jan  1 1980";
pub const __TIME__: &str = "00:00:00";

/// The text labels for PTypes
pub static PTYPE_NAMES: [&str; 2] = ["rb", "kway"];

/// The text labels for ObjTypes
pub static OBJTYPE_NAMES: [&str; 3] = ["cut", "vol", "node"];

/// The text labels for CTypes
pub static CTYPE_NAMES: [&str; 2] = ["rm", "shem"];

/// The text labels for RTypes
pub static RTYPE_NAMES: [&str; 4] = ["fm", "greedy", "2sided", "1sided"];

/// The text labels for ITypes
pub static IPTYPE_NAMES: [&str; 5] = ["grow", "random", "edge", "node", "metisrb"];

/// The text labels for GTypes
#[expect(dead_code)]
pub static GTYPE_NAMES: [&str; 2] = ["dual", "nodal"];

pub const METISTITLE: &str = if !cfg!(feature = "normalized") {
    "METIS 5.2.1 Copyright 1998-22, Regents of the University of Minnesota. Copyright 2023-25, Gavin Rohrer"
} else {
    "METIS NORMALIZED OUTPUT"
};

/// converts a user provided ufactor into a real ubfactor
pub fn i2rubfactor(ufactor: idx_t) -> real_t {
    1.0 + 0.001 * ufactor as real_t
}

pub fn vec_ptr<T>(v: &mut Option<Vec<T>>) -> *mut T {
    match v.as_mut() {
        Some(v) => v.as_mut_ptr(),
        None => std::ptr::null_mut(),
    }
}

#[expect(dead_code)]
pub fn vec_ptr2<T>(v: &mut Vec<T>) -> *mut T {
    if v.is_empty() {
        std::ptr::null_mut()
    } else {
        v.as_mut_ptr()
    }
}

pub use metis::params::params_t;

fn parse_inner<'a, 'o>(
    args: &'a [String],
    opts: &'o [getopt::Opt],
    p: &mut params_t,
) -> Result<(), Box<getopt::ParsedOptionError<'a, 'o>>> {
    let mut graphfile = None;
    let mut nparts: Option<idx_t> = None;
    for arg in getopt::parse_options(args, opts) {
        let arg = arg?;
        let (val, opt) = match arg {
            getopt::ParsedOption::Opt { val, opt } => (val, opt),
            getopt::ParsedOption::Nonopt(val) if graphfile.is_none() => {
                graphfile = Some(val);
                continue;
            }
            getopt::ParsedOption::Nonopt(val) if nparts.is_none() => {
                nparts = Some(val.parse()?);
                continue;
            }
            getopt::ParsedOption::Nonopt(_) => {
                return Err(getopt::ParsedOptionError::TooManyArgs.into());
            }
        };
        let val = val.unwrap_or("<NO ARGUMENT>");
        match opt.idx as moptions_et {
            METIS_OPTION_PTYPE => p.ptype = opt.parse_enum(val)?,
            METIS_OPTION_OBJTYPE => p.objtype = opt.parse_enum(val)?,
            METIS_OPTION_CTYPE => p.ctype = opt.parse_enum(val)?,
            METIS_OPTION_IPTYPE => p.iptype = opt.parse_enum(val)?,
            METIS_OPTION_RTYPE => p.rtype = opt.parse_enum(val)?,
            METIS_OPTION_DBGLVL => p.dbglvl = val.parse()?,
            METIS_OPTION_NIPARTS => p.niparts = val.parse()?,
            METIS_OPTION_NITER => p.niter = val.parse()?,
            METIS_OPTION_NCUTS => p.ncuts = val.parse()?,
            METIS_OPTION_SEED => p.seed = val.parse()?,
            METIS_OPTION_ONDISK => p.ondisk = true,
            METIS_OPTION_MINCONN => p.minconn = true,
            METIS_OPTION_CONTIG => p.contig = true,
            METIS_OPTION_COMPRESS => p.compress = 1,
            METIS_OPTION_CCORDER => p.ccorder = 1,
            METIS_OPTION_PFACTOR => p.pfactor = val.parse()?,
            METIS_OPTION_NSEPS => p.nseps = val.parse()?,
            METIS_OPTION_UFACTOR => p.ufactor = val.parse()?,
            METIS_OPTION_NUMBERING => unimplemented!(),
            METIS_OPTION_DROPEDGES => p.dropedges = 1,
            METIS_OPTION_NO2HOP => p.no2hop = true,
            METIS_OPTION_TWOHOP => p.no2hop = false,
            METIS_OPTION_FAST => unimplemented!(),

            METIS_OPTION_HELP => return Err(getopt::ParsedOptionError::Help.into()),
            METIS_OPTION_TPWGTS => p.tpwgtsfile = Some(val.into()),
            METIS_OPTION_NCOMMON => p.ncommon = val.parse::<std::num::NonZeroU16>()?.get() as idx_t,
            METIS_OPTION_NOOUTPUT => p.nooutput = true,
            METIS_OPTION_BALANCE => unimplemented!(),
            METIS_OPTION_GTYPE => p.gtype = opt.parse_enum(val)?,
            METIS_OPTION_UBVEC => p.ubvecstr = Some(val.into()),
            METIS_OPTION_OUTFILE => p.outfile = Some(val.into()),
            _ => panic!("unexpected option: {opt:?}"),
        }
    }

    let (Some(graphfile), Some(nparts)) = (graphfile, nparts) else {
        return Err(getopt::ParsedOptionError::NotEnoughArgs.into());
    };
    p.filename = graphfile.into();
    p.nparts = nparts.try_into().unwrap();
    Ok(())
}

pub fn parse_standard(
    args: &[String],
    help: &[&str],
    shorthelp: &[&str],
    opts: &[getopt::Opt],
    params: &mut params_t,
) -> Option<ExitCode> {
    match parse_inner(args, opts, params) {
        Ok(()) => return None,
        Err(e) if e.is_help() => {
            for line in help {
                println!("{line}");
            }
            return Some(ExitCode::SUCCESS);
        }
        Err(e) => {
            eprintln!("ERROR: {e:?}");
            for line in shorthelp {
                eprintln!("{line}");
            }
            return Some(ExitCode::FAILURE);
        }
    }
}
