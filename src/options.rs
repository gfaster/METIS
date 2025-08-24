/**
\file
\brief This file contains various routines for dealing with options and ctrl_t.

\date   Started 5/12/2011
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: options.c 17717 2014-10-03 19:09:31Z dominique $ \endverbatim
*/
use crate::*;

// this 'needs' to be a macro because some of the type coersion is somewhat annoying. This will be
// pretty easy to remove in the future, but this works for now
macro_rules! get_option {
    ($options:expr, $opt:expr, $default:expr) => {
        if $options.is_null() || *$options.add($opt as usize) == -1 {
            $default
        } else {
            *$options.add($opt as usize) as _
        }
    };
}

/*************************************************************************/
/* This function creates and sets the run parameters (ctrl_t) */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupCtrl(
    optype: moptype_et,
    options: *const idx_t,
    ncon: idx_t,
    nparts: idx_t,
    tpwgts: *const real_t,
    ubvec: *const real_t,
) -> *mut ctrl_t {
    // idx_t i, j;
    // ctrl_t *ctrl;

    let ctrl: *mut ctrl_t = gk_malloc(size_of::<ctrl_t>(), c"SetupCtrl: ctrl".as_ptr()).cast();
    *ctrl = std::mem::zeroed();
    let ctrl = &mut *ctrl;

    // memset((void *)ctrl, 0, sizeof(ctrl_t));

    ctrl.pid = std::process::id() as libc::pid_t;

    match optype {
        METIS_OP_PMETIS => {
            ctrl.objtype = get_option!(options, METIS_OPTION_OBJTYPE, METIS_OBJTYPE_CUT);
            ctrl.rtype = METIS_RTYPE_FM;
            ctrl.ncuts = get_option!(options, METIS_OPTION_NCUTS, 1);
            ctrl.niter = get_option!(options, METIS_OPTION_NITER, 10);

            if ncon == 1 {
                ctrl.iptype = get_option!(options, METIS_OPTION_IPTYPE, METIS_IPTYPE_GROW);
                ctrl.ufactor = get_option!(options, METIS_OPTION_UFACTOR, PMETIS_DEFAULT_UFACTOR);
                ctrl.CoarsenTo = 20;
            } else {
                ctrl.iptype = get_option!(options, METIS_OPTION_IPTYPE, METIS_IPTYPE_RANDOM);
                ctrl.ufactor = get_option!(options, METIS_OPTION_UFACTOR, MCPMETIS_DEFAULT_UFACTOR);
                ctrl.CoarsenTo = 100;
            }
        }

        METIS_OP_KMETIS => {
            ctrl.objtype = get_option!(options, METIS_OPTION_OBJTYPE, METIS_OBJTYPE_CUT);
            ctrl.iptype = get_option!(options, METIS_OPTION_IPTYPE, METIS_IPTYPE_METISRB);
            ctrl.rtype = METIS_RTYPE_GREEDY;
            ctrl.nIparts = get_option!(options, METIS_OPTION_NIPARTS, -1);
            ctrl.ncuts = get_option!(options, METIS_OPTION_NCUTS, 1);
            ctrl.niter = get_option!(options, METIS_OPTION_NITER, 10);
            ctrl.ufactor = get_option!(options, METIS_OPTION_UFACTOR, KMETIS_DEFAULT_UFACTOR);
            ctrl.minconn = get_option!(options, METIS_OPTION_MINCONN, 0);
            ctrl.contig = get_option!(options, METIS_OPTION_CONTIG, 0);
        }

        METIS_OP_OMETIS => {
            ctrl.objtype = get_option!(options, METIS_OPTION_OBJTYPE, METIS_OBJTYPE_NODE);
            ctrl.rtype = get_option!(options, METIS_OPTION_RTYPE, METIS_RTYPE_SEP1SIDED);
            ctrl.iptype = get_option!(options, METIS_OPTION_IPTYPE, METIS_IPTYPE_EDGE);
            ctrl.nseps = get_option!(options, METIS_OPTION_NSEPS, 1);
            ctrl.niter = get_option!(options, METIS_OPTION_NITER, 10);
            ctrl.ufactor = get_option!(options, METIS_OPTION_UFACTOR, OMETIS_DEFAULT_UFACTOR);
            ctrl.compress = get_option!(options, METIS_OPTION_COMPRESS, 1);
            ctrl.ccorder = get_option!(options, METIS_OPTION_CCORDER, 0);
            ctrl.pfactor = 0.1 * get_option!(options, METIS_OPTION_PFACTOR, 0) as real_t;

            ctrl.CoarsenTo = 100;
        }

        _ => panic!("Unknown optype of {}", optype),
    }

    /* common options */
    ctrl.ctype = get_option!(options, METIS_OPTION_CTYPE, METIS_CTYPE_SHEM);
    ctrl.no2hop = get_option!(options, METIS_OPTION_NO2HOP, 0);
    ctrl.ondisk = get_option!(options, METIS_OPTION_ONDISK, 0);
    ctrl.seed = get_option!(options, METIS_OPTION_SEED, -1);
    ctrl.dbglvl = get_option!(options, METIS_OPTION_DBGLVL, 0);
    ctrl.numflag = get_option!(options, METIS_OPTION_NUMBERING, 0);
    ctrl.dropedges = get_option!(options, METIS_OPTION_DROPEDGES, 0);

    /* set non-option information */
    ctrl.optype = optype;
    ctrl.ncon = ncon;
    ctrl.nparts = nparts;
    ctrl.maxvwgt = gk::ismalloc(ncon as usize, 0, c"SetupCtrl: maxvwgt").as_buf_ptr();

    /* setup the target partition weights */
    if ctrl.optype != METIS_OP_OMETIS {
        ctrl.tpwgts =
            gk::rsmalloc((nparts * ncon) as usize, 0.0, c"SetupCtrl: ctrl.tpwgts").as_buf_ptr();
        if !tpwgts.is_null() {
            std::ptr::copy_nonoverlapping(tpwgts, ctrl.tpwgts, (nparts * ncon) as usize);
        } else {
            mkslice_mut!(ctrl->tpwgts, nparts * ncon);
            for i in (0)..(nparts) {
                for j in (0)..(ncon) {
                    tpwgts[(i * ncon + j) as usize] = 1.0 / nparts as real_t;
                }
            }
        }
    } else {
        /* METIS_OP_OMETIS */
        /* this is required to allow the pijbm to be defined properly for
        the edge-based refinement during initial partitioning */
        ctrl.tpwgts = gk::rsmalloc(2, 0.5, c"SetupCtrl: ctrl.tpwgts").as_buf_ptr();
    }

    /* setup the ubfactors */
    ctrl.ubfactors = gk::rsmalloc(
        ctrl.ncon as usize,
        util::i2rubfactor(ctrl.ufactor) as real_t,
        c"SetupCtrl: ubfactors",
    )
    .as_buf_ptr();
    if !ubvec.is_null() {
        std::ptr::copy_nonoverlapping(ubvec, ctrl.ubfactors, ctrl.ncon as usize);
    }
    for i in (0)..(ctrl.ncon) {
        *ctrl.ubfactors.add(i as usize) += 0.0000499;
    }

    /* Allocate memory for balance multipliers.
    Note that for PMETIS/OMETIS routines the memory allocated is more
    than required as balance multipliers for 2 parts is sufficient. */
    ctrl.pijbm = gk::rmalloc((nparts * ncon) as usize, c"SetupCtrl: ctrl.pijbm").as_buf_ptr();

    util::InitRandom(ctrl.seed);

    ifset!(ctrl.dbglvl, METIS_DBG_INFO, PrintCtrl(ctrl));

    if CheckParams(ctrl) == 0 {
        let mut ctrl: *mut _ = ctrl;
        FreeCtrl(&mut ctrl);
        return std::ptr::null_mut();
    } else {
        return ctrl;
    }
}

/*************************************************************************/
/* Computes the per-partition/constraint balance multipliers */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupKWayBalMultipliers(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    let nparts = ctrl.nparts as usize;
    let ncon = graph.ncon as usize;
    mkslice!(ctrl->tpwgts, nparts * ncon);
    mkslice_mut!(ctrl->pijbm, nparts * ncon);
    get_graph_slices!(graph => invtvwgt);

    for i in (0)..(nparts) {
        for j in (0)..(ncon) {
            pijbm[i * ncon + j] = invtvwgt[j as usize] / tpwgts[i * ncon + j];
        }
    }
}

/*************************************************************************/
/* Computes the per-partition/constraint balance multipliers */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Setup2WayBalMultipliers(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    tpwgts: *mut real_t,
) {
    let ctrl = ctrl.as_ref().unwrap();
    let graph = graph.as_mut().unwrap();
    let nparts = ctrl.nparts as usize;
    let ncon = graph.ncon as usize;
    mkslice!(tpwgts, nparts * ncon);
    mkslice_mut!(ctrl->pijbm, nparts * ncon);
    get_graph_slices!(graph => invtvwgt);

    for i in (0)..(2) {
        for j in (0)..(ncon) {
            pijbm[i * ncon + j] = invtvwgt[j] / tpwgts[i * ncon + j];
        }
    }
}

/*************************************************************************/
/* This function prints the various control fields */
/*************************************************************************/
#[metis_func]
pub extern "C" fn PrintCtrl(ctrl: *const ctrl_t) {
    // idx_t i, j, modnum;
    let ctrl = ctrl.as_ref().unwrap();
    println!(" Runtime parameters:");

    print!("   Objective type: ");
    match ctrl.objtype {
        METIS_OBJTYPE_CUT => println!("METIS_OBJTYPE_CUT"),

        METIS_OBJTYPE_VOL => println!("METIS_OBJTYPE_VOL"),

        METIS_OBJTYPE_NODE => println!("METIS_OBJTYPE_NODE"),

        _ => println!("Unknown!"),
    }

    print!("   Coarsening type: ");
    match ctrl.ctype {
        METIS_CTYPE_RM => println!("METIS_CTYPE_RM"),

        METIS_CTYPE_SHEM => println!("METIS_CTYPE_SHEM"),

        _ => println!("Unknown!"),
    }

    print!("   Initial partitioning type: ");
    match ctrl.iptype {
        METIS_IPTYPE_GROW => println!("METIS_IPTYPE_GROW"),

        METIS_IPTYPE_RANDOM => println!("METIS_IPTYPE_RANDOM"),

        METIS_IPTYPE_EDGE => println!("METIS_IPTYPE_EDGE"),

        METIS_IPTYPE_NODE => println!("METIS_IPTYPE_NODE"),

        METIS_IPTYPE_METISRB => println!("METIS_IPTYPE_METISRB"),

        _ => println!("Unknown!"),
    }

    print!("   Refinement type: ");
    match ctrl.rtype {
        METIS_RTYPE_FM => println!("METIS_RTYPE_FM"),

        METIS_RTYPE_GREEDY => println!("METIS_RTYPE_GREEDY"),

        METIS_RTYPE_SEP2SIDED => println!("METIS_RTYPE_SEP2SIDED"),

        METIS_RTYPE_SEP1SIDED => println!("METIS_RTYPE_SEP1SIDED"),

        _ => println!("Unknown!"),
    }

    println!(
        "   Perform a 2-hop matching: {}",
        if ctrl.no2hop != 0 { "No" } else { "Yes" }
    );

    println!(
        "   On disk storage: {}",
        if ctrl.ondisk != 0 { "Yes" } else { "No" }
    );
    println!(
        "   Drop edges: {}",
        if ctrl.dropedges != 0 { "Yes" } else { "No" }
    );

    println!("   Number of balancing constraints: {:}", ctrl.ncon);
    println!("   Number of refinement iterations: {:}", ctrl.niter);
    println!("   Number of initial partitionings: {:}", ctrl.nIparts);
    println!("   Random number seed: {:}", ctrl.seed);

    if ctrl.optype == METIS_OP_OMETIS {
        println!("   Number of separators: {:}", ctrl.nseps);
        println!(
            "   Compress graph prior to ordering: {}",
            if ctrl.compress != 0 { "Yes" } else { "No" }
        );
        println!(
            "   Detect & order connected components separately: {}",
            if ctrl.ccorder != 0 { "Yes" } else { "No" }
        );
        println!(
            "   Prunning factor for high degree vertices: {:}",
            ctrl.pfactor
        );
    } else {
        println!("   Number of partitions: {:}", ctrl.nparts);
        println!("   Number of cuts: {:}", ctrl.ncuts);
        println!("   User-supplied ufactor: {:}", ctrl.ufactor);

        if ctrl.optype == METIS_OP_KMETIS {
            println!(
                "   Minimize connectivity: {}",
                if ctrl.minconn != 0 { "Yes" } else { "No" }
            );
            println!(
                "   Create contiguous partitions: {}",
                if ctrl.contig != 0 { "Yes" } else { "No" }
            );
        }
        let modnum = match ctrl.ncon {
            1 => 5,
            2 => 3,
            3 => 2,
            _ => 1,
        };
        print!("   Target partition weights: ");
        for i in (0)..(ctrl.nparts) {
            if i % modnum == 0 {
                print!("\n     ");
            }
            print!("{:4}=[", i);
            for j in (0)..(ctrl.ncon) {
                mkslice!(ctrl->tpwgts, ctrl.ncon * ctrl.nparts);
                // C and Rust scientific notation formatting is slightly different
                if cfg!(feature = "normalized") {
                    print!(
                        "{}{:.5}",
                        (if j == 0 { "" } else { " " }),
                        tpwgts[(i * ctrl.ncon + j) as usize]
                    );
                } else {
                    print!(
                        "{}{:.2e}",
                        (if j == 0 { "" } else { " " }),
                        tpwgts[(i * ctrl.ncon + j) as usize]
                    );
                }
            }
            print!("]");
        }
        println!();
    }

    print!("   Allowed maximum load imbalance: ");
    mkslice!(ctrl->ubfactors, ctrl.ncon);
    for i in (0)..(ctrl.ncon) {
        print!("{:.3} ", ubfactors[i as usize]);
    }
    println!();

    println!();
}

/*************************************************************************/
/* This function checks the validity of user-supplied parameters */
/*************************************************************************/
#[metis_func]
pub extern "C" fn CheckParams(ctrl: *const ctrl_t) -> libc::c_int {
    // idx_t i, j;
    // real_t sum;
    let dbglvl = METIS_DBG_INFO;
    let ctrl = ctrl.as_ref().unwrap();

    match ctrl.optype {
        METIS_OP_PMETIS => {
            if ctrl.objtype != METIS_OBJTYPE_CUT {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect objective type.\n")
                );
                return 0;
            }
            if ctrl.ctype != METIS_CTYPE_RM && ctrl.ctype != METIS_CTYPE_SHEM {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect coarsening scheme.\n")
                );
                return 0;
            }
            if ctrl.iptype != METIS_IPTYPE_GROW && ctrl.iptype != METIS_IPTYPE_RANDOM {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect initial partitioning scheme.\n")
                );
                return 0;
            }
            if ctrl.rtype != METIS_RTYPE_FM {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect refinement scheme.\n")
                );
                return 0;
            }
            if ctrl.ncuts <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ncuts.\n")
                );
                return 0;
            }
            if ctrl.niter <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect niter.\n")
                );
                return 0;
            }
            if ctrl.ufactor <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ufactor.\n")
                );
                return 0;
            }
            if ctrl.numflag != 0 && ctrl.numflag != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect numflag.\n")
                );
                return 0;
            }
            if ctrl.nparts <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect nparts.\n")
                );
                return 0;
            }
            if ctrl.ncon <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ncon.\n")
                );
                return 0;
            }

            mkslice!(ctrl->tpwgts, ctrl.nparts * ctrl.ncon);
            for i in (0)..(ctrl.ncon) {
                let sum = tpwgts[i as usize..]
                    .iter()
                    .step_by(ctrl.ncon as usize)
                    .sum::<real_t>();
                // let sum = rsum(ctrl.nparts, ctrl.tpwgts+i, ctrl.ncon);
                if !(0.99..=1.01).contains(&sum) {
                    ifset!(
                        dbglvl,
                        METIS_DBG_INFO,
                        print!(
                            "Input Error: Incorrect sum of {:} for tpwgts for constraint {:}.\n",
                            sum, i
                        )
                    );
                    return 0;
                }
            }
            for i in (0)..(ctrl.ncon) {
                for j in (0)..(ctrl.nparts) {
                    if tpwgts[(j * ctrl.ncon + i) as usize] <= 0.0 {
                        ifset!(
                            dbglvl,
                            METIS_DBG_INFO,
                            print!(
                                "Input Error: Incorrect tpwgts for partition {:} and constraint {:}.\n",
                                j, i
                            )
                        );
                        return 0;
                    }
                }
            }

            mkslice!(ctrl->ubfactors, ctrl.ncon);
            for i in (0)..(ctrl.ncon) {
                if ubfactors[i as usize] <= 1.0 {
                    ifset!(
                        dbglvl,
                        METIS_DBG_INFO,
                        print!("Input Error: Incorrect ubfactor for constraint {:}.\n", i)
                    );
                    return 0;
                }
            }
        }

        METIS_OP_KMETIS => {
            if ctrl.objtype != METIS_OBJTYPE_CUT && ctrl.objtype != METIS_OBJTYPE_VOL {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect objective type.\n")
                );
                return 0;
            }
            if ctrl.ctype != METIS_CTYPE_RM && ctrl.ctype != METIS_CTYPE_SHEM {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect coarsening scheme.\n")
                );
                return 0;
            }
            if ctrl.iptype != METIS_IPTYPE_METISRB && ctrl.iptype != METIS_IPTYPE_GROW {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect initial partitioning scheme.\n")
                );
                return 0;
            }
            if ctrl.rtype != METIS_RTYPE_GREEDY {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect refinement scheme.\n")
                );
                return 0;
            }
            if ctrl.ncuts <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ncuts.\n")
                );
                return 0;
            }
            if ctrl.niter <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect niter.\n")
                );
                return 0;
            }
            if ctrl.ufactor <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ufactor.\n")
                );
                return 0;
            }
            if ctrl.numflag != 0 && ctrl.numflag != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect numflag.\n")
                );
                return 0;
            }
            if ctrl.nparts <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect nparts.\n")
                );
                return 0;
            }
            if ctrl.ncon <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ncon.\n")
                );
                return 0;
            }
            if ctrl.contig != 0 && ctrl.contig != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect contig.\n")
                );
                return 0;
            }
            if ctrl.minconn != 0 && ctrl.minconn != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect minconn.\n")
                );
                return 0;
            }

            mkslice!(ctrl->tpwgts, ctrl.nparts * ctrl.ncon);
            for i in (0)..(ctrl.ncon) {
                // sum = rsum(ctrl.nparts, ctrl.tpwgts+i, ctrl.ncon);
                let sum = tpwgts[i as usize..]
                    .iter()
                    .step_by(ctrl.ncon as usize)
                    .sum::<real_t>();
                if !(0.99..=1.01).contains(&sum) {
                    ifset!(
                        dbglvl,
                        METIS_DBG_INFO,
                        print!(
                            "Input Error: Incorrect sum of {:} for tpwgts for constraint {:}.\n",
                            sum, i
                        )
                    );
                    return 0;
                }
            }
            for i in (0)..(ctrl.ncon) {
                for j in (0)..(ctrl.nparts) {
                    if tpwgts[(j * ctrl.ncon + i) as usize] <= 0.0 {
                        ifset!(
                            dbglvl,
                            METIS_DBG_INFO,
                            print!(
                                "Input Error: Incorrect tpwgts for partition {:} and constraint {:}.\n",
                                j, i
                            )
                        );
                        return 0;
                    }
                }
            }

            mkslice!(ctrl->ubfactors, ctrl.ncon);
            for i in (0)..(ctrl.ncon) {
                if ubfactors[i as usize] <= 1.0 {
                    ifset!(
                        dbglvl,
                        METIS_DBG_INFO,
                        print!("Input Error: Incorrect ubfactor for constraint {:}.\n", i)
                    );
                    return 0;
                }
            }
        }

        METIS_OP_OMETIS => {
            if ctrl.objtype != METIS_OBJTYPE_NODE {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect objective type.\n")
                );
                return 0;
            }
            if ctrl.ctype != METIS_CTYPE_RM && ctrl.ctype != METIS_CTYPE_SHEM {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect coarsening scheme.\n")
                );
                return 0;
            }
            if ctrl.iptype != METIS_IPTYPE_EDGE && ctrl.iptype != METIS_IPTYPE_NODE {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect initial partitioning scheme.\n")
                );
                return 0;
            }
            if ctrl.rtype != METIS_RTYPE_SEP1SIDED && ctrl.rtype != METIS_RTYPE_SEP2SIDED {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect refinement scheme.\n")
                );
                return 0;
            }
            if ctrl.nseps <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect nseps.\n")
                );
                return 0;
            }
            if ctrl.niter <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect niter.\n")
                );
                return 0;
            }
            if ctrl.ufactor <= 0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ufactor.\n")
                );
                return 0;
            }
            if ctrl.numflag != 0 && ctrl.numflag != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect numflag.\n")
                );
                return 0;
            }
            if ctrl.nparts != 3 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect nparts.\n")
                );
                return 0;
            }
            if ctrl.ncon != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ncon.\n")
                );
                return 0;
            }
            if ctrl.compress != 0 && ctrl.compress != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect compress.\n")
                );
                return 0;
            }
            if ctrl.ccorder != 0 && ctrl.ccorder != 1 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect ccorder.\n")
                );
                return 0;
            }
            if ctrl.pfactor < 0.0 {
                ifset!(
                    dbglvl,
                    METIS_DBG_INFO,
                    print!("Input Error: Incorrect pfactor.\n")
                );
                return 0;
            }

            mkslice!(ctrl->ubfactors, ctrl.ncon);
            for i in (0)..(ctrl.ncon) {
                if ubfactors[i as usize] <= 1.0 {
                    ifset!(
                        dbglvl,
                        METIS_DBG_INFO,
                        print!("Input Error: Incorrect ubfactor for constraint {:}.\n", i)
                    );
                    return 0;
                }
            }
        }

        _ => {
            ifset!(
                dbglvl,
                METIS_DBG_INFO,
                print!("Input Error: Incorrect optype\n")
            );
            return 0;
        }
    }

    return 1;
}

/*************************************************************************/
/* This function frees the memory associated with a ctrl_t */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeCtrl(ctrl: *mut *mut ctrl_t) {
    // ctrl_t *ctrl = *r_ctrl;

    wspace::FreeWorkSpace(*ctrl);
    gk::free_ref(&mut (**ctrl).tpwgts);
    gk::free_ref(&mut (**ctrl).pijbm);
    gk::free_ref(&mut (**ctrl).ubfactors);
    gk::free_ref(&mut (**ctrl).maxvwgt);
    gk::free_ref(&mut (*ctrl));

    // gk_free((void **)&ctrl.tpwgts, &ctrl.pijbm,
    //         &ctrl.ubfactors, &ctrl.maxvwgt, &ctrl, LTERM);
    // *r_ctrl = std::ptr::null_mut();
}
