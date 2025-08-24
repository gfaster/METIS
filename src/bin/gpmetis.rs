mod getopt;
mod params;
use std::process::ExitCode;

use metis::bindings::*;
use metis::timing::{gk_getcputimer, gk_startcputimer, gk_stopcputimer};
use metis::{KMETIS_DEFAULT_UFACTOR, MCPMETIS_DEFAULT_UFACTOR, PMETIS_DEFAULT_UFACTOR};

use crate::params::*;

#[rustfmt::skip]
const OPTIONS: &[getopt::Opt] = getopt::opt! {
    ("ptype",          Required,      METIS_OPTION_PTYPE, {"rb" = METIS_PTYPE_RB, "kway" = METIS_PTYPE_KWAY}),

    ("objtype",        Required,      METIS_OPTION_OBJTYPE, {"cut" = METIS_OBJTYPE_CUT, "vol" = METIS_OBJTYPE_VOL}), 

    ("ctype",          Required,      METIS_OPTION_CTYPE, {"rm" = METIS_CTYPE_RM, "shem" = METIS_CTYPE_SHEM}),

    ("iptype",         Required,      METIS_OPTION_IPTYPE, {"grow" = METIS_IPTYPE_GROW, "random" = METIS_IPTYPE_RANDOM, "rb" = METIS_IPTYPE_METISRB}),

    /* ("rtype",          Required,      METIS_OPTION_RTYPE, {"fm" = METIS_RTYPE_FM, "random" = METIS_RTYPE_RANDOM, "greedy" = METIS_RTYPE_GREEDY}), */
    /* ("balanced",       NoArg,         METIS_OPTION_BALANCE) , */

    ("no2hop",         NoArg,         METIS_OPTION_NO2HOP),
    ("minconn",        NoArg,         METIS_OPTION_MINCONN),
    ("contig",         NoArg,         METIS_OPTION_CONTIG),

    ("ondisk",         NoArg,         METIS_OPTION_ONDISK),

    ("dropedges",      NoArg,         METIS_OPTION_DROPEDGES),

    ("nooutput",       NoArg,         METIS_OPTION_NOOUTPUT),

    ("ufactor",        Required,      METIS_OPTION_UFACTOR),
    ("niter",          Required,      METIS_OPTION_NITER),
    ("ncuts",          Required,      METIS_OPTION_NCUTS),
    ("niparts",        Required,      METIS_OPTION_NIPARTS),

    ("tpwgts",         Required,      METIS_OPTION_TPWGTS),
    ("ubvec",          Required,      METIS_OPTION_UBVEC),

    ("seed",           Required,      METIS_OPTION_SEED),

    ("dbglvl",         Required,      METIS_OPTION_DBGLVL),

    ("outfile",        Required,      METIS_OPTION_OUTFILE),

    ("help",           NoArg,         METIS_OPTION_HELP),
};

// for ndmetis
// #[rustfmt::skip]
// const OPTIONS: &[getopt::Opt] = getopt::opt!{
//     ("ctype",          Required,      METIS_OPTION_CTYPE),
//     ("iptype",         Required,      METIS_OPTION_IPTYPE),
//     ("rtype",          Required,      METIS_OPTION_RTYPE),
//     ("ufactor",        Required,      METIS_OPTION_UFACTOR),
//     ("pfactor",        Required,      METIS_OPTION_PFACTOR),
//     ("nocompress",     NoArg,         METIS_OPTION_COMPRESS),
//     ("ccorder",        NoArg,         METIS_OPTION_CCORDER),
//     ("no2hop",         NoArg,         METIS_OPTION_NO2HOP),
//     ("ondisk",         NoArg,         METIS_OPTION_ONDISK),
//     ("nooutput",       NoArg,         METIS_OPTION_NOOUTPUT),
//     ("niter",          Required,      METIS_OPTION_NITER),
//     ("nseps",          Required,      METIS_OPTION_NSEPS),
//     ("seed",           Required,      METIS_OPTION_SEED),
//     ("dbglvl",         Required,      METIS_OPTION_DBGLVL),
//     ("help",           NoArg,         METIS_OPTION_HELP),
// };

#[rustfmt::skip]
const HELPSTR: &[&str] =
&[
" ",
"Usage: gpmetis [options] graphfile nparts",
" ",
" Required parameters",
"    graphfile   Stores the graph to be partitioned.",
"    nparts      The number of partitions to split the graph.",
" ",
" Optional parameters",
"  -ptype=string",
"     Specifies the scheme to be used for computing the k-way partitioning.",
"     The possible values are:",
"        rb       - Recursive bisectioning",
"        kway     - Direct k-way partitioning [default]",
" ",
"  -ctype=string",
"     Specifies the scheme to be used to match the vertices of the graph",
"     during the coarsening.",
"     The possible values are:",
"        rm       - Random matching",
"        shem     - Sorted heavy-edge matching [default]",
" ",
"  -iptype=string [applies only when -ptype=rb]",
"     Specifies the scheme to be used to compute the initial partitioning",
"     of the graph.",
"     The possible values are:",
"        grow     - Grow a bisection using a greedy scheme [default for ncon=1]",
"        random   - Compute a bisection at random [default for ncon>1]",
" ",
"  -objtype=string [applies only when -ptype=kway]",
"     Specifies the objective that the partitioning routines will optimize.",
"     The possible values are:",
"        cut      - Minimize the edgecut [default]",
"        vol      - Minimize the total communication volume",
" ",
/*
"  -rtype=string",
"     Specifies the scheme to be used for refinement.",
"     The possible values are:",
"        fm       - 2-way FM refinement [default for -ptype=rb]",
"        random   - Random k-way refinement",
"        greedy   - Greedy k-way refinement [default for -ptype=kway]",
" ",
*/
"  -no2hop",
"     Specifies that the coarsening will not perform any 2-hop matchings",
"     when the standard matching fails to sufficiently contract the graph.",
" ",
"  -contig [applies only when -ptype=kway]",
"     Specifies that the partitioning routines should try to produce",
"     partitions that are contiguous. Note that if the input graph is not",
"     connected this option is ignored.",
" ",
"  -minconn [applies only when -ptype=kway]",
"     Specifies that the partitioning routines should try to minimize the",
"     maximum degree of the subdomain graph, i.e., the graph in which each",
"     partition is a node, and edges connect subdomains with a shared",
"     interface.",
" ",
"  -tpwgts=filename",
"     Specifies the name of the file that stores the target weights for",
"     each partition. By default, all partitions are assumed to be of ",
"     the same size.",
" ",
"  -ufactor=int",
"     Specifies the maximum allowed load imbalance among the partitions.",
"     A value of x indicates that the allowed load imbalance is 1+x/1000.",
"     For ptype=rb, the load imbalance is measured as the ratio of the ",
"     2*max(left,right)/(left+right), where left and right are the sizes",
"     of the respective partitions at each bisection. ",
"     For ptype=kway, the load imbalance is measured as the ratio of ",
"     max_i(pwgts[i])/avgpwgt, where pwgts[i] is the weight of the ith",
"     partition and avgpwgt is the sum of the total vertex weights divided",
"     by the number of partitions requested.",
"     For ptype=rb, the default value is 1 (i.e., load imbalance of 1.001).",
"     For ptype=kway, the default value is 30 (i.e., load imbalance of 1.03).",
" ",
"  -ubvec=string",
"     Applies only for multi-constraint partitioning and specifies the per",
"     constraint allowed load imbalance among partitions. The required ",
"     parameter corresponds to a space separated set of floating point",
"     numbers, one for each of the constraints. For example, for three",
"     constraints, the string can be \"1.02 1.2 1.35\" indicating a ",
"     desired maximum load imbalance of 2%, 20%, and 35%, respectively.",
"     The load imbalance is defined in a way similar to ufactor.",
"     If supplied, this parameter takes priority over ufactor.",
" ",
"  -niparts=int",
"     Specifies the number of initial partitions to compute. The default",
"     value is determined by the algorithm automatically.",
" ",
"  -niter=int",
"     Specifies the number of iterations for the refinement algorithms",
"     at each stage of the uncoarsening process. Default is 10.",
" ",
"  -ncuts=int",
"     Specifies the number of different partitionings that it will compute.",
"     The final partitioning is the one that achieves the best edgecut or",
"     communication volume. Default is 1.",
" ",
"  -ondisk",
"     Specifies if graphs will be stored to disk during coarsening in order",
"     to reduce memory requirements.",
" ",
"  -nooutput",
"     Specifies that no partitioning file should be generated.",
" ",
/*
"  -balance",
"     Specifies that the final partitioning should contain nparts-1 equal",
"     size partitions with the last partition having upto nparts-1 fewer",
"     vertices.",
" ",
*/
"  -seed=int",
"     Selects the seed of the random number generator.  ",
" ",
"  -dbglvl=int      ",
"     Selects the dbglvl.  ",
" ",
"  -help",
"     Prints this message.",
];

const SHORTHELPSTR: &[&str] = &[
    " ",
    "   Usage: gpmetis [options] <filename> <nparts>",
    "          use 'gpmetis -help' for a summary of the options.",
];

#[rustfmt::skip]
fn mkopt(p: &params_t) -> [idx_t; METIS_NOPTIONS as usize] {
    let mut o = [-1; _];
    o[METIS_OPTION_OBJTYPE as usize]   = p.objtype;
    o[METIS_OPTION_CTYPE as usize]     = p.ctype;
    o[METIS_OPTION_IPTYPE as usize]    = p.iptype;
    o[METIS_OPTION_RTYPE as usize]     = p.rtype;
    o[METIS_OPTION_NO2HOP as usize]    = p.no2hop as _;
    o[METIS_OPTION_ONDISK as usize]    = p.ondisk as _;
    o[METIS_OPTION_DROPEDGES as usize] = p.dropedges;
    o[METIS_OPTION_MINCONN as usize]   = p.minconn as _;
    o[METIS_OPTION_CONTIG as usize]    = p.contig as _;
    o[METIS_OPTION_SEED as usize]      = p.seed;
    o[METIS_OPTION_NIPARTS as usize]   = p.niparts;
    o[METIS_OPTION_NITER as usize]     = p.niter;
    o[METIS_OPTION_NCUTS as usize]     = p.ncuts;
    o[METIS_OPTION_UFACTOR as usize]   = p.ufactor;
    o[METIS_OPTION_DBGLVL as usize]    = p.dbglvl;
    o
}

fn error<T: std::fmt::Display>(msg: T) -> ExitCode {
    println!("{msg}");
    ExitCode::FAILURE
}

fn main() -> ExitCode {
    metis::set_bin_overrides();

    let mut params = params_t {
        ptype: METIS_PTYPE_KWAY as _,
        objtype: METIS_OBJTYPE_CUT as _,
        ctype: METIS_CTYPE_SHEM as _,
        nparts: 1,
        dropedges: 0,
        ncuts: 1,
        niter: 10,
        balance: 0,
        ..params_t::new()
    };
    let args: Vec<_> = std::env::args().collect();
    if let Some(code) = params::parse_standard(&args, HELPSTR, SHORTHELPSTR, OPTIONS, &mut params) {
        return code;
    }

    gk_startcputimer(&mut params.iotimer);

    let graph = unsafe { &mut *metis::graphio::ReadGraph(&params) };
    metis::graphio::ReadTPwgts(&mut params, graph.ncon as usize);

    gk_stopcputimer(&mut params.iotimer);

    // in cmdline_gpmetis.c
    if params.nparts < 2 {
        return error(format_args!(
            "The number of partitions must be greater than 1!"
        ));
    }

    if params.ptype == METIS_PTYPE_RB as _ {
        params.rtype = METIS_RTYPE_FM as _;
    } else if params.ptype == METIS_PTYPE_KWAY as _ {
        if params.iptype == -1 {
            params.iptype = METIS_IPTYPE_METISRB as _;
        }
        params.iptype = if params.iptype != -1 {
            params.iptype
        } else {
            METIS_IPTYPE_METISRB as _
        };
        params.rtype = METIS_RTYPE_GREEDY as _;
    } else {
        unreachable!();
    }

    // in gpmetis.c
    if params.contig && unsafe { metis::contig::IsConnected(graph, 0) == 0 } {
        println!(
            "***The input graph is not contiguous.\n\
            ***The specified -contig option will be ignored."
        );
        params.contig = false
    }

    if let Some(ubstr) = params.ubvecstr.as_deref() {
        let mut ubvec = vec![0.0; graph.ncon as usize];
        let mut ubs = ubstr.split_whitespace();
        for (i, ubi) in ubvec.iter_mut().enumerate() {
            let Some(ub) = ubs.next() else {
                return error(format_args!("Missing entry {i} of ubvec [{ubstr}]."));
            };
            let ub: real_t = match ub.parse() {
                Ok(ub) => ub,
                Err(e) => {
                    return error(format_args!(
                        "Error parsing entry {i} of ubvec [{ubstr}]: {e}"
                    ));
                }
            };
            *ubi = ub;
        }
        params.ubvec = Some(ubvec);
    }

    if params.ufactor == -1 {
        if params.ptype == METIS_PTYPE_KWAY as _ {
            params.ufactor = KMETIS_DEFAULT_UFACTOR as _
        } else if graph.ncon == 1 {
            params.ufactor = PMETIS_DEFAULT_UFACTOR as _
        } else {
            params.ufactor = MCPMETIS_DEFAULT_UFACTOR as _
        }
    }

    if params.iptype == -1
        && params.ptype == METIS_PTYPE_RB as _ {
            if graph.ncon == 1 {
                params.iptype = METIS_IPTYPE_GROW as _
            } else {
                params.iptype = METIS_IPTYPE_RANDOM as _
            }
        }

    print_gp_info(&params, graph);

    let options = mkopt(&params);
    let mut part = vec![0; graph.nvtxs as usize];

    gk_startcputimer(&mut params.parttimer);
    let mut objval = 0;
    let status = match params.ptype as _ {
        METIS_PTYPE_RB => unsafe {
            metis::pmetis::METIS_PartGraphRecursive(
                &graph.nvtxs,
                &graph.ncon,
                graph.xadj,
                graph.adjncy,
                graph.vwgt,
                graph.vsize,
                graph.adjwgt,
                &(params.nparts as idx_t),
                params.tpwgts.as_mut_ptr(),
                vec_ptr(&mut params.ubvec),
                options.as_ptr(),
                &mut objval,
                part.as_mut_ptr(),
            )
        },

        METIS_PTYPE_KWAY => unsafe {
            metis::kmetis::METIS_PartGraphKway(
                &graph.nvtxs,
                &graph.ncon,
                graph.xadj,
                graph.adjncy,
                graph.vwgt,
                graph.vsize,
                graph.adjwgt,
                &(params.nparts as idx_t),
                params.tpwgts.as_mut_ptr(),
                vec_ptr(&mut params.ubvec),
                options.as_ptr(),
                &mut objval,
                part.as_mut_ptr(),
            )
        },
        _ => unreachable!(),
    };
    gk_stopcputimer(&mut params.parttimer);

    // when we mix-and match with printf and we output to a pipe, things are annoying
    unsafe {
        libc::fflush(std::ptr::null_mut());
    }

    if status != METIS_OK {
        println!("\n***Metis returned with an error.");
    } else {
        if !params.nooutput {
            gk_startcputimer(&mut params.iotimer);
            metis::graphio::WritePartition(params.outfile.as_deref().unwrap_or(&params.filename), &part, params.nparts as idx_t);
            gk_stopcputimer(&mut params.iotimer);
        }

        report_gp_results(&mut params, graph, &part, objval);
    }

    let mut graph: *mut graph_t = graph;
    unsafe { metis::graph::FreeGraph(&mut graph) };

    if status == METIS_OK {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

fn print_gp_info(params: &params_t, graph: &graph_t) {
    println!("******************************************************************************");
    println!("{METISTITLE}");
    println!(
        " (HEAD: {}, Built on: {}, {})",
        SVNINFO, __DATE__, __TIME__
    );
    println!(
        " size of idx_t: {}bits, real_t: {}bits, idx_t *: {}bits",
        (idx_t::BITS as usize),
        8 * size_of::<real_t>(),
        8 * size_of::<&idx_t>()
    );
    println!();
    println!("Graph Information -----------------------------------------------------------");
    println!(
        " Name: {}, #Vertices: {}, #Edges: {}, #Parts: {}",
        params.filename.display(),
        graph.nvtxs,
        graph.nedges / 2,
        params.nparts
    );
    if graph.ncon > 1 {
        println!(" Balancing constraints: {}", graph.ncon);
    }

    println!();
    println!("Options ---------------------------------------------------------------------");
    println!(
        " ptype={}, objtype={}, ctype={}, rtype={}, iptype={}",
        PTYPE_NAMES[params.ptype as usize],
        OBJTYPE_NAMES[params.objtype as usize],
        CTYPE_NAMES[params.ctype as usize],
        RTYPE_NAMES[params.rtype as usize],
        IPTYPE_NAMES[params.iptype as usize]
    );

    println!(
        " dbglvl={}, ufactor={:.3}, no2hop={}, minconn={}, contig={}",
        params.dbglvl,
        i2rubfactor(params.ufactor),
        if params.no2hop { "YES" } else { "NO" },
        if params.minconn { "YES" } else { "NO" },
        if params.contig { "YES" } else { "NO" }
    );

    println!(
        " ondisk={}, nooutput={}",
        if params.ondisk { "YES" } else { "NO" },
        if params.nooutput { "YES" } else { "NO" }
    );

    println!(
        " seed={}, niparts={}, niter={}, ncuts={}",
        params.seed, params.niparts, params.niter, params.ncuts
    );

    if let Some(ubvec) = params.ubvec.as_deref() {
        print!(" ubvec=(");
        let mut it = ubvec.iter().copied();
        if let Some(i) = it.next() {
            print!("{i:.2}");
        }
        for i in it {
            print!(" {i:.2}");
        }
        println!(")");
    }

    println!();
    match params.ptype as _ {
        METIS_PTYPE_RB => println!(
            "Recursive Partitioning ------------------------------------------------------"
        ),
        METIS_PTYPE_KWAY => println!(
            "Direct k-way Partitioning ---------------------------------------------------"
        ),
        _ => unreachable!(),
    }
}

fn report_gp_results(params: &mut params_t, graph: &graph_t, part: &[idx_t], _objval: idx_t) {
    gk_startcputimer(&mut params.reporttimer);
    unsafe {
        metis::stat::ComputePartitionInfo(params, graph, part);
    }

    gk_stopcputimer(&mut params.reporttimer);

    if !cfg!(feature = "normalized") {
        print!("\nTiming Information ----------------------------------------------------------\n");
        println!(
            "  I/O:          \t\t {:7.3} sec",
            gk_getcputimer(params.iotimer)
        );
        println!(
            "  Partitioning: \t\t {:7.3} sec   (METIS time)",
            gk_getcputimer(params.parttimer)
        );
        println!(
            "  Reporting:    \t\t {:7.3} sec",
            gk_getcputimer(params.reporttimer)
        );
    }
    // print!("\nMemory Information ----------------------------------------------------------\n");
    // print!("  Max memory used:\t\t {:7.3} MB\n", (real_t)(params.maxmemory/(1024.0*1024.0)));

    // #if !defined(MACOS) && !defined(WIN32) && !defined(__MINGW32__)
    //   {
    //     struct rusage usage;
    //     getrusage(RUSAGE_SELF, &usage);
    //     print!("  rusage.ru_maxrss:\t\t {:7.3} MB\n", (real_t)(usage.ru_maxrss/(1024.0)));
    //   }
    //   print!("  proc/self/stat/VmPeak:\t {:7.3} MB\n", (real_t)gk_GetProcVmPeak()/(1024.0*1024.0));
    // #endif

    println!("******************************************************************************");
}
