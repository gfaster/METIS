//! # NOT YET TESTED
//!
//!
#![expect(non_snake_case)]
/*
* Copyright 1997, Regents of the University of Minnesota
*
* io.c
*
* This file contains routines related to I/O
*
* Started 8/28/94
* George
*
* $Id: io.c 17513 2014-08-05 16:20:50Z dominique $
*
*/

// const MAXLINE: usize = 1280000;

use std::{
    ffi::OsString,
    fs::OpenOptions,
    io::{self, prelude::*, BufReader, Cursor},
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::{params::params_t, scanf::Scanf, *};

fn strtoidx(s: &mut &str) -> Result<idx_t, <idx_t as FromStr>::Err> {
    *s = s.trim_start();
    let mut first = true;
    let end = s
        .bytes()
        .take_while(|&b| {
            if first {
                first = false;
                if b == b'-' {
                    return true;
                }
            }
            b.is_ascii_digit()
        })
        .count();
    let ret = idx_t::from_str_radix(&s[..end], 10);
    *s = &s[end..];
    ret
}

fn strtoreal(s: &mut &str) -> Result<real_t, <real_t as FromStr>::Err> {
    *s = s.trim_start();
    let mut first = true;
    let mut decimal = false;
    let end = s
        .bytes()
        .take_while(|&b| {
            if first {
                first = false;
                if b == b'-' {
                    return true;
                }
            }
            if !decimal {
                if b == b'.' {
                    decimal = true;
                    return true;
                }
            }
            b.is_ascii_digit()
        })
        .count();
    let ret = s[..end].parse();
    *s = &s[end..];
    ret
}

fn extend_filename(p: &Path, ext: impl std::fmt::Display) -> PathBuf {
    // this is bad to preserve original semantics as much as possible
    let mut filename: Vec<u8> = Vec::new();
    filename.extend_from_slice(p.as_os_str().as_encoded_bytes());
    write!(filename, "{ext}").unwrap();
    let filename = unsafe { OsString::from_encoded_bytes_unchecked(filename) };
    PathBuf::from(filename)
}

/*************************************************************************/
/* This function reads in a sparse graph */
/*************************************************************************/
pub unsafe fn ReadGraph(params: &params_t) -> *mut graph_t {
    // idx_t i, j, k, l, fmt, ncon, nfields, readew, readvw, readvs, edge, ewgt;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *vsize;

    // char *line=std::ptr::null_mut(), fmtstr[256 as usize], *curstr, *newstr;
    // size_t lnlen=0;
    // FILE *fpin;
    // graph_t *graph;

    // if (!gk_fexists(params.filename))
    //   {
    //       errexit("File %s does not exist!", params.filename);
    //   }

    let mut fpin = match OpenOptions::new().read(true).open(&params.filename) {
        Ok(f) => BufReader::new(f),
        Err(e) => panic!("Failed to open {:?}: {e}", params.filename),
    };

    let graphp = graph::CreateGraph();
    let graph: &mut graph_t = graphp.as_mut().unwrap();

    // fpin = gk_fopen(params.filename, "r", "ReadGRaph: Graph");

    /* Skip comment lines until you get to the first valid line */

    // do {
    //   if (gk_getline(&line, &lnlen, fpin) == -1)
    //   {
    //       errexit("Premature end of input file: file: %s", params.filename);
    //   }
    // } while (line[0 as usize] == '%');
    let mut line = String::from('%');
    while line.starts_with('%') {
        line.clear();
        if 0 == fpin.read_line(&mut line).unwrap() {
            panic!("Premature end of input file: {:?}", params.filename)
        }
    }

    // fmt = ncon = 0;
    let mut fmt: idx_t = 0;
    let mut ncon: idx_t = 0;

    let mut cur = Cursor::new(line.as_bytes());
    graph.nvtxs = cur
        .scan()
        .expect("The input file does not specify the number of vertices");
    let _ = cur.scanstr(" ");
    graph.nedges = cur
        .scan()
        .expect("The input file does not specify the number of edges");
    let _ = cur.scanstr(" ");
    if let Ok(f) = cur.scan() {
        fmt = f;
    }
    let _ = cur.scanstr(" ");
    if let Ok(n) = cur.scan() {
        ncon = n;
        drop(cur)
    }

    // nfields = sscanf(line, "%"SCIDX" %"SCIDX" %"SCIDX" %"SCIDX, &(graph.nvtxs), &(graph.nedges), &fmt, &ncon);

    if graph.nvtxs <= 0 || graph.nedges <= 0 {
        panic!(
            "The supplied nvtxs:{:} and nedges:{:} must be positive.",
            graph.nvtxs, graph.nedges
        );
    }

    if fmt > 111 {
        panic!("Cannot read this type of file format [fmt={fmt:}]!");
    }

    // write!(fmtstr, "{:03}", fmt%1000);
    // readvs = (fmtstr[0 as usize] == '1');
    // readvw = (fmtstr[1 as usize] == '1');
    // readew = (fmtstr[2 as usize] == '1');
    let readvs = (fmt / 100) % 10 == 1;
    let readvw = (fmt / 10) % 10 == 1;
    let readew = (fmt / 1) % 10 == 1;

    /*println!("{:s %} {:} {:}", fmtstr, readvs, readvw, readew); */

    if ncon > 0 && !readvw {
        panic!(
            concat!(
                "------------------------------------------------------------------------------\n",
                "***  I detected an error in your input file  ***\n",
                "You specified ncon={:}, but the fmt parameter does not specify vertex weights\n",
                "Make sure that the fmt parameter is set to either 10 or 11.\n",
                "------------------------------------------------------------------------------\n"
            ),
            ncon
        );
    }

    graph.nedges *= 2;
    graph.ncon = if ncon == 0 { 1 } else { ncon };
    let ncon = graph.ncon as usize;

    graph.xadj = ismalloc(graph.nvtxs as usize + 1, 0, c"ReadGraph: xadj".as_ptr()) as *mut idx_t;
    graph.adjncy = imalloc(graph.nedges as usize, c"ReadGraph: adjncy".as_ptr()) as *mut idx_t;
    graph.vwgt =
        ismalloc(ncon * graph.nvtxs as usize, 1, c"ReadGraph: vwgt".as_ptr()) as *mut idx_t;
    graph.adjwgt = if readew {
        imalloc(graph.nedges as usize, c"ReadGraph: adjwgt".as_ptr()) as *mut idx_t
    } else {
        std::ptr::null_mut()
    };
    graph.vsize = ismalloc(graph.nvtxs as usize, 1, c"ReadGraph: vsize".as_ptr()) as *mut idx_t;
    get_graph_slices_mut!(graph => xadj adjncy vsize vwgt);
    get_graph_slices_optional!(mut graph => adjwgt);

    /*----------------------------------------------------------------------
     * Read the sparse graph file
     *---------------------------------------------------------------------*/
    xadj[0] = 0;
    let mut k = 0;

    for i in (0)..(graph.nvtxs as usize) {
        line = String::from("%");
        while line.starts_with('%') {
            line.clear();
            if 0 == fpin.read_line(&mut line).unwrap() {
                panic!("Premature end of input file while reading vertex {}", i + 1)
            }
        }
        //     do {
        //       if (gk_getline(&line, &lnlen, fpin) == -1)
        // {
        //     errexit("Premature end of input file while reading vertex {:}.", i+1);
        // }
        //     } while (line[0 as usize] == '%');

        let mut curstr = line.as_ref();

        /* Read vertex sizes */
        if readvs {
            vsize[i as usize] = match strtoidx(&mut curstr) {
                Ok(i) => i,
                Err(e) => panic!(
                    "The line for vertex {} does not have vsize information [{e}]",
                    i + 1
                ),
            };
            if vsize[i as usize] < 0 {
                panic!("The size for vertex {:} must be >= 0", i + 1);
            }
        }

        /* Read vertex weights */
        if readvw {
            for l in (0)..(ncon) {
                vwgt[(i*ncon+l) as usize] = match strtoidx(&mut curstr) {
                    Ok(i) => i,
                    Err(e) => panic!("The line for vertex {:} does not have enough weights for the {:} constraints. [{e}]", i+1, ncon)
                };
                if vwgt[(i * ncon + l) as usize] < 0 {
                    panic!(
                        "The weight vertex {:} and constraint {:} must be >= 0",
                        i + 1,
                        l
                    );
                }
            }
        }

        loop {
            let Ok(edge) = strtoidx(&mut curstr) else {
                break; /* End of line */
            };
            if edge < 1 || edge > graph.nvtxs {
                panic!("Edge {:} for vertex {:} is out of bounds", edge, i + 1);
            }

            let mut ewgt = 1;
            if readew {
                ewgt = match strtoidx(&mut curstr) {
                    Ok(i) => i,
                    Err(e) => panic!("Premature end of line for vertex {} [{e}]", i + 1),
                };
                if ewgt <= 0 {
                    panic!(
                        "The weight ({:}) for edge ({:}, {:}) must be positive.",
                        ewgt,
                        i + 1,
                        edge
                    );
                }
            }

            if k == graph.nedges {
                panic!(
                    "There are more edges in the file than the {:} specified.",
                    graph.nedges / 2
                );
            }

            adjncy[k as usize] = edge - 1;
            if readew {
                adjwgt.as_mut().unwrap()[k as usize] = ewgt
            };
            {
                k += 1;
            }
        }
        xadj[(i + 1) as usize] = k;
    }
    // gk_fclose(fpin);
    drop(fpin);

    if k != graph.nedges {
        println!("------------------------------------------------------------------------------");
        println!("***  I detected an error in your input file  ***\n");
        println!(
            "In the first line of the file, you specified that the graph contained \
                {:} edges. However, I only found {:} edges in the file.",
            graph.nedges / 2,
            k / 2
        );
        if 2 * k == graph.nedges {
            println!(
                "\n *> I detected that you specified twice the number of edges that you have in"
            );
            println!("    the file. Remember that the number of edges specified in the first line");
            println!("    counts each edge between vertices v and u only once.\n");
        }
        println!("Please specify the correct number of edges in the first line of the file.");
        println!("------------------------------------------------------------------------------");
        panic!("invalid input file")
    }

    // gk_free((void *)&line, LTERM);

    return graphp;
}

/*************************************************************************/
/* This function reads in a mesh */
/*************************************************************************/

// pub fn ReadMesh(params: &params_t) -> *mut mesh_t
// {
//   // idx_t i, j, k, l, nfields, ncon, node;
//   // idx_t *eptr, *eind, *ewgt;
//   // size_t nlines, ntokens;
//   // char *line=std::ptr::null_mut(), *curstr, *newstr;
//   // size_t lnlen=0;
//   // FILE *fpin;
//   // mesh_t *mesh;
//   if (!gk_fexists(params.filename))
//     {
//         errexit("File %s does not exist!", params.filename);
//     }
//
//   mesh = CreateMesh();
//
//   /* get some file stats */
//   gk_getfilestats(params.filename, &nlines, &ntokens, std::ptr::null_mut(), std::ptr::null_mut());
//
//   fpin = gk_fopen(params.filename, "r", __func__);
//
//   /* Skip comment lines until you get to the first valid line */
//   do {
//     if (gk_getline(&line, &lnlen, fpin) == -1)
//     {
//         errexit("Premature end of input file: file: %s", params.filename);
//     }
//   } while (line[0 as usize] == '%');
//
//
//   mesh.ncon = 0;
//   nfields = sscanf(line, "%"SCIDX" %"SCIDX, &(mesh.ne), &(mesh.ncon));
//
//   if (nfields < 1)
// {
// errexit("The input file does not specify the number of elements.");
//     }
//
//   if (mesh.ne <= 0)
// {
//     errexit("The supplied number of elements:{:} must be positive.", mesh.ne);
// }
//
//   if (mesh.ne > nlines)
//     errexit("The file has %zu lines which smaller than the number of "
//             "elements of {:} specified in the header line.", nlines, mesh.ne);
//
//   ncon = mesh.ncon;
//   eptr = mesh.eptr = ismalloc(mesh.ne+1, 0, "ReadMesh: eptr");
//   eind = mesh.eind = imalloc(ntokens, "ReadMesh: eind");
//   ewgt = mesh.ewgt = ismalloc((ncon == 0 ? 1 : ncon)*mesh.ne, 1, "ReadMesh: ewgt");
//
//
//   /*----------------------------------------------------------------------
//    * Read the mesh file
//    *---------------------------------------------------------------------*/
//   for (eptr[0 as usize]=0, k=0, i=0; i<mesh.ne; i++) {
//     do {
//       if (gk_getline(&line, &lnlen, fpin) == -1)
// {
//     errexit("Premature end of input file while reading element {:}.", i+1);
// }
//     } while (line[0 as usize] == '%');
//
//     curstr = line;
//     newstr = std::ptr::null_mut();
//
//     /* Read element weights */
//     for l in (0)..(ncon) {
//       ewgt[(i*ncon+l) as usize] = strtoidx(curstr, &newstr, 10);
//       if (newstr == curstr)
//         errexit("The line for vertex {:} does not have enough weights "
//                 "for the {:} constraints.", i+1, ncon);
//       if (ewgt[(i*ncon+l) as usize] < 0)
// {
// errexit("The weight for element {:} and constraint {:} must be >= 0", i+1, l);
//         }
//       curstr = newstr;
//     }
//
//     while (1) {
//       node = strtoidx(curstr, &newstr, 10);
//       if (newstr == curstr)
//     {
//         break; /* End of line */
//     }
//       curstr = newstr;
//
//       if (node < 1)
//     {
//         errexit("Node {:} for element {:} is out of bounds", node, i+1);
//     }
//
//       eind[k as usize] = node-1;k+=1;
//     }
//     eptr[(i+1) as usize] = k;
//   }
//   gk_fclose(fpin);
//
//   mesh.ncon = (if ncon == 0  {  1  } else {  ncon });
//   mesh.nn   = imax(eptr[mesh.ne], eind, 1)+1;
//
//   gk_free((void *)&line, LTERM);
//
//   return mesh;
// }

/*************************************************************************/
/* This function reads in the target partition weights. If no file is
specified the weights are set to 1/nparts */
/*************************************************************************/
pub unsafe fn ReadTPwgts(params: &mut params_t, ncon: usize) {
    // idx_t i, j, from, to, fromcnum, tocnum, nleft;
    // real_t awgt=0.0, twgt;
    // char *line=std::ptr::null_mut(), *curstr, *newstr;
    // size_t lnlen=0;
    // FILE *fpin;

    params.tpwgts =
        crate::rsmalloc(params.nparts * ncon, -1.0, c"ReadTPwgts: tpwgts".as_ptr()) as _;

    let Some(tpwgtsfile) = params.tpwgtsfile.as_deref() else {
        mkslice_mut!(params->tpwgts, params.nparts * ncon);
        for i in (0)..(params.nparts) {
            for j in (0)..(ncon) {
                tpwgts[(i * ncon + j) as usize] = 1.0 / (params.nparts as real_t);
            }
        }
        return;
    };

    let fpin = match OpenOptions::new().read(true).open(tpwgtsfile) {
        Ok(f) => f,
        Err(e) => panic!("Could not open graph file {:?}: {e}", tpwgtsfile),
    };
    let mut fpin = BufReader::new(fpin);
    //   if (!gk_fexists(params.tpwgtsfile))
    // {
    //     errexit("Graph file %s does not exist!", params.tpwgtsfile);
    // }

    // fpin = gk_fopen(params.tpwgtsfile, "r", "ReadTPwgts: tpwgtsfile");

    // while (gk_getline(&line, &lnlen, fpin) != -1) {
    let mut line = String::new();
    while fpin.read_line(&mut line).unwrap() != 0 {
        // gk_strchr_replace(line, " ", "");
        let line = line.replace(" ", "");

        /* start extracting the fields */

        let mut curstr = line.as_str();
        let from = match strtoidx(&mut curstr) {
            Ok(x) => x,
            Err(e) => panic!(
                "The 'from' component of line <{line}> in the tpwgts file is incorrect [{e}]"
            ),
        };
        //     from = strtoidx(curstr, &newstr, 10);
        //     if (newstr == curstr)
        // {
        //     errexit("The 'from' component of line <%s> in the tpwgts file is incorrect.", line);
        // }
        //     curstr = newstr;

        let to;
        if let Some(s) = curstr.strip_prefix('-') {
            curstr = s;
            to = match strtoidx(&mut curstr) {
                Ok(x) => x,
                Err(e) => panic!(
                    "The 'to' component of line <{line}> in the tpwgts file is incorrect. [{e}]"
                ),
            };
        } else {
            to = from;
        }

        let fromcnum;
        let tocnum;
        if let Some(s) = curstr.strip_prefix(':') {
            curstr = s;
            fromcnum = match strtoidx(&mut curstr) {
                Ok(x) => x,
                Err(e) => panic!("The 'fromcnum' component of line <{line}> in the tpwgts file is incorrect. [{e}]")
            };

            // if (curstr[0 as usize] == '-') {
            if let Some(s) = curstr.strip_prefix('-') {
                curstr = s;

                tocnum = match strtoidx(&mut curstr) {
                    Ok(x) => x,
                    Err(e) => panic!("The 'tocnum' component of line <{line}> in the tpwgts file is incorrect. [{e}]")
                };
            } else {
                tocnum = fromcnum;
            }
        } else {
            fromcnum = 0;
            tocnum = ncon as idx_t - 1;
        }

        let awgt;
        if let Some(s) = curstr.strip_prefix('=') {
            curstr = s;
            awgt = match strtoreal(&mut curstr) {
                Ok(x) => x,
                Err(e) => panic!(
                    "The 'wgt' component of line <{line}> in the tpwgts file is incorrect. [{e}]"
                ),
            };
        } else {
            panic!("The 'wgt' component of line <{line}> in the tpwgts file is missing.");
        }

        /*println!("Read: {:}-{:}:{:}-{:}={:}",
        from, to, fromcnum, tocnum, awgt);*/

        if from < 0 || to < 0 || from >= params.nparts as idx_t || to >= params.nparts as idx_t {
            panic!("Invalid partition range for {:}:{:}", from, to);
        }
        if fromcnum < 0 || tocnum < 0 || fromcnum as usize >= ncon || tocnum as usize >= ncon {
            panic!(
                "Invalid constraint number range for {:}:{:}",
                fromcnum, tocnum
            );
        }
        if awgt <= 0.0 || awgt >= 1.0 {
            panic!("Invalid partition weight of {:}", awgt);
        }
        mkslice_mut!(params->tpwgts, params.nparts * ncon);
        for i in (from)..=(to) {
            for j in (fromcnum)..=(tocnum) {
                tpwgts[(i as usize * ncon + j as usize) as usize] = awgt;
            }
        }
    }

    // gk_fclose(fpin);
    drop(fpin);

    /* Assign weight to the unspecified constraints x partitions */
    for j in (0)..(ncon) {
        /* Sum up the specified weights for the jth constraint */
        let mut twgt = 0.0;
        let mut nleft = params.nparts;
        mkslice!(params->tpwgts, params.nparts * ncon);
        for i in (0)..(params.nparts) {
            if tpwgts[(i * ncon + j) as usize] > 0.0 {
                twgt += tpwgts[(i * ncon + j) as usize];
                nleft -= 1;
            }
        }

        /* Rescale the weights to be on the safe side */
        if nleft == 0 {
            mkslice_mut!(params->tpwgts, params.nparts * ncon);
            util::rscale(params.nparts, 1.0 / twgt, &mut tpwgts[j..], ncon);
        }

        /* Assign the left-over weight to the remaining partitions */
        if nleft > 0 {
            if twgt > 1.0 {
                panic!(
                    "The total specified target partition weights for constraint #{} \
                        of {:} exceeds 1.0.",
                    j, twgt
                );
            }

            let awgt = (1.0 - twgt) / (nleft as real_t);
            for i in (0)..(params.nparts) {
                mkslice_mut!(params->tpwgts, params.nparts * ncon);
                tpwgts[(i * ncon + j) as usize] = if tpwgts[(i * ncon + j) as usize] < 0.0 {
                    awgt
                } else {
                    tpwgts[(i * ncon + j) as usize]
                };
            }
        }
    }
}

/*************************************************************************/
/* This function reads in a partition/ordering vector  */
/**************************************************************************/
pub unsafe fn ReadPOVector(graph: &mut graph_t, filename: &Path, vector: *mut idx_t) {
    // idx_t i;
    // FILE *fpin;

    // fpin = gk_fopen(filename, "r", __func__);
    let mut fpin = BufReader::new(OpenOptions::new().read(true).open(filename).unwrap());
    mkslice_mut!(vector, graph.nvtxs);
    let mut buf = String::new();
    let mut i = 0;
    while fpin.read_line(&mut buf).unwrap() != 0 {
        vector[i] = match strtoidx(&mut &*buf) {
            Ok(x) => x,
            Err(e) => panic!(
                "Failed to read {filename:?} at line {i} [(nvtxs: {})] [{e}]",
                graph.nvtxs
            ),
        };
        i += 1;
        if i >= graph.nvtxs as usize {
            break;
        }
        buf.clear();
    }
    if i < graph.nvtxs as usize {
        panic!(
            "Premature end of file {filename:?} at line {i} [(nvtxs: {})]",
            graph.nvtxs
        );
    }
    // for i in (0)..((*graph).nvtxs) {
    //   vector[i] =
    //   if (fscanf(fpin, "%"SCIDX"", vector+i) != 1)
    //     errexit("[%s] Premature end of file %s at line {} [(nvtxs: {}) as usize]",
    //         __func__, filename, i, graph.nvtxs);
    // }
}

/*************************************************************************/
/* This function writes out the partition vector */
/*************************************************************************/
pub unsafe fn WritePartition(fname: &Path, part: *mut idx_t, n: idx_t, nparts: idx_t) {
    // FILE *fpout;
    // idx_t i;
    // char filename[MAXLINE as usize];

    // sprintf(filename, "%s.part.%"PRIDX, fname, nparts);
    let filename = extend_filename(fname, format_args!(".part.{nparts}"));

    // fpout = gk_fopen(filename, "w", __func__);
    let mut fpout = OpenOptions::new().write(true).open(filename).unwrap();

    mkslice!(part, n);
    for i in (0)..(n) {
        writeln!(fpout, "{}", part[i as usize]).unwrap();
    }

    // gk_fclose(fpout);
}

/*************************************************************************/
/* This function writes out the partition vectors for a mesh */
/*************************************************************************/
pub unsafe fn WriteMeshPartition(
    fname: &Path,
    nparts: idx_t,
    ne: idx_t,
    epart: *mut idx_t,
    nn: idx_t,
    npart: *mut idx_t,
) {
    // FILE *fpout;
    // idx_t i;
    // char filename[256 as usize];

    // write!(filename,"%s.epart.%"PRIDX,fname, nparts);
    let filename = extend_filename(fname, format_args!(".epart.{nparts}"));

    let mut fpout = OpenOptions::new().write(true).open(filename).unwrap();
    // fpout = gk_fopen(filename, "w", __func__);

    mkslice!(epart, ne);
    for i in (0)..(ne) {
        writeln!(fpout, "{}", epart[i as usize]).unwrap();
    }

    // gk_fclose(fpout);
    drop(fpout);

    // write!(filename,"%s.npart.%"PRIDX,fname, nparts);
    let filename = extend_filename(fname, format_args!(".npart.{nparts}"));
    let mut fpout = OpenOptions::new().write(true).open(filename).unwrap();

    // fpout = gk_fopen(filename, "w", __func__);

    mkslice!(npart, nn);
    for i in (0)..(nn) {
        writeln!(fpout, "{}", npart[i as usize]).unwrap();
    }

    // gk_fclose(fpout);
}

/*************************************************************************/
/* This function writes out the permutation vector */
/*************************************************************************/
pub unsafe fn WritePermutation(fname: &Path, iperm: *mut idx_t, n: idx_t) {
    // FILE *fpout;
    // idx_t i;
    // char filename[MAXLINE as usize];

    // write!(filename, "%s.iperm", fname);
    let filename = extend_filename(fname, ".iperm");

    // fpout = gk_fopen(filename, "w", __func__);
    let mut fpout = OpenOptions::new().write(true).open(filename).unwrap();

    mkslice!(iperm, n);
    for i in (0)..(n) {
        writeln!(fpout, "{}", iperm[i as usize]).unwrap();
    }

    // gk_fclose(fpout);
}

/*************************************************************************/
/* This function writes a graph into a file  */
/*************************************************************************/
pub unsafe fn WriteGraph(graph: &graph_t, filename: &Path) {
    // idx_t i, j, nvtxs, ncon;
    // idx_t *xadj, *adjncy, *adjwgt, *vwgt, *vsize;
    // int hasvwgt=0, hasewgt=0, hasvsize=0;
    // FILE *fpout;
    let mut hasvwgt = false;
    let mut hasewgt = false;
    let mut hasvsize = false;
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_optional!(graph => vwgt vsize adjwgt);

    /* determine if the graph has non-unity vwgt, vsize, or adjwgt */
    if let Some(vwgt) = vwgt {
        for i in (0)..(nvtxs * ncon) {
            if vwgt[i as usize] != 1 {
                hasvwgt = true;
                break;
            }
        }
    }
    if let Some(vsize) = vsize {
        for i in (0)..(nvtxs) {
            if vsize[i as usize] != 1 {
                hasvsize = true;
                break;
            }
        }
    }
    if let Some(adjwgt) = adjwgt {
        for i in (0)..(xadj[nvtxs as usize]) {
            if adjwgt[i as usize] != 1 {
                hasewgt = true;
                break;
            }
        }
    }

    // fpout = gk_fopen(filename, "w", __func__);
    let mut fpout = std::fs::OpenOptions::new()
        .write(true)
        .open(filename)
        .unwrap();

    /* write the header line */
    write!(fpout, "{:} {}", nvtxs, xadj[nvtxs as usize] / 2).unwrap();
    if hasvwgt || hasvsize || hasewgt {
        write!(
            fpout,
            " {}{}{}",
            hasvsize as i8, hasvwgt as i8, hasewgt as i8
        )
        .unwrap();
        if hasvwgt {
            write!(fpout, " {}", graph.ncon).unwrap();
        }
    }
    // note no trailing newline

    /* write the rest of the graph */
    for i in (0)..(nvtxs) {
        writeln!(fpout).unwrap();
        if let Some(vsize) = vsize {
            writeln!(fpout, " {}", vsize[i as usize]).unwrap();
        }

        if let Some(vwgt) = vwgt {
            for j in (0)..(ncon) {
                write!(fpout, " {}", vwgt[(i * ncon + j) as usize]).unwrap();
            }
        }

        for j in (xadj[i as usize])..(xadj[i as usize]) {
            write!(fpout, " {}", adjncy[j as usize] + 1).unwrap();
            if let Some(adjwgt) = adjwgt {
                write!(fpout, " {}", adjwgt[j as usize]).unwrap();
            }
        }
    }

    // gk_fclose(fpout);
}
