use std::io::BufRead;
use std::io::Read;
use std::io::Write;
use std::os::raw::c_void;

/**
\file
\brief Functions that deal with setting up the graphs for METIS.

\date   Started 7/25/1997
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version\verbatim $Id: graph.c 15817 2013-11-25 14:58:41Z karypis $ \endverbatim
*/
use crate::*;

trait ReadTransmuted {
    fn read_idx(&mut self, buf: &mut [idx_t]) -> std::io::Result<()>;
    fn read_real(&mut self, buf: &mut [real_t]) -> std::io::Result<()>;
}

impl<T: Read> ReadTransmuted for T {
    fn read_idx(&mut self, buf: &mut [idx_t]) -> std::io::Result<()> {
        unsafe {
            let len = buf.len();
            let p = buf.as_mut_ptr() as *mut u8;
            let buf = std::slice::from_raw_parts_mut(p, len * std::mem::size_of::<idx_t>());
            self.read_exact(buf)
        }
    }

    fn read_real(&mut self, buf: &mut [real_t]) -> std::io::Result<()> {
        unsafe {
            let len = buf.len();
            let p = buf.as_mut_ptr() as *mut u8;
            let buf = std::slice::from_raw_parts_mut(p, len * std::mem::size_of::<real_t>());
            self.read_exact(buf)
        }
    }
}

trait WriteTransmuted {
    fn write_idx(&mut self, buf: &[idx_t]) -> std::io::Result<()>;
    fn write_real(&mut self, buf: &[real_t]) -> std::io::Result<()>;
}

impl<T: std::io::Write> WriteTransmuted for T {
    fn write_idx(&mut self, buf: &[idx_t]) -> std::io::Result<()> {
        unsafe {
            let len = buf.len();
            let p = buf.as_ptr() as *const u8;
            let buf = std::slice::from_raw_parts(p, len * std::mem::size_of::<idx_t>());
            self.write_all(buf)
        }
    }

    fn write_real(&mut self, buf: &[real_t]) -> std::io::Result<()> {
        unsafe {
            let len = buf.len();
            let p = buf.as_ptr() as *const u8;
            let buf = std::slice::from_raw_parts(p, len * std::mem::size_of::<real_t>());
            self.write_all(buf)
        }
    }
}

/*************************************************************************/
/* This function sets up the graph from the user input */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupGraph(
    ctrl: *mut ctrl_t,
    nvtxs: idx_t,
    ncon: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    adjwgt: *mut idx_t,
) -> *mut graph_t {
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, sum;
    // real_t *nvwgt;
    // graph_t *graph;

    /* allocate the graph and fill in the fields */
    let graph = CreateGraph();
    let mut adjwgt = adjwgt;
    let mut vsize = vsize;
    let mut vwgt = vwgt;
    let adjncy = adjncy;

    {
        let graph = graph.as_mut().unwrap();

        graph.nvtxs = nvtxs;
        graph.ncon = ncon;
        let nvtxs = nvtxs as usize;
        let ncon = ncon as usize;

        graph.nedges = *xadj.add(nvtxs);

        graph.xadj = xadj;
        graph.free_xadj = 0;

        graph.adjncy = adjncy;
        graph.free_adjncy = 0;

        graph.droppedewgt = 0;

        /* setup the vertex weights */
        if !vwgt.is_null() {
            (*graph).vwgt = vwgt;
            (*graph).free_vwgt = 0;
        } else {
            vwgt = ismalloc(ncon * nvtxs, 1, "SetupGraph: vwgt\0".as_ptr()) as *mut idx_t;
            (*graph).vwgt = vwgt;
        }

        (*graph).tvwgt = imalloc(ncon, "SetupGraph: tvwgts\0".as_ptr()) as *mut idx_t;
        (*graph).invtvwgt = rmalloc(ncon, "SetupGraph: invtvwgts\0".as_ptr()) as *mut real_t;
        for i in (0)..(ncon) {
            mkslice_mut!(vwgt, ncon * nvtxs);
            // (*graph).tvwgt[i]    = isum(nvtxs, vwgt+i, ncon);
            *(*graph).tvwgt.add(i) = vwgt.iter().skip(i).step_by(ncon).sum::<idx_t>();
            *(*graph).invtvwgt.add(i) = 1.0
                / (if *(*graph).tvwgt.add(i) > 0 {
                    *(*graph).tvwgt.add(i) as real_t
                } else {
                    1.0
                });
        }

        if ctrl.objtype == METIS_OBJTYPE_VOL {
            /* Setup the vsize */
            if !vsize.is_null() {
                (*graph).vsize = vsize;
                (*graph).free_vsize = 0;
            } else {
                vsize = ismalloc(nvtxs, 1, "SetupGraph: vsize\0".as_ptr()) as *mut idx_t;
                (*graph).vsize = vsize;
            }

            /* Allocate memory for edge weights and initialize them to the sum of the vsize */
            adjwgt = imalloc(graph.nedges as usize, "SetupGraph: adjwgt\0".as_ptr()) as *mut idx_t;
            (*graph).adjwgt = adjwgt;
            {
                get_graph_slices!(graph => xadj vsize adjncy);
                get_graph_slices_mut!(graph => adjwgt);
                for i in (0)..(nvtxs) {
                    for j in (xadj[i] as usize)..(xadj[(i + 1) as usize] as usize) {
                        adjwgt[j] = 1 + vsize[i] + vsize[adjncy[j] as usize];
                    }
                }
            }
        } else {
            /* For edgecut minimization */
            /* setup the edge weights */
            if !adjwgt.is_null() {
                graph.adjwgt = adjwgt;
                graph.free_adjwgt = 0;
            } else {
                adjwgt = ismalloc(graph.nedges as usize, 1, "SetupGraph: adjwgt\0".as_ptr())
                    as *mut idx_t;
                graph.adjwgt = adjwgt;
            }
        }

        /* setup various derived info */
        SetupGraph_tvwgt(graph);

        if ctrl.optype == METIS_OP_PMETIS || ctrl.optype == METIS_OP_OMETIS {
            SetupGraph_label(graph);
        }

        assert!(CheckGraph(graph, ctrl.numflag, 1) != 0);
    }

    return graph;
}

/*************************************************************************/
/* Set's up the tvwgt/invtvwgt info */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupGraph_tvwgt(graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    // idx_t i;

    let ncon = graph.ncon as usize;
    let _nvtxs = graph.nvtxs as usize;
    if graph.tvwgt == std::ptr::null_mut() {
        graph.tvwgt = imalloc(ncon, "SetupGraph_tvwgt: tvwgt\0".as_ptr()) as *mut idx_t;
    }
    if graph.invtvwgt == std::ptr::null_mut() {
        graph.invtvwgt = rmalloc(ncon, "SetupGraph_tvwgt: invtvwgt\0".as_ptr()) as *mut real_t;
    }
    get_graph_slices!(graph => vwgt);
    get_graph_slices_mut!(graph => tvwgt invtvwgt);

    for i in (0)..(ncon) {
        // tvwgt[i] = isum(graph.nvtxs, graph.vwgt + i, graph.ncon);
        // invtvwgt[i] = 1.0 / (if graph.tvwgt[i] > 0 { graph.tvwgt[i] } else { 1 });
        tvwgt[i] = vwgt.iter().skip(i).step_by(ncon).sum();
        invtvwgt[i] = 1.0
            / (if tvwgt[i] > 0 {
                tvwgt[i] as real_t
            } else {
                1.0
            });
    }
}

/*************************************************************************/
/* Set's up the label info */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupGraph_label(graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    // idx_t i;

    let nvtxs = graph.nvtxs as usize;
    if graph.label == std::ptr::null_mut() {
        graph.label = imalloc(nvtxs, "SetupGraph_label: label\0".as_ptr()) as *mut idx_t;
    }
    get_graph_slices_mut!(graph => label);

    for i in (0)..(nvtxs) {
        label[i] = i as idx_t;
    }
}

/*************************************************************************/
/* Setup the various arrays for the split graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupSplitGraph(
    graph: *mut graph_t,
    snvtxs: idx_t,
    snedges: idx_t,
) -> *mut graph_t {
    let graph = graph.as_mut().unwrap();
    let sgraph = CreateGraph();

    {
        let sgraph = sgraph.as_mut().unwrap();
        sgraph.nvtxs = snvtxs;
        sgraph.nedges = snedges;
        sgraph.ncon = graph.ncon;
        let snvtxs = snvtxs as usize;
        let snedges = snedges as usize;
        let sncon = sgraph.ncon as usize;

        /* Allocate memory for the split graph */
        sgraph.xadj = imalloc(snvtxs + 1, "SetupSplitGraph: xadj\0".as_ptr()) as *mut idx_t;
        sgraph.vwgt = imalloc(sncon * snvtxs, "SetupSplitGraph: vwgt\0".as_ptr()) as *mut idx_t;
        sgraph.adjncy = imalloc(snedges, "SetupSplitGraph: adjncy\0".as_ptr()) as *mut idx_t;
        sgraph.adjwgt = imalloc(snedges, "SetupSplitGraph: adjwgt\0".as_ptr()) as *mut idx_t;
        sgraph.label = imalloc(snvtxs, "SetupSplitGraph: label\0".as_ptr()) as *mut idx_t;
        sgraph.tvwgt = imalloc(sncon, "SetupSplitGraph: tvwgt\0".as_ptr()) as *mut idx_t;
        sgraph.invtvwgt = rmalloc(sncon, "SetupSplitGraph: invtvwgt\0".as_ptr()) as *mut real_t;

        if !graph.vsize.is_null() {
            sgraph.vsize = imalloc(snvtxs, "SetupSplitGraph: vsize\0".as_ptr()) as *mut idx_t;
        }
    }

    return sgraph;
}

/*************************************************************************/
/* This function creates and initializes a graph_t data structure */
/*************************************************************************/
#[metis_func]
pub extern "C" fn CreateGraph() -> *mut graph_t {
    let graph = gk_malloc(
        std::mem::size_of::<graph_t>(),
        "CreateGraph: graph\0".as_ptr(),
    ) as *mut graph_t;

    InitGraph(graph);

    return graph;
}

/*************************************************************************/
/* This function initializes a graph_t data structure */
/*************************************************************************/
#[metis_func]
pub extern "C" fn InitGraph(graph: *mut graph_t) {
    // memset((void *)graph, 0, sizeof(graph_t));
    *graph = std::mem::zeroed();
    let graph = graph.as_mut().unwrap();

    /* graph size constants */
    graph.nvtxs = -1;
    graph.nedges = -1;
    graph.ncon = -1;
    graph.mincut = -1;
    graph.minvol = -1;
    graph.nbnd = -1;

    /* memory for the graph structure */
    graph.xadj = std::ptr::null_mut();
    graph.vwgt = std::ptr::null_mut();
    graph.vsize = std::ptr::null_mut();
    graph.adjncy = std::ptr::null_mut();
    graph.adjwgt = std::ptr::null_mut();
    graph.label = std::ptr::null_mut();
    graph.cmap = std::ptr::null_mut();
    graph.tvwgt = std::ptr::null_mut();
    graph.invtvwgt = std::ptr::null_mut();

    /* by default these are set to true, but the can be explicitly changed afterwards */
    graph.free_xadj = 1;
    graph.free_vwgt = 1;
    graph.free_vsize = 1;
    graph.free_adjncy = 1;
    graph.free_adjwgt = 1;

    /* memory for the partition/refinement structure */
    graph.where_ = std::ptr::null_mut();
    graph.pwgts = std::ptr::null_mut();
    graph.id = std::ptr::null_mut();
    graph.ed = std::ptr::null_mut();
    graph.bndptr = std::ptr::null_mut();
    graph.bndind = std::ptr::null_mut();
    graph.nrinfo = std::ptr::null_mut();
    graph.ckrinfo = std::ptr::null_mut();
    graph.vkrinfo = std::ptr::null_mut();

    /* linked-list structure */
    graph.coarser = std::ptr::null_mut();
    graph.finer = std::ptr::null_mut();
}

/*************************************************************************/
/* This function frees the memory storing the structure of the graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeSData(graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    /* free graph structure */
    if graph.free_xadj != 0 {
        free_field!(graph.xadj);
    }
    if graph.free_vwgt != 0 {
        free_field!(graph.vwgt);
    }
    if graph.free_vsize != 0 {
        free_field!(graph.vsize);
    }
    if graph.free_adjncy != 0 {
        free_field!(graph.adjncy);
    }
    if graph.free_adjwgt != 0 {
        free_field!(graph.adjwgt);
    }
}

/*************************************************************************/
/* This function frees the refinement/partition memory stored in a graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeRData(graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    /* The following is for the -minconn and -contig to work properly in
    the vol-refinement routines */
    // if ((void *)graph.ckrinfo == (void *)graph.vkrinfo)
    if std::ptr::eq(graph.ckrinfo as *const (), graph.vkrinfo as *const ()) {
        graph.ckrinfo = std::ptr::null_mut();
    }

    /* free partition/refinement structure */
    free_field!(graph.where_);
    free_field!(graph.pwgts);
    free_field!(graph.id);
    free_field!(graph.ed);
    free_field!(graph.bndptr);
    free_field!(graph.bndind);
    free_field!(graph.nrinfo);
    free_field!(graph.ckrinfo);
    free_field!(graph.vkrinfo);
    // gk_free_one(&graph.where_, &graph.pwgts, &graph.id, &graph.ed,
    //     &graph.bndptr, &graph.bndind, &graph.nrinfo, &graph.ckrinfo,
    //     &graph.vkrinfo);
}

/*************************************************************************/
/* This function deallocates any memory stored in a graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeGraph(r_graph: *mut *mut graph_t) {
    let graph = *r_graph;

    if (*graph).ondisk != 0 {
        panic!("I don't thing graphs on disk should be freed");
    }

    /* free the graph structure's fields */
    FreeSData(graph);

    /* free the partition/refinement fields */
    FreeRData(graph);

    gk_free_one(&mut ((*graph).tvwgt as *mut c_void));
    gk_free_one(&mut ((*graph).invtvwgt as *mut c_void));
    gk_free_one(&mut ((*graph).label as *mut c_void));
    gk_free_one(&mut ((*graph).cmap as *mut c_void));
    gk_free_one(&mut (graph as *mut c_void));

    *r_graph = std::ptr::null_mut();
}

/*************************************************************************/
/* This function writes the key contents of the graph on disk and frees
the associated memory */
/*************************************************************************/
#[metis_func]
pub extern "C" fn graph_WriteToDisk(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t nvtxs, ncon, *xadj;
    // static int gID = 1;
    // char outfile[1024];
    // FILE *fpout;

    static GID: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(1);

    if ctrl.ondisk == 0 {
        // todo!("write tests to use disk and also get rid of sleep");
        return;
    }

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(graph => xadj);

    if std::mem::size_of::<idx_t>() * (nvtxs * (ncon + 1) + 2 * xadj[nvtxs] as usize)
        < 128 * 1024 * 1024
    {
        return;
    }

    if graph.gID > 0 {
        // sprintln!(outfile, "metis%d.%d", ctrl.pid, graph.gID);
        // gk_rmpath(outfile);
        let outfile = format!("metis{}.{}", ctrl.pid, graph.gID);
        let _ = std::fs::remove_file(outfile);
    }

    // graph.gID = gID += 1;
    graph.gID = GID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    // sprintln!(outfile, "metis%d.%d", ctrl.pid, graph.gID);
    let outfile = format!("metis{}.{}", ctrl.pid, graph.gID);
    let _ = std::fs::remove_file(&outfile);

    let ret: std::io::Result<()> = (|| {
        let fpout = std::fs::OpenOptions::new().write(true).create_new(true).open(&outfile)?;
        let mut fpout = std::io::BufWriter::new(fpout);
        // if ((fpout = fopen(outfile, "wb")) == std::ptr::null_mut()) {
        //     return;
        // }

        get_graph_slices!(graph => xadj vwgt adjncy adjwgt);
        if graph.free_xadj != 0 {
            fpout.write_idx(xadj)?;
            // if (fwrite(graph.xadj, sizeof(idx_t), nvtxs + 1, fpout) != nvtxs + 1) {
            //     todo!("goto error\0".as_ptr());
            // }
        }
        if graph.free_vwgt != 0 {
            fpout.write_idx(vwgt)?;

            // if (fwrite(graph.vwgt, sizeof(idx_t), nvtxs * ncon, fpout) != nvtxs * ncon) {
            //     todo!("goto error\0".as_ptr());
            // }
        }
        if graph.free_adjncy != 0 {
            fpout.write_idx(adjncy)?;
            // if (fwrite(graph.adjncy, sizeof(idx_t), xadj[nvtxs], fpout) != xadj[nvtxs]) {
            //     todo!("goto error\0".as_ptr());
            // }
        }
        if graph.free_adjwgt != 0 {
            fpout.write_idx(adjwgt)?;
            // if (fwrite(graph.adjwgt, sizeof(idx_t), xadj[nvtxs], fpout) != xadj[nvtxs]) {
            //     todo!("goto error\0".as_ptr());
            // }
        }
        if ctrl.objtype == METIS_OBJTYPE_VOL {
            get_graph_slices!(graph => vsize);
            if graph.free_vsize != 0 {
                fpout.write_idx(vsize)?;
                // if (fwrite(graph.vsize, sizeof(idx_t), nvtxs, fpout) != nvtxs) {
                //     todo!("goto error\0".as_ptr());
                // }
            }
        }

        // fclose(fpout);
        Ok(())
    })();

    match ret {
        Ok(_) => (),
        Err(e) => {
            let _ = std::fs::remove_file(&outfile);
            panic!("failed writing to {outfile} with error {e}")
        }
    };

    if graph.free_xadj != 0 {
        free_field!(graph.xadj);
    }
    if graph.free_vwgt != 0 {
        free_field!(graph.vwgt);
    }
    if graph.free_vsize != 0 {
        free_field!(graph.vsize);
    }
    if graph.free_adjncy != 0 {
        free_field!(graph.adjncy);
    }
    if graph.free_adjwgt != 0 {
        free_field!(graph.adjwgt);
    }

    graph.ondisk = 1;
    return;

    // error:
    // println!("Failed on writing %s\n", outfile);
    // fclose(fpout);
    // gk_rmpath(outfile);
    // graph.ondisk = 0;
}

/*************************************************************************/
/* This function reads the key contents of a graph from the disk */
/*************************************************************************/
#[metis_func]
pub extern "C" fn graph_ReadFromDisk(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t nvtxs, ncon, *xadj;
    // char infile[1024];
    // FILE *fpin;

    if graph.ondisk == 0 {
        return; /* this graph is not on the disk */
    }
    let infile = format!("metis{}.{}", ctrl.pid, graph.gID);
    let ret: std::io::Result<()> = (|| {
        //   if ((fpin = fopen(infile, "rb")) == std::ptr::null_mut())
        // {
        //   return ;
        // }
        let fpin = std::fs::OpenOptions::new().read(true).open(&infile)?;
        let mut fpin = std::io::BufReader::new(fpin);

        let nvtxs = graph.nvtxs as usize;
        let ncon = graph.ncon as usize;

        if graph.free_xadj != 0 {
            graph.xadj = imalloc(nvtxs + 1, "graph_ReadFromDisk: xadj\0".as_ptr()) as *mut idx_t;
            get_graph_slices_mut!(graph => xadj);
            fpin.read_idx(xadj)?;

            // if (fread(graph.xadj, sizeof(idx_t), nvtxs+1, fpin) != nvtxs+1)
            // {
            //   todo!("goto error");
            // }
        }
        get_graph_slices_mut!(graph => xadj);

        if graph.free_vwgt != 0 {
            graph.vwgt = imalloc(nvtxs * ncon, "graph_ReadFromDisk: vwgt\0".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->vwgt, nvtxs * ncon);
            fpin.read_idx(vwgt)?;

            // if (fread(graph.vwgt, sizeof(idx_t), nvtxs * ncon, fpin) != nvtxs * ncon) {
            //     todo!("goto error");
            // }
        }

        if graph.free_adjncy != 0 {
            let len = xadj[nvtxs] as usize;
            graph.adjncy = imalloc(len, "graph_ReadFromDisk: adjncy\0".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->adjncy, len);
            fpin.read_idx(adjncy)?;
            // if (fread(graph.adjncy, sizeof(idx_t), xadj[nvtxs], fpin) != xadj[nvtxs]) {
            //     todo!("goto error");
            // }
        }

        if graph.free_adjwgt != 0 {
            let len = xadj[nvtxs] as usize;
            graph.adjwgt = imalloc(len, "graph_ReadFromDisk: adjwgt\0".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->adjwgt, len);
            fpin.read_idx(adjwgt)?;

            // if (fread(graph.adjwgt, sizeof(idx_t), xadj[nvtxs], fpin) != xadj[nvtxs]) {
            //     todo!("goto error");
            // }
        }

        if ctrl.objtype == METIS_OBJTYPE_VOL {
            if graph.free_vsize != 0 {
                graph.vsize = imalloc(nvtxs, "graph_ReadFromDisk: vsize\0".as_ptr()) as *mut idx_t;
                mkslice_mut!(graph->vsize, nvtxs);
                fpin.read_idx(vsize)?;
                // if (fread(graph.vsize, sizeof(idx_t), nvtxs, fpin) != nvtxs) {
                //     todo!("goto error");
                // }
            }
        }

        // fclose(fpin);
        //  println!("ondisk: deleting %s\n", infile);
        // gk_rmpath(infile);
        std::mem::drop(fpin);
        graph.gID = 0;
        graph.ondisk = 0;
        return Ok(());
    })();

    let _ = std::fs::remove_file(&infile);

    match ret {
        Ok(_) => return,
        Err(e) => {
            graph.ondisk = 0;
            panic!("Failed to restore graph {infile} from the disk with error {e}")
        }
    };

    // error:
    // fclose(fpin);
    // gk_rmpath(infile);
    // graph.ondisk = 0;
    // gk_errexit(
    //     SIGERR,
    //     "Failed to restore graph %s from the disk.\n",
    //     infile,
    // );
}
