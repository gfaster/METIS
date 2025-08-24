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
    #[expect(dead_code)]
    fn read_real(&mut self, buf: &mut [real_t]) -> std::io::Result<()>;
}

impl<T: Read> ReadTransmuted for T {
    fn read_idx(&mut self, buf: &mut [idx_t]) -> std::io::Result<()> {
        unsafe {
            let p = buf.as_mut_ptr() as *mut u8;
            let buf = std::slice::from_raw_parts_mut(p, std::mem::size_of_val(buf));
            self.read_exact(buf)
        }
    }

    fn read_real(&mut self, buf: &mut [real_t]) -> std::io::Result<()> {
        unsafe {
            let p = buf.as_mut_ptr() as *mut u8;
            let buf = std::slice::from_raw_parts_mut(p, std::mem::size_of_val(buf));
            self.read_exact(buf)
        }
    }
}

trait WriteTransmuted {
    fn write_idx(&mut self, buf: &[idx_t]) -> std::io::Result<()>;
    #[expect(dead_code)]
    fn write_real(&mut self, buf: &[real_t]) -> std::io::Result<()>;
}

impl<T: std::io::Write> WriteTransmuted for T {
    fn write_idx(&mut self, buf: &[idx_t]) -> std::io::Result<()> {
        unsafe {
            let p = buf.as_ptr() as *const u8;
            let buf = std::slice::from_raw_parts(p, std::mem::size_of_val(buf));
            self.write_all(buf)
        }
    }

    fn write_real(&mut self, buf: &[real_t]) -> std::io::Result<()> {
        unsafe {
            let p = buf.as_ptr() as *const u8;
            let buf = std::slice::from_raw_parts(p, std::mem::size_of_val(buf));
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
            graph.vwgt = vwgt;
            graph.free_vwgt = 0;
        } else {
            vwgt = ismalloc(ncon * nvtxs, 1, c"SetupGraph: vwgt".as_ptr()) as *mut idx_t;
            graph.vwgt = vwgt;
        }

        graph.tvwgt = imalloc(ncon, c"SetupGraph: tvwgts".as_ptr()) as *mut idx_t;
        graph.invtvwgt = rmalloc(ncon, c"SetupGraph: invtvwgts".as_ptr()) as *mut real_t;
        for i in (0)..(ncon) {
            mkslice_mut!(vwgt, ncon * nvtxs);
            // (*graph).tvwgt[i]    = isum(nvtxs, vwgt+i, ncon);
            *graph.tvwgt.add(i) = vwgt[i..].iter().step_by(ncon).sum::<idx_t>();
            *graph.invtvwgt.add(i) = (1.0_f64
                / (if *graph.tvwgt.add(i) > 0 {
                    *graph.tvwgt.add(i) as f64
                } else {
                    1.0_f64
                })) as real_t;
        }

        if ctrl.objtype == METIS_OBJTYPE_VOL {
            /* Setup the vsize */
            if !vsize.is_null() {
                graph.vsize = vsize;
                graph.free_vsize = 0;
            } else {
                vsize = ismalloc(nvtxs, 1, c"SetupGraph: vsize".as_ptr()) as *mut idx_t;
                graph.vsize = vsize;
            }

            /* Allocate memory for edge weights and initialize them to the sum of the vsize */
            adjwgt = imalloc(graph.nedges as usize, c"SetupGraph: adjwgt".as_ptr()) as *mut idx_t;
            graph.adjwgt = adjwgt;
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
                adjwgt = ismalloc(graph.nedges as usize, 1, c"SetupGraph: adjwgt".as_ptr())
                    as *mut idx_t;
                graph.adjwgt = adjwgt;
            }
        }

        /* setup various derived info */
        SetupGraph_tvwgt(graph);

        if ctrl.optype == METIS_OP_PMETIS || ctrl.optype == METIS_OP_OMETIS {
            SetupGraph_label(graph);
        }

        assert!(checkgraph::CheckGraph(graph, ctrl.numflag, 1) != 0);
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
    if graph.tvwgt.is_null() {
        graph.tvwgt = imalloc(ncon, c"SetupGraph_tvwgt: tvwgt".as_ptr()) as *mut idx_t;
    }
    if graph.invtvwgt.is_null() {
        graph.invtvwgt = rmalloc(ncon, c"SetupGraph_tvwgt: invtvwgt".as_ptr()) as *mut real_t;
    }
    get_graph_slices!(graph => vwgt);
    get_graph_slices_mut!(graph => tvwgt invtvwgt);

    for i in (0)..(ncon) {
        // tvwgt[i] = isum(graph.nvtxs, graph.vwgt + i, graph.ncon);
        // invtvwgt[i] = 1.0 / (if graph.tvwgt[i] > 0 { graph.tvwgt[i] } else { 1 });
        // TODO: I don't think this gets specialized well
        // Using skip rather than `i..' to handle the nvtxs = 0 case
        tvwgt[i] = vwgt.iter().skip(i).step_by(ncon).sum();
        invtvwgt[i] = (1.0
            / (if tvwgt[i] > 0 {
                tvwgt[i] as double
            } else {
                1.0
            })) as real_t;
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
    if graph.label.is_null() {
        graph.label = imalloc(nvtxs, c"SetupGraph_label: label".as_ptr()) as *mut idx_t;
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
        sgraph.xadj = imalloc(snvtxs + 1, c"SetupSplitGraph: xadj".as_ptr()) as *mut idx_t;
        sgraph.vwgt = imalloc(sncon * snvtxs, c"SetupSplitGraph: vwgt".as_ptr()) as *mut idx_t;
        sgraph.adjncy = imalloc(snedges, c"SetupSplitGraph: adjncy".as_ptr()) as *mut idx_t;
        sgraph.adjwgt = imalloc(snedges, c"SetupSplitGraph: adjwgt".as_ptr()) as *mut idx_t;
        sgraph.label = imalloc(snvtxs, c"SetupSplitGraph: label".as_ptr()) as *mut idx_t;
        sgraph.tvwgt = imalloc(sncon, c"SetupSplitGraph: tvwgt".as_ptr()) as *mut idx_t;
        sgraph.invtvwgt = rmalloc(sncon, c"SetupSplitGraph: invtvwgt".as_ptr()) as *mut real_t;

        if !graph.vsize.is_null() {
            sgraph.vsize = imalloc(snvtxs, c"SetupSplitGraph: vsize".as_ptr()) as *mut idx_t;
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
        c"CreateGraph: graph".as_ptr(),
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
        let fpout = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&outfile)?;
        let mut fpout = std::io::BufWriter::new(fpout);

        get_graph_slices!(graph => xadj vwgt adjncy adjwgt);
        if graph.free_xadj != 0 {
            fpout.write_idx(xadj)?;
        }
        if graph.free_vwgt != 0 {
            fpout.write_idx(vwgt)?;
        }
        if graph.free_adjncy != 0 {
            fpout.write_idx(adjncy)?;
        }
        if graph.free_adjwgt != 0 {
            fpout.write_idx(adjwgt)?;
        }
        if ctrl.objtype == METIS_OBJTYPE_VOL
            && graph.free_vsize != 0 {
                get_graph_slices!(graph => vsize);
                fpout.write_idx(vsize)?;
            }
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
}

/*************************************************************************/
/* This function reads the key contents of a graph from the disk */
/*************************************************************************/
#[metis_func]
pub extern "C" fn graph_ReadFromDisk(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    if graph.ondisk == 0 {
        return; /* this graph is not on the disk */
    }
    let infile = format!("metis{}.{}", ctrl.pid, graph.gID);
    let ret: std::io::Result<()> = (|| {
        let fpin = std::fs::OpenOptions::new().read(true).open(&infile)?;
        let mut fpin = std::io::BufReader::new(fpin);

        let nvtxs = graph.nvtxs as usize;
        let ncon = graph.ncon as usize;

        if graph.free_xadj != 0 {
            graph.xadj = imalloc(nvtxs + 1, c"graph_ReadFromDisk: xadj".as_ptr()) as *mut idx_t;
            get_graph_slices_mut!(graph => xadj);
            fpin.read_idx(xadj)?;
        }
        get_graph_slices_mut!(graph => xadj);

        if graph.free_vwgt != 0 {
            graph.vwgt = imalloc(nvtxs * ncon, c"graph_ReadFromDisk: vwgt".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->vwgt, nvtxs * ncon);
            fpin.read_idx(vwgt)?;
        }

        if graph.free_adjncy != 0 {
            let len = xadj[nvtxs] as usize;
            graph.adjncy = imalloc(len, c"graph_ReadFromDisk: adjncy".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->adjncy, len);
            fpin.read_idx(adjncy)?;
        }

        if graph.free_adjwgt != 0 {
            let len = xadj[nvtxs] as usize;
            graph.adjwgt = imalloc(len, c"graph_ReadFromDisk: adjwgt".as_ptr()) as *mut idx_t;
            mkslice_mut!(graph->adjwgt, len);
            fpin.read_idx(adjwgt)?;
        }

        if ctrl.objtype == METIS_OBJTYPE_VOL
            && graph.free_vsize != 0 {
                graph.vsize = imalloc(nvtxs, c"graph_ReadFromDisk: vsize".as_ptr()) as *mut idx_t;
                mkslice_mut!(graph->vsize, nvtxs);
                fpin.read_idx(vsize)?;
            }
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
}
