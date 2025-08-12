/*
\file
\brief Functions dealing with memory allocation and workspace management

\date Started 2/24/96
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: wspace.c 10492 2011-07-06 09:28:42Z karypis $
*/

use crate::*;

/*************************************************************************/
/* This function allocates memory for the workspace */
/*************************************************************************/
#[metis_func]
pub extern "C" fn AllocateWorkSpace(ctrl: *mut ctrl_t, graph: *const graph_t) {
    let graph = graph.as_ref().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    let nvtxs = graph.nvtxs as usize;
    let nparts = ctrl.nparts as usize;
    let ncon = ctrl.ncon as usize;

    let coresize: usize;
    match ctrl.optype {
        METIS_OP_PMETIS => {
            coresize = 3 * (nvtxs + 1) * size_of::<idx_t>()
                + 5 * (nparts + 1) * ncon * size_of::<idx_t>()
                + 5 * (nparts + 1) * ncon * size_of::<real_t>()
        }
        _ => {
            coresize = 4 * (nvtxs + 1) * size_of::<idx_t>()
                + 5 * (nparts + 1) * ncon * size_of::<idx_t>()
                + 5 * (nparts + 1) * ncon * size_of::<real_t>()
        }
    }
    ctrl.mcore = gk_mcoreCreate(coresize);

    ctrl.nbrpoolsize = 0;
    ctrl.nbrpoolcpos = 0;
}

/*************************************************************************/
/* This function allocates refinement-specific memory for the workspace */
/*************************************************************************/
#[metis_func]
pub extern "C" fn AllocateRefinementWorkSpace(
    ctrl: *mut ctrl_t,
    nbrpoolsize_max: idx_t,
    nbrpoolsize: idx_t,
) {
    let ctrl = ctrl.as_mut().unwrap();

    ctrl.nbrpoolsize_max = nbrpoolsize_max as usize;
    ctrl.nbrpoolsize = nbrpoolsize as usize;
    ctrl.nbrpoolcpos = 0;
    ctrl.nbrpoolreallocs = 0;

    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            ctrl.cnbrpool = gk_malloc(
                ctrl.nbrpoolsize as usize * size_of::<cnbr_t>(),
                c"AllocateRefinementWorkSpace: cnbrpool".as_ptr(),
            )
            .cast()
        }

        METIS_OBJTYPE_VOL => {
            ctrl.vnbrpool = gk_malloc(
                ctrl.nbrpoolsize as usize * size_of::<vnbr_t>(),
                c"AllocateRefinementWorkSpace: vnbrpool".as_ptr(),
            )
            .cast()
        }

        _ => panic!("Unknown objtype of {}", ctrl.objtype),
    }

    /* Allocate the memory for the sparse subdomain graph */
    if ctrl.minconn != 0 {
        let nparts = ctrl.nparts as usize;
        ctrl.pvec1 = imalloc(nparts + 1, c"AllocateRefinementWorkSpace: pvec1".as_ptr()).cast();
        ctrl.pvec2 = imalloc(nparts + 1, c"AllocateRefinementWorkSpace: pvec2".as_ptr()).cast();
        ctrl.maxnads = ismalloc(
            nparts,
            INIT_MAXNAD,
            c"AllocateRefinementWorkSpace: maxnads".as_ptr(),
        )
        .cast();
        ctrl.nads = imalloc(nparts, c"AllocateRefinementWorkSpace: nads".as_ptr()).cast();
        ctrl.adids = iAllocMatrix(
            nparts,
            INIT_MAXNAD as usize,
            0,
            c"AllocateRefinementWorkSpace: adids".as_ptr(),
        )
        .cast();
        ctrl.adwgts = iAllocMatrix(
            nparts,
            INIT_MAXNAD as usize,
            0,
            c"AllocateRefinementWorkSpace: adwgts".as_ptr(),
        )
        .cast();
    }
}

/*************************************************************************/
/* This function frees the workspace */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeWorkSpace(ctrl: *mut ctrl_t) {
    let ctrl = ctrl.as_mut().unwrap();

    gk_mcoreDestroy(
        &mut ctrl.mcore,
        (ctrl.dbglvl & METIS_DBG_INFO) as std::ffi::c_int,
    );

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_INFO,
        print!(
            concat!(
                " nbrpool statistics\n",
                "        nbrpoolsize: {:12}   nbrpoolcpos: {:12}\n",
                "    nbrpoolreallocs: {:12}\n\n"
            ),
            ctrl.nbrpoolsize, ctrl.nbrpoolcpos, ctrl.nbrpoolreallocs
        )
    );

    if ctrl.cnbrpool as usize == ctrl.vnbrpool as usize {
        ctrl.vnbrpool = std::ptr::null_mut();
    }
    gk::free_ref(&mut ctrl.cnbrpool);
    gk::free_ref(&mut ctrl.vnbrpool);
    ctrl.nbrpoolsize_max = 0;
    ctrl.nbrpoolsize = 0;
    ctrl.nbrpoolcpos = 0;

    if ctrl.minconn != 0 {
        iFreeMatrix(&mut ctrl.adids, ctrl.nparts as usize, INIT_MAXNAD as usize);
        iFreeMatrix(&mut ctrl.adwgts, ctrl.nparts as usize, INIT_MAXNAD as usize);

        gk::free_ref(&mut ctrl.pvec1);
        gk::free_ref(&mut ctrl.pvec2);
        gk::free_ref(&mut ctrl.maxnads);
        gk::free_ref(&mut ctrl.nads);
        // gk_free((void **)&ctrl.pvec1, &ctrl.pvec2,
        //     &ctrl.maxnads, &ctrl.nads, LTERM);
    }
}

/*************************************************************************/
/* This function resets the cnbrpool */
/*************************************************************************/
#[metis_func]
pub extern "C" fn cnbrpoolReset(ctrl: &mut ctrl_t) {
    ctrl.nbrpoolcpos = 0;
}

/*************************************************************************/
/* This function gets the next free index from cnbrpool */
/*************************************************************************/
#[metis_func]
pub extern "C" fn cnbrpoolGetNext(ctrl: *mut ctrl_t, nnbrs: idx_t) -> idx_t {
    let ctrl = ctrl.as_mut().unwrap();
    /* add 1 because when moving vertices, an extra neighbor can be temporarily
     * needed (particularly when minconn is set) */
    let nnbrs = (ctrl.nparts).min(nnbrs) as usize + 1;
    ctrl.nbrpoolcpos += nnbrs as usize;

    if ctrl.nbrpoolcpos > ctrl.nbrpoolsize {
        ctrl.nbrpoolsize += (10 * nnbrs).max(ctrl.nbrpoolsize / 2);
        ctrl.nbrpoolsize = (ctrl.nbrpoolsize).min(ctrl.nbrpoolsize_max);

        ctrl.cnbrpool = gk_realloc(
            ctrl.cnbrpool.cast(),
            ctrl.nbrpoolsize * size_of::<cnbr_t>(),
            c"cnbrpoolGet: cnbrpool".as_ptr(),
        )
        .cast();
        ctrl.nbrpoolreallocs += 1;

        debug_assert!(ctrl.nbrpoolcpos <= ctrl.nbrpoolsize);
    }

    return (ctrl.nbrpoolcpos - nnbrs) as idx_t;
}

/*************************************************************************/
/* This function resets the vnbrpool */
/*************************************************************************/
#[metis_func]
pub extern "C" fn vnbrpoolReset(ctrl: &mut ctrl_t) {
    ctrl.nbrpoolcpos = 0;
}

/*************************************************************************/
/* This function gets the next free index from vnbrpool */
/*************************************************************************/
#[metis_func]
pub extern "C" fn vnbrpoolGetNext(ctrl: *mut ctrl_t, nnbrs: idx_t) -> idx_t {
    let ctrl = ctrl.as_mut().unwrap();

    /* add 1 because when moving vertices, an extra neighbor can be temporarily
     * needed (particularly when minconn is set) */
    let nnbrs = (ctrl.nparts).min(nnbrs) as usize + 1;
    ctrl.nbrpoolcpos += nnbrs;

    if ctrl.nbrpoolcpos > ctrl.nbrpoolsize {
        ctrl.nbrpoolsize += (10 * nnbrs).max(ctrl.nbrpoolsize / 2);
        ctrl.nbrpoolsize = (ctrl.nbrpoolsize).min(ctrl.nbrpoolsize_max);

        ctrl.vnbrpool = gk_realloc(
            ctrl.vnbrpool.cast(),
            ctrl.nbrpoolsize * size_of::<vnbr_t>(),
            c"vnbrpoolGet: vnbrpool".as_ptr(),
        )
        .cast();
        ctrl.nbrpoolreallocs += 1;

        debug_assert!(ctrl.nbrpoolcpos <= ctrl.nbrpoolsize);
    }

    return (ctrl.nbrpoolcpos - nnbrs) as idx_t;
}
