/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * meshpart.c
 *
 * This file contains routines for partitioning finite element meshes.
 *
 * Started 9/29/97
 * George
 *
 * $Id: meshpart.c 17513 2014-08-05 16:20:50Z dominique $
 *
 */

use crate::*;

/// stolen from options.rs
macro_rules! get_option {
    ($options:expr, $opt:expr, $default:expr) => {
        if $options.is_null() || *$options.add($opt as usize) == -1 {
            $default
        } else {
            *$options.add($opt as usize) as _
        }
    };
}

/*************************************************************************
* This function partitions a finite element mesh by partitioning its nodal
* graph using KMETIS and then assigning elements in a load balanced fashion.
**************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_PartMeshNodal(
    ne: *mut idx_t,
    nn: *mut idx_t,
    eptr: *mut idx_t,
    eind: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    nparts: *mut idx_t,
    tpwgts: *mut real_t,
    options: *mut idx_t,
    objval: *mut idx_t,
    epart: *mut idx_t,
    npart: *mut idx_t,
) -> std::ffi::c_int {
    // int sigrval=0, renumber=0, ptype;
    // idx_t *xadj=std::ptr::null_mut(), *adjncy=std::ptr::null_mut();
    // idx_t ncon=1, pnumflag=0;
    // int rstatus=METIS_OK;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        return METIS_ERROR_MEMORY;
    }

    //   gk_sigtrap();
    //
    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    let _renumber = get_option!(options, METIS_OPTION_NUMBERING, 0);
    let ptype = get_option!(options, METIS_OPTION_PTYPE, METIS_PTYPE_KWAY);

    /* renumber the mesh */
    // if (renumber) {
    //   ChangeMesh2CNumbering(*ne, eptr, eind);
    //   options[METIS_OPTION_NUMBERING as usize] = 0;
    // }

    /* get the nodal graph */
    let mut xadj = std::ptr::null_mut();
    let mut adjncy = std::ptr::null_mut();
    let rstatus = mesh::METIS_MeshToNodal(ne, nn, eptr, eind, &0, &mut xadj, &mut adjncy);
    assert!(rstatus == METIS_OK);
    //   if (rstatus != METIS_OK) {
    //     libc::raise(SIGERR);
    // }

    /* partition the graph */
    let rstatus = if ptype == METIS_PTYPE_KWAY {
        kmetis::METIS_PartGraphKway(
            nn,
            &1,
            xadj,
            adjncy,
            vwgt,
            vsize,
            std::ptr::null_mut(),
            nparts,
            tpwgts,
            std::ptr::null_mut(),
            options,
            objval,
            npart,
        )
    } else {
        pmetis::METIS_PartGraphRecursive(
            nn,
            &1,
            xadj,
            adjncy,
            vwgt,
            vsize,
            std::ptr::null_mut(),
            nparts,
            tpwgts,
            std::ptr::null_mut(),
            options,
            objval,
            npart,
        )
    };

    assert!(rstatus == METIS_OK);
    //   if (rstatus != METIS_OK) {
    //     raise(SIGERR);
    // }

    /* partition the other side of the mesh */
    InduceRowPartFromColumnPart(*nn, *ne, eptr, eind, epart, npart, *nparts, tpwgts);

    // SIGTHROW:
    // if (renumber) {
    //   ChangeMesh2FNumbering2(*ne, *nn, eptr, eind, epart, npart);
    //   options[METIS_OPTION_NUMBERING as usize] = 1;
    // }

    METIS_Free(xadj.cast());
    METIS_Free(adjncy.cast());

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    return METIS_OK;
}

/*************************************************************************
* This function partitions a finite element mesh by partitioning its dual
* graph using KMETIS and then assigning nodes in a load balanced fashion.
**************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_PartMeshDual(
    ne: *const idx_t,
    nn: *const idx_t,
    eptr: *mut idx_t,
    eind: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    ncommon: *mut idx_t,
    nparts: *mut idx_t,
    tpwgts: *mut real_t,
    options: *mut idx_t,
    objval: *mut idx_t,
    epart: *mut idx_t,
    npart: *mut idx_t,
) -> std::ffi::c_int {
    // int sigrval=0, renumber=0, ptype;
    // idx_t i, j;
    // idx_t *xadj=std::ptr::null_mut(), *adjncy=std::ptr::null_mut(), *nptr=std::ptr::null_mut(), *nind=std::ptr::null_mut();
    // idx_t ncon=1, pnumflag=0;
    // int rstatus = METIS_OK;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        return METIS_ERROR_MEMORY;
    }

    let ne = *ne as usize;
    let nn = *nn as usize;

    //   gk_sigtrap();
    //
    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    let _renumber = get_option!(options, METIS_OPTION_NUMBERING, 0);
    let ptype = get_option!(options, METIS_OPTION_PTYPE, METIS_PTYPE_KWAY);

    /* renumber the mesh */
    // if (renumber) {
    //   ChangeMesh2CNumbering(*ne, eptr, eind);
    //   options[METIS_OPTION_NUMBERING as usize] = 0;
    // }

    /* get the dual graph */
    let mut xadj = std::ptr::null_mut();
    let mut adjncy = std::ptr::null_mut();
    let rstatus = mesh::METIS_MeshToDual(
        &(ne as idx_t),
        &(nn as idx_t),
        eptr,
        eind,
        ncommon,
        &0,
        &mut xadj,
        &mut adjncy,
    );
    assert!(rstatus == METIS_OK);
    //   if (rstatus != METIS_OK) {
    //     raise(SIGERR);
    // }

    /* partition the graph */
    let rstatus = if ptype == METIS_PTYPE_KWAY {
        kmetis::METIS_PartGraphKway(
            &(ne as idx_t),
            &1,
            xadj,
            adjncy,
            vwgt,
            vsize,
            std::ptr::null_mut(),
            nparts,
            tpwgts,
            std::ptr::null_mut(),
            options,
            objval,
            epart,
        )
    } else {
        pmetis::METIS_PartGraphRecursive(
            &(ne as idx_t),
            &1,
            xadj,
            adjncy,
            vwgt,
            vsize,
            std::ptr::null_mut(),
            nparts,
            tpwgts,
            std::ptr::null_mut(),
            options,
            objval,
            epart,
        )
    };
    assert!(rstatus == METIS_OK);
    //   if (rstatus != METIS_OK) {
    //     raise(SIGERR);
    // }

    /* construct the node-element list */
    // nptr = ismalloc(*nn+1, 0, "METIS_PartMeshDual: nptr");
    // nind = imalloc(eptr[(*ne) as usize], "METIS_PartMeshDual: nind");
    mkslice!(eptr, ne + 1);
    mkslice!(eind, eptr[ne]);

    let nptr: &mut [idx_t] = &mut vec![0; nn + 1];
    let nind: &mut [idx_t] = &mut vec![0; eptr[ne] as usize];

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nptr[eind[j as usize] as usize] += 1;
        }
    }
    util::make_csr(nn, nptr);

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nind[nptr[eind[j as usize] as usize] as usize] = i as idx_t;
            nptr[eind[j as usize] as usize] += 1
        }
    }
    util::shift_csr(nn, nptr);

    assert_eq!(nptr[nn], eptr[ne], "idk if this is right");

    /* partition the other side of the mesh */
    InduceRowPartFromColumnPart(
        ne as idx_t,
        nn as idx_t,
        nptr.as_ptr(),
        nind.as_ptr(),
        npart,
        epart,
        *nparts,
        tpwgts,
    );

    // gk_free((void **)&nptr, &nind, LTERM);

    // SIGTHROW:
    //   if (renumber) {
    //     ChangeMesh2FNumbering2(*ne, *nn, eptr, eind, epart, npart);
    //     options[METIS_OPTION_NUMBERING as usize] = 1;
    //   }

    METIS_Free(xadj.cast());
    METIS_Free(adjncy.cast());

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    return METIS_OK;
}

/*************************************************************************/
/* Induces a partitioning of the rows based on a a partitioning of the
columns. It is used by both the Nodal and Dual routines. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn InduceRowPartFromColumnPart(
    ncols: idx_t,
    nrows: idx_t,
    rowptr: *const idx_t,
    rowind: *const idx_t,
    rpart: *mut idx_t,
    cpart: *mut idx_t,
    nparts: idx_t,
    tpwgts: *const real_t,
) {
    // idx_t i, j, k, me;
    // idx_t nnbrs, *pwgts, *nbrdom, *nbrwgt, *nbrmrk;
    // idx_t *itpwgts;

    let nparts = nparts as usize;
    let nrows = nrows as usize;
    let ncols = ncols as usize;

    mkslice!(rowptr, nrows + 1);
    mkslice!(rowind, rowptr[nrows]); // I'm like 50% sure this is right

    mkslice_mut!(cpart, ncols);
    mkslice_mut!(rpart, nrows);
    rpart.fill(-1);

    // pwgts  = ismalloc(nparts, 0, "InduceRowPartFromColumnPart: pwgts");
    // nbrdom = ismalloc(nparts, 0, "InduceRowPartFromColumnPart: nbrdom");
    // nbrwgt = ismalloc(nparts, 0, "InduceRowPartFromColumnPart: nbrwgt");
    // nbrmrk = ismalloc(nparts, -1, "InduceRowPartFromColumnPart: nbrmrk");
    let pwgts: &mut [idx_t] = &mut vec![0; nparts];
    let nbrdom: &mut [idx_t] = &mut vec![0; nparts];
    let nbrwgt: &mut [idx_t] = &mut vec![0; nparts];
    let nbrmrk: &mut [idx_t] = &mut vec![-1; nparts];

    // iset(nrows, -1, rpart);

    /* setup the integer target partition weights */
    let itpwgts: &mut [idx_t] = &mut vec![0; nparts];
    if tpwgts.is_null() {
        itpwgts.fill(1 + (nrows / nparts) as idx_t);
    } else {
        mkslice!(tpwgts, nparts);
        for (itpwgt, &tpwgt) in itpwgts.iter_mut().zip(tpwgts) {
            //     itpwgts[i as usize] = 1+nrows*tpwgts[i as usize];
            *itpwgt = (1.0 + (nrows as real_t) * tpwgt) as idx_t;
        }
    }

    /* first assign the rows consisting only of columns that belong to
    a single partition. Assign rows that are empty to -2 (un-assigned) */
    for i in (0)..(nrows) {
        if rowptr[(i + 1) as usize] - rowptr[i as usize] == 0 {
            rpart[i as usize] = -2;
            continue;
        }

        let me = cpart[rowind[rowptr[i as usize] as usize] as usize];
        // for j in (rowptr[i as usize] + 1)..(rowptr[(i + 1) as usize]) {
        //     if cpart[rowind[j as usize] as usize] != me {
        //         break;
        //     }
        // }
        // if j == rowptr[(i + 1) as usize] {
        //     rpart[i as usize] = me;
        //     pwgts[me as usize] += 1;
        // }
        if !(rowptr[i] + 1..rowptr[i + 1]).all(|j| cpart[rowind[j as usize] as usize] == me) {
            rpart[i] = me;
            pwgts[me as usize] += 1;
        }
    }

    /* next assign the rows consisting of columns belonging to multiple
    partitions in a  balanced way */
    for i in (0)..(nrows) {
        if rpart[i as usize] == -1 {
            let mut nnbrs = 0;
            for j in (rowptr[i as usize])..(rowptr[(i + 1) as usize]) {
                let me = cpart[rowind[j as usize] as usize];
                if nbrmrk[me as usize] == -1 {
                    nbrdom[nnbrs] = me;
                    nbrwgt[nnbrs] = 1;
                    nbrmrk[me as usize] = nnbrs as idx_t;
                    nnbrs += 1;
                } else {
                    nbrwgt[nbrmrk[me as usize] as usize] += 1;
                }
            }
            debug_assert!(nnbrs > 0);

            /* assign it first to the domain with most things in common */
            rpart[i as usize] = nbrdom[(util::iargmax(&nbrwgt[..nnbrs], 1)) as usize];

            /* if overweight, assign it to the light domain */
            if pwgts[rpart[i as usize] as usize] > itpwgts[rpart[i as usize] as usize] {
                for j in (0)..(nnbrs) {
                    if pwgts[nbrdom[j as usize] as usize] < itpwgts[nbrdom[j as usize] as usize]
                        || pwgts[nbrdom[j as usize] as usize] - itpwgts[nbrdom[j as usize] as usize]
                            < pwgts[rpart[i as usize] as usize]
                                - itpwgts[rpart[i as usize] as usize]
                    {
                        rpart[i as usize] = nbrdom[j as usize];
                        break;
                    }
                }
            }
            pwgts[rpart[i as usize] as usize] += 1;

            /* reset nbrmrk array */
            for j in (0)..(nnbrs) {
                nbrmrk[nbrdom[j as usize] as usize] = -1;
            }
        }
    }

    // gk_free((void **)&pwgts, &nbrdom, &nbrwgt, &nbrmrk, &itpwgts, LTERM);
}
