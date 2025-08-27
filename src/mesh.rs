/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh.c
 *
 * This file contains routines for converting 3D and 4D finite element
 * meshes into dual or nodal graphs
 *
 * Started 8/18/97
 * George
 *
 * $Id: mesh.c 13804 2013-03-04 23:49:08Z karypis $
 *
 */

use crate::*;

/*****************************************************************************/
/* This function creates a graph corresponding to the dual of a finite element
    mesh.

    \param ne is the number of elements in the mesh.
    \param nn is the number of nodes in the mesh.
    \param eptr is an array of size ne+1 used to mark the start and end
           locations in the nind array.
    \param eind is an array that stores for each element the set of node IDs
           (indices) that it is made off. The length of this array is equal
           to the total number of nodes over all the mesh elements.
    \param ncommon is the minimum number of nodes that two elements must share
           in order to be connected via an edge in the dual graph.
    \param numflag is either 0 or 1 indicating if the numbering of the nodes
           starts from 0 or 1, respectively. The same numbering is used for the
           returned graph as well.
    \param r_xadj indicates where_ the adjacency list of each vertex is stored
           in r_adjncy. The memory for this array is allocated by this routine.
           It can be freed by calling METIS_free().
    \param r_adjncy stores the adjacency list of each vertex in the generated
           dual graph. The memory for this array is allocated by this routine.
           It can be freed by calling METIS_free().

*/
/*****************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_MeshToDual(
    ne: *const idx_t,
    nn: *const idx_t,
    eptr: *mut idx_t,
    eind: *mut idx_t,
    ncommon: *const idx_t,
    _numflag: *const idx_t,
    r_xadj: *mut *mut idx_t,
    r_adjncy: *mut *mut idx_t,
) -> std::ffi::c_int {
    // int sigrval=0, renumber=0;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        return METIS_ERROR_MEMORY;
    }

    // gk_sigtrap();

    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    /* renumber the mesh */
    // if (*numflag == 1) {
    //   ChangeMesh2CNumbering(*ne, eptr, eind);
    //   renumber = 1;
    // }

    /* create dual graph */
    *r_xadj = std::ptr::null_mut();
    *r_adjncy = std::ptr::null_mut();
    CreateGraphDual(*ne, *nn, eptr, eind, *ncommon, r_xadj, r_adjncy);

    // SIGTHROW:
    //   if (renumber) {
    //     ChangeMesh2FNumbering(*ne, eptr, eind, *ne, *r_xadj, *r_adjncy);
    // }

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    return METIS_OK;
}

/*****************************************************************************/
/* This function creates a graph corresponding to (almost) the nodal of a
    finite element mesh. In the nodal graph, each node is connected to the
    nodes corresponding to the union of nodes present in all the elements
    in which that node belongs.

    \param ne is the number of elements in the mesh.
    \param nn is the number of nodes in the mesh.
    \param eptr is an array of size ne+1 used to mark the start and end
           locations in the nind array.
    \param eind is an array that stores for each element the set of node IDs
           (indices) that it is made off. The length of this array is equal
           to the total number of nodes over all the mesh elements.
    \param numflag is either 0 or 1 indicating if the numbering of the nodes
           starts from 0 or 1, respectively. The same numbering is used for the
           returned graph as well.
    \param r_xadj indicates where_ the adjacency list of each vertex is stored
           in r_adjncy. The memory for this array is allocated by this routine.
           It can be freed by calling METIS_free().
    \param r_adjncy stores the adjacency list of each vertex in the generated
           dual graph. The memory for this array is allocated by this routine.
           It can be freed by calling METIS_free().

*/
/*****************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_MeshToNodal(
    ne: *const idx_t,
    nn: *const idx_t,
    eptr: *mut idx_t,
    eind: *mut idx_t,
    _numflag: *const idx_t,
    r_xadj: *mut *mut idx_t,
    r_adjncy: *mut *mut idx_t,
) -> std::ffi::c_int {
    // int sigrval=0, renumber=0;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        return METIS_ERROR_MEMORY;
    }

    // gk_sigtrap();

    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    /* renumber the mesh */
    // if (*numflag == 1) {
    //   ChangeMesh2CNumbering(*ne, eptr, eind);
    //   renumber = 1;
    // }

    /* create nodal graph */
    *r_xadj = std::ptr::null_mut();
    *r_adjncy = std::ptr::null_mut();
    CreateGraphNodal(*ne, *nn, eptr, eind, r_xadj, r_adjncy);

    // SIGTHROW:
    //   if (renumber) {
    //     ChangeMesh2FNumbering(*ne, eptr, eind, *nn, *r_xadj, *r_adjncy);
    // }

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    //   if (sigrval != 0) {
    //     if (*r_xadj != std::ptr::null_mut()) {
    //       free(*r_xadj);
    // }
    //     if (*r_adjncy != std::ptr::null_mut()) {
    //       free(*r_adjncy);
    // }
    //     *r_xadj = *r_adjncy = std::ptr::null_mut();
    //   }

    return METIS_OK;
}

/*****************************************************************************/
/* This function creates the dual of a finite element mesh */
/*****************************************************************************/
#[metis_func]
pub extern "C" fn CreateGraphDual(
    ne: idx_t,
    nn: idx_t,
    eptr: *const idx_t,
    eind: *const idx_t,
    ncommon: idx_t,
    r_xadj: *mut *mut idx_t,
    r_adjncy: *mut *mut idx_t,
) {
    let nn = nn as usize;
    let ne = ne as usize;
    // idx_t i, j, nnbrs;
    // idx_t *nptr, *nind;
    // idx_t *xadj, *adjncy;
    // idx_t *marker, *nbrs;

    let ncommon = if ncommon < 1 {
        println!("  Increased ncommon to 1, as it was initially {:}", ncommon);
        1
    } else {
        ncommon
    };

    mkslice!(eptr, ne + 1);
    mkslice!(eind, eptr[ne]);

    /* construct the node-element list first */
    let nptr: &mut [idx_t] = &mut vec![0; nn + 1];
    let nind: &mut [idx_t] = &mut vec![0; eptr[ne] as usize];

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nptr[eind[j as usize] as usize] += 1;
        }
    }
    // MAKECSR(i, nn, nptr);
    util::make_csr(nn, nptr);

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nind[nptr[eind[j as usize] as usize] as usize] = i as idx_t;
            nptr[eind[j as usize] as usize] += 1
        }
    }
    util::shift_csr(nn, nptr);

    /* Allocate memory for xadj, since you know its size.
    These are done using standard malloc as they are returned
    to the calling function */
    let xadj = libc::calloc(ne + 1, size_of::<idx_t>()).cast::<idx_t>();
    if xadj.is_null() {
        panic!("***Failed to allocate memory for xadj.");
    }
    *r_xadj = xadj;
    mkslice_mut!(xadj, ne + 1);

    /* allocate memory for working arrays used by FindCommonElements */
    // marker = ismalloc(ne, 0, "CreateGraphDual: marker");
    let marker: &mut [idx_t] = &mut vec![0; ne];
    // nbrs   = imalloc(ne, "CreateGraphDual: nbrs");
    let nbrs: &mut [idx_t] = &mut vec![0; ne];

    for i in (0)..(ne) {
        xadj[i as usize] = FindCommonElements(
            ne as idx_t,
            nn as idx_t,
            i as idx_t,
            eptr[i + 1] - eptr[i],
            eind[eptr[i as usize] as usize..].as_ptr(),
            nptr.as_ptr(),
            nind.as_ptr(),
            eptr.as_ptr(),
            ncommon,
            marker.as_mut_ptr(),
            nbrs.as_mut_ptr(),
        );
    }
    util::make_csr(ne, xadj);

    /* Allocate memory for adjncy, since you now know its size.
    These are done using standard malloc as they are returned
    to the calling function */
    let adjncy = libc::calloc(xadj[ne as usize] as usize, size_of::<idx_t>()).cast::<idx_t>();
    if adjncy.is_null() {
        libc::free(xadj.as_mut_ptr().cast());
        *r_xadj = std::ptr::null_mut();
        panic!("***Failed to allocate memory for adjncy.");
    }
    *r_adjncy = adjncy;
    mkslice_mut!(adjncy, xadj[ne]);

    for i in (0)..(ne) {
        let nnbrs = FindCommonElements(
            ne as idx_t,
            nn as idx_t,
            i as idx_t,
            eptr[(i + 1) as usize] - eptr[i as usize],
            eind[eptr[i as usize] as usize..].as_ptr(),
            nptr.as_ptr(),
            nind.as_ptr(),
            eptr.as_ptr(),
            ncommon,
            marker.as_mut_ptr(),
            nbrs.as_mut_ptr(),
        );
        for j in (0)..(nnbrs) {
            adjncy[xadj[i as usize] as usize] = nbrs[j as usize];
            xadj[i as usize] += 1;
        }
    }
    util::shift_csr(ne, xadj);

    // gk_free((void **)&nptr, &nind, &marker, &nbrs, LTERM);
}

/*****************************************************************************/
/* This function finds all elements that share at least ncommon nodes with
    the ``query'' element.
*/
/*****************************************************************************/
#[metis_func]
pub extern "C" fn FindCommonElements(
    ne: idx_t,
    nn: idx_t,
    qid: idx_t,
    elen: idx_t,
    eind: *const idx_t,
    nptr: *const idx_t,
    nind: *const idx_t,
    eptr: *const idx_t,
    ncommon: idx_t,
    marker: *mut idx_t,
    nbrs: *mut idx_t,
) -> idx_t {
    // idx_t i, ii, j, jj, k, l, overlap;

    mkslice!(eind, elen);
    mkslice!(eptr, ne + 1);
    mkslice_mut!(marker, ne);
    mkslice_mut!(nbrs, ne);
    // FIXME: pass slices to avoid this computation
    mkslice!(nptr, nn + 1);
    mkslice!(nind, eptr[ne as usize]);

    /* find all elements that share at least one node with qid */
    let mut k = 0;
    for i in (0)..(elen as usize) {
        let j = eind[i as usize] as usize;
        for ii in (nptr[j as usize])..(nptr[(j + 1) as usize]) {
            let jj = nind[ii as usize] as usize;

            if marker[jj as usize] == 0 {
                nbrs[k] = jj as idx_t;
                k += 1;
            }
            marker[jj as usize] += 1;
        }
    }

    /* put qid into the neighbor list (in case it is not there) so that it
    will be removed in the next step */
    if marker[qid as usize] == 0 {
        nbrs[(k) as usize] = qid;
        k += 1;
    }
    marker[qid as usize] = 0;

    /* compact the list to contain only those with at least ncommon nodes */
    let mut j = 0;
    for i in (0)..(k) {
        let l = nbrs[i as usize] as usize;
        let overlap = marker[l];
        if overlap >= ncommon
            || overlap >= elen - 1
            || overlap >= eptr[(l + 1) as usize] - eptr[l as usize] - 1
        {
            nbrs[(j) as usize] = l as idx_t;
            j += 1;
        }
        marker[l as usize] = 0;
    }

    return j;
}

/*****************************************************************************/
/* This function creates the (almost) nodal of a finite element mesh */
/*****************************************************************************/
#[metis_func]
pub extern "C" fn CreateGraphNodal(
    ne: idx_t,
    nn: idx_t,
    eptr: *const idx_t,
    eind: *const idx_t,
    r_xadj: *mut *mut idx_t,
    r_adjncy: *mut *mut idx_t,
) {
    // idx_t i, j, nnbrs;
    // idx_t *nptr, *nind;
    // idx_t *xadj, *adjncy;
    // idx_t *marker, *nbrs;

    let nn = nn as usize;
    let ne = ne as usize;

    mkslice!(eptr, ne + 1);
    mkslice!(eind, ne + 1);

    /* construct the node-element list first */
    let nptr = &mut vec![0 as idx_t; nn + 1];
    let nind = &mut vec![0 as idx_t; eptr[ne] as usize];

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nptr[eind[j as usize] as usize] += 1;
        }
    }
    util::make_csr(nn, nptr);

    for i in (0)..(ne) {
        for j in (eptr[i as usize])..(eptr[(i + 1) as usize]) {
            nind[nptr[eind[j as usize] as usize] as usize] = i as idx_t;
            nptr[eind[j as usize] as usize] += 1;
        }
    }
    util::shift_csr(nn, nptr);

    /* Allocate memory for xadj, since you know its size.
    These are done using standard malloc as they are returned
    to the calling function */
    let xadj = libc::calloc(nn + 1, size_of::<idx_t>()).cast::<idx_t>();
    if xadj.is_null() {
        panic!("***Failed to allocate memory for xadj.");
    }
    *r_xadj = xadj;
    // iset(nn+1, 0, xadj);
    mkslice_mut!(xadj, nn + 1);

    /* allocate memory for working arrays used by FindCommonElements */
    let marker: &mut [idx_t] = &mut vec![0; nn];
    let nbrs: &mut [idx_t] = &mut vec![0; nn];

    for i in (0)..(nn) {
        xadj[i as usize] = FindCommonNodes(
            ne as idx_t,
            nn as idx_t,
            i as idx_t,
            nptr[(i + 1) as usize] - nptr[i as usize],
            nind[nptr[i as usize] as usize..].as_ptr(),
            eptr.as_ptr(),
            eind.as_ptr(),
            marker.as_mut_ptr(),
            nbrs.as_mut_ptr(),
        );
    }
    util::make_csr(nn, xadj);

    /* Allocate memory for adjncy, since you now know its size.
    These are done using standard malloc as they are returned
    to the calling function */
    // not originally zeroed
    let adjncy = libc::calloc(xadj[nn as usize] as usize, size_of::<idx_t>()).cast::<idx_t>();
    if adjncy.is_null() {
        libc::free(xadj.as_mut_ptr().cast());
        *r_xadj = std::ptr::null_mut();
        panic!("***Failed to allocate memory for adjncy.");
    }
    *r_adjncy = adjncy;
    mkslice_mut!(adjncy, xadj[nn]);

    for i in (0)..(nn) {
        let nnbrs = FindCommonNodes(
            ne as idx_t,
            nn as idx_t,
            i as idx_t,
            nptr[(i + 1) as usize] - nptr[i as usize],
            nind[nptr[i as usize] as usize..].as_ptr(),
            eptr.as_ptr(),
            eind.as_ptr(),
            marker.as_mut_ptr(),
            nbrs.as_mut_ptr(),
        );
        for j in (0)..(nnbrs) {
            adjncy[xadj[i as usize] as usize] = nbrs[j as usize];
            xadj[i as usize] += 1;
        }
    }
    util::shift_csr(nn, xadj);

    // gk_free((void **)&nptr, &nind, &marker, &nbrs, LTERM);
}

/*****************************************************************************/
/* This function finds the union of nodes that are in the same elements with
    the ``query'' node.
*/
/*****************************************************************************/
#[metis_func]
pub extern "C" fn FindCommonNodes(
    ne: idx_t,
    nn: idx_t,
    qid: idx_t,
    nelmnts: idx_t,
    elmntids: *const idx_t,
    eptr: *const idx_t,
    eind: *const idx_t,
    marker: *mut idx_t,
    nbrs: *mut idx_t,
) -> idx_t {
    // idx_t i, ii, j, jj, k;
    mkslice!(elmntids, nelmnts);

    mkslice_mut!(marker, nn);
    mkslice_mut!(nbrs, nn);
    mkslice!(eptr, ne + 1);
    mkslice!(eind, ne + 1);

    /* find all nodes that share at least one element with qid */
    marker[qid as usize] = 1; /* this is to prevent self-loops */
    let mut k = 0;
    for i in (0)..(nelmnts) {
        let j = elmntids[i as usize] as usize;
        for ii in (eptr[j as usize])..(eptr[(j + 1) as usize]) {
            let jj = eind[ii as usize] as usize;
            if marker[jj as usize] == 0 {
                nbrs[(k) as usize] = jj as idx_t;
                k += 1;
                marker[jj as usize] = 1;
            }
        }
    }

    /* reset the marker */
    marker[qid as usize] = 0;
    for i in (0)..(k) {
        marker[nbrs[i as usize] as usize] = 0;
    }

    return k;
}

/*************************************************************************/
/* This function creates and initializes a mesh_t structure */
/*************************************************************************/
#[metis_func]
pub extern "C" fn CreateMesh() -> *mut mesh_t {
    // mesh_t *mesh;

    let mesh = gk_malloc(size_of::<mesh_t>(), c"CreateMesh: mesh".as_ptr()).cast();

    *mesh = std::mem::zeroed();

    return mesh;
}

/*************************************************************************/
/* This function deallocates any memory stored in a mesh */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FreeMesh(r_mesh: *mut *mut mesh_t) {
    // mesh_t *mesh = *r_mesh;

    let mesh = &mut **r_mesh;
    gk::free_ref(&mut mesh.eptr);
    gk::free_ref(&mut mesh.eind);
    gk::free_ref(&mut mesh.ewgt);
    gk::free_ref(&mut *r_mesh);
    // gk_free(&mesh.eptr, &mesh.eind, &mesh.ewgt, &mesh, LTERM);

    *r_mesh = std::ptr::null_mut();
}
