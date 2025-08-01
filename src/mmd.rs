//! # Multiple Minimum External Degree (mmd)
//!
//! This was originally a FORTRAN subroutine in SPARSPAK (`genqmd` in particular). Then it was
//! ported to C. Now it has been ported to Rust. It was cursed to begin with and now it's so much
//! worse.
//!
//! The original package is somewhat hard to find, but can be a good reference. I've collected
//! links to all relevant pages that may give some insight
//!
//! - [SPARSPAK homepage] ([archive][homepage-archive]) -- not very helpful
//! - [Subroutine explanation (F77)] ([archive][f77-exp-archive])
//!     - Contains FORTRAN 77 source code link
//!     - Spells it as SPARSEPAK (I believe this is a mistake)
//! - [Subroutine explaination (F90) (archive)][f90-exp-archive] (original is gone)
//!     - Contains FORTRAN 90 source code link
//!     - Spells it as SPARSEPAK (I believe this is a mistake)
//!     - Claims future version have unfree Licensing
//! - [Julia implementation][julia]
//! - [A Fast Implementation of the Minimum Degree Algorithm Using Quotient Graphs][MMD intro]
//! - [Modification of the minimum degree algorithm by multiple elimination][MMD modification]
//! - [The evolution of the minimum degree algorithm][MMD evolution]
//!     - Paywalled :(
//!
//! [SPARSPAK homepage]: https://cs.uwaterloo.ca/~jageorge/Sparspak/sparspak.html
//! [homepage-archive]: https://web.archive.org/web/20250505112845/https://cs.uwaterloo.ca/~jageorge/Sparspak/sparspak.html
//! [Subroutine explanation (F77)]: https://people.sc.fsu.edu/~jburkardt/f77_src/sparsepak/sparsepak.html
//! [f77-exp-archive]: https://web.archive.org/web/20220628125400/https://people.sc.fsu.edu/~jburkardt/f77_src/sparsepak/sparsepak.html
//! [f90-exp-archive]: https://web.archive.org/web/20120105012708/http://people.sc.fsu.edu/~jburkardt/f_src/sparsepak/sparsepak.html
//! [julia]: https://github.com/PetrKryslUCSD/Sparspak.jl/blob/b49510aa705ffc4407e624c0ec36fe6da7517e7d/src/SparseSpdMethod/SpkMMD.jl
//! [MMD intro]: https://doi.org/10.1145/355900.355906
//! [MMD modification]: https://doi.org/10.1145/214392.214398
//! [MMD evolution]: https://doi.org/10.1137/1031049
//!
//! ## Hacks
//!
//! There are/were a number of hacks that were introduced over the years to deal with this
//! routine's history. For the sake of the optimizer, I want to get rid of as many as possible.
//!
//! ### 1-based indexing
//!
//! Since fortran is 1-indexed and C is not, the FORTRAN-to-C port simply decremented all of the
//! pointers passed to it. Strictly speaking, this causes [immediate ub] (run with miri), though as
//! of writing (2025-07-11), no compiler I could find acts on even a very simple case. That being
//! said, I'll use a slice wrapper that adjusts the indices by one. I expect that it won't be too
//! horrible to remove most of the required one-indexing using Rust features, but I may eat those
//! words.
//!
//! Outside of general ease of use (and hopefully correctness), getting rid of 1-based indexing
//! should have optimization benefits, since inclusive ranges (`a..=b`) have a somewhat rocky history
//! of optimization and lack a lot of the specialization that have been given to the "normal" half
//! open range (`a..b`).
//!
//! ### Out-of-bounds accesses
//!
//! It appears that at one point, this routine was liable to perform out-of-bounds accesses.
//!
//! Frustratingly, METIS used svn(?) when these hacks were introduced, so it's difficult to find
//! the historical context of their use. Regardless, this can be seen in [`ometis::MMDOrder`] where
//! all the working arrays are allocated with 5 extra elements. [This github discussion][jsf67]
//! was very helpful for piecing together some of the history here (The user who provided insight
//! into this bug, `jsf67` appears to have prior experience working with genmmd, but they have what
//! appears to be zero online presence outside of those comments).
//!
//! Regardless, this seems like it will make it exceptionally difficult to test, since we will
//! probably run into out-of-bounds panics that were present in the original. We could address this
//! by adopting the 5-extra-elements hack, but then we can't use slice iteration to work around
//! FORTRAN indexing.
//!
//! ### Go To Statement Considered Harmful
//!
//! We've all read [the famous letter], but it was all very abstract until I started work on
//! porting this file to Rust did I appreciate it. Thankfully, we can work around it with labeled
//! blocks, though it does look pretty ugly.
//!
//! [the famous letter]: https://doi.org/10.1145/362929.362947
//! [jsf67]: https://github.com/KarypisLab/METIS/issues/46
//! [immediate ub]: https://play.rust-lang.org/?version=stable&mode=debug&edition=2024&gist=34331f6a82be75de727ef2479302c7dd
//!
//! ## Strategy
//!
//! My actual porting strategy here is going to be different from the rest of METIS. Namely, only
//! the top-level [`genmmd`] will be marked with `#[metis_func]` and will be tested with the A/B
//! tests that use `METIS_*` functions. All other tests will be done directly on [`genmmd`]. I'll
//! exploit the fact that genmmd is only used on small graphs (nvtxs <= MMDSWITCH), is entirely
//! deterministic, and is very self-contained, to write specific test cases (and maybe eventually
//! fuzzing)
//!

use crate::*;

/*
 * mmd.c
 *
 * **************************************************************
 * The following C function was developed from a FORTRAN subroutine
 * in SPARSPAK written by Eleanor Chu, Alan George, Joseph Liu
 * and Esmond Ng.
 *
 * The FORTRAN-to-C transformation and modifications such as dynamic
 * memory allocation and deallocation were performed by Chunguang
 * Sun.
 * **************************************************************
 *
 * Taken from SMMS, George 12/13/94
 *
 * The meaning of invperm, and perm vectors is different from that
 * in genqmd_ of SparsPak
 *
 * $Id: mmd.c 22385 2019-06-03 22:08:48Z karypis $
 */

/// wrapper for slices that use FORTRAN indexing
#[repr(transparent)]
struct Fslice([idx_t]);

impl Fslice {
    #![allow(dead_code)]

    fn new(slice: &[idx_t]) -> &Fslice {
        unsafe { std::mem::transmute(slice) }
    }

    fn new_mut(slice: &mut [idx_t]) -> &mut Fslice {
        unsafe { std::mem::transmute(slice) }
    }

    unsafe fn from_raw_parts<'a>(p: *const idx_t, len: usize) -> &'a Fslice {
        std::mem::transmute(std::slice::from_raw_parts(p, len))
    }

    unsafe fn from_raw_parts_mut<'a>(p: *mut idx_t, len: usize) -> &'a mut Fslice {
        std::mem::transmute(std::slice::from_raw_parts_mut(p, len))
    }

    fn as_slice(&self) -> &[idx_t] {
        &self.0
    }

    fn as_slice_mut(&mut self) -> &mut [idx_t] {
        &mut self.0
    }
}

trait FsliceIndex {
    fn to_slice_index(self) -> Self;
}

impl<I: FsliceIndex + std::slice::SliceIndex<[idx_t]>> std::ops::Index<I> for Fslice {
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        &self.0[index.to_slice_index()]
    }
}

impl<I: FsliceIndex + std::slice::SliceIndex<[idx_t]>> std::ops::IndexMut<I> for Fslice {
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.0[index.to_slice_index()]
    }
}

impl FsliceIndex for std::ops::Range<usize> {
    fn to_slice_index(self) -> Self {
        (self.start - 1)..(self.end - 1)
    }
}

impl FsliceIndex for std::ops::RangeToInclusive<usize> {
    fn to_slice_index(self) -> Self {
        ..=(self.end - 1)
    }
}

impl FsliceIndex for std::ops::RangeInclusive<usize> {
    fn to_slice_index(self) -> Self {
        (self.start() - 1)..=(self.end() - 1)
    }
}

impl FsliceIndex for usize {
    fn to_slice_index(self) -> Self {
        self - 1
    }
}

/*************************************************************************
*  genmmd  -- multiple minimum external degree
*  purpose -- this routine implements the minimum degree
*     algorithm. it makes use of the implicit representation
*     of elimination graphs by quotient graphs, and the notion
*     of indistinguishable nodes. It also implements the modifications
*     by multiple elimination and minimum external degree.
*     Caution -- the adjacency vector adjncy will be destroyed.
*  Input parameters --
*     neqns -- number of equations.
*     (xadj, adjncy) -- the adjacency structure.
*     delta  -- tolerance value for multiple elimination.
*     maxint -- maximum machine representable (short) integer
*               (any smaller estimate will do) for marking nodes.
*  Output parameters --
*     perm -- the minimum degree ordering.
*     invp -- the inverse of perm.
*     *ncsub -- an upper bound on the number of nonzero subscripts
*               for the compressed storage scheme.
*  Working parameters --
*     head -- vector for head of degree lists.
*     invp  -- used temporarily for degree forward link.
*     perm  -- used temporarily for degree backward link.
*     qsize -- vector for size of supernodes.
*     list -- vector for temporary linked lists.
*     marker -- a temporary marker vector.
*  Subroutines used -- mmdelm, mmdint, mmdnum, mmdupd.
**************************************************************************/
/// this is the entry point to for `genmmd`. Right now it just forwards to [`genmmd_rs_entry`]
#[metis_func]
pub extern "C" fn genmmd(
    neqns: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,

    invp: *mut idx_t,
    perm: *mut idx_t,
    delta: idx_t,
    head: *mut idx_t,
    qsize: *mut idx_t,
    list: *mut idx_t,
    marker: *mut idx_t,

    maxint: idx_t,
    ncsub: *mut idx_t,
) {
    let xadj = std::slice::from_raw_parts_mut(xadj, neqns as usize + 1);
    let adjncy = std::slice::from_raw_parts_mut(adjncy, xadj[neqns as usize] as usize);

    // do fortran renumbering -- this was previously done in MMDOrder
    for i in &mut *xadj {
        *i += 1;
    }
    for i in &mut *adjncy {
        *i += 1;
    }

    let xadj = Fslice::new_mut(xadj);
    let adjncy = Fslice::new_mut(adjncy);

    debug_assert_eq!(*xadj.as_slice().last().unwrap(), xadj[neqns as usize + 1]); // sanity check


    // I might need to include the length hack provided in MMDOrder (+5 to len) for each of these
    // scratch slices
    let invp = Fslice::from_raw_parts_mut(invp, neqns as usize);
    let perm = Fslice::from_raw_parts_mut(perm, neqns as usize);
    let head = Fslice::from_raw_parts_mut(head, neqns as usize);
    let qsize = Fslice::from_raw_parts_mut(qsize, neqns as usize);
    let list = Fslice::from_raw_parts_mut(list, neqns as usize);
    let marker = Fslice::from_raw_parts_mut(marker, neqns as usize);

    assert_eq!(maxint, idx_t::MAX);

    let res = genmmd_rs_entry(
        neqns, xadj, adjncy, invp, perm, delta, head, qsize, list, marker,
    );
    *ncsub = res;

    // Relabel the vertices so that it starts from 0
    // This was previously done in MMDOrder
    for i in xadj.as_slice_mut() {
        *i -= 1;
    }
    for i in adjncy.as_slice_mut() {
        *i -= 1;
    }
}

fn genmmd_rs_entry(
    neqns: idx_t,
    xadj: &mut Fslice,
    adjncy: &mut Fslice,
    invp: &mut Fslice,
    perm: &mut Fslice,
    delta: idx_t,
    head: &mut Fslice,
    qsize: &mut Fslice,
    list: &mut Fslice,
    marker: &mut Fslice,
) -> idx_t {
    // idx_t  ehead, i, mdeg, mdlmt, mdeg_node, nextmd, num, tag;

    if neqns <= 0 {
        return 0;
    }

    /* initialization for the minimum degree algorithm */
    let mut ncsub = 0;
    mmdint(neqns, xadj, head, invp, perm, qsize, list, marker);

    /* 'num' counts the number of ordered nodes plus 1 */
    let mut num = 1;

    /* eliminate all isolated nodes */
    let mut nextmd = head[1 as usize];
    while nextmd > 0 {
        let mdeg_node = nextmd;
        nextmd = invp[mdeg_node as usize];
        marker[mdeg_node as usize] = idx_t::MAX;
        invp[mdeg_node as usize] = -num;
        num += 1;
    }

    /* search for node of the minimum degree. 'mdeg' is the current */
    /* minimum degree; 'tag' is used to facilitate marking nodes.   */
    if num > neqns {
        mmdnum(neqns, perm, invp, qsize);
        return 0;
    }

    let mut tag = 1;
    head[1 as usize] = 0;

    let mut mdeg = 2;

    /* infinite loop here */
    'outer: loop {
        while head[mdeg as usize] <= 0 {
            mdeg += 1;
        }

        /* use value of 'delta' to set up 'mdlmt', which governs */
        /* when a degree update is to be performed.              */
        //mdlmt = mdeg + delta;
        // the need for gk_min() was identified by jsf67
        // https://github.com/KarypisLab/METIS/issues/46
        let mdlmt = neqns.min(mdeg + delta);
        let mut ehead = 0;

        // merged n500 and n900
        // goto n500 => continue n500
        // goto n900 => break n500
        'n500: loop {
            let mut mdeg_node = head[mdeg as usize];
            while mdeg_node <= 0 {
                mdeg += 1;

                if mdeg > mdlmt {
                    break 'n500;
                }
                mdeg_node = head[mdeg as usize];
            }

            /* remove 'mdeg_node' from the degree structure */
            nextmd = invp[mdeg_node as usize];
            head[mdeg as usize] = nextmd;
            if nextmd > 0 {
                perm[nextmd as usize] = -mdeg;
            }
            invp[mdeg_node as usize] = -num;
            ncsub += mdeg + qsize[mdeg_node as usize] - 2;
            if num + qsize[mdeg_node as usize] > neqns {
                break 'outer;
            }

            /*  eliminate 'mdeg_node' and perform quotient graph */
            /*  transformation. reset 'tag' value if necessary.    */
            tag += 1;
            if tag >= idx_t::MAX {
                tag = 1;
                for i in (1)..=(neqns) {
                    if marker[i as usize] < idx_t::MAX {
                        marker[i as usize] = 0;
                    }
                }
            };

            mmdelm(
                mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, tag,
            );

            num += qsize[mdeg_node as usize];
            list[mdeg_node as usize] = ehead;
            ehead = mdeg_node;
            if delta >= 0 {
                continue 'n500;
            } else {
                break 'n500;
            }
        }

        // was merged with n500
        // n900:

        /* update degrees of the nodes involved in the  */
        /* minimum degree nodes elimination.            */
        if num > neqns {
            break 'outer;
        }

        mmdupd(
            ehead, neqns, xadj, adjncy, delta, &mut mdeg, head, invp, perm, qsize, list, marker, &mut tag,
        );
    }

    mmdnum(neqns, perm, invp, qsize);

    return ncsub;
}

/**************************************************************************
*           mmdelm ...... multiple minimum degree elimination
* Purpose -- This routine eliminates the node mdeg_node of minimum degree
*     from the adjacency structure, which is stored in the quotient
*     graph format. It also transforms the quotient graph representation
*     of the elimination graph.
* Input parameters --
*     mdeg_node -- node of minimum degree.
*     tag    -- tag value.
* Updated parameters --
*     (xadj, adjncy) -- updated adjacency structure.
*     (head, forward, backward) -- degree doubly linked structure.
*     qsize -- size of supernode.
*     marker -- marker vector.
*     list -- temporary linked list of eliminated nabors.
***************************************************************************/
fn mmdelm(
    mdeg_node: idx_t,
    xadj: &mut Fslice,
    adjncy: &mut Fslice,
    head: &mut Fslice,
    forward: &mut Fslice,
    backward: &mut Fslice,
    qsize: &mut Fslice,
    list: &mut Fslice,
    marker: &mut Fslice,
    tag: idx_t,
) {
    // idx_t   element, i,   istop, istart, j,
    //       jstop, jstart, link,
    //       nabor, node, npv, nqnbrs, nxnode,
    //       pvnode, rlmt, rloc, rnode, xqnbr;

    /* find the reachable set of 'mdeg_node' and */
    /* place it in the data structure.           */
    marker[mdeg_node as usize] = tag;
    let istart = xadj[mdeg_node as usize];
    let istop = xadj[(mdeg_node + 1) as usize] - 1;

    /* 'element' points to the beginning of the list of  */
    /* eliminated nabors of 'mdeg_node', and 'rloc' gives the */
    /* storage location for the next reachable node.   */
    let mut element = 0;
    let mut rloc = istart;
    let mut rlmt = istop;
    for i in istart..=istop {
        let nabor = adjncy[i as usize];
        if nabor == 0 {
            break;
        }
        if marker[nabor as usize] < tag {
            marker[nabor as usize] = tag;
            if forward[nabor as usize] < 0 {
                list[nabor as usize] = element;
                element = nabor;
            } else {
                adjncy[rloc as usize] = nabor;
                rloc += 1;
            };
        }; /* end of -- if -- */
    } /* end of -- for -- */

    /* merge with reachable nodes from generalized elements. */
    while element > 0 {
        adjncy[rlmt as usize] = -element;
        let mut link = element;

        'n400: loop {
            let jstart = xadj[link as usize];
            let jstop = xadj[(link + 1) as usize] - 1;
            for j in jstart..=jstop {
                let node = adjncy[j as usize];
                link = -node;
                if node < 0 {
                    continue 'n400; // goto
                }
                if node == 0 {
                    break;
                }
                if marker[node as usize] < tag && forward[node as usize] >= 0 {
                    marker[node as usize] = tag;
                    /*use storage from eliminated nodes if necessary.*/
                    while rloc >= rlmt {
                        link = -adjncy[rlmt as usize];
                        rloc = xadj[link as usize];
                        rlmt = xadj[(link + 1) as usize] - 1;
                    }
                    adjncy[rloc as usize] = node;
                    rloc += 1;
                };
            } /* end of -- for ( j = jstart; -- */
            break 'n400; // synthetic
        }

        element = list[element as usize];
    } /* end of -- while  element > 0  -- */

    if rloc <= rlmt {
        adjncy[rloc as usize] = 0;
    }
    /* for each node in the reachable set, do the following. */
    let mut link = mdeg_node;

    'n1100: loop {
        let istart = xadj[link as usize];
        let istop = xadj[(link + 1) as usize] - 1;
        for i in istart..=istop {
            let rnode = adjncy[i as usize];
            link = -rnode;
            if rnode < 0 {
                continue 'n1100; // goto
            }
            if rnode == 0 {
                return;
            }

            /* 'rnode' is in the degree list structure. */
            let pvnode = backward[rnode as usize];
            if pvnode != 0 && (pvnode != (-idx_t::MAX)) {
                /* then remove 'rnode' from the structure. */
                let nxnode = forward[rnode as usize];
                if nxnode > 0 {
                    backward[nxnode as usize] = pvnode;
                }
                if pvnode > 0 {
                    forward[pvnode as usize] = nxnode;
                }
                let npv = -pvnode;
                if pvnode < 0 {
                    head[npv as usize] = nxnode;
                };
            }

            /* purge inactive quotient nabors of 'rnode'. */
            let jstart = xadj[rnode as usize];
            let jstop = xadj[(rnode + 1) as usize] - 1;
            let mut xqnbr = jstart;
            for j in jstart..=jstop {
                let nabor = adjncy[j as usize];
                if nabor == 0 {
                    break;
                }
                if marker[nabor as usize] < tag {
                    adjncy[xqnbr as usize] = nabor;
                    xqnbr += 1;
                };
            }

            /* no active nabor after the purging. */
            let nqnbrs = xqnbr - jstart;
            if nqnbrs <= 0 {
                /* merge 'rnode' with 'mdeg_node'. */
                qsize[mdeg_node as usize] += qsize[rnode as usize];
                qsize[rnode as usize] = 0;
                marker[rnode as usize] = idx_t::MAX;
                forward[rnode as usize] = -mdeg_node;
                backward[rnode as usize] = -idx_t::MAX;
            } else {
                /* flag 'rnode' for degree update, and  */
                /* add 'mdeg_node' as a nabor of 'rnode'.      */
                forward[rnode as usize] = nqnbrs + 1;
                backward[rnode as usize] = 0;
                adjncy[xqnbr as usize] = mdeg_node;
                xqnbr += 1;
                if xqnbr <= jstop {
                    adjncy[xqnbr as usize] = 0;
                }
            };
        }

        break 'n1100; // synthetic
    }
    return;
}

/***************************************************************************
*    mmdint ---- mult minimum degree initialization
*    purpose -- this routine performs initialization for the
*       multiple elimination version of the minimum degree algorithm.
*    input parameters --
*       neqns  -- number of equations.
*       (xadj, adjncy) -- adjacency structure.
*    output parameters --
*       (head, dfrow, backward) -- degree doubly linked structure.
*       qsize -- size of supernode ( initialized to one).
*       list -- linked list.
*       marker -- marker vector.
****************************************************************************/
fn mmdint(
    neqns: idx_t,
    xadj: &mut Fslice,
    // adjncy: &mut Fslice,
    head: &mut Fslice,
    forward: &mut Fslice,
    backward: &mut Fslice,
    qsize: &mut Fslice,
    list: &mut Fslice,
    marker: &mut Fslice,
) -> idx_t {
    // idx_t fnode, ndeg, node;

    for node in (1)..=(neqns) {
        head[node as usize] = 0;
        qsize[node as usize] = 1;
        marker[node as usize] = 0;
        list[node as usize] = 0;
    }

    /* initialize the degree doubly linked lists. */
    for node in (1)..=(neqns) {
        let ndeg = xadj[(node + 1) as usize] - xadj[node as usize] + 1;
        let fnode = head[ndeg as usize];
        forward[node as usize] = fnode;
        head[ndeg as usize] = node;
        if fnode > 0 {
            backward[fnode as usize] = node;
        }
        backward[node as usize] = -ndeg;
    }

    return 0;
}

/****************************************************************************
* mmdnum --- multi minimum degree numbering
* purpose -- this routine performs the final step in producing
*    the permutation and inverse permutation vectors in the
*    multiple elimination version of the minimum degree
*    ordering algorithm.
* input parameters --
*     neqns -- number of equations.
*     qsize -- size of supernodes at elimination.
* updated parameters --
*     invp -- inverse permutation vector. on input,
*             if qsize[node as usize] = 0, then node has been merged
*             into the node -invp[node as usize]; otherwise,
*            -invp[node as usize] is its inverse labelling.
* output parameters --
*     perm -- the permutation vector.
****************************************************************************/
fn mmdnum(neqns: idx_t, perm: &mut Fslice, invp: &mut Fslice, qsize: &mut Fslice) {
    // idx_t father, nextf, node, nqsize, num, root;

    for node in 1..=neqns {
        let nqsize = qsize[node as usize];
        if nqsize <= 0 {
            perm[node as usize] = invp[node as usize];
        }
        if nqsize > 0 {
            perm[node as usize] = -invp[node as usize];
        }
    }

    /* for each node which has been merged, do the following. */
    for node in 1..=neqns {
        if perm[node as usize] <= 0 {
            /* trace the merged tree until one which has not */
            /* been merged, call it root.                    */
            let mut father = node;
            while perm[father as usize] <= 0 {
                father = -perm[father as usize];
            }

            /* number node after root. */
            let root = father;
            let num = perm[root as usize] + 1;
            invp[node as usize] = -num;
            perm[root as usize] = num;

            /* shorten the merged tree. */
            let mut father = node;
            let mut nextf = -perm[father as usize];
            while nextf > 0 {
                perm[father as usize] = -root;
                father = nextf;
                nextf = -perm[father as usize];
            }
        }
    }

    /* ready to compute perm. */
    for node in 1..=neqns {
        let num = -invp[node as usize];
        invp[node as usize] = num;
        perm[num as usize] = node;
    }
    return;
}

/****************************************************************************
* mmdupd ---- multiple minimum degree update
* purpose -- this routine updates the degrees of nodes after a
*            multiple elimination step.
* input parameters --
*    ehead -- the beginning of the list of eliminated nodes
*             (i.e., newly formed elements).
*    neqns -- number of equations.
*    (xadj, adjncy) -- adjacency structure.
*    delta -- tolerance value for multiple elimination.
* updated parameters --
*    mdeg -- new minimum degree after degree update.
*    (head, forward, backward) -- degree doubly linked structure.
*    qsize -- size of supernode.
*    list -- marker vector for degree update.
*    *tag   -- tag value.
****************************************************************************/
fn mmdupd(
    ehead: idx_t,
    neqns: idx_t,
    xadj: &mut Fslice,
    adjncy: &mut Fslice,
    delta: idx_t,
    mdeg: &mut idx_t,
    head: &mut Fslice,
    forward: &mut Fslice,
    backward: &mut Fslice,
    qsize: &mut Fslice,
    list: &mut Fslice,
    marker: &mut Fslice,
    tag: &mut idx_t,
) {
    // idx_t  deg, deg0, element, enode, fnode, i, iq2, istop,
    //      istart, j, jstop, jstart, link, mdeg0, mtag, nabor,
    //      node, q2head, qxhead;
    //
    let mdeg0 = *mdeg + delta;
    let mut element = ehead;

    // n100:
    loop {
        let mut mtag;

        // header (so n2300 is only ever broken out of)
        'n2300: {
            if element <= 0 {
                return;
            }

            /* for each of the newly formed element, do the following. */
            /* reset tag value if necessary.                           */
            mtag = *tag + mdeg0;
            if mtag >= idx_t::MAX {
                *tag = 1;
                for i in 1..=neqns {
                    if marker[i as usize] < idx_t::MAX {
                        marker[i as usize] = 0;
                    }
                }
                mtag = *tag + mdeg0;
            };

            /* create two linked lists from nodes associated with 'element': */
            /* one with two nabors (q2head) in the adjacency structure, and the*/
            /* other with more than two nabors (qxhead). also compute 'deg0',*/
            /* number of nodes in this element.                              */
            let mut q2head = 0;
            let mut qxhead = 0;
            let mut deg0 = 0;
            let mut link = element;

            'n400: loop {
                let istart = xadj[link as usize];
                let istop = xadj[(link + 1) as usize] - 1;
                for i in istart..=istop {
                    let enode = adjncy[i as usize];
                    link = -enode;
                    if enode < 0 {
                        continue 'n400; // goto
                    }
                    if enode == 0 {
                        break;
                    }
                    if qsize[enode as usize] != 0 {
                        deg0 += qsize[enode as usize];
                        marker[enode as usize] = mtag;

                        /*'enode' requires a degree update*/
                        if backward[enode as usize] == 0 {
                            /* place either in qxhead or q2head list. */
                            if forward[enode as usize] != 2 {
                                list[enode as usize] = qxhead;
                                qxhead = enode;
                            } else {
                                list[enode as usize] = q2head;
                                q2head = enode;
                            }
                        }
                    }
                }

                break 'n400; // synthetic fallthrough
            }

            /* for each node in q2 list, do the following. */
            let mut enode = q2head;
            let mut iq2 = 1;
            let mut deg;

            'n900: loop {
                // instead of breaking to n2200 before the n1600 loop, we duplicate the (short)
                // code of n2200, replacing its final 'goto n1600' with 'break 1600_pre'
                'n1600_pre: {
                    'n1500: {
                        'n2100_dup: {
                            if enode <= 0 {
                                break 'n1500; // goto (down)
                            }

                            if backward[enode as usize] != 0 {
                                // goto n2200 (duplicated)
                                /* get next enode in current element. */
                                enode = list[enode as usize];
                                if iq2 == 1 {
                                    continue 'n900; // goto
                                }
                                break 'n1600_pre;
                            }

                            (*tag) += 1;
                            deg = deg0;

                            /* identify the other adjacent element nabor. */
                            let istart = xadj[enode as usize];
                            let mut nabor = adjncy[istart as usize];
                            if nabor == element {
                                nabor = adjncy[(istart + 1) as usize];
                            }
                            link = nabor;
                            if forward[nabor as usize] >= 0 {
                                /* nabor is uneliminated, increase degree count. */
                                deg += qsize[nabor as usize];
                                break 'n2100_dup; // goto dup (down)
                            };

                            /* the nabor is eliminated. for each node in the 2nd element */
                            /* do the following.                                         */
                            'n1000: loop {
                                let istart = xadj[link as usize];
                                let istop = xadj[(link + 1) as usize] - 1;
                                for i in istart..=istop {
                                    let node = adjncy[i as usize];
                                    link = -node;
                                    if node != enode {
                                        if node < 0 {
                                            continue 'n1000; // goto
                                        }
                                        if node == 0 {
                                            break 'n2100_dup; // goto dup (down)
                                        }
                                        if qsize[node as usize] != 0 {
                                            if marker[node as usize] < *tag {
                                                /* 'node' is not yet considered. */
                                                marker[node as usize] = *tag;
                                                deg += qsize[node as usize];
                                            } else {
                                                if backward[node as usize] == 0 {
                                                    if forward[node as usize] == 2 {
                                                        /* 'node' is indistinguishable from 'enode'.*/
                                                        /* merge them into a new supernode.         */
                                                        qsize[enode as usize] +=
                                                            qsize[node as usize];
                                                        qsize[node as usize] = 0;
                                                        marker[node as usize] = idx_t::MAX;
                                                        forward[node as usize] = -enode;
                                                        backward[node as usize] = -idx_t::MAX;
                                                    } else {
                                                        /* 'node' is outmacthed by 'enode' */
                                                        if backward[node as usize] == 0 {
                                                            backward[node as usize] = -idx_t::MAX;
                                                        }
                                                    };
                                                }
                                            }
                                        }
                                    }
                                }

                                break 'n1000; // synthetic fallthrough
                            }
                        } // n2100 dup:
                        {
                            // n2100: duplicated and moved earlier

                            /* update external degree of 'enode' in degree structure, */
                            /* and '*mdeg' if necessary.                     */
                            deg = deg - qsize[enode as usize] + 1;
                            let fnode = head[deg as usize];
                            forward[enode as usize] = fnode;
                            backward[enode as usize] = -deg;
                            if fnode > 0 {
                                backward[fnode as usize] = enode;
                            }
                            head[deg as usize] = enode;
                            if deg < *mdeg {
                                *mdeg = deg;
                            }
                            // goto n2200 (duplicated)
                            /* get next enode in current element. */
                            enode = list[enode as usize];
                            if iq2 == 1 {
                                continue 'n900; // goto
                            }
                            break 'n1600_pre;
                        }
                    } // n1500:

                    /* for each 'enode' in the 'qx' list, do the following. */
                    enode = qxhead;
                    iq2 = 0;
                } // n1600: (pre)

                'n1600: loop {
                    if enode <= 0 {
                        break 'n2300; // goto (down)
                    }
                    if backward[enode as usize] != 0 {
                        // goto n2200 (duplicated)
                        /* get next enode in current element. */
                        enode = list[enode as usize];
                        if iq2 == 1 {
                            continue 'n900; // goto
                        }
                        continue 'n1600;
                    }
                    (*tag) += 1;
                    deg = deg0;

                    /*for each unmarked nabor of 'enode', do the following.*/
                    let istart = xadj[enode as usize];
                    let istop = xadj[(enode + 1) as usize] - 1;
                    for i in istart..=istop {
                        let nabor = adjncy[i as usize];
                        if nabor == 0 {
                            break;
                        }
                        if marker[nabor as usize] < *tag {
                            marker[nabor as usize] = *tag;
                            link = nabor;
                            if forward[nabor as usize] >= 0 {
                                /*if uneliminated, include it in deg count.*/
                                deg += qsize[nabor as usize];
                            } else {
                                'n1700: loop {
                                    /* if eliminated, include unmarked nodes in this*/
                                    /* element into the degree count.             */
                                    let jstart = xadj[link as usize];
                                    let jstop = xadj[(link + 1) as usize] - 1;
                                    for j in jstart..=jstop {
                                        let node = adjncy[j as usize];
                                        link = -node;
                                        if node < 0 {
                                            continue 'n1700; // goto
                                        }
                                        if node == 0 {
                                            break;
                                        }
                                        if marker[node as usize] < *tag {
                                            marker[node as usize] = *tag;
                                            deg += qsize[node as usize];
                                        }
                                    }

                                    break 'n1700; // synthetic fallthrough
                                }
                            }
                        }
                    }

                    // n2100: (original - now only reached by fallthrough)

                    /* update external degree of 'enode' in degree structure, */
                    /* and '*mdeg' if necessary.                     */
                    deg = deg - qsize[enode as usize] + 1;
                    let fnode = head[deg as usize];
                    forward[enode as usize] = fnode;
                    backward[enode as usize] = -deg;
                    if fnode > 0 {
                        backward[fnode as usize] = enode;
                    }
                    head[deg as usize] = enode;
                    if deg < *mdeg {
                        *mdeg = deg;
                    }
                    // n2200: (original, now fallthrough only)

                    /* get next enode in current element. */
                    enode = list[enode as usize];
                    if iq2 == 1 {
                        continue 'n900; // goto
                    }

                    continue 'n1600; // goto, superfluous (it's a loop!)
                } // bottom of n1600

                #[expect(unreachable_code)]
                continue 'n900; // goto, superfluous (it's a loop!)
            } // bottom of n900

            // this is also unreachable
            break 'n2300; // synthetic, superfluous (not a loop)
        } // n2300:

        /* get next element in the list. */
        *tag = mtag;
        element = list[element as usize];
    }
}

#[cfg(test)]
mod tests {
    use crate::{dyncall::ab_test_single_eq, graph_gen::Csr, tests::{ab_test_partition_test_graphs, TestGraph}};

    use super::*;

    #[test]
    fn ab_genmmd_through_metis() {
        ab_test_partition_test_graphs("genmmd:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_genmmd_through_metis_cc() {
        ab_test_partition_test_graphs("genmmd:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g.set_ccorder(true);
            g
        });
    }

    fn mmd_ab_test(graph: Csr, delta: idx_t) {
        ab_test_single_eq("genmmd:rs", || {
            let csr = graph.clone();
            let nvtxs = csr.nvtxs();
            let (mut xadj, mut adjncy) = csr.into_parts();
            let mut perm = vec![0; nvtxs + 5];
            let mut iperm = vec![0; nvtxs + 5];
            let mut head = vec![0; nvtxs + 5];
            let mut qsize = vec![0; nvtxs + 5];
            let mut list = vec![0; nvtxs + 5];
            let mut marker = vec![0; nvtxs + 5];
            let mut nofsub = 0;
            unsafe { genmmd(
                nvtxs as idx_t,
                xadj.as_mut_ptr(),
                adjncy.as_mut_ptr(),
                iperm.as_mut_ptr(),
                perm.as_mut_ptr(),
                delta,
                head.as_mut_ptr(),
                qsize.as_mut_ptr(),
                list.as_mut_ptr(),
                marker.as_mut_ptr(),
                idx_t::MAX,
                &mut nofsub,
            )};
            // idk if I want to compare literally everything?
            (
                nofsub,
                iperm,
                perm,
                head,
                qsize,
                list,
                marker
            )
        });
    }

    #[test]
    fn ab_mmd_direct_simple() {
        let g = Csr::test_graph(TestGraph::Webbase2004);
        mmd_ab_test(g, 1);
    }

    #[test]
    fn ab_mmd_direct_larger_delta() {
        let g = Csr::test_graph(TestGraph::Webbase2004);
        mmd_ab_test(g.clone(), 2);
        mmd_ab_test(g.clone(), 3);
        mmd_ab_test(g.clone(), 4);
    }

    #[test]
    fn ab_mmd_direct_full() {
        for tg in TestGraph::test_suite() {
            let graph = Csr::test_graph(tg);
            for i in 1..=10 {
                println!("testing {tg:?} with delta={i}");
                mmd_ab_test(graph.clone(), i);
            }
        }
    }
}
