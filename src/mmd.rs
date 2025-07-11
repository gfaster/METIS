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
    // idx_t  ehead, i, mdeg, mdlmt, mdeg_node, nextmd, num, tag;

    if (neqns <= 0) {
        return;
    }

    /* adjust from C to Fortran */
    // FIXME: this is horrible UB. It doesn't seem like any compiler currently
    // destroys this but I want to fix it
    // https://play.rust-lang.org/?version=stable&mode=debug&edition=2024&gist=34331f6a82be75de727ef2479302c7dd
    todo!("xadj--; adjncy--; invp--; perm--; head--; qsize--; list--; marker--;")

    /* initialization for the minimum degree algorithm */
    *ncsub = 0;
    mmdint(neqns, xadj, adjncy, head, invp, perm, qsize, list, marker);

    /* 'num' counts the number of ordered nodes plus 1 */
    num = 1;

    /* eliminate all isolated nodes */
    nextmd = head[1 as usize];
    while (nextmd > 0) {
        mdeg_node = nextmd;
        nextmd = invp[mdeg_node as usize];
        marker[mdeg_node as usize] = maxint;
        invp[mdeg_node as usize] = -num;
        num += 1;
    }

    /* search for node of the minimum degree. 'mdeg' is the current */
    /* minimum degree; 'tag' is used to facilitate marking nodes.   */
    if (num > neqns) {
        mmdnum(neqns, perm, invp, qsize);
        return;
    }

    tag = 1;
    head[1 as usize] = 0;
    mdeg = 2;

    /* infinite loop here */
    'outer: loop {
        while (head[mdeg as usize] <= 0) {
            mdeg += 1;
        }

        /* use value of 'delta' to set up 'mdlmt', which governs */
        /* when a degree update is to be performed.              */
        //mdlmt = mdeg + delta;
        // the need for gk_min() was identified by jsf67
        // https://github.com/KarypisLab/METIS/issues/46
        mdlmt = gk_min(neqns, mdeg + delta);
        ehead = 0;

        // merged n500 and n900
        // goto n500 => continue n500
        // goto n900 => break n500
        'n500: loop {
            mdeg_node = head[mdeg as usize];
            while (mdeg_node <= 0) {
                mdeg += 1;

                if (mdeg > mdlmt) {
                    break 'n500;
                }
                mdeg_node = head[mdeg as usize];
            }

            /* remove 'mdeg_node' from the degree structure */
            nextmd = invp[mdeg_node as usize];
            head[mdeg as usize] = nextmd;
            if (nextmd > 0) {
                perm[nextmd as usize] = -mdeg;
            }
            invp[mdeg_node as usize] = -num;
            *ncsub += mdeg + qsize[mdeg_node as usize] - 2;
            if ((num + qsize[mdeg_node as usize]) > neqns) {
                break 'outer;
            }

            /*  eliminate 'mdeg_node' and perform quotient graph */
            /*  transformation. reset 'tag' value if necessary.    */
            tag += 1;
            if (tag >= maxint) {
                tag = 1;
                for i in (1)..=(neqns) {
                    if (marker[i as usize] < maxint) {
                        marker[i as usize] = 0;
                    }
                }
            };

            mmdelm(
                mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, maxint, tag,
            );

            num += qsize[mdeg_node as usize];
            list[mdeg_node as usize] = ehead;
            ehead = mdeg_node;
            if (delta >= 0) {
                continue 'n500;
            } else {
                break 'n500;
            }
        }

        // was merged with n500
        // n900: 

        /* update degrees of the nodes involved in the  */
        /* minimum degree nodes elimination.            */
        if (num > neqns) {
            break 'outer;
        }

        mmdupd(
            ehead, neqns, xadj, adjncy, delta, &mdeg, head, invp, perm, qsize, list, marker,
            maxint, &tag,
        );
    }

    mmdnum(neqns, perm, invp, qsize);

    /* Adjust from Fortran back to C*/
    // xadj; adjncy++; invp++; perm++; head++; qsize++; list++; marker += 1;xadj+=1;
}

/**************************************************************************
*           mmdelm ...... multiple minimum degree elimination
* Purpose -- This routine eliminates the node mdeg_node of minimum degree
*     from the adjacency structure, which is stored in the quotient
*     graph format. It also transforms the quotient graph representation
*     of the elimination graph.
* Input parameters --
*     mdeg_node -- node of minimum degree.
*     maxint -- estimate of maximum representable (short) integer.
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
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    head: *mut idx_t,
    forward: *mut idx_t,
    backward: *mut idx_t,
    qsize: *mut idx_t,
    list: *mut idx_t,
    marker: *mut idx_t,
    maxint: idx_t,
    tag: idx_t,
) {
    // idx_t   element, i,   istop, istart, j,
    //       jstop, jstart, link,
    //       nabor, node, npv, nqnbrs, nxnode,
    //       pvnode, rlmt, rloc, rnode, xqnbr;

    /* find the reachable set of 'mdeg_node' and */
    /* place it in the data structure.           */
    marker[mdeg_node as usize] = tag;
    istart = xadj[mdeg_node as usize];
    istop = xadj[(mdeg_node + 1) as usize] - 1;

    /* 'element' points to the beginning of the list of  */
    /* eliminated nabors of 'mdeg_node', and 'rloc' gives the */
    /* storage location for the next reachable node.   */
    element = 0;
    rloc = istart;
    rlmt = istop;
    for i in istart..=istop {
        nabor = adjncy[i as usize];
        if (nabor == 0) {
            break;
        }
        if (marker[nabor as usize] < tag) {
            marker[nabor as usize] = tag;
            if (forward[nabor as usize] < 0) {
                list[nabor as usize] = element;
                element = nabor;
            } else {
                adjncy[rloc as usize] = nabor;
                rloc += 1;
            };
        }; /* end of -- if -- */
    } /* end of -- for -- */

    /* merge with reachable nodes from generalized elements. */
    while (element > 0) {
        adjncy[rlmt as usize] = -element;
        link = element;

        todo!("n400:"); // small loop
        jstart = xadj[link as usize];
        jstop = xadj[(link + 1) as usize] - 1;
        for j in jstart..=jstop {
            node = adjncy[j as usize];
            link = -node;
            if (node < 0) {
                todo!("goto n400");
            }
            if (node == 0) {
                break;
            }
            if (marker[node as usize] < tag && forward[node as usize] >= 0) {
                marker[node as usize] = tag;
                /*use storage from eliminated nodes if necessary.*/
                while (rloc >= rlmt) {
                    link = -adjncy[rlmt as usize];
                    rloc = xadj[link as usize];
                    rlmt = xadj[(link + 1) as usize] - 1;
                }
                adjncy[rloc as usize] = node;
                rloc += 1;
            };
        } /* end of -- for ( j = jstart; -- */

        element = list[element as usize];
    } /* end of -- while ( element > 0 ) -- */

    if (rloc <= rlmt) {
        adjncy[rloc as usize] = 0;
    }
    /* for each node in the reachable set, do the following. */
    link = mdeg_node;

    todo!("n1100:");
    istart = xadj[link as usize];
    istop = xadj[(link + 1) as usize] - 1;
    for i in istart..=istop {
        rnode = adjncy[i as usize];
        link = -rnode;
        if (rnode < 0) {
            todo!("goto n1100");
        }
        if (rnode == 0) {
            return;
        }

        /* 'rnode' is in the degree list structure. */
        pvnode = backward[rnode as usize];
        if ((pvnode != 0) && (pvnode != (-maxint))) {
            /* then remove 'rnode' from the structure. */
            nxnode = forward[rnode as usize];
            if (nxnode > 0) {
                backward[nxnode as usize] = pvnode;
            }
            if (pvnode > 0) {
                forward[pvnode as usize] = nxnode;
            }
            npv = -pvnode;
            if (pvnode < 0) {
                head[npv as usize] = nxnode;
            };
        }

        /* purge inactive quotient nabors of 'rnode'. */
        jstart = xadj[rnode as usize];
        jstop = xadj[(rnode + 1) as usize] - 1;
        xqnbr = jstart;
        for j in jstart..=jstop {
            nabor = adjncy[j as usize];
            if (nabor == 0) {
                break;
            }
            if (marker[nabor as usize] < tag) {
                adjncy[xqnbr as usize] = nabor;
                xqnbr += 1;
            };
        }

        /* no active nabor after the purging. */
        nqnbrs = xqnbr - jstart;
        if (nqnbrs <= 0) {
            /* merge 'rnode' with 'mdeg_node'. */
            qsize[mdeg_node as usize] += qsize[rnode as usize];
            qsize[rnode as usize] = 0;
            marker[rnode as usize] = maxint;
            forward[rnode as usize] = -mdeg_node;
            backward[rnode as usize] = -maxint;
        } else {
            /* flag 'rnode' for degree update, and  */
            /* add 'mdeg_node' as a nabor of 'rnode'.      */
            forward[rnode as usize] = nqnbrs + 1;
            backward[rnode as usize] = 0;
            adjncy[xqnbr as usize] = mdeg_node;
            xqnbr += 1;
            if (xqnbr <= jstop) {
                adjncy[xqnbr as usize] = 0;
            }
        };
    } /* end of -- for ( i = istart; -- */
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
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    head: *mut idx_t,
    forward: *mut idx_t,
    backward: *const idx_t,
    qsize: *const idx_t,
    list: *const idx_t,
    marker: *const idx_t,
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
        ndeg = xadj[(node + 1) as usize] - xadj[node as usize] + 1;
        fnode = head[ndeg as usize];
        forward[node as usize] = fnode;
        head[ndeg as usize] = node;
        if (fnode > 0) {
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
fn mmdnum(neqns: idx_t, perm: *mut idx_t, invp: *mut idx_t, qsize: *mut idx_t) {
    // idx_t father, nextf, node, nqsize, num, root;

    // for ( node = 1; node <= neqns; node ) {
    for node in 1..=neqns {
        nqsize = qsize[node as usize];
        if (nqsize <= 0) {
            perm[node as usize] = invp[node as usize];
        }
        if (nqsize > 0) {
            perm[node as usize] = -invp[node as usize];
        }
    }

    /* for each node which has been merged, do the following. */
    for node in 1..=neqns {
        if (perm[node as usize] <= 0) {
            /* trace the merged tree until one which has not */
            /* been merged, call it root.                    */
            father = node;
            while (perm[father as usize] <= 0) {
                father = -perm[father as usize];
            }

            /* number node after root. */
            root = father;
            num = perm[root as usize] + 1;
            invp[node as usize] = -num;
            perm[root as usize] = num;

            /* shorten the merged tree. */
            father = node;
            nextf = -perm[father as usize];
            while (nextf > 0) {
                perm[father as usize] = -root;
                father = nextf;
                nextf = -perm[father as usize];
            }
        }; /* end of -- if ( perm[node as usize] <= 0 ) -- */
    } /* end of -- for ( node = 1; -- */

    /* ready to compute perm. */
    for node in 1..=neqns {
        num = -invp[node as usize];
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
*    maxint -- maximum machine representable (short) integer.
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
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    delta: idx_t,
    mdeg: *mut idx_t,
    head: *const idx_t,
    forward: *const idx_t,
    backward: *const idx_t,
    qsize: *const idx_t,
    list: *const idx_t,
    marker: *const idx_t,
    maxint: idx_t,
    tag: *const idx_t,
) {
    // idx_t  deg, deg0, element, enode, fnode, i, iq2, istop,
    //      istart, j, jstop, jstart, link, mdeg0, mtag, nabor,
    //      node, q2head, qxhead;
    //
    //      mdeg0 = *mdeg + delta;
    //      element = ehead;

    'n100: loop {
        if (element <= 0) {
            return;
        }

        /* for each of the newly formed element, do the following. */
        /* reset tag value if necessary.                           */
        mtag = *tag + mdeg0;
        if (mtag >= maxint) {
            *tag = 1;
            for i in 1..=neqns {
                if (marker[i as usize] < maxint) {
                    marker[i as usize] = 0;
                }
            }
            mtag = *tag + mdeg0;
        };

        /* create two linked lists from nodes associated with 'element': */
        /* one with two nabors (q2head) in the adjacency structure, and the*/
        /* other with more than two nabors (qxhead). also compute 'deg0',*/
        /* number of nodes in this element.                              */
        q2head = 0;
        qxhead = 0;
        deg0 = 0;
        link = element;

        todo!("n400:");
        istart = xadj[link as usize];
        istop = xadj[(link + 1) as usize] - 1;
        for i in istart..=istop {
            enode = adjncy[i as usize];
            link = -enode;
            if (enode < 0) {
                todo!("goto n400");
            }
            if (enode == 0) {
                break;
            }
            if (qsize[enode as usize] != 0) {
                deg0 += qsize[enode as usize];
                marker[enode as usize] = mtag;

                /*'enode' requires a degree update*/
                if (backward[enode as usize] == 0) {
                    /* place either in qxhead or q2head list. */
                    if (forward[enode as usize] != 2) {
                        list[enode as usize] = qxhead;
                        qxhead = enode;
                    } else {
                        list[enode as usize] = q2head;
                        q2head = enode;
                    };
                };
            }; /* enf of -- if ( qsize[enode as usize] != 0 ) -- */
        } /* end of -- for ( i = istart; -- */

        /* for each node in q2 list, do the following. */
        enode = q2head;
        iq2 = 1;

        todo!("n900:"); // big loop (outer)
        if (enode <= 0) {
            todo!("goto n1500");
        }

        if (backward[enode as usize] != 0) {
            todo!("goto n2200");
        }

        (*tag) += 1;
        deg = deg0;

        /* identify the other adjacent element nabor. */
        istart = xadj[enode as usize];
        nabor = adjncy[istart as usize];
        if (nabor == element) {
            nabor = adjncy[(istart + 1) as usize];
        }
        link = nabor;
        if (forward[nabor as usize] >= 0) {
            /* nabor is uneliminated, increase degree count. */
            deg += qsize[nabor as usize];
            todo!("goto n2100");
        };

        /* the nabor is eliminated. for each node in the 2nd element */
        /* do the following.                                         */
        todo!("n1000:");
        istart = xadj[link as usize];
        istop = xadj[(link + 1) as usize] - 1;
        for i in istart..=istop {
            node = adjncy[i as usize];
            link = -node;
            if (node != enode) {
                if (node < 0) {
                    todo!("goto n1000");
                }
                if (node == 0) {
                    todo!("goto n2100");
                }
                if (qsize[node as usize] != 0) {
                    if (marker[node as usize] < *tag) {
                        /* 'node' is not yet considered. */
                        marker[node as usize] = *tag;
                        deg += qsize[node as usize];
                    } else {
                        if (backward[node as usize] == 0) {
                            if (forward[node as usize] == 2) {
                                /* 'node' is indistinguishable from 'enode'.*/
                                /* merge them into a new supernode.         */
                                qsize[enode as usize] += qsize[node as usize];
                                qsize[node as usize] = 0;
                                marker[node as usize] = maxint;
                                forward[node as usize] = -enode;
                                backward[node as usize] = -maxint;
                            } else {
                                /* 'node' is outmacthed by 'enode' */
                                if (backward[node as usize] == 0) {
                                    backward[node as usize] = -maxint;
                                }
                            };
                        }; /* end of -- if ( backward[node as usize] == 0 ) -- */
                    }; /* end of -- if ( marker[node as usize] < *tag ) -- */
                }; /* end of -- if ( qsize[node as usize] != 0 ) -- */
            }; /* end of -- if ( node != enode ) -- */
        } /* end of -- for ( i = istart; -- */
        todo!("goto n2100");

        todo!("n1500:");
        /* for each 'enode' in the 'qx' list, do the following. */
        enode = qxhead;
        iq2 = 0;

        todo!("n1600:");
        if (enode <= 0) {
            todo!("goto n2300");
        }
        if (backward[enode as usize] != 0) {
            todo!("goto n2200");
        }
        (*tag) += 1;
        deg = deg0;

        /*for each unmarked nabor of 'enode', do the following.*/
        istart = xadj[enode as usize];
        istop = xadj[(enode + 1) as usize] - 1;
        for i in istart..=istop {
            nabor = adjncy[i as usize];
            if (nabor == 0) {
                break;
            }
            if (marker[nabor as usize] < *tag) {
                marker[nabor as usize] = *tag;
                link = nabor;
                if (forward[nabor as usize] >= 0) {
                    /*if uneliminated, include it in deg count.*/
                    deg += qsize[nabor as usize];
                } else {
                    todo!("n1700:"); // smaller loop

                    /* if eliminated, include unmarked nodes in this*/
                    /* element into the degree count.             */
                    jstart = xadj[link as usize];
                    jstop = xadj[(link + 1) as usize] - 1;
                    for j in jstart..=jstop {
                        j += 1;
                        node = adjncy[j as usize];
                        link = -node;
                        if (node < 0) {
                            todo!("goto n1700");
                        }
                        if (node == 0) {
                            break;
                        }
                        if (marker[node as usize] < *tag) {
                            marker[node as usize] = *tag;
                            deg += qsize[node as usize];
                        };
                    } /* end of -- for ( j = jstart; -- */
                }; /* end of -- if ( forward[nabor as usize] >= 0 ) -- */
            }; /* end of -- if ( marker[nabor as usize] < *tag ) -- */
        } /* end of -- for ( i = istart; -- */

        todo!("n2100:");
        /* update external degree of 'enode' in degree structure, */
        /* and '*mdeg' if necessary.                     */
        deg = deg - qsize[enode as usize] + 1;
        fnode = head[deg as usize];
        forward[enode as usize] = fnode;
        backward[enode as usize] = -deg;
        if (fnode > 0) {
            backward[fnode as usize] = enode;
        }
        head[deg as usize] = enode;
        if (deg < *mdeg) {
            *mdeg = deg;
        }

        todo!("n2200:");
        /* get next enode in current element. */
        enode = list[enode as usize];
        if (iq2 == 1) {
            todo!("goto n900");
        }

        todo!("goto n1600");

        todo!("n2300:"); // epilogue of n100 loop - only way to continue

        /* get next element in the list. */
        *tag = mtag;
        element = list[element as usize];
    }
}
