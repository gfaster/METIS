/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * debug.c
 *
 * This file contains code that performs self debugging
 *
 * Started 7/24/97
 * George
 *
 */

use crate::*;

/// Computes the total edgecut
#[metis_func]
pub extern "C" fn ComputeCut(graph: *const graph_t, where_: *const idx_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, cut;

    get_graph_slices!(graph => xadj adjncy);
    let nvtxs = graph.nvtxs as usize;
    mkslice!(where_, nvtxs);

    let mut cut = 0;
    if graph.adjwgt.is_null() {
        // I don't think this is ever hit
        for i in (0)..(nvtxs as usize) {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                    cut += 1;
                }
            }
        }
    } else {
        get_graph_slices!(graph => adjwgt);
        for i in (0)..(nvtxs) {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                    cut += adjwgt[j as usize];
                }
            }
        }
    }

    return cut / 2;
}

/// Computes the total edgecut with unit weights
#[metis_func]
pub extern "C" fn ComputeCutUnweighted(graph: *const graph_t, where_: *const idx_t) -> idx_t {
    let graph = graph.as_ref().unwrap();

    get_graph_slices!(graph => xadj adjncy);
    let nvtxs = graph.nvtxs as usize;
    mkslice!(where_, nvtxs);

    let mut cut = 0;
    for i in (0)..(nvtxs as usize) {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                cut += 1;
            }
        }
    }

    return cut / 2;
}

/// Computes the total volume
#[metis_func]
pub extern "C" fn ComputeVolume(graph: *const graph_t, where_: *const idx_t) -> idx_t {
    let graph = graph.as_ref().unwrap();

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_optional!(graph => vsize);
    mkslice!(where_, nvtxs);

    let nparts = where_.iter().copied().max().unwrap_or_default() as usize + 1;
    let mut marker: Vec<idx_t> = vec![-1; nparts];

    let mut totalv = 0;

    for i in 0..nvtxs {
        marker[where_[i] as usize] = i as idx_t;
        for j in xadj[i]..xadj[i + 1] {
            let k = where_[adjncy[j as usize] as usize] as usize;
            if marker[k] != i as idx_t {
                marker[k] = i as idx_t;
                totalv += if let Some(vsize) = vsize {
                    vsize[i]
                } else {
                    1
                };
            }
        }
    }

    return totalv;
}

/// (*unused*) This function computes the cut given the graph and a where_ vector
#[metis_func]
pub extern "C" fn ComputeMaxCut(graph: *mut graph_t, nparts: idx_t, where_: *mut idx_t) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // idx_t i, j, maxcut;
    // idx_t *cuts;
    let nparts = nparts as usize;
    let nvtxs = graph.nvtxs as usize;
    let mut cuts = vec![0; nparts];
    mkslice!(where_, nvtxs);
    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_optional!(graph => adjwgt);

    if let Some(adjwgt) = adjwgt {
        for i in (0)..(nvtxs) {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                    cuts[where_[i as usize] as usize] += adjwgt[j as usize];
                }
            }
        }
    } else {
        for i in (0)..(nvtxs) {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                    cuts[where_[i as usize] as usize] += 1;
                }
            }
        }
    }

    let maxcut = cuts[(util::iargmax(&cuts, 1)) as usize];

    println!("{} => {:}", util::iargmax(&cuts, 1), maxcut);

    return maxcut;
}

/// This function checks whether or not the boundary information is correct
#[metis_func]
pub extern "C" fn CheckNodeBnd(graph: *mut graph_t, onbnd: idx_t) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *xadj, *adjncy, *where_, *bndptr, *bndind;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => where_ bndptr);

    let mut nbnd = 0;
    for i in (0)..(nvtxs) {
        if where_[i as usize] == 2 {
            nbnd += 1;
        }
    }

    assert!(nbnd == onbnd, "{:} {:}", nbnd, onbnd);

    for i in (0)..(nvtxs) {
        if where_[i as usize] != 2 {
            assert!(bndptr[i as usize] == -1, "{:} {:}", i, bndptr[i as usize]);
        } else {
            assert!(bndptr[i as usize] != -1, "{:} {:}", i, bndptr[i as usize]);
        }
    }

    return 1;
}

/// Checks whether or not the rinfo of a vertex is consistent
#[metis_func]
pub extern "C" fn CheckRInfo(ctrl: *const ctrl_t, rinfo: *const ckrinfo_t) -> idx_t {
    let ctrl = ctrl.as_ref().unwrap();
    let rinfo = rinfo.as_ref().unwrap();
    // idx_t i, j;
    // cnbr_t *nbrs;

    // assert!(ctrl.nbrpoolcpos >= 0);
    assert!(rinfo.nnbrs < ctrl.nparts);

    // let nbrs = ctrl.cnbrpool + rinfo.inbr;
    let nbrs =
        std::slice::from_raw_parts(ctrl.cnbrpool.add(rinfo.inbr as usize), rinfo.nnbrs as usize);

    for i in (0)..(rinfo.nnbrs) {
        for j in (i + 1)..(rinfo.nnbrs) {
            assert_ne!(
                nbrs[i as usize].pid, nbrs[j as usize].pid,
                "{:} {:} {:} {:}",
                i,
                j,
                nbrs[i as usize].pid,
                nbrs[j as usize].pid
            );
        }
    }

    return 1;
}

#[metis_func(no_c)]
pub unsafe extern "C" fn CheckRInfoExtendedFull(ctrl: &ctrl_t, graph: &graph_t) -> idx_t {
    for i in 0..graph.nvtxs {
        CheckRInfoExtended(ctrl, graph, i);
    }
    1
}

/// asserts rinfo of a vertex is consistent and that is correct
#[metis_func(no_c)]
pub unsafe extern "C" fn CheckRInfoExtended(ctrl: &ctrl_t, graph: &graph_t, vtx: idx_t) -> idx_t {
    let vtx = vtx as usize;
    if ctrl.cnbrpool.is_null() || graph.ckrinfo.is_null() {
        return 1
    }
    let (myrinfo, mynbrs) = kwayfm::crinfos(graph.ckrinfo, ctrl.cnbrpool, vtx);
    CheckRInfo(ctrl, myrinfo);
    get_graph_slices!(graph => xadj adjncy where_);


    let mut actual_eds = vec![0; ctrl.nparts as usize];
    let this_pid = where_[vtx];
    // eprintln!("checking vtx {vtx} in pid {this_pid}");

    // check that ed and id in myrinfo is correct, and compute actual ed for each neighboring
    // partition
    let mut id = 0;
    let mut ed = 0;
    let mut actual_nnbrs = 0;
    if ctrl.objtype == METIS_OBJTYPE_CUT {
        get_graph_slices!(graph => adjwgt);
        for i in xadj[vtx]..xadj[vtx + 1] {
            let i = i as usize;
            let other = adjncy[i] as usize;
            if where_[other] == this_pid {
                id += adjwgt[i];
            } else {
                if actual_eds[where_[other] as usize] == 0 {
                    actual_nnbrs += 1;
                }
                actual_eds[where_[other] as usize] += adjwgt[i];
                ed += adjwgt[i];
            }
        }
    } else {
        for i in xadj[vtx]..xadj[vtx + 1] {
            let i = i as usize;
            let other = adjncy[i] as usize;
            if where_[other] == this_pid {
                id += 1;
            } else {
                if actual_eds[where_[other] as usize] == 0 {
                    actual_nnbrs += 1;
                }
                actual_eds[where_[other] as usize] += 1;
                ed += 1;
            }
        }
    }

    assert!(id == myrinfo.id, "internal degree of vtx {vtx} does not match: {id} != {}", myrinfo.id);
    assert!(ed == myrinfo.ed, "external degree of vtx {vtx} does not match: {ed} != {}", myrinfo.ed);
    assert_eq!(ed, myrinfo.ed);

    for nbr in mynbrs {
        let actual_ed = actual_eds[nbr.pid as usize];
        assert_ne!(actual_ed, 0, "rinfo should not contain {}", nbr.pid);
        assert_eq!(actual_ed, nbr.ed);
    }

    assert!(myrinfo.nnbrs <= xadj[vtx+1]-xadj[vtx]);
    assert_eq!(myrinfo.nnbrs, actual_nnbrs as idx_t);
    assert_eq!(mynbrs.len(), actual_nnbrs);

    1
}

/// Checks the correctness of the NodeFM data structures
#[metis_func]
pub extern "C" fn CheckNodePartitionParams(graph: *mut graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, k, l, nvtxs, me, other;
    // idx_t *xadj, *adjncy, *adjwgt, *vwgt, *where_;
    // idx_t edegrees[2 as usize], pwgts[3 as usize];

    let nvtxs = graph.nvtxs as usize;
    // xadj   = graph.xadj;
    // vwgt   = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // where_  = graph.where_;

    get_graph_slices!(graph => xadj adjncy where_ vwgt nrinfo);
    mkslice!(gpwgts: graph->pwgts, 3);

    /*------------------------------------------------------------
    / Compute now the separator external degrees
    /------------------------------------------------------------*/
    let mut pwgts = [0; 3];
    for i in (0)..(nvtxs) {
        let me = where_[i as usize];
        pwgts[me as usize] += vwgt[i as usize];

        if me == 2 {
            /* If it is on the separator do some computations */
            let mut edegrees = [0; 2];

            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let other = where_[adjncy[j as usize] as usize];
                if other != 2 {
                    edegrees[other as usize] += vwgt[adjncy[j as usize] as usize];
                }
            }
            if edegrees[0 as usize] != nrinfo[i as usize].edegrees[0 as usize]
                || edegrees[1 as usize] != nrinfo[i as usize].edegrees[1 as usize]
            {
                println!(
                    "Something wrong with edegrees: {:} {:} {:} {:} {:}",
                    i,
                    edegrees[0 as usize],
                    edegrees[1 as usize],
                    nrinfo[i as usize].edegrees[0 as usize],
                    nrinfo[i as usize].edegrees[1 as usize]
                );
                return 0;
            }
        }
    }

    if pwgts[0 as usize] != gpwgts[0 as usize]
        || pwgts[1 as usize] != gpwgts[1 as usize]
        || pwgts[2 as usize] != gpwgts[2 as usize]
    {
        println!(
            "Something wrong with part-weights: {:} {:} {:} {:} {:} {:}",
            pwgts[0 as usize],
            pwgts[1 as usize],
            pwgts[2 as usize],
            gpwgts[0 as usize],
            gpwgts[1 as usize],
            gpwgts[2 as usize]
        );
        return 0;
    }

    return 1;
}

/// Checks if the separator is indeed a separator
#[metis_func]
pub extern "C" fn IsSeparable(graph: *mut graph_t) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // idx_t i, j, nvtxs, other;
    // idx_t *xadj, *adjncy, *where_;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy where_);

    for i in (0)..(nvtxs) {
        if where_[i as usize] == 2 {
            continue;
        }
        let other = (where_[i as usize] + 1) % 2;
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            assert!(
                where_[adjncy[j as usize] as usize] != other,
                "{:} {:} {:} {:} {:} {:}",
                i,
                where_[i as usize],
                adjncy[j as usize],
                where_[adjncy[j as usize] as usize],
                xadj[(i + 1) as usize] - xadj[i as usize],
                xadj[adjncy[j as usize] as usize + 1] - xadj[adjncy[j as usize] as usize]
            );
        }
    }

    return 1;
}

/// (*only used for debug*) Recomputes the `vrinfo` fields and checks them against those in the
/// `graph.vrinfo` structure
#[metis_func(no_c)]
pub unsafe extern "C" fn CheckKWayVolPartitionParams(ctrl: *const ctrl_t, graph: *const graph_t) {
    let graph = graph.as_ref().unwrap();
    let ctrl = ctrl.as_ref().unwrap();

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vsize adjncy where_);

    let mut tmpnbrs: Vec<vnbr_t> = vec![vnbr_t::default(); ctrl.nparts as usize];

    /*------------------------------------------------------------
    / Compute now the iv/ev degrees
    /------------------------------------------------------------*/
    let mut errors = 0;
    for i in (0)..(nvtxs) {
        let me = where_[i as usize];

        let (_myrinfo, mynbrs) = kwayfm::vrinfos(graph.vkrinfo, ctrl.vnbrpool, i);

        tmpnbrs[..mynbrs.len()].clone_from_slice(mynbrs);

        // tmprinfo.nnbrs = myrinfo.nnbrs;
        // tmprinfo.nid    = myrinfo.nid;
        // tmprinfo.ned    = myrinfo.ned;
        // let tmprinfo = vkrinfo_t {
        //     gv: 0,
        //     inbr: 0,
        //     ..myrinfo.clone()
        // };

        let mynbrs = &mut tmpnbrs[..mynbrs.len()];

        for nbr in &mut *mynbrs {
            nbr.gv = 0;
        }

        for j in xadj[i]..xadj[i + 1] {
            let ii = adjncy[j as usize] as usize;
            let other = where_[ii];
            let (_orinfo, onbrs) = kwayfm::vrinfos(graph.vkrinfo, ctrl.vnbrpool, ii);
            let vsize = vsize[ii];

            if me == other {
                /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                for nbr in &mut *mynbrs {
                    if onbrs.iter().all(|onbr| onbr.pid != nbr.pid) {
                        nbr.gv -= vsize;
                    }
                }
            } else {
                /* Find the orinfo[me].ed and see if I'm the only connection */
                let Some(onbr) = onbrs.iter().find(|o| o.pid == me) else {
                    panic!("missing vnbr entry")
                };

                if onbr.ned == 1 {
                    /* I'm the only connection of 'ii' in 'me' */
                    mynbrs.iter_mut().find(|nbr| nbr.pid == other).expect("missing mynbr entry").gv += vsize;

                    /* Increase the gains for all the common domains between 'i' and 'ii' */
                    for nbr in &mut *mynbrs {
                        if nbr.pid == other {
                            continue;
                        }
                        if onbrs.iter().any(|onbr| onbr.pid == nbr.pid) {
                            nbr.gv += vsize;
                        }
                    }
                } else {
                    /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                    for nbr in &mut *mynbrs {
                        if nbr.pid == other {
                            continue
                        }

                        if onbrs.iter().all(|onbr| onbr.pid != nbr.pid) {
                            nbr.gv -= vsize;
                        }
                    }
                }
            }
        }

        let (_myrinfo, mynbrs) = kwayfm::vrinfos(graph.vkrinfo, ctrl.vnbrpool, i);

        for nbr in mynbrs {
            for tnbr in &tmpnbrs[..mynbrs.len()] {
                if tnbr.pid == nbr.pid {
                    if tnbr.gv != nbr.gv {
                        errors += 1;
                        println!(
                            "[{:8} {:8} {:8} {:+8} {:+8}]",
                            i,
                            where_[i],
                            nbr.pid,
                            nbr.gv,
                            tnbr.gv
                        );
                    }
                    break;
                }
            }
        }
    }

    if errors != 0 {
        panic!("found {errors} errors")
    }
}

#[metis_func(no_c)]
pub unsafe extern "C" fn PrintRInfoVol(ctrl: *const ctrl_t, graph: *const graph_t, vtx: idx_t) {
    let ctrl = ctrl.as_ref().unwrap();
    let graph = graph.as_ref().unwrap();
    let vtx = vtx as usize;

    get_graph_slices!(graph => xadj adjncy where_ vsize);

    let (myrinfo, mynbrs) = kwayfm::vrinfos(graph.vkrinfo, ctrl.vnbrpool, vtx);

    let pid = where_[vtx];

    println!("\nRinfo of {vtx} [{pid}] [nid,ned:{:3},{:3}] [gv:{:8}]:", myrinfo.nid, myrinfo.ned, myrinfo.gv);
    for nbr in mynbrs {
        let vnbr_t { pid: nbr_pid, ned, gv } = nbr;
        println!("\t{nbr_pid:3} [ned:{ned:3}] [gv:{gv:8}]");

        println!("\t\texternal:");
        for i in xadj[vtx]..xadj[vtx + 1] {
            let i = adjncy[i as usize] as usize;
            if where_[i] != *nbr_pid {
                continue
            }
            for j in xadj[i]..xadj[i + 1] {
                let j = adjncy[j as usize] as usize;
                if where_[j] != pid {
                    println!("\t\t[{i:5}->{j:5}] [{nbr_pid:3} => {opid:3}] ({vsizei:5} <=> {vsizej:5})", vsizei=vsize[i], vsizej=vsize[j], opid=where_[j]);
                }
            }
        }

        println!("\n\t\tinternal (where == {pid}):");
        for i in xadj[vtx]..xadj[vtx + 1] {
            let i = adjncy[i as usize] as usize;
            if where_[i] != *nbr_pid {
                continue
            }
            for j in xadj[i]..xadj[i + 1] {
                let j = adjncy[j as usize] as usize;
                if where_[j] == pid {
                    println!("\t\t[{i:5}->{j:5}] ({vsizei:5} <=> {vsizej:5})", vsizei=vsize[i], vsizej=vsize[j]);
                }
            }
        }
        println!();
    }
}

#[metis_func]
pub extern "C" fn CheckBnd(graph: *const graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *xadj, *adjncy, *where, *bndptr, *bndind;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy where_ bndptr bndind);

    let mut nbnd = 0;
    for i in 0..nvtxs {
        if xadj[i + 1] - xadj[i] == 0 {
            nbnd += 1; /* Islands are considered to be boundary vertices */
        }
        for j in xadj[i]..xadj[i + 1] {
            if where_[i] != where_[adjncy[j as usize] as usize] {
                nbnd += 1;
                assert!(bndptr[i] != -1);
                assert!(bndind[bndptr[i] as usize] == i as idx_t);
                break;
            }
        }
    }

    assert!(
        nbnd == graph.nbnd,
        "calculated boundry size vs actual: {} vs {}",
        nbnd,
        graph.nbnd
    );

    1
}

/// This function checks whether or not the boundary information is correct
#[metis_func]
pub extern "C" fn CheckBnd2(graph: *const graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    let nvtxs = graph.nvtxs;
    get_graph_slices!(graph => xadj adjncy where_ bndptr bndind adjwgt);
    let mut id;
    let mut ed;

    let mut nbnd = 0;
    for i in 0..(nvtxs as usize) {
        id = 0;
        ed = 0;
        for j in xadj[i]..xadj[i + 1] {
            let j = j as usize;
            if where_[i] != where_[adjncy[j] as usize] {
                ed += adjwgt[j];
            } else {
                id += adjwgt[j];
            }
        }

        if ed - id >= 0 && xadj[i] < xadj[i + 1] {
            nbnd += 1;

            assert_ne!(bndptr[i], -1, "{i} {id} {ed}");
            assert_eq!(bndind[bndptr[i] as usize], i as idx_t);
        }
    }

    assert_eq!(nbnd, graph.nbnd);

    1
}

#[track_caller]
pub unsafe fn check_adj(graph: &graph_t) -> bool {
    get_graph_slices!(graph => xadj adjncy);

    let nvtxs = graph.nvtxs as usize;

    assert_eq!(adjncy.len(), xadj[nvtxs] as usize);
    assert_eq!(adjncy.len(), graph.nedges as usize);

    let mut seen = vec![false; nvtxs];
    let mut seenlist = vec![];
    for (vtx, [istart, iend]) in xadj.windows(2).map(|a| [a[0], a[1]]).enumerate() {
        assert!(
            istart <= iend,
            "vtx {vtx} edges start after end ({istart} > {iend})"
        );
        assert!(istart >= 0, "vtx {vtx} start is {istart}");
        assert!(
            iend <= adjncy.len() as idx_t,
            "vtx {vtx} end ({iend}) is out of adjncy ({})",
            adjncy.len()
        );
        let adjrng = (istart as usize)..(iend as usize);
        for &adj in &adjncy[adjrng] {
            assert_ne!(vtx, adj as usize, "vtx {vtx} has self loop");
            assert!(
                adj >= 0 && (adj as usize) < nvtxs,
                "vtx {vtx} has out bounds edge to {adj} (max: {nvtxs})"
            );
            assert!(!seen[adj as usize], "vtx {vtx} has duplicate edge to {adj}");
            seen[adj as usize] = true;
            seenlist.push(adj);
        }
        for i in seenlist.drain(..) {
            seen[i as usize] = false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use crate::tests::ab_test_partition_test_graphs;

    use super::*;

    #[test]
    #[cfg_attr(not(debug_assertions), ignore = "requires debug assertions")]
    fn ab_ComputeCut() {
        ab_test_partition_test_graphs("ComputeCut:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.set_objective(Objtype::Cut);
            g.set_seed(1234);
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    #[cfg_attr(not(debug_assertions), ignore = "requires debug assertions")]
    fn ab_ComputeCut_MC() {
        ab_test_partition_test_graphs("ComputeCut:rs", Optype::Kmetis, 20, 3, |mut g| {
            g.set_objective(Objtype::Cut);
            g.set_seed(1234);
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    #[cfg_attr(not(debug_assertions), ignore = "requires debug assertions")]
    fn ab_ComputeVolume() {
        ab_test_partition_test_graphs("ComputeVolume:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.set_objective(Objtype::Vol);
            g.set_seed(1234);
            g.random_vsize();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    #[cfg_attr(not(debug_assertions), ignore = "requires debug assertions")]
    fn ab_ComputeVolume_MC() {
        ab_test_partition_test_graphs("ComputeVolume:rs", Optype::Kmetis, 20, 3, |mut g| {
            g.set_objective(Objtype::Vol);
            g.set_seed(1234);
            g.random_vsize();
            g.random_tpwgts();
            g
        });
    }
}
