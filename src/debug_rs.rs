//! partial port of debug.c

use crate::*;

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
        assert!(istart <= iend, "vtx {vtx} edges start after end ({istart} > {iend})");
        assert!(istart >= 0, "vtx {vtx} start is {istart}");
        assert!(iend <= adjncy.len() as idx_t, "vtx {vtx} end ({iend}) is out of adjncy ({})", adjncy.len());
        let adjrng = (istart as usize)..(iend as usize);
        for &adj in &adjncy[adjrng] {
            assert_ne!(vtx, adj as usize, "vtx {vtx} has self loop");
            assert!(adj >= 0 && (adj as usize) < nvtxs, "vtx {vtx} has out bounds edge to {adj} (max: {nvtxs})");
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
