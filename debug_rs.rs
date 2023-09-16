use crate::*;

/// This function checks whether or not the boundary information is correct 
#[metis_func]
pub extern "C" fn CheckBnd2(graph: *const graph_t) -> idx_t {
  let graph = graph.as_ref().unwrap();
  let nvtxs  = graph.nvtxs;
  let xadj   = graph.xadj;
  let adjncy = graph.adjncy;
  let where_  = graph.where_;
  let bndptr = graph.bndptr;
  let bndind = graph.bndind;
  let mut id;
  let mut ed;
  

for i in 0..(nvtxs as usize) {
    id = 0;
    ed = 0;
    for j in xadj[i]..xadj[i+1] {
      if (where_[i] != where_[adjncy[j]]) 
      {
        ed += graph.adjwgt[j];
      } else {
        id += graph.adjwgt[j];
      }
    }

    if (ed - id >= 0 && xadj[i] < xadj[i+1]) {
      nbnd++;

      assert_ne!(bndptr[i] != -1, "{i} {id} {ed}");
      ASSERT(bndind[bndptr[i]] == i);
    }
  }

  ASSERTP(nbnd == graph.nbnd, ("%"PRIDX" %"PRIDX"\n", nbnd, graph.nbnd));

  return 1;
}
