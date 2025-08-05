/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * separator.c
 *
 * This file contains code for separator extraction
 *
 * Started 8/1/97
 * George
 *
 * $Id: separator.c 10481 2011-07-05 18:01:23Z karypis $
 *
 */

use crate::*;

/*************************************************************************
* This function takes a bisection and constructs a minimum weight vertex 
* separator out of it. It uses the node-based separator refinement for it.
**************************************************************************/
#[metis_func]
pub extern "C" fn ConstructSeparator(ctrl: *mut ctrl_t, graph: *mut graph_t) 
{
  let graph = graph.as_mut().unwrap();
  let ctrl = ctrl.as_mut().unwrap();

  // WCOREPUSH;

  get_graph_slices!(graph => xadj bndind);

  let mut where_ = Vec::from(get_graph_slice!(graph => where_));

  /* Put the nodes in the boundary into the separator */
  for i in (0)..(graph.nbnd) {
    let j = bndind[i as usize] as usize;
    // ignore islands
    if xadj[j+1]-xadj[j] <= 0 {
      continue
    }

    where_[j] = 2;
  }

  graph::FreeRData(graph);

  srefine::Allocate2WayNodePartitionMemory(ctrl, graph);
  get_graph_slice_mut!(graph => where_).copy_from_slice(&where_);

  // WCOREPOP;

  debug_assert!(debug::IsSeparable(graph) != 0);

  srefine::Compute2WayNodePartitionParams(ctrl, graph);

  debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

  sfm::FM_2WayNodeRefine2Sided(ctrl, graph, 1); 
  sfm::FM_2WayNodeRefine1Sided(ctrl, graph, 4); 

  debug_assert!(debug::IsSeparable(graph) != 0);

}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use crate::tests::ab_test_partition_test_graphs;

    #[test]
    fn ab_ConstructSeparator() {
        ab_test_partition_test_graphs("ConstructSeparator:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }
}
