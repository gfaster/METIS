use std::sync::Arc;

use fastrand::Rng;

use crate::{costs::{CaseCost, CostFactors}, graph::Graph, utils::{make_csr, RetainIndexed}};

use super::{Case, Strategy};


pub struct DeleteVtxs {
    /// inverse probability a vertex will be deleted (i.e. `prob = 3` implies about a third of
    /// vertices will be deleted)
    pub prob: usize,
    /// random number generator used for choosing vertices
    pub rng: Rng, 
    /// whether nparts should be reduced proportionally to the number of vertices deleted
    pub adjust_nparts: bool,
    pub allow_backtrack: bool,
}

impl Strategy for DeleteVtxs {
    fn cost(&self, case: &Case) -> CaseCost {
        let nvtxs = case.graph.nvtxs();
        let nedges = case.graph.nedges();
        let num_del = nvtxs / self.prob;
        let avg_deg = nedges / nvtxs;

        let new_nvtxs = nvtxs - num_del;
        let new_nedges = nedges - (avg_deg * num_del);
        let new_nparts = if self.adjust_nparts { adjust_nparts(case.nparts, nvtxs, new_nvtxs) } else { case.nparts };
        CostFactors {
            nvtxs: new_nvtxs as u32,
            nedges: new_nedges as u32,
            nparts: new_nparts as u32,
            ..CostFactors::from_case(case)
        }.into()
    }

    fn dont_backtrack(&self) -> bool {
        !self.allow_backtrack
    }

    fn is_valid(&self, case: &Case) -> bool {
        self.prob >= 1 && self.prob <= case.graph.nvtxs()
    }

    fn apply(&self, case: Case) -> Option<impl super::StrategyIter> {
        if !self.is_valid(&case) {
            return None
        }
        let mut rng = self.rng.clone();
        let nvtxs = case.graph.nvtxs();
        Some(std::iter::from_fn(move || {
            let mut case = case.clone();
            let mut rng = rng.fork();

            let to_delete: Vec<usize> = if self.prob * 16 < nvtxs {
                // will likely have few enough collisions that it'll just be best to generate a
                // bunch of random vertices
                std::iter::repeat_with(|| rng.usize(..nvtxs)).take(nvtxs.div_ceil(self.prob)).collect()
            } else {
                rng.choose_multiple(0..nvtxs, nvtxs.div_ceil(self.prob))
            };

            // it's fine if we end iteration if we end up with an empty graph -- It should be rare
            // and it implies a very sparse graph along with high probability
            case.graph = delete_vtxs(case.graph, to_delete)?;

            if self.adjust_nparts {
                case.nparts = adjust_nparts(case.nparts, nvtxs, case.graph.nvtxs())
            }

            Some(case)
        }))
    }
}

fn adjust_nparts(old_nparts: usize, old_nvtxs: usize, new_nvtxs: usize) -> usize {
    let part_size = old_nvtxs / old_nparts;
    let nparts = new_nvtxs / part_size;
    nparts.clamp(2, old_nparts)
}

/// construct a vec that maps the old vertices to the new ones, or [`usize::MAX`] if it's deleted
fn construct_deletion_mapping(nvtxs: usize, to_delete: &[usize]) -> Vec<usize> {
    let mut m = Vec::with_capacity(nvtxs);

    let mut to = 0;
    let mut idel = 0;
    for from in 0..nvtxs {
        let del = *to_delete.get(idel).unwrap_or(&usize::MAX);
        if from < del {
            m.push(to);
            to += 1;
        } else if from == del {
            m.push(usize::MAX);
            idel += 1;
        } else {
            unreachable!("to_delete isn't sorted")
        }
    }
    debug_assert_eq!(to_delete.len(), m.iter().filter(|&&v| v == usize::MAX).count());
    m
}

fn delete_vtxs(mut graph: Graph, mut to_delete: Vec<usize>) -> Option<Graph> {
    if to_delete.is_empty() {
        return None
    }
    to_delete.sort_unstable();
    to_delete.dedup();
    let nvtxs = graph.nvtxs();
    let Graph { xadj, adjncy, ncon, adjwgt, vwgt, vsize } = &mut graph;
    let ncon = *ncon;
    let xadj = Arc::make_mut(xadj);
    let adjncy = Arc::make_mut(adjncy);
    let adjwgt = adjwgt.as_mut().map(Arc::make_mut);
    let vwgt = vwgt.as_mut().map(Arc::make_mut);
    let vsize = vsize.as_mut().map(Arc::make_mut);

    let initial_mapping = construct_deletion_mapping(nvtxs, &to_delete);

    // mark deleted edges, but don't update their indices yet since we want to also delete vertices
    // that have their degree reduced to zero
    {
        // mark edges going outward
        for edge in &mut *adjncy {
            if initial_mapping[*edge] == usize::MAX {
                *edge = usize::MAX;
            }
        }
        // mark edges coming from deleted vertices
        for &vtx in &to_delete {
            adjncy[xadj[vtx]..xadj[vtx + 1]].fill(usize::MAX);
        }
    }

    let mapping = {
        // we want to also delete vertices that were orphaned
        let to_delete_more: Vec<usize> = (0..nvtxs).filter(|&vtx| {
            initial_mapping[vtx] != usize::MAX
            && adjncy[xadj[vtx]..xadj[vtx+1]].iter().all(|&e| e == usize::MAX)
        }).collect();

        // need to reconstruct mapping only if we found orphaned vertices
        if to_delete_more.is_empty() {
            initial_mapping
        } else {
            to_delete.extend(to_delete_more);
            to_delete.sort_unstable();
            debug_assert!(to_delete.windows(2).all(|w| w[0] < w[1]),
                "additional vertices added for deletion should not have been in the original set");
            drop(initial_mapping);
            construct_deletion_mapping(nvtxs, &to_delete)
        }
    };

    let new_nvtxs = nvtxs - to_delete.len();

    if new_nvtxs < 2 {
        return None
    }

    // update edges touching deleted vertices -- this shouldn't mark any new edges
    for edge in &mut *adjncy {
        // if *edge != usize::MAX || mapping[*edge] == usize::MAX {
        // }
        if *edge != usize::MAX {
            debug_assert_ne!(mapping[*edge], usize::MAX,
                "second edge pass shouldn't mark any additional edges for deletion");
            *edge = mapping[*edge];
        }
    }

    // apply deletion of adjwgts
    if let Some(adjwgt) = adjwgt {
        adjwgt.retain_indexed(|i, _| adjncy[i] != usize::MAX);
    }

    // apply deletion of vwgts -- I need to double check this is correct
    if let Some(vwgt) = vwgt {
        vwgt.retain_indexed(|i, _| mapping[i/ncon] != usize::MAX);
    }

    // apply deletion of vsize
    if let Some(vsize) = vsize {
        vsize.retain_indexed(|i, _| mapping[i] != usize::MAX);
    }

    // update xadj by first modifying 0..nvtxs to be the new degree and then constructing a CSR
    // from there.
    {
        for (iw, i) in (0..nvtxs).filter(|&i| mapping[i] != usize::MAX).enumerate() {
            let deg = adjncy[xadj[i]..xadj[i+1]].iter().filter(|&&e| e != usize::MAX).count();
            xadj[iw] = deg;
        }
        make_csr(&mut xadj[..=new_nvtxs]);
        xadj.resize(new_nvtxs + 1, usize::MAX);
    }

    // now we can finally delete the marked edges in adjncy
    adjncy.retain(|&e| e != usize::MAX);

    Some(graph)
}
