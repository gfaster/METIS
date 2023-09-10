/*
\file
\brief Driving routines for multilevel k-way refinement

\date   Started 7/28/1997
\author George 
\author  Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: kwayrefine.c 20398 2016-11-22 17:17:12Z karypis $ 
*/

use std::ffi::c_void;
use crate::*;


/*************************************************************************/
/* This function is the entry point of cut-based refinement */

/*************************************************************************/
#[metis_func]
pub extern "C" fn RefineKWay(ctrl: *mut ctrl_t, orggraph: *mut graph_t, graph: *mut graph_t) -> () {
  let mut i;
  let mut nlevels;
  let contig = ctrl.contig;
  let mut ptr: *mut graph_t;

  /* Determine how many levels are there */
  ptr = graph;
  nlevels = 0;
  while ptr != orggraph {
    nlevels += 1;
    ptr = ptr.finer;
  }

  /* Compute the parameters of the coarsest graph */
  ComputeKWayPartitionParams(ctrl, graph);

  /* Try to minimize the sub-domain connectivity */
  if ctrl.minconn != 0 {
    EliminateSubDomainEdges(ctrl, graph);
  }
  
  /* Deal with contiguity constraints at the beginning */
  if contig != 0 && FindPartitionInducedComponents(graph, graph.where_, std::ptr::null(), std::ptr::null()) > ctrl.nparts { 
    EliminateComponents(ctrl, graph);

    ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
    Greedy_KWayOptimize(ctrl, graph, 5, 0, OMODE_BALANCE); 

    ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
    Greedy_KWayOptimize(ctrl, graph, ctrl.niter, 0, OMODE_REFINE); 

    ctrl.contig = 0;
  }

  /* Refine each successively finer graph */
  for i in 0.. {
    if ctrl.minconn && i == nlevels/2 {
      EliminateSubDomainEdges(ctrl, graph);
    }

    // IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.RefTmr));

    if (2*i >= nlevels && !IsBalanced(ctrl, graph, 0.02)) {
      ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
      Greedy_KWayOptimize(ctrl, graph, 1, 0, OMODE_BALANCE); 
      ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
    }

    Greedy_KWayOptimize(ctrl, graph, ctrl.niter, 5.0, OMODE_REFINE); 

    // IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.RefTmr));

    /* Deal with contiguity constraints in the middle */
    if (contig && i == nlevels/2) {
      if (FindPartitionInducedComponents(graph, graph.where_, ptr::null_mut(), ptr::null_mut()) > ctrl.nparts) {
        EliminateComponents(ctrl, graph);

        if (!IsBalanced(ctrl, graph, 0.02)) {
          ctrl.contig = 1;
          ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
          Greedy_KWayOptimize(ctrl, graph, 5, 0, OMODE_BALANCE); 
  
          ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
          Greedy_KWayOptimize(ctrl, graph, ctrl.niter, 0, OMODE_REFINE); 
          ctrl.contig = 0;
        }
      }
    }

    if (graph == orggraph)
    {
      break;
    }

    graph = graph.finer;

    graph_ReadFromDisk(ctrl, graph);

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.ProjectTmr));
    ASSERT(graph.vwgt != ptr::null_mut());

    ProjectKWayPartition(ctrl, graph);
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.ProjectTmr));
  }

  /* Deal with contiguity requirement at the end */
  ctrl.contig = contig;
  if (contig && FindPartitionInducedComponents(graph, graph.where_, ptr::null_mut(), ptr::null_mut()) > ctrl.nparts) 
    {
      EliminateComponents(ctrl, graph);
    }

  if (!IsBalanced(ctrl, graph, 0.0)) {
    ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
    Greedy_KWayOptimize(ctrl, graph, 10, 0, OMODE_BALANCE); 

    ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
    Greedy_KWayOptimize(ctrl, graph, ctrl.niter, 0, OMODE_REFINE); 
  }

  if (ctrl.contig) 
    {
      ASSERT(FindPartitionInducedComponents(graph, graph.where_, ptr::null_mut(), ptr::null_mut()) == ctrl.nparts);
    }

  IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.UncoarsenTmr));

}


/*************************************************************************/
/* This function allocates memory for the k-way cut-based refinement */
/*************************************************************************/
#[metis_func]
fn AllocateKWayPartitionMemory(ctrl: *mut ctrl_t, graph: *mut graph_t) -> c_void
{

  graph.pwgts  = imalloc(ctrl.nparts*graph.ncon, "AllocateKWayPartitionMemory: pwgts");
  graph.where_  = imalloc(graph.nvtxs,  "AllocateKWayPartitionMemory: where");
  graph.bndptr = imalloc(graph.nvtxs,  "AllocateKWayPartitionMemory: bndptr");
  graph.bndind = imalloc(graph.nvtxs,  "AllocateKWayPartitionMemory: bndind");

  match (ctrl.objtype) {
    ETIS_OBJTYPE_CUT =>
      graph.ckrinfo  = gk_malloc(graph.nvtxs*sizeof(ckrinfo_t), 
                          "AllocateKWayPartitionMemory: ckrinfo"),

    METIS_OBJTYPE_VOL => {
      graph.vkrinfo = gk_malloc(graph.nvtxs*sizeof(vkrinfo_t), 
                          "AllocateKWayVolPartitionMemory: vkrinfo");

      /* This is to let the cut-based -minconn and -contig large-scale graph
         changes to go through */
      graph.ckrinfo = graph.vkrinfo;
},

    _ =>
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype)
  }

}


/*************************************************************************/
/* This function computes the initial id/ed  for cut-based partitioning */
/*************************************************************************/
#[metis_func]
extern "C" fn ComputeKWayPartitionParams(ctrl: *mut ctrl_t, graph: *mut graph_t)
{
  let mut i: idx_t;
  let mut j: idx_t;
  let mut k: idx_t;
  let mut l: idx_t;
  let nvtxs: idx_t;
  let ncon: idx_t;
  let nparts: idx_t;
  let nbnd: idx_t;
  let mincut: idx_t;
  let me: idx_t;
  let other: idx_t;

  let xadj: *mut idx_t;
  let vwgt: *mut idx_t;
  let adjncy: *mut idx_t;
  let adjwgt: *mut idx_t;
  let pwgts: *mut idx_t;
  let where_: *mut idx_t;
  let bndind: *mut idx_t;
  let bndptr: *mut idx_t;
  

  nparts = ctrl.nparts;

  nvtxs  = graph.nvtxs;
  ncon   = graph.ncon;
  xadj   = graph.xadj;
  vwgt   = graph.vwgt;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  where_  = graph.where_;
  pwgts  = iset(nparts*ncon, 0, graph.pwgts);
  bndind = graph.bndind;
  bndptr = iset(nvtxs, -1, graph.bndptr);

  nbnd = mincut = 0;

  /* Compute pwgts */
  if (ncon == 1) {
    for i in 0..nvtxs {
      ASSERT(where_[i] >= 0 && where_[i] < nparts);
      pwgts[where_[i]] += vwgt[i];
    }
  }
  else {
    for i in 0..nvtxs {
      me = where_[i];
      for j in 0..ncon {
        pwgts[me*ncon+j] += vwgt[i*ncon+j];
      }
    }
  }

  /* Compute the required info for refinement */
  match (ctrl.objtype) {
     METIS_OBJTYPE_CUT =>
      {
        ckrinfo_t *myrinfo;
        cnbr_t *mynbrs;

        memset(graph.ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);
        cnbrpoolReset(ctrl);

        for i in 0..nvtxs {
          me      = where_[i];
          myrinfo = graph.ckrinfo+i;

          for j in xadj[i]..xadj[i+1] {
            if (me == where_[adjncy[j]]) {
              myrinfo.id += adjwgt[j];
            } else {
              myrinfo.ed += adjwgt[j];
            }
          }

          /* Time to compute the particular external degrees */
          if (myrinfo.ed > 0) {
            mincut += myrinfo.ed;

            myrinfo.inbr = cnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]);
            mynbrs        = ctrl.cnbrpool + myrinfo.inbr;

            for j in xadj[i]..xadj[i+1] {
              other = where_[adjncy[j]];
              if (me != other) {
                for k in 0..myrinfo.nnbrs {
                  if (mynbrs[k].pid == other) {
                    mynbrs[k].ed += adjwgt[j];
                    break;
                  }
                }
                if (k == myrinfo.nnbrs) {
                  mynbrs[k].pid = other;
                  mynbrs[k].ed  = adjwgt[j];
                  myrinfo.nnbrs+=1 ;
                }
              }
            }

            ASSERT(myrinfo.nnbrs <= xadj[i+1]-xadj[i]);

            /* Only ed-id>=0 nodes are considered to be in the boundary */
            if (myrinfo.ed-myrinfo.id >= 0) {
              BNDInsert(nbnd, bndind, bndptr, i);
            }
          }
          else {
            myrinfo.inbr = -1;
          }
        }

        graph.mincut = mincut/2;
        graph.nbnd   = nbnd;

      ASSERT(CheckBnd2(graph));
      },

    METIS_OBJTYPE_VOL =>
      {
        vkrinfo_t *myrinfo;
        vnbr_t *mynbrs;

        memset(graph.vkrinfo, 0, sizeof(vkrinfo_t)*nvtxs);
        vnbrpoolReset(ctrl);

        /* Compute now the id/ed degrees */
        for i in 0..nvtxs {
          me      = where_[i];
          myrinfo = graph.vkrinfo+i;
      
          for j in xadj[i]..xadj[i+1] {
            if (me == where_[adjncy[j]]) {
              myrinfo.nid+=1;
            } else {
              myrinfo.ned+=1;
            }
          }
      
          /* Time to compute the particular external degrees */
          if (myrinfo.ned > 0) { 
            mincut += myrinfo.ned;

            myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]);
            mynbrs        = ctrl.vnbrpool + myrinfo.inbr;

            for j in xadj[0]..xadj[i+1] {
              other = where_[adjncy[j]];
              if (me != other) {
                for k in 0..myrinfo.nnbrs {
                  if (mynbrs[k].pid == other) {
                    mynbrs[k].ned+=1;
                    break;
                  }
                }
                if (k == myrinfo.nnbrs) {
                  mynbrs[k].gv  = 0;
                  mynbrs[k].pid = other;
                  mynbrs[k].ned = 1;
                  myrinfo.nnbrs+=1;
                }
              }
            }
            ASSERT(myrinfo.nnbrs <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo.inbr = -1;
          }
        }
        graph.mincut = mincut/2;
      
        ComputeKWayVolGains(ctrl, graph);
      ASSERT(graph.minvol == ComputeVolume(graph, graph.where__));
      },
      _ =>
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype)
  }

}


/*************************************************************************/
/* This function projects a partition, and at the same time computes the
 parameters for refinement. */
/*************************************************************************/

#[metis_func]
pub extern "C" fn ProjectKWayPartition(ctrl: *mut ctrl_t, graph: *mut graph_t) -> ()
{
  let i: idx_t;
  let j: idx_t;
  let k: idx_t;
  let nvtxs: idx_t;
  let nbnd: idx_t;
  let nparts: idx_t;
  let me: idx_t;
  let other: idx_t;
  let istart: idx_t;
  let iend: idx_t;
  let tid: idx_t;
  let ted: idx_t;
  

  let xadj: *mut idx_t;
  let adjncy: *mut idx_t;
  let adjwgt: *mut idx_t;
  
  let cmap: *mut idx_t;
  let where_: *mut idx_t;
  let bndptr: *mut idx_t;
  let bndind: *mut idx_t;
  let cwhere: *mut idx_t;
  let htable: *mut idx_t;
  
  let cgraph: *mut graph_t;
  let mut dropedges;


  dropedges = ctrl.dropedges;

  nparts = ctrl.nparts;

  cgraph = graph.coarser;
  cwhere = cgraph.where_;

  if (ctrl.objtype == METIS_OBJTYPE_CUT) {
    ASSERT(CheckBnd2(cgraph));
  } else {
    ASSERT(cgraph.minvol == ComputeVolume(cgraph, cgraph.where_));
  }

  /* free the coarse graph's structure (reduce maxmem) */
  FreeSData(cgraph);

  nvtxs   = graph.nvtxs;
  cmap    = graph.cmap;
  xadj    = graph.xadj;
  adjncy  = graph.adjncy;
  adjwgt  = graph.adjwgt;

  AllocateKWayPartitionMemory(ctrl, graph);

  where_  = graph.where_;
  bndind = graph.bndind;
  bndptr = iset(nvtxs, -1, graph.bndptr);

  htable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

  /* Compute the required info for refinement */
  match (ctrl.objtype) {
    METIS_OBJTYPE_CUT => {
        ckrinfo_t *myrinfo;
        cnbr_t *mynbrs;

        /* go through and project partition and compute id/ed for the nodes */
        for i in 0..nvtxs {
          k        = cmap[i];
          where_[i] = cwhere[k];
          cmap[i]  = if dropedges != 0 { 1 } else {cgraph.ckrinfo[k].ed};  /* For optimization */
        }

        memset(graph.ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);
        cnbrpoolReset(ctrl);

        nbnd = 0;
        for i in 0..nvtxs {
          istart = xadj[i];
          iend   = xadj[i+1];

          myrinfo = graph.ckrinfo+i;

          if (cmap[i] == 0) { /* Interior node. Note that cmap[i] = crinfo[cmap[i]].ed */
            tid = 0;
            for j in istart..iend {
              tid += adjwgt[j];
            }

            myrinfo.id   = tid;
            myrinfo.inbr = -1;
          }
          else { /* Potentially an interface node */
            myrinfo.inbr = cnbrpoolGetNext(ctrl, iend-istart);
            mynbrs        = ctrl.cnbrpool + myrinfo.inbr;

            me = where_[i];
            tid = 0;
            ted = 0;
            for j in istart..iend {
              other = where_[adjncy[j]];
              if (me == other) {
                tid += adjwgt[j];
              }
              else {
                ted += adjwgt[j];
                if ((k = htable[other]) == -1) {
                  htable[other]               = myrinfo.nnbrs;
                  mynbrs[myrinfo.nnbrs].pid  = other;
                  mynbrs[myrinfo.nnbrs+=1].ed = adjwgt[j];
                }
                else {
                  mynbrs[k].ed += adjwgt[j];
                }
              }
            }
            myrinfo.id = tid;
            myrinfo.ed = ted;
      
            /* Remove space for edegrees if it was interior */
            if (ted == 0) { 
              ctrl.nbrpoolcpos -= gk_min(nparts, iend-istart);
              myrinfo.inbr      = -1;
            }
            else {
              if (ted-tid >= 0) {
                BNDInsert(nbnd, bndind, bndptr, i);
              }
      
              for j in 0..myrinfo.nnbrs {
                htable[mynbrs[j].pid] = -1;
              }
            }
          }
        }
      
        graph.nbnd = nbnd;
      ASSERT(CheckBnd2(graph));
      },

    METIS_OBJTYPE_VOL =>
      {
        
        let myrinfo: *mut vkrinfo_t;
        let mynbrs: *mut vnbr_t;

        /* go through and project partition and compute id/ed for the nodes */
        for i in 0..nvtxs {
          k        = cmap[i];
          where_[i] = cwhere[k];
          cmap[i]  = if dropedges != 0 { 1 } else {cgraph.vkrinfo[k].ned};  /* For optimization */
        }

        memset(graph.vkrinfo, 0, sizeof(vkrinfo_t)*nvtxs);
        vnbrpoolReset(ctrl);

        for i in 0..nvtxs {
          istart = xadj[i];
          iend   = xadj[i+1];
          myrinfo = graph.vkrinfo+i;

          if (cmap[i] == 0) { /* Note that cmap[i] = crinfo[cmap[i]].ed */
            myrinfo.nid  = iend-istart;
            myrinfo.inbr = -1;
          }
          else { /* Potentially an interface node */
            myrinfo.inbr = vnbrpoolGetNext(ctrl, iend-istart);
            mynbrs        = ctrl.vnbrpool + myrinfo.inbr;

            me = where_[i];
            tid = 0;
            ted = 0;
            for j in istart..iend {
              other = where_[adjncy[j]];
              if (me == other) {
                tid+=1;
              }
              else {
                ted+=1;
                if ((k = htable[other]) == -1) {
                  htable[other]                = myrinfo.nnbrs;
                  mynbrs[myrinfo.nnbrs].gv    = 0;
                  mynbrs[myrinfo.nnbrs].pid   = other;
                  mynbrs[myrinfo.nnbrs+=1].ned = 1;
                }
                else {
                  mynbrs[k].ned+=1;
                }
              }
            }
            myrinfo.nid = tid;
            myrinfo.ned = ted;
      
            /* Remove space for edegrees if it was interior */
            if (ted == 0) { 
              ctrl.nbrpoolcpos -= gk_min(nparts, iend-istart);
              myrinfo.inbr = -1;
            }
            else {
              for j in 0..myrinfo.nnbrs {
                htable[mynbrs[j].pid] = -1;
              }
            }
          }
        }
      
        ComputeKWayVolGains(ctrl, graph);

        ASSERT(graph.minvol == ComputeVolume(graph, graph.where_));
      },

  _ =>
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype)
  }

  graph.mincut = if dropedges { ComputeCut(graph, where_) } else { cgraph.mincut };
  icopy(nparts*graph.ncon, cgraph.pwgts, graph.pwgts);

  FreeGraph(&graph.coarser);

}


/*************************************************************************/
/* This function computes the boundary definition for balancing. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn ComputeKWayBoundary(ctrl: *mut ctrl_t, graph: *mut graph_t, bndtype: idx_t) {
  let nvtxs;
  let nbnd;
  let bndind: *const idx_t;
  let bndptr: *const idx_t;

  nvtxs  = graph.nvtxs;
  bndind = graph.bndind;
  bndptr = iset(nvtxs, -1, graph.bndptr);

  nbnd = 0;

  match (ctrl.objtype) {
    METIS_OBJTYPE_CUT => {
      /* Compute the boundary */
      if (bndtype == BNDTYPE_REFINE) {
        for i in 0..nvtxs {
          if (graph.ckrinfo[i].ed > 0 && graph.ckrinfo[i].ed-graph.ckrinfo[i].id >= 0) {
        BNDInsert(nbnd, bndind, bndptr, i);
      }
        }
      }
      else { /* BNDTYPE_BALANCE */
        for i in 0..nvtxs {
          if (graph.ckrinfo[i].ed > 0) 
      {
        BNDInsert(nbnd, bndind, bndptr, i);
      }
        }
      }
    },

    METIS_OBJTYPE_VOL => {
      /* Compute the boundary */
      if (bndtype == BNDTYPE_REFINE) {
        for i in 0..nvtxs {
          if (graph.vkrinfo[i].gv >= 0) 
          {
            BNDInsert(nbnd, bndind, bndptr, i);
          }
        }
      }
      else { /* BNDTYPE_BALANCE */
        for i in 0..nvtxs {
          if (graph.vkrinfo[i].ned > 0) 
          {
            BNDInsert(nbnd, bndind, bndptr, i);
          }
        }
      }
    },
    _ =>
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype)
  }

  graph.nbnd = nbnd;
}


/*************************************************************************/
/* This function computes the initial gains in the communication volume */
/*************************************************************************/
#[metis_func]
pub extern "C" fn ComputeKWayVolGains(ctrl: *mut ctrl_t, graph: *mut graph_t) {
  let i: idx_t;
  let ii: idx_t;
  let j: idx_t;
  let k: idx_t;
  let l: idx_t;
  let nvtxs: idx_t;
  let nparts: idx_t;
  let me: idx_t;
  let other: idx_t;
  let pid: idx_t;
  
  let dx_t: *mut idx_t;
  let xadj: *mut idx_t;
  let vsize: *mut idx_t;
  let adjncy: *mut idx_t;
  let adjwgt: *mut idx_t;
  let where_: *mut idx_t;
  let bndind: *mut idx_t;
  let bndptr: *mut idx_t;
  let ophtable: *mut idx_t;
  
  let myrinfo: *mut vkrinfo_t;
  let orinfo: *mut vkrinfo_t;

  let mynbrs: *mut vnbr_t;
  let onbrs: *mut vnbr_t;


  nparts = ctrl.nparts;

  nvtxs  = graph.nvtxs;
  xadj   = graph.xadj;
  vsize  = graph.vsize;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  where_  = graph.where_;
  bndind = graph.bndind;
  bndptr = iset(nvtxs, -1, graph.bndptr);

  ophtable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

  /* Compute the volume gains */
  graph.minvol = graph.nbnd = 0;
  for i in 0..nvtxs {
    myrinfo     = graph.vkrinfo+i;
    myrinfo.gv = IDX_MIN;

    if (myrinfo.nnbrs > 0) {
      me     = where_[i];
      mynbrs = ctrl.vnbrpool + myrinfo.inbr;

      graph.minvol += myrinfo.nnbrs*vsize[i];

      for j in xadj[i]..xadj[i+1] {
        ii     = adjncy[j];
        other  = where_[ii];
        orinfo = graph.vkrinfo+ii;
        onbrs  = ctrl.vnbrpool + orinfo.inbr;

        for k in 0..orinfo.nnbrs {
          ophtable[onbrs[k].pid] = k;
        }
        ophtable[other] = 1;  /* this is to simplify coding */

        if (me == other) {
          /* Find which domains 'i' is connected to but 'ii' is not 
             and update their gain */
          for k in 0..myrinfo.nnbrs {
            if (ophtable[mynbrs[k].pid] == -1) {
              mynbrs[k].gv -= vsize[ii];
            }
          }
        }
        else {
          ASSERT(ophtable[me] != -1);

          if (onbrs[ophtable[me]].ned == 1) { 
            /* I'm the only connection of 'ii' in 'me' */
            /* Increase the gains for all the common domains between 'i' and 'ii' */
            for k in 0..myrinfo.nnbrs {
              if (ophtable[mynbrs[k].pid] != -1) {
                mynbrs[k].gv += vsize[ii];
              }
            }
          }
          else {
            /* Find which domains 'i' is connected to and 'ii' is not 
               and update their gain */
            for k in 0..myrinfo.nnbrs {
              if (ophtable[mynbrs[k].pid] == -1) 
              {
                mynbrs[k].gv -= vsize[ii];
              }
            }
          }
        }

        /* Reset the marker vector */
        for k in 0..orinfo.nnbrs {
          ophtable[onbrs[k].pid] = -1;
        }
        ophtable[other] = -1;
      }

      /* Compute the max vgain */
      for k in 0..myrinfo.nnbrs {
        if (mynbrs[k].gv > myrinfo.gv)
        {
          myrinfo.gv = mynbrs[k].gv;
        }
      }

      /* Add the extra gain due to id == 0 */
      if (myrinfo.ned > 0 && myrinfo.nid == 0)
      {
        myrinfo.gv += vsize[i];
      }
    }

    if (myrinfo.gv >= 0)
    {
      BNDInsert(graph.nbnd, bndind, bndptr, i);
    }
  }

}


/*************************************************************************/
/* This function checks if the partition weights are within the balance
constraints */
/*************************************************************************/
#[metis_func]
extern "C" fn IsBalanced(ctrl: *mut ctrl_t, graph: *mut graph_t, ffactor: real_t) -> std::ffi::c_int
{
  return 
    (ComputeLoadImbalanceDiff(graph, ctrl.nparts, ctrl.pijbm, ctrl.ubfactors) 
         <= ffactor);
}

