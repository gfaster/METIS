/*!
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 */

use crate::*;

/*************************************************************************/
/*! This function computes the initial bisection of the coarsest graph */
/*************************************************************************/
#[metis_func]
pub fn Init2WayPartition(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    let dbglvl;

    ASSERT(graph.tvwgt[0] >= 0);

    dbglvl = ctrl.dbglvl;
    IFSET(ctrl.dbglvl, METIS_DBG_REFINE, ctrl.dbglvl -= METIS_DBG_REFINE);
    IFSET(ctrl.dbglvl, METIS_DBG_MOVEINFO, ctrl.dbglvl -= METIS_DBG_MOVEINFO);

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.InitPartTmr));

    match (ctrl.iptype) {
        METIS_IPTYPE_RANDOM => {
        if (graph.ncon == 1)
    {
        RandomBisection(ctrl, graph, ntpwgts, niparts);
    }
        else
    {
        McRandomBisection(ctrl, graph, ntpwgts, niparts);
    }
},

        METIS_IPTYPE_GROW => {
        if (graph.nedges == 0)
            {
                if (graph.ncon == 1)
                {
                    RandomBisection(ctrl, graph, ntpwgts, niparts);
                } else {
                    McRandomBisection(ctrl, graph, ntpwgts, niparts);
                }
            }
        else
                    {
                    if (graph.ncon == 1)
                {
                    GrowBisection(ctrl, graph, ntpwgts, niparts);
                } else {
                McGrowBisection(ctrl, graph, ntpwgts, niparts);
                }
            }
        },

        _ =>
        gk_errexit(SIGERR, "Unknown initial partition type: %d\n", ctrl.iptype)
    }

    IFSET(ctrl.dbglvl, METIS_DBG_IPART, printf("Initial Cut: {}\n", graph.mincut));
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.InitPartTmr));
    ctrl.dbglvl = dbglvl;

}


/*************************************************************************/
/*! This function computes the initial separator of the coarsest graph */
/*************************************************************************/
#[metis_func]
pub fn InitSeparator(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) 
{
    let mut ntpwgts = [0.5, 0.5];
  let dbglvl: mdbglvl_et;

  dbglvl = ctrl.dbglvl;
  IFSET(ctrl.dbglvl, METIS_DBG_REFINE, ctrl.dbglvl -= METIS_DBG_REFINE);
  IFSET(ctrl.dbglvl, METIS_DBG_MOVEINFO, ctrl.dbglvl -= METIS_DBG_MOVEINFO);

  IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.InitPartTmr));

  /* this is required for the cut-based part of the refinement */
  Setup2WayBalMultipliers(ctrl, graph, ntpwgts);

  match (ctrl.iptype) {
    METIS_IPTYPE_EDGE => {
      if (graph.nedges == 0)
            {
                RandomBisection(ctrl, graph, ntpwgts, niparts);
            }
      else
            {
                GrowBisection(ctrl, graph, ntpwgts, niparts);
            }
      Compute2WayPartitionParams(ctrl, graph);
      ConstructSeparator(ctrl, graph);
},

    METIS_IPTYPE_NODE => {
      GrowBisectionNode(ctrl, graph, ntpwgts, niparts);
    },

    _ =>
      gk_errexit(SIGERR, "Unknown iptype of {}\n", ctrl.iptype)
  }

  IFSET(ctrl.dbglvl, METIS_DBG_IPART, printf("Initial Sep: {}\n", graph.mincut));
  IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.InitPartTmr));

  ctrl.dbglvl = dbglvl;

}


/*************************************************************************/
/*! This function computes a bisection of a graph by randomly assigning
    the vertices followed by a bisection refinement.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn RandomBisection(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
  // idx_t i, ii, j, k, nvtxs, pwgts[2], zeromaxpwgt, from, me, 
  //       bestcut=0, icut, mincut, inbfs;
  // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
  // idx_t *perm, *bestwhere;

let bestcut = 0;
let pwgts: [_;2];

  WCOREPUSH;

  nvtxs  = graph.nvtxs;
  xadj   = graph.xadj;
  vwgt   = graph.vwgt;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where_ = graph.where_;

  bestwhere = iwspacemalloc(ctrl, nvtxs);
  perm      = iwspacemalloc(ctrl, nvtxs);

  zeromaxpwgt = ctrl.ubfactors[0]*graph.tvwgt[0]*ntpwgts[0];

  for inbfs in 0..niparts {
    iset(nvtxs, 1, where_);

    if (inbfs > 0) {
      irandArrayPermute(nvtxs, perm, nvtxs/2, 1);
      pwgts[1] = graph.tvwgt[0];
      pwgts[0] = 0;

      for ii in 0..nvtxs {
        i = perm[ii];
        if (pwgts[0]+vwgt[i] < zeromaxpwgt) {
          where_[i] = 0;
          pwgts[0] += vwgt[i];
          pwgts[1] -= vwgt[i];
          if (pwgts[0] > zeromaxpwgt)
                    {
                        break;
                    }
        }
      }
    }

    /* Do some partition refinement  */
    Compute2WayPartitionParams(ctrl, graph);
    /* printf("IPART: %3"PRIDX" [%5"PRIDX" %5"PRIDX"] [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", graph.nvtxs, pwgts[0], pwgts[1], graph.pwgts[0], graph.pwgts[1], graph.mincut); */

    Balance2Way(ctrl, graph, ntpwgts);
    /* printf("BPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", graph.pwgts[0], graph.pwgts[1], graph.mincut); */

    FM_2WayRefine(ctrl, graph, ntpwgts, 4);
    /* printf("RPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", graph.pwgts[0], graph.pwgts[1], graph.mincut); */

    if (inbfs==0 || bestcut > graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
      if (bestcut == 0)
    {
        break;
    }
}
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}


/*************************************************************************/
/*! This function takes a graph and produces a bisection by using a region
    growing algorithm. The resulting bisection is refined using FM.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn GrowBisection(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{

  // idx_t i, j, k, nvtxs, drain, nleft, first, last, 
  //       pwgts[2], oneminpwgt, onemaxpwgt, 
  //       from, me, bestcut=0, icut, mincut, inbfs;
  // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
  // idx_t *queue, *touched, *gain, *bestwhere;

  WCOREPUSH;

  nvtxs  = graph.nvtxs;
  xadj   = graph.xadj;
  vwgt   = graph.vwgt;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where_ = graph.where_;

  bestwhere = iwspacemalloc(ctrl, nvtxs);
  queue     = iwspacemalloc(ctrl, nvtxs);
  touched   = iwspacemalloc(ctrl, nvtxs);

  onemaxpwgt = ctrl.ubfactors[0]*graph.tvwgt[0]*ntpwgts[1];
  oneminpwgt = (1.0/ctrl.ubfactors[0])*graph.tvwgt[0]*ntpwgts[1];

  for inbfs in 0..niparts {
    iset(nvtxs, 1, where_);

    iset(nvtxs, 0, touched);

    pwgts[1] = graph.tvwgt[0];
    pwgts[0] = 0;


    queue[0] = irandInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; 
    last  = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    loop {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0 || drain)
{
    break;
}
        k = irandInRange(nleft);
        for i in 0..nvtxs {
          if (touched[i] == 0) {
            if (k == 0)
            {
                break;
            }
            else
            {
                k-= 1;
            }
        }
        }

        queue[0]   = i;
        touched[i] = 1;
        first      = 0; 
        last       = 1;
        nleft-=1;
      }

      i = queue[first];
first += 1;
      if (pwgts[0] > 0 && pwgts[1]-vwgt[i] < oneminpwgt) {
        drain = 1;
        continue;
      }

      where_[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= onemaxpwgt)
{
    break;
}
      drain = 0;
      for j in xadj[i]..xadj[i+1] {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last] = k;
            last += 1;
          touched[k] = 1;
          nleft-= 1;
            
        }
      }
    }

    /* Check to see if we hit any bad limiting cases */
    if (pwgts[1] == 0) 
    {
        where_[irandInRange(nvtxs)] = 1;
    }
    if (pwgts[0] == 0) 
    {
where_[irandInRange(nvtxs)] = 0;
      }
    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    /*
    printf("IPART: %3"PRIDX" [%5"PRIDX" %5"PRIDX"] [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", 
        graph.nvtxs, pwgts[0], pwgts[1], graph.pwgts[0], graph.pwgts[1], graph.mincut); 
    */

    Balance2Way(ctrl, graph, ntpwgts);
    /*
    printf("BPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", graph.pwgts[0],
        graph.pwgts[1], graph.mincut); 
    */

    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
    /*
    printf("RPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"\n", graph.pwgts[0], 
        graph.pwgts[1], graph.mincut);
    */

    if (inbfs == 0 || bestcut > graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
      if (bestcut == 0)
            {
                break;
            }
        }
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}


/*************************************************************************/
/*! This function takes a multi-constraint graph and computes a bisection 
    by randomly assigning the vertices and then refining it. The resulting
    partition is returned in graph.where_.
*/
/**************************************************************************/
#[metis_func]
pub fn McRandomBisection(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
  // idx_t i, ii, j, k, nvtxs, ncon, from, bestcut=0, mincut, inbfs, qnum;
  // idx_t *bestwhere, *where_, *perm, *counts;
  // idx_t *vwgt;

  WCOREPUSH;

  nvtxs = graph.nvtxs;
  ncon  = graph.ncon;
  vwgt  = graph.vwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where_ = graph.where_;

  bestwhere = iwspacemalloc(ctrl, nvtxs);
  perm      = iwspacemalloc(ctrl, nvtxs);
  counts    = iwspacemalloc(ctrl, ncon);

  for inbfs in 0..2*niparts {
    irandArrayPermute(nvtxs, perm, nvtxs/2, 1);
    iset(ncon, 0, counts);

    /* partition by splitting the queues randomly */
    for ii in 0..nvtxs {
      i        = perm[ii];
      qnum     = iargmax(ncon, vwgt+i*ncon,1);
      where_[i] = (counts[qnum])%2;
counts[qnum]+=1;
    }

    Compute2WayPartitionParams(ctrl, graph);

    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
    Balance2Way(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
    Balance2Way(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

    if (inbfs == 0 || bestcut >= graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
      if (bestcut == 0)
            {
                break;
            }
        }
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}


/*************************************************************************/
/*! This function takes a multi-constraint graph and produces a bisection 
    by using a region growing algorithm. The resulting partition is 
    returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn McGrowBisection(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
  // idx_t i, j, k, nvtxs, ncon, from, bestcut=0, mincut, inbfs;
  // idx_t *bestwhere, *where_;

  WCOREPUSH;

  nvtxs = graph.nvtxs;

  Allocate2WayPartitionMemory(ctrl, graph);
  where_ = graph.where_;

  bestwhere = iwspacemalloc(ctrl, nvtxs);

  for inbfs in 0..2*niparts {
    iset(nvtxs, 1, where_);
    where_[irandInRange(nvtxs)] = 0;

    Compute2WayPartitionParams(ctrl, graph);

    Balance2Way(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
    Balance2Way(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

    if (inbfs == 0 || bestcut >= graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
      if (bestcut == 0)
            {
                break;
            }
        }
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}


/*************************************************************************/
/* This function takes a graph and produces a tri-section into left, right,
   and separator using a region growing algorithm. The resulting separator
   is refined using node FM.
   The resulting partition is returned in graph.where_.
*/
/**************************************************************************/
#[metis_func]
pub fn GrowBisectionNode(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
  // idx_t i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], oneminpwgt, 
  //       onemaxpwgt, from, me, bestcut=0, icut, mincut, inbfs;
  // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_, *bndind;
  // idx_t *queue, *touched, *gain, *bestwhere;

  WCOREPUSH;

  nvtxs  = graph.nvtxs;
  xadj   = graph.xadj;
  vwgt   = graph.vwgt;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  bestwhere = iwspacemalloc(ctrl, nvtxs);
  queue     = iwspacemalloc(ctrl, nvtxs);
  touched   = iwspacemalloc(ctrl, nvtxs);

  onemaxpwgt = ctrl.ubfactors[0]*graph.tvwgt[0]*0.5;
  oneminpwgt = (1.0/ctrl.ubfactors[0])*graph.tvwgt[0]*0.5;


  /* Allocate refinement memory. Allocate sufficient memory for both edge and node */
  graph.pwgts  = imalloc(3, "GrowBisectionNode: pwgts");
  graph.where_  = imalloc(nvtxs, "GrowBisectionNode: where_");
  graph.bndptr = imalloc(nvtxs, "GrowBisectionNode: bndptr");
  graph.bndind = imalloc(nvtxs, "GrowBisectionNode: bndind");
  graph.id     = imalloc(nvtxs, "GrowBisectionNode: id");
  graph.ed     = imalloc(nvtxs, "GrowBisectionNode: ed");
  graph.nrinfo = gk_malloc(nvtxs*sizeof(nrinfo_t), "GrowBisectionNode: nrinfo");
  
  where_  = graph.where_;
  bndind = graph.bndind;

  for inbfs in 0..niparts {
    iset(nvtxs, 1, where_);
    iset(nvtxs, 0, touched);

    pwgts[1] = graph.tvwgt[0];
    pwgts[0] = 0;

    queue[0] = irandInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    loop {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0 || drain)
{
    break;
}
        k = irandInRange(nleft);
        for i in 0..nvtxs { /* select the kth untouched vertex */
          if (touched[i] == 0) {
            if (k == 0)
                        {
                            break;
                        }
            else
                        {
                            k-=1;
                        }
                    }
        }

        queue[0]   = i;
        touched[i] = 1;
        first      = 0; 
        last       = 1;
        nleft-=1;
      }

      i = queue[first];
            first += 1;
      if (pwgts[1]-vwgt[i] < oneminpwgt) {
        drain = 1;
        continue;
      }

      where_[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= onemaxpwgt)
            {
                break;
            }
      drain = 0;
      for j in xadj[i]..xadj[i+1] {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last] = k;
        last += 1;
          touched[k] = 1;
        nleft -= 1;
        }
      }
    }

    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    Balance2Way(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, 4);

    /* Construct and refine the vertex separator */
    for i in 0..graph.nbnd {
      j = bndind[i];
      if (xadj[j+1]-xadj[j] > 0) /* ignore islands */
            {
                where_[j] = 2;
            }
        }

    Compute2WayNodePartitionParams(ctrl, graph); 
    FM_2WayNodeRefine2Sided(ctrl, graph, 1);
    FM_2WayNodeRefine1Sided(ctrl, graph, 4);

    /*
    printf("ISep: [{} %"PRIDX" %"PRIDX" %"PRIDX"] %"PRIDX"\n", 
        inbfs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2], bestcut); 
    */
    
    if (inbfs == 0 || bestcut > graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
    }
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}


/*************************************************************************/
/* This function takes a graph and produces a tri-section into left, right,
   and separator using a region growing algorithm. The resulting separator
   is refined using node FM.
   The resulting partition is returned in graph.where_.
*/
/**************************************************************************/
    #[metis_func]
pub fn GrowBisectionNode2(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
{
  // idx_t i, j, k, nvtxs, bestcut=0, mincut, inbfs;
  // idx_t *xadj, *where_, *bndind, *bestwhere;

  WCOREPUSH;

  nvtxs  = graph.nvtxs;
  xadj   = graph.xadj;

  /* Allocate refinement memory. Allocate sufficient memory for both edge and node */
  graph.pwgts  = imalloc(3, "GrowBisectionNode: pwgts");
  graph.where_  = imalloc(nvtxs, "GrowBisectionNode: where_");
  graph.bndptr = imalloc(nvtxs, "GrowBisectionNode: bndptr");
  graph.bndind = imalloc(nvtxs, "GrowBisectionNode: bndind");
  graph.id     = imalloc(nvtxs, "GrowBisectionNode: id");
  graph.ed     = imalloc(nvtxs, "GrowBisectionNode: ed");
  graph.nrinfo = gk_malloc(nvtxs*sizeof(nrinfo_t), "GrowBisectionNode: nrinfo");
  
  bestwhere = iwspacemalloc(ctrl, nvtxs);

  where_  = graph.where_;
  bndind = graph.bndind;

  for inbfs in 0..niparts {
    iset(nvtxs, 1, where_);
    if (inbfs > 0)
        {
            where_[irandInRange(nvtxs)] = 0;
        }
    Compute2WayPartitionParams(ctrl, graph);
    General2WayBalance(ctrl, graph, ntpwgts);
    FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

    /* Construct and refine the vertex separator */
    for i in 0..graph.nbnd {
      j = bndind[i];
      if (xadj[j+1]-xadj[j] > 0) /* ignore islands */
            {
                where_[j] = 2;
            }
        }

    Compute2WayNodePartitionParams(ctrl, graph); 
    FM_2WayNodeRefine2Sided(ctrl, graph, 4);

    /*
    printf("ISep: [{} %"PRIDX" %"PRIDX" %"PRIDX"] %"PRIDX"\n", 
        inbfs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2], bestcut); 
    */

    if (inbfs == 0 || bestcut > graph.mincut) {
      bestcut = graph.mincut;
      icopy(nvtxs, where_, bestwhere);
    }
  }

  graph.mincut = bestcut;
  icopy(nvtxs, bestwhere, where_);

  WCOREPOP;
}

