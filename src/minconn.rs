/*
\file 
\brief Functions that deal with prunning the number of adjacent subdomains in kmetis

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: minconn.c 17513 2014-08-05 16:20:50Z dominique $
*/

use crate::*;


/*************************************************************************/
/* This function computes the subdomain graph storing the result in the
pre-allocated worspace arrays */
/*************************************************************************/
#[metis_func]
pub extern "C" fn ComputeSubDomainGraph(ctrl: *mut ctrl_t, graph: *mut graph_t) 
{
  // idx_t i, ii, j, pid, other, nparts, nvtxs, nnbrs;
  // idx_t *xadj, *adjncy, *adjwgt, *where_;
  // idx_t *pptr, *pind;
  // idx_t nads=0, *vadids, *vadwgts;

  // WCOREPUSH;

  let   nvtxs  = graph.nvtxs as usize;
  xadj   = graph.xadj;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;
  where_  = graph.where_;

  nparts = ctrl.nparts; 

  vadids  = ctrl.pvec1;
  vadwgts = iset(nparts, 0, ctrl.pvec2);

  pptr = vec![0; nparts+1 as usize];
  pind = vec![0; nvtxs as usize];
  iarray2csr(nvtxs, nparts, where_, pptr, pind);

  for pid in (0)..(nparts) {
    match (ctrl.objtype) {
      METIS_OBJTYPE_CUT =>
      {
        ckrinfo_t *rinfo;
        cnbr_t *nbrs;

        rinfo = graph.ckrinfo;
        nads=0;
        for ii in (pptr[pid as usize])..(pptr[(pid+1) as usize]) {
          i = pind[ii as usize];
          assert!(pid == where_[i as usize]);

          if (rinfo[i as usize].ed > 0) {
            nnbrs = rinfo[i as usize].nnbrs;
            nbrs  = ctrl.cnbrpool + rinfo[i as usize].inbr;

            for j in (0)..(nnbrs) {
              other = nbrs[j as usize].pid;
              if (vadwgts[other as usize] == 0) {
                vadids[(nads) as usize] = other;nads+=1;
              }
              vadwgts[other as usize] += nbrs[j as usize].ed;
            }
          }
        }
      }

      METIS_OBJTYPE_VOL =>
      {
        vkrinfo_t *rinfo;
        vnbr_t *nbrs;

        rinfo = graph.vkrinfo;
        nads=0;
        for ii in (pptr[pid as usize])..(pptr[(pid+1) as usize]) {
          i = pind[ii as usize];
          assert!(pid == where_[i as usize]);

          if (rinfo[i as usize].ned > 0) {
            nnbrs = rinfo[i as usize].nnbrs;
            nbrs  = ctrl.vnbrpool + rinfo[i as usize].inbr;

            for j in (0)..(nnbrs) {
              other = nbrs[j as usize].pid;
              if (vadwgts[other as usize] == 0) {
                vadids[(nads) as usize] = other;nads+=1;
              }
              vadwgts[other as usize] += nbrs[j as usize].ned;
            }
          }
        }
      }

      _ => gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl.objtype)
    }

    /* See if you have enough memory to store the adjacent info for that subdomain */
    if (ctrl.maxnads[pid as usize] < nads) {
      ctrl.maxnads[pid as usize] = 2*nads;
      ctrl.adids[pid as usize]   = irealloc(ctrl.adids[pid as usize], ctrl.maxnads[pid as usize], 
        "ComputeSubDomainGraph: adids[pid as usize]");
      ctrl.adwgts[pid as usize]  = irealloc(ctrl.adwgts[pid as usize], ctrl.maxnads[pid as usize], 
        "ComputeSubDomainGraph: adids[pid as usize]");
    }

    ctrl.nads[pid as usize] = nads;
    for j in (0)..(nads) {
      ctrl.adids[pid as usize][j as usize]  = vadids[j as usize];
      ctrl.adwgts[pid as usize][j as usize] = vadwgts[vadids[j as usize] as usize];

      vadwgts[vadids[j as usize] as usize] = 0;
    }
  }

  // WCOREPOP;
}


/*************************************************************************/
/* This function updates the weight of an edge in the subdomain graph by
adding to it the value of ewgt. The update can either increase or
decrease the weight of the subdomain edge based on the value of ewgt.

\param u is the ID of one of the incident subdomains to the edge
\param v is the ID of the other incident subdomains to the edge
\param ewgt is the weight to be added to the subdomain edge
\param nparts is the number of subdomains
\param r_maxndoms is the maximum number of adjacent subdomains and is
updated as necessary. The update is skipped if a std::ptr::null_mut() value is
supplied.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn UpdateEdgeSubDomainGraph(ctrl: *mut ctrl_t, u: idx_t, v: idx_t, ewgt: idx_t, r_maxndoms: *mut idx_t) 
{
  // idx_t i, j, nads;

  if (ewgt == 0) {
    return;
  }

  for i in (0)..(2) {
    nads = ctrl.nads[u as usize];
    /* Find the edge */
    for j in (0)..(nads) {
      if (ctrl.adids[u as usize][j as usize] == v) {
        ctrl.adwgts[u as usize][j as usize] += ewgt;
        break;
      }
    }

    if (j == nads) {
      /* Deal with the case in which the edge was not found */
      assert!(ewgt > 0);
      if (ctrl.maxnads[u as usize] == nads) {
        ctrl.maxnads[u as usize] = 2*(nads+1);
        ctrl.adids[u as usize]   = irealloc(ctrl.adids[u as usize], ctrl.maxnads[u as usize], 
          "IncreaseEdgeSubDomainGraph: adids[pid as usize]");
        ctrl.adwgts[u as usize]  = irealloc(ctrl.adwgts[u as usize], ctrl.maxnads[u as usize], 
          "IncreaseEdgeSubDomainGraph: adids[pid as usize]");
      }
      ctrl.adids[u as usize][nads as usize]  = v;
      ctrl.adwgts[u as usize][nads as usize] = ewgt;
      nads += 1;
      if (r_maxndoms != std::ptr::null_mut() && nads > *r_maxndoms) {
        print!("You just increased the maxndoms: {:} {:}\n", 
          nads, *r_maxndoms);
        *r_maxndoms = nads;
      }
    }
    else {
      /* See if the updated edge becomes 0 */
      assert!(ctrl.adwgts[u as usize][j as usize] >= 0);
      if (ctrl.adwgts[u as usize][j as usize] == 0) {
        ctrl.adids[u as usize][j as usize]  = ctrl.adids[u as usize][(nads-1) as usize];
        ctrl.adwgts[u as usize][j as usize] = ctrl.adwgts[u as usize][(nads-1) as usize];
        nads;nads-=1;
        if (r_maxndoms != std::ptr::null_mut() && nads+1 == *r_maxndoms) {
          *r_maxndoms = ctrl.nads[(iargmax(ctrl.nparts, ctrl.nads,1)) as usize];
        }
      }
    }
    ctrl.nads[u as usize] = nads;

    SWAP(u, v, j);
  }
}


/*************************************************************************/
/* This function computes the subdomain graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn EliminateSubDomainEdges(ctrl: *mut ctrl_t, graph: *mut graph_t) {
  // idx_t i, ii, j, k, ncon, nparts, scheme, pid_from, pid_to, me, other, nvtxs, 
  // total, max, avg, totalout, nind=0, ncand=0, ncand2, target, target2, 
  // nadd, bestnadd=0;
  // idx_t min, move, *cpwgt;
  // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *pwgts, *where_, *maxpwgt, 
  // *mypmat, *otherpmat, *kpmat, *ind;
  // idx_t *nads, **adids, **adwgts;
  // ikv_t *cand, *cand2;
  // ipq_t queue;
  // real_t *tpwgts, badfactor=1.4;
  // idx_t *pptr, *pind;
  // idx_t *vmarker=std::ptr::null_mut(), *pmarker=std::ptr::null_mut(), *modind=std::ptr::null_mut();  /* volume specific work arrays */

  // WCOREPUSH;

  let   nvtxs  = graph.nvtxs as usize;
  let   ncon   = graph.ncon as usize;
  xadj   = graph.xadj;
  adjncy = graph.adjncy;
  vwgt   = graph.vwgt;
  adjwgt = ( if (ctrl.objtype == METIS_OBJTYPE_VOL)  { (std::ptr::null_mut()) } else { (graph.adjwgt) });

  where_ = graph.where_;
  pwgts = graph.pwgts;  /* We assume that this is properly initialized */

  nparts = ctrl.nparts;
  tpwgts = ctrl.tpwgts;

  cpwgt     = vec![0; ncon as usize];
  maxpwgt   = vec![0; nparts*ncon as usize];
  ind       = vec![0; nvtxs as usize];
  otherpmat = vec![0; nparts as usize];

  cand  = ikvwspacemalloc(ctrl, nparts);
  cand2 = ikvwspacemalloc(ctrl, nparts);

  pptr = vec![0; nparts+1 as usize];
  pind = vec![0; nvtxs as usize];
  iarray2csr(nvtxs, nparts, where_, pptr, pind);

  if (ctrl.objtype == METIS_OBJTYPE_VOL) {
    /* Vol-refinement specific working arrays */
    modind  = vec![0; nvtxs as usize];
    vmarker = vec![0; nvtxs as usize];
    pmarker = vec![-1; nparts as usize];
  }


  /* Compute the pmat matrix and ndoms */
  ComputeSubDomainGraph(ctrl, graph);

  nads   = ctrl.nads;
  adids  = ctrl.adids;
  adwgts = ctrl.adwgts;

  mypmat = iset(nparts, 0, ctrl.pvec1);
  kpmat  = iset(nparts, 0, ctrl.pvec2);

  /* Compute the maximum allowed weight for each domain */
  for i in (0)..(nparts) {
    for j in (0)..(ncon) {
      maxpwgt[(i*ncon+j) as usize] =
        (if ncon == 1 { 1.25 } else { 1.025 })
        * tpwgts[i as usize]
        * (*graph.tvwgt[j as usize])
        * (*ctrl.ubfactors[j as usize]);
    }
  }

  ipqInit(&queue, nparts);

  /* Get into the loop eliminating subdomain connections */
  loop {
    total = isum(nparts, nads, 1);
    avg   = total/nparts;
    max   = nads[(iargmax(nparts, nads,1)) as usize];

    IFSET(ctrl.dbglvl, METIS_DBG_CONNINFO, 
      print!("Adjacent Subdomain Stats: Total: {:3}, "
        "Max: {:3}[%zu], Avg: {:3}\n", 
        total, max, iargmax(nparts, nads,1), avg)); 

    if (max < badfactor*avg) {
      break;
    }

    /* Add the subdomains that you will try to reduce their connectivity */
    ipqReset(&queue);
    for i in (0)..(nparts) {
      if (nads[i as usize] >= avg + (max-avg)/2) {
        ipqInsert(&queue, i, nads[i as usize]);
      }
    }

    move_ = 0;
    while let Some(me) = queue.pop() {
      totalout = isum(nads[me as usize], adwgts[me as usize], 1);
      ncand2=0;
      for i in 0..nads[me as usize] {
        mypmat[adids[me as usize][i as usize] as usize] = adwgts[me as usize][i as usize];

        /* keep track of the weakly connected adjacent subdomains */
        if (2*nads[me as usize]*adwgts[me as usize][i as usize] < totalout) {
          cand2[ncand2 as usize].val   = adids[me as usize][i as usize];
          cand2[(ncand2) as usize].key = adwgts[me as usize][i as usize];ncand2+=1;
        }
      }

      IFSET(ctrl.dbglvl, METIS_DBG_CONNINFO, 
        print!("Me: {:}, Degree: {:4}, TotalOut: {:},\n", 
          me, nads[me as usize], totalout));

      /* Sort the connections according to their cut */
      ikvsorti(ncand2, cand2);

      /* Two schemes are used for eliminating subdomain edges.
      The first, tries to eliminate subdomain edges by moving remote groups 
      of vertices to subdomains that 'me' is already connected to.
      The second, tries to eliminate subdomain edges by moving entire sets of 
      my vertices that connect to the 'other' subdomain to a subdomain that 
      I'm already connected to.
      These two schemes are applied in sequence. */
      target = target2 = -1;
      for scheme in (0)..(2) {
        for min in (0)..(ncand2) {
          other = cand2[min as usize].val;

          /* pid_from is the subdomain from where_ the vertices will be removed.
          pid_to is the adjacent subdomain to pid_from that defines the 
          (me, other) subdomain edge that needs to be removed */
          if (scheme == 0) {
            pid_from = other;
            pid_to   = me;
          }
          else {
            pid_from  = me;
            pid_to    = other;
          }

          /* Go and find the vertices in 'other' that are connected in 'me' */
          nind=0;
          for ii in (pptr[pid_from as usize])..(pptr[(pid_from+1) as usize]) {
            i = pind[ii as usize];
            assert!(where_[i as usize] == pid_from);
            for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
              if (where_[adjncy[j as usize] as usize] == pid_to) {
                ind[(nind) as usize] = i;nind+=1;
                break;
              }
            }
          }

          /* Go and construct the otherpmat to see where_ these nind vertices are 
          connected to */
          iset(ncon, 0, cpwgt);
          ncand=0;
          for ii in (0)..(nind) {
            i = ind[ii as usize];
            iaxpy(ncon, 1, vwgt+i*ncon, 1, cpwgt, 1);

            for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
              if ((k = where_[adjncy[j as usize] as usize]) == pid_from) {
                continue;
              }
              if (otherpmat[k as usize] == 0) {
                cand[(ncand) as usize].val = k;ncand+=1;
              }
              otherpmat[k as usize] += (adjwgt ? adjwgt[j as usize] : 1);
            }
          }

          for i in (0)..(ncand) {
            cand[i as usize].key = otherpmat[cand[i as usize].val];
            assert!(cand[i as usize].key > 0);
          }

          ikvsortd(ncand, cand);

          IFSET(ctrl.dbglvl, METIS_DBG_CONNINFO, 
            print!("\tMinOut: {:4}, to: {:3}, TtlWgt: {:5}[#:{:}]\n", 
              mypmat[other as usize], other, isum(ncon, cpwgt, 1), nind));

          /* Go through and select the first domain that is common with 'me', and does
          not increase the nads[target as usize] higher than nads[me as usize], subject to the maxpwgt
          constraint. Traversal is done from the mostly connected to the least. */
          for i in (0)..(ncand) {
            k = cand[i as usize].val;

            if (mypmat[k as usize] > 0) {
              /* Check if balance will go off */
              if (!ivecaxpylez(ncon, 1, cpwgt, pwgts+k*ncon, maxpwgt+k*ncon)) {
                continue;
              }

              /* get a dense vector out of k's connectivity */
              for j in (0)..(nads[k as usize])  {
                kpmat[adids[k as usize][j as usize] as usize] = adwgts[k as usize][j as usize];
              }

              /* Check if the move to domain k will increase the nads of another
              subdomain j that the set of vertices being moved are connected
              to but domain k is not connected to. */
              for j in (0)..(nparts) {
                if (otherpmat[j as usize] > 0 && kpmat[j as usize] == 0 && nads[j as usize]+1 >= nads[me as usize])  {
                  break;
                }
              }

              /* There were no bad second level effects. See if you can find a
              subdomain to move to. */
              if (j == nparts) { 
                nadd=0;
                for j in (0)..(nparts) {
                  if (otherpmat[j as usize] > 0 && kpmat[j as usize] == 0) {
                    nadd += 1;
                  }
                }

                IFSET(ctrl.dbglvl, METIS_DBG_CONNINFO, 
                  print!("\t\tto={:}, nadd={:}, {:}\n", k, nadd, nads[k as usize]));

                if (nads[k as usize]+nadd < nads[me as usize]) {
                  if (target2 == -1 || nads[target2 as usize]+bestnadd > nads[k as usize]+nadd ||
                  (nads[target2 as usize]+bestnadd == nads[k as usize]+nadd && bestnadd > nadd)) {
                    target2  = k;
                    bestnadd = nadd;
                  }
                }

                if (nadd == 0)  {
                  target = k;
                }
              }

              /* reset kpmat for the next iteration */
              for j in (0)..(nads[k as usize])  {
                kpmat[adids[k as usize][j as usize] as usize] = 0;
              }
            }

            if (target != -1) {
              break;
            }
          }

          /* reset the otherpmat for the next iteration */
          for i in (0)..(ncand)  {
            otherpmat[cand[i as usize].val] = 0;
          }

          if (target == -1 && target2 != -1) {
            target = target2;
          }

          if (target != -1) {
            IFSET(ctrl.dbglvl, METIS_DBG_CONNINFO, 
              print!("\t\tScheme: {:}. Moving to {:}\n", scheme, target));
            move_ = 1;
            break;
          }
        }

        if (target != -1) {
          break;
        }  /* A move was found. No need to try the other scheme */
      }

      /* reset the mypmat for next iteration */
      for i in (0)..(nads[me as usize])  {
        mypmat[adids[me as usize][i as usize] as usize] = 0;
      }

      /* Note that once a target is found the above loops exit right away. So the
      following variables are valid */
      if (target != -1) {
        match (ctrl.objtype) {
          METIS_OBJTYPE_CUT =>
          MoveGroupMinConnForCut(ctrl, graph, target, nind, ind)
          METIS_OBJTYPE_VOL =>
          MoveGroupMinConnForVol(ctrl, graph, target, nind, ind, vmarker, 
            pmarker, modind)
          _ =>
          gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype)
        }

        /* Update the csr representation of the partitioning vector */
        iarray2csr(nvtxs, nparts, where_, pptr, pind);
      }
    }

    if (move_ == 0) {
      break;
    }
  }

  ipqFree(&queue);

  // WCOREPOP;
}


/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupMinConnForCut(ctrl: *mut ctrl_t, graph: *mut graph_t, to: idx_t, nind: idx_t, ind: *mut idx_t) 
{
  // idx_t i, ii, j, jj, k, l, nvtxs, nbnd, from, me;
  // idx_t *xadj, *adjncy, *adjwgt, *where_, *bndptr, *bndind;
  // ckrinfo_t *myrinfo;
  // cnbr_t *mynbrs;

  let   nvtxs  = graph.nvtxs as usize;
  xadj   = graph.xadj;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  where_  = graph.where_;
  bndptr = graph.bndptr;
  bndind = graph.bndind;

  nbnd = graph.nbnd;

  while ({nind -= 1; nind }>=0) {
    i    = ind[nind as usize];
    from = where_[i as usize];

    myrinfo = graph.ckrinfo+i;
    if (myrinfo.inbr == -1) {
      myrinfo.inbr  = cnbrpoolGetNext(ctrl, xadj[(i+1) as usize]-xadj[i as usize]);
      myrinfo.nnbrs = 0;
    }
    mynbrs = ctrl.cnbrpool + myrinfo.inbr;

    /* find the location of 'to' in myrinfo or create it if it is not there */
    for k in (0)..(myrinfo.nnbrs) {
      if (mynbrs[k as usize].pid == to) {
        break;
      }
    }
    if (k == myrinfo.nnbrs) {
      assert!(k < xadj[(i+1) as usize]-xadj[i as usize]);
      mynbrs[k as usize].pid = to;
      mynbrs[k as usize].ed  = 0;
      myrinfo.nnbrs += 1;
    }

    /* Update pwgts */
    iaxpy(graph.ncon,  1, graph.vwgt+i*graph.ncon, 1, graph.pwgts+to*graph.ncon,   1);
    iaxpy(graph.ncon, -1, graph.vwgt+i*graph.ncon, 1, graph.pwgts+from*graph.ncon, 1);

    /* Update mincut */
    graph.mincut -= mynbrs[k as usize].ed-myrinfo.id;

    /* Update subdomain connectivity graph to reflect the move of 'i' */
    UpdateEdgeSubDomainGraph(ctrl, from, to, myrinfo.id-mynbrs[k as usize].ed, std::ptr::null_mut());

    /* Update ID/ED and BND related information for the moved vertex */
    UpdateMovedVertexInfoAndBND(i, from, k, to, myrinfo, mynbrs, where_, nbnd, 
      bndptr, bndind, BNDTYPE_REFINE);

    /* Update the degrees of adjacent vertices */
    for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
      ii = adjncy[j as usize];
      me = where_[ii as usize];
      myrinfo = graph.ckrinfo+ii;

      UpdateAdjacentVertexInfoAndBND(ctrl, ii, xadj[(ii+1) as usize]-xadj[ii as usize], me,
        from, to, myrinfo, adjwgt[j as usize], nbnd, bndptr, bndind, BNDTYPE_REFINE);

      /* Update subdomain graph to reflect the move of 'i' for domains other 
      than 'from' and 'to' */
      if (me != from && me != to) {
        UpdateEdgeSubDomainGraph(ctrl, from, me, -adjwgt[j as usize], std::ptr::null_mut());
        UpdateEdgeSubDomainGraph(ctrl, to, me, adjwgt[j as usize], std::ptr::null_mut());
      }
    }
  }

  assert!(ComputeCut(graph, where_) == graph.mincut);

  graph.nbnd = nbnd;

}


/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupMinConnForVol(ctrl: *mut ctrl_t, graph: *mut graph_t, to: idx_t, nind: idx_t, ind: *mut idx_t, vmarker: *mut idx_t, pmarker: *mut idx_t, modind: *mut idx_t) 
{
  // idx_t i, ii, j, jj, k, l, nvtxs, from, me, other, xgain, ewgt;
  // idx_t *xadj, *vsize, *adjncy, *where_;
  // vkrinfo_t *myrinfo, *orinfo;
  // vnbr_t *mynbrs, *onbrs;

  let   nvtxs  = graph.nvtxs as usize;
  xadj   = graph.xadj;
  vsize  = graph.vsize;
  adjncy = graph.adjncy;
  where_  = graph.where_;

  while ({nind -= 1; nind}>=0) {
    i    = ind[nind as usize];
    from = where_[i as usize];

    myrinfo = graph.vkrinfo+i;
    if (myrinfo.inbr == -1) {
      myrinfo.inbr  = vnbrpoolGetNext(ctrl, xadj[(i+1) as usize]-xadj[i as usize]);
      myrinfo.nnbrs = 0;
    }
    mynbrs = ctrl.vnbrpool + myrinfo.inbr;

    xgain = (myrinfo.nid == 0 && myrinfo.ned > 0 ? vsize[i as usize] : 0);

    //print!("Moving {:} from {:} to {:} [(vsize: {:}) as usize] [(xgain: {:}) as usize]\n", 
    //    i, from, to, vsize[i as usize], xgain);

    /* find the location of 'to' in myrinfo or create it if it is not there */
    for k in (0)..(myrinfo.nnbrs) {
      if (mynbrs[k as usize].pid == to) {
        break;
      }
    }

    if (k == myrinfo.nnbrs) {
      //print!("Missing neighbor\n");

      if (myrinfo.nid > 0) {
        xgain -= vsize[i as usize];
      }

      /* determine the volume gain resulting from that move */
      for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
        ii     = adjncy[j as usize];
        other  = where_[ii as usize];
        orinfo = graph.vkrinfo+ii;
        onbrs  = ctrl.vnbrpool + orinfo.inbr;
        assert!(other != to);

        //print!("  %8d %8d %3d\n", (int)ii, (int)vsize[ii as usize], (int)other);

        if (from == other) {
          /* Same subdomain vertex: Decrease the gain if 'to' is a new neighbor. */
          for l in (0)..(orinfo.nnbrs) {
            if (onbrs[l as usize].pid == to) {
              break;
            }
          }
          if (l == orinfo.nnbrs)  {
            xgain -= vsize[ii as usize];
          }
        }
        else {
          /* Remote vertex: increase if 'to' is a new subdomain */
          for l in (0)..(orinfo.nnbrs) {
            if (onbrs[l as usize].pid == to) {
              break;
            }
          }
          if (l == orinfo.nnbrs)  {
            xgain -= vsize[ii as usize];
          }

          /* Remote vertex: decrease if i is the only connection to 'from' */
          for l in (0)..(orinfo.nnbrs) {
            if (onbrs[l as usize].pid == from && onbrs[l as usize].ned == 1) {
              xgain += vsize[ii as usize];
              break;
            }
          }
        }
      }
      graph.minvol -= xgain;
      graph.mincut -= -myrinfo.nid;
      ewgt = myrinfo.nid;
    }
    else {
      graph.minvol -= (xgain + mynbrs[k as usize].gv);
      graph.mincut -= mynbrs[k as usize].ned-myrinfo.nid;
      ewgt = myrinfo.nid-mynbrs[k as usize].ned;
    }

    /* Update where_ and pwgts */
    where_[i as usize] = to;
    iaxpy(graph.ncon,  1, graph.vwgt+i*graph.ncon, 1, graph.pwgts+to*graph.ncon,   1);
    iaxpy(graph.ncon, -1, graph.vwgt+i*graph.ncon, 1, graph.pwgts+from*graph.ncon, 1);

    /* Update subdomain connectivity graph to reflect the move of 'i' */
    UpdateEdgeSubDomainGraph(ctrl, from, to, ewgt, std::ptr::null_mut());

    /* Update the subdomain connectivity of the adjacent vertices */
    for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
      me = where_[adjncy[j as usize] as usize];
      if (me != from && me != to) {
        UpdateEdgeSubDomainGraph(ctrl, from, me, -1, std::ptr::null_mut());
        UpdateEdgeSubDomainGraph(ctrl, to, me, 1, std::ptr::null_mut());
      }
    }

    /* Update the id/ed/gains/bnd of potentially affected nodes */
    KWayVolUpdate(ctrl, graph, i, from, to, std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
      std::ptr::null_mut(), BNDTYPE_REFINE, vmarker, pmarker, modind);

    /*CheckKWayVolPartitionParams(ctrl, graph);*/
  }
  assert!(ComputeCut(graph, where_) == graph.mincut);
  assert!(ComputeVolume(graph, where_) == graph.minvol, "{:} {:}\n", ComputeVolume(graph, where_), graph.minvol);

}


/*************************************************************************/
/* This function computes the subdomain graph. For deubugging purposes. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn PrintSubDomainGraph(graph: *mut graph_t, nparts: idx_t, where_: *mut idx_t) 
{
  // idx_t i, j, k, me, nvtxs, total, max;
  // idx_t *xadj, *adjncy, *adjwgt, *pmat;

  let   nvtxs  = graph.nvtxs as usize;
  xadj   = graph.xadj;
  adjncy = graph.adjncy;
  adjwgt = graph.adjwgt;

  pmat = ismalloc(nparts*nparts, 0, "ComputeSubDomainGraph: pmat");

  for i in (0)..(nvtxs) {
    me = where_[i as usize];
    for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
      k = adjncy[j as usize] as usize;
      if (where_[k as usize] != me)  {
        pmat[(me*nparts+where_[k]) as usize] += adjwgt[j as usize];
      }
    }
  }

  /* print!("Subdomain Info\n"); */
  total = max = 0;
  for i in (0)..(nparts) {
    k=0;
    for j in (0)..(nparts) {
      if (pmat[(i*nparts+j) as usize] > 0) {
        k += 1;
      }
    }
    total += k;

    if (k > max) {
      max = k;
    }
    /*
    print!("{:2} . {:2}  ", i, k);
    for j in (0)..(nparts) {
    if (pmat[(i*nparts+j) as usize] > 0) {
    print!("[({:2} {:4}) as usize] ", j, pmat[(i*nparts+j) as usize]);
    }
    }
    print!("\n");
    */
  }
  print!("Total adjacent subdomains: {:}, Max: {:}\n", total, max);

  gk_free((void **)&pmat, LTERM);
}


