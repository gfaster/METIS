/*!
\file 
\brief Functions dealing with memory allocation and workspace management

\date Started 2/24/96
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: wspace.c 10492 2011-07-06 09:28:42Z karypis $
*/

#include "metislib.h"
#include "ifunc.h"

// some of the functions in the original wspace.c will remain unported and are
// instead moved to wspace_unported.c

/*************************************************************************/
/*! This function allocates memory for the workspace */
/*************************************************************************/
IFUNC(void, AllocateWorkSpace, (ctrl_t *ctrl, graph_t *graph));
void c__libmetis__AllocateWorkSpace(ctrl_t *ctrl, graph_t *graph)
{
  size_t coresize;

  switch (ctrl->optype) {
    case METIS_OP_PMETIS:
      coresize = 3*(graph->nvtxs+1)*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(real_t);
      break;
    default:
      coresize = 4*(graph->nvtxs+1)*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(real_t);
  }
  ctrl->mcore = gk_mcoreCreate(coresize);

  ctrl->nbrpoolsize = 0;
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function allocates refinement-specific memory for the workspace */
/*************************************************************************/
IFUNC(void, AllocateRefinementWorkSpace, (ctrl_t *ctrl, idx_t nbrpoolsize_max, idx_t nbrpoolsize));
void c__libmetis__AllocateRefinementWorkSpace(ctrl_t *ctrl, idx_t nbrpoolsize_max, idx_t nbrpoolsize)
{
  ctrl->nbrpoolsize_max = nbrpoolsize_max;
  ctrl->nbrpoolsize     = nbrpoolsize;
  ctrl->nbrpoolcpos     = 0;
  ctrl->nbrpoolreallocs = 0;

  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      ctrl->cnbrpool = (cnbr_t *)gk_malloc(ctrl->nbrpoolsize*sizeof(cnbr_t), 
                             "AllocateRefinementWorkSpace: cnbrpool");
      break;

    case METIS_OBJTYPE_VOL:
      ctrl->vnbrpool = (vnbr_t *)gk_malloc(ctrl->nbrpoolsize*sizeof(vnbr_t), 
                             "AllocateRefinementWorkSpace: vnbrpool");
      break;

    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }


  /* Allocate the memory for the sparse subdomain graph */
  if (ctrl->minconn) {
    ctrl->pvec1   = imalloc(ctrl->nparts+1, "AllocateRefinementWorkSpace: pvec1");
    ctrl->pvec2   = imalloc(ctrl->nparts+1, "AllocateRefinementWorkSpace: pvec2");
    ctrl->maxnads = ismalloc(ctrl->nparts, INIT_MAXNAD, "AllocateRefinementWorkSpace: maxnads");
    ctrl->nads    = imalloc(ctrl->nparts, "AllocateRefinementWorkSpace: nads");
    ctrl->adids   = iAllocMatrix(ctrl->nparts, INIT_MAXNAD, 0, "AllocateRefinementWorkSpace: adids");
    ctrl->adwgts  = iAllocMatrix(ctrl->nparts, INIT_MAXNAD, 0, "AllocateRefinementWorkSpace: adwgts");
  }
}


/*************************************************************************/
/*! This function frees the workspace */
/*************************************************************************/
IFUNC(void, FreeWorkSpace, (ctrl_t *ctrl));
void c__libmetis__FreeWorkSpace(ctrl_t *ctrl)
{
  gk_mcoreDestroy(&ctrl->mcore, ctrl->dbglvl&METIS_DBG_INFO);

  IFSET(ctrl->dbglvl, METIS_DBG_INFO,
      printf(" nbrpool statistics\n" 
             "        nbrpoolsize: %12zu   nbrpoolcpos: %12zu\n"
             "    nbrpoolreallocs: %12zu\n\n",
             ctrl->nbrpoolsize,  ctrl->nbrpoolcpos, 
             ctrl->nbrpoolreallocs));

  if (ctrl->cnbrpool == (cnbr_t *)ctrl->vnbrpool) {
    ctrl->vnbrpool = NULL;
  }
  gk_free((void **)&ctrl->cnbrpool, &ctrl->vnbrpool, LTERM);
  ctrl->nbrpoolsize_max = 0;
  ctrl->nbrpoolsize     = 0;
  ctrl->nbrpoolcpos     = 0;

  if (ctrl->minconn) {
    iFreeMatrix(&(ctrl->adids),  ctrl->nparts, INIT_MAXNAD);
    iFreeMatrix(&(ctrl->adwgts), ctrl->nparts, INIT_MAXNAD);

    gk_free((void **)&ctrl->pvec1, &ctrl->pvec2, 
        &ctrl->maxnads, &ctrl->nads, LTERM);
  }
}



/*************************************************************************/
/*! This function resets the cnbrpool */
/*************************************************************************/
IFUNC(void, cnbrpoolReset, (ctrl_t *ctrl));
void c__libmetis__cnbrpoolReset(ctrl_t *ctrl)
{
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function gets the next free index from cnbrpool */
/*************************************************************************/
IFUNC(idx_t, cnbrpoolGetNext, (ctrl_t *ctrl, idx_t nnbrs));
idx_t c__libmetis__cnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs)
{
  /* add 1 because when moving vertices, an extra neighbor can be temporarily
   * needed (particularly when minconn is set) */
  nnbrs = gk_min(ctrl->nparts, nnbrs) + 1;
  ctrl->nbrpoolcpos += nnbrs;

  if (ctrl->nbrpoolcpos > ctrl->nbrpoolsize) {
    ctrl->nbrpoolsize += gk_max(10*nnbrs, ctrl->nbrpoolsize/2);
    ctrl->nbrpoolsize = gk_min(ctrl->nbrpoolsize, ctrl->nbrpoolsize_max);

    ctrl->cnbrpool = (cnbr_t *)gk_realloc(ctrl->cnbrpool,  
                          ctrl->nbrpoolsize*sizeof(cnbr_t), "cnbrpoolGet: cnbrpool");
    ctrl->nbrpoolreallocs++;

    ASSERT(ctrl->nbrpoolcpos <= ctrl->nbrpoolsize);
  }

  return ctrl->nbrpoolcpos - nnbrs;
}


/*************************************************************************/
/*! This function resets the vnbrpool */
/*************************************************************************/
IFUNC(void, vnbrpoolReset, (ctrl_t *ctrl));
void c__libmetis__vnbrpoolReset(ctrl_t *ctrl)
{
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function gets the next free index from vnbrpool */
/*************************************************************************/
IFUNC(idx_t, vnbrpoolGetNext, (ctrl_t *ctrl, idx_t nnbrs));
idx_t c__libmetis__vnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs)
{
  /* add 1 because when moving vertices, an extra neighbor can be temporarily
   * needed (particularly when minconn is set) */
  nnbrs = gk_min(ctrl->nparts, nnbrs) + 1;
  ctrl->nbrpoolcpos += nnbrs;

  if (ctrl->nbrpoolcpos > ctrl->nbrpoolsize) {
    ctrl->nbrpoolsize += gk_max(10*nnbrs, ctrl->nbrpoolsize/2);
    ctrl->nbrpoolsize = gk_min(ctrl->nbrpoolsize, ctrl->nbrpoolsize_max);

    ctrl->vnbrpool = (vnbr_t *)gk_realloc(ctrl->vnbrpool,  
                          ctrl->nbrpoolsize*sizeof(vnbr_t), "vnbrpoolGet: vnbrpool");
    ctrl->nbrpoolreallocs++;

    ASSERT(ctrl->nbrpoolcpos <= ctrl->nbrpoolsize);
  }

  return ctrl->nbrpoolcpos - nnbrs;
}

