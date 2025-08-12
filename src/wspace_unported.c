/*!
\file 
\brief Functions dealing with memory allocation and workspace management

\date Started 2/24/96
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: wspace.c 10492 2011-07-06 09:28:42Z karypis $
*/


// None of these will be ported to Rust

#include "metislib.h"


/*************************************************************************/
/*! This function allocate space from the workspace/heap */
/*************************************************************************/
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes)
{
  return gk_mcoreMalloc(ctrl->mcore, nbytes);
}


/*************************************************************************/
/*! This function sets a marker in the stack of malloc ops to be used
    subsequently for freeing purposes */
/*************************************************************************/
void wspacepush(ctrl_t *ctrl)
{
  gk_mcorePush(ctrl->mcore);
}


/*************************************************************************/
/*! This function frees all mops since the last push */
/*************************************************************************/
void wspacepop(ctrl_t *ctrl)
{
  gk_mcorePop(ctrl->mcore);
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
idx_t *iwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (idx_t *)wspacemalloc(ctrl, n*sizeof(idx_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
real_t *rwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (real_t *)wspacemalloc(ctrl, n*sizeof(real_t));
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (ikv_t *)wspacemalloc(ctrl, n*sizeof(ikv_t));
}

