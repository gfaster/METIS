#include "metislib.h"

/*************************************************************************/
/*! This function selects the partition number and the queue from which
    we will move vertices out. */
/*************************************************************************/ 
// PORTING: since we use a totally different system for the Rust version, this is not ifunc
// we keep this as a fake unported so we can switch versions
void SelectQueue( graph_t *graph, real_t *pijbm, real_t *ubfactors,
                              rpq_t **queues, idx_t *from, idx_t *cnum)
{
  idx_t ncon, i, part;
  real_t max, tmp;

  ncon = graph->ncon;

  *from = -1;
  *cnum = -1;

  /* First determine the side and the queue, irrespective of the presence of nodes. 
     The side & queue is determined based on the most violated balancing constraint. */
  for (max=0.0, part=0; part<2; part++) {
    for (i=0; i<ncon; i++) {
      tmp = graph->pwgts[part*ncon+i]*pijbm[part*ncon+i] - ubfactors[i];
      /* the '=' in the test below is to ensure that under tight constraints
         the partition that is at the max is selected */
      if (tmp >= max) { 
        max   = tmp;
        *from = part;
        *cnum = i;
      }
    }
  }


  if (*from != -1) {
    /* in case the desired queue is empty, select a queue from the same side */
    if (rpqLength(queues[2*(*cnum)+(*from)]) == 0) {
      for (i=0; i<ncon; i++) {
        if (rpqLength(queues[2*i+(*from)]) > 0) {
          max   = graph->pwgts[(*from)*ncon+i]*pijbm[(*from)*ncon+i] - ubfactors[i];
          *cnum = i;
          break;
        }
      }

      for (i++; i<ncon; i++) {
        tmp = graph->pwgts[(*from)*ncon+i]*pijbm[(*from)*ncon+i] - ubfactors[i];
        if (tmp > max && rpqLength(queues[2*i+(*from)]) > 0) {
          max   = tmp;
          *cnum = i;
        }
      }
    }

    /*
    printf("Selected1 %"PRIDX"(%"PRIDX") -> %"PRIDX" [%5"PRREAL"]\n", 
        *from, *cnum, rpqLength(queues[2*(*cnum)+(*from)]), max); 
    */
  }
  else {
    /* the partitioning does not violate balancing constraints, in which case select 
       a queue based on cut criteria */
    for (part=0; part<2; part++) {
      for (i=0; i<ncon; i++) {
        if (rpqLength(queues[2*i+part]) > 0 && 
            (*from == -1 || rpqSeeTopKey(queues[2*i+part]) > max)) {
          max   = rpqSeeTopKey(queues[2*i+part]); 
          *from = part;
          *cnum = i;
        }
      }
    }
    /*
    printf("Selected2 %"PRIDX"(%"PRIDX") -> %"PRIDX"\n", 
        *from, *cnum, rpqLength(queues[2*(*cnum)+(*from)]), max); 
    */
  }
}
