/*
\file  gk_mkpqueue.h
\brief Templates for priority queues

\date   Started 4/09/07
\author George
\version\verbatim $Id: gk_mkpqueue.h 21742 2018-01-26 16:59:15Z karypis $ \endverbatim
*/



/*************************************************************************/
/* This function initializes the data structures of the priority queue */
/**************************************************************************/
pub fn init( maxnodes: usize) -> Self
{
  Self {
  nnodes: 0,
  maxnodes,

  heap   : KVMALLOC(maxnodes, "gk_PQInit: heap"),
  locator: gk_idxsmalloc(maxnodes, -1, "gk_PQInit: locator"),
}
}


/*************************************************************************/
/* This function resets the priority queue */
/**************************************************************************/
pub fn reset(&mut self)
{
  ssize_t i;
  ssize_t *locator=self.locator;
  KVT *heap=self.heap;

  for (i=self.nnodes-1; i>=0; i -= 1)
    locator[heap[i].val] = -1;
  self.nnodes = 0;
}


/*************************************************************************/
/* This function returns the length of the queue */
/**************************************************************************/
pub fn length(&mut self) -> usize
{
  return self.nnodes;
}


/*************************************************************************/
/* This function adds an item in the priority queue */
/**************************************************************************/
pub fn insert(&mut self, node: VT, key: KT) -> i32
{
  let locator = &mut self.locator;
  let heap = &mut self.heap;

  debug_assert!(Self::CheckHeap(queue));

  debug_assert!(locator[node] == -1);

  i = self.nnodes += 1;
  while (i > 0) {
    j = (i-1)>>1;
    if key < heap[j].key {
      heap[i] = heap[j];
      locator[heap[i].val] = i;
      i = j;
    }
    else
    {
      break;
    }
  }
  debug_assert!(i >= 0);
  heap[i].key   = key;
  heap[i].val   = node;
  locator[node] = i;

  debug_assert!(self.CheckHeap());

  return 0;
}


/*************************************************************************/
/* This function deletes an item from the priority queue */
/**************************************************************************/
pub fn delete(&mut self, node: VT) -> i32
{
  let locator=&mut self.locator;
  let heap=&mut self.heap;

  debug_assert!(locator[node] != -1);
  debug_assert!(heap[locator[node]].val == node);

  debug_assert!(Self::CheckHeap(queue));

  i = locator[node];
  locator[node] = -1;

  self.nnodes -= 1;
  if self.nnodes > 0 && heap[self.nnodes].val != node {
    node   = heap[self.nnodes].val;
    newkey = heap[self.nnodes].key;
    oldkey = heap[i].key;

    if newkey < oldkey { /* Filter-up */
      while (i > 0) {
        j = (i-1)>>1;
        if newkey < heap[j].key {
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
          break;
      }
    }
    else { /* Filter down */
      nnodes = self.nnodes;
      while ((j=(i<<1)+1) < nnodes) {
        if heap[j].key < newkey {
          if j+1 < nnodes && heap[j+1].key < heap[j].key
          {
            j += 1;
          }
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else if j+1 < nnodes && heap[j+1].key < newkey {
          j += 1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
        {
          break;
        }
      }
    }

    heap[i].key   = newkey;
    heap[i].val   = node;
    locator[node] = i;
  }

  debug_assert!(Self::CheckHeap(queue));

  return 0;
}


/*************************************************************************/
/* This function updates the key values associated for a particular item */ 
/**************************************************************************/
pub fn update(&mut self, node: VT, newkey: KT)
{
  let locator=self.locator;
  let heap=self.heap;

  let oldkey = heap[locator[node]].key;
  if !(newkey < oldkey) && !(oldkey < newkey) {return;}

  debug_assert!(locator[node] != -1);
  debug_assert!(heap[locator[node]].val == node);
  debug_assert!(self.CheckHeap());

  i = locator[node];

  if newkey < oldkey { /* Filter-up */
    while (i > 0) {
      j = (i-1)>>1;
      if newkey < heap[j].key {
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
      {
        break;
      }
    }
  }
  else { /* Filter down */
    nnodes = self.nnodes;
    while ((j=(i<<1)+1) < nnodes) {
      if heap[j].key < newkey {
        if j+1 < nnodes && heap[j+1].key < heap[j].key
        {
          j += 1;
        }
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else if j+1 < nnodes && heap[j+1].key < newkey {
        j += 1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
      {
        break;
      }
    }
  }

  heap[i].key   = newkey;
  heap[i].val   = node;
  locator[node] = i;

  debug_assert!(Self::CheckHeap(queue));

  return;
}


/*************************************************************************/
/* This function returns the item at the top of the queue and removes
    it from the priority queue */
/**************************************************************************/
pub fn get_top(&mut self) -> VT
{

  debug_assert!(Self::CheckHeap(queue));

  if self.nnodes == 0
  {
    return -1;
  }

  self.nnodes -= 1;

  let heap    = &mut self.heap;
  let locator = &mut self.locator;

  let mut vtx = heap[0].val;
  locator[vtx] = -1;

  i = self.nnodes;
  if i > 0 {
    let key  = heap[i].key;
    let node = heap[i].val;
    i = 0;
    while ((j=2*i+1) < self.nnodes) {
      if heap[j].key < key {
        if j+1 < self.nnodes && heap[j+1].key < heap[j].key
          {
            j = j+1;
          }
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else if j+1 < self.nnodes && heap[j+1].key < key {
        j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
        {
          break;
        }
      }

    heap[i].key   = key;
    heap[i].val   = node;
    locator[node] = i;
  }

  debug_assert!(self.CheckHeap());
  return vtx;
}


/*************************************************************************/
/* This function returns the item at the top of the queue. The item is not
    deleted from the queue. */
/**************************************************************************/
pub fn see_top_val(&mut self) -> VT
{
  if self.nnodes == 0 {
    -1
  } else {
    self.heap[0].val
  }
}


/*************************************************************************/
/* This function returns the key of the top item. The item is not
    deleted from the queue. */
/**************************************************************************/
pub fn see_top_key(&mut self) -> KT
{
  if self.nnodes == 0 {
    KMAX
  } else {
    self.heap[0].key
  }
}


/*************************************************************************/
/* This function returns the key of a specific item */
/**************************************************************************/
pub fn see_key(&self, node: VT) -> KT
{
  return self.heap[self.locator[node]].key;
}


/*************************************************************************/
/* This functions checks the consistency of the heap */
/**************************************************************************/
pub fn check_heap(&mut self) -> i32
{
  let heap    = &mut self.heap;
  let locator = &mut self.locator;
  let nnodes  = &mut self.nnodes;

  if nnodes == 0
{
  return 1;
}

  debug_assert!(locator[heap[0].val] == 0);
  for i in 1..nnodes {
    debug_assert!(locator[heap[i].val] == i);
    debug_assert!(!heap[i].key < heap[(i-1)/2].key);
  }
  for i in 1..nnodes
  {
    debug_assert!(!heap[i].key < heap[0].key);
  }

  let mut j = 0;
  for i in 0..self.maxnodes {
    if locator[i] != -1
    {
      j += 1;
    }
  }
  debug_assert!(j == nnodes, "{} {}", j, nnodes);

  return 1;
}
