#![allow(dead_code)]
/*
\file  gk_mkpqueue.h
\brief Templates for priority queues

\date   Started 4/09/07
\author George
\version\verbatim $Id: gk_mkpqueue.h 21742 2018-01-26 16:59:15Z karypis $ \endverbatim
*/
use std::fmt::Debug;

use crate::{idx_t, real_t};

pub type IPQueue = IndexedPriorityQueue<idx_t, idx_t>;
pub type RPQueue = IndexedPriorityQueue<real_t, idx_t>;

macro_rules! idx {
    ($key:expr) => {
        TryInto::<usize>::try_into($key).expect("valid index")
    };
}

#[derive(Default, Clone, Copy)]
struct Node<K, V>
where
    V: Copy + Default + PartialEq,
    usize: TryFrom<V>,
    <usize as TryFrom<V>>::Error: Debug,
    K: PartialOrd + Copy + Default,
{
    key: K,
    val: V,
}

/// Priority queue augmented with indices. After inserting a (key, index) pair, it can be quickly
/// accessed and/or removed using the index.
pub struct IndexedPriorityQueue<K, V>
where
    V: Copy + Default + PartialEq,
    usize: TryFrom<V>,
    <usize as TryFrom<V>>::Error: Debug,
    K: PartialOrd + Copy + Default,
{
    nnodes: usize,
    locator: Box<[isize]>,
    heap: Box<[Node<K, V>]>,
}

impl<K, V> IndexedPriorityQueue<K, V>
where
    V: Copy + Default + PartialEq,
    usize: TryFrom<V>,
    <usize as TryFrom<V>>::Error: Debug,
    K: PartialOrd + Copy + Default,
{
    /*************************************************************************/
    /* This function initializes the data structures of the priority queue */
    /**************************************************************************/
    pub fn new(maxnodes: usize) -> Self {
        Self {
            nnodes: 0,

            heap: Vec::from_iter((0..maxnodes).map(|_| Default::default())).into_boxed_slice(),
            locator: vec![-1; maxnodes].into_boxed_slice(),
        }
    }

    /// clear the queue
    pub fn reset(&mut self) {
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        // for (i=self.nnodes-1; i>=0; i -= 1)
        for i in (0..self.nnodes).rev() {
            locator[idx!(heap[i].val)] = -1;
        }
        self.nnodes = 0;
    }

    /// get the length
    pub fn length(&mut self) -> usize {
        self.nnodes
    }

    /// insert an item
    pub fn insert(&mut self, index: V, key: K) {
        debug_assert!(self.check_heap());

        let locator = &mut self.locator;
        let heap = &mut self.heap;

        debug_assert!(locator[idx!(index)] == -1);

        let mut i = self.nnodes;
        self.nnodes += 1;
        while i > 0 {
            let j = (i - 1) >> 1;
            if key < heap[j].key {
                heap[i] = heap[j];
                locator[idx!(heap[i].val)] = i as isize;
                i = j;
            } else {
                break;
            }
        }
        // debug_assert!(i >= 0);
        heap[i].key = key;
        heap[i].val = index;
        locator[idx!(index)] = i as isize;

        debug_assert!(self.check_heap());
    }

    /// delete an item
    pub fn delete(&mut self, index: V) {
        debug_assert!(self.check_heap());
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        debug_assert!(locator[idx!(index)] != -1);
        debug_assert!(heap[locator[idx!(index)] as usize].val == index);

        let mut i = locator[idx!(index)];
        locator[idx!(index)] = -1;

        self.nnodes -= 1;
        if self.nnodes > 0 && heap[self.nnodes].val != index {
            let node = heap[self.nnodes].val;
            let newkey = heap[self.nnodes].key;
            let oldkey = heap[i as usize].key;

            if newkey < oldkey {
                /* Filter-up */
                while i > 0 {
                    let j = (i - 1) >> 1;
                    if newkey < heap[j as usize].key {
                        heap[i as usize] = heap[j as usize];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val)
                            .expect("valid index")] = i;
                        i = j;
                    } else {
                        break;
                    }
                }
            } else {
                /* Filter down */
                let nnodes = self.nnodes;
                loop {
                    let mut j = ((i as usize) << 1) + 1;
                    if j >= nnodes {
                        break;
                    }
                    if heap[j].key < newkey {
                        if j + 1 < nnodes && heap[j + 1].key < heap[j].key {
                            j += 1;
                        }
                        heap[i as usize] = heap[j];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val)
                            .expect("valid index")] = i;
                        i = j as isize;
                    } else if j + 1 < nnodes && heap[j + 1].key < newkey {
                        j += 1;
                        heap[i as usize] = heap[j];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val)
                            .expect("valid index")] = i;
                        i = j as isize;
                    } else {
                        break;
                    }
                }
            }

            heap[i as usize].key = newkey;
            heap[i as usize].val = node;
            locator[idx!(node)] = i;
        }

        debug_assert!(self.check_heap());
    }

    /*************************************************************************/
    /* This function updates the key values associated for a particular item */
    /**************************************************************************/
    pub fn update(&mut self, node: V, newkey: K) {
        debug_assert!(self.check_heap());
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        let oldkey = heap[locator[idx!(node)] as usize].key;
        if (newkey >= oldkey) && (oldkey >= newkey) {
            return;
        }

        debug_assert!(locator[idx!(node)] != -1);
        debug_assert!(heap[locator[idx!(node)] as usize].val == node);

        let mut i = locator[idx!(node)];

        if newkey < oldkey {
            /* Filter-up */
            while i > 0 {
                let j = (i as usize - 1) >> 1;
                if newkey < heap[j].key {
                    heap[i as usize] = heap[j];
                    locator[idx!(heap[i as usize].val)] = i;
                    i = j as isize;
                } else {
                    break;
                }
            }
        } else {
            /* Filter down */
            let nnodes = self.nnodes;
            loop {
                let mut j = ((i as usize) << 1) + 1;
                if j >= nnodes {
                    break;
                }
                if heap[j].key < newkey {
                    if j + 1 < nnodes && heap[j + 1].key < heap[j].key {
                        j += 1;
                    }
                    heap[i as usize] = heap[j];
                    locator[idx!(heap[i as usize].val)] = i;
                    i = j as isize;
                } else if j + 1 < nnodes && heap[j + 1].key < newkey {
                    j += 1;
                    heap[i as usize] = heap[j];
                    locator[idx!(heap[i as usize].val)] = i;
                    i = j as isize;
                } else {
                    break;
                }
            }
        }

        heap[i as usize].key = newkey;
        heap[i as usize].val = node;
        locator[idx!(node)] = i;

        debug_assert!(self.check_heap());
    }

    /// This function returns the item at the top of the queue and removes it from the priority queue
    /// I have not yet convinced myself it will never return -1
    pub fn get_top(&mut self) -> Option<V> {
        debug_assert!(self.check_heap());

        if self.nnodes == 0 {
            return None;
        }

        self.nnodes -= 1;

        let heap = &mut self.heap;
        let locator = &mut self.locator;

        let vtx = heap[0].val;
        locator[idx!(vtx)] = -1;

        let mut i = self.nnodes;
        if i > 0 {
            let key = heap[i].key;
            let node = heap[i].val;
            i = 0;
            loop {
                let mut j = 2 * i + 1;
                if j >= self.nnodes {
                    break;
                }
                if heap[j].key < key {
                    if j + 1 < self.nnodes && heap[j + 1].key < heap[j].key {
                        j += 1;
                    }
                    heap[i] = heap[j];
                    locator[idx!(heap[i].val)] = i as isize;
                    i = j;
                } else if j + 1 < self.nnodes && heap[j + 1].key < key {
                    j += 1;
                    heap[i] = heap[j];
                    locator[idx!(heap[i].val)] = i as isize;
                    i = j;
                } else {
                    break;
                }
            }

            heap[i].key = key;
            heap[i].val = node;
            locator[idx!(node)] = i as isize;
        }

        debug_assert!(self.check_heap());
        Some(vtx)
    }

    /*************************************************************************/
    /* This function returns the item at the top of the queue. The item is not
    deleted from the queue. */
    /**************************************************************************/
    pub fn see_top_val(&mut self) -> Option<V> {
        if self.nnodes == 0 {
            None
        } else {
            Some(self.heap[0].val)
        }
    }

    /*************************************************************************/
    /* This function returns the key of the top item. The item is not
    deleted from the queue. */
    /**************************************************************************/
    pub fn see_top_key(&mut self) -> Option<K> {
        if self.nnodes == 0 {
            None
        } else {
            Some(self.heap[0].key)
        }
    }

    /*************************************************************************/
    /* This function returns the key of a specific item */
    /**************************************************************************/
    pub fn see_key(&self, node: V) -> K {
        self.heap[self.locator[idx!(node)] as usize].key
    }

    /*************************************************************************/
    /* This functions checks the consistency of the heap */
    /**************************************************************************/
    pub fn check_heap(&self) -> bool {
        let heap = &self.heap;
        let locator = &self.locator;
        let nnodes = self.nnodes;

        if nnodes == 0 {
            return true;
        }

        assert!(locator[idx!(heap[0].val)] == 0);
        for i in 1..nnodes {
            assert!(locator[idx!(heap[i].val)] == i as isize);
            assert!(heap[i].key >= heap[(i - 1) / 2].key);
        }
        for i in 1..nnodes {
            assert!(heap[i].key >= heap[0].key);
        }

        let mut j = 0;
        for i in 0..self.heap.len() {
            if locator[i] != -1 {
                j += 1;
            }
        }
        assert!(j == nnodes, "{} {}", j, nnodes);

        true
    }
}

/// This function selects the partition number and the queue from which we will move vertices out.
///
/// this is originally from fm.c, but I needed to rewrite it early
pub fn select_queue(
    graph: &mut crate::graph_t,
    pijbm: &[real_t],
    ubfactors: &[real_t],
    queues: &mut [RPQueue],
    from: &mut idx_t,
    cnum: &mut idx_t,
) {
    // idx_t ncon, i, part;
    // real_t max, tmp;

    let ncon = graph.ncon;
    let pwgts: &[_] = unsafe { std::slice::from_raw_parts(graph.pwgts, (ncon * 2) as usize) };

    *from = -1;
    *cnum = -1;

    /* First determine the side and the queue, irrespective of the presence of nodes.
    The side & queue is determined based on the most violated balancing constraint. */
    let mut max = 0.0;
    for part in 0..2 {
        for i in 0..ncon {
            let tmp = pwgts[(part * ncon + i) as usize] as real_t
                * pijbm[(part * ncon + i) as usize]
                - ubfactors[i as usize];
            /* the '=' in the test below is to ensure that under tight constraints
            the partition that is at the max is selected */
            if tmp >= max {
                max = tmp;
                *from = part;
                *cnum = i;
            }
        }
    }

    if *from != -1 {
        /* in case the desired queue is empty, select a queue from the same side */
        if (queues[(2 * *cnum + *from) as usize]).length() == 0 {
            let mut i = 0;
            for ii in 0..ncon {
                i = ii;
                if (queues[(2 * i + *from) as usize]).length() > 0 {
                    max = pwgts[(*from * ncon + i) as usize] as real_t
                        * pijbm[(*from * ncon + i) as usize]
                        - ubfactors[(i) as usize];
                    *cnum = i;
                    break;
                }
            }
            i += 1;
            // for (i++; i<ncon; i++) {
            while i < ncon {
                let tmp = pwgts[(*from * ncon + i) as usize] as real_t
                    * pijbm[(*from * ncon + i) as usize]
                    - ubfactors[i as usize];
                if tmp > max && (queues[(2 * i + *from) as usize]).length() > 0 {
                    max = tmp;
                    *cnum = i;
                }
                i += 1;
            }
        }

        /*
        printf("Selected1 %"PRIDX"(%"PRIDX") . %"PRIDX" [%5"PRREAL"]\n",
            *from, *cnum, rpqLength(queues[2*(*cnum)+(*from)]), max);
        */
    } else {
        /* the partitioning does not violate balancing constraints, in which case select
        a queue based on cut criteria */
        for part in 0..2 {
            for i in 0..ncon {
                if (queues[(2 * i + part) as usize]).length() > 0
                    && (*from == -1
                        || (queues[(2 * i + part) as usize].see_top_key().unwrap()) > max)
                {
                    max = queues[(2 * i + part) as usize].see_top_key().unwrap();
                    *from = part;
                    *cnum = i;
                }
            }
        }
        /*
        printf("Selected2 %"PRIDX"(%"PRIDX") . %"PRIDX"\n",
            *from, *cnum, rpqLength(queues[2*(*cnum)+(*from)]), max);
        */
    }
}

#[cfg(test)]
mod test {
    use super::IndexedPriorityQueue;
    use std::collections::BinaryHeap;

    #[test]
    fn in_order() {
        let mut heap = IndexedPriorityQueue::new(10);
        for x in 0..10 {
            heap.insert(x, x);
        }
        let mut it = 0..10;
        while let Some(k) = heap.get_top() {
            assert_eq!(Some(k), it.next());
        }
        assert_eq!(it.next(), None);
    }

    #[test]
    fn random_order() {
        let mut rand = fastrand::Rng::with_seed(1);

        let mut heap = IndexedPriorityQueue::new(100);
        let mut truth = BinaryHeap::new();

        let mut items = Vec::from_iter(0..100);
        rand.shuffle(&mut items);
        for x in 0..100 {
            let k = items[x as usize];
            heap.insert(x, k);
            truth.push((99 - k, k)); // std::heap is min-heap
        }

        for _ in 0..80 {
            eprintln!(
                "gk: {:>4?}  |  std: {:>4?}",
                heap.see_top_key(),
                truth.peek().map(|x| x.1)
            );
            assert_eq!(heap.see_top_key(), truth.peek().map(|x| x.1));
            heap.get_top();
            truth.pop();
        }
    }

    #[test]
    fn random_order_with_removes() {
        let mut rand = fastrand::Rng::with_seed(1);

        let mut heap = IndexedPriorityQueue::new(100);
        let mut truth = BinaryHeap::new();

        let mut items = Vec::from_iter(0..100);
        rand.shuffle(&mut items);
        for x in 0..35 {
            let k = items[x as usize];
            heap.insert(x, k);
        }
        for x in 35..100 {
            let k = items[x as usize];
            heap.insert(x, k);
            truth.push((99 - k, k));
        }
        for x in 0..35 {
            heap.delete(x);
        }

        for _ in 0..80 {
            eprintln!(
                "gk: {:>4?}  |  std: {:>4?}",
                heap.see_top_key(),
                truth.peek().map(|x| x.1)
            );
            assert_eq!(heap.see_top_key(), truth.peek().map(|x| x.1));
            heap.get_top();
            truth.pop();
        }
    }
}
