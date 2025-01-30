#![allow(dead_code)]
/*
\file  gk_mkpqueue.h
\brief Templates for priority queues

\date   Started 4/09/07
\author George
\version\verbatim $Id: gk_mkpqueue.h 21742 2018-01-26 16:59:15Z karypis $ \endverbatim
*/
use std::{cmp::Reverse, collections::VecDeque, fmt::Debug};

use crate::{dal, idx_t, real_t};

pub type IPQueue = IndexedPriorityQueue<idx_t, idx_t>;
pub type RPQueue = IndexedPriorityQueue<real_t, idx_t>;

pub trait Idx: Copy + Default + PartialEq + Eq {
    fn idx(self) -> usize;
    fn uidx(idx: usize) -> Self;
}

macro_rules! impl_idx {
    ($($ty:ident)*) => {
        $(impl Idx for $ty {
            #[inline(always)]
            fn idx(self) -> usize {
                self as usize
            }

            #[inline(always)]
            fn uidx(idx: usize) -> Self {
                idx as $ty
            }
        })*
    };
}

impl_idx!(usize idx_t);

#[derive(Default, Clone, Copy)]
pub struct Node<K, V: Idx>
where
    K: PartialOrd + Copy + Default,
{
    key: K,
    val: V,
}

impl<K, V: Idx> Node<K, V>
where
    K: PartialOrd + Copy + Default,
{
    /// return the val field as an index
    #[inline(always)]
    fn idx(&self) -> usize {
        self.val.idx()
    }
}

/// Priority queue augmented with indices. After inserting a (key, index) pair, it can be quickly
/// accessed and/or removed using the index.
#[derive(Clone)]
pub struct IndexedPriorityQueue<K, I>
where
    I: Idx,
    K: PartialOrd + Copy + Default,
{
    locator: Box<[idx_t]>,
    heap: Vec<Node<K, I>>,
}

impl<K, I> IndexedPriorityQueue<K, I>
where
    I: Idx,
    K: PartialOrd + Copy + Default,
{
    pub fn new(maxnodes: usize) -> Self {
        Self {
            heap: Vec::new(),
            locator: vec![-1; maxnodes].into_boxed_slice(),
        }
    }

    /// clear the queue
    pub fn reset(&mut self) {
        // for (i=self.length()-1; i>=0; i -= 1)
        for i in (0..self.length()).rev() {
            self.locator[self.heap[i].idx()] = -1;
        }
        self.heap.clear();
    }

    /// get the length
    pub fn length(&self) -> usize {
        self.heap.len()
    }

    /// insert an item
    pub fn insert(&mut self, index: I, key: K) {
        debug_assert!(self.locator[index.idx()] == -1);

        let mut i = self.length();
        self.heap.push(Node { key, val: index });
        while i > 0 {
            let j = (i - 1) >> 1;
            if key < self.heap[j].key {
                self.heap[i] = self.heap[j];
                self.locator[self.heap[i].idx()] = i as idx_t;
                i = j;
            } else {
                break;
            }
        }
        // debug_assert!(i >= 0);
        self.heap[i].key = key;
        self.heap[i].val = index;
        self.locator[index.idx()] = i as idx_t;

        debug_assert!(self.check_heap());
    }

    /// delete an item
    pub fn delete(&mut self, index: I) {
        debug_assert!(self.locator[index.idx()] != -1);
        debug_assert!(self.heap[self.locator[index.idx()].idx()].val == index);

        let mut i = self.locator[index.idx()];
        self.locator[index.idx()] = -1;

        if let Some(end) = self.heap.pop() {
            if end.val != index {
                let node = end.val;
                let newkey = end.key;
                let oldkey = self.heap[i.idx()].key;

                if newkey < oldkey {
                    /* Filter-up */
                    while i > 0 {
                        let j = (i - 1) >> 1;
                        if newkey < self.heap[j.idx()].key {
                            self.heap[i.idx()] = self.heap[j.idx()];
                            self.locator[self.heap[i.idx()].idx()] = i;
                            i = j;
                        } else {
                            break;
                        }
                    }
                } else {
                    /* Filter down */
                    let nnodes = self.length();
                    loop {
                        let mut j = ((i.idx()) << 1) + 1;
                        if j >= nnodes {
                            break;
                        }
                        if self.heap[j].key < newkey {
                            if j + 1 < nnodes && self.heap[j + 1].key < self.heap[j].key {
                                j += 1;
                            }
                            self.heap[i.idx()] = self.heap[j];
                            self.locator[self.heap[i.idx()].idx()] = i;
                            i = j as idx_t;
                        } else if j + 1 < nnodes && self.heap[j + 1].key < newkey {
                            j += 1;
                            self.heap[i.idx()] = self.heap[j];
                            self.locator[self.heap[i.idx()].idx()] = i;
                            i = j as idx_t;
                        } else {
                            break;
                        }
                    }
                }

                self.heap[i.idx()].key = newkey;
                self.heap[i.idx()].val = node;
                self.locator[node.idx()] = i;
            }
        }

        debug_assert!(self.check_heap());
    }

    /// This function updates the key values associated for a particular item
    pub fn update(&mut self, node: I, newkey: K) {
        let oldkey = self.heap[self.locator[node.idx()].idx()].key;
        if (newkey >= oldkey) && (oldkey >= newkey) {
            return;
        }

        debug_assert!(self.locator[node.idx()] != -1);
        debug_assert!(self.heap[self.locator[node.idx()].idx()].val == node);

        let mut i = self.locator[node.idx()];

        if newkey < oldkey {
            /* Filter-up */
            while i > 0 {
                let j = (i.idx() - 1) >> 1;
                if newkey < self.heap[j].key {
                    self.heap[i.idx()] = self.heap[j];
                    self.locator[self.heap[i.idx()].idx()] = i;
                    i = j as idx_t;
                } else {
                    break;
                }
            }
        } else {
            /* Filter down */
            let nnodes = self.length();
            loop {
                let mut j = ((i.idx()) << 1) + 1;
                if j >= nnodes {
                    break;
                }
                if self.heap[j].key < newkey {
                    if j + 1 < nnodes && self.heap[j + 1].key < self.heap[j].key {
                        j += 1;
                    }
                    self.heap[i.idx()] = self.heap[j];
                    self.locator[self.heap[i.idx()].idx()] = i;
                    i = j as idx_t;
                } else if j + 1 < nnodes && self.heap[j + 1].key < newkey {
                    j += 1;
                    self.heap[i.idx()] = self.heap[j];
                    self.locator[self.heap[i.idx()].idx()] = i;
                    i = j as idx_t;
                } else {
                    break;
                }
            }
        }

        self.heap[i.idx()].key = newkey;
        self.heap[i.idx()].val = node;
        self.locator[node.idx()] = i;

        debug_assert!(self.check_heap());
    }

    /// This function returns the item at the top of the queue and removes it from the priority queue
    /// This should never return -1
    ///
    /// equivalent to pqueue_get_top
    #[doc(alias = "get_top")]
    pub fn pop(&mut self) -> Option<I> {
        Some(self.pop_node()?.1)
    }

    pub fn pop_node(&mut self) -> Option<(K, I)> {
        if self.length() == 0 {
            return None;
        }

        let vtx = self.heap[0];
        self.locator[vtx.idx()] = -1;

        let Node { key, val: node } = self.heap.pop().expect("heap length nonzero");
        if self.length() == 0 {
            return Some((vtx.key, vtx.val));
        }
        let mut i = 0;
        loop {
            let mut j = 2 * i + 1;
            if j >= self.length() {
                break;
            }
            if self.heap[j].key < key {
                if j + 1 < self.length() && self.heap[j + 1].key < self.heap[j].key {
                    j += 1;
                }
                self.heap[i] = self.heap[j];
                self.locator[self.heap[i].val.idx()] = i as idx_t;
                i = j;
            } else if j + 1 < self.length() && self.heap[j + 1].key < key {
                j += 1;
                self.heap[i] = self.heap[j];
                self.locator[self.heap[i].val.idx()] = i as idx_t;
                i = j;
            } else {
                break;
            }
        }

        self.heap[i].key = key;
        self.heap[i].val = node;
        self.locator[node.idx()] = i as idx_t;

        debug_assert!(self.check_heap());
        assert_ne!(vtx.idx() as isize, -1);
        Some((vtx.key, vtx.val))
    }

    /// Returns the item at the top of the queue. The item is not deleted from the queue.
    #[doc(alias = "see_top_val")]
    pub fn peek_val(&mut self) -> Option<I> {
        self.heap.first().map(|n| n.val)
    }

    /// Returns the key of the top item. The item is not deleted from the queue.
    #[doc(alias = "see_top_key")]
    pub fn peek_key(&self) -> Option<K> {
        self.heap.first().map(|n| n.key)
    }

    /// Returns the key of a specific item
    #[doc(alias = "see_key")]
    pub fn get(&self, index: I) -> K {
        self.heap[self.locator[index.idx()].idx()].key
    }

    pub fn try_get(&self, index: I) -> Option<K> {
        let loc = self.locator[index.idx()];
        if loc == -1 {
            None
        } else {
            Some(self.heap[loc.idx()].key)
        }
    }

    /// asserts the consistency of the heap
    #[allow(unreachable_code)]
    pub fn check_heap(&self) -> bool {
        // too slow to do often
        let heap = &self.heap;
        let locator = &self.locator;
        let nnodes = self.length();

        if nnodes == 0 {
            return true;
        }

        assert!(locator[heap[0].idx()] == 0);
        for i in 1..nnodes {
            assert!(locator[heap[i].idx()] == i as idx_t);
            assert!(heap[i].key >= heap[(i - 1) / 2].key);
        }

        assert_eq!(locator.iter().filter(|&&i| i != -1).count(), nnodes);

        true
    }
}

impl IPQueue {
    pub fn to_dal(self) -> dal::DirectAccessList {
        let ptr = self.locator.to_vec();
        let ind = self.heap.into_iter().map(|Node { val, .. }| val).collect();
        // eprintln!("ptr: {ptr:?}\nind: {ind:?}");
        if cfg!(debug_assertions) {
            dal::DirectAccessList::from_vecs(ptr, ind)
        } else {
            unsafe { dal::DirectAccessList::from_vecs_unchecked(ptr, ind) }
        }
    }
}

/// This function selects the partition number and the queue from which we will move vertices out.
///
/// this is originally from fm.c, but I needed to rewrite it early
///
/// Returns (from, cnum) from original signature
pub fn select_queue(
    graph: &crate::graph_t,
    pijbm: &[real_t],
    ubfactors: &[real_t],
    queues: &[RPQueue],
) -> (idx_t, idx_t) {
    // idx_t ncon, i, part;
    // real_t max, tmp;

    let ncon = graph.ncon;
    let pwgts: &[_] = unsafe { std::slice::from_raw_parts(graph.pwgts, (ncon * 2) as usize) };

    let mut from = -1;
    let mut cnum = -1;
    // from = -1;
    // cnum = -1;

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
                from = part;
                cnum = i;
            }
        }
    }

    if from != -1 {
        /* in case the desired queue is empty, select a queue from the same side */
        if (queues[(2 * cnum + from) as usize]).length() == 0 {
            let mut i = 0;
            for ii in 0..ncon {
                i = ii;
                if (queues[(2 * i + from) as usize]).length() > 0 {
                    max = pwgts[(from * ncon + i) as usize] as real_t
                        * pijbm[(from * ncon + i) as usize]
                        - ubfactors[(i) as usize];
                    cnum = i;
                    break;
                }
            }
            i += 1;
            // for (i++; i<ncon; i++) {
            while i < ncon {
                let tmp = pwgts[(from * ncon + i) as usize] as real_t
                    * pijbm[(from * ncon + i) as usize]
                    - ubfactors[i as usize];
                if tmp > max && (queues[(2 * i + from) as usize]).length() > 0 {
                    max = tmp;
                    cnum = i;
                }
                i += 1;
            }
        }

        /*
        printf("Selected1 %"PRIDX"(%"PRIDX") . %"PRIDX" [%5"PRREAL"]\n",
            from, cnum, rpqLength(queues[2*(cnum)+(from)]), max);
        */
    } else {
        /* the partitioning does not violate balancing constraints, in which case select
        a queue based on cut criteria */
        for part in 0..2 {
            for i in 0..ncon {
                if (queues[(2 * i + part) as usize]).length() > 0
                    && (from == -1 || (queues[(2 * i + part) as usize].peek_key().unwrap()) > max)
                {
                    max = queues[(2 * i + part) as usize].peek_key().unwrap();
                    from = part;
                    cnum = i;
                }
            }
        }
        /*
        printf("Selected2 %"PRIDX"(%"PRIDX") . %"PRIDX"\n",
            *from, *cnum, rpqLength(queues[2*(*cnum)+(*from)]), max);
        */
    }
    (from, cnum)
}

// pub struct BatchedPriorityQueueUninit {
//     keys: VecDeque<idx_t>,
//     locator: Box<[idx_t]>
// }

pub struct DoublePriorityQueue<K, I>
where
    I: Idx,
    K: PartialOrd + Copy + Default,
{
    min: IndexedPriorityQueue<K, I>,
    max: IndexedPriorityQueue<Reverse<K>, I>,
}

impl<K, I> DoublePriorityQueue<K, I>
where
    I: Idx,
    K: PartialOrd + Copy + Default,
{
    pub fn new(size: usize) -> Self {
        Self {
            min: IndexedPriorityQueue::new(size),
            max: IndexedPriorityQueue::new(size),
        }
    }
}

#[cfg(test)]
mod test {
    use super::IndexedPriorityQueue;
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;

    mod heap {
        use super::*;

        #[test]
        fn in_order() {
            let mut heap = IndexedPriorityQueue::new(10);
            for x in 0..10 {
                heap.insert(x, x);
            }
            let mut it = 0..10;
            while let Some(k) = heap.pop() {
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
                truth.push((Reverse(k), k)); // std::heap is max-heap
            }

            for _ in 0..80 {
                eprintln!(
                    "gk: {:>4?}  |  std: {:>4?}",
                    heap.peek_key(),
                    truth.peek().map(|x| x.1)
                );
                assert_eq!(heap.peek_key(), truth.peek().map(|x| x.1));
                heap.pop();
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
                truth.push((Reverse(k), k));
            }
            for x in 0..35 {
                heap.delete(x);
            }

            for _ in 0..80 {
                eprintln!(
                    "gk: {:>4?}  |  std: {:>4?}",
                    heap.peek_key(),
                    truth.peek().map(|x| x.1)
                );
                assert_eq!(heap.peek_key(), truth.peek().map(|x| x.1));
                heap.pop();
                truth.pop();
            }
        }

        #[test]
        fn update_entries() {
            let mut queue = IndexedPriorityQueue::new(10);
            let v = vec![0, 3, 2, 1, 5, 9, 7, 6, 4, 8];
            for i in 0..10 {
                queue.insert(i, v[i as usize]);
            }
            for i in 0..10 {
                assert_eq!(queue.get(i), v[i as usize]);
            }
            assert!(queue.check_heap());
            for i in 0..10 {
                queue.update(i, i);
                assert!(queue.check_heap());
            }
            for i in 0..10 {
                assert_eq!(queue.get(i), i);
            }
        }

        #[test]
        fn delete_indices() {
            let mut queue = IndexedPriorityQueue::new(10);
            let v = vec![0, 3, 2, 1, 5, 9, 7, 6, 4, 8];
            let step = 3;
            for i in 0..10 {
                queue.insert(i, v[i as usize]);
            }
            let mut len = queue.length();
            for i in (0..10).step_by(step) {
                queue.delete(i);
                len -= 1;
            }
            assert_eq!(len, queue.length());
            for i in 0..10 {
                if (i as usize) % step == 0 {
                    assert_eq!(queue.locator[i as usize], -1);
                } else {
                    assert_eq!(queue.get(i), v[i as usize]);
                }
            }
        }
    }
}
