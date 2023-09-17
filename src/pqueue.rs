/*
\file  gk_mkpqueue.h
\brief Templates for priority queues

\date   Started 4/09/07
\author George
\version\verbatim $Id: gk_mkpqueue.h 21742 2018-01-26 16:59:15Z karypis $ \endverbatim
*/
use std::fmt::Debug;

#[derive(Default)]
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

pub struct Mheap<K, V>
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

impl<KT, VT> Mheap<KT, VT>
where
    VT: Copy + Default + PartialEq,
    usize: TryFrom<VT>,
    <usize as TryFrom<VT>>::Error: Debug,
    KT: PartialOrd + Copy + Default,
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

    /*************************************************************************/
    /* This function resets the priority queue */
    /**************************************************************************/
    pub fn reset(&mut self) {
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        // for (i=self.nnodes-1; i>=0; i -= 1)
        for i in (0..self.nnodes).rev() {
            locator[TryInto::<usize>::try_into(heap[i].val).expect("valid index")] = -1;
        }
        self.nnodes = 0;
    }

    /*************************************************************************/
    /* This function returns the length of the queue */
    /**************************************************************************/
    pub fn length(&mut self) -> usize {
        return self.nnodes;
    }

    /*************************************************************************/
    /* This function adds an item in the priority queue */
    /**************************************************************************/
    pub fn insert(&mut self, node: VT, key: KT) {
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        debug_assert!(self.check_heap());

        debug_assert!(locator[TryInto::<usize>::try_into(node).expect("valid index")] == -1);

        let mut i = self.nnodes;
        self.nnodes += 1;
        while i > 0 {
            let j = (i - 1) >> 1;
            if key < heap[j].key {
                heap[i] = heap[j];
                locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i as isize;
                i = j;
            } else {
                break;
            }
        }
        debug_assert!(i >= 0);
        heap[i].key = key;
        heap[i].val = node;
        locator[TryInto::<usize>::try_into(node).expect("valid index")] = i as isize;

        debug_assert!(self.check_heap());
    }

    /*************************************************************************/
    /* This function deletes an item from the priority queue */
    /**************************************************************************/
    pub fn delete(&mut self, node: VT) {
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        debug_assert!(locator[TryInto::<usize>::try_into(node).expect("valid index")] != -1);
        debug_assert!(heap[locator[TryInto::<usize>::try_into(node).expect("valid index")] as usize].val == node);

        debug_assert!(self.check_heap());

        let mut i = locator[TryInto::<usize>::try_into(node).expect("valid index")];
        locator[TryInto::<usize>::try_into(node).expect("valid index")] = -1;

        self.nnodes -= 1;
        if self.nnodes > 0 && heap[self.nnodes].val != node {
            let node = heap[self.nnodes].val;
            let mut newkey = heap[self.nnodes].key;
            let oldkey = heap[i as usize].key;

            if newkey < oldkey {
                /* Filter-up */
                while i > 0 {
                    let j = (i - 1) >> 1;
                    if newkey < heap[j as usize].key {
                        heap[i as usize] = heap[j as usize];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                        i = j;
                    } else {
                        break;
                    }
                }
            } else {
                /* Filter down */
                let nnodes = self.nnodes;
                while let mut j = ((i as usize) << 1) + 1{
                    if j >= nnodes { break; }
                    if heap[j].key < newkey {
                        if j + 1 < nnodes && heap[j + 1].key < heap[j].key {
                            j += 1;
                        }
                        heap[i as usize] = heap[j];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                        i = j as isize;
                    } else if j + 1 < nnodes && heap[j + 1].key < newkey {
                        j += 1;
                        heap[i as usize] = heap[j];
                        locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                        i = j as isize;
                    } else {
                        break;
                    }
                }
            }

            heap[i as usize].key = newkey;
            heap[i as usize].val = node;
            locator[TryInto::<usize>::try_into(node).expect("valid index")] = i;
        }

        debug_assert!(self.check_heap());
    }

    /*************************************************************************/
    /* This function updates the key values associated for a particular item */
    /**************************************************************************/
    pub fn update(&mut self, node: VT, newkey: KT) {
        let locator = &mut self.locator;
        let heap = &mut self.heap;

        let oldkey = heap[locator[TryInto::<usize>::try_into(node).expect("valid index")] as usize].key;
        if !(newkey < oldkey) && !(oldkey < newkey) {
            return;
        }

        debug_assert!(locator[TryInto::<usize>::try_into(node).expect("valid index")] != -1);
        debug_assert!(heap[locator[TryInto::<usize>::try_into(node).expect("valid index")] as usize].val == node);
        debug_assert!(self.check_heap());

        let mut i = locator[TryInto::<usize>::try_into(node).expect("valid index")];

        if newkey < oldkey {
            /* Filter-up */
            while i > 0 {
                let j = (i as usize - 1) >> 1;
                if newkey < heap[j].key {
                    heap[i as usize] = heap[j];
                    locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                    i = j as isize;
                } else {
                    break;
                }
            }
        } else {
            /* Filter down */
            let nnodes = self.nnodes;
            while let mut j = ((i as usize) << 1) + 1 {
                if !(j < nnodes) { break; }
                if heap[j].key < newkey {
                    if j + 1 < nnodes && heap[j + 1].key < heap[j].key {
                        j += 1;
                    }
                    heap[i as usize] = heap[j];
                    locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                    i = j as isize;
                } else if j + 1 < nnodes && heap[j + 1].key < newkey {
                    j += 1;
                    heap[i as usize] = heap[j];
                    locator[TryInto::<usize>::try_into(heap[i as usize].val).expect("valid index")] = i;
                    i = j as isize;
                } else {
                    break;
                }
            }
        }

        heap[i as usize].key = newkey;
        heap[i as usize].val = node;
        locator[TryInto::<usize>::try_into(node).expect("valid index")] = i;

        debug_assert!(self.check_heap());

        return;
    }

    /*************************************************************************/
    /* This function returns the item at the top of the queue and removes
    it from the priority queue */
    /**************************************************************************/
    pub fn get_top(&mut self) -> Option<VT> {
        debug_assert!(self.check_heap());

        if self.nnodes == 0 {
            return None;
        }

        self.nnodes -= 1;

        let heap = &mut self.heap;
        let locator = &mut self.locator;

        let mut vtx= heap[0].val;
        locator[TryInto::<usize>::try_into(vtx).expect("valid index")] = -1;

        let mut i = self.nnodes;
        if i > 0 {
            let key = heap[i].key;
            let node = heap[i].val.into();
            i = 0;
            while let mut j = 2 * i + 1 {
                if !(j < self.nnodes) { break; }
                if heap[j].key < key {
                    if j + 1 < self.nnodes && heap[j + 1].key < heap[j].key {
                        j = j + 1;
                    }
                    heap[i] = heap[j];
                    locator[TryInto::<usize>::try_into(heap[i].val).expect("valid index")] = i as isize;
                    i = j;
                } else if j + 1 < self.nnodes && heap[j + 1].key < key {
                    j = j + 1;
                    heap[i] = heap[j];
                    locator[TryInto::<usize>::try_into(heap[i].val).expect("valid index")] = i as isize;
                    i = j;
                } else {
                    break;
                }
            }

            heap[i].key = key;
            heap[i].val = node;
            locator[TryInto::<usize>::try_into(node).expect("valid index")] = i as isize;
        }

        debug_assert!(self.check_heap());
        return Some(vtx);
    }

    /*************************************************************************/
    /* This function returns the item at the top of the queue. The item is not
    deleted from the queue. */
    /**************************************************************************/
    pub fn see_top_val(&mut self) -> Option<VT> {
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
    pub fn see_top_key(&mut self) -> Option<KT> {
        if self.nnodes == 0 {
            None
        } else {
            Some(self.heap[0].key)
        }
    }

    /*************************************************************************/
    /* This function returns the key of a specific item */
    /**************************************************************************/
    pub fn see_key(&self, node: VT) -> KT {
        return self.heap[self.locator[TryInto::<usize>::try_into(node).expect("valid index")] as usize].key;
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

        debug_assert!(locator[TryInto::<usize>::try_into(heap[0].val).expect("valid index")] == 0);
        for i in 1..nnodes {
            debug_assert!(locator[TryInto::<usize>::try_into(heap[i].val).expect("valid index")] == i as isize);
            debug_assert!(heap[i].key >= heap[(i - 1) / 2].key);
        }
        for i in 1..nnodes {
            debug_assert!(heap[i].key >= heap[0].key);
        }

        let mut j = 0;
        for i in 0..self.heap.len() {
            if locator[i] != -1 {
                j += 1;
            }
        }
        debug_assert!(j == nnodes, "{} {}", j, nnodes);

        return true;
    }
}

#[cfg(test)]
mod test {
    use super::Mheap;

    macro_rules! basic_test {
        ($name:ident: $amt:expr, $items:expr) => {

        };
    }

    #[test]
    fn basic() {
        let mut heap = Mheap::new(10);
        for x in 0..10 {
            heap.insert(x, x);
        }
        let mut it = 0..10;
        while let Some(k) = heap.get_top() {
            assert_eq!(Some(k), it.next());
        }
        assert_eq!(it.next(), None);
    }
}
