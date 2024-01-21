//! Direct access list wrapper. This isn't used in the main routines for now, but I want to in the
//! future
#![allow(unused)]

use std::mem::MaybeUninit;

use crate::*;

#[derive(Default, Clone)]
pub struct DirectAccessList {
    ptr: Vec<idx_t>,
    ind: Vec<idx_t>,
}

impl DirectAccessList {
    pub fn from_vecs(ptr: Vec<idx_t>, ind: Vec<idx_t>) -> Self {
        let mut count = 0;
        let mut has_entry = vec![false; ind.len()];
        for (i, &k) in ptr.iter().enumerate() {
            assert!(k >= -1);
            assert!(k < ind.len() as idx_t);
            if k != -1 {
                count += 1;
                assert!(!has_entry[k as usize]);
                has_entry[k as usize] = true;
                assert_eq!(i as idx_t, ind[k as usize]);
            }
        }
        assert_eq!(count, ind.len());
        Self { ptr, ind }
    }

    pub fn clear(&mut self) {
        // Safety: we clear ind
        unsafe {
            self.clear_ptr();
        }
        self.ind.clear();
    }

    /// clears just ptr, leaving ind intact
    ///
    /// ## Safety
    /// ind must be manually cleared after
    unsafe fn clear_ptr(&mut self) {
        let max = self.ptr.len();
        let len = self.len();
        if len >= max / 8 || max < 32 {
            // do a memset if enough of the ptr array is filled enough
            // the criteria was a random guess
            self.ptr.fill(-1);
        } else {
            for &idx in self.ind.iter() {
                self.ptr[idx as usize] = -1;
            }
        }
    }

    /// sort the ind list
    //
    // TODO: optimize me
    pub fn sort(&mut self) {
        // don't do anything if we're already sorted
        if self.ind.windows(2).all(|w| w[0] <= w[1]) {
            return;
        }

        unsafe {
            self.clear_ptr();
        }
        let mut ind = std::mem::take(&mut self.ind);
        ind.sort_unstable();
        for (i, &val) in ind.iter().enumerate() {
            self.ptr[val as usize] = i as idx_t;
        }
    }

    /// return a (ptr, ind) tuple
    pub fn as_slices(&self) -> (&[idx_t], &[idx_t]) {
        (&self.ptr, &self.ind)
    }

    /// return a mutable (ptr, ind) tuple
    pub unsafe fn as_slices_mut(&mut self) -> (&mut [idx_t], &mut [idx_t]) {
        (&mut self.ptr, &mut self.ind)
    }

    pub unsafe fn from_vecs_unchecked(ptr: Vec<idx_t>, ind: Vec<idx_t>) -> Self {
        debug_assert_eq!(ind.iter().find(|&&x| x < 0), None);
        debug_assert!(ptr.len() >= ind.len());
        Self { ptr, ind }
    }

    pub fn new(cap: usize) -> Self {
        Self {
            ptr: vec![-1; cap],
            ind: Vec::with_capacity(cap),
        }
    }

    /// returns true if key was not already in set
    pub fn insert(&mut self, key: idx_t) -> bool {
        if self.contains(key) {
            return false;
        }
        let key = key as usize;

        self.ptr[key] = self.ind.len() as idx_t;
        self.ind.push(key as idx_t);
        true
    }

    /// removes key from the list, returns true if it existed
    pub fn remove(&mut self, key: idx_t) -> bool {
        if !self.contains(key) {
            return false;
        }
        let key = key as usize;

        let idx = self.ptr[key];
        // set the value of the most recent value to the slot it will be moved to. This may be the
        // same as key, so we have to do it before overwriting ptr[key]
        self.ptr[*self.ind.last().expect("ind is nonempty") as usize] = self.ptr[key];
        self.ptr[key] = -1;
        self.ind.swap_remove(idx as usize);
        true
    }

    /// returns true if list contains key
    pub fn contains(&self, key: idx_t) -> bool {
        debug_assert!(key >= 0);
        let key = key as usize;
        debug_assert!(key < self.ptr.len());
        self.ptr[key] != -1
    }

    /// iterator of the keys in the list in an unspecified order. O(n).
    pub fn iter(&self) -> impl Iterator<Item = &idx_t> {
        self.ind.iter()
    }

    /// in-order iterator of the keys in the list. O(capacity).
    pub fn iter_ordered(&self) -> impl Iterator<Item = &idx_t> {
        self.ptr
            .iter()
            .filter(|&&i| i >= 0)
            .map(|&i| &self.ind[i as usize])
    }

    pub fn len(&self) -> usize {
        self.ind.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ind.is_empty()
    }

    pub fn cap(&self) -> usize {
        self.ptr.len()
    }

    fn key_ind(&self, key: idx_t) -> Option<usize> {
        if !self.contains(key) {
            return None;
        }
        let key = key as usize;
        return Some(self.ptr[key] as usize);
    }
}

impl From<pqueue::IPQueue> for DirectAccessList {
    fn from(value: pqueue::IPQueue) -> Self {
        value.to_dal()
    }
}

impl Extend<idx_t> for DirectAccessList {
    fn extend<T: IntoIterator<Item = idx_t>>(&mut self, iter: T) {
        for key in iter {
            self.insert(key);
        }
    }
}

impl std::fmt::Debug for DirectAccessList {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_set().entries(self.iter()).finish()
    }
}

pub struct DirectAccessMap<T> {
    list: DirectAccessList,
    map: Vec<T>,
}

impl<T> DirectAccessMap<T> {
    pub fn new(cap: usize) -> Self {
        Self {
            list: DirectAccessList::new(cap),
            map: Vec::with_capacity(cap),
        }
    }

    pub fn inner(&self) -> &DirectAccessList {
        &self.list
    }

    pub fn insert(&mut self, key: idx_t, item: T) -> Option<T> {
        if let Some(key_ind) = self.list.key_ind(key) {
            Some(std::mem::replace(&mut self.map[key_ind], item))
        } else {
            self.list.insert(key);
            self.map.push(item);
            None
        }
    }

    pub fn remove(&mut self, key: idx_t) -> Option<T> {
        if let Some(key_ind) = self.list.key_ind(key) {
            self.list.remove(key);
            Some(self.map.swap_remove(key_ind))
        } else {
            None
        }
    }

    pub fn sort(&mut self) {
        let oldind = self.list.ind.clone();
        let mut newarr: Vec<MaybeUninit<T>> = Vec::with_capacity(oldind.len());
        self.list.sort();
        let ptr = &self.list.ptr;
        unsafe {
            newarr.set_len(oldind.len());
            for i in 0..oldind.len() {
                let dst = newarr[ptr[oldind[i] as usize] as usize].as_mut_ptr();
                let src = &self.map[i] as *const T;
                // eprintln!("{i} moved to {}", ptr[oldind[i] as usize]);
                std::ptr::copy_nonoverlapping(src, dst, 1);
            }
            assert!(
                core::alloc::Layout::new::<T>() == core::alloc::Layout::new::<MaybeUninit<T>>()
            );
            self.map.set_len(0);
            let p = newarr.as_mut_ptr();
            let len = newarr.len();
            let cap = newarr.capacity();
            std::mem::forget(newarr);
            self.map = Vec::from_raw_parts(p as *mut T, len, cap);
        }
    }

    pub fn sort_values(&mut self)
    where
        T: Ord,
    {
        let mut mapping: Vec<_> = (0..self.len()).collect();
        mapping.sort_unstable_by_key(|&i| &self.map[i]);
        let oldind = self.list.ind.clone();
        for i in 0..self.len() {
            self.list.ind[i] = oldind[i];
            self.list.ptr[oldind[i] as usize] = i as idx_t;
        }
        todo!("this doesn't work like at all - what was I thinking")
    }

    pub fn values_slice_mut(&mut self) -> &mut [T] {
        &mut self.map
    }
    pub fn values_slice(&self) -> &[T] {
        &self.map
    }

    /// returns (ptr, ind) of inner list
    pub fn list_slices(&self) -> (&[idx_t], &[idx_t]) {
        self.list.as_slices()
    }

    pub fn keys_slice(&self) -> &[idx_t] {
        &self.list.ind
    }

    pub fn contains(&self, key: idx_t) -> bool {
        self.list.contains(key)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&idx_t, &T)> {
        self.list.iter().zip(self.map.iter())
    }

    pub fn keys(&self) -> impl Iterator<Item = &idx_t> {
        self.list.iter()
    }

    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.map.iter()
    }

    pub fn iter_ordered(&self) -> impl Iterator<Item = (&idx_t, &T)> {
        self.list
            .ptr
            .iter()
            .enumerate()
            .filter(|(_, &k)| k >= 0)
            .map(|(i, k)| {
                (
                    &self.list.ind[self.list.ptr[i] as usize],
                    &self.map[self.list.ptr[i] as usize],
                )
            })
    }

    pub fn keys_ordered(&self) -> impl Iterator<Item = &idx_t> {
        self.list.iter_ordered()
    }

    pub fn values_ordered(&self) -> impl Iterator<Item = &T> {
        self.list
            .ptr
            .iter()
            .enumerate()
            .filter(|(_, &k)| k >= 0)
            .map(|(i, _)| &self.map[self.list.ptr[i] as usize])
    }

    pub fn len(&self) -> usize {
        self.list.len()
    }
    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    pub fn cap(&self) -> usize {
        self.list.cap()
    }

    pub fn get(&self, key: idx_t) -> Option<&T> {
        if let Some(key_ind) = self.list.key_ind(key) {
            Some(&self.map[key_ind])
        } else {
            None
        }
    }

    pub fn get_mut(&mut self, key: idx_t) -> Option<&mut T> {
        if let Some(key_ind) = self.list.key_ind(key) {
            Some(&mut self.map[key_ind])
        } else {
            None
        }
    }
}

impl<T> std::ops::Index<usize> for DirectAccessMap<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        #[inline(never)]
        fn index_fail(index: usize) -> ! {
            panic!("No entry found for key {index}");
        }
        match index.try_into().map_or(None, |idx| self.get(idx)) {
            Some(r) => r,
            None => index_fail(index),
        }
    }
}

impl<T> std::ops::Index<idx_t> for DirectAccessMap<T> {
    type Output = T;

    fn index(&self, index: idx_t) -> &Self::Output {
        #[inline(never)]
        fn index_fail(index: idx_t) -> ! {
            panic!("No entry found for key {index}");
        }
        match self.get(index) {
            Some(r) => r,
            None => index_fail(index),
        }
    }
}

impl<T> std::ops::IndexMut<usize> for DirectAccessMap<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        #[inline(never)]
        fn index_fail(index: usize) -> ! {
            panic!("No entry found for key {index}");
        }
        match index.try_into().map_or(None, |idx| self.get_mut(idx)) {
            Some(r) => r,
            None => index_fail(index),
        }
    }
}

impl<T> std::ops::IndexMut<idx_t> for DirectAccessMap<T> {
    fn index_mut(&mut self, index: idx_t) -> &mut Self::Output {
        #[inline(never)]
        fn index_fail(index: idx_t) -> ! {
            panic!("No entry found for key {index}");
        }
        match self.get_mut(index) {
            Some(r) => r,
            None => index_fail(index),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    mod dal {
        use std::collections::BTreeSet;

        use super::*;

        fn dal_equiv(dal: &DirectAccessList, bt: &BTreeSet<idx_t>) {
            assert_eq!(dal.len(), bt.len());
            assert_eq!(dal.iter().count(), dal.len());
            let normal_it_set: BTreeSet<idx_t> = dal.iter().copied().collect();
            assert_eq!(&normal_it_set, bt);
            for (dal_item, bt_item) in dal.iter_ordered().zip(bt.iter()) {
                assert_eq!(dal_item, bt_item);
            }
        }

        macro_rules! dal_test {
            ($name:ident, $cap:expr => [$($op:ident($($arg:tt),* $(,)?)),* $(,)?]) => {
                #[test]
                fn $name() {
                    #[allow(unused_mut)]
                    let mut dal = DirectAccessList::new($cap);
                    let mut bt = BTreeSet::<idx_t>::new();
                    $(
                    println!("executing dal.{}({})", stringify!($op), concat!($(stringify!($arg), ", "),*));
                    let dal_res = dal.$op($($arg,)*);
                    let bt_res = dal_test!(@coerce_bt bt, $op, $($arg),*);
                    assert_eq!(dal_res, bt_res);
                    dal_equiv(&dal, &bt);
                )*
                }
            };
            (@coerce_bt $bt:ident, contains, $arg:expr) => {
                $bt.contains(&$arg)
            };
            (@coerce_bt $bt:ident, remove, $arg:expr) => {
                $bt.remove(&$arg)
            };
            (@coerce_bt $bt:ident, $fallback:ident, $arg:expr) => {
                $bt.$fallback($arg)
            };
        }

        dal_test! {
            add_one, 5 => [insert(1)]
        }
        dal_test! {
            add_multiple, 5 => [insert(1), insert(3)]
        }
        dal_test! {
            add_til_full, 5 => [insert(1), insert(2), insert(3), insert(0), insert(4)]
        }
        dal_test! {
            add_til_dup, 5 => [insert(1), insert(2), insert(3), insert(2), insert(4)]
        }
        dal_test! {
            add_contains, 5 => [contains(2), insert(2), contains(2)]
        }
        dal_test! {
            add_remove, 5 => [insert(2), contains(2), remove(2), contains(2)]
        }
        dal_test! {
            add_remove_repeat, 5 => [
                insert(2),
                contains(2),
                remove(2),
                contains(2),
                remove(2),
                remove(2),
                insert(2),
                remove(2),
                insert(3),
                insert(2),
                remove(2),
            ]
        }
        dal_test! {
            add_many_and_remove, 5 => [remove(2), insert(2), insert(1), insert(3), remove(2), contains(2)]
        }
        dal_test! {
            extend, 5 => [extend([1, 3, 4, 2, 0]), insert(1), insert(3), remove(2), contains(2)]
        }
        dal_test! {
            extend_dups, 5 => [extend([1, 3, 4, 3, 0]), extend([0, 0])]
        }
        dal_test! {
            larger, 15 => [
                extend([1, 3, 4, 3, 0, 3, 12, 10, 5, 9, 0, 3]),
                extend([0, 0, 11, 14]),
                remove(12),
                remove(0),
                remove(3),
                extend([0, 0, 11, 14]),
                remove(14),
                contains(0),
                remove(0),
            ]
        }

        #[test]
        fn sort() {
            let v = vec![0, 1, 4, 2, 3, 9, 5, 30, 24, 29, 5, 4, 2];
            let mut dal = DirectAccessList::new(40);
            dal.extend(v.iter().copied());
            dal.sort();
            for w in dal.ind.windows(2) {
                assert!(w[0] <= w[1]);
            }
            for (i, &k) in dal.ind.iter().enumerate() {
                assert_eq!(dal.key_ind(k), Some(i));
            }
        }
    }
    mod map {
        use std::collections::BTreeMap;

        use super::*;

        fn dal_equiv(dal: &DirectAccessMap<i64>, bt: &BTreeMap<idx_t, i64>) {
            assert_eq!(dal.len(), bt.len());
            assert_eq!(dal.iter().count(), dal.len());
            let normal_it_set: BTreeMap<idx_t, i64> = dal.iter().map(|(&k, &v)| (k, v)).collect();
            assert_eq!(&normal_it_set, bt);
            for (dal_item, bt_item) in dal.iter_ordered().zip(bt.iter()) {
                assert_eq!(dal_item, bt_item);
            }
        }

        macro_rules! dal_test {
            ($name:ident, $cap:expr => [$($op:ident($($arg:tt),* $(,)?)),* $(,)?]) => {
                #[test]
                fn $name() {
                    #[allow(unused_mut)]
                    let mut dal = DirectAccessMap::<i64>::new($cap);
                    let mut bt = BTreeMap::<idx_t, i64>::new();
                    $(
                    println!("executing dal.{}({})", stringify!($op), concat!($(stringify!($arg), ", "),*));
                    let dal_res = dal.$op($($arg,)*);
                    let bt_res = dal_test!(@coerce_bt bt, $op, $($arg),*);
                    assert_eq!(dal_test!(@coerce_res $op, dal_res), dal_test!(@coerce_res $op, bt_res));
                    dal_equiv(&dal, &bt);
                )*
                }
            };
            (@coerce_res values, $var:ident) => {
                {
                    let mut v = $var.collect::<Vec<_>>();
                    v.sort();
                    v
                }
            };
            (@coerce_res iter, $var:ident) => {
                $var.collect::<Vec<_>>()
            };
            (@coerce_res keys, $var:ident) => {
                {
                    let mut v = $var.collect::<Vec<_>>();
                    v.sort();
                    v
                }
            };
            (@coerce_res $fallback:ident, $var:ident) => {
                $var
            };
            (@coerce_bt $bt:ident, contains, $arg:expr) => {
                $bt.contains_key(&$arg)
            };
            (@coerce_bt $bt:ident, remove, $arg:expr) => {
                $bt.remove(&$arg)
            };
            (@coerce_bt $bt:ident, $fallback:ident, $($arg:expr),*) => {
                $bt.$fallback($($arg),*)
            };
        }

        dal_test! {
            add_one, 5 => [insert(1, 3)]
        }
        dal_test! {
            add_multiple, 5 => [insert(1, 3), insert(3, 9)]
        }
        dal_test! {
            add_til_full, 5 => [insert(1, 1), insert(2, 0), insert(3, 3), insert(0, 2), insert(4, 7)]
        }
        dal_test! {
            add_til_dup, 5 => [insert(1, 3), insert(2, 8), insert(3, 2), insert(2, 5), insert(4, 8)]
        }
        dal_test! {
            add_contains, 5 => [contains(2), insert(2, 0), contains(2)]
        }
        dal_test! {
            add_remove, 5 => [insert(2, 3), contains(2), remove(2), contains(2)]
        }
        dal_test! {
            add_remove_repeat, 5 => [
                insert(2, 0),
                contains(2),
                remove(2),
                contains(2),
                remove(2),
                remove(2),
                insert(2, 1),
                remove(2),
                insert(3, 2),
                insert(2, 3),
                remove(2),
            ]
        }
        dal_test! {
            add_many_and_remove, 5 => [remove(2), insert(2, 0), insert(1, 1), insert(3, 2), remove(2), contains(2)]
        }
        dal_test! {
            larger, 15 => [
                keys(),
                values(),
                insert(1, 0),
                insert(4, 1),
                insert(14, 2),
                insert(13, 3),
                values(),
                insert(13, 2),
                insert(14, 4),
                insert(11, 5),
                values(),
                insert(13, 7),
                insert(5, 6),
                insert(5, 8),
                remove(12),
                keys(),
                remove(0),
                values(),
                remove(3),
                remove(14),
                contains(0),
                remove(0),
                keys(),
            ]
        }

        #[test]
        fn sort() {
            let vals = || (0..20).map(|i| Box::new(i));
            let keys = || (25..35).chain(10..20);
            let mut map = DirectAccessMap::new(40);
            let mut truth = BTreeMap::new();
            truth.extend(keys().zip(vals()));
            for (k, v) in keys().zip(vals()) {
                map.insert(k, v);
            }
            map.sort();
            let map_vals: Vec<_> = map.values().collect();
            let truth_vals: Vec<_> = truth.values().collect();
            assert_eq!(map_vals, truth_vals);
        }

        #[ignore = "does not work"]
        #[test]
        fn sort_values() {
            let keys = || ([1, 4, 2, 3, 5, 9, 7, 6, 8].into_iter()).map(|i| Box::new(i));
            let idx = || (25..35).chain(10..15);
            let mut map = DirectAccessMap::new(40);
            let mut truth = BTreeMap::new();
            truth.extend(keys().zip(idx()));
            for (k, i) in keys().zip(idx()) {
                map.insert(i, k);
            }
            map.sort_values();
            let map_vals: Vec<_> = map.values().collect();
            let truth_vals: Vec<_> = truth.keys().collect();
            assert_eq!(map_vals, truth_vals);
        }
    }
}
