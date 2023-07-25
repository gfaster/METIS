/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * bucketsort.c
 *
 * This file contains code that implement a variety of counting sorting
 * algorithms
 *
 * Started 7/25/97
 * George
 *
 */

use crate::*;

/// This function uses simple counting sort to return a permutation array
/// corresponding to the sorted order. The keys are arsumed to start from
/// 0 and they are positive.  This sorting is used during matching.
#[metis_func]
pub extern "C" fn BucketSortKeysInc(
    _ctrl: *const ctrl_t,
    n: idx_t,
    max: idx_t,
    keys: *const idx_t,
    tperm: *const idx_t,
    perm: *mut idx_t,
) -> () {
    let keys = unsafe { std::slice::from_raw_parts(keys, n as usize) };
    let tperm = unsafe { std::slice::from_raw_parts(tperm, n as usize) };
    let perm = unsafe { std::slice::from_raw_parts_mut(perm, n as usize) };

    let mut counts = vec![0; max as usize + 2];

    for i in keys {
        counts[*i as usize] += 1;
    }

    util::make_csr(max as usize + 1, &mut counts);

    for i in tperm {
        let new_cnt = &mut counts[keys[*i as usize] as usize];
        perm[*new_cnt as usize] = *i;
        *new_cnt += 1;
    }
}
