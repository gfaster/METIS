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

/// # Description
/// This function uses simple counting sort to return a permutation array
/// corresponding to the sorted order. The keys are arsumed to start from
/// 0 and they are positive.  This sorting is used during matching.
///
/// # Notes
/// I'm really not sure about the semantics of this function. It's not super complicated, I just
/// don't get it.
#[metis_func]
pub extern "C" fn BucketSortKeysInc(
    _ctrl: *const ctrl_t,
    n: idx_t,
    max: idx_t,
    keys: *const idx_t,
    tperm: *const idx_t,
    perm: *mut idx_t,
) {
    let keys = unsafe { std::slice::from_raw_parts(keys, n as usize) };
    let tperm = unsafe { std::slice::from_raw_parts(tperm, n as usize) };
    let perm = unsafe { std::slice::from_raw_parts_mut(perm, n as usize) };

    let mut counts = vec![0; max as usize + 2];

    for i in keys {
        counts[*i as usize] += 1;
    }

    util::make_csr(max as usize + 1, &mut counts);

    for &i in tperm {
        let key = keys[i as usize] as usize;
        let cnt = &mut counts[key];
        perm[*cnt as usize] = i;
        *cnt += 1;
    }
}

mod test {
    #![allow(unused_imports)]
    use crate::*;

    #[ab_test_eq(super::BucketSortKeysInc)]
    fn bucketsort_keys_inc() -> Vec<idx_t> {
        let n = 10;
        let mut perm = vec![-1; 10];
        let tperm = vec![0, 9, 3, 2, 1, 8, 4, 5, 6, 7];
        let keys = vec![0, 5, 3, 2, 0, 8, 5, 2, 8, 7];
        let mut ctrl = util::Ctrl::new_kmetis_basic();
        ctrl.init_dummy_graph(20);
        let max = 9;

        assert_eq!(perm.len(), n);
        assert_eq!(tperm.len(), n);
        assert_eq!(keys.len(), n);

        unsafe {
            BucketSortKeysInc(
                ctrl.inner,
                n as i32,
                max,
                keys.as_ptr(),
                tperm.as_ptr(),
                perm.as_mut_ptr(),
            )
        };
        perm
    }
}
