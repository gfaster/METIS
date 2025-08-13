//! Replacement routines for gklib that ought to never use the original
//!
//! IMPORTANT: Any functions added here probably require corresponding functions in
//! `normalization/functions.rs`

use crate::{idx_t, ikv_t};

#[no_mangle]
pub extern "C" fn gk_randinit(seed: u64) {
    fastrand::seed(seed);
}

#[no_mangle]
pub extern "C" fn gk_randint64() -> u64 {
    // must be positive, or you get horrible segfaults!
    fastrand::u64(..) & (u64::MAX >> 1)
}

#[no_mangle]
pub extern "C" fn gk_randint32() -> u32 {
    // must be positive, or you get horrible segfaults!
    fastrand::u32(..) & (u32::MAX >> 1)
}


// this is not declared metis_func since I want to use the nicer version
#[doc(hidden)]
#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsorti(n: usize, ikv: *mut ikv_t) {
    crate::mkslice_mut!(ikv, n);
    ikvsorti(ikv)
}

pub fn ikvsorti(ikv: &mut [ikv_t]) {
    ikv.sort_unstable_by_key(|kv| kv.key);
}

#[doc(hidden)]
#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsortd(n: usize, ikv: *mut ikv_t) {
    crate::mkslice_mut!(ikv, n);
    ikvsortd(ikv);
}

pub fn ikvsortd(ikv: &mut [ikv_t]) {
    ikv.sort_unstable_by_key(|kv| std::cmp::Reverse(kv.key));
}
