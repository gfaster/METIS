#![crate_type = "cdylib"]
#![allow(non_snake_case, non_camel_case_types)]

use std::cell::Cell;

type idx_t = i32;

#[repr(C)]
pub struct ikv_t {
    pub key: idx_t,
    pub val: idx_t,
}

// I'm using a manual implementation of Wyrand here to make sure random number generation within is
// identical to the normalized C version. For testing purposes, we want a thread_local

// for now just use 0 as default state, since that lines up with C version
thread_local! { static WYSTATE: Cell<u64> = const { Cell::new(0) }; }


#[no_mangle]
pub extern "C" fn gk_randinit(seed: u64) {
    WYSTATE.set(seed);
}

#[no_mangle]
pub extern "C" fn gk_randint64() -> u64 {
    // using fastrand's constants
    const WY_CONST_0: u64 = 0x2d358dccaa6c78a5;
    const WY_CONST_1: u64 = 0x8bb84b93962eacc9;

    let s =  WYSTATE.get().wrapping_add(WY_CONST_0);
    WYSTATE.set(s);

    // this can't overflow
    let t = (s as u128) * ((s ^ WY_CONST_1) as u128);

    t as u64 ^ (t >> 64) as u64
}

#[no_mangle]
pub extern "C" fn gk_randint32() -> u32 {
    // must be positive, or you get horrible segfaults!
    (gk_randint64() & 0x7FFFFFFF) as u32
}


// need sorting functions since Rust's functions re-order identical keys differently
//
// Why do we need to link against this even in the Rust version? Because it appears that the
// version of rustc that metis-normalized builds against (in Nix) and the current nightly I'm using
// have different implementations of `sort_unstable_by_key`.

#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsorti(n: usize, ikv: *mut ikv_t) {
    let ikv = std::slice::from_raw_parts_mut(ikv, n);
    ikv.sort_unstable_by_key(|kv| kv.key);
}

#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsortd(n: usize, ikv: *mut ikv_t) {
    let ikv = std::slice::from_raw_parts_mut(ikv, n);
    ikv.sort_unstable_by_key(|kv| std::cmp::Reverse(kv.key));
}
