//! Replacement routines for gklib that ought to never use the original
//!
//! IMPORTANT: Any functions added here probably require corresponding functions in
//! `normalization/functions.rs`

use std::cell::Cell;

use crate::{idx_t, ikv_t};


#[cfg(feature = "normalized")]
mod replace_decl {
    use crate::ikv_t;

    // Implementations here are in normalization/functions.rs

    // Why do we need to link against this even in the Rust version? Because it appears that the
    // version of rustc that metis-normalized builds against (in Nix) and the current nightly I'm
    // using have different implementations of `sort_unstable_by_key`.
    #[metis_decl]
    extern "C" {
        pub fn ikvsorti(n: usize, ikv: *mut ikv_t);
        pub fn ikvsortd(n: usize, ikv: *mut ikv_t);
    }

    unsafe extern "C" {
        pub safe fn gk_randinit(seed: u64);
        pub safe fn gk_randint64() -> u64;
        pub safe fn gk_randint32() -> u32;
    }
}

#[cfg(not(feature = "normalized"))]
mod replace_decl {
    use crate::ikv_t;

    #[no_mangle]
    pub unsafe extern "C" fn libmetis__ikvsorti(n: usize, ikv: *mut ikv_t) {
        crate::mkslice_mut!(ikv, n);
        ikv.sort_unstable_by_key(|kv| kv.key);
        ikvsorti(ikv)
    }

    #[no_mangle]
    #[cfg(not(feature = "normalized"))]
    pub unsafe extern "C" fn libmetis__ikvsortd(n: usize, ikv: *mut ikv_t) {
        crate::mkslice_mut!(ikv, n);
        ikv.sort_unstable_by_key(|kv| std::cmp::Reverse(kv.key));
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
    #[inline(never)]
    pub extern "C" fn gk_randint32() -> u32 {
        // must be positive, or you get horrible segfaults!
        (gk_randint64() & 0x7FFFFFFF) as u32
    }
}


/// Always use this over stdlib sorting functions!
pub fn ikvsorti(ikv: &mut [ikv_t]) {
    unsafe { replace_decl::ikvsorti(ikv.len(), ikv.as_mut_ptr()) }
}

/// Always use this over stdlib sorting functions!
pub fn ikvsortd(ikv: &mut [ikv_t]) {
    unsafe { replace_decl::ikvsortd(ikv.len(), ikv.as_mut_ptr()) }
}

pub fn gk_randinit(seed: u64) {
    replace_decl::gk_randinit(seed) 
}

pub fn gk_randint64() -> u64 {
    replace_decl::gk_randint64() 
}

pub fn gk_randint32() -> u32 {
    replace_decl::gk_randint32() 
}
