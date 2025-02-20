//! Replacement routines for gklib that ought to never use the original

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


#[cfg(any())]
mod for_future_miri {
    use std::{alloc::{self, Layout}, ffi::c_void};
    use crate::{idx_t, real_t};

    #[metis_func]
    pub extern "C" fn imalloc(nelm: usize, _label: *const i8) -> *mut c_void {
        assert!(nelm > 0);
        alloc::alloc_zeroed(Layout::array::<idx_t>(nelm).expect("valid size")).cast()
    }

    #[metis_func]
    pub extern "C" fn ismalloc(nelm: usize, val: idx_t, _label: *const i8) -> *mut c_void {
        assert!(nelm > 0);
        let ret = alloc::alloc(Layout::array::<idx_t>(nelm).expect("valid size")).cast::<idx_t>();
        for i in 0..nelm {
            (ret.add(i)).write(val);
        }
        ret.cast()
    }

    #[metis_func]
    pub extern "C" fn rmalloc(nelm: usize, _label: *const i8) -> *mut c_void {
        assert!(nelm > 0);
        alloc::alloc_zeroed(Layout::array::<real_t>(nelm).expect("valid size")).cast()
    }

    #[metis_func]
    pub extern "C" fn rsmalloc(nelm: usize, val: real_t, _label: *const i8) -> *mut c_void {
        assert!(nelm > 0);
        let ret = alloc::alloc(Layout::array::<real_t>(nelm).expect("valid size")).cast::<real_t>();
        for i in 0..nelm {
            (ret.add(i)).write(val);
        }
        ret.cast()
    }


    #[no_mangle]
    pub extern "C" fn gk_malloc(sz: usize, _label: *const i8) -> *mut c_void {
        assert!(sz > 0);
        unsafe { alloc::alloc_zeroed(Layout::from_size_align(sz, align_of::<usize>()).expect("valid size")) }.cast()
    }
}
