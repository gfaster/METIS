//! GKLib wrappers for using in Rust (will be removed eventually)
#![allow(dead_code)]

use crate::*;

use std::alloc::Layout;
use std::ffi::{c_char, c_void, CStr};
use std::ptr::NonNull;

type XMallocFn = unsafe fn(usize, *const c_char) -> *mut c_void;
type XSMallocFn<T> = unsafe fn(usize, T, *const c_char) -> *mut c_void;

unsafe trait GKAllocs {
    const MALLOC: XMallocFn;
    const SMALLOC: XSMallocFn<Self>;
}

unsafe impl GKAllocs for idx_t {
    const MALLOC: XMallocFn = bindings::imalloc;
    const SMALLOC: XSMallocFn<Self> = bindings::ismalloc;
}

unsafe impl GKAllocs for real_t {
    const MALLOC: XMallocFn = bindings::rmalloc;
    const SMALLOC: XSMallocFn<Self> = bindings::rsmalloc;
}

unsafe fn wrap_xmalloc<T: GKAllocs>(nmemb: usize, msg: &'static CStr) -> NonNull<[T]> {
    let layout = Layout::array::<T>(nmemb).expect("bad allocation: layout overflow");
    let p = unsafe { T::MALLOC(nmemb, msg.as_ptr()) };
    let Some(p) = NonNull::new(p) else {
        std::alloc::handle_alloc_error(layout)
    };
    NonNull::slice_from_raw_parts(p.cast::<T>(), nmemb)
}

unsafe fn wrap_xsmalloc<T: GKAllocs>(nmemb: usize, val: T, msg: &'static CStr) -> NonNull<[T]> {
    let layout = Layout::array::<T>(nmemb).expect("bad allocation: layout overflow");
    let p = unsafe { T::SMALLOC(nmemb, val, msg.as_ptr()) };
    let Some(p) = NonNull::new(p) else {
        std::alloc::handle_alloc_error(layout)
    };
    NonNull::slice_from_raw_parts(p.cast::<T>(), nmemb)
}

pub unsafe fn imalloc(nmemb: usize, msg: &'static CStr) -> NonNull<[idx_t]> {
    wrap_xmalloc(nmemb, msg)
}

pub unsafe fn ismalloc(nmemb: usize, val: idx_t, msg: &'static CStr) -> NonNull<[idx_t]> {
    wrap_xsmalloc(nmemb, val, msg)
}

pub unsafe fn rmalloc(nmemb: usize, msg: &'static CStr) -> NonNull<[real_t]> {
    wrap_xmalloc(nmemb, msg)
}

pub unsafe fn rsmalloc(nmemb: usize, val: real_t, msg: &'static CStr) -> NonNull<[real_t]> {
    wrap_xsmalloc(nmemb, val, msg)
}
