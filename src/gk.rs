//! GKLib wrappers for using in Rust (will be removed eventually)
#![allow(dead_code)]

use crate::*;

use std::alloc::Layout;
use std::ffi::{c_char, c_void, CStr};
use std::marker::PhantomData;
use std::ptr::NonNull;

type XMallocFn = unsafe fn(usize, *const c_char) -> *mut c_void;
type XSMallocFn<T> = unsafe fn(usize, T, *const c_char) -> *mut c_void;
type XReallocFn = unsafe fn(*mut c_void, usize, *const c_char) -> *mut c_void;

fn realloc_unused(_: *mut c_void, _: usize, _: *const c_char) -> *mut c_void {
    unimplemented!("Xrealloc")
}

unsafe trait GKAllocs {
    const MALLOC: XMallocFn;
    const SMALLOC: XSMallocFn<Self>;
    const REALLOC: XReallocFn;
}

unsafe impl GKAllocs for idx_t {
    const MALLOC: XMallocFn = bindings::imalloc;
    const SMALLOC: XSMallocFn<Self> = bindings::ismalloc;
    const REALLOC: XReallocFn = bindings::irealloc;
}

unsafe impl GKAllocs for real_t {
    const MALLOC: XMallocFn = bindings::rmalloc;
    const SMALLOC: XSMallocFn<Self> = bindings::rsmalloc;
    const REALLOC: XReallocFn = realloc_unused;
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
    assert!(nmemb > 0);
    let layout = Layout::array::<T>(nmemb).expect("bad allocation: layout overflow");
    let p = unsafe { T::SMALLOC(nmemb, val, msg.as_ptr()) };
    let Some(p) = NonNull::new(p) else {
        std::alloc::handle_alloc_error(layout)
    };
    NonNull::slice_from_raw_parts(p.cast::<T>(), nmemb)
}

unsafe fn wrap_xrealloc<T: GKAllocs>(old: NonNull<[T]>, nmemb: usize, msg: &'static CStr) -> NonNull<[T]> {
    assert!(nmemb > 0);
    let layout = Layout::array::<T>(nmemb).expect("bad allocation: layout overflow");
    let p = unsafe { T::REALLOC(old.as_buf_ptr().cast(), nmemb, msg.as_ptr()) };
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

pub unsafe fn irealloc(old: NonNull<[idx_t]>, nmemb: usize, msg: &'static CStr) -> NonNull<[idx_t]> {
    wrap_xrealloc(old, nmemb, msg)
}

pub unsafe fn rmalloc(nmemb: usize, msg: &'static CStr) -> NonNull<[real_t]> {
    wrap_xmalloc(nmemb, msg)
}

pub unsafe fn rsmalloc(nmemb: usize, val: real_t, msg: &'static CStr) -> NonNull<[real_t]> {
    wrap_xsmalloc(nmemb, val, msg)
}

/// wrapper over `gk_free` that is polymorphic
pub unsafe fn free_ref<T: 'static + Sized>(p: &mut *mut T) {
    assert!(size_of_val(p) == size_of::<*mut ()>());
    let p: *mut *mut T = p;
    gk_free_one(p.cast());
}

pub struct FreeGuard<'a>(PhantomData<&'a mut ()>);


/// free a pointer using gk_free that we've made a mutable reference from. This obviously still isn't safe, but
/// it should be safe enough to help me spot any potential use-after-frees that show up because we
/// no longer zero the pointer.
pub unsafe fn free_guard<T: 'static + Sized>(p: &'static mut T) -> FreeGuard<'static> {
    let mut praw: *mut T = p;
    gk_free_one((&raw mut praw).cast());
    debug_assert!(praw.is_null());
    FreeGuard(PhantomData)
}


/// like [`free_guard`], but for a slice
pub unsafe fn free_slice_guard<T: 'static + Sized>(p: &'static mut [T]) -> FreeGuard<'static> {
    let mut praw: *mut T = p.as_mut_ptr();
    gk_free_one((&raw mut praw).cast());
    debug_assert!(praw.is_null());
    FreeGuard(PhantomData)
}
