//! Function stubs only used in intermediate stages of translation, to be replaced with `rust-analyzer
//! ssr`
#![allow(dead_code)]
#![allow(unused)]

use crate::idx_t;

#[deprecated]
pub fn gk_max<T: PartialOrd>(a: T, b: T) -> T {
    unimplemented!()
}

#[deprecated]
pub fn gk_min<T: PartialOrd>(a: T, b: T) -> T {
    unimplemented!()
}

#[deprecated]
pub fn iabs(a: idx_t) -> idx_t {
    unimplemented!()
}

#[deprecated]
pub fn isum(n: usize, slice: impl AsRef<[idx_t]>, stride: usize) -> idx_t {
    unimplemented!()
}

#[deprecated]
pub fn imax(n: usize, slice: impl AsRef<[idx_t]>, stride: usize) -> idx_t {
    unimplemented!()
}

#[deprecated]
pub fn iargmax(n: usize, slice: impl AsRef<[idx_t]>, stride: usize) -> usize {
    unimplemented!()
}

#[deprecated]
pub fn iargmin(n: usize, slice: impl AsRef<[idx_t]>, stride: usize) -> usize {
    unimplemented!()
}

#[deprecated]
pub fn iargmax_strd(n: usize, slice: impl AsRef<[idx_t]>, stride: usize) -> usize {
    unimplemented!()
}

#[deprecated]
pub fn iset(n: usize, val: idx_t, slice: impl AsMut<[idx_t]>) -> usize {
    unimplemented!()
}

/// C double, declared in `gk_temp` so we can find all the weird places C inserted it when we're
/// ready to break bug-for-bug backwards compatibility
#[allow(non_camel_case_types)]
pub type double = f64;
