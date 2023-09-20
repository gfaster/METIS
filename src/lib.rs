#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::missing_safety_doc
)]

#[macro_use]
extern crate macros;

pub mod bindings;
pub use bindings::*;

pub mod bucketsort;
pub mod debug_rs;
pub use debug_rs as debug;
pub mod balance;
pub mod initpart;
pub mod kmetis;
pub mod kwayrefine;

pub(crate) mod blas;
pub(crate) mod pqueue;

#[macro_use]
pub mod util;

#[cfg(test)]
mod tests;
