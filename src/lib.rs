#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

#[macro_use]
extern crate macros;

pub mod bindings;
pub use bindings::*;

pub mod bucketsort;
pub mod debug_rs;
pub use debug_rs as debug;
pub mod kwayrefine;

#[macro_use]
pub mod util;

#[cfg(test)]
mod tests;
