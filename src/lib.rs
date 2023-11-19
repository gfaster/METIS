#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::missing_safety_doc,
    unused_imports
)]

// use logging_allocator::LoggingAllocator;
// #[global_allocator]
// static ALLOC: LoggingAllocator = LoggingAllocator::new(true);

#[macro_use]
extern crate macros;

pub mod bindings;
pub use bindings::*;

pub mod bucketsort;
pub mod debug_rs;
pub use debug_rs as debug;
pub mod balance;
pub mod coarsen;
pub mod compress;
pub mod contig;
pub mod initpart;
pub mod kmetis;
pub mod kwayrefine;
pub mod pmetis;

pub(crate) mod blas;
pub(crate) mod dal;
pub(crate) mod defs;
pub(crate) mod graph_gen;
pub(crate) mod pqueue;
use defs::*;

#[macro_use]
pub mod util;

#[cfg(test)]
mod tests;
