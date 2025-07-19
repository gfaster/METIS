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

pub mod balance;
pub mod bucketsort;
pub mod checkgraph;
pub mod coarsen;
pub mod compress;
pub mod contig;
pub mod debug;
pub mod fm;
pub mod graph;
pub mod initpart;
pub mod kmetis;
pub mod kwayrefine;
pub mod mcutil;
pub mod minconn;
pub mod pmetis;
pub mod refine;
pub mod ometis;
pub mod mmd;
pub mod separator;
pub mod kwayfm;

pub mod gklib_replace;
pub mod graphio;
pub mod params;
pub mod scanf;

mod dyncall;

pub(crate) mod blas;
pub(crate) mod dal;
pub(crate) mod defs;
pub(crate) mod gk;
pub(crate) mod graph_gen;
pub(crate) mod pqueue;
use defs::*;

#[macro_use]
pub mod util;
pub use util::{AsBufPtrExt, AsNullablePtr, AsNullablePtrMut};

#[cfg(test)]
mod tests;
