#[macro_use]
extern crate macros;

pub mod bindings;
pub use bindings::*;

pub mod bucketsort;
pub mod kwayrefine;

#[macro_use]
pub mod util;

#[cfg(test)]
mod tests;
