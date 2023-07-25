#[macro_use]
extern crate macros;

pub mod bindings;
pub use bindings::*;

pub mod bucketsort;
pub mod util;

#[cfg(test)]
mod tests;
