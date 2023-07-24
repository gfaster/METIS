
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c 10495 2011-07-06 16:04:45Z karypis $
 */


use super::macros::*;
use super::bindings::*;
use std::ffi::{c_void};



/// Initialize the random number generator
#[metis_func]
extern "C" fn InitRandom(seed: idx_t) -> c_void {
  unsafe {
    isrand(if seed == -1 { 4321 } else { seed })
  }
}


/// index of heighest weight of x[i] * y[i]
#[metis_func]
extern "C" fn iargmax_nrm(n: usize, x: *const idx_t, y: *const real_t) -> idx_t {
  let x = unsafe { std::slice::from_raw_parts(x, n) };
  let y = unsafe { std::slice::from_raw_parts(y, n) };

  let mut max = 0;
  for i in 1..n {
    if x[i] as real_t * y[i] > x[max] as real_t * y[max] {
      max = i;
    }
  }

  max as idx_t
}


/// return the index of the maximum element of a vector
#[metis_func]
extern "C" fn iargmax_strd(n: usize, x: *const idx_t, incx: idx_t) -> idx_t {
  let incx = incx as usize;
  let n = n * incx;
  let x = unsafe { std::slice::from_raw_parts(x, n) };

  let mut max = 0;
  for i in (incx..n).step_by(incx) {
    if x[i] > x[max] {
      max = i;
    }
  }
  (max / incx) as idx_t
}


/// return the index of the almost max element in a vector
#[metis_func]
extern "C" fn rargmax2(n: usize, x: *const real_t) -> idx_t {
  let x = unsafe { std::slice::from_raw_parts(x, n) };
  let mut max1;
  let mut max2;

  if x[0] > x[1] {
    max1 = 0;
    max2 = 1;
  } else {
    max1 = 0;
    max2 = 0;
  }

  for i in 2..n {
    if x[i] > x[max1] {
      max2 = max1;
      max1 = i as usize;
    } else if x[i] > x[max2] {
      max2 = i as usize;
    }
  }
  max2 as idx_t
}



/// return the index of the second largest x[i] * y[i]
#[metis_func]
extern "C" fn iargmax2_nrm(n: usize, x: *const idx_t, y: *const real_t) -> idx_t {
  let x = unsafe { std::slice::from_raw_parts(x, n) };
  let y = unsafe { std::slice::from_raw_parts(y, n) };

  let mut max1;
  let mut max2;

  if x[0] as real_t * y[0] > x[1] as real_t * y[1] {
    max1 = 0;
    max2 = 1;
  } else {
    max1 = 1;
    max2 = 0;
  }

  for i in 2..n {
    if x[i] as real_t * y[i] > x[max1] as real_t * y[max1] {
      max2 = max1;
      max1 = i;
    } else if x[i] as real_t * y[i] > x[max2] as real_t * y[max2] {
      max2 = i;
    }
  }
  max2 as idx_t
}



/// converts a signal code to a metis return code
#[metis_func]
extern "C" fn metis_rcode(sigrval: std::ffi::c_int) -> std::ffi::c_int {
  match sigrval {
    0 => METIS_OK,
    SIGMEM => METIS_ERROR_MEMORY,
    _ => METIS_ERROR
  }
}
