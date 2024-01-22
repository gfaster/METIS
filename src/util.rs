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

use super::bindings::*;
use std::io::BufRead;
use std::ptr;
use std::{error::Error, ffi::c_void};

/// Initialize the random number generator
#[metis_func]
pub extern "C" fn InitRandom(seed: idx_t) -> c_void {
    unsafe { isrand(if seed == -1 { 4321 } else { seed }) }
}

/// index of heighest weight of x[i] * y[i]
#[metis_func]
pub extern "C" fn iargmax_nrm(n: usize, x: *const idx_t, y: *const real_t) -> idx_t {
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
///
/// In reality we go by a stride
///
/// # Safety
///
/// `n * incx` must be less than the length of x
#[metis_func]
pub extern "C" fn iargmax_strd(n: usize, x: *const idx_t, incx: idx_t) -> idx_t {
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
pub extern "C" fn rargmax2(n: usize, x: *const real_t) -> idx_t {
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
            max1 = i;
        } else if x[i] > x[max2] {
            max2 = i;
        }
    }
    max2 as idx_t
}

/// return the index of the second largest x[i] * y[i]
#[metis_func]
pub extern "C" fn iargmax2_nrm(n: usize, x: *const idx_t, y: *const real_t) -> idx_t {
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
pub extern "C" fn metis_rcode(sigrval: std::ffi::c_int) -> std::ffi::c_int {
    match sigrval {
        0 => METIS_OK,
        SIGMEM => METIS_ERROR_MEMORY,
        _ => METIS_ERROR,
    }
}

#[no_mangle]
pub unsafe extern "C" fn trigger_panic(msg: *const std::ffi::c_char) -> ! {
    let s = std::ffi::CStr::from_ptr(msg);
    let s = s.to_string_lossy();
    panic!("{}", s);
}

/// returns the max element in slice by stride incx
///
/// eventually will be `metis_func`, but need to port everything from the `gk_mkblas` at once
///
/// In the original, `iargmax` takes an `n` parameter which is the number of steps that will be
/// taken. This is implied by the length of the slice
///
/// ```
/// # use metis::util::iargmax;
/// let vals = &[2, 3, 4, 3, 1];
/// assert_eq!(iargmax(vals, 1), 2);
/// assert_eq!(iargmax(vals, 2), 1);
/// assert_eq!(iargmax(vals, 3), 1);
/// assert_eq!(iargmax(vals, 4), 0);
/// ```
#[inline(always)]
pub fn iargmax(x: &[idx_t], incx: usize) -> usize {
    let mut max = 0;

    for j in (incx..x.len()).step_by(incx) {
        if x[j] > x[max] {
            max = j;
        }
    }
    max / incx
}

/// converts a user provided ufactor into a real ubfactor
#[inline(always)]
pub fn i2rubfactor(ufactor: idx_t) -> real_t {
    1.0 + 0.001 * (ufactor as f32)
}

/// Increment and decrement by the same value
///
/// ```
/// # use metis::inc_dec;
/// let mut a = 0;
/// let mut b = 0;
/// inc_dec!(a, b, 3);
/// assert_eq!(a, 3);
/// assert_eq!(b, -3);
/// ```
#[macro_export]
macro_rules! inc_dec {
    ($a:expr, $b:expr, $val:expr) => {
        $a += $val;
        $b -= $val;
    };
}

/*
* ====================================================
* below here is my additional functions
* ====================================================
*/

#[allow(unused_macros)]
macro_rules! options_match {
    ($options:ident, Cut) => {
        $options[$crate::METIS_OPTION_OBJTYPE as usize] = $crate::Objtype::Cut as idx_t;
        $options[$crate::METIS_OPTION_NCUTS as usize] = 1;
    };
    ($options:ident, Vol) => {
        $options[$crate::METIS_OPTION_OBJTYPE as usize] = $crate::Objtype::Vol as idx_t;
    };

    ($options:ident, Grow) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Grow as idx_t;
    };
    ($options:ident, Random) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Random as idx_t;
    };
    ($options:ident, Edge) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Edge as idx_t;
    };
    ($options:ident, Node) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Node as idx_t;
    };
    ($options:ident, Rb) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Rb as idx_t;
    };

    // mostly unused - I can add it back later
    // ($options:ident, Greedy) => {
    //     $options[$crate::METIS_OPTION_RTYPE as usize] = $crate::Rtype::Greedy as idx_t;
    // };
    // ($options:ident, Fm) => {
    //     $options[$crate::METIS_OPTION_RTYPE as usize] = $crate::Rtype::Fm as idx_t;
    // };
    // ($options:ident, Sep1) => {
    //     $options[$crate::METIS_OPTION_RTYPE as usize] = $crate::Rtype::Sep1Sided as idx_t;
    // };
    // ($options:ident, Sep2) => {
    //     $options[$crate::METIS_OPTION_RTYPE as usize] = $crate::Rtype::Sep2Sided as idx_t;
    // };
    ($options:ident, Rm) => {
        $options[$crate::METIS_OPTION_CTYPE as usize] = $crate::Ctype::Rm as idx_t;
    };
    ($options:ident, Shem) => {
        $options[$crate::METIS_OPTION_CTYPE as usize] = $crate::Ctype::Shem as idx_t;
    };

    ($options:ident, Minconn) => {
        $options[$crate::METIS_OPTION_MINCONN as usize] = 1;
    };

    ($options:ident, Contig) => {
        $options[$crate::METIS_OPTION_CONTIG as usize] = 1;
    };

    ($options:ident, None) => {
        // dummy option for test macros
    };
}

/// Build options array: populate it with a number of idents
///
/// `Vol`: communication volume minimization
/// `Cut`: edgecut minimization
///
/// `Grow`: initial partition via greedy strategy
/// `Edge`: make separator from edge cut
/// `Node`: initial partition via greedy node-based strategy
///
/// `Minconn`: minimize the degree of partition graph
///
/// *I need to add more - check manual*
#[macro_export]
macro_rules! make_options {
    ($($val:ident)*) => {{
        let mut options = [-1; METIS_NOPTIONS as usize];
        options[$crate::METIS_OPTION_SEED as usize] = 0;
        options[$crate::METIS_OPTION_ONDISK as usize] = 1;
        $(
            options_match!(options, $val);
        )*
        options
    }};
}

#[macro_export]
macro_rules! slice_len {
    ($ctrl:expr, $graph:expr, where_) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, xadj) => {
        $graph.nvtxs + 1
    };
    ($ctrl:expr, $graph:expr, adjncy) => {{
        let xadj = std::slice::from_raw_parts($graph.xadj, $graph.nvtxs as usize + 1);
        xadj[xadj.len() - 1]
    }};
    ($ctrl:expr, $graph:expr, adjwgt) => {{
        // let xadj = std::slice::from_raw_parts($graph.xadj, $graph.nvtxs as usize + 1);
        // xadj[xadj.len() - 1]
        $graph.nedges
    }};
    ($ctrl:expr, $graph:expr, vwgt) => {
        $graph.nvtxs * $graph.ncon
    };
    ($ctrl:expr, $graph:expr, pwgts) => {
        $ctrl.nparts * $graph.ncon
    };
    ($ctrl:expr, $graph:expr, cmap) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, id) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, ed) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, vsize) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, vkrinfo) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, ckrinfo) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, bndptr) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, bndind) => {
        $graph.nvtxs
    };
    ($ctrl:expr, $graph:expr, tvwgt) => {
        $graph.ncon
    };
    ($ctrl:expr, $graph:expr, invtvwgt) => {
        $graph.ncon
    };
    ($ctrl:expr, $graph:expr, label) => {
        $graph.nvtxs
    };
}

#[macro_export]
macro_rules! get_graph_slices_mut {
    ($ctrl:expr, $graph:expr => $($val:ident)*) => {
        $(
            assert!(!$graph.$val.is_null());
            let $val = std::slice::from_raw_parts_mut($graph.$val, slice_len!($ctrl, $graph, $val) as usize);
        )*
    };
    ($graph:expr => $($val:ident)*) => {
        $(
            assert!(!$graph.$val.is_null());
            let $val = std::slice::from_raw_parts_mut($graph.$val, slice_len!((), $graph, $val) as usize);
        )*
    };
}

#[macro_export]
macro_rules! get_graph_slices {
    ($ctrl:expr, $graph:expr => $($val:ident)*) => {
        $(
            assert!(!$graph.$val.is_null());
            let $val = std::slice::from_raw_parts($graph.$val, slice_len!($ctrl, $graph, $val) as usize);
        )*
    };
    ($graph:expr => $($val:ident)*) => {
        $(
            assert!(!$graph.$val.is_null());
            let $val = std::slice::from_raw_parts($graph.$val, slice_len!((), $graph, $val) as usize);
        )*
    };
}

#[macro_export]
macro_rules! free_field {
    ($struct:ident.$field:ident) => {{
        gk_free_one(std::ptr::addr_of_mut!((*$struct).$field) as *mut *mut std::ffi::c_void);
        assert_eq!((*$struct).$field, std::ptr::null_mut());
    }};
}

/// debug level handling from `gk_macros.h`
#[macro_export]
macro_rules! ifset {
    ($a:expr, $flag:expr, $cmd:expr $(,)?) => {
        if $a & $flag != 0 {
            $cmd;
        }
    };
}

/// Makes an index range from start and count. Also casts to usize.
///
/// ```
/// # use metis::cntrng;
/// assert_eq!(cntrng!(3, 2), 3..5);
/// ```
#[macro_export]
macro_rules! cntrng {
    ($start:expr, $cnt:expr) => {{
        let start = $start as usize;
        start..(start + $cnt as usize)
    }};
}

/// Equivalent of `MAKECSR` in `gk_macros.h`
///
/// n is the length of the slice, but it's often used shorter in METIS
///
/// Converts a slice of lengths (vertex degrees) to be the xadj in a CSR representation. The last
/// index is ignored.
///
/// ```rust
/// # use metis::util::make_csr;
/// let mut a = [3, 2, 4, 3, 1, -999];
///
/// make_csr(a.len() - 1, &mut a);
///
/// assert_eq!(a, [0, 3, 5, 9, 12, 13]);
/// ```
///
#[inline(always)]
pub fn make_csr(n: usize, a: &mut [idx_t]) {
    assert!(
        n < a.len(),
        "making a csr indexes up to n: {n} (a.len = {})",
        a.len()
    );
    debug_assert_eq!(n, a.len() - 1, "I want to see if this ever happens - this assert can be removed. If it never triggers, then we can remove n as an argument");
    if n == 0 {
        return;
    }

    for i in 1..n {
        a[i] += a[i - 1];
    }
    for i in (1..=n).rev() {
        a[i] = a[i - 1];
    }
    a[0] = 0;
}

/// Equivalent of `SHIFTCSR` in `gk_macros.h`
///
/// Shifts every element in the slice over by one, discarding the last element.
///
/// ```rust
/// # use metis::util::shift_csr;
/// let mut a = [3, 1, 2];
///
/// shift_csr(a.len() - 1, &mut a);
///
/// assert_eq!(a, [0, 3, 1]);
#[inline(always)]
pub fn shift_csr(n: usize, a: &mut [idx_t]) {
    assert!(n < a.len(), "making a csr indexes up to n");
    debug_assert_eq!(n, a.len() - 1, "I want to see if this ever happens - this assert can be removed. If it never triggers, then we can remove n as an argument");

    if n == 0 {
        return;
    }
    for i in (1..=n).rev() {
        a[i] = a[i - 1];
    }
    a[0] = 0;
}

/// from mcutil.c
///
/// returns true if forall `i`, `a * x[i] + y[i] <= z[i]`
/// original took length argument at beginning
pub fn ivecaxpylez(a: idx_t, x: &[idx_t], y: &[idx_t], z: &[idx_t]) -> bool {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), z.len());
    x.into_iter()
        .zip(y)
        .map(|(&xi, &yi)| a * xi + yi)
        .zip(z)
        .all(|(li, &zi)| li <= zi)
}

/// from mcutil.c
///
/// returns true if `x[i] <= z[i]` forall `i`
/// original took length argument at beginning
pub fn ivecle(x: &[idx_t], z: &[idx_t]) -> bool {
    assert_eq!(x.len(), z.len());
    x.into_iter().zip(z).all(|(xi, zi)| xi <= zi)
}

#[macro_export]
macro_rules! BNDInsert {
    ($n:expr, $bndind:expr, $bndptr:expr, $i:expr) => {
        assert_eq!($bndptr[$i], -1);
        $bndind[$n as usize] = $i as _;
        $bndptr[$i as usize] = $n as _;
        $n += 1;
    };
}

#[macro_export]
macro_rules! BNDDelete {
    ($nbnd:expr, $bndind:expr, $bndptr:expr, $i:expr) => {
        assert_ne!($bndptr[$i], -1);
        $nbnd -= 1;
        $bndind[$bndptr[$i as usize] as usize] = $bndind[$nbnd as usize];
        $bndptr[$bndind[$nbnd as usize] as usize] = $bndptr[$i as usize];
        $bndptr[$i] = -1;
    };
}

/// note the original mangles 'j'
#[macro_export]
macro_rules! UpdateMovedVertexInfoAndBND {
    ($i:expr, $from:expr, $k:expr, $to:expr, $myrinfo:expr, $mynbrs:expr, $where:expr, $nbnd:expr,
    $bndptr:expr, $bndind:expr, $bndtype:expr) => {
        $where[$i] = $to;
        $myrinfo.ed += $myrinfo.id - $mynbrs[$k].ed;
        std::mem::swap(&mut $myrinfo.id, &mut $mynbrs[$k].ed);
        if ($mynbrs[$k].ed == 0) {
            $myrinfo.nnbrs -= 1;
            $mynbrs[$k as usize] = $mynbrs[$myrinfo.nnbrs as usize];
        } else {
            $mynbrs[$k].pid = $from;
        }

        /* Update the boundary information. Both deletion and addition is
        allowed as this routine can be used for moving arbitrary nodes. */
        if ($bndtype == BNDTYPE_REFINE) {
            if ($bndptr[$i] != -1 && $myrinfo.ed - $myrinfo.id < 0) {
                BNDDelete!($nbnd, $bndind, $bndptr, $i);
            }
            if ($bndptr[$i] == -1 && $myrinfo.ed - $myrinfo.id >= 0) {
                BNDInsert!($nbnd, $bndind, $bndptr, $i);
            }
        } else {
            if ($bndptr[$i] != -1 && $myrinfo.ed <= 0) {
                BNDDelete!($nbnd, $bndind, $bndptr, $i);
            }
            if ($bndptr[$i] == -1 && $myrinfo.ed > 0) {
                BNDInsert!($nbnd, $bndind, $bndptr, $i);
            }
        }
    };
}

/// note the original mangles 'j'
#[macro_export]
macro_rules! UpdateAdjacentVertexInfoAndBND {
    ($ctrl:expr, $vid:expr, $adjlen:expr, $me:expr, $from:expr, $to:expr, $myrinfo:expr, $ewgt:expr, $nbnd:expr, $bndptr:expr, $bndind:expr, $bndtype:expr) => {
        // idx_t k;
        // cnbr_t *mynbrs;

        if ($myrinfo.inbr == -1) {
            $myrinfo.inbr = cnbrpoolGetNext($ctrl, $adjlen);
            $myrinfo.nnbrs = 0;
        }
        assert!(CheckRInfo($ctrl, $myrinfo) != 0);

        let mynbrs = std::slice::from_raw_parts_mut(
            $ctrl.cnbrpool.add($myrinfo.inbr as usize),
            $myrinfo.nnbrs as usize + 1,
        );

        /* Update global ID/ED and boundary */
        if ($me == $from) {
            inc_dec!($myrinfo.ed, $myrinfo.id, ($ewgt));
            if ($bndtype == BNDTYPE_REFINE) {
                if ($myrinfo.ed - $myrinfo.id >= 0 && $bndptr[($vid)] == -1) {
                    BNDInsert!($nbnd, $bndind, $bndptr, ($vid));
                }
            } else {
                if ($myrinfo.ed > 0 && $bndptr[($vid)] == -1) {
                    BNDInsert!($nbnd, $bndind, $bndptr, ($vid));
                }
            }
        } else if ($me == $to) {
            inc_dec!($myrinfo.id, $myrinfo.ed, ($ewgt));
            if ($bndtype == BNDTYPE_REFINE) {
                if ($myrinfo.ed - $myrinfo.id < 0 && $bndptr[($vid)] != -1) {
                    BNDDelete!($nbnd, $bndind, $bndptr, ($vid));
                }
            } else {
                if ($myrinfo.ed <= 0 && $bndptr[($vid)] != -1) {
                    BNDDelete!($nbnd, $bndind, $bndptr, ($vid));
                }
            }
        }

        /* Remove contribution $from the .ed of '$from' */
        if ($me != $from) {
            for k in 0..($myrinfo.nnbrs as usize) {
                if (mynbrs[k].pid == $from) {
                    if (mynbrs[k].ed == ($ewgt)) {
                        $myrinfo.nnbrs -= 1;
                        mynbrs[k] = mynbrs[$myrinfo.nnbrs as usize];
                    } else {
                        mynbrs[k].ed -= ($ewgt);
                    }
                    break;
                }
            }
        }

        /* Add contribution $to the .ed of '$to' */
        if ($me != $to) {
            let mut k = 0;
            while k < $myrinfo.nnbrs as usize {
                if (mynbrs[k].pid == $to) {
                    mynbrs[k].ed += ($ewgt);
                    break;
                }
                k += 1;
            }
            if (k == $myrinfo.nnbrs as usize) {
                mynbrs[k].pid = $to;
                mynbrs[k].ed = ($ewgt);
                $myrinfo.nnbrs += 1;
            }
        }

        assert!(CheckRInfo($ctrl, $myrinfo) != 0);
    };
}

#[macro_export]
macro_rules! mkslice_mut {
    ($struct:ident->$var:ident, $len:expr) => {
        debug_assert!(!(*$struct).$var.is_null());
        let $var: &mut [_] = std::slice::from_raw_parts_mut((*$struct).$var, $len as usize);
    };
    ($newvar:ident: $struct:ident->$var:ident, $len:expr) => {
        debug_assert!(!(*$struct).$var.is_null());
        let $newvar: &mut [_] = std::slice::from_raw_parts_mut((*$struct).$var, $len as usize);
    };
    ($var:ident, $len:expr) => {
        debug_assert!(!$var.is_null());
        let $var: &mut [_] = std::slice::from_raw_parts_mut($var, $len as usize);
    };
    ($newvar:ident: $var:expr, $len:expr) => {
        debug_assert!(!$var.is_null());
        let $newvar: &mut [_] = std::slice::from_raw_parts_mut($var, $len as usize);
    };
}

#[macro_export]
macro_rules! mkslice {
    ($struct:ident->$var:ident, $len:expr) => {
        let $var: &[_] = std::slice::from_raw_parts((*$struct).$var, $len as usize);
    };
    ($newvar:ident: $struct:ident->$var:ident, $len:expr) => {
        let $newvar: &[_] = std::slice::from_raw_parts((*$struct).$var, $len as usize);
    };
    ($var:ident, $len:expr) => {
        let $var: &[_] = std::slice::from_raw_parts($var, $len as usize);
    };
    ($newvar:ident: $var:ident, $len:expr) => {
        let $newvar: &[_] = std::slice::from_raw_parts($var, $len as usize);
    };
}

/// makes a slice or initialize a vec with a default value
/// ```
/// # use metis::slice_default;
/// # unsafe {
/// let v = vec![0, 1, 2];
/// let len = v.len();
/// let arr_v: *const i32 = v.as_ptr();
/// let arr: *mut i32 = std::ptr::null_mut();
/// slice_default!(arr_v, [0; len]);
/// slice_default!(mut arr, [1; len]);
///
/// assert_eq!(arr_v, &[0, 1, 2]);
/// assert_eq!(arr, &mut [1, 1, 1]);
/// # }
/// ```
#[macro_export]
macro_rules! slice_default {
    ($arr:ident, [$fill:expr; $len:expr]) => {
        let arr_v;
        let $arr = if $arr.is_null() {
            arr_v = vec![0; $len as usize];
            &arr_v[..]
        } else {
            std::slice::from_raw_parts($arr, $len as usize)
        };
    };

    ($arr:ident, [$fill:expr; $len:expr]) => {
        let arr_v;
        let $arr = if $arr.is_null() {
            arr_v = vec![$fill; $len as usize];
            &arr_v[..]
        } else {
            std::slice::from_raw_parts($arr, $len as usize)
        };
    };

    (mut $arr:ident, [$fill:expr; $len:expr]) => {
        let mut arr_v;
        let $arr = if $arr.is_null() {
            arr_v = vec![$fill; $len as usize];
            &mut arr_v[..]
        } else {
            std::slice::from_raw_parts_mut($arr, $len as usize)
        };
    };
}
/*

/// makes a slice or initialize a vec with a default value
/// ```
/// let v = vec![0, 1, 2];
/// let arr_v: *mut i32 = v.as_ptr_mut();
/// let arr: *mut i32 = std::ptr::null_mut();
/// let len = v.len();
///
/// slice_default_mut!(arr_v, len => 0);
/// slice_default_mut!(arr, len => 1);
///
/// assert_eq!(arr_v, &mut [0, 1, 2]);
/// assert_eq!(arr, &mut [1, 1, 1]);
/// ```
#[macro_export]
macro_rules! slice_default_mut {
    // ($arr:ident, $len:expr) => {
    //     let mut arr_v;
    //     let $arr = if $arr.is_null() {
    //         arr_v = vec![0; $len as usize];
    //         &mut arr_v[..]
    //     } else {
    //         std::slice::from_raw_parts_mut($arr, $len as usize)
    //     };
    // };

    ($arr:ident, $len:expr => $fill:expr) => {
        let mut arr_v;
        let $arr = if $arr.is_null() {
            arr_v = vec![$fill; $len as usize];
            &mut arr_v[..]
        } else {
            std::slice::from_raw_parts_mut($arr, $len as usize)
        };
    };
}
    */

/// inverse of `util::make_csr`. Last element is unspecified.
///
/// ```rust
/// # use metis::util::from_csr;
///
/// let mut a = [0, 2, 4, 7, 8];
///
/// from_csr(&mut a);
///
/// assert_eq!(&a[..a.len() - 1], &[2, 2, 3, 1]);
/// ```
pub fn from_csr(a: &mut [idx_t]) {
    if a.is_empty() {
        return;
    }
    for i in 0..(a.len() - 1) {
        a[i] = a[i + 1] - a[i];
        assert!(a[i] >= 0, "a[{i}] = {} has negative degree", a[i]);
    }
}

/// inverse of `util::make_csr`, but returns an iterator
///
/// ```rust
/// # use metis::util::from_csr_iter;
///
/// let a = [0, 2, 4, 7, 8];
///
/// let mut it = from_csr_iter(&a);
/// assert_eq!(it.next(), Some(2));
/// assert_eq!(it.next(), Some(2));
/// assert_eq!(it.next(), Some(3));
/// assert_eq!(it.next(), Some(1));
/// assert_eq!(it.next(), None);
/// ```
pub fn from_csr_iter(a: &[idx_t]) -> impl Iterator<Item = idx_t> + '_ {
    a.windows(2).map(|s| {
        let deg = s[1] - s[0];
        debug_assert!(deg >= 0, "deg = {deg} is negative");
        deg
    })
}

#[cfg(test)]
mod test {
    #[allow(unused_imports)]
    use crate::idx_t;

    #[ab_test_eq(super::iargmax_nrm)]
    fn iargmax_nrm() -> idx_t {
        let x = &[1, 3, 2, 4];
        let y = &[5.0, 3.0, 3.0, 1.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(super::iargmax_nrm)]
    fn iargmax_nrm_same() -> idx_t {
        let x = &[1, 2, 3, 4, 6];
        let y = &[12.0, 6.0, 4.0, 3.0, 2.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(super::iargmax2_nrm)]
    fn iargmax2_nrm() -> idx_t {
        let x = &[1, 3, 2, 4];
        let y = &[5.0, 3.0, 3.0, 1.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax2_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(super::iargmax2_nrm)]
    fn iargmax2_nrm_same() -> idx_t {
        let x = &[1, 2, 3, 4, 6];
        let y = &[12.0, 6.0, 4.0, 3.0, 2.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax2_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(super::iargmax_strd)]
    fn iargmax_strd_basic() -> i32 {
        let x = &[1, 3, 2, 4];
        let n = x.len();
        let incx = 1;

        assert!((n * incx as usize) <= x.len());

        unsafe { iargmax_strd(n, x.as_ptr(), incx) }
    }

    #[ab_test_eq(super::iargmax_strd)]
    fn iargmax_strd_stride() -> i32 {
        let x = &[1, 3, 2, 4, 3, 2, 1, 5, 1];
        let n = 3;
        let incx = 3;

        assert!((n * incx as usize) <= x.len());

        unsafe { iargmax_strd(n, x.as_ptr(), incx) }
    }

    #[ab_test_basic(super::iargmax_strd)]
    fn iargmax_strd_stride_expl() {
        let x = &[1, 3, 2, 4, 3, 2, 5, 2, 1];
        let n = 3;
        let incx = 3;
        assert!((n * incx as usize) == x.len());
        assert_eq!(unsafe { iargmax_strd(n, x.as_ptr(), incx) }, 2);

        let n = 4;
        let incx = 2;
        assert!((n * incx as usize) <= x.len());
        assert_eq!(unsafe { iargmax_strd(n, x.as_ptr(), incx) }, 3);
    }

    #[ab_test_eq(super::rargmax2)]
    fn rargmax2() -> idx_t {
        let x = &[12.0, 6.0, 4.0, 3.0, 2.0];
        let n = x.len();

        unsafe { rargmax2(n, x.as_ptr()) }
    }
}
