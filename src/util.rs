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

macro_rules! options_match {
    ($options:ident, Cut) => {
        $options[$crate::METIS_OPTION_PTYPE as usize] = $crate::Objtype::Cut as idx_t;
        $options[$crate::METIS_OPTION_NCUTS as usize] = 1;
    };
    ($options:ident, Vol) => {
        $options[$crate::METIS_OPTION_PTYPE as usize] = $crate::Objtype::Vol as idx_t;
    };
    ($options:ident, Node) => {
        $options[$crate::METIS_OPTION_PTYPE as usize] = $crate::Objtype::Node as idx_t;
    };
    ($options:ident, Grow) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Grow as idx_t;
    };
    ($options:ident, Random) => {
        $options[$crate::METIS_OPTION_IPTYPE as usize] = $crate::Iptype::Random as idx_t;
    };
    ($options:ident, Minconn) => {
        $options[$crate::METIS_OPTION_MINCONN as usize] = 1;
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
}

#[macro_export]
macro_rules! get_graph_slices_mut {
    ($ctrl:expr, $graph:expr => $($val:ident)*) => {
        $(
            let $val = std::slice::from_raw_parts_mut($graph.$val, slice_len!($ctrl, $graph, $val) as usize);
        )*
    };
    ($graph:expr => $($val:ident)*) => {
        $(
            let $val = std::slice::from_raw_parts_mut($graph.$val, slice_len!((), $graph, $val) as usize);
        )*
    };
}

#[macro_export]
macro_rules! get_graph_slices {
    ($ctrl:expr, $graph:expr => $($val:ident)*) => {
        $(
            let $val = std::slice::from_raw_parts($graph.$val, slice_len!($ctrl, $graph, $val) as usize);
        )*
    };
    ($graph:expr => $($val:ident)*) => {
        $(
            let $val = std::slice::from_raw_parts($graph.$val, slice_len!((), $graph, $val) as usize);
        )*
    };
}

/// debug level handling from gk_macros.h
#[macro_export]
macro_rules! ifset {
    ($a:expr, $flag:expr, $cmd:expr) => {
        if $a & $flag != 0 {
            $cmd;
        }
    };
    ($a:expr, $flag:expr, $cmd:expr,) => {
        if $a & $flag != 0 {
            $cmd;
        }
    };
}

/// Equivalent of MAKECSR in gk_macros.h
///
/// n is the length of the slice, but it's often used shorter in METIS
#[inline(always)]
pub fn make_csr(n: usize, a: &mut [idx_t]) {
    assert!(n < a.len(), "making a csr indexes up to n");
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

/// Equivalent of SHIFTCSR in gk_macros.h
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

#[macro_export]
macro_rules! mkslice_mut {
    ($struct:ident->$var:ident, $len:expr) => {
        let $var: &mut [_] = std::slice::from_raw_parts_mut((*$struct).$var, $len as usize);
    };
    ($newvar:ident: $struct:ident->$var:ident, $len:expr) => {
        let $newvar: &mut [_] = std::slice::from_raw_parts_mut((*$struct).$var, $len as usize);
    };
    ($var:ident, $len:expr) => {
        let $var: &mut [_] = std::slice::from_raw_parts_mut($var, $len as usize);
    };
    ($newvar:ident: $var:expr, $len:expr) => {
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

/// read a graph from a file in a simplified version of the format specified in the manual
///
/// returns (xadj, adjncy)
pub fn read_graph(f: &mut impl BufRead) -> Result<(Vec<idx_t>, Vec<idx_t>), Box<dyn Error>> {
    let mut buf = String::with_capacity(80);
    f.read_line(&mut buf)?;

    let mut xadj = Vec::<idx_t>::new();
    let mut adjncy = Vec::<idx_t>::new();

    let mut x = 0;
    buf.clear();
    while 0 != f.read_line(&mut buf)? {
        if buf.trim_start().starts_with('#') || buf.trim().is_empty() {
            continue;
        }

        xadj.push(x);

        let split = buf.split_whitespace();
        for a in split {
            adjncy.push(a.parse()?);
            x += 1;
        }

        buf.clear();
    }
    xadj.push(x);

    for (i, (start, end)) in xadj.windows(2).map(|w| (w[0], w[1])).enumerate() {
        assert!((start as usize) < adjncy.len());
        assert!((end as usize) <= adjncy.len());
        assert!(start < end);
        for j in &adjncy[(start as usize)..(end as usize)] {
            assert!(j >= &0, "no negatives");
            assert_ne!(i, *j as usize, "no self loops");
            assert!(*j < xadj.len() as idx_t - 1, "adj in bounds");
        }
    }

    Ok((xadj, adjncy))
}

/// adapted from mtest.c: VerifyPart
pub fn verify_part(
    xadj: &[idx_t],
    adjncy: &[idx_t],
    vwgt: Option<&[idx_t]>,
    adjwgt: Option<&[idx_t]>,
    objval: idx_t,
    part: &[idx_t],
    nparts: idx_t,
) {
    assert!(nparts > 0);
    assert_eq!(
        xadj.len() - 1,
        part.len(),
        "part is the partition that each vertex goes to"
    );

    let nvtxs = xadj.len() - 1;
    let mut pwgts = vec![0; nparts as usize];

    assert_eq!(
        *part.iter().max().unwrap_or(&0),
        nparts - 1,
        "total number of partitions eq to nparts"
    );

    let mut cut = 0;
    for i in 0..nvtxs {
        pwgts[part[i] as usize] += vwgt.map(|v| v[i]).unwrap_or(1);
        for j in xadj[i]..xadj[i + 1] {
            if part[i] != part[adjncy[j as usize] as usize] {
                cut += adjwgt.map(|v| v[j as usize]).unwrap_or(1);
            }
        }
    }

    assert_eq!(
        cut,
        2 * objval,
        "objval should be edgecut, and the calculated cut should be double it"
    );
    let actual = (nparts * pwgts.iter().max().unwrap()) as f64;
    let expected = 1.10 * pwgts.iter().sum::<idx_t>() as f64;
    assert!(
        // (nparts * pwgts[iargmax(nparts, pwgts)]) as f64
        actual <= expected,
        "actual: {actual:.1}, expected: {expected:.1}.\n\tThis assert spuriously fails sometimes - rerun tests."
    );
}

/// creates a set of dummy weights for partition testing (vwgt, adjwgt)
#[cfg(debug_assertions)]
pub fn create_dummy_weights(
    ncon: usize,
    xadj: &[idx_t],
    adjncy: &[idx_t],
) -> (Vec<idx_t>, Vec<idx_t>) {
    let mut rng = fastrand::Rng::new();
    let nvtxs = xadj.len() - 1;

    let vwgt = Vec::from_iter((0..(nvtxs * ncon)).map(|_| rng.i32(0..10)));

    let mut adjwgt = vec![0; adjncy.len()];
    for i in 0..nvtxs {
        for j in xadj[i]..xadj[i + 1] {
            let k = adjncy[j as usize] as usize;
            if i < k {
                adjwgt[j as usize] = rng.i32(1..=5);
                for jj in xadj[k]..xadj[k + 1] {
                    if adjncy[jj as usize] as usize == i {
                        adjwgt[jj as usize] = adjwgt[j as usize];
                        break;
                    }
                }
            }
        }
    }
    (vwgt, adjwgt)
}

/// returns (xadj, adjncy)
pub fn create_dummy_graph(nvtxs: idx_t) -> (Vec<idx_t>, Vec<idx_t>) {
    let mut xadj = Vec::with_capacity(nvtxs as usize + 1);
    let mut adjncy = Vec::with_capacity(nvtxs as usize * 2);
    for x in 0..nvtxs {
        xadj.push(adjncy.len() as idx_t);
        adjncy.push((x + 1) % nvtxs);
        adjncy.push((x + nvtxs - 1) % nvtxs);
    }
    xadj.push(adjncy.len() as idx_t);

    (xadj, adjncy)
}

/// wrapper for [`ctrl_t`] that frees it properly on drop
pub struct Ctrl {
    pub inner: *mut ctrl_t,

    /// for testing purposes, we often need to initialize a [`graph_t`] associated with ctrl. This
    /// is because all allocation that METIS does relies on the memory arena created for a specific
    /// graph. For convienence purposes (for now), I want to be able to dump those structures in
    /// along with [`Ctrl`] so they can be freed properly.
    graph_parts: Vec<Box<dyn std::any::Any>>,
}

impl Ctrl {
    /// calls [`SetupCtrl`] to create
    pub fn new(
        optype: Optype,
        options: &mut [idx_t; METIS_NOPTIONS as usize],
        ncon: idx_t,
        nparts: idx_t,
        tpwgts: Option<&[real_t]>,
        ubvec: Option<&[real_t]>,
    ) -> Self {
        if unsafe { gk_malloc_init() } == 0 {
            panic!("gk_malloc_init failed")
        }

        let tpwgts = tpwgts.map_or(ptr::null(), |x| x.as_ptr());
        let ubvec = ubvec.map_or(ptr::null(), |x| x.as_ptr());

        let inner = unsafe {
            SetupCtrl(
                optype as u32,
                options.as_mut_ptr(),
                ncon,
                nparts,
                tpwgts,
                ubvec,
            )
        };
        if inner.is_null() {
            panic!("setup ctrl failed")
        }
        Ctrl {
            inner,
            graph_parts: Vec::new(),
        }
    }

    pub fn new_kmetis_basic() -> Self {
        let mut options = make_options!(Cut Grow);
        let ncon = 1;
        let nparts = 2;
        let tpwgts = None;
        let ubvec = None;
        let optype = Optype::Kmetis;
        Self::new(optype, &mut options, ncon, nparts, tpwgts, ubvec)
    }

    pub fn init_dummy_graph(&mut self, nvtxs: idx_t) {
        let (mut xadj, mut adjncy) = create_dummy_graph(nvtxs);
        let graph = unsafe {
            SetupGraph(
                self.inner,
                nvtxs,
                1,
                xadj.as_mut_ptr(),
                adjncy.as_mut_ptr(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
            )
        };
        self.graph_parts.push(Box::new(adjncy));
        self.graph_parts.push(Box::new(xadj));
        unsafe { AllocateWorkSpace(self.inner, graph) };
    }
}

impl Drop for Ctrl {
    fn drop(&mut self) {
        // core (worksapce) is freed if it has been set: noop on null

        if self.inner.is_null() {
            panic!("dropped null ctrl")
        } else {
            unsafe { FreeCtrl(&mut self.inner as *mut _) };
        }
    }
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
