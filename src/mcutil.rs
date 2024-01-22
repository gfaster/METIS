/*
 * mutil.c
 *
 * This file contains various utility functions for the MOC portion of the
 * code
 *
 * Started 2/15/98
 * George
 *
 * $Id: mcutil.c 13901 2013-03-24 16:17:03Z karypis $
 *
 */

use crate::*;

/*************************************************************************/
/* This function compares two vectors x & y and returns true
    if \forall i, x[i] <= y[i].
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn rvecle(n: idx_t, x: *const real_t, y: *const real_t) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(y, n);
    x.iter().zip(y).rev().all(|(x, y)| x <= y) as std::ffi::c_int
}

/*************************************************************************/
/* This function compares two vectors x & y and returns true
    if \forall i, x[i] >= y[i].
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn rvecge(n: idx_t, x: *const real_t, y: *const real_t) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(y, n);
    x.iter().zip(y).rev().all(|(x, y)| x >= y) as std::ffi::c_int
}

/*************************************************************************/
/* This function compares vectors x1+x2 against y and returns true
    if \forall i, x1[i]+x2[i] <= y[i].
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn rvecsumle(
    n: idx_t,
    x1: *const real_t,
    x2: *const real_t,
    y: *const real_t,
) -> std::ffi::c_int {
    mkslice!(x1, n);
    mkslice!(x2, n);
    mkslice!(y, n);
    x1.iter()
        .zip(x2)
        .map(|(x1, x2)| x1 + x2)
        .zip(y)
        .all(|(x, &y)| x <= y) as _
}

/*************************************************************************/
/* This function returns max_i(x[i]-y[i]) */
/**************************************************************************/
#[metis_func]
pub extern "C" fn rvecmaxdiff(n: idx_t, x: *const real_t, y: *const real_t) -> real_t {
    // real_t max;
    assert!(n > 0);
    mkslice!(x, n);
    mkslice!(y, n);
    x.iter()
        .zip(y)
        .map(|(x, &y)| x - y)
        .max_by(real_t::total_cmp)
        .unwrap()
}

/*************************************************************************/
/* This function returns true if \forall i, x[i] <= z[i]. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ivecle(n: idx_t, x: *const idx_t, z: *const idx_t) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(z, n);
    x.iter().zip(z).rev().all(|(x, z)| x <= z) as std::ffi::c_int
}

/*************************************************************************/
/* This function returns true if \forall i, x[i] >= z[i]. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ivecge(n: idx_t, x: *const idx_t, z: *const idx_t) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(z, n);
    x.iter().zip(z).rev().all(|(x, z)| x >= z) as std::ffi::c_int
}

/*************************************************************************/
/* This function returns true if \forall i, a*x[i]+y[i] <= z[i]. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ivecaxpylez(
    n: idx_t,
    a: idx_t,
    x: *const idx_t,
    y: *const idx_t,
    z: *const idx_t,
) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(y, n);
    mkslice!(z, n);
    x.iter()
        .zip(y)
        .map(|(&x, &y)| a * x + y)
        .zip(z)
        .all(|(l, &z)| l <= z) as _
}

/*************************************************************************/
/* This function returns true if \forall i, a*x[i]+y[i] >= z[i]. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ivecaxpygez(
    n: idx_t,
    a: idx_t,
    x: *const idx_t,
    y: *const idx_t,
    z: *const idx_t,
) -> std::ffi::c_int {
    mkslice!(x, n);
    mkslice!(y, n);
    mkslice!(z, n);
    x.iter()
        .zip(y)
        .map(|(&x, &y)| a * x + y)
        .zip(z)
        .all(|(l, &z)| l >= z) as _
}

/*************************************************************************/
/* This function checks if v+u2 provides a better balance in the weight
vector that v+u1 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn BetterVBalance(
    ncon: idx_t,
    invtvwgt: *const real_t,
    v_vwgt: *const idx_t,
    u1_vwgt: *const idx_t,
    u2_vwgt: *const idx_t,
) -> std::ffi::c_int {
    // idx_t i;
    // real_t sum1=0.0, sum2=0.0, diff1=0.0, diff2=0.0;
    let ncon = ncon as usize;
    mkslice!(v_vwgt, ncon);
    mkslice!(u1_vwgt, ncon);
    mkslice!(u2_vwgt, ncon);
    mkslice!(invtvwgt, ncon);
    let mut sum1 = 0.0;
    let mut sum2 = 0.0;
    for i in (0)..(ncon) {
        sum1 += (v_vwgt[i] + u1_vwgt[i]) as real_t * invtvwgt[i];
        sum2 += (v_vwgt[i] + u2_vwgt[i]) as real_t * invtvwgt[i];
    }
    sum1 = sum1 / ncon as real_t;
    sum2 = sum2 / ncon as real_t;

    let mut diff1 = 0.0;
    let mut diff2 = 0.0;
    for i in (0)..(ncon) {
        diff1 += (sum1 - (v_vwgt[i] + u1_vwgt[i]) as real_t * invtvwgt[i]).abs();
        diff2 += (sum2 - (v_vwgt[i] + u2_vwgt[i]) as real_t * invtvwgt[i]).abs();
    }

    return (diff1 - diff2 >= 0.0) as _;
}

/*************************************************************************/
/* This function takes two ubfactor-centered load imbalance vectors x & y,
and returns true if y is better balanced than x. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn BetterBalance2Way(n: idx_t, x: *const real_t, y: *const real_t) -> std::ffi::c_int {
    // real_t nrm1=0.0, nrm2=0.0;
    mkslice!(x, n);
    mkslice!(y, n);
    let nrm1: real_t = x.iter().filter(|&&x| x > 0.0).sum();
    let nrm2: real_t = y.iter().filter(|&&y| y > 0.0).sum();
    (nrm1 < nrm2) as _
    // for (--n; n>=0; n--) {
    //   if (x[n] > 0) nrm1 += x[n]*x[n];
    //   {
    //     if (y[n] > 0) nrm2 += y[n]*y[n];
    //   }
    // }
    // return nrm2 < nrm1;
}

/*************************************************************************/
/* Given a vertex and two weights, this function returns 1, if the second
    partition will be more balanced than the first after the weighted
    additional of that vertex.
    The balance determination takes into account the ideal target weights
    of the two partitions.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn BetterBalanceKWay(
    ncon: idx_t,
    vwgt: *const idx_t,
    ubvec: *const real_t,
    a1: idx_t,
    pt1: *const idx_t,
    bm1: *const real_t,
    a2: idx_t,
    pt2: *const idx_t,
    bm2: *const real_t,
) -> std::ffi::c_int {
    // idx_t i;
    // real_t tmp, nrm1=0.0, nrm2=0.0, max1=0.0, max2=0.0;
    let ncon = ncon as usize;
    mkslice!(bm1, ncon);
    mkslice!(pt1, ncon);
    mkslice!(bm2, ncon);
    mkslice!(pt2, ncon);
    mkslice!(vwgt, ncon);
    mkslice!(ubvec, ncon);
    let mut nrm1 = 0.0;
    let mut nrm2 = 0.0;
    let mut max1 = 0.0;
    let mut max2 = 0.0;
    for i in (0)..(ncon) {
        let tmp = bm1[i] * (pt1[i] + a1 * vwgt[i]) as real_t - ubvec[i];
        //printf("BB: %d %+.4f ", (int)i, (float)tmp);
        nrm1 += tmp * tmp;
        max1 = if tmp > max1 { tmp } else { max1 };

        let tmp = bm2[i] * (pt2[i] + a2 * vwgt[i]) as real_t - ubvec[i];
        //printf("%+.4f ", (float)tmp);
        nrm2 += tmp * tmp;
        max2 = if tmp > max2 { tmp } else { max2 };

        //printf("%4d %4d %4d %4d %4d %4d %4d %.2f\n",
        //    (int)vwgt[i],
        //    (int)a1, (int)pt1[i], (int)tpt1[i],
        //    (int)a2, (int)pt2[i], (int)tpt2[i], ubvec[i]);
    }
    //printf("   %.3f %.3f %.3f %.3f\n", (float)max1, (float)nrm1, (float)max2, (float)nrm2);

    if max2 < max1 {
        return 1;
    }

    if max2 == max1 && nrm2 < nrm1 {
        return 1;
    }

    return 0;
}

/*************************************************************************/
/* Computes the maximum load imbalance of a partitioning solution over
all the constraints. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ComputeLoadImbalance(
    graph: *const graph_t,
    nparts: idx_t,
    pijbm: *const real_t,
) -> real_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, ncon, *pwgts;
    // real_t max, cur;

    let ncon = graph.ncon as usize;
    let nparts = nparts as usize;
    mkslice!(graph->pwgts, ncon * nparts);
    mkslice!(pijbm, ncon * nparts);

    let mut max = 1.0;
    for i in (0)..(ncon) {
        for j in (0)..(nparts) {
            let cur = pwgts[j * ncon + i] as real_t * pijbm[j * ncon + i];
            if cur > max {
                max = cur;
            }
        }
    }

    return max;
}

/*************************************************************************/
/* Computes the maximum load imbalance difference of a partitioning
   solution over all the constraints.
   The difference is defined with respect to the allowed maximum
   unbalance for the respective constraint.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn ComputeLoadImbalanceDiff(
    graph: *const graph_t,
    nparts: idx_t,
    pijbm: *const real_t,
    ubvec: *const real_t,
) -> real_t {
    // idx_t i, j, ncon, *pwgts;
    // real_t max, cur;

    let graph = graph.as_ref().unwrap();
    let ncon = graph.ncon as usize;
    let nparts = nparts as usize;
    mkslice!(graph->pwgts, ncon * nparts);
    mkslice!(pijbm, ncon * nparts);
    mkslice!(ubvec, ncon);

    let mut max = -1.0;
    for i in (0)..(ncon) {
        for j in (0)..(nparts) {
            let cur = pwgts[j * ncon + i] as real_t * pijbm[j * ncon + i] - ubvec[i];
            if cur > max {
                max = cur;
            }
        }
    }

    return max;
}

/*************************************************************************/
/* Computes the difference between load imbalance of each constraint across
the partitions minus the desired upper bound on the load imabalnce.
It also returns the maximum load imbalance across the partitions &
constraints. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ComputeLoadImbalanceDiffVec(
    graph: *const graph_t,
    nparts: idx_t,
    pijbm: *const real_t,
    ubfactors: *const real_t,
    diffvec: *mut real_t,
) -> real_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, ncon, *pwgts;
    // real_t cur, max;

    let ncon = graph.ncon as usize;
    let nparts = nparts as usize;
    mkslice_mut!(diffvec, ncon);
    mkslice!(graph->pwgts, ncon * nparts);
    mkslice!(pijbm, ncon * nparts);
    mkslice!(ubfactors, ncon);

    let mut max = -1.0;
    for i in (0)..(ncon) {
        diffvec[i] = pwgts[i] as real_t * pijbm[i] - ubfactors[i];
        for j in (1)..(nparts) {
            let cur = pwgts[j * ncon + i] as real_t * pijbm[j * ncon + i] - ubfactors[i];
            if cur > diffvec[i] {
                diffvec[i] = cur;
            }
        }
        if max < diffvec[i] {
            max = diffvec[i];
        }
    }

    return max;
}

/*************************************************************************/
/* Computes the load imbalance of each constraint across the partitions. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn ComputeLoadImbalanceVec(
    graph: *const graph_t,
    nparts: idx_t,
    pijbm: *const real_t,
    lbvec: *mut real_t,
) {
    // idx_t i, j, ncon, *pwgts;
    // real_t cur;

    let graph = graph.as_ref().unwrap();
    let ncon = graph.ncon as usize;
    let nparts = nparts as usize;
    mkslice_mut!(lbvec, ncon);
    mkslice!(graph->pwgts, ncon * nparts);
    mkslice!(pijbm, ncon * nparts);

    for i in (0)..(ncon) {
        lbvec[i] = pwgts[i] as real_t * pijbm[i];
        for j in (1)..(nparts) {
            let cur = pwgts[j * ncon + i] as real_t * pijbm[j * ncon + i];
            if cur > lbvec[i] {
                lbvec[i] = cur;
            }
        }
    }
}
