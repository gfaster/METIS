use crate::idx_t;

/// ?axpy routine - standard to Level 1 BLAS routines
///
/// computes: y + alpha * x where x,y are vectors
#[inline]
pub fn iaxpy<'a>(
    n: usize,
    alpha: idx_t,
    x: &'_ [idx_t],
    incx: idx_t,
    y: &'a mut [idx_t],
    incy: idx_t,
) -> &'a mut [idx_t] {
    assert!(n <= x.len(), "x len: {}, n: {n}", x.len());
    assert!(n <= y.len(), "y len: {}, n: {n}", y.len());

    for (x, y) in x
        .iter()
        .step_by(incx as usize)
        .zip(y.iter_mut().step_by(incy as usize))
        .take(n)
    {
        *y += *x * alpha;
    }

    y
}
