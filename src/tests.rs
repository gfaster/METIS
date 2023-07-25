#![cfg(test)]
#![allow(unused_mut, nonstandard_style)]

use crate::bindings::{
    idx_t, METIS_PartGraphKway, METIS_PartGraphRecursive, METIS_NOPTIONS, METIS_OK,
};

#[test]
fn basic_part_graph_recursive() {
    let mut xadj = &mut [0, 1, 2];
    let mut adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let mut vwgt = std::ptr::null_mut();
    let mut vsize = std::ptr::null_mut();
    let mut adjwgt = std::ptr::null_mut();
    let mut nparts: idx_t = 2;
    let mut tpwgts = std::ptr::null_mut();
    let mut ubvec = std::ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        METIS_PartGraphRecursive(
            &mut nvtxs as *mut _,
            &mut ncon as *mut _,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            vsize,
            adjwgt,
            &mut nparts as *mut _,
            tpwgts,
            ubvec,
            options.as_mut_ptr(),
            &mut objval as *mut _,
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}

#[test]
fn basic_part_graph_kway() {
    let mut xadj = &mut [0, 1, 2];
    let mut adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let mut vwgt = std::ptr::null_mut();
    let mut vsize = std::ptr::null_mut();
    let mut adjwgt = std::ptr::null_mut();
    let mut nparts: idx_t = 2;
    let mut tpwgts = std::ptr::null_mut();
    let mut ubvec = std::ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        METIS_PartGraphKway(
            &mut nvtxs as *mut _,
            &mut ncon as *mut _,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            vsize,
            adjwgt,
            &mut nparts as *mut _,
            tpwgts,
            ubvec,
            options.as_mut_ptr(),
            &mut objval as *mut _,
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}

mod util {
    use crate::util;

    #[ab_test_eq(util::iargmax_nrm)]
    fn iargmax_nrm() -> i32 {
        let x = &[1, 3, 2, 4];
        let y = &[5.0, 3.0, 3.0, 1.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(util::iargmax_nrm)]
    fn iargmax_nrm_same() -> i32 {
        let x = &[1, 2, 3, 4, 6];
        let y = &[12.0, 6.0, 4.0, 3.0, 2.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(util::iargmax2_nrm)]
    fn iargmax2_nrm() -> i32 {
        let x = &[1, 3, 2, 4];
        let y = &[5.0, 3.0, 3.0, 1.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax2_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(util::iargmax2_nrm)]
    fn iargmax2_nrm_same() -> i32 {
        let x = &[1, 2, 3, 4, 6];
        let y = &[12.0, 6.0, 4.0, 3.0, 2.0];
        let n = x.len();

        assert_eq!(x.len(), y.len());

        unsafe { iargmax2_nrm(n, x.as_ptr(), y.as_ptr()) }
    }

    #[ab_test_eq(util::iargmax_strd)]
    fn iargmax_strd_basic() -> i32 {
        let x = &[1, 3, 2, 4];
        let n = x.len();
        let incx = 1;

        assert!((n * incx as usize) <= x.len());

        unsafe { iargmax_strd(n, x.as_ptr(), incx) }
    }

    #[ab_test_eq(util::iargmax_strd)]
    fn iargmax_strd_stride() -> i32 {
        let x = &[1, 3, 2, 4, 3, 2, 1, 5, 1];
        let n = 3;
        let incx = 3;

        assert!((n * incx as usize) <= x.len());

        unsafe { iargmax_strd(n, x.as_ptr(), incx) }
    }

    #[ab_test_basic(util::iargmax_strd)]
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
}
