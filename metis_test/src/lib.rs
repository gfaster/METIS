#![allow(unused, nonstandard_style)]
mod metis {
    include!("../../libmetis/metis-sys.rs");
}
use metis::*;

#[test]
fn it_works() {
    let mut adjncy: [idx_t; 2] = [1, 0];
    let mut xadj: [idx_t; 3] = [0, 1, 2];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon = 1;
    let mut vwgt = std::ptr::null_mut();
    let mut vsize = std::ptr::null_mut();
    let mut adjwgt = std::ptr::null_mut();
    let mut nparts = 2;
    let mut tpwgts = std::ptr::null_mut();
    let mut ubvec = std::ptr::null_mut();
    let mut options = std::ptr::null_mut();
    let mut objval = 0;
    let mut part = [0; 2];

    let res = unsafe {
        METIS_PartGraphKway(
            std::ptr::addr_of_mut!(nvtxs),
            std::ptr::addr_of_mut!(ncon),
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            vsize,
            adjwgt,
            std::ptr::addr_of_mut!(nparts),
            tpwgts,
            ubvec,
            options,
            std::ptr::addr_of_mut!(objval),
            part.as_mut_ptr()
            )
    };

    assert_eq!(res, rstatus_et_METIS_OK);
    assert_ne!(part[0], part[1]);
}
