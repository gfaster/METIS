#[allow(unused_imports)]
// use metis::util::create_dummy_weights;
use metis::*;

// macro_rules! options_arr {
//     ($objtype:ident $iptype:ident) => {{
//         let mut options = [-1; METIS_NOPTIONS as usize];
//         options[crate::METIS_OPTION_PTYPE as usize] = crate::Objtype::$objtype as i32;
//         options[crate::METIS_OPTION_IPTYPE as usize] = crate::Iptype::$iptype as i32;
//
//         options
//     }};
// }

fn main() {
    // let n = 10;
    // let mut perm = vec![-1; 10];
    // let tperm = vec![0, 9, 3, 2, 1, 8, 4, 5, 6, 7];
    // let keys = vec![0, 5, 3, 2, 0, 8, 5, 2, 8, 7];
    // let ctrl = util::Ctrl::new_kmetis_basic();
    // let max = 10;
    //
    // assert_eq!(perm.len(), n);
    // assert_eq!(tperm.len(), n);
    // assert_eq!(keys.len(), n);
    //
    // unsafe { bucketsort::libmetis__BucketSortKeysInc(ctrl.inner, n as i32, max, keys.as_ptr(), tperm.as_ptr(), perm.as_mut_ptr()) };
    // dbg!(perm);

    let xadj = &mut [0, 1, 2];
    let adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let vwgt = std::ptr::null_mut();
    let vsize = std::ptr::null_mut();
    let adjwgt = std::ptr::null_mut();
    let mut nparts: idx_t = 2;
    let tpwgts = std::ptr::null_mut();
    let ubvec = std::ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let _res = unsafe {
        pmetis::METIS_PartGraphRecursive(
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
}

#[test]
fn playground_can_utilize_lib() {
    let xadj = &mut [0, 1, 2];
    let adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let vwgt = std::ptr::null_mut();
    let vsize = std::ptr::null_mut();
    let adjwgt = std::ptr::null_mut();
    let mut nparts: idx_t = 2;
    let tpwgts = std::ptr::null_mut();
    let ubvec = std::ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        pmetis::METIS_PartGraphRecursive(
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
