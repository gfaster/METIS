#![crate_type = "staticlib"]
#![allow(non_snake_case, non_camel_case_types)]

type idx_t = i32;

#[repr(C)]
pub struct ikv_t {
    pub key: idx_t,
    pub val: idx_t,
}


// need sorting functions since Rust's functions re-order identical keys differently

#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsorti(n: usize, ikv: *mut ikv_t) {
    let ikv = std::slice::from_raw_parts_mut(ikv, n);
    ikv.sort_unstable_by_key(|kv| kv.key);
}

#[no_mangle]
pub unsafe extern "C" fn libmetis__ikvsortd(n: usize, ikv: *mut ikv_t) {
    let ikv = std::slice::from_raw_parts_mut(ikv, n);
    ikv.sort_unstable_by_key(|kv| std::cmp::Reverse(kv.key));
}
