//! Replacement routines for gklib that ought to never use the original

#[no_mangle]
pub extern "C" fn gk_randinit(seed: u64) {
    fastrand::seed(seed);
}


#[no_mangle]
pub extern "C" fn gk_randint64() -> u64 {
    // must be positive, or you get horrible segfaults!
    fastrand::u64(..) & (u64::MAX >> 1)
}

#[no_mangle]
pub extern "C" fn gk_randint32() -> u32 {
    // must be positive, or you get horrible segfaults!
    fastrand::u32(..) & (u32::MAX >> 1)
}
