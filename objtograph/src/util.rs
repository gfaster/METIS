use std::mem::MaybeUninit;

pub const fn empty_arr<T, const N: usize>() -> [T; N] {
    // don't want to const assert this
    assert!(N == 0);
    assert!(size_of::<[T; N]>() == 0);
    let ret: [T; 0] = [];
    unsafe { std::mem::transmute_copy(&ret) }
}

pub struct PartialArray<T, const N: usize> {
    arr: [MaybeUninit<T>; N],
    len: usize,
}

impl<T, const N: usize> PartialArray<T, N> {
    pub const fn new() -> Self {
        Self {
            arr: [const { MaybeUninit::uninit() }; N],
            len: 0
        }
    }

    pub fn push(&mut self, item: T) {
        assert!(self.len < self.arr.len());
        self.arr[self.len].write(item);
        self.len += 1;
    }

    pub fn to_array(self) -> [T; N] {
        assert_eq!(self.len, self.arr.len());
        const { assert!(size_of::<[T; N]>() == size_of::<[MaybeUninit<T>; N]>()); }
        unsafe { std::mem::transmute_copy(&self.arr) }
    }

    pub fn is_complete(&self) -> bool {
        self.len == self.arr.len()
    }
}

impl<T, const N: usize> Drop for PartialArray<T, N> {
    fn drop(&mut self) {
        if std::mem::needs_drop::<T>() {
            // use Vec implementation
            unsafe {
                std::ptr::drop_in_place(std::ptr::slice_from_raw_parts_mut(self.arr.as_mut_ptr(), self.len));
            }
        }
    }
}
