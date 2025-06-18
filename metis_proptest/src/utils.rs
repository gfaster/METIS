

/// zips two iterators, but makes the second item optional and has the length of the first
pub fn zip_first<I1, I2>(i1: I1, i2: I2) -> impl Iterator<Item = (I1::Item, Option<I2::Item>)> 
where
I1: IntoIterator,
I2: IntoIterator,
{
    let mut i2 = i2.into_iter();
    i1.into_iter().zip(std::iter::from_fn(move || Some(i2.next())))
}


pub trait RetainIndexed {
    type Item;
    fn retain_indexed<F>(&mut self, f: F) where F: FnMut(usize, &Self::Item) -> bool;
}

impl<T> RetainIndexed for Vec<T> {
    type Item = T;

    fn retain_indexed<F>(&mut self, mut f: F) where F: FnMut(usize, &Self::Item) -> bool {
        let mut i = 0;
        self.retain(|item| {
            let ret = f(i, item);
            i += 1;
            ret
        });
    }
}

/// construct a CSR from degrees specified in `degs[..(degs.len()-1)]`
pub fn make_csr(degs: &mut [usize]) {
    let mut i = 0;
    for d in degs {
        let inc = *d;
        *d = i;
        i = i.wrapping_add(inc);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_retain_indexed() {
        let keep = [true, false, true, true, false];
        let mut v: Vec<u8> = vec![0, 1, 2, 3, 4];
        v.retain_indexed(|i, _| keep[i]);
        assert_eq!(v, &[0, 2, 3][..])
    }

    #[test]
    fn test_make_csr() {
        let mut xadj = [2, 3, 1, 4, usize::MAX];
        make_csr(&mut xadj);
        assert_eq!(xadj, [0, 2, 5, 6, 10]);
    }
}
