#![allow(dead_code)]
use std::{any::Any, io, ptr::NonNull};

use crate::{idx_t, real_t};

type Error = Box<dyn std::error::Error>;
type Result = std::result::Result<(), Error>;

pub trait Scanf<'a> {
    fn scanstr(&mut self, s: &str) -> Result;
    fn scan<R: ReadItem>(&mut self) -> std::result::Result<R, Error>;
}

impl<'a> Scanf<'a> for io::Cursor<&'a [u8]> {
    fn scanstr(&mut self, s: &str) -> Result {
        read_str(self, s)
    }

    fn scan<R: ReadItem>(&mut self) -> std::result::Result<R, Error> {
        R::read_item(self)
    }
}

fn take_while_limited<'a>(
    cur: &mut io::Cursor<&[u8]>,
    buf: &'a mut [u8],
    mut f: impl FnMut(u8) -> bool,
) -> &'a mut [u8] {
    let &backing = cur.get_ref();
    let mut len = 0;
    let start = cur.position() as usize;
    for (&src, dst) in backing[start..].iter().zip(&mut *buf) {
        if !f(src) {
            break;
        }
        *dst = src;
        len += 1
    }
    cur.set_position((start + len) as u64);
    &mut buf[..len]
}

type ReadItemFn = unsafe fn(&mut io::Cursor<&[u8]>, Option<NonNull<u8>>) -> Result;

pub trait ReadItem {
    fn read_item(cur: &mut io::Cursor<&[u8]>) -> std::result::Result<Self, Error>
    where
        Self: Sized;
}

impl ReadItem for idx_t {
    fn read_item(cur: &mut io::Cursor<&[u8]>) -> std::result::Result<Self, Error> {
        let mut buf = [0_u8; Self::MAX.ilog10() as usize + 3];
        let mut first = true;
        let &mut ref buf = take_while_limited(cur, &mut buf, |b| {
            if first {
                if b == b'-' {
                    return true;
                }
                first = false;
            }
            b.is_ascii_digit()
        });
        // Safety: buf contains only '-' and ascii digits
        let s = unsafe { std::str::from_utf8_unchecked(buf) };
        Self::from_str_radix(s, 10).map_err(|e| e.into())
    }
}

impl ReadItem for real_t {
    fn read_item(cur: &mut io::Cursor<&[u8]>) -> std::result::Result<Self, Error> {
        let mut buf = [0_u8; 64];
        let mut first = true;
        let mut found_decimal = false;
        let &mut ref buf = take_while_limited(cur, &mut buf, |b| {
            if first {
                first = false;
                if b == b'-' {
                    return true;
                }
            }
            if !found_decimal && b == b'.' {
                found_decimal = false;
                return true;
            }
            b.is_ascii_digit()
        });
        // Safety: buf contains only '-' and ascii digits
        let s = unsafe { std::str::from_utf8_unchecked(buf) };
        Ok(s.parse()?)
    }
}

pub enum Arg {
    F(ReadItemFn),
    Str(&'static str),
}

pub struct Arguments<'a> {
    args: &'a [Arg],
}

impl Arguments<'_> {
    unsafe fn read(self, source: &mut io::Cursor<&[u8]>, dests: &[Option<NonNull<u8>>]) -> Result {
        let mut di = 0;
        for arg in self.args {
            match arg {
                Arg::F(f) => {
                    debug_assert!(dests[di].is_some());
                    let ret = f(&mut *source, dests[di]);
                    di += 1;
                    ret
                }
                Arg::Str(s) => {
                    debug_assert!(dests[di].is_none());
                    di += 1;
                    read_str(&mut *source, s)
                }
            }?
        }
        Ok(())
    }
}

fn read_ws(source: &mut io::Cursor<&[u8]>) -> Result {
    let &backing = source.get_ref();
    let mut bi = source.position() as usize;
    let rem = backing[bi..].trim_ascii_start();
    let diff = backing[bi..].len() - rem.len();
    if diff == 0 {
        Err("No whitespace")?;
    }
    bi += diff;
    source.set_position(bi as u64);
    Ok(())
}

fn read_str(source: &mut io::Cursor<&[u8]>, s: &str) -> Result {
    if s.chars().all(|ch| ch.is_ascii_whitespace()) {
        return read_ws(source)
    }
    let mut bi = source.position() as usize;
    // debug_assert!(bi < source.get_ref().len());
    let &backing = source.get_ref();
    for s in s.bytes() {
        if s.is_ascii_whitespace() {
            bi += backing[bi..]
                .iter()
                .take_while(|&&b| b.is_ascii_whitespace())
                .count();
        } else {
            if Some(&s) != backing.get(bi) {
                Err("string literal did not match")?;
            }
            bi += 1;
        }
    }
    todo!()
}
