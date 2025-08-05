use std::io::{self, BufRead};

use crate::{idx_t, invalid_format};

pub fn edges(mut r: impl BufRead) -> io::Result<Vec<(idx_t, idx_t)>> {
    let mut buf = String::with_capacity(80);

    let mut ret = Vec::new();

    loop {
        buf.clear();
        if r.read_line(&mut buf)? == 0 {
            break
        };
        let s = buf.trim_ascii_start();
        if s.starts_with(|c: char| c == '%' || c == '#') || s.is_empty() {
            continue
        }

        let mut it = s.split_whitespace();

        let from = it.next().ok_or_else(invalid_format)?.parse::<idx_t>().map_err(|_| invalid_format())?.checked_sub(1).expect("not 1-indexed");
        let to = it.next().ok_or_else(invalid_format)?.parse::<idx_t>().map_err(|_| invalid_format())?.checked_sub(1).expect("not 1-indexed");
        if to != from {
            ret.push((to, from));
            ret.push((from, to));
        }
    }


    Ok(ret)
}
