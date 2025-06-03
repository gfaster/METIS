//! Parsing of the Matrix Market format
//!
//! See: https://networkrepository.com/format-info.php
//!
//! See: https://math.nist.gov/MatrixMarket/formats.html#MMformat

use std::io::{self, BufRead};

use crate::{idx_t, invalid_format};
use crate::parser as p;
use p::parse_str;

enum Symmetry {
    Symm,
    General
}

enum Field {
    Int,
    Real,
    Pattern,
}

enum Format {
    Coord,
    Array
}

pub fn edges(mut r: impl BufRead) -> io::Result<Vec<(idx_t, idx_t)>> {
    let mut buf = String::with_capacity(80);

    // read header
    if r.read_line(&mut buf)? == 0 {
        return Err(invalid_format())
    };
    let Ok((format, field, symm)) = parse_header(&buf) else {
        return Err(invalid_format())
    };

    let ret = match format {
        Format::Coord => parse_coordinate(r, symm),
        Format::Array => parse_array(r, field),
    }?;


    Ok(ret)
}

fn parse_array(mut r: impl BufRead, field: Field) -> io::Result<Vec<(idx_t, idx_t)>> {
    let mut buf = String::with_capacity(80);

    // read size line
    if r.read_line(&mut buf)? == 0 {
        return Err(invalid_format())
    };
    let Some([m, n]) = parse_str(p::whitespace_sep_array(p::usize()), &buf) else {
        return Err(invalid_format())
    };

    fn process(mut r: impl BufRead, m: usize, n: usize, f: impl Fn(&str) -> bool) 
        -> io::Result<Vec<(idx_t, idx_t)>>
    {
        // times 2 since it's two sided - we merge later
        let mut ret = Vec::with_capacity(m * n * 2);
        let mut buf = String::with_capacity(80);
        let mut idx = 0; 
        while ({buf.clear(); r.read_line(&mut buf)?}) != 0 {
            let s = buf.trim();
            if s.starts_with('%') {
                continue
            }
            if f(&s) {
                let r = (idx % m) as idx_t;
                let c = (idx / m) as idx_t;
                ret.push((r, c));
                ret.push((c, r));
            }
            idx += 1;
        }
        Ok(ret)
    }

    match field {
        Field::Int => {
            process(r, m, n, |val| val != "0")
        },
        Field::Real => {
            process(r, m, n, |val| {
                let Ok(r) = val.parse::<f64>() else {
                    return false;
                };
                r < 0.00001
            })
        },
        Field::Pattern => Err(invalid_format()),
    }
}

fn parse_coordinate(mut r: impl BufRead, symm: Symmetry) -> io::Result<Vec<(idx_t, idx_t)>> {
    let mut buf = String::with_capacity(80);

    // read size line
    let [_n, _m, nonzeros] = {
        while r.read_line(&mut buf)? != 0 {
            if !buf.starts_with('%') { break }
            buf.clear()
        }
        parse_str(p::whitespace_sep_array(p::usize()), &buf).ok_or_else(invalid_format)?
    };

    // times 2 since it's two sided -- in the general symmetry case, there could be a symmetric
    // pattern but we won't find out until later
    let cap = match symm {
        Symmetry::Symm => nonzeros * 2,
        Symmetry::General => nonzeros * 4,
    };
    let mut ret = Vec::with_capacity(cap);
    let mut buf = String::with_capacity(80);
    while ({buf.clear(); r.read_line(&mut buf)?}) != 0 {
        let s = buf.trim_start();
        if s.starts_with('%') || s.is_empty() {
            continue
        }
        let ([v0, v1], _) = parse_str((p::whitespace_sep_array(idx_t()), p::everything()), s).ok_or_else(invalid_format)?;
        let v0 = v0.checked_sub(1).expect("graph is zero indexed!");
        let v1 = v1.checked_sub(1).expect("graph is zero indexed!");

        if v0 == v1 {
            continue
        }

        ret.push((v0, v1));
        ret.push((v1, v0));
    }
    Ok(ret)

}

fn parse_header(h: &str) -> Result<(Format, Field, Symmetry), ()> {
    let mut it = h.split_whitespace();
    let mut fields = [None; 5];
    fields.iter_mut().for_each(|f| *f = it.next());
    let fields = fields;
    let format = match fields[2] {
        Some("coordinate") => Format::Coord,
        Some("array") => Format::Array,
        _ => return Err(())
    };
    let field = match fields[3] {
        Some("real") | Some("double") => Field::Real,
        Some("integer") => Field::Int,
        Some("pattern") => Field::Pattern,
        Some("complex") => {
            eprintln!("Complex field are unsupported");
            return Err(())
        },
        _ => return Err(())
    };
    let symm = match fields[4] {
        Some("symmetric") => Symmetry::Symm,
        Some("general") => Symmetry::General,
        _ => return Err(())
    };
    let valid = matches!(fields, [
        Some("%MatrixMarket"|"%%MatrixMarket"),
        Some("matrix"),
        _,
        Some("real" | "double" | "complex" | "integer" | "pattern"),
        _,
    ]);
    valid.then_some((format, field, symm)).ok_or(())
}
