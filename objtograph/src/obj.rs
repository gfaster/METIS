use std::io::{BufRead, self};
use crate::idx_t;

/// get the first component of `s.split_once(delim)` or `s` if not found
fn split_first(s: &[u8], delim: u8) -> &[u8] {
    if let Some(pos) = s.iter().position(|&b| b == delim) {
        &s[..pos]
    } else {
        s
    }
}


/// graph of the edges
pub fn edges(mut r: impl BufRead) -> io::Result<Vec<(idx_t, idx_t)>> {
    let mut buf = Vec::with_capacity(80);
    let mut ret = Vec::new();
    while r.read_until(b'\n', &mut buf)? != 0 {
        buf.pop(); // read_until includes pattern
        let Some(edge) = edge_line(&buf) else { continue };
        for (v1, v2) in edge {
            if v1 == v2 {
                continue
            }
            ret.push((v1, v2));
            ret.push((v2, v1));
        }
        buf.clear();
    }

    Ok(ret)
}

fn edge_line(s: &[u8]) -> Option<impl Iterator<Item = (idx_t, idx_t)> + use<'_>> {
    let s = split_first(s, b'#'); // strip comment

    let s = s.strip_prefix(b"f ")?; // ignore anything that isn't a face

    let mut it = s.split(|&b| b.is_ascii_whitespace()).filter(|b| !b.is_empty()).flat_map(parse_face_triplet);
    let first = it.next()?;
    let mut prev = first;
    let mut fin = false;

    Some(std::iter::from_fn(move || {
        if let Some(next) = it.next() {
            let ret = (prev, next);
            prev = next;
            Some(ret)
        } else if !fin {
            // emit looped iteration
            fin = true;
            Some((prev, first))
        } else {
            None
        }
    }))
}

fn parse_face_triplet(s: &[u8]) -> Option<idx_t> {
    debug_assert_eq!(s.trim_ascii(), s);
    let s = split_first(s, b'/');

    if s.is_empty() {
        return None
    }
    let s = std::str::from_utf8(s).ok()?;
    let s: idx_t = s.parse().ok()?;
    s.checked_sub(1)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_face_triplets() {
        assert_eq!(parse_face_triplet(b"2/3/4"), Some(1));
        assert_eq!(parse_face_triplet(b"2//4"), Some(1));
        assert_eq!(parse_face_triplet(b"2//"), Some(1));
        assert_eq!(parse_face_triplet(b"2"), Some(1));
        assert_eq!(parse_face_triplet(b"1"), Some(0));
        assert_eq!(parse_face_triplet(b"0"), None);
        assert_eq!(parse_face_triplet(b""), None);
        assert_eq!(parse_face_triplet(b"-1"), None);
    }

    #[test]
    fn test_edge_line() {
        let s = b"f 1 2 3 4 # first face\r\n";
        let es: Vec<_> = edge_line(s).unwrap().collect();
        assert_eq!(es, [(0, 1), (1, 2), (2, 3), (3, 0)]);
    }
}
