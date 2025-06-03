//! Based on the Helix editor's 
#![allow(dead_code)]

use crate::util::{empty_arr, PartialArray};

type ParseResult<'a, Out> = Result<(&'a str, Out), &'a str>;

pub trait Parser<'a> {
    type Out;

    fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out>;
}

pub fn parse_str<'a, P>(p: P, s: &'a str) -> Option<P::Out> 
where P: Parser<'a>
{
    match p.parse(s) {
        Ok((rem, out)) if rem.trim_end().is_empty() => Some(out),
        _ => None,
    }
}

impl<'a, F, T> Parser<'a> for F
    where 
    F: Fn(&'a str) -> ParseResult<'a, T> 
{
    type Out = T;

    fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out> {
        self(input)
    }
}

impl<'a> Parser<'a> for str {
    type Out = &'a str;

    fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out> {
        if input.starts_with(self) {
            Ok((&input[..self.len()], &input[self.len()..]))
        } else {
            Err(input)
        }
    }
}


impl<'a> Parser<'a> for &str {
    type Out = &'a str;

    fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out> {
        if input.starts_with(self) {
            Ok((&input[self.len()..], &input[..self.len()]))
        } else {
            Err(input)
        }
    }
}

pub fn tok<'a>(input: &str) -> impl Parser<'a> {
    input
}

pub fn tok_ignore_case<'a>(tok: &str) -> impl Parser<'a> {
    move |input: &'a str| {
        match input.get(..tok.len()) {
            Some(val) if tok.eq_ignore_ascii_case(val) => {
                Ok((&input[tok.len()..], &input[..tok.len()]))
            },
            _ => Err(input),
        }
    }
}

pub fn single_char_fn<'a, F>(f: F) -> impl Parser<'a, Out = char> 
where
F: Fn(char) -> bool
{
    move |input: &'a str| {
        let Some(ch) = input.chars().next() else {
            return Err(input)
        };
        if f(ch) {
            Ok((&input[ch.len_utf8()..], ch))
        } else {
            return Err(input)
        }
    }
}

impl<'a> Parser<'a> for char {
    type Out = char;

    fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out> {
        single_char_fn(|ch| *self == ch).parse(input)
    }
}


pub fn or<'a, P1, P2, T>(p1: P1, p2: P2) -> impl Parser<'a, Out = T> 
where 
    P1: Parser<'a, Out = T>,
    P2: Parser<'a, Out = T>
{
    move |input| {
        p1.parse(input).or_else(|_| p2.parse(input)).or(Err(input))
    }
}

pub fn or_hetero<'a, P1, P2>(p1: P1, p2: P2) -> impl Parser<'a, Out = &'a str> 
where 
    P1: Parser<'a>,
    P2: Parser<'a>
{
    move |input| {
        capture(reff(&p1)).parse(input).or_else(|_| capture(reff(&p2)).parse(input)).or(Err(input))
    }
}

pub fn discard_or<'a, D, P>(discard: D, p: P) -> impl Parser<'a, Out = Option<P::Out>> 
where 
    D: Parser<'a>,
    P: Parser<'a>,
{
    move |input| {
        if let Ok((rem, _)) = discard.parse(input) {
            Ok((rem, None))
        } else {
            map(reff(&p), Some).parse(input)
        }
    }
}

pub fn maybe<'a, P, T>(p: P) -> impl Parser<'a, Out = Option<T>> 
where 
    P: Parser<'a, Out = T>,
{
    move |input| {
        match p.parse(input) {
            Ok((rem, out)) => Ok((rem, Some(out))),
            Err(_) => Ok((input, None)),
        }
    }
}

pub fn some<'a, P, T>(p: P) -> impl Parser<'a, Out = T> 
where 
    P: Parser<'a, Out = Option<T>>,
{
    move |input| {
        match p.parse(input) {
            Ok((rem, Some(out))) => Ok((rem, out)),
            _ => Err(input),
        }
    }
}

pub fn take_while<'a, F>(pat: F) -> impl Parser<'a, Out = &'a str>
where F: Fn(char) -> bool
{
    move |input: &'a str| {
        if let Some((i, _)) = input.char_indices().skip_while(|&(_i, c)| pat(c)).next() {
            Ok((&input[i..], &input[..i]))
        } else if !input.is_empty() {
            Ok((&input[input.len()..], input))
        } else {
            Err(input)
        }
    }
}

pub fn take_until<'a, F>(pat: F) -> impl Parser<'a, Out = &'a str>
where F: Fn(char) -> bool
{
    take_while(move |ch| !pat(ch))
}


/// never fails, 
pub fn trim_while<'a, F>(pat: F) -> impl Parser<'a, Out = &'a str>
where F: Fn(char) -> bool
{
    move |input: &'a str| {
        if let Some((i, _)) = input.char_indices().skip_while(|&(_i, c)| pat(c)).next() {
            Ok((&input[i..], &input[..i]))
        } else {
            Ok((&input[input.len()..], input))
        }
    }
}

pub fn everything<'a>() -> impl Parser<'a, Out = &'a str>
{
    trim_while(|_| true)
}

pub fn filter<'a, P, F>(p: P, filter: F) -> impl Parser<'a, Out = P::Out>
where 
    P: Parser<'a>,
    F: Fn(&P::Out) -> bool
{
    move |input: &'a str| {
        match p.parse(input) {
            Ok(ok) if filter(&ok.1) => Ok(ok),
            _ => Err(input),
        }
    }
}

pub fn sep<'a, P, S> (item: P, sep: S) -> impl Parser<'a, Out = Vec<P::Out>> 
where 
    P: Parser<'a>,
    S: Parser<'a>,
{
    move |input: &'a str| {
        let (mut input, first) = item.parse(input).or(Err(input))?;
        let mut out = vec![first];
        loop {
            let Ok((next, _)) = sep.parse(input) else {
                break
            };
            input = next;

            let Ok((next, res)) = item.parse(input) else {
                break
            };
            out.push(res);
            input = next;
        }

        Ok((input, out))
    }
}

pub fn whitespace_sep<'a, P> (item: P) -> impl Parser<'a, Out = Vec<P::Out>> 
where 
    P: Parser<'a>,
{
    sep(item, take_while(|c| c.is_whitespace()))
}

pub fn sep_array<'a, P, S, const N: usize>(item: P, sep: S) -> impl Parser<'a, Out = [P::Out; N]>
where 
    P: Parser<'a>,
    S: Parser<'a>,
{
    // TODO: don't allocate
    move |input| {
        if N == 0 {
            return Ok((input, empty_arr()))
        }
        let mut ret = PartialArray::new();
        let (mut input, first) = item.parse(input).or(Err(input))?;
        ret.push(first);
        for _ in 1..N {
            let Ok((next, _)) = sep.parse(input) else {
                break
            };
            input = next;

            let Ok((next, res)) = item.parse(input) else {
                break
            };
            ret.push(res);
            input = next;
        }

        if ret.is_complete() {
            Ok((input, ret.to_array()))
        } else {
            Err(input)
        }
    }
}

pub fn whitespace_sep_array<'a, P, const N: usize>(item: P) -> impl Parser<'a, Out = [P::Out; N]> 
where 
    P: Parser<'a>,
{
    sep_array(item, whitespace())
}

pub fn zero_plus<'a, P> (item: P) -> impl Parser<'a, Out = Vec<P::Out>> 
where 
    P: Parser<'a>,
{
    move |mut input: &'a str| {
        let mut out = vec![];
        loop {
            let Ok((next, res)) = item.parse(input) else {
                break
            };
            out.push(res);
            input = next;
        }

        Ok((input, out))
    }
}

/// like [`zero_plus`], but it just captures the string and doesn't allocatore
pub fn skip_zero_plus<'a, P> (item: P) -> impl Parser<'a, Out = &'a str> 
where 
    P: Parser<'a>,
{
    capture(move |mut input: &'a str| {
        loop {
            let Ok((next, _)) = item.parse(input) else {
                break
            };
            input = next;
        }
        Ok((input, ()))
    })
}

pub fn one_plus<'a, P> (item: P) -> impl Parser<'a, Out = Vec<P::Out>>
where 
    P: Parser<'a>,
{
    move |input: &'a str| {
        let (mut input, first) = item.parse(input).or(Err(input))?;
        let mut out = vec![first];
        loop {
            let Ok((next, res)) = item.parse(input) else {
                break
            };
            out.push(res);
            input = next;
        }

        Ok((input, out))
    }
}

/// like [`one_plus`], but it just captures the string and doesn't allocatore
pub fn skip_one_plus<'a, P> (item: P) -> impl Parser<'a, Out = &'a str> 
where 
    P: Parser<'a>,
{
    capture(move |input: &'a str| {
        let (mut input, _first) = item.parse(input).or(Err(input))?;
        loop {
            let Ok((next, _)) = item.parse(input) else {
                break
            };
            input = next;
        }
        Ok((input, ()))
    })
}

pub fn whitespace<'a>() -> impl Parser<'a, Out = ()> {
    map(take_while(|ch| ch.is_whitespace()), |_| ())
}

pub fn line_sep<'a, P> (item: P) -> impl Parser<'a, Out = Vec<P::Out>> 
where 
    P: Parser<'a>,
{
    sep(item, and(take_while(|c| c.is_whitespace() && c != '\n'), maybe("\n")))
}


fn integer_base<'a, I>() -> impl Parser<'a, Out = I> 
where I: std::str::FromStr
{
    move |input| {
        let (rem, int) = capture((maybe(or("+", "-")), skip_one_plus(single_char_fn(|ch| ch.is_ascii_digit())))).parse(input)?;
        let int = I::from_str(int).or(Err(input))?;
        Ok((rem, int))
    }
}

fn float_base<'a, I>() -> impl Parser<'a, Out = I> 
where I: std::str::FromStr
{
    move |input| {
        let sign = || maybe(or("+", "-"));
        let digit = || single_char_fn(|ch| ch.is_ascii_digit());
        // this is liberal about what it accepts, but will never go further than the grammar of
        // fXX::from_str. i.e. this will accept just `'.'`
        let (rem, int) = capture((
            sign(),
            or_hetero(
                (
                    skip_zero_plus(digit()),
                    ".",
                    skip_zero_plus(digit()),
                    maybe((
                        or("e", "E"),
                        sign(),
                        skip_one_plus(digit())
                    ))
                ),
                or(
                    tok_ignore_case("nan"),
                    or(
                        tok_ignore_case("inf"),
                        tok_ignore_case("infinity")
                    ),
                )
            ),
        )).parse(input)?;
        let int = I::from_str(int).or(Err(input))?;
        Ok((rem, int))
    }
}

macro_rules! from_str_impl {
    ($fn:ident => $($ty:ident)*) => {
        $(
        #[allow(dead_code)]
        #[doc = concat!("parse a single ", stringify!($ty))]
        pub fn $ty<'a>() -> impl Parser<'a, Out = $ty> {
            $fn()
        })*
    };
}

from_str_impl!(integer_base => u8 u16 u32 u64 usize i8 i16 i32 i64 isize);
from_str_impl!(float_base => f32 f64);


macro_rules! tuple_impl {
    ($($p:ident)*) => {
        impl<'a, $($p),*> Parser<'a> for ($($p,)*)
        where
            $($p: Parser<'a>,)*
        {
            type Out = ($($p::Out,)*);
            
            fn parse(&self, input: &'a str) -> ParseResult<'a, Self::Out> {
                #![allow(non_snake_case)]
                let original = input;
                let ($($p,)*) = self;
                $(
                let (input, $p) = $p.parse(input).or(Err(original))?;
                )*
                Ok((input, ($($p,)*)))
            }
        }
    };
}

tuple_impl!(T U);
tuple_impl!(T U V);
tuple_impl!(T U V W);
tuple_impl!(T U V W X);
tuple_impl!(T U V W X Y);
tuple_impl!(T U V W X Y Z);
tuple_impl!(T U V W X Y Z A);
tuple_impl!(T U V W X Y Z A B);
tuple_impl!(T U V W X Y Z A B C);
tuple_impl!(T U V W X Y Z A B C D);


pub fn reff<'a, P>(p: &P) -> impl Parser<'a, Out = P::Out> 
where
P: Parser<'a>
{
    move |input| p.parse(input)
}

pub fn and<'a, P1, P2>(p1: P1, p2: P2) -> impl Parser<'a, Out = (P1::Out, P2::Out)> 
where 
    P1: Parser<'a>,
    P2: Parser<'a>
{
    move |input| {
        (|i| p1.parse(i), |i| p2.parse(i)).parse(input)
    }
}

pub fn then<'a, P1, P2>(p1: P1, p2: P2) -> impl Parser<'a, Out = P2::Out> 
where 
    P1: Parser<'a>,
    P2: Parser<'a>
{
    move |input| {
        (reff(&p1), reff(&p2)).parse(input).map(|(rem, (_r1, r2))| (rem, r2))
    }
}

pub fn map<'a, P, F, T>(p: P, f: F) -> impl Parser<'a, Out = T> 
where
    P: Parser<'a>,
    F: Fn(P::Out) -> T
{
    move |input| p.parse(input).map(|(rem, x)| (rem, f(x)))
}

pub fn capture<'a, P>(p: P) -> impl Parser<'a, Out = &'a str>
where P: Parser<'a>
{
    move |input| {
        let (rem, _) = p.parse(input).or(Err(input))?;
        let len = input.len().saturating_sub(rem.len());
        Ok((rem, &input[..len]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_array() {
        let res: [i32; 4] = parse_str(whitespace_sep_array(i32()), "-12  4 0 121").unwrap();
        assert_eq!(res, [-12, 4, 0, 121]);

        let res: [i32; 3] = parse_str(whitespace_sep_array(i32()), "-12 40 121").unwrap();
        assert_eq!(res, [-12, 40, 121]);

        let res: [i32; 2] = parse_str(whitespace_sep_array(i32()), "-12 121").unwrap();
        assert_eq!(res, [-12, 121]);

        let res: [i32; 1] = parse_str(whitespace_sep_array(i32()), "121").unwrap();
        assert_eq!(res, [121]);

        let res: [i32; 0] = parse_str(whitespace_sep_array(i32()), "").unwrap();
        assert_eq!(res, []);
    }

    #[test]
    fn basic_rust() {
        fn trws<'a, P: Parser<'a>>(p: P) -> impl Parser<'a, Out = P::Out> {
            map((p, maybe(whitespace())), |(x, _)| x)
        }
        fn surround<'a, P: Parser<'a>>(l: char, r: char, p: P) -> impl Parser<'a, Out = P::Out> {
            map((trws(l), trws(p), r), |(_, x, _)| x)
        }
        fn paren<'a, P: Parser<'a>>(p: P) -> impl Parser<'a, Out = P::Out> {
            surround('(', ')', p)
        }
        let ws = whitespace;
        let sep_ch = |ch: char| !ch.is_ascii_alphanumeric() && ch != '_';
        let ident = || trws(capture((single_char_fn(|ch| ch.is_ascii_alphabetic()), take_until(sep_ch))));
        let ty = || trws(capture((ident(), maybe(surround('<', '>', ident())))));
        let arg = || trws((ident(), trws(':'), ty()));
        let fnitem = || (
            "fn", 
            ws(),
            ident(),
            trws(paren(sep(arg(), trws(',')))),
            maybe((trws("->"), trws(ty()))),
            surround('{', '}', maybe(ws()))
        );
        let input = "fn foo  (arg1: Vec<String>, arg2: usize) -> Vec<usize> { }";
        let res = parse_str(fnitem(), input).unwrap();
        assert_eq!(res, (
            "fn",
            (),
            "foo",
            vec![
                ("arg1", ':', "Vec<String>"),
                ("arg2", ':', "usize")
            ],
            Some(("->", "Vec<usize>")),
            Some(())
        ));
    }
}
