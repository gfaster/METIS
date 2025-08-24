//! Rust one-to-one (functionality-wise) implementation of GKlib's getopt. The goal is to be almost
//! indistinguishable to the user whether we're using this or GKlib's getopt. For the most part, it
//! seems that GKlib's getopt is lifted from (probably an older) glibc, so I will black-box
//! re-implement the POSIX behavior that is used by the CLI programs.
//!
//! I will add some minor improvements for error reporting though, as I found the original to be
//! quite lacking
//!
//! In the future I want to have a feature switch to use this or clap

use std::collections::VecDeque;

use metis::idx_t;

#[derive(Debug, Clone, Copy)]
pub enum ArgPresence {
    NoArg,
    Required,
}

impl ArgPresence {
    /// Returns `true` if the arg presence is [`Required`].
    ///
    /// [`Required`]: ArgPresence::Required
    #[must_use]
    pub fn is_required(&self) -> bool {
        matches!(self, Self::Required)
    }
}

#[derive(Debug)]
pub struct Opt {
    pub name: &'static str,
    pub has_arg: ArgPresence,
    pub idx: usize,
    pub table: Option<&'static [OptEnumEntry]>,
}

impl Opt {
    pub fn parse_enum<'a, 'o>(
        &'o self,
        variant: &'a str,
    ) -> Result<idx_t, Box<ParsedOptionError<'a, 'o>>> {
        let table = self.table.expect("not an enum option");
        if let Some(e) = table.iter().find(|e| e.name == variant) {
            Ok(e.val)
        } else {
            Err(Box::new(ParsedOptionError::InvalidEnumVariant(
                variant, self,
            )))
        }
    }
}

#[derive(Debug)]
pub struct OptEnumEntry {
    pub name: &'static str,
    pub val: idx_t,
}

macro_rules! opt_enum {
    ($($name:literal = $val:expr),* $(,)?) => {
        &[
            $( getopt::OptEnumEntry { name: $name, val: $val as idx_t } ),*
        ]
    };
}
pub(crate) use opt_enum;

macro_rules! opt {
    ($(($($tt:tt)*)),*$(,)?) => {
        &[$(getopt::opt!(@one $($tt)*)),*]
    };
    (@one $s:literal, $has:ident, $idx:expr $(,)?) => {
        getopt::Opt {
            name: $s,
            has_arg: getopt::ArgPresence::$has,
            idx: $idx as usize,
            table: None,
        }
    };
    (@one $s:literal, $has:ident, $idx:expr, {$($tt:tt)*} $(,)?) => {
        getopt::Opt {
            name: $s,
            has_arg: getopt::ArgPresence::$has,
            idx: $idx as usize,
            table: Some(getopt::opt!(@table $($tt)*)),
        }
    };
    (@table $($name:literal = $val:expr),* $(,)?) => {
        getopt::opt_enum!($($name = $val),*)
    };
}
pub(crate) use opt;

fn valid_optstr(s: &str) -> bool {
    !s.is_empty()
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"-_".contains(&b))
        && !s.starts_with(['-', '_'])
}

pub enum ParsedOption<'a, 'o> {
    Opt { val: Option<&'a str>, opt: &'o Opt },
    Nonopt(&'a str),
}

#[derive(Debug, Clone)]
pub enum ParsedOptionError<'a, 'o> {
    Help,
    MissingArg(&'o Opt),
    SuperfluousArg(&'o Opt),
    InvalidName,
    TooManyArgs,
    NotEnoughArgs,
    AmbiguousAbbreviation(&'a str),
    InvalidEnumVariant(&'a str, &'o Opt),
    IntParse(std::num::ParseIntError),
    FloatParse(std::num::ParseFloatError),
}

impl<'a, 'o> ParsedOptionError<'a, 'o> {
    /// Returns `true` if the parsed option error is [`Help`].
    ///
    /// [`Help`]: ParsedOptionError::Help
    #[must_use]
    pub fn is_help(&self) -> bool {
        matches!(self, Self::Help)
    }
}

impl<'a, 'o> From<std::num::ParseFloatError> for Box<ParsedOptionError<'a, 'o>> {
    fn from(v: std::num::ParseFloatError) -> Self {
        Box::new(ParsedOptionError::FloatParse(v))
    }
}

impl<'a, 'o> From<std::num::ParseIntError> for Box<ParsedOptionError<'a, 'o>> {
    fn from(v: std::num::ParseIntError) -> Self {
        Box::new(ParsedOptionError::IntParse(v))
    }
}

impl<'a, 'o> From<std::num::ParseFloatError> for ParsedOptionError<'a, 'o> {
    fn from(v: std::num::ParseFloatError) -> Self {
        Self::FloatParse(v)
    }
}

impl<'a, 'o> From<std::num::ParseIntError> for ParsedOptionError<'a, 'o> {
    fn from(v: std::num::ParseIntError) -> Self {
        Self::IntParse(v)
    }
}

pub type ParseResult<'a, 'o> = Result<ParsedOption<'a, 'o>, Box<ParsedOptionError<'a, 'o>>>;

// TODO: switch to OsString
/// Parses command line options. Doesn't validate enums or parse numbers. Can be called with
/// `&std::env::args()`
pub fn parse_options<'a, 'o>(
    args: &'a [String],
    opts: &'o [Opt],
) -> impl Iterator<Item = ParseResult<'a, 'o>> {
    if cfg!(debug_assertions) {
        for opt in opts {
            assert!(valid_optstr(opt.name));
            for entry in opt.table.into_iter().flatten() {
                assert!(valid_optstr(entry.name))
            }
        }

        // make sure there's no duplicates
        let len = opts.len();
        let mut opts: Vec<_> = opts.iter().map(|o| o.name).collect();
        opts.sort_unstable();
        opts.dedup();
        assert_eq!(opts.len(), len, "duplicate arguments");
    }

    // don't panic if missing argv[0]
    let mut args = args.get(1..).unwrap_or(args);

    let mut error = false;

    let mut nonopts: VecDeque<&str> = VecDeque::with_capacity(args.len());

    macro_rules! err {
        ($error:expr) => {{
            error = true;
            return Some(Err(Box::new($error)));
        }};
    }
    std::iter::from_fn(move || {
        loop {
            if error {
                return None;
            }
            let Some(next) = args.split_off_first() else {
                // no more options, return the non-option arguments
                let nonopt = nonopts.pop_front()?;
                return Some(Ok(ParsedOption::Nonopt(nonopt)));
            };
            let next: &str = next;

            // '--' means all remaining arguments are non-options
            if next == "--" {
                nonopts.extend(args.iter().map(std::ops::Deref::deref));
                args = &[];
                continue;
            }

            let Some(next) = next.strip_prefix('-') else {
                // non-option argument, push it to nonopts
                nonopts.push_back(next);
                continue;
            };

            // remove possible remaining leading '-'
            let next = next.strip_prefix('-').unwrap_or(next);

            // if the argument is in the form --name=val, we already have the value in this
            // argument
            let (name, val) = next
                .split_once('=')
                .map_or((next, None), |(n, v)| (n, Some(v)));

            // option argument, find the option it corresponds to
            let mut abbrev = None;
            let mut ambiguous = false;
            let mut exact = None;
            for opt in opts {
                if opt.name.starts_with(name) {
                    // if lengths match, then we have exact match
                    if opt.name.len() == name.len() {
                        debug_assert_eq!(opt.name, name);
                        exact = Some(opt);
                        break;
                    }
                    ambiguous |= abbrev.is_some();
                    abbrev = Some(opt)
                }
            }

            let opt = match (exact, abbrev) {
                (Some(exact), _) => {
                    // exact match
                    exact
                }
                (None, None) => {
                    // found nothing
                    err!(ParsedOptionError::InvalidName);
                }
                (None, Some(_)) if ambiguous => {
                    // no exact match but multiple possible abbreviations
                    err!(ParsedOptionError::AmbiguousAbbreviation(name));
                }
                (None, Some(abbrev)) => {
                    // no exact match but only one abbreviation
                    abbrev
                }
            };

            // process --name=val arg, if we had it
            match opt.has_arg {
                ArgPresence::NoArg if val.is_some() => {
                    // had value assigned in --name=val form, but this option does not accept
                    // a value
                    err!(ParsedOptionError::SuperfluousArg(opt));
                }
                ArgPresence::NoArg => return Some(Ok(ParsedOption::Opt { val: None, opt })),
                ArgPresence::Required => {
                    if let Some(val) = val {
                        return Some(Ok(ParsedOption::Opt {
                            val: Some(val),
                            opt,
                        }));
                    }
                }
            }

            // next value should be an argument
            let Some(val) = args.split_off_first() else {
                err!(ParsedOptionError::MissingArg(opt))
            };

            return Some(Ok(ParsedOption::Opt {
                val: Some(val),
                opt,
            }));
        }
    })
}
