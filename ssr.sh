#!/usr/bin/env bash

set -ex

rust-analyzer ssr '(iargmax_strd($a, $b, $c)) as usize ==>> iargmax_strd($a, $b, $c)'
rust-analyzer ssr 'iargmax_strd($n, $slice + $off, $stride) ==>> util::iargmax(&$slice[$off..], $stride)'

rust-analyzer ssr '(iargmax($a, $b, $c)) as usize ==>> iargmax($a, $b, $c)'
rust-analyzer ssr 'iargmax($n, $slice + $off, $stride) ==>> util::iargmax(&$slice[cntrng!($off, $n)], $stride)'
rust-analyzer ssr 'iargmax($n, $slice, $stride) ==>> util::iargmax(&$slice[..$n], $stride)'

rust-analyzer ssr '(iargmin($a, $b, $c)) as usize ==>> iargmin($a, $b, $c)'
rust-analyzer ssr 'iargmin($n, $slice + $off, $stride) ==>> util::iargmin(&$slice[cntrng!($off, $n)], $stride)'
rust-analyzer ssr 'iargmin($n, $slice, $stride) ==>> util::iargmin(&$slice[..$n], $stride)'

rust-analyzer ssr 'isum($n, $slice + $off, 1) ==>> $slice[cntrng!($off, $n)].iter().sum::<idx_t>()'
rust-analyzer ssr 'isum($n, $slice + $off, $stride) ==>> $slice[cntrng!($off, $n)].iter().step_by($stride).sum::<idx_t>()'
# rust-analyzer ssr 'isum($n, $slice, 1) ==>> $slice[..$n].iter().sum::<idx_t>()'
rust-analyzer ssr 'isum($n, $slice, $stride) ==>> $slice[..$n].iter().step_by($stride).sum::<idx_t>()'

rust-analyzer ssr 'iset($n, $val, $slice) ==>> $slice[..$n].fill($val)'

rust-analyzer ssr 'gk_min($a, $b) ==>> ($a).min($b)'
rust-analyzer ssr 'gk_max($a, $b) ==>> ($a).max($b)'
