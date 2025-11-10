# METIS 

METIS is a set of serial programs for partitioning graphs, partitioning finite
element meshes, and producing fill reducing orderings for sparse matrices. The
algorithms implemented in METIS are based on the multilevel
recursive-bisection, multilevel k-way, and multi-constraint partitioning
schemes developed in the Karypis Lab at UMN.

This repository is my attempt of porting METIS to Rust file-by-file, line-by-line.

As of 2025-11-10, nearly 100% of the C code has been converted to Rust, with
just a little bit more testing we'll be able to make it usable as a drop-in
replacement.

Check [`translation.md`](./translation.md) for my porting process (mostly
historical).

Check [`appendix.md`](./appendix.md) for some info on common functions used
(mostly historical).

## Disclaimer

This project is in a perpetually broken state until everything is ported. In
order to make porting easier, we flaunt many rules of Rust. This means that **we
cause tons of Undefined Behavior**. Once fully ported, I expect METIS to
require very little, if any, uses of `unsafe`. However, for now, expect
enabling optimizations to Cause Problems.

To mitigate *some* of this harm, we use unstable rustc flags to disable
`noalias` annotations on pointers.

Furthermore, we rely on a GNU toolchain and unportable (even unspecified)
behaviors.

This should all change *eventually*.

## Roadmap

This is a rough roadmap of the major milestones

- [ ] Literal C port (~90%). Nearly every function is literally interchangeable.
    - Part of the testing process is literally swapping out individual
    functions
- [ ] Break internal ABI to allow more use of Rust types and less flagrant
Undefined Behavior
    - I expect this to become mostly safe
    - Still largely going to be unidiomatic Rust
    - Still going to have identical observable behavior to upstream
    - Likely first nonzero version
- [ ] Optimize into fast and idiomatic Rust
    - There are a bunch of idioms that make sense in C but not in Rust, let's
    change those
    - To make porting easier, I removed a number of optimizations, most notably
    around arena allocation. I may try to reintroduce that here.
    - I would not expect the Rust version to be faster than the C version until
    this step nears completion.
- [ ] Add a decent Rust API and prepare packaging
    - Designing a usable Rust API has been a challenge that's been on my mind
    throughout the porting process, but I'm still unsure of how to make it
    idiomatic and match the overhead of the original API.
    - I'll also have to figure out how to handle Metis's hard-coded `idx_t`
    type


## Building

Building requires a recent-ish GNU/Linux system, with gcc, glibc, cargo, and a
recent rustc nightly version. With that, `cargo b` should be sufficient to
build. Note that during the porting process, assumptions are made about the
directory structure of the output, so installation most likely will not work.

It is recommended that you build under the Nix shell.

Again, *this software is not yet suitable for production*. Do not use it
because it will break and you will be sad.

## Copyright & License Notice

Copyright 1998-2020, Regents of the University of Minnesota
Copyright 2023-2025, Gavin Rohrer

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

