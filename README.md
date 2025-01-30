# METIS 

METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, 
and producing fill reducing orderings for sparse matrices. The algorithms implemented in 
METIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
partitioning schemes developed in our lab.

This repository is my attempt of porting METIS to Rust file-by-file, line-by-line.

As of 2024-01-20, 5,210 of 15,654 lines of C code have been ported
(not including header files, using `wc`) - that means 10,444 to go. It has
taken 11,397 lines of Rust code to achieve this.

Check [`translation.md`](./translation.md) for my porting process.

Check [`appendix.md`](./appendix.md) for some info on common functions used.

## Disclaimer

This project is in a perpetually broken state until everything is ported. In
order to make porting easier, we flaunt many rules of Rust. This means that **we
cause tons of Undefined Behavior**. Once fully ported, I expect METIS to
require very little, if any, uses of `unsafe`. However, for now, expect
enabling optimizations to Cause Problems.

To mitigate *some* of this harm, build with:

```bash
RUSTFLAGS="-Zbox-noalias=false -Zmutable-noalias=false" cargo +nightly b
```

Furthermore, we rely on a GNU toolchain and unportable (even unspecified)
behaviors.

This should all change *eventually*.

## Goals

My immediate goal is to get METIS ported to idiomatic Rust. This is a highly
mechanical process that is fairly uninteresting. Once that is done, I want to
try to add support for partitioning euclidian meshes.


## Building standalone METIS binaries and library

To build METIS you can follow the instructions below:

### Dependencies

General dependencies for building this METIS refactor are: gcc, Rust, and Cargo. 
In Ubuntu systems these can be obtained from the apt package manager (e.g., apt-get install make, etc) 

```
sudo apt-get install build-essential
```

Rust and Cargo should be installed via [`rustup`](https://rustup.rs/).

Note that during the porting process, we rely on a lot of horrible hacks, which
means, among other things, that the final build relies on the directory
structure of `/target`.

### Building

1. run `cargo build` and/or `cargo test`

## Testing ported functions


## Copyright & License Notice
Copyright 1998-2020, Regents of the University of Minnesota

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

