# METIS 

METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, 
and producing fill reducing orderings for sparse matrices. The algorithms implemented in 
METIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
partitioning schemes developed in our lab.

This repository is my attempt of porting METIS to Rust file-by-file, line-by-line.

As of 2023-09-26, 3,869 of 15,645 lines of C code have been ported
(not including header files, using `wc`) - that means 11,776 to go. It has
taken 8,376 lines of Rust code to achieve this.

Check [`translation.md`](./translation.md) for my porting process.

Check [`appendix.md`](./appendix.md) for some info on common functions used.

## Goals

My immediate goal is to get METIS ported to idiomatic Rust. This is a highly
mechanical process that is fairly uninteresting. Once that is done, I want to
try to add support for partitioning euclidian meshes.

##  Downloading METIS

You can download METIS by simply cloning it using the command:
```
git clone https://github.com/gfaster/METIS.git
```

## Building standalone METIS binaries and library

To build METIS you can follow the instructions below:

### Dependencies

General dependencies for building this METIS refactor are: gcc, Rust, and Cargo. 
In Ubuntu systems these can be obtained from the apt package manager (e.g., apt-get install make, etc) 

```
sudo apt-get install build-essential
```

Rust and Cargo should be installed via [`rustup`](https://rustup.rs/)

### Building

1. run `cargo build` and/or `cargo test`


## Copyright & License Notice
Copyright 1998-2020, Regents of the University of Minnesota

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

