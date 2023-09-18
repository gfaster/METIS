# METIS 

METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, 
and producing fill reducing orderings for sparse matrices. The algorithms implemented in 
METIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
partitioning schemes developed in our lab.

This repository is my attempt of porting METIS to Rust file-by-file, line-by-line.

As of 2023-09-17 (`8cde126`), 2,056 of 15,642 lines of C code have been ported
(not including header files, using `wc`) - that means 13,586 to go. It has
taken 4,900 lines of Rust code to achieve this.

Check [`translation.md`](./translation.md) for my porting process.

Check [`appendix.md`](./appendix.md) for some info on common functions used.

##  Downloading METIS

You can download METIS by simply cloning it using the command:
```
git clone https://github.com/gfaster/METIS.git
```

## Building standalone METIS binaries and library

To build METIS you can follow the instructions below:

### Dependencies

General dependencies for building this METIS refactor are: gcc, GNUmake, cmake, build-essential, Rust, and Cargo. 
In Ubuntu systems these can be obtained from the apt package manager (e.g., apt-get install make, etc) 

```
sudo apt-get install build-essential
sudo apt-get install make
sudo apt-get install cmake
```

Rust and Cargo should be installed via [`rustup`](https://rustup.rs/)

### Building

3. `cd` to `METIS/GKlib` and run `make config` and then `make install`
4. return to the parent directory (`METIS`)
5. run `cargo build` and/or `cargo test`

full commands, assuming `gcc`, `cmake`, `build-essential`, `rustc`, and `cargo` are already installed:
```sh
git clone https://github.com/gfaster/METIS
cd METIS

# should be changed properly, but just in case
cd GKlib && make config prefix=../install 

make install
cd ..
cargo test
```


## Copyright & License Notice
Copyright 1998-2020, Regents of the University of Minnesota

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

