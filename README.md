# METIS 

METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, 
and producing fill reducing orderings for sparse matrices. The algorithms implemented in 
METIS are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
partitioning schemes developed in our lab.

##  Downloading METIS

You can download METIS by simply cloning it using the command:
```
git clone https://github.com/KarypisLab/METIS.git
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

In addition, you need to clone [GKlib](https://github.com/KarypisLab/GKlib) as a subdirectory of METIS.

### Building

1. Clone Gklib
    - Should be as a subdirectory in this repo: `METIS/GKlib`
2. in the newly cloned `METIS/GKlib/Makefile`, set `prefix` to `../install`
3. `cd` to `METIS/GKlib` and run `make config` and then `make install`
4. return to the parent directory (`METIS`)
5. run `make lib` and then `make test`

full commands, assuming `gcc`, `cmake`, `build-essential`, `rustc`, and `cargo` are already installed:
```sh
git clone https://github.com/gfaster/METIS
cd METIS
git clone https://github.com/KarypisLab/GKlib
cd GKlib

# should be changed properly, but this still works
make config prefix=../install 

make install
cd ..
make lib
make test
```

### Advanced debugging related options:

    gdb=1           - Build with support for GDB [off by default]
    debug=1         - Enable debugging support [off by default]
    assert=1        - Enable asserts [off by default]
    assert2=1       - Enable very expensive asserts [off by default]

### Other make commands

    make uninstall
         Removes all files installed by 'make install'.

    make clean
         Removes all object files but retains the configuration options.

    make distclean
         Performs clean and completely removes the build directory.


## Copyright & License Notice
Copyright 1998-2020, Regents of the University of Minnesota

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

