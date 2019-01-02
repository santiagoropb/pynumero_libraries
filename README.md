PyNumero libraries
==================

This repository contains a cmake project for building shared
libraries for the python numerical optimization package PyNumero.

The current version of PyNumero relies on two share libraries for
efficient numerical calculations while solving optimization problems

1. libpynumero_ASL.so
2. libpynumero_HSL.so

The first library, libpynumero_ASL, provides functionality for fast automatic
differentiation (AD) subroutines required for computing first- and second-
oder derivatives necessary for gradient-based optimization algorithms. The
subroutines for AD are the same ones used by AMPL nonlinear optimization
solvers since libpynumero_ASL links with the ampl solver library ASL.

The second library, libptnumero_HSL, provides functionality to interface
with the MA27 linear solver that is available in the Harwell subroutine
library. The linear solvers in HSL have proofed to be excellent direct
solver for finding solutions of saddle point problems.

A conda-forge recipe has been created to facilitate the distribution of
libpynumero_ASL without users requiring to compile any code. To get
libpynumero_ASL library via conda-forge users may run:

```
conda install -c conda-forge pynumero_libraries
```

If conda is not available, compilation of libpynumero_ASL library can be
achieve with the following steps:
```
cd /third_party/ASL
./get.ASL.sh
cd solvers
./configurehere
```
if compiling from linux:
```
find . -name "makefile" -exec sed -i "s/CFLAGS = -DNo_dtoa -fPIC -O/CFLAGS = -fPIC -O/g" {} \;
make
cd ../../../
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_ASL=ON
make
```	
if compiling from mac:
```
find . -name "makefile" -exec sed -ie 's/CFLAGS = -DNo_dtoa -fPIC -O/CFLAGS = -fPIC -O/g' {} \;
make
cd ../../../
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_ASL=ON
make
```
If ASL (AMPL-MP) is installed already, users can compile and link libpynumero_ASL as follows:

```
cmake .. -DMP_PATH=<PATH_TO_ASL> -DBUILD_SHARED_LIBS=ON -DBUILD_ASL=ON
make
```

To compile the HSL library, download the HSL code from http://hsl.rl.ac.uk/ipopt. Unpack the archive, then move and rename the resulting directory so that it becomes /third_party/HSL/coinhsl. 

```
cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_HSL=ON
```

Once the libraries have been compiled can be moved to pyomo/contrib/pynumero/extension/lib/<OS>
