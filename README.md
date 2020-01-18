PyNumero libraries
==================

This repository contains a cmake project for building shared
libraries for the python numerical optimization package PyNumero.

The current version of PyNumero relies on share libraries for
efficient numerical calculations while solving optimization problems

1. libpynumero_ASL.so
2. libpynumero_MA27.so
3. libpynumero_MA57.so (only available if linked against an existing ipopt installation)

The first library, libpynumero_ASL, provides functionality for fast automatic
differentiation (AD) subroutines required for computing first- and second-
oder derivatives necessary for gradient-based optimization algorithms. The
subroutines for AD are the same ones used by AMPL nonlinear optimization
solvers since libpynumero_ASL links with the ampl solver library ASL.

The second and third libraries provide functionality to interface
with HSL linear solvers available in the Harwell subroutine
library.

A conda-forge recipe has been created to facilitate the distribution of
pynumero libraries without users requiring to compile any code. To get
libpynumero_ASL library via conda-forge users may run:

```
conda install -c conda-forge pynumero_libraries
```

For open source linear solver we recommend installing

```
conda install -c conda-forge pymumps
```

If conda is not available, users that have a compiled version of IPOPT can use the following steps to build the libraries:

```
mkdir build
cmake .. -DBUILD_ASL=ON -DIPOPT_PATH=<your-path-to-ipopt-install> -DCMAKE_INSTALL_PREFIX=../pynumero-install
make
make install
```

The libraries are then compiled and copied into ../pynumero-install. If the user has a compiled version of ipopt with HSL solver and would like to compile MA27 or MA57 libraries, the options -DBUILD_MA27=ON and -DMBUILD_MA57=ON can be used. 

If a compiled version of ipopt is not available, users will need to get the source code for ASL and MA27

```
cd /third_party/ASL
./get.ASL.sh
cd solvers
./configurehere
```

If compiling from linux:
```
find . -name "makefile" -exec sed -i "s/CFLAGS = -DNo_dtoa -fPIC -O/CFLAGS = -fPIC -O/g" {} \; 
```
If compiling from mac:
```
find . -name "makefile" -exec sed -i 's/CFLAGS = -DNo_dtoa -fPIC -O/CFLAGS = -fPIC -O/g' {} \; 
```

Then

```
make
cd ../../../
mkdir build
cd build
cmake .. -DBUILD_ASL=ON -DCMAKE_INSTALL_PREFIX=../pynumero-install
make
make install
```

If ASL (AMPL-MP) is installed already, users can compile and link libpynumero_ASL as follows:

```
cmake .. -DMP_PATH=<PATH_TO_AMPL-MP>  -DBUILD_ASL=ON -DCMAKE_INSTALL_PREFIX=../pynumero-install
make
make install
```

To compile the MA27 library, download the HSL code from http://hsl.rl.ac.uk/ipopt. Unpack the archive, then move and rename the resulting directory so that it becomes /third_party/HSL/coinhsl. 

```
cmake .. -DBUILD_MA27=ON -DCMAKE_INSTALL_PREFIX=../pynumero-install
make
make install
```

Once the libraries have been compiled they can be moved to pyomo/contrib/pynumero/extension/lib/<OS>
