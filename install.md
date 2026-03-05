# Installation
This part will be replaced with an automatical installation scheme. 

You can install as follows: 
```
$ cd src
$ mkdir ../lib
$ gfortran -O3 -fPIC -c parallax_core.f90
$ gfortran -O3 -fPIC -c c_interface.f90
$ gfortran -shared parallax_core.o c_interface.o -o ../lib/libparallax.so
```
