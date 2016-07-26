#!/bin/bash

#remove previous and re-generate signature file
rm -f external.pyf
f2py -m external -h external.pyf external.f90

#remove previous and re-generate the object files (-fPIC option is mandatory)
#in marenostrum you will need ifort instead of gfortran
rm -f external.o demo_subroutine.o
gfortran -c -fPIC external.f90 demo_subroutine.f90
#ifort -c -fPIC external.f90 demo_subroutine.f90

#link all
rm -f external.so
f2py -c external.pyf *.o
