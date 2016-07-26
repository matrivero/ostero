#!/bin/bash

#generate signature file
f2py -m external -h external.pyf external.f90
#generate the object files (-fPIC option is mandatory)
gfortran -c -fPIC external.f90 demo_subroutine.f90
#link all
f2py -c external.pyf *.o
