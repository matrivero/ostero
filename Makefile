#Makefile for ostero code

external.so: external.f90
	f2py -c -m external external.f90
