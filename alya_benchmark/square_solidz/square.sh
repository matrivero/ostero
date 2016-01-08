#!/bin/bash

#BSUB -n 1
#BSUB -oo output.txt
#BSUB -eo error.txt
#BSUB -J benchmark
#BSUB -W 01:30

/home/bsc21/bsc21494/AlyaContact/Executables/unix/Alya.x benchmark
