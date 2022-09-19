#!/bin/bash

#PBS -l walltime=120:00:00 
#PBS -l nodes=1:ppn=8
#PBS -l pmem=1200mb 

cd preprocess
make clean
make
./a.out

cd ../evolver
make clean
make
./muse_evolver.out

cd ../
python3 plot.py