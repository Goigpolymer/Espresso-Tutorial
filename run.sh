#!/bin/bash

mpirun -n 2 ./pypresso LJ3.py
mkdir Run
cp Diff.py Run
cp pylj_liquid.* Run
cp  *.txt  Run
cp *.vtf Run
cp *.png Run
./clean.sh
cd Run
./Diff.py
cd ..

