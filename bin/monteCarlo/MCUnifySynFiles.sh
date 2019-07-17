#!/bin/bash
# This script unify all Monte Carlo files
base==$(cd ..; pwd)"/"
cd $base"data_files/MonteCarlo/ecoli60Tv"

rm allSyn.csv
for i in $(ls proc*Syn.csv); do 
    tail -n +2 $i >> allSyn.csv; 
done

wc -l allSyn.csv

cd -
