#!/bin/bash
base="/home/clovis/Doutorado/Projetos/Ecoli60"
cd $base"/data_files/MonteCarlo/ecoli60Tv"

rm allSyn.csv
for i in $(ls proc*Syn.csv); do 
    tail -n +2 $i >> allSyn.csv; 
done

wc -l allSyn.csv

cd -
