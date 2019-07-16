#!/bin/bash
#source activate python2.7

binDir=$(cd ..; pwd)
mkdir -p ${binDir}"/data_files/AuxFiles"
mkdir -p ${binDir}"/data_files/base"

#python2.7 0.1GeraHaplotypeTimecourse.py 
#All scripts labeled as Original use the modified version of parse_file.py
python2.7 0.1calculate_clade_hmm_wrapper.py #Original
python2.7 0.2calculate_well_mixed_hmm_wrapper.py #Original
python2.7 0.3calculate_convergence_matrices.py #Original
python2.7 0.4calculate_parallel_genes_table.py #Original
python2.7 0.5calculate_survival_sum_test.py #Original
python2.7 1.1GenerateFrequency.py # Modified from plot_full_trafic_plot_figure.py
python2.7 1.2CalculateCodonUsage.py 
Rscript   1.3CreateTopGenesList.R
python2.7 1.4CalculateDeltaW.py
python2.7 1.5createCodonAll.py
Rscript   2.1CorrectTGCError.R #Correct the error in Cysteine synonymous in parse_file.py from 'TGC':'D' to 'TGC':'C'
Rscript   2.2totalizaCodons.R #create aux files with codons totalizations and relation of MutT mutations
# 2.3enrichDepletion.R
# 2.4plotaComparaCodon.R
source deactivate 
