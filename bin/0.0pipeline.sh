#!/bin/bash
#source activate python2.7

baseDir=$(cd ..; pwd)"/"
mkdir -p ${baseDir}"data_files/AuxFiles"
mkdir -p ${baseDir}"data_files/base"
mkdir -p ${baseDir}"figures"

#All scripts labeled as Original are adapted from the original pipeline and
# use our modified version of parse_file.py

# Final operations in the original pipeline 
python2.7 0.1calculate_clade_hmm_wrapper.py #Original
python2.7 0.2calculate_well_mixed_hmm_wrapper.py #Original
python2.7 0.3calculate_convergence_matrices.py #Original
python2.7 0.4calculate_parallel_genes_table.py #Original
python2.7 0.5calculate_survival_sum_test.py #Original

# Generate files needed in our pipeline
python2.7 1.1GenerateFrequency.py # Modified from plot_full_trafic_plot_figure.py
python2.7 1.2CalculateCodonUsage.py 
Rscript   1.3CreateTopGenesList.R $baseDir
python2.7 1.4CalculateDeltaWCAI.py 86
python2.7 1.5createCodonAll.py

#Create the R enviroment needed
Rscript   0.6base.R $baseDir

#Correct the error in Cysteine synonymous generated in the 
#original pipeline parse_file.py file, where 'TGC':'D' is corrected to 'TGC':'C'
Rscript   2.1CorrectTGCError.R $baseDir 

#just delete some garbage
rm ${baseDir}"figures/supplemental_gene_parallelism_pvalue.pdf"

# Our main pipeline script
Rscript main.R $baseDir

