#!/bin/bash
#source activate python2.7

baseDir=$(cd ..; pwd)"/"
mkdir -p ${baseDir}"data_files/AuxFiles"
mkdir -p ${baseDir}"data_files/base"
mkdir -p ${baseDir}"figures"

#All scripts labeled as Original are adapted from the original pipeline and
# use our modified version of parse_file.py

# Final operations in the original pipeline
echo "######################################"
echo "Calculating clade hmm wrapper"
echo "######################################"
python2.7 0.1calculate_clade_hmm_wrapper.py #Original
echo "######################################"
echo "Calculating well mixed hmm wrapper"
echo "######################################"
python2.7 0.2calculate_well_mixed_hmm_wrapper.py #Original
echo "######################################"
echo "Calculating convergence matrices"
echo "######################################"
python2.7 0.3calculate_convergence_matrices.py #Original
echo "######################################"
echo "Calculating parallel gene tables"
echo "######################################"
python2.7 0.4calculate_parallel_genes_table.py #Original
echo "######################################"
echo "Calculating survival"
echo "######################################"
python2.7 0.5calculate_survival_sum_test.py #Original

# Generate files needed in our pipeline
echo "######################################"
echo "Calculating mutation frequencies"
echo "######################################"
python2.7 1.1GenerateFrequency.py # Modified from plot_full_trafic_plot_figure.py
echo "######################################"
echo "Calculating codon usage"
echo "######################################"
python2.7 1.2CalculateCodonUsage.py 
echo "######################################"
echo "Creating top expressed genes list"
echo "######################################"
Rscript   1.3CreateTopGenesList.R $baseDir
echo "######################################"
echo "Calculating Delta w"
echo "######################################"
python2.7 1.4CalculateDeltaWCAI.py 86
echo "######################################"
echo "Creating aux files"
echo "######################################"
python2.7 1.5createCodonAll.py

#Create the R enviroment needed
echo "######################################"
echo "Seting up enviroment"
echo "######################################"
Rscript   0.6base.R $baseDir

#Correct the error in Cysteine synonymous generated in the 
#original pipeline parse_file.py file, where 'TGC':'D' is corrected to 'TGC':'C'
echo "######################################"
echo "Correctin TGC error"
echo "######################################"
Rscript   2.1CorrectTGCError.R $baseDir 

#just delete some garbage
rm ${baseDir}"figures/supplemental_gene_parallelism_pvalue.pdf"

# Our main pipeline script
echo "######################################"
echo "Running Main function"
echo "######################################"
Rscript main.R $baseDir

