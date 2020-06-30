# Main routines

This folder contains the scripts used to generate and process the data used by our team after pre-processing all sequencing files. This pre-processing can be done following the instructions found in [Time-resolved metagenomic sequencing of Lenski's long-term evolution experiment with Escherichia coli ](https://github.com/benjaminhgood/LTEE-metagenomic)\[1\]. 

After that, please place the tragectory files in *data_files/preProcess* folder.

Alternatively, you can run our pipeline using pre-processed sequencing files provided [here](https://github.com/LabBiosystemUFRN/Ecoli60/tree/master/data_files/preProcess/). To do that, execute:
> 00pipeline.sh


## Description of files

### Adapted from the original pipeline

#### 0.1calculate_clade_hmm_wrapper.py \*

>Calculates the clade hmm wrapper (see the [original pipeline](https://github.com/benjaminhgood/LTEE-metagenomic) for details)

#### 0.2calculate_well_mixed_hmm_wrapper.py \*

>Calculates well mixed hmm wrapper (see the [original pipeline](https://github.com/benjaminhgood/LTEE-metagenomic) for details)

#### 0.3calculate_convergence_matrices.py \*

>Calculates the convergence matrices (see the [original pipeline](https://github.com/benjaminhgood/LTEE-metagenomic) for details)

#### 0.4calculate_parallel_genes_table.py \*

>Calculates parallel gene tables (see the [original pipeline](https://github.com/benjaminhgood/LTEE-metagenomic) for details)

#### 0.5calculate_survival_sum_test.py \*

>Calculates the time of survival (see the [original pipeline](https://github.com/benjaminhgood/LTEE-metagenomic) for details)

### Our routines

#### 1.1GenerateFrequency.py \*\*

>Calculates the mutation frequencies

#### 1.2CalculateCodonUsage.py 

>Calculates the codon usage

#### 1.3CreateTopGenesList.R

>R script to create top expressed genes list

#### 1.4CalculateDeltaWCAI.py 86

>R script to calculate all Delta w

#### 1.5createCodonAll.py

>Creates some auxiliary files needed

#### 0.6base.R 

>R script to create and setting up the R enviroment needed

#### 2.1CorrectTGCError.R 

>This script correct an error found in the original script (parse_file.py) where the where 'TGC' codon is assinged to Aspartic acid instead Cysteine.


#### main.R 

>Run our main pipeline


<sup>\* Those scripts are adapted from the original pipeline and use our modified version of parse_file.py</sup>

<sup>\*\* Modified from the original plot_full_trafic_plot_figure.py script</sup>

