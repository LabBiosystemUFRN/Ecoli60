### Main routines

This folder contains the scripts used to generate and process the data used by our team after pre-processing. This pre-processing can be done folowing the instructions find in [Time-resolved metagenomic sequencing of Lenski's long-term evolution experiment with Escherichia coli ](https://github.com/benjaminhgood/LTEE-metagenomic)

The first step is create de symulated data . It will take a long time.
> monteCarlo.sh

After the process finish the files must be unified and a MySQL database will be created, runing the scripts:
>MCUnifySynFiles.sh
createDatabase.sql


Finaly create the needed data frames.
>justCreateMonteCarloDF.R

Or save a lot of time and use the data provided in files that could be found at *data_files/base/* folder:
>MonteCarloPEmp.csv
>aHighMutSynFreqDW.csv
>aMutTMutSynFreqDW.csv
