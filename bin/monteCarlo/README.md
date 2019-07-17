### Monte Carlo symulation

This folder contains the scripts to generate the Monte Carlo symulation to check if the unobserved codons are relevant. 

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
