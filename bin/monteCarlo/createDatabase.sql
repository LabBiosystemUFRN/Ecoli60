#Create and populate the Monte Carlo database - mySQL

#drop table MonteCarloBaseTv;

#create table for results
create table MonteCarloBaseTv (simulation int, run int, position int, symbol varchar(20), 
                    type varchar(15), ref varchar(3), mut varchar(3), nucRef varchar(1),
                    nucMut varchar(1), allele varchar(4), strand varchar(7));

#Load consolidated data
load data local infile '../../data_files/MonteCarlo/ecoli60Tv/allSyn.csv' 
    into table MonteCarloBaseTv
    FIELDS TERMINATED by ','
    LINES TERMINATED BY '\n';

#create indexes
create index IdxSimulation ON MonteCarloBaseTv (simulation);
create index IdxRun ON MonteCarloBaseTv (run);
create index IdxRef ON MonteCarloBaseTv (ref);
create index IdxMut ON MonteCarloBaseTv (mut);

#summary table
create table MCMutationsTv (simulation int, run int, ref varchar(3), mut varchar(3), freq int);

#populate sumary table
insert into MCMutationsTv (select simulation, run, ref, mut, count(*) 
                            from MonteCarloBaseTv
                            group by simulation, run, ref, mut);

#and it indexes
create index IdxSimulationTot ON MCMutationsTv (simulation);
create index IdxRunTot ON MCMutationsTv (run);
create index IdxRefTot ON MCMutationsTv (ref);
create index IdxMutTot ON MCMutationsTv (mut);

    


