#!/bin/bash
#SBATCH --time=6-0:0
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-user=cfreis230@gmail.com
#SBATCH --mail-type=ALL


module load softwares/python/2.7-anaconda-5.0.1

if [ $# -ne 2 ]; then
	echo 1>&2 Sintaxe: $0 \<Qtd Synonymous\> \<ID da Simulação\>
	exit 999
fi

synMut=$1
tipo=$2 #"ecoli60Tv"

hora=`date`
echo "Inciado em: $hora"


dir="./Resultados/$tipo"

mkdir -p $dir

#realiza sempre 20000 simulações 
let proces=32 #125 #$1
let runs=625 #160 #(64/$1) #100 #$2

procRef=0;
for i in `seq 1 $proces`; do
    let procRef=$procRef+1;
    arq="$dir/proc$procRef"

    echo #$proces $repete
    echo 
    echo "Iniciando processo $procRef"
    echo 
    
    python2.7 ./monteCarlo.py $arq $procRef $runs $synMut & 
    echo $arq $procRef $runs $synMut    
done;
wait

hora=`date`
echo "Finalizado em: $hora"

