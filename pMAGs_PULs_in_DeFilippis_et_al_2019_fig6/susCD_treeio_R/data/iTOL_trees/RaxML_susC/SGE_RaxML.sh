#!/bin/bash
#$ -N RaxML
#$ -l arch=linux-x64
#$ -pe multislot 20
#$ -l vf=5G
#$ -b n
## -q all.q
#$ -i /dev/null
#$ -e /vol/projects/egalvez/SusDC_phylo/RaxML_Logfiles
#$ -o /vol/projects/egalvez/SusDC_phylo/RaxML_Logfiles
#$ -cwd


# setup pathes
export MY_INPUT=${1}
export MY_OUTPUT=${MY_INPUT%.*}_Rml
export TRIM_log=${MY_INPUT%.*}_RaxML.log
export RAXMLHPC=/vol/cluster-data/egalvez/egalvez_tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3 
export USED_CPUS=20

cd $PWD

$RAXMLHPC -f a -m PROTGAMMAVT -p 12345 -x 12345 -# 100 -s $MY_INPUT -n $MY_OUTPUT -T $USED_CPUS 
