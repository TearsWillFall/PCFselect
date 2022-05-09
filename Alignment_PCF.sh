#!/bin/bash -l


#$ -l h_rt=48:0:0
#$ -l mem=8G

module unload compilers mpi
module load r/recommended
export PATH=$PATH:/home/regmova/sratoolkit.2.10.8-ubuntu64/bin
export R_LIBS=/home/regmova/MyRlibs:$R_LIBS

if [ "${11}" = "${12}" ]; then
	Rscript ~/Scripts/PCF/Alignment_PCF.R -F $1  -f $2 -o  $3  -g  $4  -q  $5  -Q $6  -T $7 -t  $8  -c  $9 -v  ${10} -L "" -l ""
else
    Rscript ~/Scripts/PCF/Alignment_PCF.R -F $1  -f $2 -o  $3  -g  $4  -q  $5  -Q $6  -T $7 -t  $8  -c  $9 -v  ${10} -L ${11} -l ${12}
fi


