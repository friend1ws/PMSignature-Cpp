#! /bin/sh
#$ -S /bin/sh
#$ -cwd

PATH=/home/yshira/bin/BEDTools-Version-2.14.3/bin:$PATH
REF=/home/yshira/ngs/ref/hg19/hg19.fasta

INPUT=$1
OUTPUT=$2
SIZE=$3

perl makeBed.pl ${INPUT} ${SIZE} > ${OUTPUT}.temp

fastaFromBed -fi ${REF} -bed ${OUTPUT}.temp -fo ${OUTPUT}.fasta -tab -name

perl checkAndformat.pl ${OUTPUT}.fasta mutation ${SIZE} > ${OUTPUT} 
