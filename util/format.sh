#! /bin/sh
#$ -S /bin/sh
#$ -cwd

PATH=/home/yshira/bin/BEDTools-Version-2.14.3/bin:$PATH
REF=/home/yshira/ngs/ref/hg19/hg19.fasta
EXON=/home/yshira/code/PMSignature/util/exon.bed

INPUT=$1
OUTPUT=$2
SIZE=$3

echo "perl makeBed.pl ${INPUT} ${SIZE} > ${OUTPUT}.temp"
perl makeBed.pl ${INPUT} ${SIZE} > ${OUTPUT}.temp

echo "fastaFromBed -fi ${REF} -bed ${OUTPUT}.temp -fo ${OUTPUT}.fasta -tab -name"
fastaFromBed -fi ${REF} -bed ${OUTPUT}.temp -fo ${OUTPUT}.fasta -tab -name

echo "intersectBed -a ${OUTPUT}.temp -b ${EXON} -wa -wb > ${OUTPUT}.exon"
intersectBed -a ${OUTPUT}.temp -b ${EXON} -wa -wb > ${OUTPUT}.exon

echo "perl checkAndformat.pl ${OUTPUT}.fasta ${OUTPUT}.exon mutation ${SIZE} > ${OUTPUT}"
perl checkAndformat.pl ${OUTPUT}.fasta ${OUTPUT}.exon mutation ${SIZE} > ${OUTPUT} 

