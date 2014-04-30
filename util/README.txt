##########
GOAL:
Create PMSignature input files from mutation data.

##########
USAGE:
sh format ${MUTATION} ${OUTPUT} ${SIZE}

${MUTATION}
TAB-delimited file for mutation result.
First four columns should be "sample name", "chromosome", "position", "reference allele", "alternated allele".
Example;


sample1       chr5    134296290       C       T
sample1       chr7    143141292       C       A
sample1       chr1    31194287        C       G
sample2       chr1    114165566       A       G
sample2       chr11   5255224 G       C
....

${OUTPUT}
PATH for the input file for PMSignature

${SIZE}
The number of flanking bases around the mutations.


