#!/bin/sh

# Script to download macaque genome, Mmul10

module load kent/385
module load samtools/intel/1.11

mkdir -p genomes
cd genomes

GENOME_FA=Mmul.fa

wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz' \
    -O ${GENOME_FA}.gz

gunzip -c ${GENOME_FA}.gz | \
    sed -e "s/^>.*chromosome \([^,]*\),.*/>Mmulchr\\1/" > \
    $GENOME_FA  # Rename Chromosomes to Mmul

rm -f ${GENOME_FA/.fa/.chr.fa}
for CHR in `seq 1 20` X Y; do
    samtools faidx $GENOME_FA Mmulchr$CHR >> ${GENOME_FA/.fa/.chr.fa}
done  # Create a new file with extracted chromosomes (not including unplaced scaffolds) 1-20 X and Y

mv ${GENOME_FA/.fa/.chr.fa} $GENOME_FA

cd ..

exit
