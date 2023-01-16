#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Macaque annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Macaca mulatta genome, mmul 10

cd genomes/ 

MMUL_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765

ANNO1=GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf
ANNO2=GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $MMUL_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

sed '
    s/^\NC_041754.1/Mmulchr1/ 
    s/^\NC_041755.1/Mmulchr2/
    s/^\NC_041756.1/Mmulchr3/
    s/^\NC_041757.1/Mmulchr4/
    s/^\NC_041758.1/Mmulchr5/
    s/^\NC_041759.1/Mmulchr6/
    s/^\NC_041760.1/Mmulchr7/
    s/^\NC_041761.1/Mmulchr8/
    s/^\NC_041762.1/Mmulchr9/
    s/^\NC_041763.1/Mmulchr10/
    s/^\NC_041764.1/Mmulchr11/
    s/^\NC_041765.1/Mmulchr12/
    s/^\NC_041766.1/Mmulchr13/
    s/^\NC_041767.1/Mmulchr14/
    s/^\NC_041768.1/Mmulchr15/
    s/^\NC_041769.1/Mmulchr16/
    s/^\NC_041770.1/Mmulchr17/
    s/^\NC_041771.1/Mmulchr18/
    s/^\NC_041772.1/Mmulchr19/
    s/^\NC_041773.1/Mmulchr20/
    s/^\NC_041774.1/MmulchrX/
    s/^\NC_027914.1/MmulchrY/g' GCF_003339765.1_Mmul_10_genomic.gtf > GCF_003339765.1_Mmul_10_genomic.fix.gtf


mv GCF_003339765.1_Mmul_10_genomic.fix.gtf Mmul.fix.gtf #rename gtf file


cd ..
