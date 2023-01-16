#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium coatneyi annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium coatneyi genome, P. coatneyi


cd genomes/ 

PLASMODIUM_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/680/005

ANNO1=GCF_001680005.1_ASM168000v1/GCF_001680005.1_ASM168000v1_genomic.gtf
ANNO2=GCF_001680005.1_ASM168000v1/GCF_001680005.1_ASM168000v1_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $PLASMODIUM_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

sed '
    s/^\NC_033556.1/Plaschr1/ 
    s/^\NC_033557.1/Plaschr2/
    s/^\NC_033558.1/Plaschr3/
    s/^\NC_033559.1/Plaschr4/
    s/^\NC_033560.1/Plaschr5/
    s/^\NC_033561.1/Plaschr6/
    s/^\NC_033562.1/Plaschr7/
    s/^\NC_033563.1/Plaschr8/
    s/^\NC_033564.1/Plaschr9/
    s/^\NC_033565.1/Plaschr10/
    s/^\NC_033566.1/Plaschr11/
    s/^\NC_033567.1/Plaschr12/
    s/^\NC_033568.1/Plaschr13/
    s/^\NC_033569.1/Plaschr14/g' GCF_001680005.1_ASM168000v1_genomic.gtf > GCF_001680005.1_ASM168000v1_genomic.fix.gtf

mv GCF_001680005.1_ASM168000v1_genomic.fix.gtf pcoat.fix.gtf #rename gtf file

cd ..

