#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium cynomolgi annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium cynomolgi genome


cd genomes/ 
CYNOMOLGI_URL=ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Plasmodium_cynomolgi/latest_assembly_versions/

ANNO1=GCF_000321355.1_PcynB_1.0/GCF_000321355.1_PcynB_1.0_genomic.gtf
ANNO2=GCF_000321355.1_PcynB_1.0/GCF_000321355.1_PcynB_1.0_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $CYNOMOLGI_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

sed '
    s/^\NC_020396.1/Plaschr1/ 
    s/^\NC_020395.1/Plaschr2/
    s/^\NC_020397.1/Plaschr3/
    s/^\NC_020408.1/Plaschr4/
    s/^\NC_020398.1/Plaschr5/
    s/^\NC_020399.1/Plaschr6/
    s/^\NC_020400.1/Plaschr7/
    s/^\NC_020401.1/Plaschr8/
    s/^\NC_020402.1/Plaschr9/
    s/^\NC_020403.1/Plaschr10/
    s/^\NC_020404.1/Plaschr11/
    s/^\NC_020405.1/Plaschr12/
    s/^\NC_020406.1/Plaschr13/
    s/^\NC_020407.1/Plaschr14/g' GCF_000321355.1_PcynB_1.0_genomic.gtf > GCF_000321355.1_PcynB_1.0_genomic.fix.gtf

mv GCF_000321355.1_PcynB_1.0_genomic.fix.gtf pcyno.fix.gtf #rename gtf file

cd ..
