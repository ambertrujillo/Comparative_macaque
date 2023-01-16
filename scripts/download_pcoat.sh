#!/bin/sh

# Script to download Plasmodium coatneyi genome

module load kent/385

mkdir -p genomes/
cd genomes/

GENOME_FA=plas.fa

PLASMODIUM_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/680/005
PLASMODIUM_URL=$PLASMODIUM_URL/GCF_001680005.1_ASM168000v1/GCF_001680005.1_ASM168000v1_genomic.fna.gz  

wget $PLASMODIUM_URL \
    -O ${GENOME_FA}.gz
gunzip -c ${GENOME_FA}.gz | \
    sed -e "s/^>.*chromosome \([^,]*\),.*/>Plaschr\\1/" > \
    $GENOME_FA

echo "Getting rid of unassembled stuff..." >&2

LAST_OK_LINE=$((`grep -n "^>[^c]" $GENOME_FA | head -n 1 | cut -d":" -f 1` - 1))
if [ $LAST_OK_LINE -gt 0 ]; then
    mv $GENOME_FA ${GENOME_FA}.backup
    head -n $LAST_OK_LINE ${GENOME_FA}.backup > ${GENOME_FA}
    rm ${GENOME_FA}.backup   
fi

mkdir tmp_for_sort
faSplit byname ${GENOME_FA} tmp_for_sort/
cd tmp_for_sort/;
ls -v | xargs cat > ../${GENOME_FA}
cd ..
rm -r tmp_for_sort

cd ..
