#!/bin/sh

# Script to download Plasmodium cynomolgi genome

module load kent/385

cd genomes/

GENOME_FA=pcyno.fa

CYNOMOLGI_URL=ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Plasmodium_cynomolgi/latest_assembly_versions/
CYNOMOLGI_URL=$CYNOMOLGI_URL/GCF_000321355.1_PcynB_1.0/GCF_000321355.1_PcynB_1.0_genomic.fna.gz

wget $CYNOMOLGI_URL \
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
