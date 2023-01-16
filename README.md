# Comparative_macaque


# Pipeline
## Prepare Reference Genome and Macaque Reads for Analysis

1. Download Reference Genomes and Annotation files (Host and each pathogen)
> Necessary module(s): kent/385
  * _Macaca mulatta_
```bash
scripts/download_Mmul.sh
```
```bash
scripts/download_Mmul_annotations.sh
```
  * _Plasmodium cynomolgi_
```bash
scripts/download_pcyno.sh
```
```bash
scripts/download_pcyno_annotations.sh
```
  * _Plasmodium coatneyi_
```bash
scripts/download_pcoat.sh
```
```bash
scripts/download_pcoat_annotations.sh
```
2. Concatenate Host-Pathogen Reference Genomes and Annotation files
  * _P.cynomoligi_ and _Macaca mulatta_
```bash
# Reference Genomes
cd genomes
cat plas.fa Mmul.fa > combined.fa
#Annotation files
cat plas.fix.gtf Mmul.fix.gtf > combined.gtf
```
