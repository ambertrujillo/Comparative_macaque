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
  * _P. cynomolgi_ and _M. mulatta_
```bash
mkdir genomes/cynomolgi
cd genomes/cynomolgi

# Reference Genomes
cat pcyno.fa Mmul.fa > pcyno_combined.fa

# Annotation files
cat pcyno.fix.gtf Mmul.fix.gtf > pcyno_combined.gtf

cd ..
```
 * _P. coatneyi_ and _M. mulatta_
 ```bash
 mkdir coatneyi
 cd coatneyi
 
 # Reference Genomes
 cat pcoat.fa Mmul.fa > pcoat_combined.fa
 
 # Annotation files
 cat pcoat.fix.gtf Mmul.fix.gtf > pcoat_combined.gtf
 
 cd ..
 ```
3. Index concatenated Reference Genome (SBATCH JOB)
> Necessary module(s): gcc/10.2.0, star/intel/2.7.6a

> To submit sbatch job: `sbatch sbatch/name_of_job.sbatch`

```sbatch
# P. cynomolgi and M. mulatta
sbatch/index_pcyno_genomes.sbatch

# P. coatneyi and M. mulatta
sbatch/index_pcoat_genomes.sbatch
```
