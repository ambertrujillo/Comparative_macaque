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
4. Download reads from each bioproject
> Necessary module(s): edirect/20210122, sra-tools/2.10.9, parallel/20201022
```bash
mkdir data/cynomolgi
mkdir data/coatneyi
```
 * _P. cynomolgi_ and _M. mulatta_
```bash
esearch -db sra -query PRJNA388645 | efetch -format runinfo | cut -d "," -f 1 > pcyno_SRR.numbers
						
# prefetch
sbatch sbatch/prefetch_pcyno.sbatch

# fasqdump
sbatch sbatch/dump_pcyno.sbatch
```
 * _P. coatneyi_ and _M. mulatta_
```bash
esearch -db sra -query PRJNA400695 | efetch -format runinfo | cut -d "," -f 1 > pcoat_SRR.number

# prefetch
sbatch sbatch/prefetch_pcoat.sbatch

# fastqdump
sbatch sbatch/dump_pcoat.sbatch
```
5. Trim _Macaque_ reads (ARRAY JOB)
> Necessary module(s): trimmomatic/0.36

> To submit an sbatch array job: `sbatch --array=1-[number of individuals] sbatch/name_of_job.sbatch`

```bash
sbatch/trim_pcyno_reads.sbatch
sbatch/trim_pcoat_reads.sbatch
```
## 5a. EdgeR pipeline
### Align Macaque reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
> Necessary module(s): gcc/10.2.0 star/intel/2.7.6a
 * _P. cynomolgi_ and _M. mulatta_
```bash
gunzip $SCRATCH/data/cynomolig/*TRIM_*.fastq.gz
```
```bash
sbatch/edgeR_pcyno_align_genome.sbatch
```
 * _P. coatneyi_ and _M. mulatta_
```bash
gunzip $SCRATCH/data/coatneyi/*TRIM_*.fastq.gz
```
```bash
sbatch/edgeR_pcoat_align_genome.sbatch
```

