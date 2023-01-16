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
6. Align Macaque reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
> Necessary module(s): gcc/10.2.0 star/intel/2.7.6a
 * _P. cynomolgi_ and _M. mulatta_
```bash
gunzip $SCRATCH/data/cynomolig/*TRIM_*.fastq.gz
```
```bash
sbatch/align_pcyno_genome.sbatch
```
 * _P. coatneyi_ and _M. mulatta_
```bash
gunzip $SCRATCH/data/coatneyi/*TRIM_*.fastq.gz
```
```bash
sbatch/align_pcoat_genome.sbatch
```
```bash
# Remove lines in gtf file with empty gene_id field
# P. cynomolgi
grep 'gene_id ""' genomes/cynomolgi/pcyno_combined.gtf # to look at them
grep -v 'gene_id ""' genomes/cynomolgi/pcyno_combined.gtf > genomes/cynomolgi/pcyno_combined_FIXED.gtf

# P. coatneyi
grep 'gene_id ""' genomes/coatneyi/pcoat_combined.gtf # to look at them
grep -v 'gene_id ""' genomes/coatneyi/pcoat_combined.gtf > genomes/coatneyi/pcoat_combined_FIXED.gtf
```
7. Obtain "Unique" Host and Pathogen Data (ARRAY JOB)
> Necessary module(s): bamtools/intel/2.5.1, samtools/intel/1.11
```bash
sbatch/extract_pcyno_reads.sbatch
sbatch/extract_pcoat_reads.sbatch
```













8. Obtain Read Count Matrix and Calculate Percent Parasitemia
> Necessary module(s): r/intel/4.0.4
> Necessary R package(s): BiocManager, Rsubread
```bash
module load r/intel/4.0.4
R
```
 * _P. cynomolgi_ and _M. mulatta_

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

module load r/gcc/4.2.0

R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

# Macaque
bams = list.files(path = "results/cynomolgi/mapped_reads/macaque", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/cynomolgi/pcyno_combined_FIXED.gtf"
MmulCyno_exonfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    GTF.featureType = "exon",
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("MmulCyno_exonfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Macaque_Reads_Mapped):
 > Macaque_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "results/cynomolgi/mapped_reads/plasmodium", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/pcyno_combined_FIXED.gtf"
PlasCyno_exonfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    GTF.featureType = "exon",
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("PlasCyno_exonfc.Rdata")```

  * As file is running, create percent_parasitemia table in excel (enter Hepatocystis_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Macaque_Reads_Mapped, Plasmodium_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Hepatocystis_Reads_Mapped / Total_Reads
 
  * _P. coatneyi_ and _M. mulatta_

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

module load r/gcc/4.2.0

R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

# Macaque
bams = list.files(path = "results/coatneyi/mapped_reads/macaque", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/coatneyi/pcoat_combined_FIXED.gtf"
MmulCyno_exonfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    GTF.featureType = "exon",
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("exon/Mmul_exonfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Macaque_Reads_Mapped):
 > Macaque_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "results/coatneyi/mapped_reads/plasmodium", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/pcoat_combined_FIXED.gtf"
PlasCyno_exonfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    GTF.featureType = "exon",
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("exon/PlasCoat_exonfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Plasmodium_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Macaque_Reads_Mapped, Plasmodium_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Plasmodium_Reads_Mapped / Total_Reads
```
