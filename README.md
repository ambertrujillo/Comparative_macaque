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
## Prepare reads from each bioproject for analysis
1. Download reads from NCBI
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
2. Trim _Macaque_ reads (ARRAY JOB)
> Necessary module(s): trimmomatic/0.36

> To submit an sbatch array job: `sbatch --array=1-[number of individuals] sbatch/name_of_job.sbatch`

```bash
sbatch/trim_pcyno_reads.sbatch
sbatch/trim_pcoat_reads.sbatch
```
## Align Macaque reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
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
## Obtain "Unique" Host and Pathogen Data (ARRAY JOB)
> Necessary module(s): bamtools/intel/2.5.1, samtools/intel/1.11
```bash
sbatch/extract_pcyno_reads.sbatch
sbatch/extract_pcoat_reads.sbatch
```
## Obtain Read Count Matrix and Calculate Percent Parasitemia
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
 ```
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

save.image("MmulCoat_exonfc.Rdata")
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

save.image("PlasCoat_exonfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Plasmodium_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Macaque_Reads_Mapped, Plasmodium_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Plasmodium_Reads_Mapped / Total_Reads

## Calculate parasitemia proxies
```R
### --> Combine macaque .Rdata files
load("Mmul_exonfc.Rdata")
Mmul_coatneyi = Mmul_exonfc$counts
load("PlasCoat_exonfc.Rdata")
Plas_coatneyi = PlasCoat_exonfc$counts

load("MmulCyno_exonfc.Rdata")
Mmul_cyno = MmulCyno_exonfc$counts 
load("PlasCyno_exonfc.Rdata")
Plas_cyno = PlasCyno_exonfc$counts

#Coatneyi
#Load info about individuals
coatneyi_info = read.csv("parasitemias/Pcoat_proportions_parasitemia.csv")
coatneyi_info <- coatneyi_info[ -c(4,11,13,17,19) ] # Getting rid of uncessary info

#Make read count dataframe
Mmul_coatneyi = data.frame(t(Mmul_coatneyi))
rownames(Mmul_coatneyi) <- coatneyi_info$Individual_SRR
Plas_coatneyi = data.frame(t(Plas_coatneyi))
rownames(Plas_coatneyi) <- coatneyi_info$Individual_SRR

#Remove or keep PCOAH (plasmodium) columns in macaque mapping
Mmul_coatneyi = Mmul_coatneyi[,!grepl("^PCOAH_",names(Mmul_coatneyi))]
Plas_coatneyi = Plas_coatneyi[,grepl("^PCOAH_",names(Plas_coatneyi))]

#Sum rows
Mmul_coatneyi_readsMapped = data.frame(rowSums(Mmul_coatneyi))
names(Mmul_coatneyi_readsMapped)[1] = "Macaque_reads"
Plas_coatneyi_readsMapped = data.frame(rowSums(Plas_coatneyi))
names(Plas_coatneyi_readsMapped)[1] = "Plasmodium_reads"
#Include the measure of parasites/ul on day of RNA collection and Rsubread proxy (gene and exon for comparison)
Plas_coatneyi_readsMapped$Individual_ID = coatneyi_info$Individual_ID
Plas_coatneyi_readsMapped$Time_point = coatneyi_info$Time_point
Plas_coatneyi_readsMapped$Parasites_ul_actual_o = coatneyi_info$Parasites_ul_actual_o
Plas_coatneyi_readsMapped$Parasites_ul_actual_b_o_a = coatneyi_info$Parasites_ul_actual_b_o_a
Plas_coatneyi_readsMapped$Parasites_ul_actual_b_o = coatneyi_info$Parasites_ul_actual_b_o
Plas_coatneyi_readsMapped$Parasites_ul_actual_2b_b_o = coatneyi_info$Parasites_ul_actual_2b_b_o
Plas_coatneyi_readsMapped$Percent_Parasitemia_gene = coatneyi_info$Percent_Parasitemia_gene
Plas_coatneyi_readsMapped$Percent_Parasitemia_exon = coatneyi_info$Percent_Parasitemia_exon
Plas_coatneyi_readsMapped$Plas_species = "coatneyi"

coatneyi_reads = merge(Mmul_coatneyi_readsMapped, Plas_coatneyi_readsMapped, by="row.names")
rownames(coatneyi_reads) = coatneyi_reads$Row.names
coatneyi_reads = coatneyi_reads[-c(1)]

#Cynomolgi
#Load info about individuals
cynomolgi_info = read.csv("parasitemias/Pcyno_proportions_parasitemia.csv")
cynomolgi_info <- cynomolgi_info[ -c(10,12,16,18) ] # getting rid of unnecessary info

#Make read count data frame
Mmul_cyno = data.frame(t(Mmul_cyno))
rownames(Mmul_cyno) <- cynomolgi_info$Individual_SRR
Plas_cyno = data.frame(t(Plas_cyno))
rownames(Plas_cyno) <- cynomolgi_info$Individual_SRR

#Remove or keep PCOAH (plasmodium) columns in macaque mapping
Mmul_cyno = Mmul_cyno[,!grepl("^PCYB_",names(Mmul_cyno))]
Plas_cyno = Plas_cyno[,grepl("^PCYB_",names(Plas_cyno))]

#Sum rows
Mmul_cyno_readsMapped = data.frame(rowSums(Mmul_cyno))
names(Mmul_cyno_readsMapped)[1] = "Macaque_reads"
Plas_cyno_readsMapped = data.frame(rowSums(Plas_cyno))
names(Plas_cyno_readsMapped)[1] = "Plasmodium_reads"

# Include other measures of parasite load
Plas_cyno_readsMapped$Individual_ID = cynomolgi_info$Individual_ID
Plas_cyno_readsMapped$Time_point = cynomolgi_info$Time_point
Plas_cyno_readsMapped$Parasites_ul_actual_o = cynomolgi_info$Parasites_ul_actual_o
Plas_cyno_readsMapped$Parasites_ul_actual_b_o_a = cynomolgi_info$Parasites_ul_actual_b_o_a
Plas_cyno_readsMapped$Parasites_ul_actual_b_o = cynomolgi_info$Parasites_ul_actual_b_o
Plas_cyno_readsMapped$Parasites_ul_actual_2b_b_o = cynomolgi_info$Parasites_ul_actual_2b_b_o
Plas_cyno_readsMapped$Percent_Parasitemia_gene = cynomolgi_info$Percent_Parasitemia_gene
Plas_cyno_readsMapped$Percent_Parasitemia_exon = cynomolgi_info$Percent_Parasitemia_exon
Plas_cyno_readsMapped$Plas_species = "cynomolgi"

cynomolgi_reads = merge(Mmul_cyno_readsMapped, Plas_cyno_readsMapped, by="row.names")
rownames(cynomolgi_reads) = cynomolgi_reads$Row.names
cynomolgi_reads = cynomolgi_reads[-c(1)]

# Combine both parasitemia tables
parasitemia = rbind(coatneyi_reads, cynomolgi_reads)
# Calculate parasitemia
parasitemia$parasitemia_proxy = parasitemia$Plasmodium_reads / parasitemia$Macaque_reads
parasitemia$SRRs = rownames(parasitemia)

# For model 1
parasitemia_outler = parasitemia[!(parasitemia$Individual_ID %in% c("CF97_donor")), ] # Only dropped transfusion sample

# For model 2
parasitemia_NOoutlier = parasitemia[!(parasitemia$Individual_ID %in% c("RFa14", "RMe14", "CF97_donor")), ]

save.image(file="exon/Parasitemia_calculation.Rdata")
```
## Combine read count matrices
```R
# load orthologs
require(data.table)

orthologs.coat = read.csv("orthologs/orthologs_starting_with_coat.csv") #4465
orthologs.cyno = read.csv("orthologs/orthologs_starting_with_cyno.csv") #4509

# Get rid of all rows that do not have 1 to 1 ortholog
coat.1ortho = orthologs.coat[unlist(lapply(orthologs.coat$Input.Ortholog.s., nchar)) == 14,] #3927
cyno.1ortho = orthologs.cyno[unlist(lapply(orthologs.cyno$Input.Ortholog.s., nchar)) == 11,] #3918

# make sure they are in each other's table
reciprocal.ortho = coat.1ortho[coat.1ortho$Input.Ortholog.s. %in% cyno.1ortho$Gene.ID,] #3747
reciprocal.ortho = reciprocal.ortho[reciprocal.ortho$Gene.ID %in% coat.1ortho$Gene.ID,] #3747
reciprocal.ortho = subset(reciprocal.ortho, select=c("Input.Ortholog.s.", "Gene.ID"))
names(reciprocal.ortho)[names(reciprocal.ortho) == "Input.Ortholog.s."] <- "Coatneyi_genes"
names(reciprocal.ortho)[names(reciprocal.ortho) == "Gene.ID"] <- "Cynomolgi_genes"

# Combine read count matrices 
library(dplyr)
library(tidyr)

# For plasmodium analysis
Plas_coat_matrix = data.frame(t(Plas_coatneyi))
Plas_coat_matrix$Coatneyi_genes = rownames(Plas_coat_matrix)
Plas_coat_matrix = full_join(Plas_coat_matrix, reciprocal.ortho, by="Coatneyi_genes")
Plas_coat_matrix = drop_na(Plas_coat_matrix)
Plas_coat_matrix = Plas_coat_matrix[,-c(37)]

Plas_cyno_matrix = data.frame(t(Plas_cyno))
Plas_cyno_matrix$Cynomolgi_genes = rownames(Plas_cyno_matrix)
Plas_cyno_matrix = full_join(Plas_cyno_matrix, reciprocal.ortho, by="Cynomolgi_genes")
Plas_cyno_matrix = drop_na(Plas_cyno_matrix)
Plas_cyno_matrix = Plas_cyno_matrix[,-c(31)]

Plas_matrix = merge(Plas_coat_matrix, Plas_cyno_matrix, by="Coatneyi_genes")
rownames(Plas_matrix) = Plas_matrix$Coatneyi_genes
Plas_matrix = Plas_matrix[,-c(1)]

save.image(file="exon/Ready_for_pathogen_DE.Rdata")
```
