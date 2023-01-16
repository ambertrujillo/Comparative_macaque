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
## DE analysis
1. Host DE analysis (parasite load)
```R
# Remove outliers
Mac_matrix_DE <- Mac_matrix[,-65] # remove transfusion individual
Mac_matrix_DE <- Mac_matrix_DE[,-c(39, 43, 48, 52, 56, 59, 63, 36, 41, 45, 50, 54, 57, 61)]

# Create DGEList
d0.species <- DGEList(Mac_matrix_DE)

# Calculate normalization factors for voom TTM normalization
d0.species <- calcNormFactors(d0.species)

# Filter lowly expressed genes
cutoff.species <- 1
drop.species <- which(apply(cpm(d0.species), 1, max) < cutoff.species)
d.species <- d0.species[-drop.species,] 
dim(d.species) # number of genes left
#16464    50

# change sample and rownames (i.e., sample ids) in DGE list
rownames(d.species$samples) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")
snames.species = colnames(d.species$counts) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")

# Get ID and time variables to create a group 
plas_species <- substr(snames.species, 8, nchar(snames.species)) 
time.species <- substr(snames.species, 6, 6)
individual.species <- substr(snames.species, 1, 5)
group.species <- factor(plas_species)
#group <- interaction(sample_ID, time)

# MDS plot
#by time
plotMDS(d.species, col = as.numeric(time.species))
#by individual
plotMDS(d.species, col = as.numeric(individual.species))
#by group
plotMDS(d.species, col = as.numeric(group.species))

# Running limma-voom for Individual being a random effect since there are more than one data points per individual (http://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/variancePartition/inst/doc/dream.html)

# Try to include sample_ID as block for each model
design_species = model.matrix( ~ parasitemia_NOoutlier$parasitemia_proxy + group.species + parasitemia_NOoutlier$parasitemia_proxy:group.species)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
# first model
y.species <- voom(d.species, design_species, plot = T)
dim(y.species$E) #50
head(y.species$E)
substr(colnames(y.species$E),1,5) == parasitemia_NOoutlier$Individual_ID # Make sure they are in the correct order
dupcor.species <- duplicateCorrelation(y.species,design_species,block=parasitemia_NOoutlier$Individual_ID) # make individual a blocking variable (random effect)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj.species = voom(d.species, design_species, plot=TRUE, block=parasitemia_NOoutlier$Individual_ID, correlation=dupcor.species$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor.species <- duplicateCorrelation(vobj.species, design_species, block=parasitemia_NOoutlier$Individual_ID)

# But this step uses only the genome-wide average for the random effect
fitDupCor.species <- lmFit(vobj.species, design_species, block=parasitemia_NOoutlier$Individual_ID, correlation=dupcor.species$consensus)

# Specify for each coefficient of interest
tmp.species_1 = contrasts.fit(fitDupCor.species, coef = "parasitemia_NOoutlier$parasitemia_proxy")
tmp.species_2 = contrasts.fit(fitDupCor.species, coef = "parasitemia_NOoutlier$parasitemia_proxy:group.speciescynomolgi")
tmp.species_3 = contrasts.fit(fitDupCor.species, coef = "group.speciescynomolgi")
# Fit Empirical Bayes for moderated t-statistics
fitDupCor2.species_1 <- eBayes(tmp.species_1)
fitDupCor2.species_2 <- eBayes(tmp.species_2)
fitDupCor2.species_3 <- eBayes(tmp.species_3)

top.table.species_coefParasitemia <- topTable(fitDupCor2.species_1, coef="parasitemia_NOoutlier$parasitemia_proxy", adjust.method="BH", sort.by = "P", n = Inf) # sort by adj. pval
top.table.species_coefParasitemia$GeneID = rownames(top.table.species_coefParasitemia)

top.table.species_coefParasitemiaSpecies <- topTable(fitDupCor2.species_2, coef="parasitemia_NOoutlier$parasitemia_proxy:group.speciescynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.species_coefParasitemiaSpecies$GeneID = rownames(top.table.species_coefParasitemiaSpecies)

top.table.species_coefSpecies <- topTable(fitDupCor2.species_3, coef="group.speciescynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.species_coefSpecies$GeneID = rownames(top.table.species_coefSpecies)

# Write results tables
all.genes.species_coefParasitemia = subset(top.table.species_coefParasitemia, adj.P.Val < 1, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefParasitemia) #1399

all.genes.species_coefParsitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefParsitemiaSpecies) #127

all.genes.species_coefSpecies = subset(top.table.species_coefSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefSpecies) #224

write.table(all.genes.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_sig_genes_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.species_coefParsitemiaSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_sig_genes_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.species_coefSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_sig_genes_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)


sig.genes.up.species_coefParasitemia = subset(top.table.species_coefParasitemia, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefParasitemia) #1078

sig.genes.up.species_coefParasitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefParasitemiaSpecies) #116

sig.genes.up.species_coefSpecies = subset(top.table.species_coefSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefSpecies) #154

write.table(sig.genes.up.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_upreg_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.species_coefParasitemiaSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_upreg_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.species_coefSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_upreg_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.genes.dn.species_coefParasitemia = subset(top.table.species_coefParasitemia, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefParasitemia) #321

sig.genes.dn.species_coefParasitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefParasitemiaSpecies) #11

sig.genes.dn.species_coefSpecies = subset(top.table.species_coefSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefSpecies) #70

write.table(sig.genes.dn.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_dn_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.species_coefParasitemiaSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_dn_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.species_coefSpecies,
            file=paste0("Species_specific_results/exon/NoOUTmacaque_dn_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
```
2. Host DE analysis (parasitemia)
```R
# Remove outliers
Mac_matrix_DE <- Mac_matrix[,-65] # remove transfusion individual
Mac_matrix_DE <- Mac_matrix_DE[,-c(39, 43, 48, 52, 56, 59, 63, 36, 41, 45, 50, 54, 57, 61)]

parasitemia_binary <- parasitemia_NOoutlier

# Identify these individuals as "uninfected" and the rest as "infected" in "status" 
parasitemia_binary$status = "infected"
parasitemia_binary[c(6, 2, 4, 16, 3, 5, 31, 26, 21, 34, 36, 50, 19, 47, 1, 25, 32, 38, 49, 37, 11), "status"] = "uninfected"

# Create DGEList
d0.binary <- DGEList(Mac_matrix_DE)

# Calculate normalization factors for voom TTM normalization
d0.binary <- calcNormFactors(d0.binary)

# Filter lowly expressed genes
cutoff.binary <- 1
drop.binary <- which(apply(cpm(d0.binary), 1, max) < cutoff.binary)
d.binary <- d0.binary[-drop.binary,] 
dim(d.binary) # number of genes left
#16464    50

# change sample and rownames (i.e., sample ids) in DGE list
rownames(d.binary$samples) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")
snames.binary = colnames(d.binary$counts) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")

# Get ID and time variables to create a group 
plas_binary <- substr(snames.binary, 8, nchar(snames.binary)) 
time.binary <- substr(snames.binary, 6, 6)
individual.binary <- substr(snames.binary, 1, 5)
group.binary <- factor(plas_binary)
#group <- interaction(sample_ID, time)

# MDS plot
#by time
plotMDS(d.binary, col = as.numeric(time.binary))
#by individual
plotMDS(d.binary, col = as.numeric(individual.binary))
#by group
plotMDS(d.binary, col = as.numeric(group.binary))

# Running limma-voom for Individual being a random effect since there are more than one data points per individual (http://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/variancePartition/inst/doc/dream.html)

# Try to include sample_ID as block for each model
parasitemia_binary$status = factor(parasitemia_binary$status, levels=c("uninfected", "infected"))
design_binary = model.matrix( ~ parasitemia_binary$status + group.binary + parasitemia_binary$status:group.binary)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
# first model
y.binary <- voom(d.binary, design_binary, plot = T)
dim(y.binary$E) #50
head(y.binary$E)
substr(colnames(y.binary$E),1,5) == parasitemia_binary$Individual_ID # Make sure they are in the correct order
dupcor.binary <- duplicateCorrelation(y.binary,design_binary,block=parasitemia_binary$Individual_ID) # make individual a blocking variable (random effect)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj.binary = voom(d.binary, design_binary, plot=TRUE, block=parasitemia_binary$Individual_ID, correlation=dupcor.binary$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor.binary <- duplicateCorrelation(vobj.binary, design_binary, block=parasitemia_binary$Individual_ID)

# But this step uses only the genome-wide average for the random effect
fitDupCor.binary <- lmFit(vobj.binary, design_binary, block=parasitemia_binary$Individual_ID, correlation=dupcor.binary$consensus)

# Specify for each coefficient of interest
tmp.binary_1 = contrasts.fit(fitDupCor.binary, coef = "parasitemia_binary$statusinfected")
tmp.binary_2 = contrasts.fit(fitDupCor.binary, coef = "parasitemia_binary$statusinfected:group.binarycynomolgi")
tmp.binary_3 = contrasts.fit(fitDupCor.binary, coef = "group.binarycynomolgi")

# Fit Empirical Bayes for moderated t-statistics
fitDupCor2.binary_1 <- eBayes(tmp.binary_1)
fitDupCor2.binary_2 <- eBayes(tmp.binary_2)
fitDupCor2.binary_3 <- eBayes(tmp.binary_3)

summary(decideTests(fitDupCor2.binary_1))
summary(decideTests(fitDupCor2.binary_2))
summary(decideTests(fitDupCor2.binary_3))

top.table.binary_coefStatus <- topTable(fitDupCor2.binary_1, coef="parasitemia_binary$statusinfected", adjust.method="BH", sort.by = "P", n = Inf) # sort by adj. pval
top.table.binary_coefStatus$GeneID = rownames(top.table.binary_coefStatus)

top.table.binary_coefStatusSpecies <- topTable(fitDupCor2.binary_2, coef="parasitemia_binary$statusinfected:group.binarycynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.binary_coefStatusSpecies$GeneID = rownames(top.table.binary_coefStatusSpecies)

top.table.binary_coefSpecies <- topTable(fitDupCor2.binary_3, coef="group.binarycynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.binary_coefSpecies$GeneID = rownames(top.table.binary_coefSpecies)

# Write results tables
all.genes.binary_coefStatus = subset(top.table.binary_coefStatus, adj.P.Val < 1, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.binary_coefStatus) #2418

all.genes.binary_coefStatusSpecies = subset(top.table.binary_coefStatusSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.binary_coefStatusSpecies) #0

all.genes.binary_coefSpecies = subset(top.table.binary_coefSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.binary_coefSpecies) #255

write.table(all.genes.binary_coefStatus,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_sig_genes_Status", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.binary_coefStatusSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_sig_genes_StatusSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.binary_coefSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_sig_genes_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)


sig.genes.up.binary_coefStatus = subset(top.table.binary_coefStatus, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.binary_coefStatus) #1080

sig.genes.up.binary_coefStatusSpecies = subset(top.table.binary_coefStatusSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.binary_coefStatusSpecies) #0

sig.genes.up.binary_coefSpecies = subset(top.table.binary_coefSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.binary_coefSpecies) #189

write.table(sig.genes.up.binary_coefStatus,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_upreg_Status", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.binary_coefStatusSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_upreg_StatusSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.binary_coefSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_upreg_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.genes.dn.binary_coefStatus = subset(top.table.binary_coefStatus, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.binary_coefStatus) #1388

sig.genes.dn.binary_StatusSpecies = subset(top.table.binary_coefStatusSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.binary_StatusSpecies) #0

sig.genes.dn.binary_coefSpecies = subset(top.table.binary_coefSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.binary_coefSpecies) #66

write.table(sig.genes.dn.binary_coefStatus,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_dn_Status", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.binary_StatusSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_dn_StatusSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.binary_coefSpecies,
            file=paste0("Species_specific_results/exon/binary_3/NoOUTmacaque_dn_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
```
3. Pathogen DE analysis
```R
# Drop outliers, individual that has no microscopy data for time point 3 (Fv13_3) transfusions
Plas_matrix_DE <- Plas_matrix[,-65] # remove transfusion individual
Plas_matrix_DE <- Plas_matrix_DE[,-c(39, 43, 48, 52, 56, 59, 63, 36, 41, 45, 50, 54, 57, 61)] # remove outliers
parasitemia_general <- parasitemia_NOoutlier 

# Create DGEList for inferred parasitemia 
d0_plas <- DGEList(Plas_matrix_DE)

# Calculate normalization factors for voom TTM normalization
d0_plas <- calcNormFactors(d0_plas)
dim(d0_plas)
#3747   50

# have at least min.count reads in a worthwhile number samples. based on design matrix
group.species <- factor(parasitemia_general$Plas_species)
design1 = model.matrix( ~ parasitemia_general$parasitemia_proxy + group.species + parasitemia_general$parasitemia_proxy:group.species)

keep <- filterByExpr(d0_plas, design1)
d_plas = d0_plas[keep, , keep.lib.sizes=FALSE]
dim(d_plas)
#2048   50

# change sample and rownames (i.e., sample ids) in DGE list
rownames(d_plas$samples) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")
snames = colnames(d_plas$counts) <- c("RCs131.coatneyi", "RTi131.coatneyi", "RUn131.coatneyi", "RWr131.coatneyi", "RZe131.coatneyi", "RCs132.coatneyi", "RTi132.coatneyi", "RUn132.coatneyi", "RWr132.coatneyi", "RZe132.coatneyi", "RCs133.coatneyi", "RTi133.coatneyi", "RUn133.coatneyi", "RWr133.coatneyi", "RZe133.coatneyi", "RCs134.coatneyi", "RTi134.coatneyi", "RUn134.coatneyi", "RWr134.coatneyi", "RZe134.coatneyi", "RCs135.coatneyi", "RTi135.coatneyi", "RUn135.coatneyi", "RWr135.coatneyi", "RZe135.coatneyi", "RCs136.coatneyi", "RTi136.coatneyi", "RUn136.coatneyi", "RWr136.coatneyi", "RZe136.coatneyi", "RCs137.coatneyi", "RTi137.coatneyi", "RUn137.coatneyi", "RWr137.coatneyi", "RZe137.coatneyi", "RFv131.cynomolgi", "RIc141.cynomolgi", "RSb141.cynomolgi", "RIc142.cynomolgi", "RSb142.cynomolgi", "RFv133.cynomolgi", "RIc143.cynomolgi", "RSb143.cynomolgi", "RIc144.cynomolgi", "RSb144.cynomolgi", "RIc145.cynomolgi", "RIc146.cynomolgi", "RSb146.cynomolgi", "RIc147.cynomolgi", "RSb147.cynomolgi")

# Get ID and time variables to create a group 
sample_ID <- substr(snames, 1, nchar(snames) - 1) 
time <- substr(snames, nchar(snames), nchar(snames))
individual <- substr(snames, 1, 5)
#group <- interaction(sample_ID, time)
#group

# MDS plot
#by time
plotMDS(d_plas, col = as.numeric(time))
#by individual
plotMDS(d_plas, col = as.numeric(individual))

# Running limma-voom for Individual being a random effect since there are more than one data points per individual (http://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/variancePartition/inst/doc/dream.html)
# Try to include sample_ID as block for each model

# first model
vobj_tmp = voom(d_plas, design1, plot=TRUE)
dim(vobj_tmp$E)
head(vobj_tmp$E)
substr(colnames(vobj_tmp$E),1,5) == parasitemia_general$Individual_ID
dupcor <- duplicateCorrelation(vobj_tmp, design1, block=parasitemia_general$Individual_ID)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
vobj_tmp = voom(d_plas, design1, plot=TRUE, block=parasitemia_general$Individual_ID, correlation=dupcor$consensus)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj_tmp, design1, block=parasitemia_general$Individual_ID)

# But this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj_tmp, design1, block=parasitemia_general$Individual_ID, correlation=dupcor$consensus)

# Specify for each coefficient of interest
tmp.plas.1 = contrasts.fit(fitDupCor, coef = "parasitemia_general$parasitemia_proxy")
tmp.plas.2 = contrasts.fit(fitDupCor, coef = "parasitemia_general$parasitemia_proxy:group.speciescynomolgi")
tmp.plas.3 = contrasts.fit(fitDupCor, coef = "group.speciescynomolgi")

# Fit Empirical Bayes for moderated t-statistics
fitDupCor.plas.1 <- eBayes(tmp.plas.1)
fitDupCor.plas.2 <- eBayes(tmp.plas.2)
fitDupCor.plas.3 <- eBayes(tmp.plas.3)

top.table.species_coefParasitemia <- topTable(fitDupCor.plas.1, coef="parasitemia_general$parasitemia_proxy", adjust.method="BH", sort.by = "P", n = Inf) # sort by adj. pval
top.table.species_coefParasitemia$GeneID = rownames(top.table.species_coefParasitemia)

top.table.species_coefParasitemiaSpecies <- topTable(fitDupCor.plas.2, coef="parasitemia_general$parasitemia_proxy:group.speciescynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.species_coefParasitemiaSpecies$GeneID = rownames(top.table.species_coefParasitemiaSpecies)

top.table.species_coefSpecies <- topTable(fitDupCor.plas.3, coef="group.speciescynomolgi", adjust.method="BH", sort.by="p", n = Inf)
top.table.species_coefSpecies$GeneID = rownames(top.table.species_coefSpecies)

# Write results tables
all.genes.species_coefParasitemia = subset(top.table.species_coefParasitemia, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefParasitemia) #151

all.genes.species_coefParsitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefParsitemiaSpecies) #20

all.genes.species_coefSpecies = subset(top.table.species_coefSpecies, adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(all.genes.species_coefSpecies) #661

write.table(all.genes.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Parasitemia/NoOUTplas_sig_genes_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.species_coefParsitemiaSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_ParasitemiaSpecies/NoOUTplas_sig_genes_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all.genes.species_coefSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Species/NoOUTplas_sig_genes_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.genes.up.species_coefParasitemia = subset(top.table.species_coefParasitemia, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefParasitemia) #0

sig.genes.up.species_coefParasitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefParasitemiaSpecies) #0

sig.genes.up.species_coefSpecies = subset(top.table.species_coefSpecies, logFC > 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.up.species_coefSpecies) #351

write.table(sig.genes.up.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Parasitemia/NoOUTplas_upreg_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.species_coefParasitemiaSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_ParasitemiaSpecies/NoOUTplas_upreg_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.up.species_coefSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Species/NoOUTplas_upreg_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)

sig.genes.dn.species_coefParasitemia = subset(top.table.species_coefParasitemia, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefParasitemia) #151

sig.genes.dn.species_coefParasitemiaSpecies = subset(top.table.species_coefParasitemiaSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefParasitemiaSpecies) #20

sig.genes.dn.species_coefSpecies = subset(top.table.species_coefSpecies, logFC < 0 & adj.P.Val < 0.05, select=c(GeneID, logFC, adj.P.Val))
dim(sig.genes.dn.species_coefSpecies) #310

write.table(sig.genes.dn.species_coefParasitemia,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Parasitemia/NoOUTplas_dn_Parasitemia", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.species_coefParasitemiaSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_ParasitemiaSpecies/NoOUTplas_dn_ParasitemiaSpecies", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig.genes.dn.species_coefSpecies,
            file=paste0("Species_specific_results/exon/traditional_1/coef_Species/NoOUTplas_dn_Species", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)
```


