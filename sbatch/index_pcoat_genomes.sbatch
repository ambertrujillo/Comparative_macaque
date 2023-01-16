#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=index_genome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=starindex_%j.out


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

module load gcc/10.2.0 star/intel/2.7.6a

mkdir -p genomes/coatneyi/index

STAR --runThreadN 2 \
--runMode genomeGenerate \
--genomeDir genomes/coatneyi/index \
--genomeFastaFiles genomes/coatneyi/pcoat_combined.fa \
--sjdbGTFfile genomes/coatneyi/pcoat_combined.gtf \
--sjdbOverhang 99 \
--outFileNamePrefix genomes/coatneyi/index/combinedindex \
--limitGenomeGenerateRAM=40712737493
