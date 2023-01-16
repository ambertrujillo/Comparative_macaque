#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=124GB
#SBATCH --job-name=extract_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=extract_%A_%a.out

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

module load bamtools/intel/2.5.1
module load samtools/intel/1.11

DATA_DIR=$SCRATCH/data/coatneyi
THIS_IND=`ls $DATA_DIR/*_1.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_1.*//"`
READ_DIR=$SCRATCH/results/coatneyi

echo "Extracting primary reads..."

#Only output primary reads
samtools view -F 256 -b $READ_DIR/${THIS_IND}Aligned.sortedByCoord.out.bam > $READ_DIR/${THIS_IND}.primary.bam


echo "Splitting bam files..."

#Splits bam files into chromosomes
bamtools split -in $READ_DIR/${THIS_IND}.primary.bam -reference 

echo "Removing unmapped reads and unplaced scaffolds..."

#Remove unmapped and unplaced scaffolds
rm $READ_DIR/${THIS_IND}.primary.REF_NW_004* #Remove unplaced macaque stuff
rm $READ_DIR/${THIS_IND}.primary.REF_unmapped.bam

echo "Merging bams..."

#-->Macaque
ls $READ_DIR/${THIS_IND}.primary.REF_Mmulchr*.bam > $READ_DIR/${THIS_IND}.Mmulbamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.Mmulbamlist.txt -out $READ_DIR/${THIS_IND}.Mmul.bam

#-->Plasmodium
ls $READ_DIR/${THIS_IND}.primary.REF_Plaschr*.bam > $READ_DIR/${THIS_IND}.Plasbamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.Plasbamlist.txt -out $READ_DIR/${THIS_IND}.Plas.bam

echo "Removing mess you made ..."

#Remove mess you made
rm $READ_DIR/${THIS_IND}.primary.REF_Mmulch*.bam
rm $READ_DIR/${THIS_IND}.primary.REF_Plaschr*.bam

echo "Moving unique BAMS..."

#Move unique colobus reads to mapped directory
mv $READ_DIR/${THIS_IND}.Mmul.bam $READ_DIR/mapped_reads/macaque
mv $READ_DIR/${THIS_IND}.Plas.bam $READ_DIR/mapped_reads/plasmodium


echo "DONE"
