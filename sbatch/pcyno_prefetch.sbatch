#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --job-name=prefetch_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=prefetch_%j.out
  
module load edirect/20210122
module load sra-tools/2.10.9
module load parallel/20201022
  
cd $SLURM_SUBMIT_DIR

parallel --verbose -j 20 prefetch --output-directory $SCRATCH/data/cynomolgi {} ::: $(grep -v Run SRR.numbers )
wait 

exit 
