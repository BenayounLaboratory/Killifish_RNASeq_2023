#!/bin/bash
#SBATCH --mail-user=alanxu@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --mem=60GB
#SBATCH --account=bbenayou_34

# Author: Alan Xu
# This script uses star/2.7.0e to generate genome index

module purge
module load gcc/8.3.0
module load intel/18.0.4
module load star/2.7.0e

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./genome_index --limitGenomeGenerateRAM 60000000000 --genomeFastaFiles Soft_Masked_Genome.fa --sjdbGTFfile GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf
