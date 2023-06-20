#!/bin/bash
#SBATCH --mail-user=alanxu@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --account=bbenayou_34
#SBATCH --partition=oneweek

# Author: Alan Xu
# This script uses fastqc/0.11.9 to quality check the filtered reads

module purge
module load gcc/11.3.0
module load fastqc/0.11.9

WORKDIR=/scratch2/alanxu/2019_killi_rna_raw/brain/
cd $WORKDIR
for f in */
do
    cd $f
    read_1=*_1_trimmed.fq.gz
    read_2=*_2_trimmed.fq.gz
    fastqc --outdir /scratch2/alanxu/2019_killi_rna_raw/fastqc/ --threads 12 $read_1 $read_2
    cd $WORKDIR
done
