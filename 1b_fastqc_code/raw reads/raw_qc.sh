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

WORKDIR_B=/scratch2/alanxu/2019_killi_rna_raw/brain/
WORKDIR_H=/scratch2/alanxu/2019_killi_rna_raw/heart/
WORKDIR_M=/scratch2/alanxu/2019_killi_rna_raw/muscle/
WORKDIR_S=/scratch2/alanxu/2019_killi_rna_raw/spleen/

cd $WORKDIR_B
for f in */
do
    cd $f
    read_1=*_1.fq.gz
    read_2=*_2.fq.gz
    fastqc --outdir /scratch2/alanxu/2019_killi_rna_raw/fastqc/raw/ --threads 12 $read_1 $read_2
    cd $WORKDIR_B
done

cd $WORKDIR_H
for f in */
do
    cd $f
    read_1=*_1.fq.gz
    read_2=*_2.fq.gz
    fastqc --outdir /scratch2/alanxu/2019_killi_rna_raw/fastqc/raw/ --threads 12 $read_1 $read_2
    cd $WORKDIR_H
done

cd $WORKDIR_M
for f in */
do
    cd $f
    read_1=*_1.fq.gz
    read_2=*_2.fq.gz
    fastqc --outdir /scratch2/alanxu/2019_killi_rna_raw/fastqc/raw/ --threads 12 $read_1 $read_2
    cd $WORKDIR_M
done

cd $WORKDIR_S
for f in */
do
    cd $f
    read_1=*_1.fq.gz
    read_2=*_2.fq.gz
    fastqc --outdir /scratch2/alanxu/2019_killi_rna_raw/fastqc/raw/ --threads 12 $read_1 $read_2
    cd $WORKDIR_S
done
