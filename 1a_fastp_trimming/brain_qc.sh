#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --mem=60GB
#SBATCH --account=bbenayou_34
#SBATCH --partition=epyc-64
#SBATCH --mail-type=all
#SBATCH --mail-user=alanxu@usc.edu

# Author: Alan Xu
# Software Used: fastp version 0.23.2

WORKDIR=/scratch2/alanxu/2019_killi_rna_raw/brain/
cd $WORKDIR
for f in */
do
    cd $f
    read_1=*_1.fq.gz
    read_2=*_2.fq.gz
    read_1_out="$(echo $read_1 | sed 's/.fq.gz/_trimmed.fq.gz/')"
    read_2_out="$(echo $read_2 | sed 's/.fq.gz/_trimmed.fq.gz/')"
    ~/fastp --in1 $read_1 --in2 $read_2 --out1 $read_1_out --out2 $read_2_out --failed_out fail_reads.out --detect_adapter_for_pe 
    cd $WORKDIR
done
