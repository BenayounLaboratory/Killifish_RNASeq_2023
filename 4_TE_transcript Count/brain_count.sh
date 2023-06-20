#!/bin/bash
#SBATCH --mail-user=alanxu@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mem=60GB
#SBATCH --account=bbenayou_34
#SBATCH --partition=oneweek

# Author: Alan Xu
# This script uses TEtranscripts 2.2.1 to generate the count matrix

module purge

TEtranscripts -c brain_1.bam brain_2.bam brain_3.bam brain_4.bam brain_5.bam brain_6.bam brain_7.bam brain_8.bam brain_9.bam brain_10.bam \
-t brain_11.bam brain_12.bam brain_13.bam brain_14.bam brain_15.bam brain_16.bam brain_17.bam brain_18.bam brain_19.bam brain_20.bam \
--GTF /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf \
--TE /scratch2/alanxu/2019_killi_rna_raw/TE_gtf_for_TETranscripts.gtf \
--sortByPos --project brain
