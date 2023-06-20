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

TEtranscripts -c Heart_1.bam Heart_2.bam Heart_3.bam Heart_4.bam Heart_5.bam Heart_6.bam Heart_7.bam Heart_8.bam Heart_9.bam Heart_10.bam \
-t Heart_11.bam Heart_12.bam Heart_13.bam Heart_14.bam Heart_15.bam Heart_16.bam Heart_17.bam Heart_18.bam Heart_19.bam Heart_20.bam \
--GTF /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf \
--TE /scratch2/alanxu/2019_killi_rna_raw/TE_gtf_for_TETranscripts.gtf \
--sortByPos --project heart
