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

TEtranscripts -c muscle_1.bam muscle_2.bam muscle_3.bam muscle_4.bam muscle_5.bam muscle_6.bam muscle_7.bam muscle_8.bam muscle_9.bam \
-t muscle_10.bam muscle_11.bam muscle_12.bam muscle_13.bam muscle_14.bam muscle_15.bam muscle_16.bam muscle_17.bam muscle_18.bam muscle_19.bam \
--GTF /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf \
--TE /scratch2/alanxu/2019_killi_rna_raw/TE_gtf_for_TETranscripts.gtf \
--sortByPos --project muscle
