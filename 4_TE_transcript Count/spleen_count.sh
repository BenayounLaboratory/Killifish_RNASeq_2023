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

TEtranscripts -c spleen_1.bam spleen_2.bam spleen_3.bam spleen_4.bam spleen_5.bam spleen_6.bam spleen_7.bam spleen_8.bam spleen_9.bam \
-t spleen_10.bam spleen_11.bam spleen_12.bam spleen_13.bam spleen_14.bam spleen_15.bam spleen_16.bam spleen_17.bam spleen_18.bam spleen_19.bam \
--GTF /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf \
--TE /scratch2/alanxu/2019_killi_rna_raw/TE_gtf_for_TETranscripts.gtf \
--sortByPos --project spleen
