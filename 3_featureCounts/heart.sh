#!/bin/bash
#SBATCH --mail-user=alanxu@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=05:00:00
#SBATCH --mem=60GB
#SBATCH --account=bbenayou_34
#SBATCH --partition=epyc-64

# Summarizing at the gene level

WORKDIR=/scratch2/alanxu/2019_killi_rna_raw/bam_sorted/heart
cd $WORKDIR
/home1/alanxu/subread-2.0.4-Linux-x86_64/bin/featureCounts -t gene -a /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/GCA_014300015.1_MPIA_NFZ_2.0_genomic.gtf \
-F GTF -O -T 20 -p --countReadPairs -o /scratch2/alanxu/2019_killi_rna_raw/featurecount/heart_gene.txt \
*.bam
