#!/bin/bash
#SBATCH --mail-user=alanxu@usc.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=72:00:00
#SBATCH --mem=60GB
#SBATCH --account=bbenayou_34
#SBATCH --partition=oneweek

# Author: Alan Xu
# This script uses star/2.7.0e to align polished reads

module purge
module load gcc/8.3.0
module load intel/18.0.4
module load star/2.7.0e

WORKDIR=/scratch2/alanxu/2019_killi_rna_raw/spleen
cd $WORKDIR
for f in */
do
    cd $f
    read_1=*_1_trimmed.fq.gz
    read_2=*_2_trimmed.fq.gz
    outname="$(echo $read_1 | sed 's/_1_trimmed.fq.gz//')"
    STAR --genomeDir /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/genome_index --readFilesIn $read_1 $read_2 --readFilesCommand zcat --runThreadN 6 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignEndsProtrude 10 ConcordantPair --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /scratch2/alanxu/2019_killi_rna_raw/STAR_Align/polished_read_align/spleen_reads/$outname
    cd $WORKDIR
done
