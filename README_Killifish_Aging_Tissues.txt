###############################
##########  README   ########## 

**Transcriptional profiling of aging tissues from female and male African turquoise killifish**

Alan Xu, Bryan B. Teefy, Ryan J. Lu, Séverine Nozownik, Alexandra M. Tyers, Dario R. Valenzano, & Bérénice A. Benayoun

This code was used to process, quality control and analyse bulk ribosomal RNA-depleted RNA-seq libraries of aging turquoise killifish tissues.

**********************************
*****   Standalone software  *****
R 4.3.0
fastp 0.23.2 
Fastqc 0.11.9
Multiqc 1.15
STAR 2.7.0e
TEtranscripts 2.2.1
featureCounts 2.0.4

**********************************

**********************************
*****   R Package versions   *****
ggpubr 0.6.0
DESeq2 1.40.1
variancePartition 1.30.0
clusterProfiler 4.8.1
org.Hs.eg.db 3.17.0
**********************************


     - 1_Cell_Ranger   : contains the code used to run cellranger from fastq files (RNA quantification)
     - 2_CITEseq_tools : contains the code used to run CITEseq tools from fastq files (HTO quantification)
     - 3_R_analysis    : contains the code to process and annotate the data 
	 - 1a_fastp_trimming      : contains the code used to run fastp trimming
	 - 1b_fastqc_code         : contains the code used to run fastqc quality trimming
	 - 2_star_alignment       : contains the code used to run STAR alignment
	 - 3_featureCounts        : contains the code used to count exonic and intronic reads 
	 - 4_TE_transcript Count  : contains the code used to run TEtranscripts
	 - 5_R_analysis           : contains the code used to run downstream analysis in R
	 - 6_nf_interactive_db    : contains the code used to deploy shiny app of differential expression results (https://alanxu-usc.shinyapps.io/nf_interactive_db/) 
            
The corresponding raw sequencing data has been deposited under PRJNA952180.
