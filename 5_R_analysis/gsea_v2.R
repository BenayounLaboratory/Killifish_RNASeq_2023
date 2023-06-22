# Set working directory
setwd("/Users/alanx/OneDrive/Desktop/2019_killi_rna_data/2019_killi_rna_data_descriptor/count_table")

# Import library
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2) 
library(scales)
library(ggthemes)

# Import DEG results
brain <- read.table('2023-06-01_brain_age.txt', sep = '\t', header = TRUE)
heart <- read.table('2023-06-01_heart_age.txt', sep = '\t', header = TRUE)
muscle <- read.table('2023-06-01_muscle_age.txt', sep = '\t', header = TRUE)
spleen <- read.table('2023-06-01_spleen_age.txt', sep = '\t', header = TRUE)

# Read in the human homology table
convert.table <- read.table('2022-08-29_Killi_human_homology_table.txt', sep = '\t', header = TRUE)
colnames(convert.table)[2] <- "locus_tag"

# Making a List of Results
deseq.res.list.genes <- list(length = 4)
deseq.res.list.genes[[1]] <- brain

deseq.res.list.genes[[2]] <- heart

deseq.res.list.genes[[3]] <- muscle

deseq.res.list.genes[[4]] <- spleen

go.results <- list(length = 4)
t.name <- c("brain", "heart", "muscle", "spleen")
# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.genes)) {
  
  # Nfur Data
  tissue<- deseq.res.list.genes[[i]]
  tissue$locus_tag <- rownames(tissue)
  
  # Get Human homolog information symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  merged<-merge(unique(convert.table[,c("locus_tag","ensembl_gene_id")]), tissue, by = "locus_tag") 
  
  # There can be duplicate values because of paralogs, I take average of those for log2FoldChange
  data.unique = aggregate(merged[,"log2FoldChange"], list(merged$ensembl_gene_id), mean)
  colnames(data.unique) = c("human", "log2FoldChange")
  data.unique <- data.unique[data.unique != "",] # discard if no human homolog
  
  # generate and sort the gene list based on log2FoldChange in decreasing order. 
  geneList = data.unique[,2]  # gene list for GO 
  names(geneList) = as.character(data.unique[,1]) # with human gene symbols as names
  geneList = sort(geneList, decreasing = TRUE)
  
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  go.bp.gsea <- gseGO(geneList     = geneList,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = 'ENSEMBL',
                      ont          = "BP",
                      minGSSize    = 25,
                      maxGSSize    = 5000)
  # Saving the Results
  write.table(go.bp.gsea@result, file = paste0(t.name[i],"_GSEA_Analysis.txt"), quote = F, sep = "\t")
  
  # save results into object
  go.results[[i]] <- go.bp.gsea@result
  
}
names(go.results) <- c("brain", "heart", "muscle", "spleen")

go.up <- list(length = 4)
go.down <- list(length = 4)

for (i in 1:length(go.results)){
  go.up[[i]] <- go.results[[i]][go.results[[i]]$NES > 0 ,]
  go.up[[i]] <- go.up[[i]][order(go.up[[i]]$NES, decreasing = TRUE), ]
  go.down[[i]] <- go.results[[i]][go.results[[i]]$NES < 0 ,]
  go.down[[i]] <- go.down[[i]][order(go.down[[i]]$NES), ]
  go.up[[i]]$pathway <- paste0(go.up[[i]]$ID, " ", go.up[[i]]$Description)
  go.down[[i]]$pathway <- paste0(go.down[[i]]$ID, " ", go.down[[i]]$Description)
  go.up[[i]]$logq <- -log10(go.up[[i]]$qvalue)
  go.down[[i]]$logq <- -log10(go.down[[i]]$qvalue)
  go.up[[i]]$pathway <- factor(go.up[[i]]$pathway, levels = rev(go.up[[i]]$pathway))
  go.down[[i]]$pathway <- factor(go.down[[i]]$pathway, levels = rev(go.down[[i]]$pathway))
  write.csv(go.up[[i]], paste0(t.name[i], "_age_up.csv"))
  write.csv(go.down[[i]], paste0(t.name[i], "_age_down.csv"))
}

# For Sex Diffrentailly Expressed Genes

# Import DEG results
brain <- read.table('2023-05-26_brain_sex.txt', sep = '\t', header = TRUE)
heart <- read.table('2023-05-26_heart_sex.txt', sep = '\t', header = TRUE)
muscle <- read.table('2023-05-26_muscle_sex.txt', sep = '\t', header = TRUE)
spleen <- read.table('2023-05-26_spleen_sex.txt', sep = '\t', header = TRUE)

# Making a List of Results
deseq.res.list.genes <- list(length = 4)
deseq.res.list.genes[[1]] <- brain

deseq.res.list.genes[[2]] <- heart

deseq.res.list.genes[[3]] <- muscle

deseq.res.list.genes[[4]] <- spleen

go.results <- list(length = 4)
t.name <- c("brain", "heart", "muscle", "spleen")
# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.genes)) {
  
  # Nfur Data
  tissue<- deseq.res.list.genes[[i]]
  tissue$locus_tag <- rownames(tissue)
  
  # Get Human homolog information symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  merged<-merge(unique(convert.table[,c("locus_tag","ensembl_gene_id")]), tissue, by = "locus_tag") 
  
  # There can be duplicate values because of paralogs, I take average of those for log2FoldChange
  data.unique = aggregate(merged[,"log2FoldChange"], list(merged$ensembl_gene_id), mean)
  colnames(data.unique) = c("human", "log2FoldChange")
  data.unique <- data.unique[data.unique != "",] # discard if no human homolog
  
  # generate and sort the gene list based on log2FoldChange in decreasing order. 
  geneList = data.unique[,2]  # gene list for GO 
  names(geneList) = as.character(data.unique[,1]) # with human gene symbols as names
  geneList = sort(geneList, decreasing = TRUE)
  
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  go.bp.gsea <- gseGO(geneList     = geneList,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = 'ENSEMBL',
                      ont          = "BP",
                      minGSSize    = 25,
                      maxGSSize    = 5000)
  # Saving the Results
  write.table(go.bp.gsea@result, file = paste0(t.name[i],"_GSEA_Analysis.txt"), quote = F, sep = "\t")
  
  # save results into object
  go.results[[i]] <- go.bp.gsea@result
  
}
names(go.results) <- c("brain", "heart", "muscle", "spleen")

go.up <- list(length = 4)
go.down <- list(length = 4)

for (i in 1:length(go.results)){
  go.up[[i]] <- go.results[[i]][go.results[[i]]$NES > 0 ,]
  go.up[[i]] <- go.up[[i]][order(go.up[[i]]$NES, decreasing = TRUE), ]
  go.down[[i]] <- go.results[[i]][go.results[[i]]$NES < 0 ,]
  go.down[[i]] <- go.down[[i]][order(go.down[[i]]$NES), ]
  go.up[[i]]$pathway <- paste0(go.up[[i]]$ID, " ", go.up[[i]]$Description)
  go.down[[i]]$pathway <- paste0(go.down[[i]]$ID, " ", go.down[[i]]$Description)
  go.up[[i]]$logq <- -log10(go.up[[i]]$qvalue)
  go.down[[i]]$logq <- -log10(go.down[[i]]$qvalue)
  go.up[[i]]$pathway <- factor(go.up[[i]]$pathway, levels = rev(go.up[[i]]$pathway))
  go.down[[i]]$pathway <- factor(go.down[[i]]$pathway, levels = rev(go.down[[i]]$pathway))
  write.csv(go.up[[i]], paste0(t.name[i], "_sex_up.csv"))
  write.csv(go.down[[i]], paste0(t.name[i], "_sex_down.csv"))
}