# Set working directory
setwd("/Users/alanx/OneDrive/Desktop/2019_killi_rna_data/2019_killi_rna_data_descriptor/count_table")

# Import library
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2) 
library(scales)
library(ggthemes)

# Import DEG results
brain <- read.table('2023-05-26_brain_age.txt', sep = '\t', header = TRUE)
heart <- read.table('2023-05-26_heart_age.txt', sep = '\t', header = TRUE)
muscle <- read.table('2023-05-26_muscle_age.txt', sep = '\t', header = TRUE)
spleen <- read.table('2023-05-26_spleen_age.txt', sep = '\t', header = TRUE)

# Read in the zebrafish homology table
convert.table <- read.csv('2023_2020_killi_zebra.txt', sep = '\t')

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
  
  # Get Zebrafish homolog information symbols based on the BLAST results file using org.Dr.eg.db package
  # Some Ids will fail to map and will be ignored
  merged<-merge(unique(convert.table[,c("locus_tag","DanRer_Symbol")]), tissue, by = "locus_tag") 
  
  # There can be duplicate values because of paralogs, I take average of those for log2FoldChange
  data.unique = aggregate(merged[,"log2FoldChange"], list(merged$DanRer_Symbol), mean)
  colnames(data.unique) = c("zebrafish", "log2FoldChange")
  data.unique <- data.unique[data.unique != "",] # discard if no zebrafish homolog
  
  # generate and sort the gene list based on log2FoldChange in decreasing order. 
  geneList = data.unique[,2]  # gene list for GO 
  names(geneList) = as.character(data.unique[,1]) # with zebrafish gene symbols as names
  geneList = sort(geneList, decreasing = TRUE)
  
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  go.bp.gsea <- gseGO(geneList     = geneList,
                      OrgDb        = org.Dr.eg.db,
                      keyType      = 'SYMBOL',
                      ont          = "BP",
                      nPerm        = 1000,
                      minGSSize    = 25,
                      maxGSSize    = 5000)
  
  # Saving the Results
  write.table(go.bp.gsea@result, file = paste0(t.name[i],"_GSEA_Analysis.txt"), quote = F, sep = "\t")
  
  # save results into object
  go.results[[i]] <- go.bp.gsea@result
  
}
names(go.results) <- c("brain", "heart", "muscle", "spleen")
# Spleen Yielded Empty Results
# We will ignore spleen for the following analysis
go.results <- go.results[c(1,2,3)]


ggplot(go.results[[1]],aes(x=,y=Description,colour=NES,size=q_value))+ 
  theme_clean() + 
  geom_point(shape = 16) +
  scale_color_gradient(low = "deepskyblue4", high = "deeppink4")






my.tissues.go <- vector(length=length(go.results), mode="list")
names(my.tissues.go) <- names(go.results)
my.pathways <- c()

for ( i in 1:length(go.results)) {
  my.GO.terms <- paste(go.results[[i]]$ID, go.results[[i]]$Description)
  go.results[[i]]$GO_Term <- my.GO.terms
  my.pathways <- unique(c(my.pathways,my.GO.terms))
}

my.matrix <- matrix(0,length(my.pathways),length(go.results)) # default: -log10(1) pval == 0 no enrichment

# Enrichment matrix
my.matrix2 <- matrix(0,length(my.pathways),length(go.results)) # initialize with Enrichment = 0 if no enrich

# matrix with record of significance
my.matrix3 <- matrix(0,length(my.pathways),length(go.results)) # to get sigificant pathways

colnames(my.matrix)  <- names(go.results)
colnames(my.matrix2) <- names(go.results)
colnames(my.matrix3) <- names(go.results)
rownames(my.matrix)  <- my.pathways
rownames(my.matrix2) <- my.pathways
rownames(my.matrix3) <- my.pathways

# collect data from clusterprofiler run
for (i in 1:length(my.pathways)) {
  for (j in 1:length(go.results)) { # tissues 
    # determine position
    my.id <- which(go.results[[j]]$GO_Term %in% my.pathways[i])
    
    if(length(my.id) == 1) { # if was present in this tissue
      # extract stats: FDR, NES and presence
      my.matrix[i,j] <- -log10(go.results[[j]]$p.adjust[my.id]) # log(0) is undefined
      my.matrix2[i,j] <- go.results[[j]]$NES[my.id]
      my.matrix3[i,j] <- 1
    }
  }
}


# find pathways significant in all 3 tissue types [18 passed threshold]
my.sigs         <- apply(my.matrix3,1,sum) == 3
sum(my.sigs) # 18
my.res.enrich   <- data.frame(my.matrix2[my.sigs,])
my.pval.enrich  <- data.frame(my.matrix[my.sigs,])

# sort by average change
my.average     <- apply(my.res.enrich,1,mean)
my.sorted      <- sort(my.average,index.return=T,decreasing=T)
my.res.enrich2 <- my.res.enrich[my.sorted$ix,]

my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])
my.res.enrich2$Pathnames <- rownames(my.res.enrich2)

# format for ggplot
my.res.enrich3 <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[1])],rep(names(go.results)[1],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[1]])
colnames(my.res.enrich3) <- c('Pathnames','NES','condition','minusLog10Pval')
for ( h in 2:length(names(go.results))) {
  my.new <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[h])],rep(names(go.results)[h],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[h]])
  colnames(my.new) <- colnames(my.res.enrich3)
  my.res.enrich3 <- rbind(my.res.enrich3, 
                          my.new)
  
}

ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=NES,size=minusLog10Pval))+ 
  theme_clean() + 
  geom_point(shape = 16) +
  labs(x = "Tissue", y = "Gene Set") +
  scale_color_gradient(low = "deepskyblue4", high = "deeppink4")
  


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
  
  # Get Zebrafish homolog information symbols based on the BLAST results file using org.Dr.eg.db package
  # Some Ids will fail to map and will be ignored
  merged<-merge(unique(convert.table[,c("locus_tag","DanRer_Symbol")]), tissue, by = "locus_tag") 
  
  # There can be duplicate values because of paralogs, I take average of those for log2FoldChange
  data.unique = aggregate(merged[,"log2FoldChange"], list(merged$DanRer_Symbol), mean)
  colnames(data.unique) = c("zebrafish", "log2FoldChange")
  data.unique <- data.unique[data.unique != "",] # discard if no zebrafish homolog
  
  # generate and sort the gene list based on log2FoldChange in decreasing order. 
  geneList = data.unique[,2]  # gene list for GO 
  names(geneList) = as.character(data.unique[,1]) # with zebrafish gene symbols as names
  geneList = sort(geneList, decreasing = TRUE)
  
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  go.bp.gsea <- gseGO(geneList     = geneList,
                      OrgDb        = org.Dr.eg.db,
                      keyType      = 'SYMBOL',
                      ont          = "BP",
                      nPerm        = 1000,
                      minGSSize    = 25,
                      maxGSSize    = 5000)
  
  # Saving the Results
  write.table(go.bp.gsea@result, file = paste0(t.name[i],"_sex_GSEA_Analysis.txt"), quote = F, sep = "\t")
  
  # save results into object
  go.results[[i]] <- go.bp.gsea@result
  
}
names(go.results) <- c("brain", "heart", "muscle", "spleen")
# Muscle and Spleen Yielded Empty Results
# We will ignore muscle and spleen for the following analysis
go.results <- go.results[c(1,2)]

my.tissues.go <- vector(length=length(go.results), mode="list")
names(my.tissues.go) <- names(go.results)
my.pathways <- c()

for ( i in 1:length(go.results)) {
  my.GO.terms <- paste(go.results[[i]]$ID, go.results[[i]]$Description)
  go.results[[i]]$GO_Term <- my.GO.terms
  my.pathways <- unique(c(my.pathways,my.GO.terms))
}

my.matrix <- matrix(0,length(my.pathways),length(go.results)) # default: -log10(1) pval == 0 no enrichment

# Enrichment matrix
my.matrix2 <- matrix(0,length(my.pathways),length(go.results)) # initialize with Enrichment = 0 if no enrich

# matrix with record of significance
my.matrix3 <- matrix(0,length(my.pathways),length(go.results)) # to get sigificant pathways

colnames(my.matrix)  <- names(go.results)
colnames(my.matrix2) <- names(go.results)
colnames(my.matrix3) <- names(go.results)
rownames(my.matrix)  <- my.pathways
rownames(my.matrix2) <- my.pathways
rownames(my.matrix3) <- my.pathways

# collect data from clusterprofiler run
for (i in 1:length(my.pathways)) {
  for (j in 1:length(go.results)) { # tissues 
    # determine position
    my.id <- which(go.results[[j]]$GO_Term %in% my.pathways[i])
    
    if(length(my.id) == 1) { # if was present in this tissue
      # extract stats: FDR, NES and presence
      my.matrix[i,j] <- -log10(go.results[[j]]$p.adjust[my.id]) # log(0) is undefined
      my.matrix2[i,j] <- go.results[[j]]$NES[my.id]
      my.matrix3[i,j] <- 1
    }
  }
}


# find pathways significant in all 2 tissue types [18 passed threshold]
my.sigs         <- apply(my.matrix3,1,sum) == 2
sum(my.sigs) # 18
my.res.enrich   <- data.frame(my.matrix2[my.sigs,])
my.pval.enrich  <- data.frame(my.matrix[my.sigs,])

# sort by average change
my.average     <- apply(my.res.enrich,1,mean)
my.sorted      <- sort(my.average,index.return=T,decreasing=T)
my.res.enrich2 <- my.res.enrich[my.sorted$ix,]

my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])
my.res.enrich2$Pathnames <- rownames(my.res.enrich2)

# format for ggplot
my.res.enrich3 <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[1])],rep(names(go.results)[1],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[1]])
colnames(my.res.enrich3) <- c('Pathnames','NES','condition','minusLog10Pval')
for ( h in 2:length(names(go.results))) {
  my.new <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[h])],rep(names(go.results)[h],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[h]])
  colnames(my.new) <- colnames(my.res.enrich3)
  my.res.enrich3 <- rbind(my.res.enrich3, 
                          my.new)
  
}

ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=NES,size=minusLog10Pval))+ 
  theme_clean() + 
  geom_point(shape = 16) +
  labs(x = "Tissue", y = "Gene Set") +
  scale_color_gradient(low = "deepskyblue4", high = "deeppink4")