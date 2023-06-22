# Setting Working Directory
setwd("/Users/alanx/OneDrive/Desktop/2019_killi_rna_data/2019_killi_rna_data_descriptor/count_table")

# Library Loading 
library(DESeq2)  # Version 1.40.1
library(ggplot2) # Version 3.4.2
library(clusterProfiler) # Version 4.8.1
library(org.Dr.eg.db) # Version 3.17.0
library(AnnotationDbi) # Version 1.62.1
library(gtools) # Version 3.9.4
library(janitor) # Version 2.2.0
library(ggthemes) # Version 4.2.4
library(plyr) # Version 1.8.8
library(reshape2) # Version 1.4.4
library(tidyr) # Version 1.3.0
library(forcats) # Version 1.0.0
library(ggpubr) # Version 0.6.0
library(rstatix) # Version 0.7.2
library(variancePartition) # Version 1.30.0

# Loading files
brain <- as.data.frame(read.delim(file = "brain.cntTable",header = TRUE, sep = "\t"))
rownames(brain) <- brain$gene.TE
brain <- brain[,2:ncol(brain)]
heart <- as.data.frame(read.delim(file = "heart.cntTable",header = TRUE, sep = "\t"))
rownames(heart) <- heart$gene.TE
heart <- heart[,2:ncol(heart)]
muscle <- as.data.frame(read.delim(file = "muscle.cntTable",header = TRUE, sep = "\t"))
rownames(muscle) <- muscle$gene.TE
muscle <- muscle[,2:ncol(muscle)]
spleen <- as.data.frame(read.delim(file = "spleen.cntTable",header = TRUE, sep = "\t"))
rownames(spleen) <- spleen$gene.TE
spleen <- spleen[,2:ncol(spleen)]
convert_table <- as.data.frame(read.delim(file = "2022-08-29_Killi_human_homology_table.txt",header = TRUE,sep = "\t"))

### there's blank space in the convert table?
convert_table <- convert_table[convert_table$Human_Gene_Symbol != "", ]

# Function to order samples based on age/lexicographical order
reorder <- function(mat){
  new_order <- mixedorder(colnames(mat))
  return(mat[,new_order])
}

# Reodering the samples
brain_ordered <- reorder(brain)
heart_ordered <- reorder(heart)
muscle_ordered <- reorder(muscle)
spleen_ordered <- reorder(spleen)

# Creating Column Information for DESeq2 (Age + Sex)
brain_coldata <- data.frame(sample = colnames(brain_ordered), age = c(rep("young",10),rep("old",10)), row.names = "sample")
brain_coldata$age <- as.factor(brain_coldata$age)
brain_coldata <- cbind(brain_coldata, data.frame(sample = colnames(brain_ordered), 
                                                 sex = c(rep("male",5), rep("female",5), rep("male",5), rep("female",5)), row.names = "sample") )
brain_coldata$sex <- as.factor(brain_coldata$sex)

heart_coldata <- data.frame(sample = colnames(heart_ordered), age = c(rep("young",10),rep("old",10)), row.names = "sample")
heart_coldata$age <- as.factor(heart_coldata$age)
heart_coldata <- cbind(heart_coldata, data.frame(sample = colnames(heart_ordered), 
                                                 sex = c(rep("male",5), rep("female",5), rep("male",5), rep("female",5)), row.names = "sample") )
heart_coldata$sex <- as.factor(heart_coldata$sex)

muscle_coldata <- data.frame(sample = colnames(muscle_ordered), age = c(rep("young",9),rep("old",10)), row.names = "sample")
muscle_coldata$age <- as.factor(muscle_coldata$age)
muscle_coldata <- cbind(muscle_coldata, data.frame(sample = colnames(muscle_ordered), 
                                                   sex = c(rep("male",5), rep("female",4), rep("male",5), rep("female",5)), row.names = "sample") )
muscle_coldata$sex <- as.factor(muscle_coldata$sex)

spleen_coldata <- data.frame(sample = colnames(spleen_ordered), age = c(rep("young",9),rep("old",10)), row.names = "sample")
spleen_coldata$age <- as.factor(spleen_coldata$age)
spleen_coldata <- cbind(spleen_coldata, data.frame(sample = colnames(spleen_ordered), 
                                                   sex = c(rep("male",5), rep("female",4), rep("male",5), rep("female",5)), row.names = "sample") )
spleen_coldata$sex <- as.factor(spleen_coldata$sex)

# Plot for % of reads mapped to TE
# TE starts with NotFur
TE.filter <- function(df){
  mask <- grepl("NotFur", rownames(df))
  return(df[mask,])
}
brain.TE <- TE.filter(brain_ordered)
heart.TE <- TE.filter(heart_ordered)
muscle.TE <- TE.filter(muscle_ordered)
spleen.TE <- TE.filter(spleen_ordered)

brain.TE.nreads <- colSums(brain.TE)
heart.TE.nreads <- colSums(heart.TE)
muscle.TE.nreads <- colSums(muscle.TE)
spleen.TE.nreads <- colSums(spleen.TE)

brain.nreads <- colSums(brain_ordered)
heart.nreads <- colSums(heart_ordered)
muscle.nreads <- colSums(muscle_ordered)
spleen.nreads <- colSums(spleen_ordered)

brain.r <- data.frame(brain.TE.nreads / brain.nreads)
heart.r <- data.frame(heart.TE.nreads / heart.nreads)
muscle.r <- data.frame(muscle.TE.nreads / muscle.nreads)
spleen.r <- data.frame(spleen.TE.nreads / spleen.nreads)
colnames(brain.r) <- "ratio"
colnames(heart.r) <- "ratio"
colnames(muscle.r) <- "ratio"
colnames(spleen.r) <- "ratio"
te.r <- rbind(brain.r, heart.r, muscle.r, spleen.r)
te.r$tissue <- "male young brain"
te.r$tissue[6:10] <- "female young brain"
te.r$tissue[11:15] <- "male old brain"
te.r$tissue[16:20] <- "female old brain"
te.r$tissue[21:25] <- "male young heart"
te.r$tissue[26:30] <- "female young heart"
te.r$tissue[31:35] <- "male old heart"
te.r$tissue[36:40] <- "female old heart"
te.r$tissue[41:45] <- "male young muscle"
te.r$tissue[46:49] <- "female young muscle"
te.r$tissue[50:54] <- "male old muscle"
te.r$tissue[55:59] <- "female old muscle"
te.r$tissue[60:64] <- "male young spleen"
te.r$tissue[65:68] <- "female young spleen"
te.r$tissue[69:73] <- "male old spleen"
te.r$tissue[74:78] <- "female old spleen"

te.r$group <- "Young Male"
te.r$group[grepl("male old", te.r$tissue)] <- "Old Male"
te.r$group[grepl("female young", te.r$tissue)] <- "Young Female"
te.r$group[grepl("female old", te.r$tissue)] <- "Old Female"

stat <- compare_means(ratio ~ tissue, data = te.r, method = "wilcox.test", p.adjust.method = "fdr")
stat <- stat[(stat$group1 == "male young brain" & stat$group2 == "male old brain") |
               (stat$group1 == "female young brain" & stat$group2 == "female old brain")  |
               (stat$group1 == "male young heart" & stat$group2 == "male old heart") |
               (stat$group1 == "female young heart" & stat$group2 == "female old heart") |
               (stat$group1 == "male young muscle" & stat$group2 == "male old muscle") |
               (stat$group1 == "female young muscle" & stat$group2 == "female old muscle") |
               (stat$group1 == "male young spleen" & stat$group2 == "male old spleen") |
               (stat$group1 == "female young spleen" & stat$group2 == "female old spleen"), ]
stat$y.position <- 0.75

ggplot(data = te.r, aes(x = factor(tissue, level = c("male young brain", 
                                                     "male old brain",
                                                     "female young brain",
                                                     "female old brain",
                                                     "male young heart",
                                                     "male old heart",
                                                     "female young heart",
                                                     "female old heart",
                                                     "male young muscle",
                                                     "male old muscle",
                                                     "female young muscle",
                                                     "female old muscle",
                                                     "male young spleen",
                                                     "male old spleen",
                                                     "female young spleen",
                                                     "female old spleen")), y = ratio)) +
  geom_boxplot(fill = "azure4") +
  theme_clean() +
  xlab("Tissue") +
  ylab("Ratio of Reads Mapped to TE") +
  geom_point(size = 2.5, aes(color = factor(group,levels = c("Young Male", "Young Female",
                                                            "Old Male", "Old Female")))) +
  scale_color_manual("Groups", values = c("deepskyblue", "deeppink", "deepskyblue4", "deeppink4")) +
  stat_pvalue_manual(stat, label = "p.adj")


ggsave(paste0(Sys.Date(), "_te_ratio.pdf"), width = 20, height = 12)

# DESEQ2 Analysis
dds_brain <- DESeqDataSetFromMatrix(countData = brain_ordered, colData = brain_coldata, design = ~ sex + age)
dds_brain <- estimateSizeFactors(dds_brain)
norm_brain <- counts(dds_brain, normalized=TRUE)
dds_brain <- dds_brain[rowSums(counts(dds_brain)) >= 6,]
dds_brain <- DESeq(dds_brain)
stabilized_brain <- as.data.frame(getVarianceStabilizedData(dds_brain))
result.brain.age <- as.data.frame(results(dds_brain, contrast = c("age","old","young"),
                        format = "DataFrame"))
result.brain.sex <- as.data.frame(results(dds_brain, contrast = c("sex","female","male"),
                                          format = "DataFrame"))

dds_heart <- DESeqDataSetFromMatrix(countData = heart_ordered, colData = heart_coldata, design = ~ sex + age)
dds_heart <- estimateSizeFactors(dds_heart)
norm_heart <- counts(dds_heart, normalized=TRUE)
dds_heart <- dds_heart[rowSums(counts(dds_heart)) >= 6,]
dds_heart <- DESeq(dds_heart)
stabilized_heart <- as.data.frame(getVarianceStabilizedData(dds_heart))
result.heart.age <- as.data.frame(results(dds_heart, contrast = c("age","old","young"),
                                      format = "DataFrame"))
result.heart.sex <- as.data.frame(results(dds_heart, contrast = c("sex","female","male"),
                                          format = "DataFrame"))

dds_muscle <- DESeqDataSetFromMatrix(countData = muscle_ordered, colData = muscle_coldata, design = ~ sex + age)
dds_muscle <- estimateSizeFactors(dds_muscle)
norm_muscle <- counts(dds_muscle, normalized=TRUE)
dds_muscle <- dds_muscle[rowSums(counts(dds_muscle)) >= 6,]
dds_muscle <- DESeq(dds_muscle)
stabilized_muscle <- as.data.frame(getVarianceStabilizedData(dds_muscle))
result.muscle.age <- as.data.frame(results(dds_muscle, contrast = c("age","old","young"),
                                      format = "DataFrame"))
result.muscle.sex <- as.data.frame(results(dds_muscle, contrast = c("sex","female","male"),
                                           format = "DataFrame"))

dds_spleen <- DESeqDataSetFromMatrix(countData = spleen_ordered, colData = spleen_coldata, design = ~ sex + age)
dds_spleen <- estimateSizeFactors(dds_spleen)
norm_spleen <- counts(dds_spleen, normalized=TRUE)
dds_spleen <- dds_spleen[rowSums(counts(dds_spleen)) >= 6,]
dds_spleen <- DESeq(dds_spleen)
stabilized_spleen <- as.data.frame(getVarianceStabilizedData(dds_spleen))
result.spleen.age <- as.data.frame(results(dds_spleen, contrast = c("age","old","young"),
                                       format = "DataFrame"))
result.spleen.sex <- as.data.frame(results(dds_spleen, contrast = c("sex","female","male"),
                                           format = "DataFrame"))

# Creating a correlation matrix with normalized read counts
cor.mat <- cbind(norm_brain, norm_heart, norm_muscle, norm_spleen)
brain.nam <- paste0(brain_coldata$age, " ", brain_coldata$sex)
heart.nam <- paste0(heart_coldata$age, " ", heart_coldata$sex)
muscle.nam <- paste0(muscle_coldata$age, " ", muscle_coldata$sex)
spleen.nam <- paste0(spleen_coldata$age, " ", spleen_coldata$sex)
comb.nam <- c(brain.nam,heart.nam,muscle.nam,spleen.nam)
colnames(cor.mat) <- paste0(comb.nam[], " ", colnames(cor.mat)[])
cor.res <- cor(cor.mat[,],method="spearman")
melt.res <- melt(cor.res)
melt.res$Var1 <- gsub(pattern = "\\..*", replacement = "", x = melt.res$Var1)
melt.res$Var1 <- gsub(pattern = "_", replacement = " ", x = melt.res$Var1)
melt.res$Var1 <- tolower(melt.res$Var1)
melt.res$Var2 <- gsub(pattern = "\\..*", replacement = "", x = melt.res$Var2)
melt.res$Var2 <- gsub(pattern = "_", replacement = " ", x = melt.res$Var2)
melt.res$Var2 <- tolower(melt.res$Var2)
colnames(melt.res) <- c("Var1", "Var2", "Correlation")
order <- unique(melt.res$Var1)
ggplot(data = melt.res, aes(x=Var1, 
                            y=Var2, fill = Correlation)) +
  geom_tile() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits = order) +
  scale_y_discrete(limits = order) +
  scale_fill_gradientn(colors = c("seashell","red4"),limits = c(0.5, 1.0),
                       breaks = c(0.5,0.6,0.7,0.8,0.9,1.0)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
ggsave(paste0(Sys.Date(),"_cor.pdf"), width = 15, height = 15)

# Exporting the results
write.table(result.brain.age, file = paste0(Sys.Date(), "_brain_age.txt"), sep = "\t", quote = FALSE)
write.table(result.brain.sex, file = paste0(Sys.Date(), "_brain_sex.txt"), sep = "\t", quote = FALSE)
write.table(result.heart.age, file = paste0(Sys.Date(), "_heart_age.txt"), sep = "\t", quote = FALSE)
write.table(result.heart.sex, file = paste0(Sys.Date(), "_heart_sex.txt"), sep = "\t", quote = FALSE)
write.table(result.muscle.age, file = paste0(Sys.Date(), "_muscle_age.txt"), sep = "\t", quote = FALSE)
write.table(result.muscle.sex, file = paste0(Sys.Date(), "_muscle_sex.txt"), sep = "\t", quote = FALSE)
write.table(result.spleen.age, file = paste0(Sys.Date(), "_spleen_age.txt"), sep = "\t", quote = FALSE)
write.table(result.spleen.sex, file = paste0(Sys.Date(), "_spleen_sex.txt"), sep = "\t", quote = FALSE)

# PCA Analysis
make_PCA <- function(count.mat) {
  my.pos.var <- apply(count.mat,1,var) > 0
  my.pca <- prcomp(t(count.mat[my.pos.var,]),scale = TRUE)
  x <- my.pca$x[,1]
  y <- my.pca$x[,2]
  my.summary <- summary(my.pca)
  return(list(x, y, my.summary))
}

# define plotting function
plot_function <- function(x,y,z, pt.shp) {
  par(mar = c(5, 5, 1, 1))
  plot(x[[1]], x[[2]],
       cex=2, pch = pt.shp, col = y,
       xlab = paste('PC1 (', round(100*x[[3]]$importance[,1][2],1),"%)", sep=""),
       ylab = paste('PC2 (', round(100*x[[3]]$importance[,2][2],1),"%)", sep=""),
       xlim = z,
       ylim = z,
       cex.lab = 1.5,
       cex.axis = 1.5)
  legend("bottomright", legend=c("Young Male", "Young Female", "Old Male", "Old Female"),
         col = c("deepskyblue1", "deeppink1", "deepskyblue4", "deeppink4"), 
         pch = c(19,19,15,15))
}
color.1 <- c(rep("deepskyblue1", 5), rep("deeppink1", 5), rep("deepskyblue4", 5),
             rep("deeppink4", 5))
shape.1 <- c(rep(19, 5), rep(19, 5), rep(15, 5),
             rep(15, 5))
color.2 <- c(rep("deepskyblue1", 5), rep("deeppink1", 4), rep("deepskyblue4", 5),
             rep("deeppink4", 5))
shape.2 <- c(rep(19, 5), rep(19, 4), rep(15, 5),
             rep(15, 5))
pdf(paste0(Sys.Date(),"_brain_pca.pdf"))
plot_function(make_PCA(stabilized_brain), color.1, c(-100,150), shape.1)
dev.off()
pdf(paste0(Sys.Date(), "_heart_pca.pdf"))
plot_function(make_PCA(stabilized_heart), color.1, c(-125,100), shape.1)
dev.off()
pdf(paste0(Sys.Date(), "_muscle_pca.pdf"))
plot_function(make_PCA(stabilized_muscle), color.2, c(-200,100), shape.2)
dev.off()
pdf(paste0(Sys.Date(), "_spleen_pca.pdf"))
plot_function(make_PCA(stabilized_spleen), color.2, c(-175,225), shape.2)
dev.off()

# Calculate Variance Explained by Sex and Age
# In this model, we will treat sex as categorical
# age as categorical as well since there is only 2 groups: Young & Old
# The following segment analyzes canonical genes
brain.var   <- fitExtractVarPartModel(stabilized_brain[!grepl("NotFur", rownames(stabilized_brain)), ], 
                                      ~ (1|age) + (1|sex), brain_coldata)
pdf(paste0(Sys.Date(), "_brain_var.pdf"))
plotVarPart(brain.var)
dev.off()

heart.var   <- fitExtractVarPartModel(stabilized_heart[!grepl("NotFur", rownames(stabilized_heart)), ], 
                                      ~ (1|age) + (1|sex), heart_coldata)
pdf(paste0(Sys.Date(), "_heart_var.pdf"))
plotVarPart(heart.var)
dev.off()

muscle.var   <- fitExtractVarPartModel(stabilized_muscle[!grepl("NotFur", rownames(stabilized_muscle)), ], 
                                      ~ (1|age) + (1|sex), muscle_coldata)
pdf(paste0(Sys.Date(), "_muscle_var.pdf"))
plotVarPart(muscle.var)
dev.off()

spleen.var   <- fitExtractVarPartModel(stabilized_spleen[!grepl("NotFur", rownames(stabilized_spleen)), ], 
                                       ~ (1|age) + (1|sex), spleen_coldata)
pdf(paste0(Sys.Date(), "_spleen_var.pdf"))
plotVarPart(spleen.var)
dev.off()

# The following segment analyzes TEs
brain.te.var   <- fitExtractVarPartModel(stabilized_brain[grepl("NotFur", rownames(stabilized_brain)), ], 
                                      ~ (1|age) + (1|sex), brain_coldata)
pdf(paste0(Sys.Date(), "_brain_te_var.pdf"))
plotVarPart(brain.te.var)
dev.off()

heart.te.var   <- fitExtractVarPartModel(stabilized_heart[grepl("NotFur", rownames(stabilized_heart)), ], 
                                      ~ (1|age) + (1|sex), heart_coldata)
pdf(paste0(Sys.Date(), "_heart_te_var.pdf"))
plotVarPart(heart.te.var)
dev.off()

muscle.te.var   <- fitExtractVarPartModel(stabilized_muscle[grepl("NotFur", rownames(stabilized_muscle)), ], 
                                       ~ (1|age) + (1|sex), muscle_coldata)
pdf(paste0(Sys.Date(), "_muscle_te_var.pdf"))
plotVarPart(muscle.te.var)
dev.off()

spleen.te.var   <- fitExtractVarPartModel(stabilized_spleen[grepl("NotFur", rownames(stabilized_spleen)), ], 
                                       ~ (1|age) + (1|sex), spleen_coldata)
pdf(paste0(Sys.Date(), "_spleen_te_var.pdf"))
plotVarPart(spleen.te.var)
dev.off()

# stripplot for diffrentially expressed canonical genes with age
deseq.res.list <- vector(mode = "list", length = 4)
deseq.res.list[[1]] <- result.brain.age[!grepl("NotFur", rownames(result.brain.age)), ]
deseq.res.list[[2]] <- result.heart.age[!grepl("NotFur", rownames(result.heart.age)), ]
deseq.res.list[[3]] <- result.muscle.age[!grepl("NotFur", rownames(result.muscle.age)), ]
deseq.res.list[[4]] <- result.spleen.age[!grepl("NotFur", rownames(result.spleen.age)), ]

age.results <- lapply(deseq.res.list,function(x) {x[order(x$padj),]})
names(age.results) <- c("Brain Age Diffrentially \nExpressed Canonical Genes", 
                        "Heart Age Diffrentially\nExpressed Canonical Genes",
                        "Muscle Age Diffrentially\nExpressed Canonical Genes",
                        "Spleen Age Diffrentially\nExpressed Canonical Genes")
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

col.scheme <- c("goldenrod4", "goldenrod1")
cols <- list(length = 4)
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 30), n[i]) # grey60
  sig.i <- age.results[[i]]$padj <= 0.05
  log.i <- age.results[[i]]$log2FoldChange > 0
  sig.i[is.na(sig.i)] <- FALSE
  for (j in 1:length(sig.i)){
    if (sig.i[j] & log.i[j]){
        cols[[i]][j] <- col.scheme[1]
    } else if (sig.i[j]){
        cols[[i]][j] <- col.scheme[2]
    }
  }
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(cols[[i]][] != rgb(153, 153, 153, maxColorValue = 255, alpha = 30)), " sig.)", sep = "")
}
names(cols) <- names(age.results)

# Making stripplot of canonical genes
pdf(paste0(Sys.Date(), "_age_canonical.pdf"), width = 20, height = 15)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 7.5),
     ylim = c(-25, 25),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-25, 25, by = 5),
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(2*i - 1, nrow(age.results[[i]])), amount = 0.3),
         y = rev(age.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
}
axis(1,
     at = seq(1,7,by=2),
     tick = FALSE,
     las = 1.5,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7,
     line = 1)
axis(2,
     las = 1,
     at = seq(-25, 25, 5))
box()
dev.off()

# stripplot for diffrentially expressed canonical genes with sex
deseq.res.list <- vector(mode = "list", length = 4)
deseq.res.list[[1]] <- result.brain.sex[!grepl("NotFur", rownames(result.brain.sex)), ]
deseq.res.list[[2]] <- result.heart.sex[!grepl("NotFur", rownames(result.heart.sex)), ]
deseq.res.list[[3]] <- result.muscle.sex[!grepl("NotFur", rownames(result.muscle.sex)), ]
deseq.res.list[[4]] <- result.spleen.sex[!grepl("NotFur", rownames(result.spleen.sex)), ]

sex.results <- lapply(deseq.res.list,function(x) {x[order(x$padj),]})
names(sex.results) <- c("Brain Sex Diffrentially\nExpressed Canonical Genes", 
                        "Heart Sex Diffrentially\nExpressed Canonical Genes",
                        "Muscle Sex Diffrentially\nExpressed Canonical Genes",
                        "Spleen Sex Diffrentially\nExpressed Canonical Genes")
n        <- sapply(sex.results, nrow)
names(n) <- names(sex.results)

col.scheme <- c("deeppink4", "deepskyblue4")
cols <- list(length = 4)
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 30), n[i]) # grey60
  sig.i <- sex.results[[i]]$padj <= 0.05
  log.i <- sex.results[[i]]$log2FoldChange > 0
  sig.i[is.na(sig.i)] <- FALSE
  for (j in 1:length(sig.i)){
    if (sig.i[j] & log.i[j]){
        cols[[i]][j] <- col.scheme[1]
    } else if (sig.i[j]){
        cols[[i]][j] <- col.scheme[2]
    }
  }
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(cols[[i]][] != rgb(153, 153, 153, maxColorValue = 255, alpha = 30)), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

pdf(paste0(Sys.Date(), "_sex_canonical.pdf"), width = 20, height = 15)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 7.5),
     ylim = c(-25, 25),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-25, 25, by = 5),
       lty = "dotted",
       col = "grey")
for(i in 1:length(sex.results)){
  set.seed(1234)
  points(x = jitter(rep(2*i - 1, nrow(sex.results[[i]])), amount = 0.3),
         y = rev(sex.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
}
axis(1,
     at = seq(1,7,by=2),
     tick = FALSE,
     las = 1.5,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7,
     line = 1)
axis(2,
     las = 1,
     at = seq(-25, 25, 5))
box()
dev.off()

# Make stripplot of TE with respect to age
deseq.res.list <- vector(mode = "list", length = 4)
deseq.res.list[[1]] <- result.brain.age[grepl("NotFur", rownames(result.brain.age)), ]
deseq.res.list[[2]] <- result.heart.age[grepl("NotFur", rownames(result.heart.age)), ]
deseq.res.list[[3]] <- result.muscle.age[grepl("NotFur", rownames(result.muscle.age)), ]
deseq.res.list[[4]] <- result.spleen.age[grepl("NotFur", rownames(result.spleen.age)), ]

age.results <- lapply(deseq.res.list,function(x) {x[order(x$padj),]})
names(age.results) <- c("Brain Age Diffrentially \nExpressed TEs",
                        "Heart Age Diffrentially\nExpressed TEs",
                        "Muscle Age Diffrentially\nExpressed TEs",
                        "Spleen Age Diffrentially\nExpressed TEs")
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

col.scheme <- c("goldenrod4", "goldenrod1")
cols <- list(length = 4)
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 30), n[i]) # grey60
  sig.i <- age.results[[i]]$padj <= 0.05
  log.i <- age.results[[i]]$log2FoldChange > 0
  sig.i[is.na(sig.i)] <- FALSE
  for (j in 1:length(sig.i)){
    if (sig.i[j] & log.i[j]){
        cols[[i]][j] <- col.scheme[1]
    } else if (sig.i[j]){
        cols[[i]][j] <- col.scheme[2]
    }
  }
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(cols[[i]][] != rgb(153, 153, 153, maxColorValue = 255, alpha = 30)), " sig.)", sep = "")
}
names(cols) <- names(age.results)

# Making stripplot
pdf(paste0(Sys.Date(), "_age_TE.pdf"), width = 20, height = 15)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 7.5),
     ylim = c(-10, 10),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-10, 10, by = 5),
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(2*i - 1, nrow(age.results[[i]])), amount = 0.3),
         y = rev(age.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
}
axis(1,
     at = seq(1,7,by=2),
     tick = FALSE,
     las = 1.5,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7,
     line = 1)
axis(2,
     las = 1,
     at = seq(-10, 10, 5))
box()
dev.off()

# Make stripplot of TE with respect to sex
deseq.res.list <- vector(mode = "list", length = 4)
deseq.res.list[[1]] <- result.brain.sex[grepl("NotFur", rownames(result.brain.sex)), ]
deseq.res.list[[2]] <- result.heart.sex[grepl("NotFur", rownames(result.heart.sex)), ]
deseq.res.list[[3]] <- result.muscle.sex[grepl("NotFur", rownames(result.muscle.sex)), ]
deseq.res.list[[4]] <- result.spleen.sex[grepl("NotFur", rownames(result.spleen.sex)), ]
sex.results <- lapply(deseq.res.list,function(x) {x[order(x$padj),]})
names(sex.results) <- c("Brain Sex Diffrentially\nExpressed TEs",
                        "Heart Sex Diffrentially\nExpressed TEs",
                        "Muscle Sex Diffrentially\nExpressed TEs",
                        "Spleen Sex Diffrentially\nExpressed TEs")
n        <- sapply(sex.results, nrow)
names(n) <- names(sex.results)

col.scheme <- c("deeppink4", "deepskyblue4")
cols <- list(length = 4)
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 30), n[i]) # grey60
  sig.i <- sex.results[[i]]$padj <= 0.05
  log.i <- sex.results[[i]]$log2FoldChange > 0
  sig.i[is.na(sig.i)] <- FALSE
  for (j in 1:length(sig.i)){
    if (sig.i[j] & log.i[j]){
      cols[[i]][j] <- col.scheme[1]
    } else if (sig.i[j]){
      cols[[i]][j] <- col.scheme[2]
    }
  }
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(cols[[i]][] != rgb(153, 153, 153, maxColorValue = 255, alpha = 30)), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

pdf(paste0(Sys.Date(), "_sex_TE.pdf"), width = 20, height = 15)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 7.5),
     ylim = c(-5, 5),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-5, 5, by = 1),
       lty = "dotted",
       col = "grey")
for(i in 1:length(sex.results)){
  set.seed(1234)
  points(x = jitter(rep(2*i - 1, nrow(sex.results[[i]])), amount = 0.3),
         y = rev(sex.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
}
axis(1,
     at = seq(1,7,by=2),
     tick = FALSE,
     las = 1.5,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7,
     line = 1)
axis(2,
     las = 1,
     at = seq(-5, 5, 1))
box()
dev.off()