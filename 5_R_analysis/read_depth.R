# set working directory
setwd("/Users/alanx/Downloads/fastqc_raw/multiqc_data")

# load library
library(ggplot2)
library(ggthemes)
library(scales)
library(gtools)

# Read in Multiqc Report
df <- read.table("multiqc_fastqc.txt", sep = "\t", header = TRUE)
df <- df[mixedorder(df$Sample),]
df <- df[seq(1,nrow(df), by = 2), ]
df$tissue <- "brain"
df$tissue[grepl("Heart", df$Sample)] <- "heart"
df$tissue[grepl("muscle", df$Sample)] <- "muscle"
df$tissue[grepl("spleen", df$Sample)] <- "spleen"
df$group <- c(rep("Young Male", 5), rep("Young Female", 5), rep("Old Male", 5), rep("Old Female", 5),
              rep("Young Male", 5), rep("Young Female", 5), rep("Old Male", 5), rep("Old Female", 5),
              rep("Young Male", 5), rep("Young Female", 4), rep("Old Male", 5), rep("Old Female", 5),
              rep("Young Male", 5), rep("Young Female", 4), rep("Old Male", 5), rep("Old Female", 5))

ggplot(df, aes(x= tissue, y= Total.Sequences)) +
  geom_boxplot(fill = "azure4") +
  theme_clean() +
  xlab("Tissue") +
  ylab("Total Reads") +
  geom_point(aes(color = factor(df$group, levels = c("Young Male", "Young Female",
                                                     "Old Male", "Old Female"))),
             size = 2.5) +
  scale_color_manual("Groups", values = c("deepskyblue", "deeppink", "deepskyblue4", "deeppink4"))
  

ggsave(paste0(Sys.Date(), "_length.pdf"), width = 15, height = 12)
