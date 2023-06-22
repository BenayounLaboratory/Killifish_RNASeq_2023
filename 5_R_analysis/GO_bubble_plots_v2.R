#Script to take GSEA results and create GO bubble plots

#set workding directory
setwd("/Users/bryanteefy/Dropbox/2023_Killifish_atlas_RNAseq_Data_Descriptor/Figures")

#import libraries
library(dplyr)
library(ggplot2)
library(grid)


#1.Import all GSEA results
#age down
brain_age_down <- read.csv("/Users/bryanteefy/Downloads/brain_age_down.csv")
heart_age_down <- read.csv("/Users/bryanteefy/Downloads/heart_age_down.csv")
muscle_age_down <- read.csv("/Users/bryanteefy/Downloads/muscle_age_down.csv")
spleen_age_down <- read.csv("/Users/bryanteefy/Downloads/spleen_age_down.csv")
#age up
brain_age_up <- read.csv("/Users/bryanteefy/Downloads/brain_age_up.csv")
heart_age_up <- read.csv("/Users/bryanteefy/Downloads/heart_age_up.csv")
muscle_age_up <- read.csv("/Users/bryanteefy/Downloads/muscle_age_up.csv")
spleen_age_up <- read.csv("/Users/bryanteefy/Downloads/spleen_age_up.csv")
#male biased
brain_sex_down <- read.csv("/Users/bryanteefy/Downloads/brain_sex_down.csv")
heart_sex_down <- read.csv("/Users/bryanteefy/Downloads/heart_sex_down.csv")
muscle_sex_down <- read.csv("/Users/bryanteefy/Downloads/muscle_sex_down.csv")
spleen_sex_down <- read.csv("/Users/bryanteefy/Downloads/spleen_sex_down.csv")
#female biased
brain_sex_up <- read.csv("/Users/bryanteefy/Downloads/brain_sex_up.csv")
heart_sex_up <- read.csv("/Users/bryanteefy/Downloads/heart_sex_up.csv")
muscle_sex_up <- read.csv("/Users/bryanteefy/Downloads/muscle_sex_up.csv")
spleen_sex_up <- read.csv("/Users/bryanteefy/Downloads/spleen_sex_up.csv")

#2. Define plotting functions
#plotting by age
bubbles <- function(up, down, low_color, high_color, tissue) {
  all <-  rbind((up %>%                                      
                   arrange(p.adjust) %>% 
                   slice(1:5)), (down %>%                                      
                                   arrange(p.adjust) %>% 
                                   slice(1:5)))
  all$tissue <- tissue
  all$minlog10fdr <- -log10(all$p.adjust)
  
  all$Names <- paste(all$ID, all$Description, sep = "\n")
  
  my.max <- max(all$NES)
  my.min <- min(all$NES)
  
  ggplot(all, aes(x = tissue, y = reorder(Names, NES) ,colour = NES, size = minlog10fdr)) +
    theme_bw() + geom_point(shape = 16) + scale_size_continuous(limits=c(1,8)) +
    scale_color_gradient2(low= low_color, high= high_color, limits = c(-3, 3), oob = scales::squish)
}

# Get the grobs aging
brain_age_plot <- ggplotGrob(bubbles(brain_age_up, brain_age_down, "darkblue", "firebrick4", "brain"))
heart_age_plot <- ggplotGrob(bubbles(heart_age_up, heart_age_down, "darkblue", "firebrick4", "heart"))
muscle_age_plot <- ggplotGrob(bubbles(muscle_age_up, muscle_age_down, "darkblue", "firebrick4", "muscle"))
spleen_age_plot <- ggplotGrob(bubbles(spleen_age_up, spleen_age_down, "darkblue", "firebrick4", "spleen"))

brain_sex_plot <- ggplotGrob(bubbles(brain_sex_up, brain_sex_down, "deepskyblue", "deeppink", "brain"))
heart_sex_plot <- ggplotGrob(bubbles(heart_sex_up, heart_sex_down, "deepskyblue", "deeppink", "heart"))
muscle_sex_plot <- ggplotGrob(bubbles(muscle_sex_up, muscle_sex_down, "deepskyblue", "deeppink", "muscle"))
spleen_sex_plot <- ggplotGrob(bubbles(spleen_sex_up, spleen_sex_down, "deepskyblue", "deeppink", "spleen"))

# Combine the plots
age_plots = cbind(brain_age_plot, heart_age_plot, muscle_age_plot, spleen_age_plot, size = "last")
sex_plots = cbind(brain_sex_plot, heart_sex_plot, muscle_sex_plot, spleen_sex_plot, size = "last")
both <- rbind(age_plots, sex_plots, size = "last")

# Draw it
grid.newpage()
pdf(paste0(Sys.Date(), "both_GO_plots.pdf"), height = 25, width = 18)
grid.draw(both)
dev.off()


#######################
sink(file = paste(Sys.Date(),"_scRNAseq_PseudoBulk_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()