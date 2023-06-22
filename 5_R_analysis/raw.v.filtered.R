setwd("/Users/alanx/Downloads/fastqc_raw/multiqc_data")
library(gtools)

raw <- read.table("multiqc_fastqc.txt", sep = "\t", header = TRUE)
raw <- raw[mixedorder(raw$Sample),]
raw <- raw[seq(1,nrow(raw), by = 2), ]

setwd("/Users/alanx/Downloads/Star/multiqc_data")
star <- read.table("multiqc_star.txt", header = TRUE, sep = "\t")
star <- star[mixedorder(star$Sample),]

output <- data.frame(cbind(raw$Sample, raw$Total.Sequences, star$total_reads, 
                           star$uniquely_mapped, star$uniquely_mapped_percent))
colnames(output) <- c("sample", "raw", "filtered",
                      "uniquely_mapped", "percent_uniquely_mapped")

write.csv(output, "raw_v_filtered.csv")
