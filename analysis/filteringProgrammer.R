install.packages('tidyverse')
library('tidyverse')
library('dplyr')
library('readr')
library('stringr')
library('ggplot2')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library('DESeq2')
BiocManager::install("apeglm")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(version='devel')
BiocManager::install("S4Vectors")
library('S4Vectors')

#for counts only
counts_tsv <- function(filename){
  tsv <- as_tibble(read.table(filename))
  filt_tsv <- tsv[c(1,7)]
  names(filt_tsv) <- c('geneid','counts')
  class(filt_tsv[[2]]) <- 'numeric'
  return(filt_tsv[-1,])
}



#For all 9 counts: gene id and counts
counts7963 <- counts_tsv('data/7963counts.txt')
counts7964 <- counts_tsv('data/7964counts.txt')                     
counts7965 <- counts_tsv('data/7965counts.txt')
counts7997 <- counts_tsv('data/7997counts.txt')
counts7999 <- counts_tsv('data/7999counts.txt')
counts8002 <- counts_tsv('data/8002counts.txt')
counts8014 <- counts_tsv('data/8014counts.txt')
counts8021 <- counts_tsv('data/8021counts.txt')
counts8047 <- counts_tsv('data/8047counts.txt')

#merging all counts together
cts_list <- list(counts7963, counts7964, counts7965, counts7997, counts7999, counts8002, counts8014, counts8021, counts8047)
cts_combined <- cts_list %>% reduce(full_join, by='geneid') 
names(cts_combined) <- (c('gene_id', 'SRR1177963', 'SRR1177964', 'SRR1177965', 'SRR1177997', 'SRR1177999', 'SRR1178002', 'SRR1178014', 'SRR1178021', 'SRR1178047'))
cts_combined <- column_to_rownames(cts_combined, 'gene_id')
write.csv(cts_combined, 'counts_matrix.csv')

treatments <- as_tibble(read.csv('counts_matrix.csv')[,-1]) %>%
  column_to_rownames('gene_id')


par(cex.axis=0.75)
logcounts <- log2(treatments + 1) #add one so there's no log(0), the log better describes this
boxplot(logcounts,
        ylab='Log2(counts)',
        las=2,
        main='Log Count Distributions Per Sample')
abline(h=mean(as.matrix(logcounts)), col = 'blue')


#The distribution of counts in each sample were plotted on a log scale, and looked similar
#across all samples, as we should see.  A sample that has a significantly different distribution
#would require further attention.