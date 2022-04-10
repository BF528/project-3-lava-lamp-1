library(ggplot2)

#controls
controls <- as_tibble(read.table('data/control_counts.csv', header = TRUE, sep = ',')[,-1]) 
sub_control_samples <- function(data, subset){
  return(data[subset])
}

controls <- sub_control_samples(controls, c('SRR1178030', 'SRR1178040', 'SRR1178056', 'SRR1178064', 'SRR1178074', 'SRR1178075')) 

#combine controls with treatments
cts_treat_control <- bind_cols(treatments, controls)

#run deseq
# filter out rows that have any zeros for funzies
cts_treat_control <- subset(cts_treat_control,rowSums(cts_treat_control==0)==0)

# sample information
info <- read.csv('data/group_6_rna_info.csv') %>% as_tibble()

info1 <- info[info$vehicle == 'CMC_.5_%',] %>% 
  subset(mode_of_action == 'AhR' | mode_of_action == 'Control')

info2 <- info[info$vehicle == 'CORN_OIL_100_%',] %>%
  subset(mode_of_action == 'CAR/PXR' | mode_of_action == 'Control')

info3 <- info[info$vehicle == 'CMC_.5_%',] %>% 
  subset(mode_of_action == 'PPARA' | mode_of_action == 'Control')

cts_tc_1 <- cts_treat_control %>% select(dput(as.character(info1$Run))) #AhR treatment
cts_tc_2 <- cts_treat_control %>% select(dput(as.character(info2$Run))) #CAR/PXR treatment                                                   
cts_tc_3 <- cts_treat_control %>% select(dput(as.character(info3$Run))) #PPARA treatment

# create the DESeq object with info1 
##CMC_.5_% vehicle
dds1 <- DESeqDataSetFromMatrix(
  countData = cts_tc_1,
  colData = info1,
  design= ~ mode_of_action
)

dds2 <- DESeqDataSetFromMatrix(
  countData = cts_tc_2,
  colData = info2,
  design= ~ mode_of_action
)

dds3 <- DESeqDataSetFromMatrix(
  countData = cts_tc_3,
  colData = info3,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')
dds2$mode_of_action <- relevel(dds2$mode_of_action, ref='Control')
dds3$mode_of_action <- relevel(dds3$mode_of_action, ref='Control')


# run DESeq
dds1 <- DESeq(dds1) #CMC_.5_% vehicle
dds2 <- DESeq(dds2) #CORN_OIL_100_% vehicle                       
dds3 <- DESeq(dds3)
#lfcShrink(dds, 2)

#against Ahr mode of action
res_Ahr <- results(dds1, contrast=c('mode_of_action','AhR','Control'))
res_Ahr <- lfcShrink(dds1, 2)
write.csv(res_Ahr, 'DESeqAhRresults.csv')
write.csv(res_Ahr[order(res_Ahr$padj),], 'AhRresultsSorted.csv') #padj sorted
length(which(res_Ahr$padj < 0.05)) #number of significant genes at threshold
#328
AHRsig <- subset(as_tibble(res_Ahr, rownames = 'genes'), padj < 0.05) %>%
  column_to_rownames('genes')
#Fold change histogram
ggplot(AHRsig, aes(x=log2FoldChange)) + geom_histogram(color='black', fill='blue') +
  ggtitle('AHR Fold Change Values from Significantly Expressed Genes') +
  theme(plot.title = element_text(hjust = 0.5))

#pval scatter plot
ggplot(AHRsig, aes(log2FoldChange, -log10(pvalue))) + geom_point(color = 'blue')+
  ggtitle('AHR Fold Change vs P-Values')+
  theme(plot.title = element_text(hjust = 0.5))

#top ten AHR genes
write.csv(top_n(as_tibble(res_Ahr, rownames = 'genes'), -10, pvalue), 'AHRtopTen.csv')

#against CAR/PXR mode of action
res_CARPXR <- results(dds2, contrast=c('mode_of_action','CAR/PXR','Control'))
res_CARPXR <- lfcShrink(dds2, 2)
write.csv(res_CARPXR, 'DESeqCAR-PXRresults.csv')
write.csv(res_CARPXR[order(res_CARPXR$padj),], 'CAR-PXRresultsSorted.csv') #padj sorted

length(which(res_CARPXR$padj < 0.05))
#3766
CARPXRsig <- subset(as_tibble(res_CARPXR, rownames = 'genes'), padj < 0.05) %>%
  column_to_rownames('genes')

#Top ten CARPXR genes
write.csv(top_n(as_tibble(res_CARPXR, rownames = 'genes'), -10, pvalue), 'CARPXRtopTen.csv')

ggplot(CARPXRsig, aes(x=log2FoldChange)) + geom_histogram(color='black', fill='blue') +
  ggtitle('CARPXR Fold Change Values from Significantly Expressed Genes') +
  theme(plot.title = element_text(hjust = 0.5))

#pval scatter plot
ggplot(CARPXRsig, aes(log2FoldChange, -log10(pvalue))) + geom_point(color = 'blue')+
  ggtitle('CARPXR Fold Change vs P-Values')+
  theme(plot.title = element_text(hjust = 0.5))

#against PPARA mode of action
res_PPARA <- results(dds3, contrast=c('mode_of_action','PPARA','Control'))
res_PPARA <- lfcShrink(dds3, 2)
#res_PPARA <- na.omit(res_PPARA)
write.csv(res_PPARA, 'DESeqPPARAresults.csv')
write.csv(res_PPARA[order(res_PPARA$padj),], 'PPARAresultsSorted.csv')
length(which(res_PPARA$padj < 0.05))
#2814
PPARAsig <- subset(as_tibble(res_PPARA, rownames = 'genes'), padj < 0.05) %>%
  column_to_rownames('genes')

#Top ten PPARA genes
write.csv(top_n(as_tibble(res_PPARA, rownames = 'genes'), -10, pvalue), 'PPARAtopTen.csv')

ggplot(PPARAsig, aes(x=log2FoldChange)) + geom_histogram(color='black', fill='blue') +
  ggtitle('PPARA Fold Change Values from Significantly Expressed Genes') +
  theme(plot.title = element_text(hjust = 0.5))

#pval scatter plot
ggplot(PPARAsig, aes(log2FoldChange, -log10(pvalue))) + geom_point(color = 'blue')+
  ggtitle('AHR Fold Change vs P-Values')+
  theme(plot.title = element_text(hjust = 0.5))

#res_CARPXRDE <- list('res' = res_CARPXR, 'dds' = dds)
#CARPXR_DDS <- estimateSizeFactors(res_CARPXRDE$dds)
#norm_CARPXR_counts <- counts(CARPXR_DDS, normalized=TRUE)
#normalize
#dds<- estimateSizeFactors(cts_treat_control)
#normdds <- counts(dds, normalized=TRUE)
#write.csv(normalizedDESeq, 'normalizedCounts.csv')

write.csv(counts(dds1,normalized=TRUE),'AhR_deseq_norm_counts.csv')
write.csv(counts(dds2,normalized=TRUE),'CARPXR_deseq_norm_counts.csv')
write.csv(counts(dds3,normalized=TRUE),'PPARA_deseq_norm_counts.csv')
