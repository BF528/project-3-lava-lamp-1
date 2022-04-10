library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

setwd('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis')
# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_6_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

#get the chemical list
chemical_list<-as.list(names(table(samples$chemical)))
chemical_list[[2]]=NULL

#write csv files
for (i in 1:3) {
  #print(i)
  #Subset samples and rma for one chemical each time
  chemical=chemical_list[[i]]
  temp_samples<-samples[samples$chemical %in% c(chemical,'Control'),] 
  rma.subset <- rma[paste0('X',temp_samples$array_id)]
  # construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(
    ~factor(
      temp_samples$chemical,
      levels=c('Control',chemical)
    )
  )
  colnames(design) <- c('Intercept',chemical)

  # run limma
  fit <- lmFit(rma.subset, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH',sort.by = 'p',p.value=0.05)
  #print(t)
  # write out the results to file
  write.csv(t,paste(chemical_list[[i]],'limma_results.csv',sep='_'))
}

limma1<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/3-METHYLCHOLANTHRENE_limma_results.csv')
limma2<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/FLUCONAZOLE_limma_results.csv')
limma3<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/PIRINIXIC_ACID_limma_results.csv')

#Histogram
h1<-ggplot(limma1)+
  geom_histogram(aes(logFC),color='black', fill='blue')
h1

h2<-ggplot(limma2)+
  geom_histogram(aes(logFC),color='black', fill='blue')
h2

h3<-ggplot(limma3)+
  geom_histogram(aes(logFC),color='black', fill='blue')
h3

s1<- ggplot(limma1) +
  geom_point(aes(x=logFC, y=P.Value),shape=1,size = 1,color='red')
s2<- ggplot(limma2) +
  geom_point(aes(x=logFC, y=P.Value),shape=1,size = 1,color='red')
s3<- ggplot(limma3) +
  geom_point(aes(x=logFC, y=P.Value),shape=1,size = 1,color='red')

#####Concordance#####
id_map<-read.csv('/project/bf528/project_3/refseq_affy_map.csv')
id_map<-na.omit(id_map)
#Read the RNA-SEQ DESEQ results
des.results.1<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqAhRresults.csv')
des.results.1<-des.results.1[des.results.1$pvalue<0.05,]
des.results.2<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqCAR-PXRresults.csv')
des.results.2<-des.results.2[des.results.2$pvalue<0.05,]
des.results.3<-read.csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqPPARAresults.csv')
des.results.3<-des.results.3[des.results.3$pvalue<0.05,]
#change names so can use id_map to merge results
colnames(limma1)[1] <- 'PROBEID'
colnames(limma2)[1] <- 'PROBEID'
colnames(limma3)[1] <- 'PROBEID'
colnames(des.results.1)[1]<-'REFSEQ'
colnames(des.results.2)[1]<-'REFSEQ'
colnames(des.results.3)[1]<-'REFSEQ'


df1 = merge(x=limma1,y=id_map,by="PROBEID")
df1 = merge(x=df1,y=des.results.1,by="REFSEQ")

df2 = merge(x=limma2,y=id_map,by="PROBEID")
df2 = merge(x=df2,y=des.results.2,by="REFSEQ")

df3 = merge(x=limma3,y=id_map,by="PROBEID")
df3 = merge(x=df3,y=des.results.3,by="REFSEQ")

N=54879

n0.1<-dim(df1)[1]
n1.1<-dim(limma1)[1]
n2.1<-dim(des.results.1)[1]

n0.2<-dim(df2)[1]
n1.2<-dim(limma2)[1]
n2.2<-dim(des.results.2)[1]

n0.3<-dim(df3)[1]
n1.3<-dim(limma3)[1]
n2.3<-dim(des.results.3)[1]

#Calculate nx to find the background corrected number of “true” overlapping genes.
nx<-function(N,n0,n1,n2){
  x<-(N*n0-n1*n2)/(n0+N-n1-n2)
  return(x)
}

#USING 25000 for N according to reference
nx.1<-nx(N,n0.1,n1.1,n2.1)
nx.2<-nx(N,n0.2,n1.2,n2.2)
nx.3<-nx(N,n0.3,n1.3,n2.3)

concordance<-function(nx,n1,n2){
  2*(nx)/(n1+n2)
}

c1<- concordance(nx.1,n1.1,n2.1)
c2<- concordance(nx.2,n1.2,n2.2)
c3<- concordance(nx.3,n1.3,n2.3)


#plot the concordance vs the number of DE genes from the microarray analysis
microarray<-data.frame("Concordance"= c(c1,c2,c3), "DEG_Number"=c(n1.1,n1.2,n1.3))
pr<-ggplot(microarray,aes(x=DEG_Number,  y=Concordance))+
  geom_point(shape=1,size = 1,color='red')+
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, aes(label = ..rr.label..))
pr

#plot the concordance vs the number of DE genes from the RNA-seq analysis
rna_seq<-data.frame("Concordance"= c(c1,c2,c3), "DEG_Number"=c(n2.1,n2.2,n2.3))
pm<-ggplot(rna_seq,aes(x=DEG_Number,  y=Concordance))+
  geom_point(shape=1,size = 1,color='red')+
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y = 1, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.9, aes(label = ..rr.label..))
pm


#Above and below median
#chemical 1 deseq
med1.des<-median(des.results.1$baseMean,na.rm = TRUE)
above.des.1<-subset(des.results.1,baseMean>=med1.des)
below.des.1<-subset(des.results.1,baseMean<med1.des)

#chemical1 limma
med1.limma<-median(limma1$AveExpr,na.rm = TRUE)
above.limma1<-subset(limma1,AveExpr>=med1.limma)
below.limma1<-subset(limma1,AveExpr<med1.limma)

#chemical 2 deseq
med2.des<-median(des.results.2$baseMean,na.rm = TRUE)
above.des.2<-subset(des.results.2,baseMean>=med2.des)
below.des.2<-subset(des.results.2,baseMean<med2.des)

#chemical 2 limma
med2.limma<-median(limma2$AveExpr,na.rm = TRUE)
above.limma2<-subset(limma2,AveExpr>=med2.limma)
below.limma2<-subset(limma2,AveExpr<med2.limma)

#chemical 3 deseq
med3.des<-median(des.results.3$baseMean,na.rm = TRUE)
above.des.3<-subset(des.results.3,baseMean>=med3.des)
below.des.3<-subset(des.results.3,baseMean<med3.des)

#chemical 3 limma
med3.limma<-median(limma3$AveExpr,na.rm = TRUE)
above.limma3<-subset(limma3,AveExpr>=med3.limma)
below.limma3<-subset(limma3,AveExpr<med3.limma)

#Caculate concordance for above median groups for all three chemicals

df1.above = merge(x=above.limma1,y=id_map,by="PROBEID")
df1.above = merge(x=df1.above,y=above.des.1,by="REFSEQ")

df2.above = merge(x=above.limma2,y=id_map,by="PROBEID")
df2.above = merge(x=df2.above,y=above.des.2,by="REFSEQ")

df3.above = merge(x=above.limma3,y=id_map,by="PROBEID")
df3.above = merge(x=df3.above,y=above.des.3,by="REFSEQ")

n0.1.above<-dim(df1.above)[1]
n1.1.above<-dim(above.limma1)[1]
n2.1.above<-dim(above.des.1)[1]

n0.2.above<-dim(df2.above)[1]
n1.2.above<-dim(above.limma2)[1]
n2.2.above<-dim(above.des.2)[1]

n0.3.above<-dim(df3.above)[1]
n1.3.above<-dim(above.limma3)[1]
n2.3.above<-dim(above.des.3)[1]

nx.1.above<-nx(N,n0.1.above,n1.1.above,n2.1.above)
nx.2.above<-nx(N,n0.2.above,n1.2.above,n2.2.above)
nx.3.above<-nx(N,n0.3.above,n1.3.above,n2.3.above)

c1.above<- concordance(nx.1.above,n1.1.above,n2.1.above)
c2.above<- concordance(nx.2.above,n1.2.above,n2.2.above)
c3.above<- concordance(nx.3.above,n1.3.above,n2.3.above)

#Caculate concordance for above median groups for all three chemicals

df1.below = merge(x=below.limma1,y=id_map,by="PROBEID")
df1.below = merge(x=df1.below,y=below.des.1,by="REFSEQ")

df2.below = merge(x=below.limma2,y=id_map,by="PROBEID")
df2.below = merge(x=df2.below,y=below.des.2,by="REFSEQ")

df3.below = merge(x=below.limma3,y=id_map,by="PROBEID")
df3.below = merge(x=df3.below,y=below.des.3,by="REFSEQ")

n0.1.below<-dim(df1.below)[1]
n1.1.below<-dim(below.limma1)[1]
n2.1.below<-dim(below.des.1)[1]

n0.2.below<-dim(df2.below)[1]
n1.2.below<-dim(below.limma2)[1]
n2.2.below<-dim(below.des.2)[1]

n0.3.below<-dim(df3.below)[1]
n1.3.below<-dim(below.limma3)[1]
n2.3.below<-dim(below.des.3)[1]

nx.1.below<-nx(N,n0.1.below,n1.1.below,n2.1.below)
nx.2.below<-nx(N,n0.2.below,n1.2.below,n2.2.below)
nx.3.below<-nx(N,n0.3.below,n1.3.below,n2.3.below)

c1.below<- concordance(nx.1.below,n1.1.below,n2.1.below)
c2.below<- concordance(nx.2.below,n1.2.below,n2.2.below)
c3.below<- concordance(nx.3.below,n1.3.below,n2.3.below)

#create concordance table in long format
concordance.table<-data.frame("chemical"=c(rep(chemical_list[[1]],3),rep(chemical_list[[2]],3),rep(chemical_list[[3]],3)),
                              "group"=rep(c('overall','above','below'),3),"concordance"=c(c1,c1.above,c1.below,c2,c2.above,c2.below,c3,c3.above,c3.below))
                              

#bar chart
bc<-ggplot(concordance.table)+
  geom_bar(aes(x=chemical,y=concordance,fill=group),position="dodge", stat = "identity")
bc