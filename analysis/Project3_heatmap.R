library(tidyverse)
hmdata <- read.csv(file='normalizedCounts.csv', sep = ',')
names(hmdata)[1] <- 'ID'
genelabel <- hmdata$ID
hmdata$ID <- NULL

#Add Average Count and Coefficient of Variation metric 
hmdata$avgcounts <- rowMeans(hmdata, na.rm = TRUE)
sd <- apply(hmdata, 1, sd)
hmdata$sd <- sd
hmdata <- mutate(hmdata, CV = sd/avgcounts)
hmdata <- subset(hmdata, select = -sd)

#Filter 
hmdata <- hmdata %>%
  filter(CV < 0.2 & avgcounts > 100)

hmdata <- subset(hmdata, select = -avgcounts)
hmdata <- subset(hmdata, select = -CV)


hmmatrix <- as.matrix(hmdata)
cols <- colnames(hmdata)
colcolors <- character(length(cols))
for(i in 1:length(cols)) {
  if(cols[i] == "SRR1177997" | cols[i] == "SRR1177999" | cols[i] == "SRR1178002") {
    colcolors[i] = "Red"
  } 
  else if(cols[i] == "SRR1178014" | cols[i] == "SRR1178021" | cols[i] == "SRR1178047") {
    colcolors[i] = "Blue"
  } 
  else if(cols[i] == "SRR1177963" | cols[i] == "SRR1177964" | cols[i] == "SRR1177965") {
    colcolors[i] = "Green"
  } else {
    colcolors[i] = "Yellow"
  }
}
print(colcolors)
heatmap(hmmatrix, ColSideColors = colcolors, labRow = genelabel, margins=c(13,5))


