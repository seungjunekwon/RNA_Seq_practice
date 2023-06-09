# MDS plot
# https://enjoybioinfo.blogspot.com/2021/06/deg-r-count-data-mds-normalization.html

Data <- read.csv("cm.csv", header = T)
CountData <- read.csv("cm.csv", header = T, row.names = 1)
GeneIDs <- Data[, 1]

library(edgeR) # Google "edgeR install"

Condition <- c("Ctrl", "Ctrl", "Ctrl", "Lin", "Lin", "Lin")
SampleID <- colnames(CountData)

MetaData <- data.frame(SampleID, Condition)
targets <- MetaData

targets <- data.frame(SampleID, Condition)
targets$Condition <- factor(targets$Condition)

levels(targets$Condition)
Condition <- factor(Condition, levels = c("Ctrl", "Lin"))
targets$Condition <- factor(targets$Condition, levels = c("Ctrl","Lin"))

design <- model.matrix(~Condition, data = targets)
y <- DGEList(counts = CountData, gene = GeneIDs)
y <- calcNormFactors(y)
y <- estimateGLMRobustDisp(y, design)

MDS_data <- plotMDS(y, top = 20000) # using top 20000 genes
