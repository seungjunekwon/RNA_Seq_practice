library(DESeq2) # 使うメインパッケージ

# データを読み込む

dat <- read.csv("en_gg6_cm_large_O-modify.csv", header = T, row.names = 1)
info <- read.table("colData.txt", header = T, sep = '\t')

######## DEG（Differentially Expressed Genes; 差次的発現遺伝子）解析 ########

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

# リード数が少ない遺伝子をフィルタリング（除く）
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]


ddsDE <- DESeq(dds)
vsd <- vst(ddsDE, blind=FALSE)

a <- plotPCA(vsd, intgroup=c("condition"))

print(a)


rld <- rlog(ddsDE, blind=FALSE)

plotPCA(rld, intgroup=c("condition"))
