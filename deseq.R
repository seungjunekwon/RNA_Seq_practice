# DeSeq2
# https://www.youtube.com/watch?v=wPzeea1Do18

# パッケージを読み込む（ないものはGoogleに " パッケージ名 install " と調べて出てくるコードをコピペ・実行することでインストール）

library(DESeq2) # 使うメインパッケージ
library(org.Gg.eg.db) # 遺伝子IDに遺伝子名をつけるためのパッケージ
library(EnhancedVolcano) # ボルケーノプロットを作るパッケージ

# データを読み込む

dat <- read.csv("cm.csv", header = T, row.names = 1)
info <- read.table("colData.txt", header = T, sep = '\t')

######## DEG（Differentially Expressed Genes; 差次的発現遺伝子）解析 ########

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

# リード数が少ない遺伝子をフィルタリング（除く）
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

ddsDE <- DESeq(dds)

# リード数を正規化
normCounts <- counts(ddsDE, normalized = T)

# DEG解析の結果ファイル（DESeq.csv）を作成
res <- results(ddsDE, alpha = 0.05)
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(res.df),keytype = "ENSEMBL", column = "SYMBOL")
res.dfO <- res.df[order(res.df$padj),]
write.csv(res.dfO, "DESeq.csv")

# 正規化の結果ファイル（NormCounts.csv）を作成
normCounts.df <- as.data.frame(normCounts)
normCounts.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(normCounts.df),keytype = "ENSEMBL", column = "SYMBOL")
write.csv(normCounts.df, "NormCounts.csv")

# 基本的なボルケーノプロットを作成（色やサイズなど、詳しい調整の仕方は Google " EnhancedVolcano "）
EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1)
