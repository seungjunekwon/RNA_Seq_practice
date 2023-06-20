# DeSeq2
# https://www.youtube.com/watch?v=wPzeea1Do18

# パッケージを読み込む（ないものはGoogleに " パッケージ名 install " と調べて出てくるコードをコピペ・実行することでインストール）

library(DESeq2) # 使うメインパッケージ
library(org.Gg.eg.db) # 遺伝子IDに遺伝子名をつけるためのパッケージ; ~galGal6をRef.Genomeで使ったらとりあえずこれはなしで
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
res.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(res.df),keytype = "ENSEMBL", column = "SYMBOL") # ~galGal6をRef.Genomeで使ったらとりあえずこれはなしで
res.dfO <- res.df[order(res.df$padj),]
write.csv(res.dfO, "DESeq.csv")

# 正規化の結果ファイル（NormCounts.csv）を作成
normCounts.df <- as.data.frame(normCounts)
normCounts.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(normCounts.df),keytype = "ENSEMBL", column = "SYMBOL") # ~galGal6をRef.Genomeで使ったらとりあえずこれはなしで
write.csv(normCounts.df, "NormCounts.csv")

#### 以下はgalGal7~をRef.Genomeで使った場合；~galGal6を使った人はvolcano.Rで！

# 基本的なボルケーノプロットを作成（色やサイズなど、詳しい調整の仕方は Google " EnhancedVolcano "）
EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1) # 閾値の設定

# 1つの遺伝子のみを見る

EGR1 <- c("EGR1")

EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 0, FCcutoff = 0, selectLab = EGR1, drawConnectors = TRUE, max.overlaps = Inf)

# 複数の遺伝子をピックアップして見る

pickup <- c("SALL4", "LIN28A", "TRIM71", "IRX3", "FGF10", "HOXD9", "MEIS2", "TBX5", "SALL1",
            "MSX1", "TSHZ2", "ZBTB16", "PRDM16", "ETV4", "TFAP2A", "MSX2", "LHX9", "DUSP6", "AXIN2",
            "ACTC1", "GATA4", "FOXF1", "HAND1", "KRT8", "BMP5", "TAGLN")

EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 0, FCcutoff = 0, selectLab = pickup, drawConnectors = TRUE, max.overlaps = Inf)
