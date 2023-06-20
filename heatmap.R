# https://enjoybioinfo.blogspot.com/search?q=heatmap

library(pheatmap)
dat <- read.csv("dat.csv", header=T, row.names = 1)
pdf("Heatmap.pdf", width =7, height=14, onefile=FALSE) #open blank PDF
pheatmap(dat, fontsize_row = 1)
dev.off() #close PDF
