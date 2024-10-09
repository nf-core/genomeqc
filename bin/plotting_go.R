#!/usr/bin/Rscript

library(pheatmap)

erefd<-read.table("Go_summary_neg.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_neg.pdf", width=9, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

mat_clean <- df2[, !grepl("Node", colnames(df2))]
pdf("Go_summary_neg_noNode.pdf", width=9, height=5)
pheatmap::pheatmap(log10(head(mat_clean, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

erefd<-read.table("Go_summary_pos.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_pos.pdf", width=9, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()


mat_clean2 <- df2[, !grepl("Node", colnames(df2))]
pdf("Go_summary_pos_noNode.pdf", width=9, height=5)
pheatmap::pheatmap(log10(head(mat_clean2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()


