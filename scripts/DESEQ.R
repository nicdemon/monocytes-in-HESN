#Environment
setwd("D:/Google Drive Laurence/Maitrise/LABO/R_stats/RNA-seq")
#setwd("C:/Users/lblon/Google Drive/Maitrise/LABO/R_stats/RNA-seq")

library(DESeq2)
library(readr)

read_count <- read_csv("summary_star_genes_readcount.stranded.annotated.csv")

read_count = read_count[5:58223,]
read_count = read_count[,1:14]
names = read_count$X1
read_count$X1 = NULL
read_count = as.matrix(read_count)
rownames(read_count) = names

coldata = data.frame(
  condition = c("VIH","VIH","VIH","TS","HESN","TS","HESN","HESN","HESN","TS","HLA","HLA","HLA"),
  row.names = colnames(read_count)
  )

dds <- DESeqDataSetFromMatrix(countData = read_count,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- factor(dds$condition, levels = c("HLA","HESN","TS","VIH"))

dds_normal <- DESeq(dds)
res_normal <- results(dds)

dds_lrt = DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt = results(dds_lrt)

write.csv(as.data.frame(res_lrt), file="DESeq2_LRT.csv")
