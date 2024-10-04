#Environment
setwd("H:/My Drive/Maitrise/LABO/R_stats/RNA-seq")
#setwd("C:/Users/lblon/Google Drive/Maitrise/LABO/R_stats/RNA-seq")
library(readr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(NMF)
library(gridExtra)
library(gtable)
library(digest)
library(ggpubr)
library(pheatmap)

#Load data
DESEQ_norm_log_readcounts <- read_csv("DESEQ_norm_log_readcounts.csv")
DESEQ_pval <- read_csv("DESeq2_LRT.csv") 
Liste_molecules_PCA_Johanne_ISG <- read_csv("Liste_molecules_PCA_Johanne_ISG.csv")
Liste_molecules_abondance <- read_csv("Liste_molecules_RNA-seq_abondances.csv")
effector_function <- read_csv("effector_function.csv")
immune_regulation <- read_csv("immune_regulation.csv")
HIV_protection <- read_csv("SHIV_HIV_protection.csv")
#
DESEQ_HESN_HLA <- read_csv("DESEQ_HESNs_vs_Ctrl (1).csv")
DESEQ_HESN_TS <- read_csv("DESEQ_HESNs_vs_TSpVIHn (1).csv")
DESEQ_HESN_VIH <- read_csv("DESEQ_HESNs_vs_TSpVIHp (1).csv")
DESEQ_TS_HLA <- read_csv("DESEQ_TSpVIHn_vs_Ctrl (1).csv")
DESEQ_VIH_HLA <- read_csv("DESEQ_TSpVIHp_vs_Ctrl (1).csv")
DESEQ_VIH_TS <- read_csv("DESEQ_TSpVIHp_vs_TSpVIHn (1).csv")

#Filter DESEQ by molecules lists
DESEQ_norm_log_readcounts = right_join(DESEQ_norm_log_readcounts, DESEQ_pval, by = "EnsemblID")
DESEQ_Johanne_ISG = right_join(DESEQ_norm_log_readcounts,Liste_molecules_PCA_Johanne_ISG,by='Gene')
DESEQ_abondance = right_join(DESEQ_norm_log_readcounts,Liste_molecules_abondance, by='Gene')
DESEQ_effector_function = right_join(DESEQ_norm_log_readcounts, effector_function, by='Gene')
DESEQ_immune_regulation = right_join(DESEQ_norm_log_readcounts, immune_regulation, by='Gene')
DESEQ_HIV_protection = right_join(DESEQ_norm_log_readcounts, HIV_protection, by='Gene')
#
DESEQ_HESN_HLA = right_join(DESEQ_norm_log_readcounts, DESEQ_HESN_HLA, by = "EnsemblID")
colnames(DESEQ_HESN_HLA)[2] = "Gene"
DESEQ_HESN_TS = right_join(DESEQ_norm_log_readcounts, DESEQ_HESN_TS, by = "EnsemblID")
colnames(DESEQ_HESN_TS)[2] = "Gene"
DESEQ_HESN_VIH = right_join(DESEQ_norm_log_readcounts, DESEQ_HESN_VIH, by = "EnsemblID")
colnames(DESEQ_HESN_VIH)[2] = "Gene"
DESEQ_TS_HLA = right_join(DESEQ_norm_log_readcounts, DESEQ_TS_HLA, by = "EnsemblID")
colnames(DESEQ_TS_HLA)[2] = "Gene"
DESEQ_VIH_HLA = right_join(DESEQ_norm_log_readcounts, DESEQ_VIH_HLA, by = "EnsemblID")
colnames(DESEQ_VIH_HLA)[2] = "Gene"
DESEQ_VIH_TS = right_join(DESEQ_norm_log_readcounts, DESEQ_VIH_TS, by = "EnsemblID")
colnames(DESEQ_VIH_TS)[2] = "Gene"
#
DESEQ_HESN_HLA = right_join(DESEQ_HESN_HLA, Liste_molecules_PCA_Johanne_ISG, by = "Gene")
DESEQ_HESN_TS = right_join(DESEQ_HESN_TS, Liste_molecules_PCA_Johanne_ISG, by = "Gene")
DESEQ_HESN_VIH = right_join(DESEQ_HESN_VIH, Liste_molecules_PCA_Johanne_ISG, by = "Gene")
DESEQ_TS_HLA = right_join(DESEQ_TS_HLA, Liste_molecules_PCA_Johanne_ISG, by = "Gene")
DESEQ_VIH_HLA = right_join(DESEQ_VIH_HLA, Liste_molecules_PCA_Johanne_ISG, by = "Gene")
DESEQ_VIH_TS = right_join(DESEQ_VIH_TS, Liste_molecules_PCA_Johanne_ISG, by = "Gene")

#Tableau p-value pour Michel
# Anova
df_mich = data.frame(DESEQ_Johanne_ISG$EnsemblID, DESEQ_Johanne_ISG$Gene, DESEQ_Johanne_ISG$log2FoldChange, DESEQ_Johanne_ISG$pvalue)
colnames(df_mich) = c("EnsemblID","Gene","log2FoldChange","pvalue")
write.csv(df_mich, "pvalues_anova.csv")
# T-test
df_mich_t = data.frame(DESEQ_HESN_HLA$EnsemblID,DESEQ_HESN_HLA$Gene,DESEQ_HESN_HLA$pvalue,DESEQ_HESN_TS$pvalue,DESEQ_HESN_VIH$pvalue,DESEQ_TS_HLA$pvalue,DESEQ_VIH_HLA$pvalue,DESEQ_VIH_TS$pvalue)
colnames(df_mich_t) = c("EnsemblID","Gene","p-value_HESN_vs_HLA","p-value_HESN_vs_TS+VIH-","p-value_HESN_vs_TS+VIH+","p-value_TS+VIH-_vs_HLA","p-value_TS+VIH+_vs_HLA","p-value_TS+VIH+_vs_TS+VIH-")
write.csv(df_mich_t, "pvalues_t-test.csv")
#Johanne + ISG
#Heatmap
df_ISG = DESEQ_Johanne_ISG[,1:15]
genes = df_ISG$Gene
df_ISG$EnsemblID = NULL
df_ISG$Gene = NULL
rownames(df_ISG) = genes
df_ISG$"Non CSWs HIV-" = rowMeans(df_ISG[,11:13])
tmp = data.frame(df_ISG$D619TM24, df_ISG$D714TM27, df_ISG$D752TM12,row.names = genes)
df_ISG$"CSWs HIV-" = rowMeans(tmp)
tmp = data.frame(df_ISG$D704TM30,df_ISG$D715TM24,df_ISG$D722TM45,df_ISG$D749TM24,row.names = genes)
df_ISG$HESN = rowMeans(tmp)
df_ISG$"CSWs HIV+" = rowMeans(df_ISG[,1:3])
df_ISG[,1:13] = NULL
row.names(df_ISG) = genes
df_ISG = df_ISG %>% arrange(factor(rownames(df_ISG), levels = Liste_molecules_PCA_Johanne_ISG$Gene))
df_ISG = as.matrix(df_ISG)
aheatmap(df_ISG, scale = "row", layout = '_*', border_color = "grey", fontsize = 10, cexRow = 0.5, Rowv = NA, Colv = NA, annRow = data.frame("p-value" = DESEQ_Johanne_ISG$pvalue), cexAnn = 1, annColors = "-RdYlBu2", labAnn = "row", annLegend = FALSE, main = "", filename = "heatmap_Johanne_ISG.png", width = 9, height = 11)

#Small heatmaps/molecule types
df_eff = DESEQ_effector_function[,1:15]
df_imm = DESEQ_immune_regulation[,1:15]
df_prot = DESEQ_HIV_protection[,1:15]
genes_eff = df_eff$Gene
genes_imm = df_imm$Gene
genes_prot = df_prot$Gene
df_eff$EnsemblID = NULL
df_imm$EnsemblID = NULL
df_prot$EnsemblID = NULL
df_eff$Gene = NULL
df_imm$Gene = NULL
df_prot$Gene = NULL
rownames(df_eff) = genes_eff
rownames(df_imm) = genes_imm
rownames(df_prot) = genes_prot
df_eff$"Non CSWs HIV-" = rowMeans(df_eff[,11:13])
df_imm$"Non CSWs HIV-" = rowMeans(df_imm[,11:13])
df_prot$"Non CSWs HIV-" = rowMeans(df_prot[,11:13])
tmp_eff = data.frame(df_eff$D619TM24, df_eff$D714TM27, df_eff$D752TM12,row.names = genes_eff)
tmp_imm = data.frame(df_imm$D619TM24, df_imm$D714TM27, df_imm$D752TM12,row.names = genes_imm)
tmp_prot = data.frame(df_prot$D619TM24, df_prot$D714TM27, df_prot$D752TM12,row.names = genes_prot)
df_eff$"CSWs HIV-" = rowMeans(tmp_eff)
df_imm$"CSWs HIV-" = rowMeans(tmp_imm)
df_prot$"CSWs HIV-" = rowMeans(tmp_prot)
tmp_eff = data.frame(df_eff$D704TM30,df_eff$D715TM24,df_eff$D722TM45,df_eff$D749TM24,row.names = genes_eff)
tmp_imm = data.frame(df_imm$D704TM30,df_imm$D715TM24,df_imm$D722TM45,df_imm$D749TM24,row.names = genes_imm)
tmp_prot = data.frame(df_prot$D704TM30,df_prot$D715TM24,df_prot$D722TM45,df_prot$D749TM24,row.names = genes_prot)
df_eff$HESN = rowMeans(tmp_eff)
df_imm$HESN = rowMeans(tmp_imm)
df_prot$HESN = rowMeans(tmp_prot)
df_eff$"CSWs HIV+" = rowMeans(df_eff[,1:3])
df_imm$"CSWs HIV+" = rowMeans(df_imm[,1:3])
df_prot$"CSWs HIV+" = rowMeans(df_prot[,1:3])
df_eff[,1:13] = NULL
df_imm[,1:13] = NULL
df_prot[,1:13] = NULL
rownames(df_eff) = genes_eff
rownames(df_imm) = genes_imm
rownames(df_prot) = genes_prot
df_eff = df_eff %>% arrange(factor(rownames(df_eff), levels = effector_function$Gene))
df_imm = df_imm %>% arrange(factor(rownames(df_imm), levels = immune_regulation$Gene))
df_prot = df_prot %>% arrange(factor(rownames(df_prot), levels = HIV_protection$Gene))
df_eff = as.matrix(df_eff)
df_imm = as.matrix(df_imm)
df_prot = as.matrix(df_prot)
aheatmap(df_eff, scale = "row", legend = F, border_color = "grey", fontsize = 10, Colv = NA, annRow = data.frame("p-value" = DESEQ_effector_function$pvalue), cexAnn = 1, annColors = "-RdYlBu2", labAnn = "row", annLegend = FALSE, filename = "heatmap_effector_functions.png", main = "Effector functions")
aheatmap(df_imm, scale = "row", legend = F, border_color = "grey", fontsize = 10, Colv = NA, annRow = data.frame("p-value" = DESEQ_immune_regulation$pvalue), cexAnn = 1, annColors = "-RdYlBu2", labAnn = "row", annLegend = FALSE, filename = "heatmap_immune_regulation.png", main = "Immune regulation")
aheatmap(df_prot, scale = "row", legend = F, border_color = "grey", fontsize = 10, Colv = NA, annRow = data.frame("p-value" = DESEQ_HIV_protection$pvalue), cexAnn = 1, annColors = "-RdYlBu2", labAnn = "row", annLegend = FALSE, filename = "heatmap_SHIV_HIV_protection.png", main = "SHIV/HIV protection")

#Layout 4 heatmaps
png("layout_heatmap.png", width = 9, height = 11, units = "in", res = 600)
layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = T), widths = c(1,1))
layout.show(4)
aheatmap(df_ISG, scale = "row", layout = '_*', border_color = "grey", fontsize = 9, cexRow = 0.6, cexCol = 0.85, Rowv = NA, Colv = NA)
title(main = "A", adj = 0,line = 2.5)
aheatmap(df_eff, scale = "row", legend = F, border_color = "grey", fontsize = 9, cexRow = 0.8, cexCol = 0.48, Colv = NA)
title(main = "B", adj = 0, line = 2.5)
aheatmap(df_imm, scale = "row", legend = F, border_color = "grey", fontsize = 9, cexRow = 0.8, cexCol = 0.48, Colv = NA)
title(main = "C", adj = 0, line = 2.5)
aheatmap(df_prot, scale = "row", legend = F, border_color = "grey", fontsize = 9, cexRow = 0.8, cexCol = 0.48, Colv = NA)
title(main = "D", adj = 0, line = 2.5)
par(mfrow = c(1,1))
dev.off()
