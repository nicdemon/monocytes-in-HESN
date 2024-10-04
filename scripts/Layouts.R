# Environment
setwd("../data")
library(readr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(gtable)
library(digest)
library(ggpubr)
library(data.table)
library(tidyverse)
library(reshape2)
library(cowplot)

### Functions ###

# RNA-seq Fonction de graph
moustache_RNA = function(df,i,min = 7.5, max = 15, bars = c(14,13,12,11,10,9)){
  tmp = data.frame()
  tmp[1:4,1] = "HESN"
  tmp[5:7,1] = "Non CSWs HIV-"
  tmp[8:10,1] = "CSWs HIV-"
  tmp[11:13,1] = "CSWs HIV+"
  tmp[,2] = df[,i]
  tmp[1:4,3] = mean(tmp[1:4,2],na.rm=TRUE)
  tmp[5:7,3] = mean(tmp[5:7,2],na.rm=TRUE)
  tmp[8:10,3] = mean(tmp[8:10,2],na.rm=TRUE)
  tmp[11:13,3] = mean(tmp[11:13,2],na.rm=TRUE)
  tmp[1:9,4] = sd(tmp[1:9,2],na.rm=TRUE)
  tmp[5:7,4] = sd(tmp[5:7,2],na.rm=TRUE)
  tmp[8:10,4] = sd(tmp[8:10,2],na.rm=TRUE)
  tmp[11:13,4] = sd(tmp[11:13,2],na.rm=TRUE)
  name = colnames(df)[i]
  colnames(tmp) = c("populations","marqueur","moy","std")
  matrice_tmp = compare_means(marqueur ~ populations, data=tmp, method = "t.test")
  comparison = list()
  for (i in seq(1,length(matrice_tmp$p.signif))){  
    if (matrice_tmp$p.signif[i] != "ns"){
      if (matrice_tmp$group1[i] == "HESN"){
        comparison[[length(comparison)+1]] = c(matrice_tmp$group1[i],matrice_tmp$group2[i])
      }
    }
  }
  tmp %>%
    arrange(marqueur) %>%
    mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x=populations, y=marqueur, group=populations)) + 
    geom_boxplot(aes(x=populations, y=marqueur)) +
    scale_shape_manual(values = c(15:18)) +
    scale_size_area(max_size=2) +
    stat_compare_means(comparisons = comparison, method = "t.test", label = "p.signif", hide.ns = TRUE, na.rm = TRUE, label.y = bars) +
    coord_cartesian(ylim = c(min, max)) +
    xlab("") +
    ylab(paste(name, "expression (log2 readcount)", sep = " ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5), axis.title.y = element_text(size = 10))
}

# FACS Identify / Remove Outliers
outliers = function(df){
  for(i in rownames(df)){
    for(j in colnames(df)){
      if(df[i,j] == 0){
        df[i,j] = NA
      }
    }
  }
  for(i in colnames(df)){
    df[,i] = remove_outliers(df[,i])
  }
  return(df)
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

#FACS fonction de graph
moustache_FACS = function(df,spop,i,min = 0, max = 1000, bars = c(0,0,0,0,0,0)){
  tmp = data.frame()
  tmp[1:9,1] = "HESN"
  tmp[10:20,1] = "CSWs HIV-"
  tmp[21:29,1] = "CSWs HIV+"
  tmp[30:36,1] = "Non CSWs HIV-"
  tmp[,2] = df[,i]
  tmp[1:9,3] = mean(tmp[1:9,2],na.rm=TRUE)
  tmp[10:20,3] = mean(tmp[10:20,2],na.rm=TRUE)
  tmp[21:29,3] = mean(tmp[21:29,2],na.rm=TRUE)
  tmp[30:36,3] = mean(tmp[30:36,2],na.rm=TRUE)
  tmp[1:9,4] = sd(tmp[1:9,2],na.rm=TRUE)
  tmp[10:20,4] = sd(tmp[10:20,2],na.rm=TRUE)
  tmp[21:29,4] = sd(tmp[21:29,2],na.rm=TRUE)
  tmp[30:36,4] = sd(tmp[30:36,2],na.rm=TRUE)
  name = colnames(df)[i]
  colnames(tmp) = c("populations","marqueur","moy","std")
  matrice_tmp = compare_means(marqueur ~ populations, data=tmp, method = "t.test")
  comparison = list()
  for (i in seq(1,length(matrice_tmp$p.signif))){  
    if (matrice_tmp$p.signif[i] != "ns"){
      if (matrice_tmp$group1[i] == "HESN"){
        comparison[[length(comparison)+1]] = c(matrice_tmp$group1[i],matrice_tmp$group2[i])
      }
    }
  }
  tmp %>%
    arrange(marqueur) %>%
    mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x=populations, y=marqueur, group=populations)) + 
    geom_point(aes(x=populations, y=marqueur, shape = populations, size = 1), binaxis = "y", stackdir = "center") +
    scale_shape_manual(values = c(15:18)) +
    scale_size_area(max_size=2) +
    geom_errorbar(aes(ymin=moy-std,ymax=moy+std), width=0.2, position=position_dodge(0.1), size = 0.5) +
    geom_errorbar(aes(ymin=moy,ymax=moy), width=0.2, position=position_dodge(0.1), size = 0.5) +
    stat_compare_means(comparisons = comparison, method = "t.test", label = "p.signif", hide.ns = TRUE, na.rm = TRUE, label.y = bars) +
    ggtitle(as.character(spop)) +
    coord_cartesian(ylim = c(min, max)) +
    xlab("") +
    ylab(paste(name, "(GeoMFI)", sep = " ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, vjust = 0.5), axis.title.y = element_text(size=10))
}

temps_travail = function(df, name, i){
  x = na.omit(df$temps_travail)
  y = df[,i]
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  eq = as.character(as.expression(eq))
  ggplot(df, aes(x=x, y=y)) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x) +
  geom_text(x = Inf, y = Inf, label = eq, parse = T, size = 3, aes(hjust = 1, vjust = 1)) +
  coord_cartesian(xlim = c(0, 20)) +
  xlab("Time of sex work") +
  ylab(paste(name, "(GeoMFI)", sep = " ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size=10))
}

### RNA-seq###

# Load RNA-seq data
DESEQ_norm_log_readcounts <- read_csv("../RNA-seq/DESEQ_norm_log_readcounts.csv")
Liste_molecules_PCA_Johanne_ISG <- read_csv("../RNA-seq/Liste_molecules_PCA_Johanne_ISG.csv")
Liste_molecules_PCA_Johanne <- read_csv("../RNA-seq/Liste_molecules_PCA_Johanne.csv")
Liste_molecules_PCA_Cocktail_Plex <- read_csv("../RNA-seq/Liste_molecules_PCA_Cocktail_Plex.csv")
Liste_molecules_FACS <- read_csv("../RNA-seq/Liste_molecules_PCA_FACS.csv")
Liste_molecules_abondance <- read_csv("../RNA-seq/Liste_molecules_RNA-seq_abondances.csv")


# RNA-seq Filter DESEQ by molecules lists
DESEQ_Johanne_ISG = right_join(DESEQ_norm_log_readcounts,Liste_molecules_PCA_Johanne_ISG,by='Gene')
DESEQ_Johanne = right_join(DESEQ_norm_log_readcounts,Liste_molecules_PCA_Johanne,by='Gene')
DESEQ_Cocktail_Plex = right_join(DESEQ_norm_log_readcounts,Liste_molecules_PCA_Cocktail_Plex,by='Gene')
DESEQ_FACS = right_join(DESEQ_norm_log_readcounts,Liste_molecules_FACS, by='Gene')
DESEQ_abondance = right_join(DESEQ_norm_log_readcounts,Liste_molecules_abondance, by='Gene')

# RNA-seq Join DESeq + liste molécules + log2
DESEQ_abondance = right_join(DESEQ_norm_log_readcounts,Liste_molecules_abondance, by='Gene')
DESEQ_abondance$EnsemblID = NULL
genes = DESEQ_abondance$Gene
DESEQ_abondance$Gene = NULL
rownames(DESEQ_abondance) = genes
DESEQ_abondance = data.frame(t(DESEQ_abondance))
DESEQ_abondance$group = c("VIH","VIH","VIH","TS","HESN","TS","HESN","HESN","HESN","TS","HLA","HLA","HLA")
DESEQ_abondance = with(DESEQ_abondance, DESEQ_abondance[order(group),])
DESEQ_abondance$group = NULL
DESEQ_abondance = data.frame(sapply(DESEQ_abondance, FUN = "log2"))
for (i in seq(1,length(rownames(DESEQ_abondance)))){
  for(j in seq(1,length(colnames(DESEQ_abondance)))){
    if (DESEQ_abondance[i,j] == -Inf){
      DESEQ_abondance[i,j] = 0
    }
  }
}

# RNA-seq Graphs expression log2
RNA_CX3CR1 = moustache_RNA(DESEQ_abondance, 2)#CX3CR1
RNA_NR4A3 = moustache_RNA(DESEQ_abondance, 3)#NR4A3
RNA_CCL3 = moustache_RNA(DESEQ_abondance, 4)#CCL3
RNA_FCGR3A = moustache_RNA(DESEQ_abondance, 6)#FCGR3A
RNA_IL10 = moustache_RNA(DESEQ_abondance, 7)#IL10
RNA_NR4A2 = moustache_RNA(DESEQ_abondance, 8)#NR4A2
RNA_CCR2 = moustache_RNA(DESEQ_abondance, 9)#CCR2
RNA_CD14 = moustache_RNA(DESEQ_abondance, 11)#CD14
RNA_CD83 = moustache_RNA(DESEQ_abondance, 13)#CD83
RNA_HLAG = moustache_RNA(DESEQ_abondance, 14)#HLA-G
RNA_HLADRA = moustache_RNA(DESEQ_abondance, 15)#HLA-DRA
RNA_HLADRB9 = moustache_RNA(DESEQ_abondance, 16)#HLA-DRB9
RNA_HLADRB5 = moustache_RNA(DESEQ_abondance, 17)#HLA-DRB5
RNA_HLADRB6 = moustache_RNA(DESEQ_abondance, 18)#HLA-DRB6
RNA_HLADRB1 = moustache_RNA(DESEQ_abondance, 19)#HLA-DRB1
RNA_AhR = moustache_RNA(DESEQ_abondance, 20)#AhR
RNA_NR4A1 = moustache_RNA(DESEQ_abondance, 21)#NR4A1
RNA_ITGAM = moustache_RNA(DESEQ_abondance, 23)#ITGAM
RNA_ITGAX = moustache_RNA(DESEQ_abondance, 24)#ITGAX
RNA_LILRB2 = moustache_RNA(DESEQ_abondance, 27)#LILRB2
RNA_IFNAR2 = moustache_RNA(DESEQ_abondance, 28)#IFNAR2
RNA_IFNAR1 = moustache_RNA(DESEQ_abondance, 29)#IFNAR1
RNA_TLR7 = moustache_RNA(DESEQ_abondance, 31)#TLR7
RNA_BAFF = moustache_RNA(DESEQ_abondance, 22)#BLys/BAFF
RNA_APRIL = moustache_RNA(DESEQ_abondance, 25)#APRIL
RNA_SEMA4A = moustache_RNA(DESEQ_abondance, 5)#SEMA4A
RNA_SLC36A1 = moustache_RNA(DESEQ_abondance, 12)#SLC36A1
RNA_SERINC5 = moustache_RNA(DESEQ_abondance, 10)#SERINC5
RNA_IL17RA = moustache_RNA(DESEQ_abondance, 30)#IL17RA
RNA_CTSD = moustache_RNA(DESEQ_abondance, 1)#CTSD
RNA_GAA = moustache_RNA(DESEQ_abondance, 26)#GAA

### FACS cellulaire ###

# Load FACS cellulaire data
cell_totaux <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_totaux.csv"))
cell_classiques <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_classiques.csv"))
cell_intermediaires <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_intermediaires.csv"))
cell_non_classiques <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_non-classiques.csv"))
cell_pDC <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_pDC.csv"))
cell_abondances_relatives <- read_csv("../FACS/cellulaire/abondances_relatives.csv")
cell_totaux_temps_travail <- read_csv("../FACS/cellulaire/cellulaire_totaux_temps_travail.csv")
cell_classiques_temps_travail <- read_csv("../FACS/cellulaire/cellulaire_classiques_temps_travail.csv")
cell_intermediaires_temps_travail <- read_csv("../FACS/cellulaire/cellulaire_intermediaires_temps_travail.csv")
cell_non_classiques_temps_travail <- read_csv("../FACS/cellulaire/cellulaire_non-classiques_temps_travail.csv")
cell_pDC_temps_travail <- read_csv("../FACS/cellulaire/cellulaire_pDC_temps_travail.csv")

# FACS cell Rownames
row.names(cell_totaux) = cell_totaux$Participantes
cell_totaux$Participantes = NULL
row.names(cell_classiques) = cell_classiques$Participantes
cell_classiques$Participantes = NULL
row.names(cell_intermediaires) = cell_intermediaires$Participantes
cell_intermediaires$Participantes = NULL
row.names(cell_non_classiques) = cell_non_classiques$Participantes
cell_non_classiques$Participantes = NULL
row.names(cell_pDC) = cell_pDC$Participantes
cell_pDC$Participantes = NULL

# FACS cell Petits df par pops
cell_HESN_tot = as.data.frame(cell_totaux[1:9,2:13])
cell_HESN_classiques = as.data.frame(cell_classiques[1:9,2:13])
cell_HESN_intermediaires = as.data.frame(cell_intermediaires[1:9,2:13])
cell_HESN_non_c = as.data.frame(cell_non_classiques[1:9,2:13])
cell_HESN_pDC = as.data.frame(cell_pDC[1:9,2:11])
row.names(cell_HESN_tot) = row.names(cell_totaux)[1:9]
row.names(cell_HESN_classiques) = row.names(cell_classiques)[1:9]
row.names(cell_HESN_intermediaires) = row.names(cell_intermediaires)[1:9]
row.names(cell_HESN_non_c) = row.names(cell_non_classiques)[1:9]
row.names(cell_HESN_pDC) = row.names(cell_pDC)[1:9]

cell_TS_tot = as.data.frame(cell_totaux[10:20,2:13])
cell_TS_classiques = as.data.frame(cell_classiques[10:20,2:13])
cell_TS_intermediaires = as.data.frame(cell_intermediaires[10:20,2:13])
cell_TS_non_c = as.data.frame(cell_non_classiques[10:20,2:13])
cell_TS_pDC = as.data.frame(cell_pDC[10:20,2:11])
row.names(cell_TS_tot) = row.names(cell_totaux)[10:20]
row.names(cell_TS_classiques) = row.names(cell_classiques)[10:20]
row.names(cell_TS_intermediaires) = row.names(cell_intermediaires)[10:20]
row.names(cell_TS_non_c) = row.names(cell_non_classiques)[10:20]
row.names(cell_TS_pDC) = row.names(cell_pDC)[10:20]

cell_VIH_tot = as.data.frame(cell_totaux[21:29,2:13])
cell_VIH_classiques = as.data.frame(cell_classiques[21:29,2:13])
cell_VIH_intermediaires = as.data.frame(cell_intermediaires[21:29,2:13])
cell_VIH_non_c = as.data.frame(cell_non_classiques[21:29,2:13])
cell_VIH_pDC = as.data.frame(cell_pDC[21:29,2:11])
row.names(cell_VIH_tot) = row.names(cell_totaux)[21:29]
row.names(cell_VIH_classiques) = row.names(cell_classiques)[21:29]
row.names(cell_VIH_intermediaires) = row.names(cell_intermediaires)[21:29]
row.names(cell_VIH_non_c) = row.names(cell_non_classiques)[21:29]
row.names(cell_VIH_pDC) = row.names(cell_pDC)[21:29]

cell_HLA_tot = as.data.frame(cell_totaux[30:36,2:13])
cell_HLA_classiques = as.data.frame(cell_classiques[30:36,2:13])
cell_HLA_intermediaires = as.data.frame(cell_intermediaires[30:36,2:13])
cell_HLA_non_c = as.data.frame(cell_non_classiques[30:36,2:13])
cell_HLA_pDC = as.data.frame(cell_pDC[30:36,2:11])
row.names(cell_HLA_tot) = row.names(cell_totaux)[30:36]
row.names(cell_HLA_classiques) = row.names(cell_classiques)[30:36]
row.names(cell_HLA_intermediaires) = row.names(cell_intermediaires)[30:36]
row.names(cell_HLA_non_c) = row.names(cell_non_classiques)[30:36]
row.names(cell_HLA_pDC) = row.names(cell_pDC)[30:36]

# FACS cell Outliers
cell_HESN_tot = outliers(cell_HESN_tot)
cell_HESN_classiques = outliers(cell_HESN_classiques)
cell_HESN_intermediaires = outliers(cell_HESN_intermediaires)
cell_HESN_non_c = outliers(cell_HESN_non_c)
cell_HESN_pDC = outliers(cell_HESN_pDC)
cell_TS_tot = outliers(cell_TS_tot)
cell_TS_classiques = outliers(cell_TS_classiques)
cell_TS_intermediaires = outliers(cell_TS_intermediaires)
cell_TS_non_c = outliers(cell_TS_non_c)
cell_TS_pDC = outliers(cell_TS_pDC)
cell_VIH_tot = outliers(cell_VIH_tot)
cell_VIH_classiques = outliers(cell_VIH_classiques)
cell_VIH_intermediaires = outliers(cell_VIH_intermediaires)
cell_VIH_non_c = outliers(cell_VIH_non_c)
cell_VIH_pDC = outliers(cell_VIH_pDC)
cell_HLA_tot = outliers(cell_HLA_tot)
cell_HLA_classiques = outliers(cell_HLA_classiques)
cell_HLA_intermediaires = outliers(cell_HLA_intermediaires)
cell_HLA_non_c = outliers(cell_HLA_non_c)
cell_HLA_pDC = outliers(cell_HLA_pDC)

# FACS cell graphs
# Totaux
cell_totaux = rbind(cell_HESN_tot[1:9,],cell_TS_tot[1:11,],cell_VIH_tot[1:9,],cell_HLA_tot[1:7,])
cell_totaux_CD14 = moustache_FACS(cell_totaux,"Total",1,max = 1500) #CD14
cell_totaux_CD16 = moustache_FACS(cell_totaux,"Total",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500)) #CD16
cell_totaux_HLADR = moustache_FACS(cell_totaux,"Total",3,max = 2250, bars = c(2100,2000,0,0,0,0)) #HLA-DR
cell_totaux_HLAG = moustache_FACS(cell_totaux,"Total",4,min = -10,max = 250, bars = c(225,0,0,0,0,0)) #HLA-G
cell_totaux_CX3CR1 = moustache_FACS(cell_totaux,"Total",5,max = 1450, bars = c(1400,1350,1250,1150,0,0)) #CX3CR1
cell_totaux_ILT4 = moustache_FACS(cell_totaux,"Total",6,max = 1500) #ILT4
cell_totaux_CD11c = moustache_FACS(cell_totaux,"Total",7,max = 7250, bars = c(7000,6500,6000,0,0,0)) #CD11c
cell_totaux_CD83 = moustache_FACS(cell_totaux, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0)) #CD83
cell_totaux_CD123 = moustache_FACS(cell_totaux,"Total",9,min = -100,max = 1750) #CD123
cell_totaux_IL10 = moustache_FACS(cell_totaux, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0)) #IL-10
cell_totaux_IFNa = moustache_FACS(cell_totaux,"Total",11,max = 2550, bars = c(2450,2300,2150,0,0,0)) #IFN-a
cell_totaux_CCL3 = moustache_FACS(cell_totaux,"Total",12,min = -50,max = 300) #CCL3
#classiques
cell_classiques = rbind(cell_HESN_classiques[1:9,],cell_TS_classiques[1:11,],cell_VIH_classiques[1:9,],cell_HLA_classiques[1:7,])
cell_classiques_CD14 = moustache_FACS(cell_classiques,"Classical",1,max = 1500)
cell_classiques_CD16 = moustache_FACS(cell_classiques,"Classical",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_classiques_HLADR = moustache_FACS(cell_classiques,"Classical",3,max = 2250, bars = c(2100,2000,0,0,0,0))
cell_classiques_HLAG = moustache_FACS(cell_classiques,"Classical",4,min = -10,max = 250, bars = c(225,0,0,0,0,0))
cell_classiques_CX3CR1 = moustache_FACS(cell_classiques,"Classical",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_classiques_ILT4 = moustache_FACS(cell_classiques,"Classical",6,max = 1500)
cell_classiques_CD11c = moustache_FACS(cell_classiques,"Classical",7,max = 7250, bars = c(7000,6500,6000,0,0,0))
cell_classiques_CD83 = moustache_FACS(cell_classiques, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_classiques_CD123 = moustache_FACS(cell_classiques,"Classical",9,min = -100,max = 1750)
cell_classiques_IL10 = moustache_FACS(cell_classiques, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_classiques_IFNa = moustache_FACS(cell_classiques,"Classical",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
cell_classiques_CCL3 = moustache_FACS(cell_classiques,"Classical",12,min = -50,max = 300)
#intermédiaires
cell_intermediaires = rbind(cell_HESN_intermediaires[1:9,],cell_TS_intermediaires[1:11,],cell_VIH_intermediaires[1:9,],cell_HLA_intermediaires[1:7,])
cell_intermediaires_CD14 = moustache_FACS(cell_intermediaires,"Intermediate",1,max = 1500)
cell_intermediaires_CD16 = moustache_FACS(cell_intermediaires,"Intermediate",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_intermediaires_HLADR = moustache_FACS(cell_intermediaires,"Intermediate",3,max = 2250, bars = c(2100,2000,0,0,0,0))
cell_intermediaires_HLAG = moustache_FACS(cell_intermediaires,"Intermediate",4,min = -10,max = 250, bars = c(225,0,0,0,0,0))
cell_intermediaires_CX3CR1 = moustache_FACS(cell_intermediaires,"Intermediate",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_intermediaires_ILT4 = moustache_FACS(cell_intermediaires,"Intermediate",6,max = 1500)
cell_intermediaires_CD11c = moustache_FACS(cell_intermediaires,"Intermediate",7,max = 7250, bars = c(7000,6500,6000,0,0,0))
cell_intermediaires_CD83 = moustache_FACS(cell_intermediaires, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_intermediaires_CD123 = moustache_FACS(cell_intermediaires,"Intermediate",9,min = -100,max = 1750)
cell_intermediaires_IL10 = moustache_FACS(cell_intermediaires, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_intermediaires_IFNa = moustache_FACS(cell_intermediaires,"Intermediate",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
cell_intermediaires_CCL3 = moustache_FACS(cell_intermediaires,"Intermediate",12,min = -50,max = 300)
#non-classiques
cell_non_classiques = rbind(cell_HESN_non_c[1:9,],cell_TS_non_c[1:11,],cell_VIH_non_c[1:9,],cell_HLA_non_c[1:7,])
cell_non_classiques_CD14 = moustache_FACS(cell_non_classiques,"Non-Classical",1,max = 1500)
cell_non_classiques_CD16 = moustache_FACS(cell_non_classiques,"Non-Classical",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_non_classiques_HLADR = moustache_FACS(cell_non_classiques,"Non-Classical",3,max = 2250, bars = c(2100,2000,0,0,0,0))
cell_non_classiques_HLAG = moustache_FACS(cell_non_classiques,"Non-Classical",4,min = -10,max = 250, bars = c(225,0,0,0,0,0))
cell_non_classiques_CX3CR1 = moustache_FACS(cell_non_classiques,"Non-Classical",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_non_classiques_ILT4 = moustache_FACS(cell_non_classiques,"Non-Classical",6,max = 1500)
cell_non_classiques_CD11c = moustache_FACS(cell_non_classiques,"Non-Classical",7,max = 7250, bars = c(7000,6500,6000,0,0,0))
cell_non_classiques_CD83 = moustache_FACS(cell_non_classiques, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_non_classiques_CD123 = moustache_FACS(cell_non_classiques,"Non-Classical",9,min = -100,max = 1750)
cell_non_classiques_IL10 = moustache_FACS(cell_non_classiques, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_non_classiques_IFNa = moustache_FACS(cell_non_classiques,"Non-Classical",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
cell_non_classiques_CCL3 = moustache_FACS(cell_non_classiques,"Non-Classical",12,min = -50,max = 300)
#pDC
cell_pDC = rbind(cell_HESN_pDC[1:9,],cell_TS_pDC[1:11,],cell_VIH_pDC[1:9,],cell_HLA_pDC[1:7,])
cell_pDC_HLADR = moustache_FACS(cell_pDC, spop = NULL,1,min = 200,max = 650)
cell_pDC_HLAG = moustache_FACS(cell_pDC,"pDC",2,min = 10, max = 25)
cell_pDC_CX3CR1 = moustache_FACS(cell_pDC,"pDC",3,min = 10,max = 35)
cell_pDC_ILT4 = moustache_FACS(cell_pDC,"pDC",4,min = 50,max = 300, bars = c(290,0,0,0,0,0))
cell_pDC_CD11c = moustache_FACS(cell_pDC,"pDC",5,min = -30,max = 10)
cell_pDC_CD83 = moustache_FACS(cell_pDC, spop = NULL,6,min = 10,max = 300)
cell_pDC_CD123 = moustache_FACS(cell_pDC,"pDC",7,max = 5100)
cell_pDC_IL10 = moustache_FACS(cell_pDC,"pDC",8,min = 100,max = 480, bars = c(450,0,0,0,0,0))
cell_pDC_IFNa = moustache_FACS(cell_pDC, spop = NULL,9,max = 550, bars = c(500,470,440,0,0,0))
cell_pDC_CCL3 = moustache_FACS(cell_pDC,"pDC",10,max = 400)

# FACS cell abondances relatives
#Rownames
names = cell_abondances_relatives$participantes
cell_abondances_relatives$participantes = NULL
row.names(cell_abondances_relatives) = names
#Petits df par pops
cell_abond_mean = data.frame(populations = c("HESN","HESN","HESN","HESN","Non CSWs HIV-","Non CSWs HIV-","Non CSWs HIV-","Non CSWs HIV-","CSWs HIV-","CSWs HIV-","CSWs HIV-","CSWs HIV-","CSWs HIV+","CSWs HIV+","CSWs HIV+","CSWs HIV+"), pop_cell  = c("totaux","classiques","intermediaires","non-classiques","totaux","classiques","intermediaires","non-classiques","totaux","classiques","intermediaires","non-classiques","totaux","classiques","intermediaires","non-classiques"))
cell_abond_mean$mean = NA
agg = t(aggregate(cell_abondances_relatives[,2:5], list(cell_abondances_relatives$populations), mean))
cell_abond_mean[1:4,3] = agg[2:5,1]
cell_abond_mean[5:8,3] = agg[2:5,2]
cell_abond_mean[9:12,3] = agg[2:5,3]
cell_abond_mean[13:16,3] = agg[2:5,4]

for (i in seq(1,length(cell_abondances_relatives$populations))){
  if (cell_abondances_relatives[i,1] == "HLA"){
    cell_abondances_relatives[i,1] = "Non CSWs HIV-"
  } else if(cell_abondances_relatives[i,1] == "TS"){
    cell_abondances_relatives[i,1] = "CSWs HIV-"
  } else if (cell_abondances_relatives[i,1] == "VIH"){
    cell_abondances_relatives[i,1] = "CSWs HIV+"
  }
}

cell_df = melt(cell_abondances_relatives[,1:5], id = "populations")

cell_df_pDC = melt(cell_abondances_relatives, id = "populations")[145:180,]

# FACS cell abondances
cell_abondances_graph = cell_df %>%
  mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x = factor(populations), y = value, group = variable, color = variable)) +
    geom_point(aes(x=populations, shape = variable, size = 1), position = position_dodge(width = 0.4)) +
    scale_color_manual(breaks = c('Total','Classical','Intermediate','Non-Classical'), values = c('red3','forestgreen','blue3','orange')) +
    scale_shape_manual(values = c(15:19)) +
    scale_size_area(max_size=2) +
    xlab("") + 
    ylab("Relative percentages (%)") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    guides(size = "none")
    

cell_abondance_pDC = cell_df_pDC %>%
  mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x = factor(populations), y = value, group = variable)) +
    geom_point(aes(x=populations, shape = variable, size = 1), position = position_dodge(width = 0.4)) +
    scale_shape_manual(values = c(16)) +
    scale_size_area(max_size=2) +
    xlab("") +
    ylab("Relative percentages (%)") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    guides(size = "none", shape = "none")

# FACS cell temps travail
# Colnames
cell_totaux_temps_travail = data.frame(cell_totaux_temps_travail)
cell_classiques_temps_travail = data.frame(cell_classiques_temps_travail)
cell_intermediaires_temps_travail = data.frame(cell_intermediaires_temps_travail)
cell_non_classiques_temps_travail = data.frame(cell_non_classiques_temps_travail)
cell_pDC_temps_travail = data.frame(cell_pDC_temps_travail)
colnames(cell_totaux_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CX3CR1","ILT4","CD11c","CD83","CD123","IL10","IFNa","CCL3")
colnames(cell_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CX3CR1","ILT4","CD11c","CD83","CD123","IL10","IFNa","CCL3")
colnames(cell_intermediaires_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CX3CR1","ILT4","CD11c","CD83","CD123","IL10","IFNa","CCL3")
colnames(cell_non_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CX3CR1","ILT4","CD11c","CD83","CD123","IL10","IFNa","CCL3")
colnames(cell_pDC_temps_travail) = c("participantes","temps_travail","HLADR","HLAG","CX3CR1","ILT4","CD11c","CD83","CD123","IL10","IFNa","CCL3")
cell_totaux_temps_travail = outliers(cell_totaux_temps_travail[,2:14])
cell_classiques_temps_travail = outliers(cell_classiques_temps_travail[,2:14])
cell_intermediaires_temps_travail = outliers(cell_intermediaires_temps_travail[,2:14])
cell_non_classiques_temps_travail = outliers(cell_non_classiques_temps_travail[,2:14])
cell_pDC_temps_travail = outliers(cell_pDC_temps_travail[,2:12])

# Totaux
cell_temps_totaux_CD14 = temps_travail(cell_totaux_temps_travail, "CD14", 2)
cell_temps_totaux_CD16 = temps_travail(cell_totaux_temps_travail, "CD16", 3)
cell_temps_totaux_HLADR = temps_travail(cell_totaux_temps_travail, "HLA-DR", 4)
cell_temps_totaux_HLAG = temps_travail(cell_totaux_temps_travail, "HLA-G", 5)
cell_temps_totaux_CX3CR1 = temps_travail(cell_totaux_temps_travail, "CX3CR1", 6)
cell_temps_totaux_ILT4 = temps_travail(cell_totaux_temps_travail, "ILT4", 7)
cell_temps_totaux_CD11c = temps_travail(cell_totaux_temps_travail, "CD11c", 8)
cell_temps_totaux_CD83 = temps_travail(cell_totaux_temps_travail, "CD83", 9)
cell_temps_totaux_CD123 = temps_travail(cell_totaux_temps_travail, "CD123", 10)
cell_temps_totaux_IL10 = temps_travail(cell_totaux_temps_travail, "IL10", 11)
cell_temps_totaux_IFNa = temps_travail(cell_totaux_temps_travail, "IFNa", 12)
cell_temps_totaux_CCL3 = temps_travail(cell_totaux_temps_travail, "CCL3", 13)
# Classiques
cell_temps_classiques_CD14 = temps_travail(cell_classiques_temps_travail, "CD14", 2)
cell_temps_classiques_CD16 = temps_travail(cell_classiques_temps_travail, "CD16", 3)
cell_temps_classiques_HLADR = temps_travail(cell_classiques_temps_travail, "HLA-DR", 4)
cell_temps_classiques_HLAG = temps_travail(cell_classiques_temps_travail, "HLA-G", 5)
cell_temps_classiques_CX3CR1 = temps_travail(cell_classiques_temps_travail, "CX3CR1", 6)
cell_temps_classiques_ILT4 = temps_travail(cell_classiques_temps_travail, "ILT4", 7)
cell_temps_classiques_CD11c = temps_travail(cell_classiques_temps_travail, "CD11c", 8)
cell_temps_classiques_CD83 = temps_travail(cell_classiques_temps_travail, "CD83", 9)
cell_temps_classiques_CD123 = temps_travail(cell_classiques_temps_travail, "CD123", 10)
cell_temps_classiques_IL10 = temps_travail(cell_classiques_temps_travail, "IL10", 11)
cell_temps_classiques_IFNa = temps_travail(cell_classiques_temps_travail, "IFNa", 12)
cell_temps_classiques_CCL3 = temps_travail(cell_classiques_temps_travail, "CCL3", 13)
# Intermédiaires
cell_temps_intermediaires_CD14 = temps_travail(cell_intermediaires_temps_travail, "CD14", 2)
cell_temps_intermediaires_CD16 = temps_travail(cell_intermediaires_temps_travail, "CD16", 3)
cell_temps_intermediaires_HLADR = temps_travail(cell_intermediaires_temps_travail, "HLA-DR", 4)
cell_temps_intermediaires_HLAG = temps_travail(cell_intermediaires_temps_travail, "HLA-G", 5)
cell_temps_intermediaires_CX3CR1 = temps_travail(cell_intermediaires_temps_travail, "CX3CR1", 6)
cell_temps_intermediaires_ILT4 = temps_travail(cell_intermediaires_temps_travail, "ILT4", 7)
cell_temps_intermediaires_CD11c = temps_travail(cell_intermediaires_temps_travail, "CD11c", 8)
cell_temps_intermediaires_CD83 = temps_travail(cell_intermediaires_temps_travail, "CD83", 9)
cell_temps_intermediaires_CD123 = temps_travail(cell_intermediaires_temps_travail, "CD123", 10)
cell_temps_intermediaires_IL10 = temps_travail(cell_intermediaires_temps_travail, "IL10", 11)
cell_temps_intermediaires_IFNa = temps_travail(cell_intermediaires_temps_travail, "IFNa", 12)
cell_temps_intermediaires_CCL3 = temps_travail(cell_intermediaires_temps_travail, "CCL3", 13)
# Non-classiques
cell_temps_non_classiques_CD14 = temps_travail(cell_non_classiques_temps_travail, "CD14", 2)
cell_temps_non_classiques_CD16 = temps_travail(cell_non_classiques_temps_travail, "CD16", 3)
cell_temps_non_classiques_HLADR = temps_travail(cell_non_classiques_temps_travail, "HLA-DR", 4)
cell_temps_non_classiques_HLAG = temps_travail(cell_non_classiques_temps_travail, "HLA-G", 5)
cell_temps_non_classiques_CX3CR1 = temps_travail(cell_non_classiques_temps_travail, "CX3CR1", 6)
cell_temps_non_classiques_ILT4 = temps_travail(cell_non_classiques_temps_travail, "ILT4", 7)
cell_temps_non_classiques_CD11c = temps_travail(cell_non_classiques_temps_travail, "CD11c", 8)
cell_temps_non_classiques_CD83 = temps_travail(cell_non_classiques_temps_travail, "CD83", 9)
cell_temps_non_classiques_CD123 = temps_travail(cell_non_classiques_temps_travail, "CD123", 10)
cell_temps_non_classiques_IL10 = temps_travail(cell_non_classiques_temps_travail, "IL10", 11)
cell_temps_non_classiques_IFNa = temps_travail(cell_non_classiques_temps_travail, "IFNa", 12)
cell_temps_non_classiques_CCL3 = temps_travail(cell_non_classiques_temps_travail, "CCL3", 13)
# pDC
cell_temps_pDC_HLADR = temps_travail(cell_pDC_temps_travail, "HLA-DR", 2)
cell_temps_pDC_HLAG = temps_travail(cell_pDC_temps_travail, "HLA-G", 3)
cell_temps_pDC_CX3CR1 = temps_travail(cell_pDC_temps_travail, "CX3CR1", 4)
cell_temps_pDC_ILT4 = temps_travail(cell_pDC_temps_travail, "ILT4", 5)
cell_temps_pDC_CD11c = temps_travail(cell_pDC_temps_travail, "CD11c", 6)
cell_temps_pDC_CD83 = temps_travail(cell_pDC_temps_travail, "CD83", 7)
cell_temps_pDC_CD123 = temps_travail(cell_pDC_temps_travail, "CD123", 8)
cell_temps_pDC_IL10 = temps_travail(cell_pDC_temps_travail, "IL10", 9)
cell_temps_pDC_IFNa = temps_travail(cell_pDC_temps_travail, "IFNa", 10)
cell_temps_pDC_CCL3 = temps_travail(cell_pDC_temps_travail, "CCL3", 11)

### FACS nucléaire ###

# Load FACS nucléaire data
nuc_totaux <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_totaux.csv"))
nuc_classiques <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_classiques.csv"))
nuc_intermediaires <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_intermediaires.csv"))
nuc_non_classiques <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_non-classiques.csv"))
nuc_totaux_temps_travail <- read_csv("../FACS/nucleaire/nucleaire_totaux_temps_travail.csv")
nuc_classiques_temps_travail <- read_csv("../FACS/nucleaire/nucleaire_classiques_temps_travail.csv")
nuc_intermediaires_temps_travail <- read_csv("../FACS/nucleaire/nucleaire_intermediaires_temps_travail.csv")
nuc_non_classiques_temps_travail <- read_csv("../FACS/nucleaire/nucleaire_non-classiques_temps_travail.csv")

# FACS nuc Rownames
row.names(nuc_totaux) = nuc_totaux$Participantes
nuc_totaux$Participantes = NULL
row.names(nuc_classiques) = nuc_classiques$Participantes
nuc_classiques$Participantes = NULL
row.names(nuc_intermediaires) = nuc_intermediaires$Participantes
nuc_intermediaires$Participantes = NULL
row.names(nuc_non_classiques) = nuc_non_classiques$Participantes
nuc_non_classiques$Participantes = NULL

# FACS nuc Petits df par pops
nuc_HESN_tot = as.data.frame(nuc_totaux[1:9,2:10])
nuc_HESN_classiques = as.data.frame(nuc_classiques[1:9,2:10])
nuc_HESN_intermediaires = as.data.frame(nuc_intermediaires[1:9,2:10])
nuc_HESN_non_c = as.data.frame(nuc_non_classiques[1:9,2:10])
row.names(nuc_HESN_tot) = row.names(nuc_totaux)[1:9]
row.names(nuc_HESN_classiques) = row.names(nuc_classiques)[1:9]
row.names(nuc_HESN_intermediaires) = row.names(nuc_intermediaires)[1:9]
row.names(nuc_HESN_non_c) = row.names(nuc_non_classiques)[1:9]

nuc_TS_tot = as.data.frame(nuc_totaux[10:20,2:10])
nuc_TS_classiques = as.data.frame(nuc_classiques[10:20,2:10])
nuc_TS_intermediaires = as.data.frame(nuc_intermediaires[10:20,2:10])
nuc_TS_non_c = as.data.frame(nuc_non_classiques[10:20,2:10])
row.names(nuc_TS_tot) = row.names(nuc_totaux)[10:20]
row.names(nuc_TS_classiques) = row.names(nuc_classiques)[10:20]
row.names(nuc_TS_intermediaires) = row.names(nuc_intermediaires)[10:20]
row.names(nuc_TS_non_c) = row.names(nuc_non_classiques)[10:20]

nuc_VIH_tot = as.data.frame(nuc_totaux[21:29,2:10])
nuc_VIH_classiques = as.data.frame(nuc_classiques[21:29,2:10])
nuc_VIH_intermediaires = as.data.frame(nuc_intermediaires[21:29,2:10])
nuc_VIH_non_c = as.data.frame(nuc_non_classiques[21:29,2:10])
row.names(nuc_VIH_tot) = row.names(nuc_totaux)[21:29]
row.names(nuc_VIH_classiques) = row.names(nuc_classiques)[21:29]
row.names(nuc_VIH_intermediaires) = row.names(nuc_intermediaires)[21:29]
row.names(nuc_VIH_non_c) = row.names(nuc_non_classiques)[21:29]

nuc_HLA_tot = as.data.frame(nuc_totaux[30:36,2:10])
nuc_HLA_classiques = as.data.frame(nuc_classiques[30:36,2:10])
nuc_HLA_intermediaires = as.data.frame(nuc_intermediaires[30:36,2:10])
nuc_HLA_non_c = as.data.frame(nuc_non_classiques[30:36,2:10])
row.names(nuc_HLA_tot) = row.names(nuc_totaux)[30:36]
row.names(nuc_HLA_classiques) = row.names(nuc_classiques)[30:36]
row.names(nuc_HLA_intermediaires) = row.names(nuc_intermediaires)[30:36]
row.names(nuc_HLA_non_c) = row.names(nuc_non_classiques)[30:36]

# FACS nuc Outliers
nuc_HESN_tot = outliers(nuc_HESN_tot)
nuc_HESN_classiques = outliers(nuc_HESN_classiques)
nuc_HESN_intermediaires = outliers(nuc_HESN_intermediaires)
nuc_HESN_non_c = outliers(nuc_HESN_non_c)
nuc_TS_tot = outliers(nuc_TS_tot)
nuc_TS_classiques = outliers(nuc_TS_classiques)
nuc_TS_intermediaires = outliers(nuc_TS_intermediaires)
nuc_TS_non_c = outliers(nuc_TS_non_c)
nuc_VIH_tot = outliers(nuc_VIH_tot)
nuc_VIH_classiques = outliers(nuc_VIH_classiques)
nuc_VIH_intermediaires = outliers(nuc_VIH_intermediaires)
nuc_VIH_non_c = outliers(nuc_VIH_non_c)
nuc_HLA_tot = outliers(nuc_HLA_tot)
nuc_HLA_classiques = outliers(nuc_HLA_classiques)
nuc_HLA_intermediaires = outliers(nuc_HLA_intermediaires)
nuc_HLA_non_c = outliers(nuc_HLA_non_c)

# FACS nuc Graphs
#Totaux
nuc_totaux = rbind(nuc_HESN_tot[1:9,],nuc_TS_tot[1:11,],nuc_VIH_tot[1:9,],nuc_HLA_tot[1:7,])
nuc_totaux_CD14 = moustache_FACS(nuc_totaux,"Total",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
nuc_totaux_CD16 = moustache_FACS(nuc_totaux,"Total",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
nuc_totaux_HLADR = moustache_FACS(nuc_totaux, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_totaux_HLAG = moustache_FACS(nuc_totaux,"Total",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_totaux_CCR2 = moustache_FACS(nuc_totaux, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_totaux_ILT4 = moustache_FACS(nuc_totaux,"Total",6,max = 1200)
nuc_totaux_CD11b = moustache_FACS(nuc_totaux,"Total",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_totaux_CD11c = moustache_FACS(nuc_totaux, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_totaux_NR4A1 = moustache_FACS(nuc_totaux,"Total",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#classiques
nuc_classiques = rbind(nuc_HESN_classiques[1:9,],nuc_TS_classiques[1:11,],nuc_VIH_classiques[1:9,],nuc_HLA_classiques[1:7,])
nuc_classiques_CD14 = moustache_FACS(nuc_classiques,"Classical",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
nuc_classiques_CD16 = moustache_FACS(nuc_classiques,"Classical",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
nuc_classiques_HLADR = moustache_FACS(nuc_classiques, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_classiques_HLAG = moustache_FACS(nuc_classiques,"Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_classiques_CCR2 = moustache_FACS(nuc_classiques, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_classiques_ILT4 = moustache_FACS(nuc_classiques,"Classical",6,max = 1200)
nuc_classiques_CD11b = moustache_FACS(nuc_classiques,"Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_classiques_CD11c = moustache_FACS(nuc_classiques, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_classiques_NR4A1 = moustache_FACS(nuc_classiques,"Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#intermédiaires
nuc_intermediaires = rbind(nuc_HESN_intermediaires[1:9,],nuc_TS_intermediaires[1:11,],nuc_VIH_intermediaires[1:9,],nuc_HLA_intermediaires[1:7,])
nuc_intermediaires_CD14 = moustache_FACS(nuc_intermediaires,"Intermediate",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
nuc_intermediaires_CD16 = moustache_FACS(nuc_intermediaires,"Intermediate",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
nuc_intermediaires_HLADR = moustache_FACS(nuc_intermediaires, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_intermediaires_HLAG = moustache_FACS(nuc_intermediaires,"Intermediate",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_intermediaires_CCR2 = moustache_FACS(nuc_intermediaires, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_intermediaires_ILT4 = moustache_FACS(nuc_intermediaires,"Intermediate",6,max = 1200)
nuc_intermediaires_CD11b = moustache_FACS(nuc_intermediaires,"Intermediate",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_intermediaires_CD11c = moustache_FACS(nuc_intermediaires, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_intermediaires_NR4A1 = moustache_FACS(nuc_intermediaires,"Intermediate",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#non-classiques
nuc_non_classiques = rbind(nuc_HESN_non_c[1:9,],nuc_TS_non_c[1:11,],nuc_VIH_non_c[1:9,],nuc_HLA_non_c[1:7,])
nuc_non_classiques_CD14 = moustache_FACS(nuc_non_classiques,"Non-Classical",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
nuc_non_classiques_CD16 = moustache_FACS(nuc_non_classiques,"Non-Classical",2,max = 12000, bars = c(11500,11000,10500,0,0,0))
nuc_non_classiques_HLADR = moustache_FACS(nuc_non_classiques, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_non_classiques_HLAG = moustache_FACS(nuc_non_classiques,"Non-Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_non_classiques_CCR2 = moustache_FACS(nuc_non_classiques, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_non_classiques_ILT4 = moustache_FACS(nuc_non_classiques,"Non-Classical",6,max = 1200)
nuc_non_classiques_CD11b = moustache_FACS(nuc_non_classiques,"Non-Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_non_classiques_CD11c = moustache_FACS(nuc_non_classiques, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_non_classiques_NR4A1 = moustache_FACS(nuc_non_classiques,"Non-Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))

# FACS nuc temps travail
#Colnames
nuc_totaux_temps_travail = data.frame(nuc_totaux_temps_travail)
nuc_classiques_temps_travail = data.frame(nuc_classiques_temps_travail)
nuc_intermediaires_temps_travail = data.frame(nuc_intermediaires_temps_travail)
nuc_non_classiques_temps_travail = data.frame(nuc_non_classiques_temps_travail)
colnames(nuc_totaux_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nuc_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nuc_intermediaires_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nuc_non_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
nuc_totaux_temps_travail = outliers(nuc_totaux_temps_travail[,2:11])
nuc_classiques_temps_travail = outliers(nuc_classiques_temps_travail[,2:11])
nuc_intermediaires_temps_travail = outliers(nuc_intermediaires_temps_travail[,2:11])
nuc_non_classiques_temps_travail = outliers(nuc_non_classiques_temps_travail[,2:11])

# Totaux
nuc_temps_totaux_CD14 = temps_travail(nuc_totaux_temps_travail, "CD14", 2)
nuc_temps_totaux_CD16 = temps_travail(nuc_totaux_temps_travail, "CD16", 3)
nuc_temps_totaux_HLADR = temps_travail(nuc_totaux_temps_travail, "HLA-DR", 4)
nuc_temps_totaux_HLAG = temps_travail(nuc_totaux_temps_travail, "HLA-G", 5)
nuc_temps_totaux_CCR2 = temps_travail(nuc_totaux_temps_travail, "CCR2", 6)
nuc_temps_totaux_ILT4 = temps_travail(nuc_totaux_temps_travail, "ILT4", 7)
nuc_temps_totaux_CD11b = temps_travail(nuc_totaux_temps_travail, "CD11B", 8)
nuc_temps_totaux_CD11c = temps_travail(nuc_totaux_temps_travail, "CD11c", 9)
nuc_temps_totaux_NR4A1 = temps_travail(nuc_totaux_temps_travail, "NR4A1", 10)
# Classiques
nuc_temps_classiques_CD14 = temps_travail(nuc_classiques_temps_travail, "CD14", 2)
nuc_temps_classiques_CD16 = temps_travail(nuc_classiques_temps_travail, "CD16", 3)
nuc_temps_classiques_HLADR = temps_travail(nuc_classiques_temps_travail, "HLA-DR", 4)
nuc_temps_classiques_HLAG = temps_travail(nuc_classiques_temps_travail, "HLA-G", 5)
nuc_temps_classiques_CCR2 = temps_travail(nuc_classiques_temps_travail, "CCR2", 6)
nuc_temps_classiques_ILT4 = temps_travail(nuc_classiques_temps_travail, "ILT4", 7)
nuc_temps_classiques_CD11b = temps_travail(nuc_classiques_temps_travail, "CD11b", 8)
nuc_temps_classiques_CD11c = temps_travail(nuc_classiques_temps_travail, "CD11c", 9)
nuc_temps_classiques_NR4A1 = temps_travail(nuc_classiques_temps_travail, "NR4A1", 10)
# Intermédiaires
nuc_temps_intermediaires_CD14 = temps_travail(nuc_intermediaires_temps_travail, "CD14", 2)
nuc_temps_intermediaires_CD16 = temps_travail(nuc_intermediaires_temps_travail, "CD16", 3)
nuc_temps_intermediaires_HLADR = temps_travail(nuc_intermediaires_temps_travail, "HLA-DR", 4)
nuc_temps_intermediaires_HLAG = temps_travail(nuc_intermediaires_temps_travail, "HLA-G", 5)
nuc_temps_intermediaires_CCR2 = temps_travail(nuc_intermediaires_temps_travail, "CCR2", 6)
nuc_temps_intermediaires_ILT4 = temps_travail(nuc_intermediaires_temps_travail, "ILT4", 7)
nuc_temps_intermediaires_CD11b = temps_travail(nuc_intermediaires_temps_travail, "CD11b", 8)
nuc_temps_intermediaires_CD11c = temps_travail(nuc_intermediaires_temps_travail, "CD11c", 9)
nuc_temps_intermediaires_NR4A1 = temps_travail(nuc_intermediaires_temps_travail, "NR4A1", 10)
# Non-classiques
nuc_temps_non_classiques_CD14 = temps_travail(nuc_non_classiques_temps_travail, "CD14", 2)
nuc_temps_non_classiques_CD16 = temps_travail(nuc_non_classiques_temps_travail, "CD16", 3)
nuc_temps_non_classiques_HLADR = temps_travail(nuc_non_classiques_temps_travail, "HLA-DR", 4)
nuc_temps_non_classiques_HLAG = temps_travail(nuc_non_classiques_temps_travail, "HLA-G", 5)
nuc_temps_non_classiques_CCR2 = temps_travail(nuc_non_classiques_temps_travail, "CCR2", 6)
nuc_temps_non_classiques_ILT4 = temps_travail(nuc_non_classiques_temps_travail, "ILT4", 7)
nuc_temps_non_classiques_CD11b = temps_travail(nuc_non_classiques_temps_travail, "CD11b", 8)
nuc_temps_non_classiques_CD11c = temps_travail(nuc_non_classiques_temps_travail, "CD11c", 9)
nuc_temps_non_classiques_NR4A1 = temps_travail(nuc_non_classiques_temps_travail, "NR4A1", 10)

### Layouts article
layout_abondance_NR4A = ggdraw() +
  draw_plot(cell_abondances_graph,0.2,0.66,0.6,0.3) +
  draw_plot(nuc_totaux_NR4A1,0.01,0.31,0.24,0.34) +
  draw_plot(nuc_classiques_NR4A1,0.25,0.31,0.24,0.34) +
  draw_plot(nuc_intermediaires_NR4A1,0.5,0.31,0.24,0.34) +
  draw_plot(nuc_non_classiques_NR4A1,0.75,0.31,0.24,0.34) +
  draw_plot(nuc_temps_totaux_NR4A1,0.01,0,0.24,0.3) +
  draw_plot(nuc_temps_classiques_NR4A1,0.25,0,0.24,0.3) +
  draw_plot(nuc_temps_intermediaires_NR4A1,0.5,0,0.24,0.3) +
  draw_plot(nuc_temps_non_classiques_NR4A1,0.75,0,0.24,0.3) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_plot_label(c("A","B","C","D","E","F","G","H","I"), c(0.19,0.015,0.25,0.5,0.75,0.015,0.26,0.51,0.76), c(0.96,0.65,0.65,0.65,0.65,0.31,0.31,0.31,0.31), size = 15)
ggsave("NR4A_abondance.png", width = 15, height = 15, units = "in", dpi = 600)

layout_CD16_HLADR = plot_grid(
  plot_grid(cell_totaux_CD16,cell_classiques_CD16,cell_intermediaires_CD16,cell_non_classiques_CD16,ncol = 4, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(cell_temps_totaux_CD16,cell_temps_classiques_CD16,cell_temps_intermediaires_CD16,cell_temps_non_classiques_CD16,ncol = 4, nrow = 1, labels = c("E","F","G","H"), align = "h"),
  plot_grid(nuc_totaux_HLADR,nuc_classiques_HLADR,nuc_intermediaires_HLADR,nuc_non_classiques_HLADR,ncol = 4, nrow = 1, labels = c("I","J","K","L"), align = "h"),
  plot_grid(nuc_temps_totaux_HLADR,nuc_temps_classiques_HLADR,nuc_temps_intermediaires_HLADR,nuc_temps_non_classiques_HLADR,ncol = 4, nrow = 1, labels = c("M","N","O","P"), align = "h"),
  nrow = 4, ncol = 1, aligns = "v")
ggsave("CD16_HLADR.png",layout_CD16_HLADR, width = 12, height = 18, units = "in", dpi = 600)

layout_CD11c_CD11b = plot_grid(
  plot_grid(nuc_totaux_CD11b,nuc_classiques_CD11b,nuc_intermediaires_CD11b,nuc_non_classiques_CD11b,ncol = 4, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(nuc_temps_totaux_CD11b,nuc_temps_classiques_CD11b,nuc_temps_intermediaires_CD11b,nuc_temps_non_classiques_CD11b,ncol = 4, nrow = 1, labels = c("E","F","G","H"), align = "h"),
  plot_grid(nuc_totaux_CD11c,nuc_classiques_CD11c,nuc_intermediaires_CD11c,nuc_non_classiques_CD11c,ncol = 4, nrow = 1, labels = c("I","J","K","L"), align = "h"),
  plot_grid(nuc_temps_totaux_CD11c,nuc_temps_classiques_CD11c,nuc_temps_intermediaires_CD11c,nuc_temps_non_classiques_CD11c,ncol = 4, nrow = 1, labels = c("M","N","O","P"), align = "h"),
  nrow = 4, ncol = 1, aligns = "v")
ggsave("CD11c_CD11b.png",layout_CD11c_CD11b, width = 12, height = 18, units = "in", dpi = 600)

layout_HLAG_IL10_AHR = plot_grid(
  plot_grid(nuc_totaux_HLAG,nuc_classiques_HLAG,nuc_intermediaires_HLAG,nuc_non_classiques_HLAG,ncol = 4, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(nuc_temps_totaux_HLAG,nuc_temps_classiques_HLAG,nuc_temps_intermediaires_HLAG,nuc_temps_non_classiques_HLAG,ncol = 4, nrow = 1, labels = c("E","F","G","H"), align = "h"),
  plot_grid(cell_totaux_IL10,cell_classiques_IL10,cell_intermediaires_IL10,cell_non_classiques_IL10,ncol = 4, nrow = 1, labels = c("I","J","K","L"), align = "h"),
  plot_grid(cell_temps_totaux_IL10,cell_temps_classiques_IL10,cell_temps_intermediaires_IL10,cell_temps_non_classiques_IL10,ncol = 4, nrow = 1, labels = c("M","N","O","P"), align = "h"),
  nrow = 4, ncol = 1, aligns = "v")
ggsave("HLAG_IL10_AHR.png",layout_HLAG_IL10_AHR, width = 12, height = 18, units = "in", dpi = 600)

layout_IFNa_TLR7 = plot_grid(
  plot_grid(cell_totaux_IFNa,cell_classiques_IFNa,cell_intermediaires_IFNa,cell_non_classiques_IFNa,ncol = 4, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(cell_temps_totaux_IFNa,cell_temps_classiques_IFNa,cell_temps_intermediaires_IFNa,cell_temps_non_classiques_IFNa,ncol = 4, nrow = 1, labels = c("E","F","G","H"), align = "h"),
  nrow = 2, ncol = 1, aligns = "v")
ggsave("IFNa_TLR7.png",layout_IFNa_TLR7, width = 10, height = 11, units = "in", dpi = 600)

layout_pDC = ggdraw() +
  draw_plot(cell_abondance_pDC,0.2,0.66,0.6,0.3) +
  draw_plot(cell_pDC_HLADR,0.04,0.31,0.3,0.34) +
  draw_plot(cell_pDC_CD83,0.35,0.31,0.3,0.34) +
  draw_plot(cell_pDC_IFNa,0.66,0.31,0.3,0.34) +
  draw_plot(cell_temps_pDC_HLADR,0.04,0,0.3,0.3) +
  draw_plot(cell_temps_pDC_CD83,0.35,0,0.3,0.3) +
  draw_plot(cell_temps_pDC_IFNa,0.66,0,0.3,0.3) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_plot_label(c("A","B","C","D","E","F","G"), c(0.2,0.04,0.35,0.66,0.04,0.35,0.66), c(0.98,0.66,0.66,0.66,0.31,0.31,0.31), size = 15)
ggsave("pDC.png", width = 15, height = 15, units = "in", dpi = 600)

layout_supplementaire = plot_grid(
  plot_grid(cell_totaux_CX3CR1, cell_classiques_CX3CR1, cell_intermediaires_CX3CR1, cell_non_classiques_CX3CR1, ncol = 4, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(nuc_totaux_CCR2, nuc_classiques_CCR2, nuc_intermediaires_CCR2, nuc_non_classiques_CCR2, ncol = 4, nrow = 1, labels = c("E","F","G","H"), align = "h"),
  plot_grid(cell_totaux_CD83,cell_classiques_CD83,cell_intermediaires_CD83,cell_non_classiques_CD83, ncol = 4, nrow = 1, labels = c("I","J","K","L"), align = "h"),
  nrow = 3, ncol = 1, align = "v")
ggsave("supplementaire.png", layout_supplementaire, width = 12, height = 18, units = "in", dpi = 600)

layout_supp_RNA = plot_grid(
  plot_grid(RNA_SEMA4A, RNA_SLC36A1, RNA_SERINC5, ncol = 3, nrow = 1, labels = "AUTO", align = "h"),
  plot_grid(RNA_IL17RA, RNA_CTSD, RNA_GAA, ncol = 3, nrow = 1, labels = c("D","E","F"), align = "h"),
  nrow = 2, ncol = 1, align = "v")
ggsave("supp_RNA.png", layout_supp_RNA, width = 12, height = 8, units = "in", dpi = 600)
