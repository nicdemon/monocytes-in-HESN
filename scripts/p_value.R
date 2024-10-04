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
library(readr)

### Functions ###

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

tableau = function(i,df_t){
  df = data.frame(matrix(nrow = 3, ncol = 1))
  row.names(df) = c("T-test HESN vs HLA", "T-test HESN vs TS+VIH-", "T-test HESN vs TS+VIH+")
  colnames(df) = c("p_value")
  df$p_value[1:3] = df_t[1:3,i]
  return(df)
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

### Load data

cell_totaux <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_totaux.csv"))
cell_classiques <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_classiques.csv"))
cell_intermediaires <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_intermediaires.csv"))
cell_non_classiques <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_non-classiques.csv"))
cell_pDC <- as.data.frame(read_csv("../FACS/cellulaire/cellulaire_pDC.csv"))
nuc_totaux <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_totaux.csv"))
nuc_classiques <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_classiques.csv"))
nuc_intermediaires <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_intermediaires.csv"))
nuc_non_classiques <- as.data.frame(read_csv("../FACS/nucleaire/nucleaire_non-classiques.csv"))

## FACS Cellulaire

# Rownames
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
# Petits df par pops
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

# Outliers
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

## FACS nucléaire

# Rownames
row.names(nuc_totaux) = nuc_totaux$Participantes
nuc_totaux$Participantes = NULL
row.names(nuc_classiques) = nuc_classiques$Participantes
nuc_classiques$Participantes = NULL
row.names(nuc_intermediaires) = nuc_intermediaires$Participantes
nuc_intermediaires$Participantes = NULL
row.names(nuc_non_classiques) = nuc_non_classiques$Participantes
nuc_non_classiques$Participantes = NULL

# Petits df par pops
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

# Outliers
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


### Graphs de l'article

## Cellulaire

# Totaux
cell_totaux = rbind(cell_HESN_tot[1:9,],cell_TS_tot[1:11,],cell_VIH_tot[1:9,],cell_HLA_tot[1:7,])
cell_totaux_CD16 = moustache_FACS(cell_totaux,"Total",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500)) #CD16
cell_totaux_CX3CR1 = moustache_FACS(cell_totaux,"Total",5,max = 1450, bars = c(1400,1350,1250,1150,0,0)) #CX3CR1
cell_totaux_CD83 = moustache_FACS(cell_totaux, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0)) #CD83
cell_totaux_IL10 = moustache_FACS(cell_totaux, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0)) #IL-10
cell_totaux_IFNa = moustache_FACS(cell_totaux,"Total",11,max = 2550, bars = c(2450,2300,2150,0,0,0)) #IFN-a
#classiques
cell_classiques = rbind(cell_HESN_classiques[1:9,],cell_TS_classiques[1:11,],cell_VIH_classiques[1:9,],cell_HLA_classiques[1:7,])
cell_classiques_CD16 = moustache_FACS(cell_classiques,"Classical",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_classiques_CX3CR1 = moustache_FACS(cell_classiques,"Classical",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_classiques_CD83 = moustache_FACS(cell_classiques, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_classiques_IL10 = moustache_FACS(cell_classiques, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_classiques_IFNa = moustache_FACS(cell_classiques,"Classical",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
#intermédiaires
cell_intermediaires = rbind(cell_HESN_intermediaires[1:9,],cell_TS_intermediaires[1:11,],cell_VIH_intermediaires[1:9,],cell_HLA_intermediaires[1:7,])
cell_intermediaires_CD16 = moustache_FACS(cell_intermediaires,"Intermediate",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_intermediaires_CX3CR1 = moustache_FACS(cell_intermediaires,"Intermediate",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_intermediaires_CD83 = moustache_FACS(cell_intermediaires, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_intermediaires_IL10 = moustache_FACS(cell_intermediaires, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_intermediaires_IFNa = moustache_FACS(cell_intermediaires,"Intermediate",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
#non-classiques
cell_non_classiques = rbind(cell_HESN_non_c[1:9,],cell_TS_non_c[1:11,],cell_VIH_non_c[1:9,],cell_HLA_non_c[1:7,])
cell_non_classiques_CD16 = moustache_FACS(cell_non_classiques,"Non-Classical",2,max = 17000, bars = c(16500, 15500, 14500, 13500, 12500, 11500))
cell_non_classiques_CX3CR1 = moustache_FACS(cell_non_classiques,"Non-Classical",5,max = 1450, bars = c(1400,1350,1250,1150,0,0))
cell_non_classiques_CD83 = moustache_FACS(cell_non_classiques, spop = NULL,i = 8,max = 500, bars = c(475,0,0,0,0,0))
cell_non_classiques_IL10 = moustache_FACS(cell_non_classiques, spop = NULL,i = 10,max = 680, bars = c(650,630,0,0,0,0))
cell_non_classiques_IFNa = moustache_FACS(cell_non_classiques,"Non-Classical",11,max = 2550, bars = c(2450,2300,2150,0,0,0))
#pDC
cell_pDC = rbind(cell_HESN_pDC[1:9,],cell_TS_pDC[1:11,],cell_VIH_pDC[1:9,],cell_HLA_pDC[1:7,])
cell_pDC_HLADR = moustache_FACS(cell_pDC, spop = NULL,1,min = 200,max = 650)
cell_pDC_CD83 = moustache_FACS(cell_pDC, spop = NULL,6,min = 10,max = 300)
cell_pDC_IFNa = moustache_FACS(cell_pDC, spop = NULL,9,max = 550, bars = c(500,470,440,0,0,0))

## Nucléaire

#Totaux
nuc_totaux = rbind(nuc_HESN_tot[1:9,],nuc_TS_tot[1:11,],nuc_VIH_tot[1:9,],nuc_HLA_tot[1:7,])
nuc_totaux_HLADR = moustache_FACS(nuc_totaux, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_totaux_HLAG = moustache_FACS(nuc_totaux,"Total",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_totaux_CCR2 = moustache_FACS(nuc_totaux, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_totaux_CD11b = moustache_FACS(nuc_totaux,"Total",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_totaux_CD11c = moustache_FACS(nuc_totaux, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_totaux_NR4A1 = moustache_FACS(nuc_totaux,"Total",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#classiques
nuc_classiques = rbind(nuc_HESN_classiques[1:9,],nuc_TS_classiques[1:11,],nuc_VIH_classiques[1:9,],nuc_HLA_classiques[1:7,])
nuc_classiques_HLADR = moustache_FACS(nuc_classiques, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_classiques_HLAG = moustache_FACS(nuc_classiques,"Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_classiques_CCR2 = moustache_FACS(nuc_classiques, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_classiques_CD11b = moustache_FACS(nuc_classiques,"Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_classiques_CD11c = moustache_FACS(nuc_classiques, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_classiques_NR4A1 = moustache_FACS(nuc_classiques,"Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#intermédiaires
nuc_intermediaires = rbind(nuc_HESN_intermediaires[1:9,],nuc_TS_intermediaires[1:11,],nuc_VIH_intermediaires[1:9,],nuc_HLA_intermediaires[1:7,])
nuc_intermediaires_HLADR = moustache_FACS(nuc_intermediaires, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_intermediaires_HLAG = moustache_FACS(nuc_intermediaires,"Intermediate",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_intermediaires_CCR2 = moustache_FACS(nuc_intermediaires, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_intermediaires_CD11b = moustache_FACS(nuc_intermediaires,"Intermediate",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_intermediaires_CD11c = moustache_FACS(nuc_intermediaires, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_intermediaires_NR4A1 = moustache_FACS(nuc_intermediaires,"Intermediate",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
#non-classiques
nuc_non_classiques = rbind(nuc_HESN_non_c[1:9,],nuc_TS_non_c[1:11,],nuc_VIH_non_c[1:9,],nuc_HLA_non_c[1:7,])
nuc_non_classiques_HLADR = moustache_FACS(nuc_non_classiques, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_non_classiques_HLAG = moustache_FACS(nuc_non_classiques,"Non-Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_non_classiques_CCR2 = moustache_FACS(nuc_non_classiques, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_non_classiques_CD11b = moustache_FACS(nuc_non_classiques,"Non-Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_non_classiques_CD11c = moustache_FACS(nuc_non_classiques, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_non_classiques_NR4A1 = moustache_FACS(nuc_non_classiques,"Non-Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))

### Tableaux de p-values

## T-test cellulaire
df_t_tot = data.frame(matrix(ncol = 12, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_tot) = colnames(cell_totaux)[1:12]

for(i in colnames(df_t_tot)){
  tryCatch({df_t_tot[1,i] = t.test(x = cell_HESN_tot[1:9,i], y = cell_HLA_tot[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[2,i] = t.test(x = cell_HESN_tot[1:9,i], y = cell_TS_tot[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[3,i] = t.test(x = cell_HESN_tot[1:9,i], y = cell_VIH_tot[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_classiques = data.frame(matrix(ncol = 12, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_classiques) = colnames(cell_classiques)[1:12]

for(i in colnames(df_t_classiques)){
  tryCatch({df_t_classiques[1,i] = t.test(x = cell_HESN_classiques[1:9,i], y = cell_HLA_classiques[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[2,i] = t.test(x = cell_HESN_classiques[1:9,i], y = cell_TS_classiques[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[3,i] = t.test(x = cell_HESN_classiques[1:9,i], y = cell_VIH_classiques[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_intermediaires = data.frame(matrix(ncol = 12, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_intermediaires) = colnames(cell_intermediaires)[1:12]

for(i in colnames(df_t_intermediaires)){
  tryCatch({df_t_intermediaires[1,i] = t.test(x = cell_HESN_intermediaires[1:9,i], y = cell_HLA_intermediaires[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[2,i] = t.test(x = cell_HESN_intermediaires[1:9,i], y = cell_TS_intermediaires[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[3,i] = t.test(x = cell_HESN_intermediaires[1:9,i], y = cell_VIH_intermediaires[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_non_c = data.frame(matrix(ncol = 12, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_non_c) = colnames(cell_non_classiques)[1:12]

for(i in colnames(df_t_non_c)){
  tryCatch({df_t_non_c[1,i] = t.test(x = cell_HESN_non_c[1:9,i], y = cell_HLA_non_c[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[2,i] = t.test(x = cell_HESN_non_c[1:9,i], y = cell_TS_non_c[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[3,i] = t.test(x = cell_HESN_non_c[1:9,i], y = cell_VIH_non_c[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_pDC = data.frame(matrix(ncol = 10, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_pDC) = colnames(cell_pDC)[1:10]

for(i in colnames(df_t_pDC)){
  tryCatch({df_t_pDC[1,i] = t.test(x = cell_HESN_pDC[1:9,i], y = cell_HLA_pDC[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_pDC[2,i] = t.test(x = cell_HESN_pDC[1:9,i], y = cell_TS_pDC[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_pDC[3,i] = t.test(x = cell_HESN_pDC[1:9,i], y = cell_VIH_pDC[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

# Save tableaux recap
write.csv(df_t_tot, "Cellulaire_T-test_totaux.csv")
write.csv(df_t_classiques, "Cellulaire_T-test_classiques.csv")
write.csv(df_t_intermediaires, "Cellulaire_T-test_intermediaires.csv")
write.csv(df_t_non_c, "Cellulaire_T-test_non-classiques.csv")
write.csv(df_t_pDC, "Cellulaire_T-test_pDC.csv")
write.csv(cell_HESN_tot,"Cellulaire_cell_HESN_totaux_normalite.csv")
write.csv(cell_HESN_classiques,"Cellulaire_cell_HESN_classiques_normalite.csv")
write.csv(cell_HESN_intermediaires,"Cellulaire_cell_HESN_intermediaires_normalite.csv")
write.csv(cell_HESN_non_c,"Cellulaire_cell_HESN_non-classiques_normalite.csv")
write.csv(cell_HESN_pDC,"Cellulaire_cell_HESN_pDC_normalite.csv")
write.csv(cell_TS_tot,"Cellulaire_cell_TS_totaux_normalite.csv")
write.csv(cell_TS_classiques,"Cellulaire_cell_TS_classiques_normalite.csv")
write.csv(cell_TS_intermediaires,"Cellulaire_cell_TS_intermediaires_normalite.csv")
write.csv(cell_TS_non_c,"Cellulaire_cell_TS_non-classiques_normalite.csv")
write.csv(cell_TS_pDC,"Cellulaire_cell_TS_pDC_normalite.csv")
write.csv(cell_VIH_tot,"Cellulaire_cell_VIH_totaux_normalite.csv")
write.csv(cell_VIH_classiques,"Cellulaire_cell_VIH_classiques_normalite.csv")
write.csv(cell_VIH_intermediaires,"Cellulaire_cell_VIH_intermediaires_normalite.csv")
write.csv(cell_VIH_non_c,"Cellulaire_cell_VIH_non-classiques_normalite.csv")
write.csv(cell_VIH_pDC,"Cellulaire_cell_VIH_pDC_normalite.csv")
write.csv(cell_HLA_tot,"Cellulaire_cell_HLA_totaux_normalite.csv")
write.csv(cell_HLA_classiques,"Cellulaire_cell_HLA_classiques_normalite.csv")
write.csv(cell_HLA_intermediaires,"Cellulaire_cell_HLA_intermediaires_normalite.csv")
write.csv(cell_HLA_non_c,"Cellulaire_cell_HLA_non-classiques_normalite.csv")
write.csv(cell_HLA_pDC,"Cellulaire_cell_HLA_pDC_normalite.csv")

## Tableaux t-test pour layouts

# Totaux
cell_tt_CD16 = tableau(2,df_t_tot)
cell_tt_CX3CR1 = tableau(5,df_t_tot)
cell_tt_CD83 = tableau(8,df_t_tot)
cell_tt_IL10 = tableau(10,df_t_tot)
cell_tt_IFNa = tableau(11,df_t_tot)
# Classiques
cell_ct_CD16 = tableau(2,df_t_classiques)
cell_ct_CX3CR1 = tableau(5,df_t_classiques)
cell_ct_CD83 = tableau(8,df_t_classiques)
cell_ct_IL10 = tableau(10,df_t_classiques)
cell_ct_IFNa = tableau(11,df_t_classiques)
# Intermédiaires
cell_it_CD16 = tableau(2,df_t_intermediaires)
cell_it_CX3CR1 = tableau(5,df_t_intermediaires)
cell_it_CD83 = tableau(8,df_t_intermediaires)
cell_it_IL10 = tableau(10,df_t_intermediaires)
cell_it_IFNa = tableau(11,df_t_intermediaires)
# Non-classiques
cell_nt_CD16 = tableau(2,df_t_non_c)
cell_nt_CX3CR1 = tableau(5,df_t_non_c)
cell_nt_CD83 = tableau(8,df_t_non_c)
cell_nt_IL10 = tableau(10,df_t_non_c)
cell_nt_IFNa = tableau(111,df_t_non_c)
# pDC
cell_pt_HLADR = tableau(1,df_t_pDC)
cell_pt_CD83 = tableau(6,df_t_pDC)
cell_pt_IFNa = tableau(9,df_t_pDC)

## T-test nucléaire
df_t_tot = data.frame(matrix(ncol = 9, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_tot) = colnames(nuc_totaux)[1:9]

for(i in colnames(df_t_tot)){
  tryCatch({df_t_tot[1,i] = t.test(x = nuc_HESN_tot[1:9,i], y = nuc_HLA_tot[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[2,i] = t.test(x = nuc_HESN_tot[1:9,i], y = nuc_TS_tot[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[3,i] = t.test(x = nuc_HESN_tot[1:9,i], y = cell_VIH_tot[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_classiques = data.frame(matrix(ncol = 9, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_classiques) = colnames(nuc_classiques)[1:9]

for(i in colnames(df_t_classiques)){
  tryCatch({df_t_classiques[1,i] = t.test(x = nuc_HESN_classiques[1:9,i], y = nuc_HLA_classiques[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[2,i] = t.test(x = nuc_HESN_classiques[1:9,i], y = nuc_TS_classiques[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[3,i] = t.test(x = nuc_HESN_classiques[1:9,i], y = cell_VIH_classiques[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_intermediaires = data.frame(matrix(ncol = 9, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_intermediaires) = colnames(nuc_intermediaires)[1:9]

for(i in colnames(df_t_intermediaires)){
  tryCatch({df_t_intermediaires[1,i] = t.test(x = nuc_HESN_intermediaires[1:9,i], y = nuc_HLA_intermediaires[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[2,i] = t.test(x = nuc_HESN_intermediaires[1:9,i], y = nuc_TS_intermediaires[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[3,i] = t.test(x = nuc_HESN_intermediaires[1:9,i], y = cell_VIH_intermediaires[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_non_c = data.frame(matrix(ncol = 9, nrow = 3),row.names = c("HESN_HLA","HESN_TS","HESN_VIH"))
colnames(df_t_non_c) = colnames(nuc_non_classiques)[1:9]

for(i in colnames(df_t_non_c)){
  tryCatch({df_t_non_c[1,i] = t.test(x = nuc_HESN_non_c[1:9,i], y = nuc_HLA_non_c[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[2,i] = t.test(x = nuc_HESN_non_c[1:9,i], y = nuc_TS_non_c[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[3,i] = t.test(x = nuc_HESN_non_c[1:9,i], y = cell_VIH_non_c[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

# Save tableaux recap
write.csv(df_t_tot, "Nucleaire_T-test_totaux.csv")
write.csv(df_t_classiques, "Nucleaire_T-test_classiques.csv")
write.csv(df_t_intermediaires, "Nucleaire_T-test_intermediaires.csv")
write.csv(df_t_non_c, "Nucleaire_T-test_non-classiques.csv")
write.csv(nuc_HESN_tot,"Nucleaire_nuc_HESN_totaux_normalite.csv")
write.csv(nuc_HESN_classiques,"Nucleaire_nuc_HESN_classiques_normalite.csv")
write.csv(nuc_HESN_intermediaires,"Nucleaire_nuc_HESN_intermediaires_normalite.csv")
write.csv(nuc_HESN_non_c,"Nucleaire_nuc_HESN_non-classiques_normalite.csv")
write.csv(nuc_TS_tot,"Nucleaire_nuc_TS_totaux_normalite.csv")
write.csv(nuc_TS_classiques,"Nucleaire_nuc_TS_classiques_normalite.csv")
write.csv(nuc_TS_intermediaires,"Nucleaire_nuc_TS_intermediaires_normalite.csv")
write.csv(nuc_TS_non_c,"Nucleaire_nuc_TS_non-classiques_normalite.csv")
write.csv(cell_VIH_tot,"Nucleaire_cell_VIH_totaux_normalite.csv")
write.csv(cell_VIH_classiques,"Nucleaire_cell_VIH_classiques_normalite.csv")
write.csv(cell_VIH_intermediaires,"Nucleaire_cell_VIH_intermediaires_normalite.csv")
write.csv(cell_VIH_non_c,"Nucleaire_cell_VIH_non-classiques_normalite.csv")
write.csv(nuc_HLA_tot,"Nucleaire_nuc_HLA_totaux_normalite.csv")
write.csv(nuc_HLA_classiques,"Nucleaire_nuc_HLA_classiques_normalite.csv")
write.csv(nuc_HLA_intermediaires,"Nucleaire_nuc_HLA_intermediaires_normalite.csv")
write.csv(nuc_HLA_non_c,"Nucleaire_nuc_HLA_non-classiques_normalite.csv")

# Tableaux t-test pour layouts
# Totaux
nuc_tt_HLADR = tableau(3,df_t_tot)
nuc_tt_HLAG = tableau(4,df_t_tot)
nuc_tt_CCR2 = tableau(5,df_t_tot)
nuc_tt_CD11b = tableau(7,df_t_tot)
nuc_tt_CD11c = tableau(8,df_t_tot)
nuc_tt_NR4A1 = tableau(9,df_t_tot)
# Classiques
nuc_ct_HLADR = tableau(3,df_t_classiques)
nuc_ct_HLAG = tableau(4,df_t_classiques)
nuc_ct_CCR2 = tableau(5,df_t_classiques)
nuc_ct_CD11b = tableau(7,df_t_classiques)
nuc_ct_CD11c = tableau(8,df_t_classiques)
nuc_ct_NR4A1 = tableau(9,df_t_classiques)
# Intermédiaires
nuc_it_HLADR = tableau(3,df_t_intermediaires)
nuc_it_HLAG = tableau(4,df_t_intermediaires)
nuc_it_CCR2 = tableau(5,df_t_intermediaires)
nuc_it_CD11b = tableau(7,df_t_intermediaires)
nuc_it_CD11c = tableau(8,df_t_intermediaires)
nuc_it_NR4A1 = tableau(9,df_t_intermediaires)
# Non-classiques
nuc_nt_HLADR = tableau(3,df_t_non_c)
nuc_nt_HLAG = tableau(4,df_t_non_c)
nuc_nt_CCR2 = tableau(5,df_t_non_c)
nuc_nt_CD11b = tableau(7,df_t_non_c)
nuc_nt_CD11c = tableau(8,df_t_non_c)
nuc_nt_NR4A1 = tableau(9,df_t_non_c)

### Layouts

## Monocytes cellulaire

# CD16
layout_cell_CD16 = grid.arrange(cell_totaux_CD16, cell_classiques_CD16, cell_intermediaires_CD16, cell_non_classiques_CD16, tableGrob(cell_tt_CD16), tableGrob(cell_ct_CD16), tableGrob(cell_it_CD16), tableGrob(cell_nt_CD16), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD16 Monocytes cellulaires")
# IL-10
layout_cell_IL10 = grid.arrange(cell_totaux_IL10, cell_classiques_IL10, cell_intermediaires_IL10, cell_non_classiques_IL10, tableGrob(cell_tt_IL10), tableGrob(cell_ct_IL10), tableGrob(cell_it_IL10), tableGrob(cell_nt_IL10), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "IL10 Monocytes cellulaires")
# CX3CR1
layout_cell_CX3CR1 = grid.arrange(cell_totaux_CX3CR1, cell_classiques_CX3CR1, cell_intermediaires_CX3CR1, cell_non_classiques_CX3CR1, tableGrob(cell_tt_CX3CR1), tableGrob(cell_ct_CX3CR1), tableGrob(cell_it_CX3CR1), tableGrob(cell_nt_CX3CR1), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CX3CR1 Monocytes cellulaires")
# IFN-a
layout_cell_IFNa = grid.arrange(cell_totaux_IFNa, cell_classiques_IFNa, cell_intermediaires_IFNa, cell_non_classiques_IFNa, tableGrob(cell_tt_IFNa), tableGrob(cell_ct_IFNa), tableGrob(cell_it_IFNa), tableGrob(cell_nt_IFNa), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "IFNa Monocytes cellulaires")
# CD83
layout_cell_CD83 = grid.arrange(cell_totaux_CD83, cell_classiques_CD83, cell_intermediaires_CD83, cell_non_classiques_CD83, tableGrob(cell_tt_CD83), tableGrob(cell_ct_CD83), tableGrob(cell_it_CD83), tableGrob(cell_nt_CD83), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD83 Monocytes cellulaires")

## pDC cellulaire

layout_pDC = grid.arrange(cell_pDC_CD83, cell_pDC_HLADR, cell_pDC_IFNa, tableGrob(cell_pt_CD83), tableGrob(cell_pt_HLADR), tableGrob(cell_pt_IFNa), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "pDC cellulaires")

## Monocytes nucléaire

# NR4A1
layout_nuc_NR4A1 = grid.arrange(nuc_totaux_NR4A1, nuc_classiques_NR4A1, nuc_intermediaires_NR4A1, nuc_non_classiques_NR4A1, tableGrob(nuc_tt_NR4A1), tableGrob(nuc_ct_NR4A1), tableGrob(nuc_it_NR4A1), tableGrob(nuc_nt_NR4A1), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "NR4A1 Monocytes nucléaires")
# HLADR
layout_nuc_HLADR = grid.arrange(nuc_totaux_HLADR, nuc_classiques_HLADR, nuc_intermediaires_HLADR, nuc_non_classiques_HLADR, tableGrob(nuc_tt_HLADR), tableGrob(nuc_ct_HLADR), tableGrob(nuc_it_HLADR), tableGrob(nuc_nt_HLADR), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "HLADR Monocytes nucléaires")
# CD11b
layout_nuc_CD11b = grid.arrange(nuc_totaux_CD11b, nuc_classiques_CD11b, nuc_intermediaires_CD11b, nuc_non_classiques_CD11b, tableGrob(nuc_tt_CD11b), tableGrob(nuc_ct_CD11b), tableGrob(nuc_it_CD11b), tableGrob(nuc_nt_CD11b), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD11b Monocytes nucléaires")
# CD11c
layout_nuc_CD11c = grid.arrange(nuc_totaux_CD11c, nuc_classiques_CD11c, nuc_intermediaires_CD11c, nuc_non_classiques_CD11c, tableGrob(nuc_tt_CD11c), tableGrob(nuc_ct_CD11c), tableGrob(nuc_it_CD11c), tableGrob(nuc_nt_CD11c), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD11c Monocytes nucléaires")
# HLA-G
layout_nuc_HLAG = grid.arrange(nuc_totaux_HLAG, nuc_classiques_HLAG, nuc_intermediaires_HLAG, nuc_non_classiques_HLAG, tableGrob(nuc_tt_HLAG), tableGrob(nuc_ct_HLAG), tableGrob(nuc_it_HLAG), tableGrob(nuc_nt_HLAG), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "HLAG Monocytes nucléaires")
# CCR2
layout_nuc_CCR2 = grid.arrange(nuc_totaux_CCR2, nuc_classiques_CCR2, nuc_intermediaires_CCR2, nuc_non_classiques_CCR2, tableGrob(nuc_tt_CCR2), tableGrob(nuc_ct_CCR2), tableGrob(nuc_it_CCR2), tableGrob(nuc_nt_CCR2), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CCR2 Monocytes nucléaires")

## Save all in one pdf

nuc_non_classiques = rbind(nuc_HESN_non_c[1:9,],nuc_TS_non_c[1:11,],nuc_VIH_non_c[1:9,],nuc_HLA_non_c[1:7,])
nuc_non_classiques_HLADR = moustache_FACS(nuc_non_classiques, spop = NULL,i = 3,max = 2000, bars = c(1800,1600,1400,0,0,0))
nuc_non_classiques_HLAG = moustache_FACS(nuc_non_classiques,"Non-Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
nuc_non_classiques_CCR2 = moustache_FACS(nuc_non_classiques, spop = NULL,i = 5,max = 1900, bars = c(1800,1600,0,0,0,0))
nuc_non_classiques_CD11b = moustache_FACS(nuc_non_classiques,"Non-Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
nuc_non_classiques_CD11c = moustache_FACS(nuc_non_classiques, spop = NULL,i = 8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
nuc_non_classiques_NR4A1 = moustache_FACS(nuc_non_classiques,"Non-Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))

ggsave("layout_pvalues.pdf", marrangeGrob(list(layout_cell_CD16, layout_cell_CX3CR1, layout_cell_CD83, layout_cell_IL10, layout_cell_IFNa, layout_pDC, layout_nuc_HLADR, layout_nuc_HLAG, layout_nuc_CCR2, layout_nuc_CD11b, layout_nuc_CD11c, layout_nuc_NR4A1), nrow = 1, ncol = 1), width = 20, height = 8.5, units = "in")
