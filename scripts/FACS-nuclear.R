#WORK_DIR + dependencies
setwd("../data/FACS/nuclear")
library(data.table)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(gtable)
library(digest)
library(ggpubr)

###Import data
totaux <- as.data.frame(read_csv("nucleaire_totaux.csv"))
classiques <- as.data.frame(read_csv("nucleaire_classiques.csv"))
intermediaires <- as.data.frame(read_csv("nucleaire_intermediaires.csv"))
non_classiques <- as.data.frame(read_csv("nucleaire_non-classiques.csv"))
#pDC <- as.data.frame(read_csv("nucleaire_pDC.csv"))

###Visualisation

#Rownames
row.names(totaux) = totaux$Participantes
totaux$Participantes = NULL
row.names(classiques) = classiques$Participantes
classiques$Participantes = NULL
row.names(intermediaires) = intermediaires$Participantes
intermediaires$Participantes = NULL
row.names(non_classiques) = non_classiques$Participantes
non_classiques$Participantes = NULL
#row.names(pDC) = pDC$Participantes
#pDC$Participantes = NULL

#Dotplot
#par(mfrow=c(3,3))
#for (i in c(2:10)){
#  plot(c(1:76),totaux[,i], ylab = colnames(totaux)[i], xlab = "")
#}
#mtext("Fluorescence du cocktail nucléaire des monocytes totaux pré-traitement", side=3, outer=TRUE, line=-3)
#par(mfrow=c(1,1))

#par(mfrow=c(2,4))
#for (i in c(2:18)){
#  plot(c(1:22),pDC[,i], ylab = colnames(pDC)[i], xlab = "")
#}
#mtext("Fluorescence du cocktail nucléaire des pDC pré-traitement", side=3, outer=TRUE, line=-3)
#par(mfrow=c(1,1))

###Split par populations

#Petits df par pops
#HESN[1:9], TS[10:20], VIH[21:29], HLA[30:36]
HESN_tot = as.data.frame(totaux[1:9,2:10])
HESN_classiques = as.data.frame(classiques[1:9,2:10])
HESN_intermediaires = as.data.frame(intermediaires[1:9,2:10])
HESN_non_c = as.data.frame(non_classiques[1:9,2:10])
#HESN_pDC = as.data.frame(pDC[1:9,2:8])
row.names(HESN_tot) = row.names(totaux)[1:9]
row.names(HESN_classiques) = row.names(classiques)[1:9]
row.names(HESN_intermediaires) = row.names(intermediaires)[1:9]
row.names(HESN_non_c) = row.names(non_classiques)[1:9]
#row.names(HESN_pDC) = row.names(pDC)[1:9]

TS_tot = as.data.frame(totaux[10:20,2:10])
TS_classiques = as.data.frame(classiques[10:20,2:10])
TS_intermediaires = as.data.frame(intermediaires[10:20,2:10])
TS_non_c = as.data.frame(non_classiques[10:20,2:10])
#TS_pDC = as.data.frame(pDC[10:20,2:8])
row.names(TS_tot) = row.names(totaux)[10:20]
row.names(TS_classiques) = row.names(classiques)[10:20]
row.names(TS_intermediaires) = row.names(intermediaires)[10:20]
row.names(TS_non_c) = row.names(non_classiques)[10:20]
#row.names(TS_pDC) = row.names(pDC)[10:20]

VIH_tot = as.data.frame(totaux[21:29,2:10])
VIH_classiques = as.data.frame(classiques[21:29,2:10])
VIH_intermediaires = as.data.frame(intermediaires[21:29,2:10])
VIH_non_c = as.data.frame(non_classiques[21:29,2:10])
#VIH_pDC = as.data.frame(pDC[21:29,2:8])
row.names(VIH_tot) = row.names(totaux)[21:29]
row.names(VIH_classiques) = row.names(classiques)[21:29]
row.names(VIH_intermediaires) = row.names(intermediaires)[21:29]
row.names(VIH_non_c) = row.names(non_classiques)[21:29]
#row.names(VIH_pDC) = row.names(pDC)[21:29]

HLA_tot = as.data.frame(totaux[30:36,2:10])
HLA_classiques = as.data.frame(classiques[30:36,2:10])
HLA_intermediaires = as.data.frame(intermediaires[30:36,2:10])
HLA_non_c = as.data.frame(non_classiques[30:36,2:10])
#HLA_pDC = as.data.frame(pDC[30:36,2:8])
row.names(HLA_tot) = row.names(totaux)[30:36]
row.names(HLA_classiques) = row.names(classiques)[30:36]
row.names(HLA_intermediaires) = row.names(intermediaires)[30:36]
row.names(HLA_non_c) = row.names(non_classiques)[30:36]
#row.names(HLA_pDC) = row.names(pDC)[30:36]

###Stats

#Outliers
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

HESN_tot = outliers(HESN_tot)
HESN_classiques = outliers(HESN_classiques)
HESN_intermediaires = outliers(HESN_intermediaires)
HESN_non_c = outliers(HESN_non_c)
#HESN_pDC = outliers(HESN_pDC)
TS_tot = outliers(TS_tot)
TS_classiques = outliers(TS_classiques)
TS_intermediaires = outliers(TS_intermediaires)
TS_non_c = outliers(TS_non_c)
#TS_pDC = outliers(TS_pDC)
VIH_tot = outliers(VIH_tot)
VIH_classiques = outliers(VIH_classiques)
VIH_intermediaires = outliers(VIH_intermediaires)
VIH_non_c = outliers(VIH_non_c)
#VIH_pDC = outliers(VIH_pDC)
HLA_tot = outliers(HLA_tot)
HLA_classiques = outliers(HLA_classiques)
HLA_intermediaires = outliers(HLA_intermediaires)
HLA_non_c = outliers(HLA_non_c)
#HLA_pDC = outliers(HLA_pDC)

#Normalité
normalite = function(df){
  maxrow = length(rownames(df))+1
  for(i in colnames(df)){
    tryCatch({df[maxrow,i] = shapiro.test(na.omit(df[,i]))[2]},error=function(cond){return(NULL)})
  }
  rownames(df)[length(rownames(df))] = 'Normalite_p-value'
  return(df)
}

HESN_tot = normalite(HESN_tot)
HESN_classiques = normalite(HESN_classiques)
HESN_intermediaires = normalite(HESN_intermediaires)
HESN_non_c = normalite(HESN_non_c)
HESN_pDC = normalite(HESN_pDC)
TS_tot = normalite(TS_tot)
TS_classiques = normalite(TS_classiques)
TS_intermediaires = normalite(TS_intermediaires)
TS_non_c = normalite(TS_non_c)
TS_pDC = normalite(TS_pDC)
VIH_tot = normalite(VIH_tot)
VIH_classiques = normalite(VIH_classiques)
VIH_intermediaires = normalite(VIH_intermediaires)
VIH_non_c = normalite(VIH_non_c)
VIH_pDC = normalite(VIH_pDC)
HLA_tot = normalite(HLA_tot)
HLA_classiques = normalite(HLA_classiques)
HLA_intermediaires = normalite(HLA_intermediaires)
HLA_non_c = normalite(HLA_non_c)
HLA_pDC = normalite(HLA_pDC)

#T-test entre pops par marqueur
df_t_tot = data.frame(matrix(ncol = 9, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
colnames(df_t_tot) = colnames(totaux)[2:10]

for(i in colnames(df_t_tot)){
  tryCatch({df_t_tot[1,i] = t.test(x = HESN_tot[1:9,i], y = HLA_tot[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[2,i] = t.test(x = HESN_tot[1:9,i], y = TS_tot[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[3,i] = t.test(x = HESN_tot[1:9,i], y = VIH_tot[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[4,i] = t.test(x = TS_tot[1:11,i], y = HLA_tot[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[5,i] = t.test(x = TS_tot[1:11,i], y = VIH_tot[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_tot[6,i] = t.test(x = VIH_tot[1:9,i], y = HLA_tot[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_classiques = data.frame(matrix(ncol = 9, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
colnames(df_t_classiques) = colnames(classiques)[2:10]

for(i in colnames(df_t_classiques)){
  tryCatch({df_t_classiques[1,i] = t.test(x = HESN_classiques[1:9,i], y = HLA_classiques[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[2,i] = t.test(x = HESN_classiques[1:9,i], y = TS_classiques[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[3,i] = t.test(x = HESN_classiques[1:9,i], y = VIH_classiques[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[4,i] = t.test(x = TS_classiques[1:11,i], y = HLA_classiques[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[5,i] = t.test(x = TS_classiques[1:11,i], y = VIH_classiques[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_classiques[6,i] = t.test(x = VIH_classiques[1:9,i], y = HLA_classiques[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_intermediaires = data.frame(matrix(ncol = 9, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
colnames(df_t_intermediaires) = colnames(intermediaires)[2:10]

for(i in colnames(df_t_intermediaires)){
  tryCatch({df_t_intermediaires[1,i] = t.test(x = HESN_intermediaires[1:9,i], y = HLA_intermediaires[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[2,i] = t.test(x = HESN_intermediaires[1:9,i], y = TS_intermediaires[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[3,i] = t.test(x = HESN_intermediaires[1:9,i], y = VIH_intermediaires[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[4,i] = t.test(x = TS_intermediaires[1:11,i], y = HLA_intermediaires[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[5,i] = t.test(x = TS_intermediaires[1:11,i], y = VIH_intermediaires[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_intermediaires[6,i] = t.test(x = VIH_intermediaires[1:9,i], y = HLA_intermediaires[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

df_t_non_c = data.frame(matrix(ncol = 9, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
colnames(df_t_non_c) = colnames(non_classiques)[2:10]

for(i in colnames(df_t_non_c)){
  tryCatch({df_t_non_c[1,i] = t.test(x = HESN_non_c[1:9,i], y = HLA_non_c[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[2,i] = t.test(x = HESN_non_c[1:9,i], y = TS_non_c[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[3,i] = t.test(x = HESN_non_c[1:9,i], y = VIH_non_c[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[4,i] = t.test(x = TS_non_c[1:11,i], y = HLA_non_c[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[5,i] = t.test(x = TS_non_c[1:11,i], y = VIH_non_c[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t_non_c[6,i] = t.test(x = VIH_non_c[1:9,i], y = HLA_non_c[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

#df_t_pDC = data.frame(matrix(ncol = 7, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
#colnames(df_t_pDC) = colnames(pDC)[2:8]

#for(i in colnames(df_t_pDC)){
#  tryCatch({df_t_pDC[1,i] = t.test(x = HESN_pDC[1:9,i], y = HLA_pDC[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#  tryCatch({df_t_pDC[2,i] = t.test(x = HESN_pDC[1:9,i], y = TS_pDC[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#  tryCatch({df_t_pDC[3,i] = t.test(x = HESN_pDC[1:9,i], y = VIH_pDC[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#  tryCatch({df_t_pDC[4,i] = t.test(x = TS_pDC[1:11,i], y = HLA_pDC[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#  tryCatch({df_t_pDC[5,i] = t.test(x = TS_pDC[1:11,i], y = VIH_pDC[1:9,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#  tryCatch({df_t_pDC[6,i] = t.test(x = VIH_pDC[1:9,i], y = HLA_pDC[1:7,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
#}

#Tableau recap
write.csv(df_t_tot, "Nucleaire_T-test_totaux.csv")
write.csv(df_t_classiques, "Nucleaire_T-test_classiques.csv")
write.csv(df_t_intermediaires, "Nucleaire_T-test_intermediaires.csv")
write.csv(df_t_non_c, "Nucleaire_T-test_non-classiques.csv")
#write.csv(df_t_pDC, "Nucleaire_T-test_pDC.csv")
write.csv(HESN_tot,"Nucleaire_HESN_totaux_normalite.csv")
write.csv(HESN_classiques,"Nucleaire_HESN_classiques_normalite.csv")
write.csv(HESN_intermediaires,"Nucleaire_HESN_intermediaires_normalite.csv")
write.csv(HESN_non_c,"Nucleaire_HESN_non-classiques_normalite.csv")
#write.csv(HESN_pDC,"Nucleaire_HESN_pDC_normalite.csv")
write.csv(TS_tot,"Nucleaire_TS_totaux_normalite.csv")
write.csv(TS_classiques,"Nucleaire_TS_classiques_normalite.csv")
write.csv(TS_intermediaires,"Nucleaire_TS_intermediaires_normalite.csv")
write.csv(TS_non_c,"Nucleaire_TS_non-classiques_normalite.csv")
#write.csv(TS_pDC,"Nucleaire_TS_pDC_normalite.csv")
write.csv(VIH_tot,"Nucleaire_VIH_totaux_normalite.csv")
write.csv(VIH_classiques,"Nucleaire_VIH_classiques_normalite.csv")
write.csv(VIH_intermediaires,"Nucleaire_VIH_intermediaires_normalite.csv")
write.csv(VIH_non_c,"Nucleaire_VIH_non-classiques_normalite.csv")
#write.csv(VIH_pDC,"Nucleaire_VIH_pDC_normalite.csv")
write.csv(HLA_tot,"Nucleaire_HLA_totaux_normalite.csv")
write.csv(HLA_classiques,"Nucleaire_HLA_classiques_normalite.csv")
write.csv(HLA_intermediaires,"Nucleaire_HLA_intermediaires_normalite.csv")
write.csv(HLA_non_c,"Nucleaire_HLA_non-classiques_normalite.csv")
#write.csv(HLA_pDC,"Nucleaire_HLA_pDC_normalite.csv")

###Graphs marqueurs par populations

#Dotplot par population
dotplot = function(df, name="Inconnu"){
  par(mfrow=c(3,3))
  for(i in c(1:length(colnames(df)))){
    tmp = as.data.frame(df)
    tmp = tmp[-length(row.names(tmp)),]
    tmp = tmp[complete.cases(tmp[row.names(tmp), i]),]
    len = length(row.names(tmp))
    if(len > 0){
      plot(c(1:len), tmp[,i], ylab = colnames(tmp)[i], xlab = "")
    }else {
      plot.new()
    }
  }
  mtext(paste("Fluorescence des marqueurs pour la population", name), side=3, outer=TRUE, line=-3)
  par(mfrow=c(1,1))
}

dotplot(HESN_tot,"Nucleaire_HESN_totaux")
dotplot(HESN_classiques,"Nucleaire_HESN_classiques")
dotplot(HESN_intermediaires,"Nucleaire_HESN_intermediaires")
dotplot(HESN_non_c,"Nucleaire_HESN_non-classiques")
#dotplot(HESN_pDC,"Nucleaire_HESN_pDC")
dotplot(TS_tot,"Nucleaire_TS_totaux")
dotplot(TS_classiques,"Nucleaire_TS_classiques")
dotplot(TS_intermediaires,"Nucleaire_TS_intermediaires")
dotplot(TS_non_c,"Nucleaire_TS_non-classiques")
#dotplot(TS_pDC,"Nucleaire_TS_pDC")
dotplot(VIH_tot,"Nucleaire_VIH_totaux")
dotplot(VIH_classiques,"Nucleaire_VIH_classiques")
dotplot(VIH_intermediaires,"Nucleaire_VIH_intermediaires")
dotplot(VIH_non_c,"Nucleaire_VIH_non-classiques")
#dotplot(VIH_pDC,"Nucleaire_VIH_pDC")
dotplot(HLA_tot,"Nucleaire_HLA_totaux")
dotplot(HLA_classiques,"Nucleaire_HLA_classiques")
dotplot(HLA_intermediaires,"Nucleaire_HLA_intermediaires")
dotplot(HLA_non_c,"Nucleaire_HLA_non-classiques")
#dotplot(HLA_pDC,"Nucleaire_HLA_pDC")

#Boxplot par marqueur
moustache = function(df,spop,i,min = 0, max = 1000, bars = c(0,0,0,0,0,0)){
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
      comparison[[length(comparison)+1]] = c(matrice_tmp$group1[i],matrice_tmp$group2[i])
    }
  }
  tmp %>%
    arrange(marqueur) %>%
    mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x=populations, y=marqueur, group=populations)) + 
    geom_point(aes(x=populations, y=marqueur, shape = populations, size = 5), binaxis = "y", stackdir = "center") +
    scale_shape_manual(values = c(15:18)) +
    geom_errorbar(aes(ymin=moy-std,ymax=moy+std), width=0.2, position=position_dodge(0.1),size = 0.5) +
    geom_errorbar(aes(ymin=moy,ymax=moy), width=0.2, position=position_dodge(0.1),size = 0.5) +
    stat_compare_means(comparisons = comparison, method = "t.test", label = "p.signif", hide.ns = TRUE, na.rm = TRUE, label.y = bars) +
    ggtitle(as.character(spop)) +
    coord_cartesian(ylim = c(min, max)) +
    xlab("") +
    ylab(paste(name, "expressing cells (GeoMFI)", sep = " ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, vjust = 0.5))
}

tableau = function(i,HESN, TS, VIH, HLA, df_t){
  df = data.frame(matrix(nrow = 10, ncol = 1))
  row.names(df) = c("Normalité HESN", "Normalité TS+VIH-", "Normalité TS+VIH+", "Normalité HLA", "T-test HESN vs HLA", "T-test HESN vs TS+VIH-", "T-test HESN vs TS+VIH+", "T-test TS+VIH- vs HLA", "T-test TS+VIH- vs TS+VIH+", "T-test TS+VIH+ vs HLA")
  colnames(df) = c("p_value")
  len = length(row.names(HESN))
  df$p_value[1] = HESN[len,i]
  len = length(row.names(TS))
  df$p_value[2] = TS[len,i]
  len = length(row.names(VIH))
  df$p_value[3] = VIH[len,i]
  len = length(row.names(HLA))
  df$p_value[4] = HLA[len,i]
  df$p_value[5:10] = df_t[1:6,i]
  return(df)
}

#Totaux
#To change scale to uniformize in fx moustache : min/max = x
totaux = rbind(HESN_tot[1:9,],TS_tot[1:11,],VIH_tot[1:9,],HLA_tot[1:7,])
tp1 = moustache(totaux,"Total",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
tp2 = moustache(totaux,"Total",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
tp3 = moustache(totaux,"Total",3,max = 2000, bars = c(1800,1600,1400,0,0,0))
tp4 = moustache(totaux,"Total",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
tp5 = moustache(totaux,"Total",5,max = 1900, bars = c(1800,1600,0,0,0,0))
tp6 = moustache(totaux,"Total",6,max = 1200)
tp7 = moustache(totaux,"Total",7,max = 1500, bars = c(1450,1350,0,0,0,0))
tp8 = moustache(totaux,"Total",8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
tp9 = moustache(totaux,"Total",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
tt1 = tableau(1,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt2 = tableau(2,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt3 = tableau(3,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt4 = tableau(4,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt5 = tableau(5,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt6 = tableau(6,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt7 = tableau(7,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt8 = tableau(8,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tt9 = tableau(9,HESN = HESN_tot, TS = TS_tot, VIH = VIH_tot, HLA = HLA_tot, df_t = df_t_tot)
tg1 = grid.arrange(tp1,tableGrob(tt1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD14")
tg2 = grid.arrange(tp2,tableGrob(tt2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD16")
tg3 = grid.arrange(tp3,tableGrob(tt3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-DR")
tg4 = grid.arrange(tp4,tableGrob(tt4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-G")
tg5 = grid.arrange(tp5,tableGrob(tt5), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CCR2")
tg6 = grid.arrange(tp6,tableGrob(tt6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "ILT4")
tg7 = grid.arrange(tp7,tableGrob(tt7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11b")
tg8 = grid.arrange(tp8,tableGrob(tt8), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11c")
tg9 = grid.arrange(tp9,tableGrob(tt9), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "NR4A1")

#classiques
classiques = rbind(HESN_classiques[1:9,],TS_classiques[1:11,],VIH_classiques[1:9,],HLA_classiques[1:7,])
cp1 = moustache(classiques,"Classical",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
cp2 = moustache(classiques,"Classical",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
cp3 = moustache(classiques,"Classical",3,max = 2000, bars = c(1800,1600,1400,0,0,0))
cp4 = moustache(classiques,"Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
cp5 = moustache(classiques,"Classical",5,max = 1900, bars = c(1800,1600,0,0,0,0))
cp6 = moustache(classiques,"Classical",6,max = 1200)
cp7 = moustache(classiques,"Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
cp8 = moustache(classiques,"Classical",8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
cp9 = moustache(classiques,"Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
ct1 = tableau(1,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct2 = tableau(2,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct3 = tableau(3,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct4 = tableau(4,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct5 = tableau(5,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct6 = tableau(6,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct7 = tableau(7,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct8 = tableau(8,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
ct9 = tableau(9,HESN = HESN_classiques, TS = TS_classiques, VIH = VIH_classiques, HLA = HLA_classiques, df_t = df_t_classiques)
cg1 = grid.arrange(cp1,tableGrob(ct1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD14")
cg2 = grid.arrange(cp2,tableGrob(ct2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD16")
cg3 = grid.arrange(cp3,tableGrob(ct3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-DR")
cg4 = grid.arrange(cp4,tableGrob(ct4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-G")
cg5 = grid.arrange(cp5,tableGrob(ct5), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CCR2")
cg6 = grid.arrange(cp6,tableGrob(ct6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "ILT4")
cg7 = grid.arrange(cp7,tableGrob(ct7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11b")
cg8 = grid.arrange(cp8,tableGrob(ct8), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11c")
cg9 = grid.arrange(cp9,tableGrob(ct9), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "NR4A1")

#intermédiaires
intermediaires = rbind(HESN_intermediaires[1:9,],TS_intermediaires[1:11,],VIH_intermediaires[1:9,],HLA_intermediaires[1:7,])
ip1 = moustache(intermediaires,"Intermediate",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
ip2 = moustache(intermediaires,"Intermediate",2,max = 7500, bars = c(11500,11000,10500,0,0,0))
ip3 = moustache(intermediaires,"Intermediate",3,max = 2000, bars = c(1800,1600,1400,0,0,0))
ip4 = moustache(intermediaires,"Intermediate",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
ip5 = moustache(intermediaires,"Intermediate",5,max = 1900, bars = c(1800,1600,0,0,0,0))
ip6 = moustache(intermediaires,"Intermediate",6,max = 1200)
ip7 = moustache(intermediaires,"Intermediate",7,max = 1500, bars = c(1450,1350,0,0,0,0))
ip8 = moustache(intermediaires,"Intermediate",8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
ip9 = moustache(intermediaires,"Intermediate",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
it1 = tableau(1,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it2 = tableau(2,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it3 = tableau(3,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it4 = tableau(4,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it5 = tableau(5,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it6 = tableau(6,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it7 = tableau(7,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it8 = tableau(8,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
it9 = tableau(9,HESN = HESN_intermediaires, TS = TS_intermediaires, VIH = VIH_intermediaires, HLA = HLA_intermediaires, df_t = df_t_intermediaires)
ig1 = grid.arrange(ip1,tableGrob(it1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD14")
ig2 = grid.arrange(ip2,tableGrob(it2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD16")
ig3 = grid.arrange(ip3,tableGrob(it3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-DR")
ig4 = grid.arrange(ip4,tableGrob(it4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-G")
ig5 = grid.arrange(ip5,tableGrob(it5), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CCR2")
ig6 = grid.arrange(ip6,tableGrob(it6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "ILT4")
ig7 = grid.arrange(ip7,tableGrob(it7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11b")
ig8 = grid.arrange(ip8,tableGrob(it8), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11c")
ig9 = grid.arrange(ip9,tableGrob(it9), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "NR4A1")

#non-classiques
non_classiques = rbind(HESN_non_c[1:9,],TS_non_c[1:11,],VIH_non_c[1:9,],HLA_non_c[1:7,])
np1 = moustache(non_classiques,"Non-Classical",1,max = 1600, bars = c(1500,1400,1300,0,0,0))
np2 = moustache(non_classiques,"Non-Classical",2,max = 12000, bars = c(11500,11000,10500,0,0,0))
np3 = moustache(non_classiques,"Non-Classical",3,max = 2000, bars = c(1800,1600,1400,0,0,0))
np4 = moustache(non_classiques,"Non-Classical",4,min = -50,max = 260, bars = c(240,220,200,180,0,0))
np5 = moustache(non_classiques,"Non-Classical",5,max = 1900, bars = c(1800,1600,0,0,0,0))
np6 = moustache(non_classiques,"Non-Classical",6,max = 1200)
np7 = moustache(non_classiques,"Non-Classical",7,max = 1500, bars = c(1450,1350,0,0,0,0))
np8 = moustache(non_classiques,"Non-Classical",8,min = 1000,max = 7500, bars = c(7200,7000,6800,0,0,0))
np9 = moustache(non_classiques,"Non-Classical",9,min = 500,max = 2900, bars = c(2800,2700,2600,0,0,0))
nt1 = tableau(1,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt2 = tableau(2,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt3 = tableau(3,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt4 = tableau(4,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt5 = tableau(5,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt6 = tableau(6,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt7 = tableau(7,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt8 = tableau(8,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
nt9 = tableau(9,HESN = HESN_non_c, TS = TS_non_c, VIH = VIH_non_c, HLA = HLA_non_c, df_t = df_t_non_c)
ng1 = grid.arrange(np1,tableGrob(nt1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD14")
ng2 = grid.arrange(np2,tableGrob(nt2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD16")
ng3 = grid.arrange(np3,tableGrob(nt3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-DR")
ng4 = grid.arrange(np4,tableGrob(nt4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-G")
ng5 = grid.arrange(np5,tableGrob(nt5), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CCR2")
ng6 = grid.arrange(np6,tableGrob(nt6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "ILT4")
ng7 = grid.arrange(np7,tableGrob(nt7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11b")
ng8 = grid.arrange(np8,tableGrob(nt8), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11c")
ng9 = grid.arrange(np9,tableGrob(nt9), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "NR4A1")

#pDC
#pDC = rbind(HESN_pDC[1:9,],TS_pDC[1:11,],VIH_pDC[1:9,],HLA_pDC[1:7,])
#pp1 = moustache(pDC,"pDC",1,min = 300,max = 475)
#pp2 = moustache(pDC,"pDC",2,min = 20,max = 80)
#pp3 = moustache(pDC,"pDC",3,min = 25,max = 125)
#pp4 = moustache(pDC,"pDC",4,max = 350)
#pp5 = moustache(pDC,"pDC",5,min = 25,max = 45)
#pp6 = moustache(pDC,"pDC",6,min = 1000,max = 1800)
#pp7 = moustache(pDC,"pDC",7,min = 300,max = 750)
#pt1 = tableau(1,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt2 = tableau(2,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt3 = tableau(3,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt4 = tableau(4,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt5 = tableau(5,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt6 = tableau(6,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pt7 = tableau(7,HESN = HESN_pDC, TS = TS_pDC, VIH = VIH_pDC, HLA = HLA_pDC, df_t = df_t_pDC)
#pg1 = grid.arrange(pp1,tableGrob(pt1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-DR")
#pg2 = grid.arrange(pp2,tableGrob(pt2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "HLA-G")
#pg3 = grid.arrange(pp3,tableGrob(pt3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CCR2")
#pg4 = grid.arrange(pp4,tableGrob(pt4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "ILT4")
#pg5 = grid.arrange(pp5,tableGrob(pt5), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11b")
#pg6 = grid.arrange(pp6,tableGrob(pt6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "CD11c")
#pg7 = grid.arrange(pp7,tableGrob(pt7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "NR4A1")

#Layout tous dans un
ggsave("nucleaire_totaux_comparison_fluorescence.png", grid.arrange(tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9, nrow = 3, top = "Comparaison des fluorescence des marqueurs de monocytes totaux par populations"), width = 15, height = 8.5, units = "in")
ggsave("nucleaire_classiques_comparison_fluorescence.png", grid.arrange(cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp9, nrow = 3, top = "Comparaison des fluorescence des marqueurs de monocytes classiques par populations"), width = 15, height = 8.5, units = "in")
ggsave("nucleaire_intermediaires_comparison_fluorescence.png", grid.arrange(ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8,ip9, nrow = 3, top = "Comparaison des fluorescence des marqueurs de monocytes intermediaires par populations"), width = 15, height = 8.5, units = "in")
ggsave("nucleaire_non-classiques_comparison_fluorescence.png", grid.arrange(np1,np2,np3,np4,np5,np6,np7,np8,np9, nrow = 3, top = "Comparaison des fluorescence des marqueurs de monocytes non-classiques par populations"), width = 15, height = 8.5, units = "in")
#ggsave("nucleaire_pDC_comparison_fluorescence.png", grid.arrange(pp1,pp2,pp3,pp4,pp5,pp6,pp7, nrow = 2, top = "Comparaison des fluorescence des marqueurs de pDC par populations"), width = 15, height = 8.5, units = "in")

#Layout graph + tableau par page
ggsave("nucleaire_totaux_fluorescence_marqueurs_p-values.pdf", marrangeGrob(list(tg1,tg2,tg3,tg4,tg5,tg6,tg7,tg8,tg9), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")
ggsave("nucleaire_classiques_fluorescence_marqueurs_p-values.pdf", marrangeGrob(list(cg1,cg2,cg3,cg4,cg5,cg6,cg7,cg8,cg9), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")
ggsave("nucleaire_intermediaires_fluorescence_marqueurs_p-values.pdf", marrangeGrob(list(ig1,ig2,ig3,ig4,ig5,ig6,ig7,ig8,ig9), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")
ggsave("nucleaire_non-classiques_fluorescence_marqueurs_p-values.pdf", marrangeGrob(list(ng1,ng2,ng3,ng4,ng5,ng6,ng7,ng8,ng9), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")
#ggsave("nucleaire_pDC_fluorescence_marqueurs_p-values.pdf", marrangeGrob(list(pg1,pg2,pg3,pg4,pg5,pg6,pg7), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")

#Layout 3 graph + tableau sous-populations de monocytes par marqueur
lg1 = grid.arrange(cp1,ip1,np1, tableGrob(ct1), tableGrob(it1), tableGrob(nt1), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "CD14")
lg2 = grid.arrange(cp2,ip2,np2, tableGrob(ct2), tableGrob(it2), tableGrob(nt2), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "CD16")
lg3 = grid.arrange(cp3,ip3,np3, tableGrob(ct3), tableGrob(it3), tableGrob(nt3), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "HLA-DR")
lg4 = grid.arrange(cp4,ip4,np4, tableGrob(ct4), tableGrob(it4), tableGrob(nt4), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "HLA-G")
lg5 = grid.arrange(cp5,ip5,np5, tableGrob(ct5), tableGrob(it5), tableGrob(nt5), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "CCR2")
lg6 = grid.arrange(cp6,ip6,np6, tableGrob(ct6), tableGrob(it6), tableGrob(nt6), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "ILT4")
lg7 = grid.arrange(cp7,ip7,np7, tableGrob(ct7), tableGrob(it7), tableGrob(nt7), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "CD11b")
lg8 = grid.arrange(cp8,ip8,np8, tableGrob(ct8), tableGrob(it8), tableGrob(nt8), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "CD11c")
lg9 = grid.arrange(cp9,ip9,np9, tableGrob(ct9), tableGrob(it9), tableGrob(nt9), nrow = 2, ncol = 3, widths = c(1,1,1), clip = F, top = "NR4A1")
ggsave("nucleaire_comparaison_sous-pops.pdf", marrangeGrob(list(lg1,lg2,lg3,lg4,lg5,lg6,lg7,lg8,lg9), nrow = 1, ncol = 1), width = 15, height = 8.5, units = "in")

#Layout Totaux + sous-pops + stats ds tableau
lg1 = grid.arrange(tp1,cp1,ip1,np1, tableGrob(tt1), tableGrob(ct1), tableGrob(it1), tableGrob(nt1), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD14")
lg2 = grid.arrange(tp2,cp2,ip2,np2, tableGrob(tt2), tableGrob(ct2), tableGrob(it2), tableGrob(nt2), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD16")
lg3 = grid.arrange(tp3,cp3,ip3,np3, tableGrob(tt3), tableGrob(ct3), tableGrob(it3), tableGrob(nt3), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "HLA-DR")
lg4 = grid.arrange(tp4,cp4,ip4,np4, tableGrob(tt4), tableGrob(ct4), tableGrob(it4), tableGrob(nt4), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "HLA-G")
lg5 = grid.arrange(tp5,cp5,ip5,np5, tableGrob(tt5), tableGrob(ct5), tableGrob(it5), tableGrob(nt5), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CCR2")
lg6 = grid.arrange(tp6,cp6,ip6,np6, tableGrob(tt6), tableGrob(ct6), tableGrob(it6), tableGrob(nt6), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "ILT4")
lg7 = grid.arrange(tp7,cp7,ip7,np7, tableGrob(tt7), tableGrob(ct7), tableGrob(it7), tableGrob(nt7), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD11b")
lg8 = grid.arrange(tp8,cp8,ip8,np8, tableGrob(tt8), tableGrob(ct8), tableGrob(it8), tableGrob(nt8), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "CD11c")
lg9 = grid.arrange(tp9,cp9,ip9,np9, tableGrob(tt9), tableGrob(ct9), tableGrob(it9), tableGrob(nt9), nrow = 2, ncol = 4, widths = c(1,1,1,1), clip = F, top = "NR4A1")
ggsave("nucleaire_layout_monocytes.pdf", marrangeGrob(list(lg1,lg2,lg3,lg4,lg5,lg6,lg7,lg8,lg9), nrow = 1, ncol = 1), width = 20, height = 8.5, units = "in")


### Régression âge des TS/HESN

#Import data
nucleaire_totaux_temps_travail <- read_csv("nucleaire_totaux_temps_travail.csv")
nucleaire_classiques_temps_travail <- read_csv("nucleaire_classiques_temps_travail.csv")
nucleaire_intermediaires_temps_travail <- read_csv("nucleaire_intermediaires_temps_travail.csv")
nucleaire_non_classiques_temps_travail <- read_csv("nucleaire_non-classiques_temps_travail.csv")

#Colnames
colnames(nucleaire_totaux_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nucleaire_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nucleaire_intermediaires_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
colnames(nucleaire_non_classiques_temps_travail) = c("participantes","temps_travail","CD14","CD16","HLADR","HLAG","CCR2","ILT4","CD11b","CD11c","NR4A1")
#Graphs + layouts

#Total
tg1 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=CD14)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD14", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg2 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=CD16)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD16", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg3 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=HLADR)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-DR", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg4 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=HLAG)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-G", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg5 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=CCR2)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CCR2", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg6 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=ILT4)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("ILT4", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg7 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=CD11b)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11b", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg8 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=CD11c)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11c", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
tg9 = ggplot(nucleaire_totaux_temps_travail, aes(x=temps_travail, y=NR4A1)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("NR4A1", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
layout_tot = grid.arrange(tg1,tg2,tg3,tg4,tg5,tg6,tg7,tg8,tg9, nrow = 3, ncol = 3, widths = c(1,1,1), clip = F, top = "Total")

#Classiques
cg1 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=CD14)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD14", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg2 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=CD16)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD16", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg3 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=HLADR)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-DR", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg4 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=HLAG)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-G", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg5 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=CCR2)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CCR2", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg6 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=ILT4)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("ILT4", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg7 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=CD11b)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11b", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg8 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=CD11c)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11c", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cg9 = ggplot(nucleaire_classiques_temps_travail, aes(x=temps_travail, y=NR4A1)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("NR4A1", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
layout_classiques = grid.arrange(cg1,cg2,cg3,cg4,cg5,cg6,cg7,cg8,cg9, nrow = 3, ncol = 3, widths = c(1,1,1), clip = F, top = "Classical")

#Intermédiaires
ig1 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=CD14)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD14", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig2 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=CD16)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD16", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig3 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=HLADR)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-DR", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig4 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=HLAG)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-G", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig5 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=CCR2)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CCR2", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig6 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=ILT4)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("ILT4", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig7 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=CD11b)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11b", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig8 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=CD11c)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11c", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ig9 = ggplot(nucleaire_intermediaires_temps_travail, aes(x=temps_travail, y=NR4A1)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("NR4A1", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
layout_intermediaires = grid.arrange(ig1,ig2,ig3,ig4,ig5,ig6,ig7,ig8,ig9, nrow = 3, ncol = 3, widths = c(1,1,1), clip = F, top = "Intermediate")

#Non-classiques
ng1 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=CD14)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD14", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng2 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=CD16)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD16", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng3 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=HLADR)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-DR", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng4 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=HLAG)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("HLA-G", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng5 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=CCR2)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CCR2", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng6 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=ILT4)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("ILT4", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng7 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=CD11b)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11b", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng8 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=CD11c)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("CD11c", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ng9 = ggplot(nucleaire_non_classiques_temps_travail, aes(x=temps_travail, y=NR4A1)) + geom_point() + geom_smooth(method = "lm") + coord_cartesian(xlim = c(0, 20)) + xlab("Time of sex work") + ylab(paste("NR4A1", "expressing cells (GeoMFI)", sep = " ")) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
layout_non_classiques = grid.arrange(ng1,ng2,ng3,ng4,ng5,ng6,ng7,ng8,ng9, nrow = 3, ncol = 3, widths = c(1,1,1), clip = F, top = "Non-classical")

ggsave("nucleaire_totaux_age_regression.png", layout_tot,width = 15, height = 8.5, units = "in")
ggsave("nucleaire_classiques_age_regression.png", layout_classiques,width = 15, height = 8.5, units = "in")
ggsave("nucleaire_intermediaires_age_regression.png", layout_intermediaires,width = 15, height = 8.5, units = "in")
ggsave("nucleaire_non_classiques_age_regression.png", layout_non_classiques,width = 15, height = 8.5, units = "in")

ggsave("nucleaire_age_regression.pdf", marrangeGrob(list(layout_tot,layout_classiques,layout_intermediaires,layout_non_classiques), nrow = 1, ncol = 1), width = 10, height = 8.5, units = "in")
