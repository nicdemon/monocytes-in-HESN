#WORK_DIR + dependencies
setwd("../data/cytokines")
library(data.table)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(gtable)
library(digest)
library(ggpubr)

#Import data
Cytokines <- as.data.frame(read_csv("Cytokines.csv"))
rownames(Cytokines) = Cytokines$Participantes
Cytokines$Participantes = NULL

###Visualisation
#Dotplot
par(mfrow=c(5,5))
for (i in c(2:26)){
  plot(c(1:40),Cytokines[,i], ylab = colnames(Cytokines)[i], xlab = "")
}
mtext("Concentrations de cytokines pré-traitement", side=3, outer=TRUE, line=-3)
par(mfrow=c(1,1))

#Split par populations
#Petits df par pops
HESN = as.data.frame(Cytokines[1:12,])
rownames(HESN) = rownames(Cytokines)[1:12]
HESN$Participantes = NULL
HESN$Status_VIH = NULL
TS = as.data.frame(Cytokines[13:24,])
rownames(TS) = rownames(Cytokines)[13:24]
TS$Participantes = NULL
TS$Status_VIH = NULL
VIH = as.data.frame(Cytokines[25:35,])
rownames(VIH) = rownames(Cytokines)[25:35]
VIH$Participantes = NULL
VIH$Status_VIH = NULL
HLA = as.data.frame(Cytokines[36:40,])
rownames(HLA) = rownames(Cytokines)[36:40]
HLA$Participantes = NULL
HLA$Status_VIH = NULL

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

remove_outliers <- function(x, na.rm = TRUE) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

HESN = outliers(HESN)
TS = outliers(TS)
VIH = outliers(VIH)
HLA = outliers(HLA)

#Normalité
normalite = function(df){
  maxrow = length(rownames(df))+1
  for(i in colnames(df)){
    tryCatch({df[maxrow,i] = shapiro.test(na.omit(df[,i]))[2]},error=function(cond){return(NULL)})
  }
  rownames(df)[length(rownames(df))] = 'Normalite_p-value'
  return(df)
}

HESN = normalite(HESN)
TS = normalite(TS)
VIH = normalite(VIH)
HLA = normalite(HLA)

#T-test entre pops par cytokine
df_t = data.frame(matrix(ncol = 25, nrow = 6),row.names = c("HESN_HLA","HESN_TS","HESN_VIH","TS_HLA","TS_VIH","VIH_HLA"))
colnames(df_t) = colnames(HESN)

for(i in colnames(df_t)){
  tryCatch({df_t[1,i] = t.test(x = HESN[1:12,i], y = HLA[1:5,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t[2,i] = t.test(x = HESN[1:12,i], y = TS[1:12,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t[3,i] = t.test(x = HESN[1:12,i], y = VIH[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t[4,i] = t.test(x = TS[1:12,i], y = HLA[1:5,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t[5,i] = t.test(x = TS[1:12,i], y = VIH[1:11,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
  tryCatch({df_t[6,i] = t.test(x = VIH[1:11,i], y = HLA[1:5,i], na.action = na.omit)$p.value},error=function(cond){return(NULL)})
}

#Tableau recap
write.csv(df_t, "T-test_results.csv")
write.csv(HESN, "HESN_normalite.csv")
write.csv(TS, "TS_normalite.csv")
write.csv(VIH, "VIH_normalite.csv")
write.csv(HLA, "HLA_normalite.csv")

###Graphs cytokines par populations
#Dotplot par population
dotplot = function(df, name="Inconnu"){
  par(mfrow=c(5,5))
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
  mtext(paste("Concentration des cytokines pour la population", name), side=3, outer=TRUE, line=-3)
  par(mfrow=c(1,1))
}

dotplot(HESN, "HESN")
dotplot(TS, "TS+ VIH-")
dotplot(VIH, "TS+ VIH+")
dotplot(HLA, "HLA")

#Boxplot par cytokine
moustache = function(df,i,min = 0, max = 1000, bars = c(0,0,0,0,0,0)){
  tmp = data.frame()
  tmp[1:12,1] = "HESN"
  tmp[13:24,1] = "CSWs HIV-"
  tmp[25:35,1] = "CSWs HIV+"
  tmp[36:40,1] = "Non CSWs HIV-"
  tmp[,2] = df[,i]
  tmp[1:12,3] = mean(tmp[1:12,2],na.rm=TRUE)
  tmp[13:24,3] = mean(tmp[13:24,2],na.rm=TRUE)
  tmp[25:35,3] = mean(tmp[25:35,2],na.rm=TRUE)
  tmp[36:40,3] = mean(tmp[36:40,2],na.rm=TRUE)
  tmp[1:12,4] = sd(tmp[1:12,2],na.rm=TRUE)
  tmp[13:24,4] = sd(tmp[13:24,2],na.rm=TRUE)
  tmp[25:35,4] = sd(tmp[25:35,2],na.rm=TRUE)
  tmp[36:40,4] = sd(tmp[36:40,2],na.rm=TRUE)
  name = colnames(df)[i]
  colnames(tmp) = c("populations","cytokine","moy","std")
  matrice_tmp = compare_means(cytokine ~ populations, data=na.omit(tmp), method = "t.test")
  comparison = list()
  for (i in seq(1,length(matrice_tmp$p.signif))){
    if (matrice_tmp$p.signif[i] != "ns"){
      comparison[[length(comparison)+1]] = c(matrice_tmp$group1[i],matrice_tmp$group2[i])
    }
  }
  tmp %>%
    arrange(cytokine) %>%
    mutate(populations = factor(populations, levels = c("Non CSWs HIV-","CSWs HIV-","HESN","CSWs HIV+"))) %>%
    ggplot(aes(x=populations, y=cytokine, group=populations)) + 
    geom_point(aes(x=populations, y=cytokine, shape = populations, size = 5), binaxis = "y", stackdir = "center") +
    scale_shape_manual(values = c(15:18)) +
    geom_errorbar(aes(ymin=moy-std,ymax=moy+std), width=0.2, position=position_dodge(0.1),size = 0.5) +
    geom_errorbar(aes(ymin=moy,ymax=moy), width=0.2, position=position_dodge(0.1),size = 0.5) +
    stat_compare_means(comparisons = comparison, method = "t.test", label = "p.signif", hide.ns = TRUE, na.rm = TRUE, label.y = bars) +
    coord_cartesian(ylim = c(min, max)) +
    xlab("") +
    ylab(paste(name, "concentration", sep = " ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, vjust = 0.5))
}

tableau = function(i){
  df = data.frame(matrix(nrow = 10, ncol = 1))
  row.names(df) = c("Normalité HESN", "Normalité TS+VIH-", "Normalité TS+VIH+", "Normalité HLA", "T-test HESN vs HLA", "T-test HESN vs TS+VIH-", "T-test HESN vs TS+VIH+", "T-test TS+VIH- vs HLA", "T-test TS+VIH- vs TS+VIH+", "T-test TS+VIH+ vs HLA")
  colnames(df) = c("p_value")
  df$p_value[1] = HESN[13,i]
  df$p_value[2] = TS[13,i]
  df$p_value[3] = VIH[12,i]
  df$p_value[4] = HLA[6,i]
  df$p_value[5:10] = df_t[1:6,i]
  return(df)
}

Cytokines = rbind(HESN[1:12,], TS[1:12,],VIH[1:11,], HLA[1:5,])
p1 = moustache(Cytokines,1, max = 30, bars = c(29,27,25,0,0,0))
p2 = moustache(Cytokines,2, max = 0.8, bars = c(0.78,0.74,0.7,0.66,0,0))
p3 = moustache(Cytokines,3, max = 40)
p4 = moustache(Cytokines,4, max = 4.5)
p6 = moustache(Cytokines,6, max = 150, bars = c(145,135,0,0,0,0))
p7 = moustache(Cytokines,7, max = 3.5, bars = c(3.2,3,2.8,2.6,0,0))
p8 = moustache(Cytokines,8, max = 160, bars = c(150,140,130,120,0,0))
p9 = moustache(Cytokines,9, max = 1.8, bars = c(1.7,1.65,1.6,1.55,0,0))
p12 = moustache(Cytokines,12, max = 250)
p13 = moustache(Cytokines,13, max = 45)
p14 = moustache(Cytokines,14, max = 12, bars = c(11,0,0,0,0,0))
p15 = moustache(Cytokines,15, max = 6)
p16 = moustache(Cytokines,16, max = 25, bars = c(24,0,0,0,0,0))
p17 = moustache(Cytokines,17, max = 260)
p18 = moustache(Cytokines,18, max = 3.8)
p19 = moustache(Cytokines,19, max = 180, bars = c(170,160,150,0,0,0))
p20 = moustache(Cytokines,20, max = 180, bars = c(170,160,150,0,0,0))
p21 = moustache(Cytokines,21, max = 600)
p22 = moustache(Cytokines,22, max = 300, bars = c(290,270,250,230,0,0))
p23 = moustache(Cytokines,23, max = 350, bars = c(340,320,300,0,0,0))
p24 = moustache(Cytokines,24, max = 45)
p25 = moustache(Cytokines,25, max = 45, bars = c(44,42,40,0,0,0))
t1 = tableau(1)
t2 = tableau(2)
t3 = tableau(3)
t4 = tableau(4)
t6 = tableau(6)
t7 = tableau(7)
t8 = tableau(8)
t9 = tableau(9)
t12 = tableau(12)
t13 = tableau(13)
t14 = tableau(14)
t15 = tableau(15)
t16 = tableau(16)
t17 = tableau(17)
t18 = tableau(18)
t19 = tableau(19)
t20 = tableau(20)
t21 = tableau(21)
t22 = tableau(22)
t23 = tableau(23)
t24 = tableau(24)
t25 = tableau(25)
g1 = grid.arrange(p1,tableGrob(t1), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "Eotaxin")
g2 = grid.arrange(p2,tableGrob(t2), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "GM-CSF")
g3 = grid.arrange(p3,tableGrob(t3), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IFN-a")
g4 = grid.arrange(p4,tableGrob(t4), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IFN-G")
g6 = grid.arrange(p6,tableGrob(t6), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-1RA")
g7 = grid.arrange(p7,tableGrob(t7), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-10")
g8 = grid.arrange(p8,tableGrob(t8), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-12")
g9 = grid.arrange(p9,tableGrob(t9), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-13")
g12 = grid.arrange(p12,tableGrob(t12), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-2")
g13 = grid.arrange(p13,tableGrob(t13), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-2R")
g14 = grid.arrange(p14,tableGrob(t14), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-4")
g15 = grid.arrange(p15,tableGrob(t15), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-6")
g16 = grid.arrange(p16,tableGrob(t16), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-7")
g17 = grid.arrange(p17,tableGrob(t17), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IL-8")
g18 = grid.arrange(p18,tableGrob(t18), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "IP-10")
g19 = grid.arrange(p19,tableGrob(t19), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "MCP-1")
g20 = grid.arrange(p20,tableGrob(t20), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "MIG")
g21 = grid.arrange(p21,tableGrob(t21), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "MIP-1a")
g22 = grid.arrange(p22,tableGrob(t22), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "MIP-1b")
g23 = grid.arrange(p23,tableGrob(t23), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "RANTES")
g24 = grid.arrange(p24,tableGrob(t24), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "TNF-a")
g25 = grid.arrange(p25,tableGrob(t25), nrow = 1, ncol = 2, widths = c(1,1.5), clip = FALSE, top = "APRIL")
grid.arrange(p1,p2,p3,p4,p6,p7,p8,p9,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25, nrow = 5, top = "Comparaison des concentrations de cytokines par populations")
ml = marrangeGrob(list(g1,g2,g3,g4,g6,g7,g8,g9,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g25), nrow = 1, ncol = 1)
ggsave("Concentrations_cytokines_p-values.pdf",ml)
