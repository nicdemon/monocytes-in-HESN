# WORK_DIR + dependencies
setwd("D:/Google Drive Laurence/Maitrise/LABO/R_stats/Séroconversion")
library(data.table)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(gtable)
library(digest)
library(ggpubr)
library(cowplot)

## Functions
# Barplot par cytokine
bars_serum = function(df,i,min = 0, max = 1000, labels_height = -1.5, label_size = 2){
  name = colnames(df[,i])
  if(i == 6){
    label_y = paste(paste("IFN-","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 7){
    label_y = paste(paste("IFN-","\U03B3", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 8){
    label_y = paste(paste("IL-1","\U03B2", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 14){
    label_y = paste(paste("IL-17","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 25){
    label_y = paste(paste("MIP-1","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 26){
    label_y = paste(paste("MIP-1","\U03B2", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 28){
    label_y = paste(paste("TNF-","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 29){
    label_y = paste(paste("MIP-3","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 31){
    label_y = "IgG1 gp41 (MFI)"
  }else{
    label_y = paste(name, "(pg/ml)", sep = " ")
  }
  
  #Early HESN
  tmp_early = data.frame(df)[1:4,1:3]
  tmp_early[,4] = df[1:4,i]
  tmp_early = data.frame(na.omit(tmp_early))
  colnames(tmp_early) = c("Participantes","Visite","Population","Cytokine")
  tmp_early$Visites = as.character(seq(1,length(tmp_early$Participantes),1))
  early = tmp_early %>%
    arrange(Visites) %>%
    mutate(Visite = factor(Visite, levels = c("0.5","1.5","0.25","1.25"))) %>%
    ggplot(aes(fill = Visites, x = Participantes, y = Cytokine, color = Participantes)) +
    geom_bar(position = position_dodge2(width = 0.9, preserve = "single"), stat = "identity") +
    scale_fill_manual(values = c("#009933","#009933","#006633","#006633")) +
    scale_color_manual(values = c("#000000","#000000","#000000","#000000")) +
    geom_text(aes(label = Visite, y = labels_height), position = position_dodge2(width = 0.9, preserve = "single"), size = label_size) +
    xlab("Early converter") +
    ylab(label_y) +
    coord_cartesian(ylim = c(min, max)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 10, color = "black"))
  
  #HESN
  tmp_HESN = data.frame(df)[5:15,1:3]
  tmp_HESN[,4] = df[5:15,i]
  tmp_HESN = data.frame(na.omit(tmp_HESN))
  colnames(tmp_HESN) = c("Participantes","Visite","Population","Cytokine")
  tmp_HESN$Visites = as.character(seq(5, 4+length(tmp_HESN$Participantes),1))
  HESN = tmp_HESN %>%
    arrange(Visites) %>%
    mutate(Visite = factor(Visite, levels = c("8.25","8.5","8.75","6","9.25","9.5","10","5","8","8.75*","9"))) %>%
    ggplot(aes(fill = Visite, x = Participantes, y = Cytokine, color = Participantes)) +
    geom_bar(position = position_dodge2(width = 0.9, preserve = "single"), stat = "identity") +
    scale_fill_manual(values = c("#66CCFF","#66CCFF","#66CCFF","#3399CC","#3399CC","#3399CC","#3399CC","#003366","#003366","#003366","#003366")) +
    scale_color_manual(values = c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")) +
    geom_text(aes(label = Visite, y = labels_height), position = position_dodge2(width = 0.9, preserve = "single"), size = label_size) +
    xlab("HESN converter") +
    coord_cartesian(ylim = c(min, max)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text = element_text(size = 10, color = "black"))
  
  #Graph global
  graph = plot_grid(early, HESN, nrow = 1, ncol = 2, align = "h", rel_widths = c(2,3))
  return(graph)
}

bars_CVL = function(df,i,min = 0, max = 1000, labels_height = -1.5, label_size = 2){
  name = colnames(df[,i])
  if(i == 6){
    label_y = paste(paste("IFN-","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 7){
    label_y = paste(paste("IFN-","\U03B3", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 8){
    label_y = paste(paste("IL-1","\U03B2", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 14){
    label_y = paste(paste("IL-17","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 25){
    label_y = paste(paste("MIP-1","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 26){
    label_y = paste(paste("MIP-1","\U03B2", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 28){
    label_y = paste(paste("TNF-","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 29){
    label_y = paste(paste("MIP-3","\U03B1", sep = ""), "(pg/ml)", sep = " ")
  }else if(i == 31){
    label_y = "IgG1 gp41 (MFI)"
  }else{
    label_y = paste(name, "(pg/ml)", sep = " ")
  }
  
  #Early HESN
  tmp_early = data.frame(df)[1:6,1:3]
  tmp_early[,4] = df[1:6,i]
  tmp_early = data.frame(na.omit(tmp_early))
  colnames(tmp_early) = c("Participantes","Visite","Population","Cytokine")
  tmp_early$Visites = as.character(seq(1,length(tmp_early$Participantes),1))
  fill_early = c("#009933","#009933","#009933","#006633","#006633","#006633")
  color_early = c("#000000","#000000","#000000","#000000","#000000","#000000")
  early = tmp_early %>%
    arrange(Visites) %>%
    mutate(Visite = factor(Visite, levels = c("0.5","1","1.5","0.25","0.75","1.25"))) %>%
    ggplot(aes(fill = Visites, x = Participantes, y = Cytokine, color = Participantes)) +
    geom_bar(position = position_dodge2(width = 1, preserve = "single"), stat = "identity") +
    scale_fill_manual(values = fill_early) +
    scale_color_manual(values = color_early) +
    geom_text(aes(label = Visite, y = labels_height), position = position_dodge2(width = 0.9, preserve = "single"), size = label_size) +
    xlab("Early converter") +
    ylab(label_y) +
    coord_cartesian(ylim = c(min, max)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 10, color = "black"))
  
  #HESN
  tmp_HESN = data.frame(df)[7:18,1:3]
  tmp_HESN[,4] = df[7:18,i]
  tmp_HESN = data.frame(na.omit(tmp_HESN))
  colnames(tmp_HESN) = c("Participantes","Visite","Population","Cytokine")
  tmp_HESN$Visites = as.character(seq(7,(6+length(tmp_HESN$Participantes)),1))
  fill_HESN = c("#66CCFF","#66CCFF","#66CCFF","#66CCFF","#3399CC","#3399CC","#3399CC","#3399CC","#003366","#003366","#003366","#003366")
  color_HESN = c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")
  HESN = tmp_HESN %>%
    arrange(Visites) %>%
    mutate(Visite = factor(Visite, levels = c("5","8.25","8.5","8.75","6","9.25","9.5","10","5*","8","8.75*","9"))) %>%
    ggplot(aes(fill = Visite, x = Participantes, y = Cytokine, color = Participantes)) +
    geom_bar(position = position_dodge2(width = 1, preserve = "single"), stat = "identity") +
    scale_fill_manual(values = fill_HESN) +
    scale_color_manual(values = color_HESN) +
    geom_text(aes(label = Visite, y = labels_height), position = position_dodge2(width = 0.9, preserve = "single"), size = label_size) +
    xlab("HESN converter") +
    coord_cartesian(ylim = c(min, max)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text = element_text(size = 10, color = "black"))
  
  #Graph global
  graph = plot_grid(early, HESN, nrow = 1, ncol = 2, align = "h", axis = "t", rel_widths = c(2,3))
  return(graph)
}

## Import data
serum <- read_csv("D:/Google Drive Laurence/Maitrise/LABO/R_stats/Séroconversion/serum.csv", col_types = cols(Participantes = col_character()))
CVL <- read_csv("D:/Google Drive Laurence/Maitrise/LABO/R_stats/Séroconversion/CVL.csv", col_types = cols(Participantes = col_character()))

## Visualisation pre-stats
# Dotplot
pairs(serum[,4:31])
pairs(CVL[,4:31])

# Bandes de données
# Timeline par participante en mois
serum_eotaxin = bars_serum(serum,4, max = 65) #Eotaxin
serum_GMCSF = bars_serum(serum,5, max = 1)#GM-CSF
serum_IFNa = bars_serum(serum,6, max = 10)#IFN-a
serum_IFNg = bars_serum(serum,7, max = 1)#IFN-g
serum_IL1b = bars_serum(serum,8, max = 1.75, labels_height = -0.1)#IL-1b
serum_IL1RA = bars_serum(serum,9, max = 400)#IL-1RA
serum_IL10 = bars_serum(serum,10, max = 25)#IL-10
serum_IL12 = bars_serum(serum,11, max = 300)#IL-12
serum_IL13 = bars_serum(serum,12, max = 15)#IL-13
serum_IL15 = bars_serum(serum,13, max = 45)#IL-15
serum_IL17a = bars_serum(serum,14)#IL-17a
serum_IL2 = bars_serum(serum,15, max = 15)#IL-2
serum_IL2R = bars_serum(serum,16, max = 275)#IL-2R
serum_IL4 = bars_serum(serum,17, max = 30)#IL-4
serum_IL5 = bars_serum(serum,18)#IL-5
serum_IL6 = bars_serum(serum,19, max = 150)#IL-6
serum_IL7 = bars_serum(serum,20, max = 20)#IL-7
serum_IL8 = bars_serum(serum,21, max = 200)#IL-8
serum_IP10 = bars_serum(serum,22, max = 15, labels_height = -0.3, label_size = 2.5)#IP-10
serum_MCP1 = bars_serum(serum,23, max = 350, labels_height = -7, label_size = 2.5)#MCP-1
serum_MIG = bars_serum(serum,24, max = 900, labels_height = -20, label_size = 2.5)#MIG
serum_MIP1a = bars_serum(serum,25, max = 120, labels_height = -2.5, label_size = 2.5)#MIP-1a
serum_MIP1b = bars_serum(serum,26, max = 250, labels_height = -5, label_size = 2.5)#MIP-1b
serum_RANTES = bars_serum(serum,27, max = 1500)#RANTES
serum_TNFa = bars_serum(serum,28, max = 60)#TNF-a
serum_MIP3a = bars_serum(serum,29, max = 4500)#MIP-3a
serum_BAFF = bars_serum(serum,30, max = 36, labels_height = -0.7, label_size = 3)#BAFF
serum_gp41 = bars_serum(serum,31, max = 250000000, labels_height = -5000000, label_size = 3)#gp41

CVL_eotaxin = bars_CVL(CVL,4, max = 5) #Eotaxin
CVL_GMCSF = bars_CVL(CVL,5)#GM-CSF
CVL_IFNa = bars_CVL(CVL,6, max = 35, labels_height = -0.5)#IFN-a
ggsave("CVL_IFNa.png", width = 5, height = 5, units = "in")
CVL_IFNg = bars_CVL(CVL,7, max = 20, labels_height = -0.1)#IFN-g
CVL_IL1b = bars_CVL(CVL,8, max = 250)#IL-1b
CVL_IL1RA = bars_CVL(CVL,9, max = 670000)#IL-1RA
CVL_IL10 = bars_CVL(CVL,10, max = 2, labels_height = -0.1)#IL-10
CVL_IL12 = bars_CVL(CVL,11, max = 30)#IL-12
CVL_IL13 = bars_CVL(CVL,12)#IL-13
CVL_IL15 = bars_CVL(CVL,13)#IL-15
CVL_IL17a = bars_CVL(CVL,14)#IL-17a
CVL_IL2 = bars_CVL(CVL,15)#IL-2
CVL_IL2R = bars_CVL(CVL,16)#IL-2R
CVL_IL4 = bars_CVL(CVL,17, max = 30)#IL-4
CVL_IL5 = bars_CVL(CVL,18)#IL-5
CVL_IL6 = bars_CVL(CVL,19, max = 100)#IL-6
CVL_IL7 = bars_CVL(CVL,20)#IL-7
CVL_IL8 = bars_CVL(CVL,21, max = 20000)#IL-8
CVL_IP10 = bars_CVL(CVL,22, max = 35, labels_height = -0.8, label_size = 2.5)#IP-10
CVL_MCP1 = bars_CVL(CVL,23, max = 205, labels_height = -4, label_size = 2.5)#MCP-1
CVL_MIG = bars_CVL(CVL,24, max = 350, labels_height = -7, label_size = 2.5)#MIG
CVL_MIP1a = bars_CVL(CVL,25, max = 150, labels_height = -3, label_size = 2.5)#MIP-1a
CVL_MIP1b = bars_CVL(CVL,26, max = 510, labels_height = -10, label_size = 2.5)#MIP-1b
CVL_RANTES = bars_CVL(CVL,27, max = 1700)#RANTES
CVL_TNFa = bars_CVL(CVL,28, max = 20, labels_height = -0.1)#TNF-a
CVL_MIP3a = bars_CVL(CVL,29, max = 11000)#MIP-3a
CVL_BAFF = bars_CVL(CVL,30, max = 90, labels_height = -1.5, label_size = 3)#BAFF
CVL_gp41 = bars_CVL(CVL,31, max = 120000, labels_height = -2400, label_size = 3)#gp41


### Layouts

## CVL/Serum seuls

# Chemokines
#MIP1a
#MCP1
#MIP1b
#MIG
#IP10
#RANTES
chimiokines_CVL = plot_grid(CVL_MIP1a,CVL_MCP1,CVL_MIP1b,CVL_MIG,CVL_IP10,CVL_RANTES, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("chimiokines_CVL.png", width = 8.5, height = 11, units = "in")
chimiokines_serum = plot_grid(serum_MIP1a,serum_MCP1,serum_MIP1b,serum_MIG,serum_IP10,serum_RANTES, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("chimiokines_serum.png", width = 8.5, height = 11, units = "in")
# Cytokines
#IL6
#IL10
#TNFa
#IFNg
#IFNa
cytokines_CVL = plot_grid(CVL_IL6,CVL_IL10,CVL_TNFa,CVL_IFNg,CVL_IFNa, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("cytokines_CVL.png", width = 8.5, height = 11, units = "in")
cytokines_serum = plot_grid(serum_IL6,serum_IL10,serum_TNFa,serum_IFNg,serum_IFNa, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("cytokines_serum.png", width = 8.5, height = 11, units = "in")
# Others
#Eotaxin
#IL12
#IL1b
others_CVL = plot_grid(CVL_eotaxin,CVL_IL12,CVL_IL1b, nrow = 1, ncol = 3, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("others_CVL.png", width = 8.5, height = 11, units = "in")
others_serum = plot_grid(serum_eotaxin,serum_IL12,serum_IL1b, nrow = 1, ncol = 3, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("others_serum.png", width = 8.5, height = 11, units = "in")
# BAFF + Ig CVL/Sérum
# Ig insert cutoff de moyenne pop normale
#BAFF
#gp41
baff_ig = plot_grid(CVL_BAFF, CVL_gp41, serum_BAFF, serum_gp41, nrow = 2, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("baff_ig.png", width = 11, height = 11, units = "in")

## Layouts Johanne
layout_CVL = plot_grid(CVL_MIP1a,CVL_MIP1b,CVL_MCP1,CVL_MIG,CVL_IP10, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("layout_CVL.png", width = 10.5, height = 11, units = "in")
layout_serum = plot_grid(serum_MIP1a,serum_MIP1b,serum_MCP1,serum_MIG,serum_IP10, nrow = 3, ncol = 2, labels = "AUTO", align = "hv") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave("layout_serum.png", width = 10.5, height = 11, units = "in")
