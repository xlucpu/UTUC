library(TeachingDemos)
library(ChAMP)
library(ggplot2)
library(ggpubr)
library(plotly)
library(reshape2)
library(plyr)
library(clusterProfiler)
library(car)
library(rgl)
library(ClassDiscovery)
library(GSVA)
library(dendsort)
library(NMF)
require(grid)
library(RColorBrewer)
library(cowplot)
#set color
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
cherry <- "#700353"
lightgrey <- "#dcddde"

nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun <- "#E53435"
peach <- "#E43889"
violet <- "#89439B"
soil <- "#EC7D21"
lightgreen <- "#54B642"
darkblue <- "#21498D"
darkgreen <- "#009047"
brown <- "#874118"
seagreen <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")

workdir <- "H:/Jan2020/UTUC_2020/UTUC_BLCA/ChAMP"; setwd(workdir)
idat.path <- "H:/Jan2020/UTUC_2020/UTUC_BLCA/1493_Malouf EPIC DNA Methylation data package/IDAT Files"
############################
txtStart("F:/Project/UTUC_BLCA/ChAMP/CONSOLE_OUTPUT_NOFILTER_C1vsC2.txt")
#########################################################################
# c1 VS c2 (derive from nofilter clustering with 35 tumor samples only) #
#########################################################################
iDatDir.NOFILTER.C1vsC2 <- file.path(workdir,"IDAT_NOFILTER_C1vsC2")
resDir.NOFILTER.C1vsC2 <- file.path(iDatDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2")
champ.process(directory = iDatDir.NOFILTER.C1vsC2,
              resultsDir = resDir.NOFILTER.C1vsC2,
              arraytype = "EPIC",
              DMRmethod = "DMRcate",
              cores=1,
              runCNA = F,
              runBlock = F,
              runEpiMod = F)
###########################
txtStop()


#txtStart("F:/Project/UTUC_BLCA/ChAMP/CONSOLE_OUTPUT_FILTER_C1vsC2.txt")
#########################################################################
# c1 VS c2 (derive from filter clustering with 35 tumor samples only) #
#########################################################################
#iDatDir.FILTER.C1vsC2 <- file.path(workdir,"IDAT_FILTER_C1vsC2")
#resDir.FILTER.C1vsC2 <- file.path(iDatDir.FILTER.C1vsC2,"ChAMP_RES_FILTER_C1vsC2")
#champ.process(directory = iDatDir.FILTER.C1vsC2,
#              resultsDir = resDir.FILTER.C1vsC2,
#              arraytype = "EPIC",
#              DMRmethod = "DMRcate",
#              cores=1,
#              runCNA = F,
#              runBlock = F,
#              runEpiMod = F)
###########################
#txtStop()


txtStart("F:/Project/UTUC_BLCA/ChAMP/CONSOLE_OUTPUT_NOFILTER_TvsN.txt")
#########################################################################
# Tumor vs Normal (derive from nofilter clustering with 43 samples) #####
#########################################################################
iDatDir.NOFILTER.TvsN <- file.path(workdir,"IDAT_NOFILTER_TvsN")
resDir.NOFILTER.TvsN <- file.path(iDatDir.NOFILTER.TvsN,"ChAMP_RES_NOFILTER_TvsN")
champ.process(directory = iDatDir.NOFILTER.TvsN,
              resultsDir = resDir.NOFILTER.TvsN,
              arraytype = "EPIC",
              cores=1,
              DMRmethod = "DMRcate",
              runCNA = F,
              runBlock = F,
              runEpiMod = F)
############################
txtStop()
#Filtering probes with a detection p-value above 0.01.
#Removing 14365 probes.
#If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

#Filtering BeadCount Start
#Filtering probes with a beadcount <3 in at least 5% of samples.
#Removing 20495 probes

#Filtering NoCG Start
#Only Keep CpGs, removing 2766 probes from the analysis.

#Filtering SNPs Start
#Using general EPIC SNP list for filtering.
#Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#Removing 77022 probes from the analysis.

#Filtering MultiHit Start
#Filtering probes that align to multiple locations as identified in Nordlund et al
#Removing 49 probes from the analysis.

#Filtering XY Start
#Filtering probes located on X,Y chromosome, removing 16101 probes from the analysis.

#do GSEA by using DMP or DMR
myGSEA <- champ.GSEA(beta = myNorm,DMP = myDMP[[1]],DMR = myDMR,arraytype = "EPIC",adjPval = 0.05,method = "fisher")

#do ebayes GSEA by just using phonotype and universal beta value
myebayGSEA <- champ.ebayGSEA(beta = myNorm,pheno = myLoad$pd$Sample_Group,arraytype = "EPIC")

#do differential methylated interation hotspots -- network
myEpiMod <- champ.EpiMod(beta = myNorm,pheno = myLoad$pd$Sample_Group)

#do CNA : should determine controlgroup
myCNA <- champ.CNA(intensity = myLoad$intensity,pheno = myLoad$pd$Sample_Group,arraytype = "EPIC",controlGroup = "Normal")

###############################
# New analysis plan 8/16/2018 #
###############################
#transmit file for detailed Tumor vs Normal
#Tumor3 vs Normal1
Sinfo <- read.csv("F:/Project/UTUC_BLCA/Results/Sinfo of detailed Tumor Normal for RnBeads.csv",header = T,stringsAsFactors = F,check.names = F)
ref <- read.csv(file.path(iDatDir.NOFILTER.TvsN,"pd_43TumorNormalSamples_nofilter_Tumor vs Normal.csv"),header = T,stringsAsFactors = F,check.names = F)
samples <- Sinfo[which(Sinfo$Sample_Group %in% c("Tumor3","Normal1")),"Sample_Name"]
barcode <- ref[which(ref$Sample_Name %in% samples),"Complete Barcode"]
iDatDir.NOFILTER.T3vsN1 <- file.path(workdir,"IDAT_NOFILTER_T3vsN1")
resDir.NOFILTER.T3vsN1 <- file.path(iDatDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1")
#locate file
toDir <- iDatDir.NOFILTER.T3vsN1
file.path <- iDatDir.NOFILTER.TvsN
file.pattern <- paste0(barcode,".*.idat$")
for (pattern in file.pattern) {
  file <- dir(file.path,pattern = pattern)
  for (n in 1: length(file)){
    file.copy(file.path(file.path,file[n]), toDir, overwrite = T)
  }
}
if (!file.exists(resDir.NOFILTER.T3vsN1)) { dir.create(resDir.NOFILTER.T3vsN1) }
pd <- ref[which(ref$Sample_Name %in% samples),]
tmp <- Sinfo; rownames(tmp) <- tmp[,1]
pd$Sample_Group <- tmp[pd$Sample_Name,"Sample_Group"]
pd$Sentrix_ID <- factor(as.character(pd$Sentrix_ID))
pd <- as.data.frame(pd)
write.csv(pd,file.path(iDatDir.NOFILTER.T3vsN1,"pd_19TumorNormalSamples_nofilter_Tumor3 vs Normal1.csv"),row.names = F)

txtStart("F:/Project/UTUC_BLCA/ChAMP/CONSOLE_OUTPUT_NOFILTER_T3vsN1.txt")
#########################################################################
# Tumor3 vs Normal1 (derive from nofilter clustering with 43 samples) ###
#########################################################################
iDatDir.NOFILTER.T3vsN1 <- file.path(workdir,"IDAT_NOFILTER_T3vsN1")
resDir.NOFILTER.T3vsN1 <- file.path(iDatDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1")
#champ.process(directory = iDatDir.NOFILTER.T3vsN1,
#              resultsDir = resDir.NOFILTER.T3vsN1,
#              arraytype = "EPIC",
#              cores=1,
#              DMRmethod = "DMRcate",
#              runCNA = F,
#              runBlock = F,
#              runEpiMod = F)
myLoad <- champ.load(iDatDir.NOFILTER.T3vsN1,arraytype = "EPIC")
champ.QC(resultsDir = file.path(resDir.NOFILTER.T3vsN1,"CHAMP_QCimages"))
myNorm <- champ.norm(beta = myLoad$beta,arraytype = "EPIC",cores = 1,resultsDir = file.path(resDir.NOFILTER.T3vsN1,"CHAMP_Normalization")); save(myNorm,file = file.path(resDir.NOFILTER.T3vsN1,"myNorm.rda"))
myLoad$pd$Slide <- as.character(myLoad$pd$Slide); save(myLoad,file = file.path(resDir.NOFILTER.T3vsN1,"myLoad.rda"))
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide")); save(myCombat,file = file.path(resDir.NOFILTER.T3vsN1,"myCombat.rda"))
myDMP <- champ.DMP(beta=myCombat,pheno=myLoad$pd$Sample_Group,arraytype = "EPIC"); save(myDMP,file = file.path(resDir.NOFILTER.T3vsN1,"myDMP.rda"))
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate",arraytype = "EPIC",cores = 1); save(myDMR,file = file.path(resDir.NOFILTER.T3vsN1,"myDMR.rda"))
myGSEA <- champ.GSEA(beta=myCombat,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher"); save(myGSEA,file = file.path(resDir.NOFILTER.T3vsN1,"myGSEA.rda"))

############################
txtStop()

#Filtering probes with a detection p-value above 0.01.
#Removing 4099 probes.
#If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

#Filtering BeadCount Start
#Filtering probes with a beadcount <3 in at least 5% of samples.
#Removing 58717 probes

#Filtering NoCG Start
#Only Keep CpGs, removing 2672 probes from the analysis.

#Filtering SNPs Start
#Using general EPIC SNP list for filtering.
#Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#Removing 74861 probes from the analysis.

#Filtering MultiHit Start
#Filtering probes that align to multiple locations as identified in Nordlund et al
#Removing 47 probes from the analysis.

#Filtering XY Start
#Filtering probes located on X,Y chromosome, removing 15825 probes from the analysis.

#Updating PD file

#Fixing Outliers Start
#Replacing all value smaller/equal to 0 with smallest positive value.
#Replacing all value greater/equal to 1 with largest value below 1..
#[ Section 2: Filtering Done ]

#All filterings are Done, now you have 709697 probes and 19 samples.

#Tumor24 vs Normal1
Sinfo <- read.csv("F:/Project/UTUC_BLCA/Results/Sinfo of detailed Tumor Normal for RnBeads.csv",header = T,stringsAsFactors = F,check.names = F)
ref <- read.csv(file.path(iDatDir.NOFILTER.TvsN,"pd_43TumorNormalSamples_nofilter_Tumor vs Normal.csv"),header = T,stringsAsFactors = F,check.names = F)
samples <- Sinfo[which(Sinfo$Sample_Group %in% c("Tumor2","Tumor4","Normal1")),"Sample_Name"]
barcode <- ref[which(ref$Sample_Name %in% samples),"Complete Barcode"]
iDatDir.NOFILTER.T24vsN1 <- file.path(workdir,"IDAT_NOFILTER_T24vsN1")
resDir.NOFILTER.T24vsN1 <- file.path(iDatDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1")
#locate file
toDir <- iDatDir.NOFILTER.T24vsN1
file.path <- iDatDir.NOFILTER.TvsN
file.pattern <- paste0(barcode,".*.idat$")
for (pattern in file.pattern) {
  file <- dir(file.path,pattern = pattern)
  for (n in 1: length(file)){
    file.copy(file.path(file.path,file[n]), toDir, overwrite = T)
  }
}
if (!file.exists(resDir.NOFILTER.T24vsN1)) { dir.create(resDir.NOFILTER.T24vsN1) }
pd <- ref[which(ref$Sample_Name %in% samples),]
tmp <- Sinfo; rownames(tmp) <- tmp[,1]
pd$Sample_Group <- tmp[pd$Sample_Name,"Sample_Group"]
pd[which(pd$Sample_Group %in% c("Tumor2","Tumor4")),"Sample_Group"] <- "Tumor24"
pd$Sentrix_ID <- as.character(pd$Sentrix_ID)
write.csv(pd,file.path(iDatDir.NOFILTER.T24vsN1,"pd_19TumorNormalSamples_nofilter_Tumor24 vs Normal1.csv"),row.names = F)

txtStart("F:/Project/UTUC_BLCA/ChAMP/CONSOLE_OUTPUT_NOFILTER_T24vsN1.txt")
#########################################################################
# Tumor24vs Normal1 (derive from nofilter clustering with 43 samples) ###
#########################################################################
iDatDir.NOFILTER.T24vsN1 <- file.path(workdir,"IDAT_NOFILTER_T24vsN1")
resDir.NOFILTER.T24vsN1 <- file.path(iDatDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1")
myLoad <- champ.load(iDatDir.NOFILTER.T24vsN1,arraytype = "EPIC")
champ.QC(resultsDir = file.path(resDir.NOFILTER.T24vsN1,"CHAMP_QCimages"))
myNorm <- champ.norm(beta = myLoad$beta,arraytype = "EPIC",cores = 1,resultsDir = file.path(resDir.NOFILTER.T24vsN1,"CHAMP_Normalization")); save(myNorm,file = file.path(resDir.NOFILTER.T24vsN1,"myNorm.rda"))
myLoad$pd$Slide <- as.character(myLoad$pd$Slide); save(myLoad,file = file.path(resDir.NOFILTER.T24vsN1,"myLoad.rda"))
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide")); save(myCombat,file = file.path(resDir.NOFILTER.T24vsN1,"myCombat.rda"))
myDMP <- champ.DMP(beta=myCombat,pheno=myLoad$pd$Sample_Group,arraytype = "EPIC"); save(myDMP,file = file.path(resDir.NOFILTER.T24vsN1,"myDMP.rda"))
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate",arraytype = "EPIC",cores = 1); save(myDMR,file = file.path(resDir.NOFILTER.T24vsN1,"myDMR.rda"))
myGSEA <- champ.GSEA(beta=myCombat,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher"); save(myGSEA,file = file.path(resDir.NOFILTER.T24vsN1,"myGSEA.rda"))

############################
txtStop()

#Filtering probes with a detection p-value above 0.01.
#Removing 13489 probes.
#If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

#Filtering BeadCount Start
#Filtering probes with a beadcount <3 in at least 5% of samples.
#Removing 27585 probes

#Filtering NoCG Start
#Only Keep CpGs, removing 2786 probes from the analysis.

#Filtering SNPs Start
#Using general EPIC SNP list for filtering.
#Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#    Removing 76604 probes from the analysis.

#  Filtering MultiHit Start
#    Filtering probes that align to multiple locations as identified in Nordlund et al
#    Removing 49 probes from the analysis.

#  Filtering XY Start
#    Filtering probes located on X,Y chromosome, removing 15920 probes from the analysis.

#  Updating PD file

#  Fixing Outliers Start
#    Replacing all value smaller/equal to 0 with smallest positive value.
#    Replacing all value greater/equal to 1 with largest value below 1..
#[ Section 2: Filtering Done ]

# All filterings are Done, now you have 729485 probes and 32 samples.

######################################################################
# Tumor vs Normal (derive from nofilter clustering with 43 samples) ##
######################################################################
iDatDir.NOFILTER.TvsN <- file.path(workdir,"IDAT_NOFILTER_TvsN")
resDir.NOFILTER.TvsN <- file.path(iDatDir.NOFILTER.TvsN,"ChAMP_RES_NOFILTER_TvsN")
myLoad <- champ.load(iDatDir.NOFILTER.TvsN,arraytype = "EPIC")
champ.QC(resultsDir = file.path(resDir.NOFILTER.TvsN,"CHAMP_QCimages"))
myNorm <- champ.norm(beta = myLoad$beta,arraytype = "EPIC",cores = 1,resultsDir = file.path(resDir.NOFILTER.TvsN,"CHAMP_Normalization")); save(myNorm,file = file.path(resDir.NOFILTER.TvsN,"myNorm.rda"))
myLoad$pd$Slide <- as.character(myLoad$pd$Slide); save(myLoad,file = file.path(resDir.NOFILTER.TvsN,"myLoad.rda"))
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide")); save(myCombat,file = file.path(resDir.NOFILTER.TvsN,"myCombat.rda"))
myDMP <- champ.DMP(beta=myCombat,pheno=myLoad$pd$Sample_Group,arraytype = "EPIC"); save(myDMP,file = file.path(resDir.NOFILTER.TvsN,"myDMP.rda"))
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate",arraytype = "EPIC",cores = 1); save(myDMR,file = file.path(resDir.NOFILTER.TvsN,"myDMR.rda"))
myGSEA <- champ.GSEA(beta=myCombat,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher"); save(myGSEA,file = file.path(resDir.NOFILTER.TvsN,"myGSEA.rda"))

# Filtering probes with a detection p-value above 0.01.
# Removing 14365 probes.
# If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
# 
# Filtering BeadCount Start
# Filtering probes with a beadcount <3 in at least 5% of samples.
# Removing 20495 probes
# 
# Filtering NoCG Start
# Only Keep CpGs, removing 2766 probes from the analysis.
# 
# Filtering SNPs Start
# Using general EPIC SNP list for filtering.
# Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
# Removing 77022 probes from the analysis.
# 
# Filtering MultiHit Start
# Filtering probes that align to multiple locations as identified in Nordlund et al
# Removing 49 probes from the analysis.
# 
# Filtering XY Start
# Filtering probes located on X,Y chromosome, removing 16101 probes from the analysis.
# 
# Updating PD file
# 
# Fixing Outliers Start
# Replacing all value smaller/equal to 0 with smallest positive value.
# Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]
# 
# All filterings are Done, now you have 735120 probes and 43 samples.


######################################################################
# C12 vs Normal (derive from nofilter clustering with 43 samples) ####
######################################################################
iDatDir.NOFILTER.C12vsN <- file.path(workdir,"IDAT_NOFILTER_C12vsN")
resDir.NOFILTER.C12vsN <- file.path(iDatDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER_C12vsN")
myLoad <- champ.load(iDatDir.NOFILTER.C12vsN,arraytype = "EPIC")
champ.QC(resultsDir = file.path(resDir.NOFILTER.C12vsN,"CHAMP_QCimages"))
myNorm <- champ.norm(beta = myLoad$beta,arraytype = "EPIC",cores = 1,resultsDir = file.path(resDir.NOFILTER.C12vsN,"CHAMP_Normalization")); save(myNorm,file = file.path(resDir.NOFILTER.C12vsN,"myNorm.rda"))
myLoad$pd$Slide <- as.character(myLoad$pd$Slide); save(myLoad,file = file.path(resDir.NOFILTER.C12vsN,"myLoad.rda"))
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide")); save(myCombat,file = file.path(resDir.NOFILTER.C12vsN,"myCombat.rda"))
myDMP <- champ.DMP(beta=myCombat,pheno=myLoad$pd$Sample_Group,arraytype = "EPIC",adjPVal = 1); save(myDMP,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.rda"))

myCombat.C12 <- myCombat[,1:35]
save(myCombat.C12,file = file.path(resDir.NOFILTER.C12vsN,"myCombat.C12.rda"))

myCombat.C1N <- myCombat[,c(1:3,6:11,14:15,19:21,23:25,27,29:30,33:35,36:43)]
save(myCombat.C1N,file = file.path(resDir.NOFILTER.C12vsN,"myCombat.C1N.rda"))

myCombat.C2N <- myCombat[,c(4:5,12:13,16:18,22,26,28,31:32,36:43)]
save(myCombat.C2N,file = file.path(resDir.NOFILTER.C12vsN,"myCombat.C2N.rda"))

pheno.C12 <- myLoad$pd$Sample_Group[1:35]
pheno.C1N <- myLoad$pd$Sample_Group[c(1:3,6:11,14:15,19:21,23:25,27,29:30,33:35,36:43)]
pheno.C2N <- myLoad$pd$Sample_Group[c(4:5,12:13,16:18,22,26,28,31:32,36:43)]
myDMP.C1vsC2 <- champ.DMP(beta=myCombat.C12,pheno=pheno.C12,arraytype = "EPIC",adjPVal = 0.05); save(myDMP.C1vsC2,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.C1vsC2.rda"))
myDMP.C1vsN <- champ.DMP(beta=myCombat.C1N,pheno=pheno.C1N,arraytype = "EPIC",adjPVal = 1); save(myDMP.C1vsN,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.C1vsN.rda"))
myDMP.C2vsN <- champ.DMP(beta=myCombat.C2N,pheno=pheno.C2N,arraytype = "EPIC",adjPVal = 1); save(myDMP.C2vsN,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.C2vsN.rda"))

myDMR.C1vsC2 <- champ.DMR(beta=myCombat.C12,pheno=pheno.C12,method="DMRcate",arraytype = "EPIC",cores = 1); save(myDMR.C1vsC2,file = file.path(resDir.NOFILTER.C12vsN,"myDMR.C1vsC2.rda"))
myGSEA.C1vsC2 <- champ.GSEA(beta=myCombat.C12,DMP=myDMP.C1vsC2[[1]], DMR=myDMR.C1vsC2,pheno = pheno.C12 , arraytype="EPIC",adjPval=0.05, method="fisher"); save(myGSEA.C1vsC2,file = file.path(resDir.NOFILTER.C12vsN,"myGSEA.C1vsC2.rda"))

# Filtering probes with a detection p-value above 0.01.
# Removing 14365 probes.
# If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
# 
# Filtering BeadCount Start
# Filtering probes with a beadcount <3 in at least 5% of samples.
# Removing 20495 probes
# 
# Filtering NoCG Start
# Only Keep CpGs, removing 2766 probes from the analysis.
# 
# Filtering SNPs Start
# Using general EPIC SNP list for filtering.
# Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#     Removing 77022 probes from the analysis.
# 
#   Filtering MultiHit Start
#     Filtering probes that align to multiple locations as identified in Nordlund et al
#     Removing 49 probes from the analysis.
# 
#   Filtering XY Start
#     Filtering probes located on X,Y chromosome, removing 16101 probes from the analysis.
# 
#   Updating PD file
# 
#   Fixing Outliers Start
#     Replacing all value smaller/equal to 0 with smallest positive value.
#     Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]
# 
#  All filterings are Done, now you have 735120 probes and 43 samples.

# create DMP information about C1 vs C2 + C1 vs Normal + C2 vs Normal
tmp <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInC1 <- rownames(tmp[which(tmp$C1_AVG >= 0.4 & tmp$C2_AVG <= 0.2 & tmp$adj.P.Val < 0.05),])
tmp <- tmp[gainMethInC1,]; colnames(tmp)[1:9] <- paste0("C1vsC2_",colnames(tmp)[1:9])

tmp1 <- read.table(file.path(resDir.NOFILTER.C12vsN,"CHAMP_RES_NOFILTER.C12N_C1vsN_DMP.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
tmp1 <- tmp1[gainMethInC1,]; colnames(tmp1)[1:9] <- paste0("C1vsN_",colnames(tmp1)[1:9])

tmp2 <- read.table(file.path(resDir.NOFILTER.C12vsN,"CHAMP_RES_NOFILTER.C12N_C2vsN_DMP.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
tmp2 <- tmp2[gainMethInC1,]; colnames(tmp2)[1:9] <- paste0("C2vsN_",colnames(tmp2)[1:9])

rs <- cbind.data.frame(tmp,tmp1[,c(1,2,4,5,7,8,9)],tmp2[,c(1,2,4,5,7,8,9)])
write.table(rs,file.path(resDir.NOFILTER.C12vsN,"CHAMP_RES_NOFILTER.C12N_C1vsC2_DMP_add_C1vsN_C2vsN_addFeatures.txt"),sep = "\t",row.names = T,col.names = NA)

tmp3 <- read.table("F:/Project/UTUC_BLCA/Results/mRNA_rmBatch_edgeR_test_result.C1_vs_C2.txt",sep = "\t",check.names = F,stringsAsFactors = F,row.names = NULL,header = T)
tmp4 <- read.table(file.path(resDir.NOFILTER.C12vsN,"CHAMP_RES_NOFILTER.C12N_C1vsC2_DMP_add_C1vsN_C2vsN_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = NULL,header = T)
tmp4$gene <- toupper(tmp4$gene); colnames(tmp4)[1] <- "probe"
tmp3$id <- toupper(tmp3$id); colnames(tmp3) <- c("gene","ExpressionlogFC","ExpressionlogCPM","ExpressionLR","ExpressionPvalue","ExpressionFDR") 
tmp5 <- merge(tmp4,tmp3,by.x = "gene",all.x=T)
tmp5$gene <- gsub(tmp5$gene,pattern = "^$",replacement = "NA")
write.table(tmp5,file.path(resDir.NOFILTER.C12vsN,"CHAMP_RES_NOFILTER.C12N_C1vsC2_DMP_add_C1vsN_C2vsN_DE.txt"),sep = "\t",row.names = F)
###############################
## Analysis plan : 8/16/2018 ##
###############################

# loose criterion with just consideration of diff
# 1. for NOfiltered clustering C1 vs C2
# 1.1 Calculate which genomic regions are enriched for gain of methylation in C1 versus C2 (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInC1 <- rownames(dmp.res[which(dmp.res$deltaBeta <= (-0.2) & dmp.res$adj.P.Val < 0.05),])

tmp <- dmp.res[gainMethInC1,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInC1)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInC1),")"))
ggsave(file.path(resDir.NOFILTER.C1vsC2,"Stacked bar of enrichment of genomic regions for gainMeth in nofiltered C1 vs C2.pdf"))

pbForSpCluster <- rownames(dmp.res[which(abs(dmp.res$deltaBeta) >= 0.2 & dmp.res$adj.P.Val < 0.05),])
pbForSpCluster_diff0.3_fdr0.01 <- rownames(dmp.res[which(abs(dmp.res$deltaBeta) >= 0.3 & dmp.res$adj.P.Val < 0.01),])

save(pbForSpCluster_diff0.3_fdr0.01,file = file.path(resDir.NOFILTER.C1vsC2,"pbForSpCluster_diff0.3_fdr0.01.rda"))
#ggplot(tmp, aes(GenomicRegions, Counts)) +   
#  geom_bar(aes(fill = LogicResults), position = "dodge", stat="identity") + 
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 1.2 Calculate which genomic regions are enriched for gain of methylation in C2 versus C1 (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
gainMethInC2 <- rownames(dmp.res[which(dmp.res$deltaBeta >= 0.2 & dmp.res$adj.P.Val < 0.05),])

tmp <- dmp.res[gainMethInC2,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInC2)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInC2),")"))
ggsave(file.path(resDir.NOFILTER.C1vsC2,"Stacked bar of enrichment of genomic regions for gainMeth in nofiltered C2 vs C1.pdf"))

# 2. for clustering Tumor vs Normal
# 2.1 Calculate which genomic regions are enriched for gain of methylation in Tumor3 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInT3 <- rownames(dmp.res[which(dmp.res$deltaBeta >= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInT3,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInT3)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInT3),")"))
ggsave(file.path(resDir.NOFILTER.T3vsN1,"Stacked bar of enrichment of genomic regions for gainMeth in nofiltered T3 vs N1.pdf"))

# 2.2 Calculate which genomic regions are enriched for gain of methylation in Tumor24 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInT24 <- rownames(dmp.res[which(dmp.res$deltaBeta >= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInT24,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInT24)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInT24),")"))
ggsave(file.path(resDir.NOFILTER.T24vsN1,"Stacked bar of enrichment of genomic regions for gainMeth in nofiltered T24 vs N1.pdf"))

# 2.3 Calculate which genomic regions are enriched for loss of methylation in Tumor3 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
lossMethInT3 <- rownames(dmp.res[which(dmp.res$deltaBeta <= (-0.2) & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[lossMethInT3,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(lossMethInT3)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(lossMethInT3),")"))
ggsave(file.path(resDir.NOFILTER.T3vsN1,"Stacked bar of enrichment of genomic regions for lossMeth in nofiltered T3 vs N1.pdf"))

# 2.4 Calculate which genomic regions are enriched for loss of methylation in Tumor24 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
lossMethInT24 <- rownames(dmp.res[which(dmp.res$deltaBeta <= (-0.2) & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[lossMethInT24,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(lossMethInT24)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(lossMethInT24),")"))
ggsave(file.path(resDir.NOFILTER.T24vsN1,"Stacked bar of enrichment of genomic regions for lossMeth in nofiltered T24 vs N1.pdf"))

# strict criterion with one methylation <= 0.2 and another one >= 0.4 which makes diff >= 0.2
# 3. for NOfiltered clustering C1 vs C2
# 3.1 Calculate which genomic regions are enriched for strict gain of methylation in C1 versus C2 (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInC1,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInC1)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInC1),")"))
ggsave(file.path(resDir.NOFILTER.C1vsC2,"Stacked bar of enrichment of genomic regions for strict gainMeth in nofiltered C1 vs C2.pdf"))

# 3.2 Calculate which genomic regions are enriched for strict gain of methylation in C2 versus C1 (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
gainMethInC2 <- rownames(dmp.res[which(dmp.res$C1_AVG <= 0.2 & dmp.res$C2_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])

tmp <- dmp.res[gainMethInC2,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInC2)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInC2),")"))
ggsave(file.path(resDir.NOFILTER.C1vsC2,"Stacked bar of enrichment of genomic regions for strict gainMeth in nofiltered C2 vs C1.pdf"))

# 4. for clustering Tumor vs Normal
# 4.1 Calculate which genomic regions are enriched for strict gain of methylation in Tumor3 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInT3 <- rownames(dmp.res[which(dmp.res$Tumor3_AVG >= 0.4 & dmp.res$Normal1_AVG <=0.2  & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInT3,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInT3)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInT3),")"))
ggsave(file.path(resDir.NOFILTER.T3vsN1,"Stacked bar of enrichment of genomic regions for strict gainMeth in nofiltered T3 vs N1.pdf"))

# 4.2 Calculate which genomic regions are enriched for strict gain of methylation in Tumor24 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInT24 <- rownames(dmp.res[which(dmp.res$Tumor24_AVG >= 0.4 & dmp.res$Normal1_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInT24,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(gainMethInT24)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(gainMethInT24),")"))
ggsave(file.path(resDir.NOFILTER.T24vsN1,"Stacked bar of enrichment of genomic regions for strict gainMeth in nofiltered T24 vs N1.pdf"))

# 4.3 Calculate which genomic regions are enriched for strict loss of methylation in Tumor3 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
lossMethInT3 <- rownames(dmp.res[which(dmp.res$Tumor3_AVG <= 0.2 & dmp.res$Normal1_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[lossMethInT3,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(lossMethInT3)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(lossMethInT3),")"))
ggsave(file.path(resDir.NOFILTER.T3vsN1,"Stacked bar of enrichment of genomic regions for strict lossMeth in nofiltered T3 vs N1.pdf"))

# 4.4 Calculate which genomic regions are enriched for strict loss of methylation in Tumor24 versus Normal (CGI, shore, shelf, enhancers, promoter, body, open sea, etc¡­)
dmp.res <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
lossMethInT24 <- rownames(dmp.res[which(dmp.res$Tumor24_AVG <= 0.2 & dmp.res$Normal1_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[lossMethInT24,21:30]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:30])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp <- melt(tmp,id.vars="Names"); tmp$value <- as.numeric(tmp$value); colnames(tmp) <- c("GenomicRegions","LogicResults","Counts")
tmp$Percentage <- round(tmp$Counts/length(lossMethInT24)*100,0)
tmp <- tmp[order(tmp$LogicResults,decreasing = T),]
tmp <- ddply(tmp, .(GenomicRegions),
             transform, pos = cumsum(Percentage) - (0.5 * Percentage))
ggplot() + 
  geom_bar(aes(y = Percentage, x = GenomicRegions, fill = LogicResults), data = tmp,stat="identity") + scale_fill_manual(values = c(red,blue)) +
  geom_text(data=tmp, aes(x = GenomicRegions, y = pos, label = paste0(Percentage,"%")),size=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank()) + labs(y = paste0("Percentage (N = ",length(lossMethInT24),")"))
ggsave(file.path(resDir.NOFILTER.T24vsN1,"Stacked bar of enrichment of genomic regions for strict lossMeth in nofiltered T24 vs N1.pdf"))

gainMethInC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
gainMethInC2 <- rownames(dmp.res[which(dmp.res$C1_AVG <= 0.2 & dmp.res$C2_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])
pbForSpCluster_strict_diff0.2_fdr0.05 <- c(gainMethInC1,gainMethInC2)
save(pbForSpCluster_strict_diff0.2_fdr0.05,file = file.path(resDir.NOFILTER.C1vsC2,"pbForSpCluster_strict_diff0.2_fdr0.05.rda"))


# 5 perform GSEA
# 5.1 C1 vs C2
dmp.res <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
dmp.gainC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
dmp.Enhancer <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isEnhancer == TRUE),]))
dmp.PromCGI <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isPromCGI == TRUE),]))
dmp.TFBS <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isTFBS == TRUE),]))

#save(dmp.Enhancer,file = file.path(resDir.NOFILTER.C1vsC2,"pbEnhancerForSpClusterTCGA.rda"))

load(file = file.path(resDir.NOFILTER.C1vsC2,"myNorm.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMP.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMR.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myLoad.rda"))
CpGlist.C1vsC2 <- list()
CpGlist.C1vsC2[["gainMethC1"]] <- dmp.gainC1
CpGlist.C1vsC2[["Enhancer"]] <- dmp.Enhancer
CpGlist.C1vsC2[["PromCGI"]] <- dmp.PromCGI
CpGlist.C1vsC2[["TFBS"]] <- dmp.TFBS

myGSEA <- champ.GSEA(CpGlist =  CpGlist.C1vsC2,arraytype = "EPIC")
for (i in 1:length(myGSEA)) {
  label <- names(myGSEA)[i]
  if(label %in% c("DMP","DMR")) {
    tmp <- myGSEA[[i]]
    write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C1vsC2_entire_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
    next()}
  tmp <- myGSEA[[i]]
  write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C1vsC2_strict_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}

# 5.2 C2 vs C1
dmp.gainC2 <- rownames(dmp.res[which(dmp.res$C1_AVG <= 0.2 & dmp.res$C2_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])
dmp.Enhancer <- intersect(dmp.gainC2,rownames(dmp.res[which(dmp.res$isEnhancer == TRUE),]))
dmp.PromCGI <- intersect(dmp.gainC2,rownames(dmp.res[which(dmp.res$deltaBeta >= 0.2 & dmp.res$adj.P.Val < 0.05 & dmp.res$isPromCGI == TRUE),]))
dmp.TFBS <- intersect(dmp.gainC2,rownames(dmp.res[which(dmp.res$deltaBeta >= 0.2 & dmp.res$adj.P.Val < 0.05 & dmp.res$isTFBS == TRUE),]))

CpGlist.C2vsC1 <- list()
CpGlist.C2vsC1[["gainMethC2"]] <- dmp.gainC2
CpGlist.C2vsC1[["Enhancer"]] <- dmp.Enhancer
CpGlist.C2vsC1[["PromCGI"]] <- dmp.PromCGI
CpGlist.C2vsC1[["TBFS"]] <- dmp.TBFS
myGSEA <- champ.GSEA(CpGlist =  CpGlist.C2vsC1,arraytype = "EPIC")
for (i in 1:length(myGSEA)) {
  label <- names(myGSEA)[i]
  if(label %in% c("DMP","DMR")) {next()}
  tmp <- myGSEA[[i]]
  write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C2vsC1_strict_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}

load(file = file.path(resDir.NOFILTER.C1vsC2,"myNorm.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMP.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMR.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myLoad.rda"))
CpGlist.C1vsC2 <- list()
CpGlist.C1vsC2[["gainMethC1"]] <- dmp.gainC1
CpGlist.C1vsC2[["Enhancer"]] <- dmp.Enhancer
CpGlist.C1vsC2[["PromCGI"]] <- dmp.PromCGI
CpGlist.C1vsC2[["TFBS"]] <- dmp.TFBS

myGSEA <- champ.GSEA(CpGlist =  CpGlist.C1vsC2,arraytype = "EPIC")
for (i in 1:length(myGSEA)) {
  label <- names(myGSEA)[i]
  if(label %in% c("DMP","DMR")) {
    tmp <- myGSEA[[i]]
    write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C1vsC2_entire_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
    next()}
  tmp <- myGSEA[[i]]
  write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C1vsC2_strict_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}
###

# 5.3 C1 vs C2 in C12N dataset
dmp.res <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
dmp.gainC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
dmp.Enhancer <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isEnhancer == TRUE),]))
dmp.PromCGI <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isPromCGI == TRUE),]))
dmp.Shelf <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isShelf == TRUE),]))
dmp.Shore <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isShore == TRUE),]))
dmp.UTRCGI <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isUTRCGI == TRUE),]))
dmp.BodyCGI <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isBodyCGI == TRUE),]))
dmp.OpenChr <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isOpenChromatin == TRUE),]))
dmp.DMR <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isDMR == TRUE),]))
dmp.DNase <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isDNase == TRUE),]))
dmp.CGI <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isCGI == TRUE),]))
dmp.TFBS <- intersect(dmp.gainC1,rownames(dmp.res[which(dmp.res$isTFBS == TRUE),]))

#pbForSpClusterTCGA <- dmp.gainC1
#save(pbForSpClusterTCGA,file = file.path(resDir.NOFILTER.C12vsN,"pbForSpClusterTCGA_14209.rda"))

load(file = file.path(resDir.NOFILTER.C12vsN,"myDMP.C1vsC2.rda"))
load(file = file.path(resDir.NOFILTER.C12vsN,"myDMR.C1vsC2.rda"))
load(file = file.path(resDir.NOFILTER.C12vsN,"myCombat.C12.rda"))

CpGlist.C1vsC2.all <- list()
CpGlist.C1vsC2.all[["gainMethC1"]] <- dmp.gainC1
CpGlist.C1vsC2.all[["Enhancer"]] <- dmp.Enhancer
CpGlist.C1vsC2.all[["PromCGI"]] <- dmp.PromCGI
CpGlist.C1vsC2.all[["Shelf"]] <- dmp.Shelf
CpGlist.C1vsC2.all[["Shore"]] <- dmp.Shore
CpGlist.C1vsC2.all[["TFBS"]] <- dmp.TFBS
CpGlist.C1vsC2.all[["UTRCGI"]] <- dmp.UTRCGI
CpGlist.C1vsC2.all[["BodyCGI"]] <- dmp.BodyCGI
CpGlist.C1vsC2.all[["OpenChr"]] <- dmp.OpenChr
CpGlist.C1vsC2.all[["annoDMR"]] <- dmp.DMR # DMR annotated in feature annotation, not calculated DMR
CpGlist.C1vsC2.all[["DNase"]] <- dmp.DNase
CpGlist.C1vsC2.all[["CGI"]] <- dmp.CGI

myGSEA <- champ.GSEA(beta=myCombat.C12,DMP=myDMP.C1vsC2[[1]],DMR=myDMR.C1vsC2,CpGlist = CpGlist.C1vsC2.all,pheno = pheno.C12,arraytype="EPIC",adjPval=0.05, method="fisher")

for (i in 1:length(myGSEA)) {
  label <- names(myGSEA)[i]
  if(label %in% c("DMP","DMR")) {
    tmp <- myGSEA[[i]]
    write.table(tmp,file.path(resDir.NOFILTER.C12vsN,paste0("ChAMP_RES_NOFILTER.C12N_C1vsC2_entire_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
    next()}
  tmp <- myGSEA[[i]]
  write.table(tmp,file.path(resDir.NOFILTER.C12vsN,paste0("ChAMP_RES_NOFILTER.C12N_C1vsC2_strict_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}



##################################
# Analysis plan dated 09/03/2018 #
##################################

# CNA for Tumor vs Normal
load(file = file.path(resDir.NOFILTER.TvsN,"myNorm.rda"))
load(file = file.path(resDir.NOFILTER.TvsN,"myDMP.rda"))
load(file = file.path(resDir.NOFILTER.TvsN,"myDMR.rda"))
load(file = file.path(resDir.NOFILTER.TvsN,"myLoad.rda"))
colnames(myLoad$intensity) <- paste0("UTUC_",colnames(myLoad$intensity))
myCNA <- champ.CNA(intensity = myLoad$intensity,pheno = myLoad$pd$Sample_Group,arraytype = "EPIC",controlGroup = "Normal",resultsDir=file.path(resDir.NOFILTER.TvsN,"CHAMP_CNA"))
write.table(myCNA[[2]]$Tumor,file = file.path(resDir.NOFILTER.TvsN,"CNA groupResult.txt"),sep = "\t",row.names = F)
feature <- read.table("F:/Project/UTUC_BLCA/Annotation/Methylation_probes_features_logical2.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
enhancer.pb <- rownames(feature[which(feature$isEnhancer==T),])

# ) the first group (G1) are related to DMR and DMP differently methylated between C1 and C2 from one side and normal samples from one side
# if use TvsN as reference, there is no results, thus we use T24vsN1 and T3vsN1
# 1. read diff res data of C1vsC2 and TvsN
tmp1 <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp2 <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp3 <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp4 <- read.table(file.path(resDir.NOFILTER.TvsN,"ChAMP_RES_NOFILTER_TvsN_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# 2. hyper in C1 and hyper in Tumor24 or Tumor3
tmp5 <- tmp1[which(tmp1$C1_AVG >= 0.4 & tmp1$C2_AVG <= 0.2 & tmp1$adj.P.Val < 0.05),]
tmp6 <- tmp2[which(tmp2$Tumor24_AVG >= 0.4 & tmp2$Normal1_AVG <= 0.2 & tmp2$adj.P.Val < 0.05),]
tmp7 <- tmp3[which(tmp3$Tumor3_AVG >= 0.4 & tmp3$Normal1_AVG <= 0.2 & tmp3$adj.P.Val < 0.05),]
tmp8 <- tmp4[which(tmp4$Tumor_AVG >= 0.4 & tmp4$Normal_AVG <= 0.2 & tmp4$adj.P.Val < 0.05),]

hyperC1_hyperT24 <- intersect(rownames(tmp5),rownames(tmp6)) #5
# "cg23187316" "cg11025705" "cg17124224" "cg19516515" "cg23947326"
hyperC1_hyperT3 <- intersect(rownames(tmp5),rownames(tmp7)) #81
hyperC1_hyperT <- intersect(rownames(tmp5),rownames(tmp8)) #7

# 3. hyper in C1 and hyper in Normal1.24(refer to T24vsN1) or Normal1.3(refer to T3vsN1)
tmp5 <- tmp1[which(tmp1$C1_AVG >= 0.4 & tmp1$C2_AVG <= 0.2 & tmp1$adj.P.Val < 0.05),]
tmp6 <- tmp2[which(tmp2$Tumor24_AVG <= 0.2 & tmp2$Normal1_AVG >= 0.4 & tmp2$adj.P.Val < 0.05),]
tmp7 <- tmp3[which(tmp3$Tumor3_AVG <= 0.2 & tmp3$Normal1_AVG >= 0.4 & tmp3$adj.P.Val < 0.05),]
tmp8 <- tmp4[which(tmp4$Tumor_AVG <= 0.2 & tmp4$Normal_AVG >= 0.4 & tmp4$adj.P.Val < 0.05),]

hyperC1_hyperN1.24 <- intersect(rownames(tmp5),rownames(tmp6)) #2713
hyperC1_hyperN1.3 <- intersect(rownames(tmp5),rownames(tmp7)) #0
hyperC1_hyperN <- intersect(rownames(tmp5),rownames(tmp8)) #0

hyperC1_hyperN1.24.enhancer <- intersect(hyperC1_hyperN1.24,enhancer.pb) #551
hyperC1_hyperT24.enhancer <- intersect(hyperC1_hyperT24,enhancer.pb) #2
hyperC1_hyperT3.enhancer <- intersect(hyperC1_hyperT3,enhancer.pb) #27
hyperC1_hyperT.enhancer <- intersect(hyperC1_hyperT,enhancer.pb) #4

#the second group (G2) are related to DMR and DMP differently methylated between C1 and C2 from one side but no difference is observed between those and  normal samples.
# 1. read diff res data of C1vsC2 and TvsN
tmp1 <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp2 <- read.table(file.path(resDir.NOFILTER.T24vsN1,"ChAMP_RES_NOFILTER_T24vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp3 <- read.table(file.path(resDir.NOFILTER.T3vsN1,"ChAMP_RES_NOFILTER_T3vsN1_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp4 <- read.table(file.path(resDir.NOFILTER.TvsN,"ChAMP_RES_NOFILTER_TvsN_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

# 2. get exclude probes that are not differentially methylated
tmp5 <- tmp1[which(tmp1$C1_AVG >= 0.4 & tmp1$C2_AVG <= 0.2 & tmp1$adj.P.Val < 0.05),]
tmp6 <- setdiff(rownames(feature),rownames(tmp2)) #non-DMPs between T24vsN1
tmp7 <- setdiff(rownames(feature),rownames(tmp3)) #non-DMPs between T3vsN1
tmp8 <- setdiff(rownames(feature),rownames(tmp4)) #non-DMPs between TvsN

hyperC1_nodiff.T24vsN1 <- intersect(rownames(tmp5),tmp6) #1024
hyperC1_nodiff.T3vsN1 <- intersect(rownames(tmp5),tmp7) #3481
hyperC1_nodiff.TvsN <- intersect(rownames(tmp5),tmp8) #2494

hyperC1_nodiff.T24vsN1.enhancer <- intersect(hyperC1_nodiff.T24vsN1,enhancer.pb) #294
hyperC1_nodiff.T3vsN1.enhancer <- intersect(hyperC1_nodiff.T3vsN1,enhancer.pb) #871
hyperC1_nodiff.TvsN.enhancer <- intersect(hyperC1_nodiff.TvsN,enhancer.pb) #687

#perfrom GSEA to these enhancer probes (load C1vsC2 results)
load(file = file.path(resDir.NOFILTER.C1vsC2,"myNorm.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMP.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myDMR.rda"))
load(file = file.path(resDir.NOFILTER.C1vsC2,"myLoad.rda"))
CpGlist.Glist <- list()
CpGlist.Glist[["hyperC1_hyperT3"]] <- hyperC1_hyperT3
CpGlist.Glist[["hyperC1_hyperT3.enhancer"]] <- hyperC1_hyperT3.enhancer
CpGlist.Glist[["hyperC1_hyperN1.24"]] <- hyperC1_hyperN1.24
CpGlist.Glist[["hyperC1_nodiff.T24vsN1"]] <- hyperC1_nodiff.T24vsN1
CpGlist.Glist[["hyperC1_nodiff.T24vsN1.enhancer"]] <- hyperC1_nodiff.T24vsN1.enhancer
CpGlist.Glist[["hyperC1_nodiff.T3vsN1"]] <- hyperC1_nodiff.T3vsN1
CpGlist.Glist[["hyperC1_nodiff.T3vsN1.enhancer"]] <- hyperC1_nodiff.T3vsN1.enhancer
CpGlist.Glist[["hyperC1_nodiff.TvsN"]] <- hyperC1_nodiff.TvsN
CpGlist.Glist[["hyperC1_nodiff.TvsN.enhancer"]] <- hyperC1_nodiff.TvsN.enhancer
save(CpGlist.Glist,file=file.path(resDir.NOFILTER.C1vsC2,"CpGlist.Glist.rda"))
myGSEA2 <- champ.GSEA(CpGlist =  CpGlist.Glist,arraytype = "EPIC")
for (i in 1:length(myGSEA2)) {
  label <- names(myGSEA2)[i]
  if(label %in% c("DMP","DMR")) {
     next()}
  tmp <- myGSEA2[[i]]
  write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,paste0("ChAMP_RES_NOFILTER_C1vsC2_",label,"_GSEA.txt")),sep = "\t",row.names = T,col.names = NA)
}

dmp.res <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
gainMethInC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
write.table(dmp.res[gainMethInC1,],file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeature_strict.txt"),sep = "\t",row.names = T,col.names = NA)

tmp <- as.data.frame(myNorm)
tmp <- tmp[gainMethInC1,]
tmp <- cbind.data.frame(tmp,dmp.res[gainMethInC1,21:30])
write.table(tmp,file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeature_strict_normBeta.txt"),sep = "\t",row.names = T,col.names = NA)

#############################################################
# compare percentage of each region to the whole EPIC array #
#dmp.res <- read.table(file.path(resDir.NOFILTER.C1vsC2,"ChAMP_RES_NOFILTER_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
dmp.res <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)

gainMethInC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInC1,21:31]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:31])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/14209*100,0)

tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage")
tmp3$Class <- rep(c("C1 vs C2","EPIC"),each=11)
ggplot(tmp3,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_continuous(breaks=seq(0, 60, 5), limits=c(0, 60))
ggsave(file.path(resDir.NOFILTER.C12vsN,"pair-wise histogram for C1vsC2 and EPIC array.pdf"))

gainMethInC2 <- rownames(dmp.res[which(dmp.res$C1_AVG <= 0.2 & dmp.res$C2_AVG >= 0.4 & dmp.res$adj.P.Val < 0.05),])

#############################################################
# Histogram for all the CGI and outside CGI
gainMethInC1 <- rownames(dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),])
tmp <- dmp.res[gainMethInC1,21:31]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:31])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/14209*100,0)
tmp$CGI <- ifelse(grepl("CGI",tmp$Names),"InCGI","OutsideCGI")
tmp1 <- tmp %>% group_by(CGI) %>% summarise(count = sum(percentage))

tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)
tmp2$CGI <- ifelse(grepl("CGI",tmp2$Names),"InCGI","OutsideCGI")

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage","CGI")
tmp3$Class <- rep(c("C1 vs C2","EPIC"),each=11)

p1 <- ggplot(tmp3, aes(x = as.numeric(interaction(Class,CGI)), y = Percentage, fill = Region, width=0.8)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=c(nake,gold,cyan,sun,peach,violet,soil,lightgreen,darkblue,darkgreen,brown)) +
  scale_x_continuous(breaks=c(1,2,3,4),labels=c("CGI: C1 vs C2","CGI: EPIC","non-CGI: c1 vs c2","non-CGI: EPIC")) +
  theme(axis.text.x=element_text(angle=45, hjust=1),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1.1)) +
  theme(axis.title.x = element_blank())+ theme(legend.text=element_text(size=8),legend.title = element_blank())
ggsave(file.path(resDir.NOFILTER.C12vsN,"Cumulative percentage for CGI and outsideCGI.pdf"))
#############################################################
# Histogram for except CGI
tmp4 <- tmp3 %>% filter(CGI=="OutsideCGI")
p2 <- ggplot(tmp4,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1.1)) +
  scale_y_continuous(breaks=seq(0, 60, 5), limits=c(0, 60)) + theme(axis.title.x = element_blank()) + theme(legend.text=element_text(size=8),legend.title = element_blank())
ggsave(file.path(resDir.NOFILTER.C12vsN,"Percentage for outsideCGI.pdf"))

library(cowplot)
p3 <- plot_grid(p1, p2, labels = c("a", "b"))
ggsave(file.path(resDir.NOFILTER.C12vsN,"Pairwise percentage for region.pdf"),width = 9,height = 6)

#############################################################
# Histogram for pb distribution of H3K27ME3
annoPb <- read.table("F:/Project/UTUC_BLCA/Annotation/Methylation_probes_annotation_simplified.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
annoPb$Simple_UCSC_RefGene_Name <- toupper(annoPb$Simple_UCSC_RefGene_Name)
MSigDB=clusterProfiler::read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
# tmp <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_entire_DMP_GSEA.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

#"BENPORATH_ES_WITH_H3K27ME3"                 "MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"
Gaby.pathway.gene <- list()
# Gaby.pathway.gene[["BENPORATH_ES_WITH_H3K27ME3"]] <- toupper(strsplit(Gaby.pathway[1,"Genes"]," ")[[1]])
# Gaby.pathway.gene[["MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"]] <- toupper(strsplit(Gaby.pathway[2,"Genes"]," ")[[1]])

Gaby.pathway.gene[["BENPORATH_ES_WITH_H3K27ME3"]] <- toupper(MSigDB[which(MSigDB$ont %in% "BENPORATH_ES_WITH_H3K27ME3"),"gene"])
Gaby.pathway.gene[["MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"]] <- toupper(MSigDB[which(MSigDB$ont %in% "MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"),"gene"])

Gaby.pathway.pb <- list()
Gaby.pathway.pb[["BENPORATH_ES_WITH_H3K27ME3"]] <- feature[rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% Gaby.pathway.gene[[1]]),]),]
Gaby.pathway.pb[["MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"]] <- feature[rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% Gaby.pathway.gene[[2]]),]),]

tmp <- Gaby.pathway.pb[[1]]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- colnames(Gaby.pathway.pb[[1]])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/nrow(Gaby.pathway.pb[[1]])*100,0)

tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage")
tmp3$Class <- factor(rep(c("H3K27ME3","EPIC Array"),each=11),levels = c("H3K27ME3","EPIC Array"))
p3 <- ggplot(tmp3,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_manual(values=c("#0068B5","#F09300")) +
  scale_y_continuous(breaks=seq(0, 70, 5), limits=c(0, 70))
ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"Histograme for pb distribution of BENPORATH_ES_WITH_H3K27ME3.pdf"))

tmp <- Gaby.pathway.pb[[2]]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- colnames(Gaby.pathway.pb[[2]])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/nrow(Gaby.pathway.pb[[2]])*100,0)

tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage")
tmp3$Class <- factor(rep(c("H3K27ME3","EPIC Array"),each=11),levels = c("H3K27ME3","EPIC Array"))
p4 <- ggplot(tmp3,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_manual(values=c("#0068B5","#F09300")) +
  scale_y_continuous(breaks=seq(0, 65, 5), limits=c(0, 65))
ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"Histograme for pb distribution of MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3.pdf"))

####################################################
# 3D plot for DMPs in C1 vs C2
dmp.res <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
df <- dmp.res
df$logPvalue <- -log10(df$adj.P.Val)
df$Class <- rep("Normal",nrow(df))

df[which(df$C1_AVG>=0.4 & df$C2_AVG<=0.2),"Class"] <- "gainMethC1"
df[which(df$C1_AVG<=0.2 & df$C2_AVG>=0.4),"Class"] <- "gainMethC2"

df <- df[,c(7,8,32,33)]
probeC1 <- sample(rownames(df[which(df$Class=="gainMethC1"),]),500)
probeC2 <- rownames(df[which(df$Class=="gainMethC2"),])
Normal <- sample(rownames(df[which(df$Class=="Normal"),]),500)
df <- df[c(probeC1,probeC2,Normal),]

plot_ly(df, x = ~C1_AVG, y = ~C2_AVG, z = ~logPvalue, color = ~Class, colors = c(red,blue,grey),size = I(3)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'C1 betaAVG',showgrid=T,gridcolor="black",gridwidth=3,showline=T,linewidth=5,autotick=F,tick0=0.4,zeroline=F),
                      yaxis = list(title = 'C2 betaAVG',showgrid=T,gridcolor="black",gridwidth=3,showline=T,linewidth=5,autotick=F,tick0=0.2,zeroline=F),
                      zaxis = list(title = 'log10Pvalue',showgrid=T,gridcolor="black",gridwidth=3,showline=T,linewidth=5,zeroline=F)))

scatter3d(x = df$C1_AVG, y = df$C2_AVG, z = df$logPvalue, groups = df$Class, 
          surface.col = c(red,blue,grey),surface=F, ellipsoid = TRUE,
          xlab = "C1 betaAVG", ylab = "C2 betaAVG",
          zlab = "log(adj.pValue)",axis.col = c("black", "black", "black"))

# heatmap for gainMethC1 and gainMethC2
# add annotation for 9p21.3 CDKN2A
tmp <- read.table(file.path("D:/UTUC_BLCA","CNA/GISTIC/1732580/CNA_segment_forGISTIC2.0.all_thresholded.by_genes.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- as.data.frame(t(tmp["CDKN2A",3:37]))
rownames(tmp) <- substr(rownames(tmp),start = 6,stop = 8)
annCol$`focal deletion in CDKN2A` <- as.character(tmp[rownames(annCol),"CDKN2A"])
annCol$`focal deletion in CDKN2A` <- gsub("0",".",annCol$`focal deletion in CDKN2A`)
annCol$`focal deletion in CDKN2A` <- gsub("-1","Loss",annCol$`focal deletion in CDKN2A`)
annCol$`focal deletion in CDKN2A` <- gsub("-2","Del",annCol$`focal deletion in CDKN2A`)
annColors[["focal deletion in CDKN2A"]] <- c("."=white,"Loss"=blue,"Del"=darkblue)
table(annCol$`focal deletion in CDKN2A`,annCol$DMBCluster)
# C1 C2
# .    18  8
# Del   3  3
# Loss  2  1
fisher.test(matrix(c(18,8,5,4),byrow = T,ncol = 2)) #0.6855

load("F:/Project/UTUC_BLCA/Results/hcs.DMBcluster.rda")
index <- c("02T","03T","04T","05T","06T","08T","09T","12T","13T","14T","16T","15T","17T","18T","19T","21T","22T","25T","26T","27T","28T","29T","30T","31T","32T","33T","34T","35T","36T","37T","38T","39T","41T","42T","43T")
df <- c(gainMethInC1,gainMethInC2)
df <- myCombat.C12[df,rownames(annCol)]
pdf(file.path(resDir.NOFILTER.C12vsN,"heatmap_gainMethInC1andC2_14243pbs.pdf"))
aheatmap(as.matrix(df[,index]), Rowv=NA, Colv=dendsort(as.dendrogram(hcs.DMBcluster)), annCol=annCol[index,],annColors = annColors, color=c("blue", "blue", "green", "yellow", "red", "red"), revC=TRUE, fontsize=5,labRow = NA,labCol = NA,cexAnn = 0.5)
invisible(dev.off())

df <- myCombat.C12[gainMethInC2,rownames(annCol)]
hcg <- hclust(distanceMatrix(as.matrix(t(df)), "euclidean"), "ward.D")
pdf(file.path(resDir.NOFILTER.C12vsN,"heatmap_gainMethInC2_34pbs.pdf"),height = 2.5,width = 7)
aheatmap(as.matrix(df[,index]), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs.DMBcluster)), annCol=annCol[index,],annColors = annColors, color=colorRampPalette(brewer.pal(11,"Spectral"))(11)[11:1], revC=TRUE, fontsize=5,labRow = NA,labCol = NA,cexAnn = 0.5)
invisible(dev.off())

df <- myCombat.C12[gainMethInC1,rownames(annCol)]
hcg <- hclust(distanceMatrix(as.matrix(t(df)), "euclidean"), "ward.D")
pdf(file.path(resDir.NOFILTER.C12vsN,"heatmap_gainMethInC1_14209pbs.pdf"),height = 3.5,width = 7)
aheatmap(as.matrix(df[,index]), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs.DMBcluster)), annCol=annCol[index,],annColors = annColors, color=colorRampPalette(brewer.pal(11,"Spectral"))(11)[11:1], revC=TRUE, fontsize=5,labRow = NA,labCol = NA,cexAnn = 0.5)
invisible(dev.off())

rm(df)
rm(probeC1)
rm(probeC2)
rm(Normal)
####################################################
# calculate ssGSEA foe 293 Gaby interested pathway #

load(file.path(resDir.NOFILTER.C12vsN,"myCombat.C12.rda"))
load(file.path(resDir.NOFILTER.C12vsN,"myLoad.rda"))

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
MeTIL <- myCombat.C12[MeTIL.marker,rownames(annCol)[1:35]]
MeTIL <- t(scale(t(MeTIL)))
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL.score <- annTrackScale(indata = pca.MeTIL$rotation[,1], halfwidth = 2, poolsd = F)

tmp <- as.data.frame(pca.MeTIL$rotation[,1]); colnames(tmp) <- "pca.MeTIL"
write.table(tmp,"H:/UTUC_BLCA/Results/pca.MeTIL.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

annCol <- data.frame("SampleName"=myLoad$pd$Sample_Name,"DMBCluster"=myLoad$pd$Sample_Group,row.names = myLoad$pd$Sample_Name,stringsAsFactors = F)
annCol <- annCol[order(annCol$DMBCluster,decreasing = F),]
annCol <- annCol[1:35,]
annCol$MeTIL.score <- MeTIL.score[rownames(annCol)]
annCol <- annCol[,-1]

annColors <- list()
annColors[["DMBCluster"]] <- c("C1"=red,"C2"=blue)
annColors[["MeTIL.score"]] <- bluered(64)
write.table(annCol,file.path(resDir.NOFILTER.C12vsN,"MeTIL in DMBcluster C1 and C2.txt"),sep = "\t",row.names = T,col.names = NA)

tmp1 <- annCol[which(annCol$DMBCluster == "C1"),"MeTIL.score"]
tmp2 <- annCol[which(annCol$DMBCluster == "C2"),"MeTIL.score"]
p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,3) # 0.001262
df <- data.frame("MeTIL"=c(tmp1,tmp2),"DMBCluster"=rep(c("C1","C2"),c(length(tmp1),length(tmp2))))
p <- ggplot(df,aes(x=DMBCluster,y=MeTIL,fill=DMBCluster)) + geom_boxplot() + scale_fill_manual(values = c(red,blue)) + ggtitle("pValue = 0.001") + theme(legend.position="none")
ggsave(file.path(resDir.NOFILTER.C12vsN,"boxplot of MeTIL score between UTUC DMBcluster C1 and C2 .pdf"),width = 3,height = 4)

gaby.interested.pathway <- read.table("F:/Project/UTUC_BLCA/InputData/Gaby_interested_pathway.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
gaby.interested.pathway.list <- list()
for (i in 1:nrow(gaby.interested.pathway)) {
  tmp <- gaby.interested.pathway$PathwayID[i]
  tmp1 <- toupper(MSigDB[which(MSigDB$ont %in% tmp),"gene"])
  tmp2 <- intersect(rownames(dmp.res),rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% tmp1),]))
  gaby.interested.pathway.list[[tmp]] <- tmp2
}
gaby.enrichscore <- gsva(as.matrix(myCombat.C12),gaby.interested.pathway.list,method="ssgsea")
gaby.enrichscore <- as.data.frame(gaby.enrichscore[,rownames(annCol)])
gaby.enrichscore.backup <- gaby.enrichscore

p <- vector()
for (i in 1:nrow(gaby.enrichscore)) {
  tmp <- round(wilcox.test(as.numeric(gaby.enrichscore[i,1:23]),as.numeric(gaby.enrichscore[i,24:35]),"greater")$p.value,4)
  p <- append(p, tmp)
}
rownames(gaby.enrichscore) <- paste0(rownames(gaby.enrichscore),"_p_",p)

plotscore <- as.data.frame(na.omit(standarize.fun(indata = gaby.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotscore)), "euclidean"), "ward.D")
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_DNA_enrichmentscore293_for annotation.pdf"), height=12)
hv <- aheatmap(as.matrix(plotscore), Rowv=dendsort(as.dendrogram(hcg)), Colv=NA, annCol=annCol, annColors=annColors, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=5, cexCol = 0.8, cexRow = 0.4,cexAnn = 0.6)
invisible(dev.off())
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_DNA_enrichmentscore293.pdf"), height=6)
hv <- aheatmap(as.matrix(plotscore), Rowv=dendsort(as.dendrogram(hcg)), Colv=NA, annCol=annCol, annColors=annColors, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=5, cexCol = 0.8, cexRow = 0.4,cexAnn = 0.6,labRow = NA)
invisible(dev.off())

# gaby methylation ssgsea pathways to plot
gaby.plot.pathway <- read.table("F:/Project/UTUC_BLCA/InputData/Gaby_methylation_ssgsea_pathways_to_plot.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
gaby.plot.pathway.list <- list()
for (i in 1:nrow(gaby.plot.pathway)) {
  tmp <- gaby.plot.pathway$PathwayID[i]
  tmp1 <- toupper(MSigDB[which(MSigDB$ont %in% tmp),"gene"])
  tmp2 <- intersect(rownames(dmp.res),rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% tmp1),]))
  gaby.plot.pathway.list[[tmp]] <- tmp2
}
gaby.plot.enrichscore <- gsva(as.matrix(myCombat.C12),gaby.plot.pathway.list,method="ssgsea")
gaby.plot.enrichscore <- as.data.frame(gaby.plot.enrichscore[,rownames(annCol)])
gaby.plot.enrichscore.backup <- gaby.plot.enrichscore

p <- vector()
for (i in 1:nrow(gaby.plot.enrichscore)) {
  tmp <- wilcox.test(as.numeric(gaby.plot.enrichscore[i,1:23]),as.numeric(gaby.plot.enrichscore[i,24:35]),"greater")$p.value
  p <- append(p, tmp)
}
rownames(gaby.plot.enrichscore) <- paste0(rownames(gaby.plot.enrichscore),"_p_",p)

plotscore <- as.data.frame(na.omit(standarize.fun(indata = gaby.plot.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotscore)), "euclidean"), "ward.D")
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_DNA_enrichmentscore_to_plot.pdf"), width = 15,height = 3.5)
mycol <- colorpanel(64,low=cyan,mid = "black",high=peach)
hv <- pheatmap(as.matrix(plotscore),border_color = NA,
               color = mycol,cluster_rows = dendsort(hcg),treeheight_row = 0,cluster_cols = F,
               annotation_col = annCol, annotation_colors = annColors,
               show_rownames = T,show_colnames = F,legend = T)
invisible(dev.off())
############################################################
## Association between the mut.sig (n=30) and MeTIL.score ##
#mut.sig.res <- read.table("F:/Project/UTUC_BLCA/exome UTUC/Results/Mut.sig.snp cluster3 results of DMBcluster MWBcluster RNAcluster.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
mut.sig.res <- read.table("H:/UTUC_BLCA/exome UTUC/Results/Mut.sig.snp cluster3 results of DMBcluster MWBcluster RNAcluster.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

rownames(mut.sig.res) <- gsub("UTUC_","",rownames(mut.sig.res))
mut.sig.res$MeTIL <- annCol[rownames(mut.sig.res),"MeTIL.score"]
mut.sig.res$pcaMeTIL <- pca.MeTIL$rotation[,1][rownames(mut.sig.res)]
mut.sig.res <- as.data.frame(na.omit(mut.sig.res))
mut.sig.res$MWBcluster <- factor(mut.sig.res$MWBcluster,levels = c("C1","C2","C3"))
my_comparisons <- list( c("C1", "C2"), c("C1", "C3"), c("C2", "C3") )
ggviolin(mut.sig.res, x = "MWBcluster", y = "pcaMeTIL", fill = "MWBcluster",
         palette = c("#F09300",lightgreen,"#0068B5"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",hide.ns = F)+ # Add significance levels
  stat_compare_means(label.y = 1)  
#ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"violin plot for association between the mut.sig and MeTIL.score in MWBcluster.pdf"))
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the mut.sig and MeTIL.score in MWBcluster.pdf")

# annCol$NMFbasis <- mut.sig.res[rownames(annCol),"NMFbasis"]
# annCol$NMFbasis[is.na(annCol$NMFbasis)] <- "N/A"
# annColors[["NMFbasis"]] <- c("1"="#9FA0FF","2"="#EE7BFF","3"="#3AB1FF","N/A"=grey)
annCol$MWBcluster <- as.character(mut.sig.res[rownames(annCol),"MWBcluster"])
annCol$MWBcluster[is.na(annCol$MWBcluster)] <- "N/A"
annColors[["MWBcluster"]] <- c("C1"="#F09300","C2"=lightgreen,"C3"="#0068B5","N/A"=grey)

mut.sig.res$MWBSigCluster <- ifelse(mut.sig.res$MWBcluster == "C1","Sig.16","Others")
mut.sig.res$MWBSigCluster <- factor(mut.sig.res$MWBSigCluster,levels = c("Sig.16","Others"))
annColors[["NMFSigCluster"]] <- annColors[["MWBSigCluster"]] <- c("Sig.16"=jco[2],"Others"=jco[1])
ggviolin(mut.sig.res, x = "MWBSigCluster", y = "pcaMeTIL", fill = "MWBSigCluster",
         palette = c(jco[2],jco[1]),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Sig.16","Others")), label = "p.signif",hide.ns = F) # Add significance levels
#ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"violin plot for association between the Sig16.Others and MeTIL.score in MWBcluster.pdf"))
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the Sig16.Others and MeTIL.score in MWBcluster.pdf")

mut.sig.res$MWBSigCluster <- ifelse(mut.sig.res$MWBcluster == "C2","Sig.1","Others")
mut.sig.res$MWBSigCluster <- factor(mut.sig.res$MWBSigCluster,levels = c("Sig.1","Others"))
annColors[["NMFSigCluster"]] <- annColors[["MWBSigCluster"]] <- c("Sig.1"=jco[2],"Others"=jco[1])
ggviolin(mut.sig.res, x = "MWBSigCluster", y = "pcaMeTIL", fill = "MWBSigCluster",
         palette = c(jco[2],jco[1]),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Sig.1","Others")), label = "p.signif",hide.ns = F) # Add significance levels
#ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"violin plot for association between the Sig1.Others and MeTIL.score in MWBcluster.pdf"))
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the Sig1.Others and MeTIL.score in MWBcluster.pdf")

mut.sig.res$MWBSigCluster <- ifelse(mut.sig.res$MWBcluster == "C3","Sig.13","Others")
mut.sig.res$MWBSigCluster <- factor(mut.sig.res$MWBSigCluster,levels = c("Sig.13","Others"))
annColors[["NMFSigCluster"]] <- annColors[["MWBSigCluster"]] <- c("Sig.13"=jco[2],"Others"=jco[1])
ggviolin(mut.sig.res, x = "MWBSigCluster", y = "pcaMeTIL", fill = "MWBSigCluster",
         palette = c(jco[2],jco[1]),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Sig.13","Others")), label = "p.signif",hide.ns = F) # Add significance levels
#ggsave(filename = file.path(resDir.NOFILTER.C12vsN,"violin plot for association between the Sig13.Others and MeTIL.score in MWBcluster.pdf"))
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the Sig13.Others and MeTIL.score in MWBcluster.pdf")

fisher.test(table(mut.sig.res$MWBSigCluster,mut.sig.res$MeTIL>0))
# FALSE TRUE
# Sig.13     7    7
# Others     8    6





mut.wt <- read.table("F:/Project/UTUC_BLCA/exome UTUC/Results/mutation.snp.signature.weightMatrix.bydeconstructSigs.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
#mut.wt <- read.table("H:/UTUC_BLCA/exome UTUC/Results/mutation.snp.signature.weightMatrix.bydeconstructSigs.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

rownames(mut.wt) <- paste0(0,sapply(strsplit(rownames(mut.wt),"_"),"[",2))
rownames(mut.wt) <- substr(rownames(mut.wt),start = nchar(rownames(mut.wt))-2,stop = nchar(rownames(mut.wt)))
mut.sig.res$sig13 <- mut.wt[rownames(mut.sig.res),"Signature.13"]
mut.sig.res$sig16 <- mut.wt[rownames(mut.sig.res),"Signature.16"]
mut.sig.res$sig1 <- mut.wt[rownames(mut.sig.res),"Signature.1"]

cor.test(mut.sig.res$pcaMeTIL,mut.sig.res$sig13,method = "spearman")
cor.test(mut.sig.res$pcaMeTIL,mut.sig.res$sig16,method = "spearman") #0.07122
cor.test(mut.sig.res$pcaMeTIL,mut.sig.res$sig1)

# METIL score bwteen ZNF36_family mutation
#UTUC.annotation <- read.table("F:/Project/UTUC_BLCA/InputData/UTUC_annotation.txt",sep = "\t",check.names = F,header = T,row.names = 1,stringsAsFactors = F)
UTUC.annotation <- read.table("H:/UTUC_BLCA/InputData/UTUC_annotation.txt",sep = "\t",check.names = F,header = T,row.names = 1,stringsAsFactors = F)
tmp <- UTUC.annotation$ZFP36_family; names(tmp) <- rownames(UTUC.annotation)
wilcox.test(MeTIL.score[names(tmp[tmp == "1"])],MeTIL.score[names(tmp[tmp == "0"])],na.rm=T) #0.724
# ##############################################
# # try customized GSEA for methylation probes #
# MSigDB <- read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
# pb2gene <- data.frame("probe"=annoPb[1:865859,"Name"],"gene"=annoPb[1:865859,"Simple_UCSC_RefGene_Name"])
# pb2gene <- na.omit(pb2gene)
# save(pb2gene,file = "C:/Users/Sugus/Desktop/pb2gene.rda")
# save(MSigDB,file = "C:/Users/Sugus/Desktop/MSigDB.rda")
# #MSigDB.pb <- merge(MSigDB,pb2gene,by="gene",all.x=T) #DONE
# 
# dmp.res <- read.table(file.path(resDir.NOFILTER.C12vsN,"ChAMP_RES_NOFILTER.C12N_C1vsC2_DMP_addFeatures.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1,header = T)
# pbList <- -dmp.res[which(dmp.res$C1_AVG >= 0.4 & dmp.res$C2_AVG <= 0.2 & dmp.res$adj.P.Val < 0.05),]$deltaBeta
# names(pbList) <- gainMethInC1
# pbList <- sort(pbList,decreasing = T)
# save(pbList,file="C:/Users/Sugus/Desktop/pbList.rda")
# GSEA.DMPs.C1vsC2 <- GSEA(geneList = pbList,TERM2GENE=MSigDB.pb,seed = T,verbose=F)

# pcaMetil score between FGFR3 mutation
tmp <- UTUC.annotation$FGFR3_mut; names(tmp) <- rownames(UTUC.annotation)
tmp1 <- na.omit(pca.MeTIL$rotation[,1][names(tmp[tmp == "1"])])
tmp2 <- na.omit(pca.MeTIL$rotation[,1][names(tmp[tmp == "0"])])
tmp <- data.frame(FGFR3 = rep(c(1,0),c(length(tmp1),length(tmp2))),
                  pcaMeTIL = c(as.numeric(tmp1),as.numeric(tmp2)),
                  stringsAsFactors = F)
tmp$FGFR3 <- factor(ifelse(tmp$FGFR3 == 0,"Wild","Mutated"))
ggviolin(tmp, x = "FGFR3", y = "pcaMeTIL", fill = "FGFR3",
         palette = c(cherry,lightgrey),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Mutated","Wild")), label = "p.signif",hide.ns = F) # Add significance levels
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the curated FGFR3 mutation and MeTIL.score p=9.118e-05.pdf")

wilcox.test(as.numeric(tmp1),as.numeric(tmp2)) #9.118e-05

# create violine plot for MeTIL score between DMBcluster C1 and C2
tmp <- UTUC.annotation$`Cluster DNA methylation no filter without normal`; names(tmp) <- rownames(UTUC.annotation)
tmp1 <- pca.MeTIL$rotation[,1][names(tmp[tmp == "C1"])]
tmp2 <- pca.MeTIL$rotation[,1][names(tmp[tmp == "C2"])]
tmp <- data.frame(DMBCluster = rep(c("C1","C2"),c(length(tmp1),length(tmp2))),
                  pcaMeTIL = c(as.numeric(tmp1),as.numeric(tmp2)),
                  stringsAsFactors = F)
ggviolin(tmp, x = "DMBCluster", y = "pcaMeTIL", fill = "DMBCluster",
         palette = c(red,blue),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("C1","C2")), label = "p.signif",hide.ns = F) # Add significance levels
ggsave("H:/UTUC_BLCA/Figures/violin plot for association between the DMBCluster and MeTIL.score.pdf")

wilcox.test(as.numeric(tmp1),as.numeric(tmp2)) #0.009174

############################################################
### DMP analysis between FGFR3 curated mutation and wild ###
myLoad$pd$FGFR3_curated <- UTUC.annotation[myLoad$pd$Sample_Name,"FGFR3_mut"]
samples.FGFR3_curated <- rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut %in% c("1","0") & UTUC.annotation$`DNA Methylation data` == "yes"),])
pheno.FGFR3_curated <- ifelse(UTUC.annotation[samples.FGFR3_curated,"FGFR3_mut"] == "1","Mutated","Wild")
myDMP.FGFR3_curated <-   champ.DMP(beta=myCombat[,samples.FGFR3_curated],
                         pheno=pheno.FGFR3_curated,
                         arraytype = "EPIC",adjPVal = 0.05); 
write.table(myDMP.FGFR3_curated$Wild_to_Mutated,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.FGFR3_curated.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
save(myDMP.FGFR3_curated,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.FGFR3_curated.rda"))

myDMR.FGFR3_curated <- champ.DMR(beta=myCombat[,samples.FGFR3_curated],
                                 pheno=pheno.FGFR3_curated,
                                 method="DMRcate",
                                 arraytype = "EPIC",
                                 cores = 1);
write.table(myDMR.FGFR3_curated$DMRcateDMR,file = file.path(resDir.NOFILTER.C12vsN,"myDMR.FGFR3_curated.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
save(myDMR.FGFR3_curated,file = file.path(resDir.NOFILTER.C12vsN,"myDMR.FGFR3_curated.rda"))

# draw DMR
phen.col <- ifelse(pheno.FGFR3_curated == "Mutated",cherry,grey); names(phen.col) <- pheno.FGFR3_curated
pdf(file.path(resDir.NOFILTER.C12vsN,"DMR_FGFR3.pdf"),width = 10,height = 8)
DMR.plot(ranges=myDMR.FGFR3_curated$DMR.ranges, 
         dmr=49, CpGs=myCombat[,samples.FGFR3_curated], what="Beta", arraytype = "EPIC",
         phen.col=phen.col, 
         genome="hg19", samps=1:31)
invisible(dev.off())

# pick up probes within DMR
RR <-as.data.frame (myDMR.FGFR3_curated$DMR.ranges)
RR$DMRID <-rownames(RR)
row.names(RR) = NULL
RR$DMRNO <- rownames (RR)
row.names(RR) = RR$DMRID
RR$DMRID = NULL
RR <-RR[order(RR$minfdr), , drop = FALSE]
cgID <- as.data.frame(myDMR.FGFR3_curated$dmrcoutput$input)
DMRNUM <- readline(prompt = "What is your DMR Number:")
assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))
assign(paste0("DMR_",DMRNUM,"_probelist"), subset(cgID, cgID$CHR==assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))$seqnames & cgID$pos>assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))$start-1 & cgID$pos<assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))$end+1))                   
                   
# GSEA                   
myGSEA.FGFR3_curated <- champ.GSEA(beta=myCombat[,samples.FGFR3_curated],
                                   DMP=myDMP.FGFR3_curated[[1]], 
                                   DMR=myDMR.FGFR3_curated,pheno = pheno.FGFR3_curated , 
                                   arraytype="EPIC",adjPval=0.05, method="fisher");
write.table(myGSEA.FGFR3_curated$DMP,file = file.path(resDir.NOFILTER.C12vsN,"GSEA.DMP.FGFR3_curated.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(myGSEA.FGFR3_curated$DMR,file = file.path(resDir.NOFILTER.C12vsN,"GSEA.DMR.FGFR3_curated.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
save(myGSEA.FGFR3_curated,file = file.path(resDir.NOFILTER.C12vsN,"myGSEA.FGFR3_curated.rda"))

# heatmap of REACTOME_NEGATIVE_REGULATION_OF_FGFR_SIGNALING and REACTOME_FGFR_LIGAND_BINDING_AND_ACTIVATION
FGFR3.dmr.list <- list("REACTOME_NEGATIVE_REGULATION_OF_FGFR_SIGNALING" = c("SPRY2", "FRS2", "KLB", "FGF1", "FGF3", "FGF6", "FGF7", "FGF8", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "FGF20", "FGF22", "GRB2", "PPP2CB", "RPS27A", "SRC", "BRAF", "FGF23", "FGF17", "KL", "FGF19"),
                       "REACTOME_FGFR_LIGAND_BINDING_AND_ACTIVATION" = c("KLB", "FGF1", "FGF3", "FGF6", "FGF7", "FGF8", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "FGF20", "FGF22", "FGF23", "FGF17", "KL", "FGF19"),
                       "BENPORATH_PRC2_TARGETS"=c("MAL","NPAS4","TP73","STXBP6","SLC10A4","FOXL1","SCN4B","FBLN7","KIAA1199",
                                                   "INA","IKZF3","SIX1","GRIA2","B4GALNT2","CCDC140","KL","PAX1","DCHS2",
                                                   "PMP22","CYP27B1","TRIM9","PTGER4","NKX6-2","FAM150A","TBX21","ZADH2",
                                                   "OLFML2B","PDX1","TBX5","SOX7","SIX2","CYP24A1","MT1B","NEFL","SLIT2",
                                                   "SRD5A2","INSRR","FLI1","ST8SIA2","HOXC5","HOXD3","PROK2","LTK","FAM5C",
                                                   "DKK1","USH1G","BARX2","COLEC12","IGSF21","OSR1","SLC30A2","NPY1R",
                                                   "FGF3","EFNA3","NPR3","SEMA6D","TLX2","ZMYND15","SYT12","WT1","CNNM1",
                                                   "ICAM5","PAX6","CDH23","PDGFRA","FAM81A","GHR","ROBO3","GATA6","PTGER3",
                                                   "VDR","MCOLN3","MT1A","KLF4","BHLHE23","NELL1","ALX3","DGKG","NDUFA4L2","ZNF436",
                                                   "HOXD12","FGF20","PITX3","RYR3","ADRB3","HOXB2","FOXD2","NAGS","BNC1",
                                                   "WNT10B","EGR3","ZIC1","GIMAP5","CHST8","BHLHE22","GSX2","PARM1",
                                                   "CRTAC1","SERTM1","CIDEA","EFNA1","SLC35F3","KCNAB1","DGKI","LYSMD2",
                                                   "CBR3","TRIM36","IL7","ADCY8","HMX3","DMRT2","PTPRT","NEFM","SLC24A4",
                                                   "DSC3","FOXA2","GATA2","FBP1","RSPO1","HOXC11","DLX4","PAX2","RBP4",
                                                   "CYP26A1","NOL4","ZNF503","LRFN5","ZEB2","FRMD3","SHH","NEUROG3",
                                                   "TMEM132E","SLITRK1","PTH2","RIPK3","CNRIP1","OTX1","NEUROD2","VAX1",
                                                   "TMEM30B","GPR88","TET2","CTNND2","ASTN1","HS6ST3","PAPPA","HTR7",
                                                   "EN1","PKP1","NKX2-1","GRM7","GBX2","PLEC","FEZ1","ADARB2","MKX",
                                                   "ARL9","PAX9","GDF7","MESP1","TBR1","LGR5","ASCL2","DCC","SPOCK3",
                                                   "SLC1A2","CACNA1G","TLX1","CA10","ECEL1","RASSF5","NPNT","SLC6A3",
                                                   "RASGRF1","SLC27A2","FAM43B","GRIK1","DHH","IRX3","TBX2","HAND2",
                                                   "TRH","SLC30A4","GRID1","FAM5B","UCN","ADCYAP1","HNF1B","HHIP","DDAH1",
                                                   "WNT6","GPC5","BTG2","NT5C1A","KCNK12","COL25A1","PTF1A","SLCO5A1",
                                                   "CACNA1D","LMX1B","NTRK2","HOXD1","DIO3","ZFHX3","TBX1","IRX4","GDNF",
                                                   "CRYBA2","KCNMA1","WNT1","B4GALNT1","RGS9BP","PGR","WNT3A","EOMES",
                                                   "SIM2","ATOH1","DSCAML1","CD8A","IRX5","SSTR2","LAYN","HLX","METRNL",
                                                   "ADAP2","GABRA2","GSC","SLC9A3","MSC","NFIX","SIX6","ITPKA","NKX3-2",
                                                   "KIRREL3","COMP","FEZF2","ADAMTS18","KCNK4","TPPP3","SOX14","HOXD8",
                                                   "BARX1","GJB2","ADAMTS15","DNAJC22","CSMD3","PTPRU","NTN1","LHX5",
                                                   "GRIK3","GABRA4","PAX3","CHRD","UCP1","PRKG1","IGF2-AS","WNT16",
                                                   "UNC5C","DUOX1","KLHL35","NPTX1","GPR12","TSLP","KCNC2","SFRP5",
                                                   "KCND3","KCNH1","SLITRK3","SHOX2","CASZ1","LBX1","NKX2-3","LRRTM1",
                                                   "OTOP3","FAM19A4","PCDH17","CALCA","GRIN3A","HPSE2","FAM84A","SHISA6",
                                                   "PRAC","GATA3","GSC2","VAX2","SLC6A5","AQP5","SIX3","TMEM59L",
                                                   "NKX2-8","ZIC4","SSTR1","SLCO2A1","FOXF1","HOXC4","FEV","NAV2",
                                                   "HOXB3","MT1H","FOXL2","OPRD1","HES2","HHEX","SCD5","OTOP2","CXCL16",
                                                   "TFAP2E","C2CD4A","ERBB4","ADRB1","LHX6","HS3ST3B1","HOXB8","DUOXA2",
                                                   "FIGLA","ANKRD27","FLJ32063","EPAS1","NEUROG1","ADRA2A","CSMD1",
                                                   "TMOD2","MSX1","KAZALD1","LPPR1","PITX1","KCNH3","FLJ45983","FAM163A",
                                                   "EN2","ESPN","PRDM12","HOXC6","NR2F2","TTYH1","COL9A2","TLL1","MNX1",
                                                   "CDH7","DUOXA1","NRG1","KCNK2","OCA2","GUCY1A3","KY","CYP26C1","FLRT2",
                                                   "SFRP1","GALR2","HPCAL4","EPHB1","EPHA5","HOXB6","LPL","HOXD4","THBD",
                                                   "PYY","VSX1","BHLHE41","WNT7A","NTRK1","TBX3","DLX3","TRIM67","LRRC71",
                                                   "DUOX2","ALX4","FOXG1","KCNV1","MAPK4","NCAM1","SORCS1"),
                       BENPORATH_EED_TARGETS=c("SELS","HOXC5","SLC30A4","HOXD4","HAP1","DUOX2","NAV2",
                                               "SLC9A3","MSX2","GALR2","TBX21","TP73","FLJ45983","MAF",
                                               "RNF152","PTGER4","MNX1","RIPK3","FBP1","DKK1","NEUROG1",
                                               "SIX6","HOPX","PLXNC1","ADAP2","EVX1","ADM","NKX2-1",
                                               "HIST1H3I","KY","DNAJC22","INTU","TNFRSF25","NEUROD2",
                                               "KCNAB1","GSTA4","SHISA2","KCNQ1","FEZ1","IL7","SIX5",
                                               "PYY","FAM150A","NRG1","GPR12","CHST8","HOXD3","FGF20",
                                               "CCDC140","FOXD2","DLX5","NTRK1","DCC","NELL1","BHLHE41",
                                               "SNAP91","ZADH2","CTNND2","GRM3","OPRD1","ADCYAP1","GATA6",
                                               "IRX3","FLI1","KAZALD1","LAYN","BHLHE23","DHRS13","DDAH1",
                                               "TRH","VWA3B","COL12A1","DIO3","GPR88","KCNH3","PAX3","HNF1B",
                                               "KCNK2","GUCY1A3","MKX","BTG4","ALX3","GRM7","HMX3","IRX4",
                                               "RGS9BP","PLXND1","CXCL12","NOL4","KIRREL3","EMX2","SIM1",
                                               "PLA2G7","CRTAC1","ASCL2","GPR6","KCNH1","FOXA2","KIAA1199","PAPOLA",
                                               "B4GALNT1","CCBE1","EN2","TPPP3","CEBPD","EOMES","NTN1","FAM19A4","KL",
                                               "NKX3-2","ZNF503","GDF7","SSTR1","NEUROG3","NPTX1","FAM84A","HOXB6","GFI1",
                                               "CASZ1","NEFM","ANKRD27","FEZF2","UCP1","SLC6A5","LRFN5","DLK1","SLC22A3","PGR",
                                               "RFX6","SOX7","HOXC6","TBX1","NAGS","IGF2-AS","CSMD3","DGKG","GRID1","DLX4",
                                               "CACNA1D","IKZF3","SLC35F3","MED31","LRFN2","NKPD1","EFNA3","TWIST1","KCNK4","HOXB3",
                                               "TGFB2","INSRR","B4GALNT2","FAM5B","SSTR2","FOS","GBX2","GSX2","SIX2","MT1A","LMO1",
                                               "CLIC5","HOXB8","SLCO2A1","CLIP4","FGF3","ADRB1","SHOX2","GSC","HTR1B","FLRT2","VAX2",
                                               "SFRP5","PTF1A","LRRC71","SOX21","GATA3","HOXD11","SHISA6","NCAM1","EPHB1","ADCY8",
                                               "CHRD","MSX1","DSCAML1","HOXC9","FBXO39","COL19A1","GABRA4","BNC1","HTR7","CYP26C1",
                                               "IGSF21","NDUFA4L2","ARL9","ITPKA","TMEM132E","FAM163A","EBF1","PRAC","SLC30A2","KCNC2",
                                               "BHLHE22","WNT1","SORCS2","PARM1","SOX1","COL9A2","NPY1R","DUOXA1","FOXF2","RYR3",
                                               "PTH2","WT1","VDR","CSNK1E","INHBB","ALX4","BTG2","HOXA3","CTGF","SEMA6D",
                                               "TAP1","EYA4","NDUFA13","FOXL1","BMPR1B","GDA","SRRM2","NPNT","CUEDC1","NR2F2","TSSK6",
                                               "OTOP2","SLITRK1","DSC3","CRYBA2","EPHA5","HOXA9","SFRP1","HOXD8","P2RX5","TMOD2",
                                               "MRPS18B","NPAS4","GATA2","ESPN","FOXF1","LHFPL3","MT1B","MT1H","QRFPR",
                                               "VSX1","ZEB2","GRIK1","SERTM1","HCG9","HOXC11","FEV","UNC5C","ESYT3","WNT7A",
                                               "GRIK3","PLEC","EGR3","VSTM4","TFAP2B","NT5C1A","NKX2-8","ADAMTS15","SLITRK3",
                                               "TAC1","PRKG1","KLF4","EFNA1","ZNF436","ETV7","HTR1E","ASTN1","PCDH10",
                                               "NPR3","PROK2","LTK","EPS8","TLL1","ZCCHC14","STEAP4","LOX","LGR5","EGFR",
                                               "HES2","SLC18A3","CNRIP1","IL20RA","PAPPA","SLCO5A1","VSIG2","GPC5",
                                               "MAPK4","CALCA","PARD3","OLIG1","GHR","PRDM12","DMRT2","WNT3A","CDH7",
                                               "CYP26A1","HOXD12","HOXB2","SPOCK3","ERBB4","RAB32","SF3A2","LYPD1","BAIAP2",
                                               "MICB","NXPH4","TBX2","TRIM9","GPR126","ICAM5","TBR1","PAX2","CIDEA","HIST1H4L",
                                               "COMP","PDGFRA","FOXQ1","PTGER3","INA","WNT6","FAM81A","TTPA","SRD5A2","RBP4",
                                               "HOXA7","OCA2","GINS1","GDNF","SLCO4A1","SORCS1","FOXG1","CTSD","ADRB3","B4GALT6",
                                               "DLX3","RASSF5","DCHS2","EPHX3","SCN4B","LHX5","GRIN3A","LRRTM1","SYT10","FZD8",
                                               "ZIC2","PPP1R10","CREG2","RGS6","LYSMD2","TUBA4B","EPAS1","COL25A1","NR2E1","DSP",
                                               "TGFA","IRX5","GIMAP5","HIST1H2AL","KCNMA1","ZIC4","HHEX","PDE10A","GFRA1","NKX6-2",
                                               "PTPRT","OTX1","FNDC1","PMP22","DHH","SDC2","USH1G","NEFL","LBX1","PITX1","IGFBP3",
                                               "GABRB2","PLCB1","DFNA5","EN1","NKAIN2","OTOP3","DCLK1","NGB","CYP1B1","ST8SIA2",
                                               "AIM1","KCNV1","SYT12","TBX3","TRIM67","ZIC1","GABRA2","MOB2","CHN2","PAX6","ROBO3",
                                               "RBM24","EYA2","CNNM1","KLHL35","SHROOM1","MSC","CD8A","MTSS1","MESP1","BCL2L10",
                                               "DUOX1","ACVR1C","FBXL14","SCD5","ADRA2A","OSR1","FUCA2","HOXD1","MCOLN3","CBR3",
                                               "ASAH1","BARX1","THBD","ADARB2","CXCL16","HHIP","GRIA2","GAD1","SLIT2","MAGI2",
                                               "CYP27B1","WNT10B","SOD3","HIC1","SLC24A4","PAPL","SOX5","SLC27A2","MEOX2",
                                               "CYP24A1","TCF21","HAND1","RASGRP1","TMEM30B","PAX1","BACH2","RASGRF1","SHH",
                                               "SCIN","TBX5","DUOXA2","LOC283392","UCN","HLX","SLC6A3","HIST1H4C","HPCAL4","PKP1",
                                               "LMX1B","NKX2-3","SFRP4","TRIM36","FLJ32063","CSMD1","LPL","PAX9","MAL","SOX14",
                                               "HOXC13","HAND2","HOXA4","DGKI","CACNA1G","NRIP3","FAM43B","STXBP6","TMEM59L","TLX1",
                                               "ATOH1","KCND3","LPPR1","NFIX","TUBA4A","ALOX5","COLEC12","FAM5C","PTPRU","HS6ST3","SLC1A2","ECEL1",
                                               "METRNL","ZFHX3","WNT16","LHX6","FOXL2","IGF2","TCF15","KCNK12","PCDH17","NEURL","SIM2",
                                               "GALNT10","RSPO1","HPSE2","PDX1","STC2","TET2","TMEM26","AHR","HS3ST3B1","PITX3","TRHDE","CALCR",
                                               "SIX1","ANKRD9","FBLN7","ADAMTS18","AQP5","BARX2","GJB2","SLC35D3","TTYH1","TFAP2E","HTRA1","VAX1",
                                               "SIX3","SLC8A3","CA10","HOXC4","PPP1R14C","OLFML2B","ZMYND15","FRMD3",
                                               "C2CD4A","FERD3L","CDH23","NTRK2","TRIM7","FIGLA","TSLP","SLC10A4","TLX2","GSC2"))
                       
                       # "KONDO_EZH2_TARGETS"=c("ZEB1","ADAM12","FTSJD1", "EPHA5", "MCPH1","WWTR1","APBB2","NCF4",
                       #                      "LMCD1","COL1A1","KLF3","HAS3","HPGD","TMTC1","EGFR","MMP1","EPB41L3","SGK1",
                       #                      "ITGBL1","CCDC50","UBXN10","SPOCK1","ARHGEF3","ADCY7","CTGF","IRX5","MATN3",
                       #                      "HSD17B2","GRK4","LOC100294362","EDIL3","SHISA3","MAP4K4","RAB3B","GBA2",
                       #                      "FAR2","NAV2","ANO1","AP3B2","APP","GPC6","COMMD10","OSR1","GLP1R","IQCA1",
                       #                      "ELMO1","DNER","ST3GAL6","PCK1","TACR1","FLI1","SERPINE1","NEK3","NPAS2",
                       #                      "KCNMA1","MYLK","GULP1","ZNF766","RNF144B","NT5E","NNMT","CD59","CPE",
                       #                      "LOC100128788","FAM43A","IGFBP3","NIP7","MSX2","CCDC56","NEBL","ABI3BP",
                       #                      "YPEL5","C1S","TPM3","CTSB","ATP8B1","ATXN3","RGMB","KAAG1","PRSS35","SCN2A",
                       #                      "TNFRSF21","TRIM2","NRP1","GAL3ST1","ZNF567","TMEM45A","AOX1","SORL1",
                       #                      "SPRY4","PCLO","ITGA6","RBP4","CAP2","SYT6","CNIH3","SLC2A10"),
                       #"BENPORATH_SUZ12_TARGETS"=c("MAL","NPAS4","TP73","STXBP6","SLC10A4","FOXL1","SCN4B","FBLN7","KIAA1199","INA","IKZF3","SIX1","GRIA2","B4GALNT2","CCDC140","KL","PAX1","DCHS2","PMP22","CYP27B1","TRIM9","PTGER4","NKX6-2","FAM150A","TBX21","ZADH2","OLFML2B","PDX1","TBX5","SOX7","SIX2","CYP24A1","MT1B","NEFL","SLIT2","SRD5A2","INSRR","FLI1","ST8SIA2","HOXC5","HOXD3","PROK2","LTK","FAM5C","DKK1","USH1G","BARX2","COLEC12","IGSF21","OSR1","SLC30A2","NPY1R","FGF3","EFNA3","NPR3","SEMA6D","TLX2","ZMYND15","SYT12","WT1","CNNM1","ICAM5","PAX6","CDH23","PDGFRA","FAM81A","GHR","ROBO3","GATA6","PTGER3","VDR","MCOLN3","MT1A","KLF4","BHLHE23","NELL1","ALX3","DGKG","NDUFA4L2","ZNF436","HOXD12","FGF20","PITX3","RYR3","ADRB3","HOXB2","FOXD2","NAGS","BNC1","WNT10B","EGR3","ZIC1","GIMAP5","CHST8","BHLHE22","GSX2","PARM1","CRTAC1","SERTM1","CIDEA","EFNA1","SLC35F3","KCNAB1","DGKI","LYSMD2","CBR3","TRIM36","IL7","ADCY8","HMX3","DMRT2","PTPRT","NEFM","SLC24A4","DSC3","FOXA2","GATA2","FBP1","RSPO1","HOXC11","DLX4","PAX2","RBP4","CYP26A1","NOL4","ZNF503","LRFN5","ZEB2","FRMD3","SHH","NEUROG3","TMEM132E","SLITRK1","PTH2","RIPK3","CNRIP1","OTX1","NEUROD2","VAX1","TMEM30B","GPR88","TET2","CTNND2","ASTN1","HS6ST3","PAPPA","HTR7","EN1","PKP1","NKX2-1","GRM7","GBX2","PLEC","FEZ1","ADARB2","MKX","ARL9","PAX9","GDF7","MESP1","TBR1","LGR5","ASCL2","DCC","SPOCK3","SLC1A2","CACNA1G","TLX1","CA10","ECEL1","RASSF5","NPNT","SLC6A3","RASGRF1","SLC27A2","FAM43B","GRIK1","DHH","IRX3","TBX2","HAND2","TRH","SLC30A4","GRID1","FAM5B","UCN","ADCYAP1","HNF1B","HHIP","DDAH1","WNT6","GPC5","BTG2","NT5C1A","KCNK12","COL25A1","PTF1A","SLCO5A1","CACNA1D","LMX1B","NTRK2","HOXD1","DIO3","ZFHX3","TBX1","IRX4","GDNF","CRYBA2","KCNMA1","WNT1","B4GALNT1","RGS9BP","PGR","WNT3A","EOMES","SIM2","ATOH1","DSCAML1","CD8A","IRX5","SSTR2","LAYN","HLX","METRNL","ADAP2","GABRA2","GSC","SLC9A3","MSC","NFIX","SIX6","ITPKA","NKX3-2","KIRREL3","COMP","FEZF2","ADAMTS18","KCNK4","TPPP3","SOX14","HOXD8","BARX1","GJB2","ADAMTS15","DNAJC22","CSMD3","PTPRU","NTN1","LHX5","GRIK3","GABRA4","PAX3","CHRD","UCP1","PRKG1","IGF2-AS","WNT16","UNC5C","DUOX1","KLHL35","NPTX1","GPR12","TSLP","KCNC2","SFRP5","KCND3","KCNH1","SLITRK3","SHOX2","CASZ1","LBX1","NKX2-3","LRRTM1","OTOP3","FAM19A4","PCDH17","CALCA","GRIN3A","HPSE2","FAM84A","SHISA6","PRAC","GATA3","GSC2","VAX2","SLC6A5","AQP5","SIX3","TMEM59L","NKX2-8","ZIC4","SSTR1","SLCO2A1","FOXF1","HOXC4","FEV","NAV2","HOXB3","MT1H","FOXL2","OPRD1","HES2","HHEX","SCD5","OTOP2","CXCL16","TFAP2E","C2CD4A","ERBB4","ADRB1","LHX6","HS3ST3B1","HOXB8","DUOXA2","FIGLA","ANKRD27","FLJ32063","EPAS1","NEUROG1","ADRA2A","CSMD1","TMOD2","MSX1","KAZALD1","LPPR1","PITX1","KCNH3","FLJ45983","FAM163A","EN2","ESPN","PRDM12","HOXC6","NR2F2","TTYH1","COL9A2","TLL1","MNX1","CDH7","DUOXA1","NRG1","KCNK2","OCA2","GUCY1A3","KY","CYP26C1","FLRT2","SFRP1","GALR2","HPCAL4","EPHB1","EPHA5","HOXB6","LPL","HOXD4","THBD","PYY","VSX1","BHLHE41","WNT7A","NTRK1","TBX3","DLX3","TRIM67","LRRC71","DUOX2","ALX4","FOXG1","KCNV1","MAPK4","NCAM1","SORCS1"),
for (i in 1:4) {
  tmp <- FGFR3.dmr.list[[i]]
  tmp <- toupper(tmp)
  tmp <- rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% tmp),])
  FGFR3.dmr.list[[i]] <- tmp
}             
FGFR3.dmr.enrichscore <- gsva(as.matrix(myCombat[,FGFR3.anno$SampleID]),FGFR3.dmr.list,method="gsva") # better than ssgsea
mycol <- colorpanel(64,low=darkblue,mid = "black",high=gold)
FGFR3.dmr.enrichscore <- standarize.fun(FGFR3.dmr.enrichscore,halfwidth = 1)
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_FGFR3_dmr_enrichmentscore_to_plot.pdf"), width = 10,height = 2.5)
hv <- pheatmap(as.matrix(FGFR3.dmr.enrichscore),border_color = NA,
               color = mycol,cluster_rows = F,treeheight_row = 0,cluster_cols = F,
               annotation_col = annCol[colnames(FGFR3.dmr.enrichscore),c("FGFR3","MeTIL.score")], annotation_colors = annColors,
               show_rownames = T,show_colnames = F,legend = T,
               gaps_row = 2)
invisible(dev.off())

wilcox.test(FGFR3.dmr.enrichscore[3,rownames(FGFR3.anno[which(FGFR3.anno$FGFR3 == "Mutated"),])],
            FGFR3.dmr.enrichscore[3,rownames(FGFR3.anno[which(FGFR3.anno$FGFR3 == "Wild"),])])

wilcox.test(FGFR3.dmr.enrichscore[4,rownames(FGFR3.anno[which(FGFR3.anno$FGFR3 == "Mutated"),])],
            FGFR3.dmr.enrichscore[4,rownames(FGFR3.anno[which(FGFR3.anno$FGFR3 == "Wild"),])])

# heatmap of promoter methylation of FGFR3, GATA3, MIR200A, KRT15 and KRT
# tmp <- read.delim(file.path(resDir.NOFILTER.C12vsN,"promoter_probes_5cancer_genes.txt"),sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
# tmp <- tmp[intersect(rownames(tmp),rownames(myCombat)),]
# tmp <- tmp[order(tmp$Simple_UCSC_RefGene_Name),]
# plotdata <- as.data.frame(myCombat[rownames(tmp),FGFR3.anno$SampleID])
# plotdata[plotdata >0.3] <- 1
# plotdata[plotdata >0.3 & plotdata<=0.7] <- 0
# plotdata[plotdata <=0.3] <- 0
# 
# #plotdata$gene <- tmp$Simple_UCSC_RefGene_Name
# #plotdata <- apply(plotdata[,setdiff(colnames(plotdata),"gene")], 2, function(x) tapply(x, INDEX=factor(plotdata$gene), FUN=mean, na.rm=TRUE))
# plotdata <- standarize.fun(plotdata,halfwidth = 1)
# annColors[["Gene"]] <- c("FGFR3" = sun,
#                          "KRT15" = blue,
#                          "KRT5" = yellow,
#                          "MIR200A" = seagreen)
# pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_FGFR3_cancer genes.pdf"), width = 7,height = 5)
# hv <- pheatmap(plotdata,border_color = NA,
#                color = greenred(64),cluster_rows = F,treeheight_row = 0,cluster_cols = F,
#                annotation_col = annCol[colnames(plotdata),c("FGFR3","MeTIL.score")], annotation_colors = annColors,
#                annotation_row = data.frame(Gene = tmp$Simple_UCSC_RefGene_Name,row.names = rownames(tmp)),
#                show_rownames = T,show_colnames = F,legend = T)
# invisible(dev.off())

tmp <- read.delim(file.path(resDir.NOFILTER.C12vsN,"myDMP.FGFR3_curated.txt"),sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
tmp <- tmp[which((tmp$deltaBeta >=0.2|tmp$deltaBeta <=-0.2)& tmp$adj.P.Val < 0.05), ]
tmp <- tmp[which(tmp$gene %in% c("FGFR3","KRT15","KRT5","MIR200A","GATA3")),]
tmp <- tmp[order(tmp$gene),]
plotdata <- as.data.frame(myCombat[rownames(tmp),FGFR3.anno$SampleID])
annColors[["Gene"]] <- c("FGFR3" = sun,
                         "KRT15" = blue,
                         "KRT5" = yellow,
                         #"MIR200A" = brown,
                         "GATA3" = seagreen)
library(viridis)
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_FGFR3_cancer genes.pdf"), width = 8,height = 7)
hv <- pheatmap(plotdata,border_color = NA,
               color = inferno(64),cluster_rows = F,treeheight_row = 0,cluster_cols = F,
               annotation_col = annCol[colnames(plotdata),c("FGFR3","MeTIL.score")], annotation_colors = annColors,
               annotation_row = data.frame(Gene = tmp$gene,row.names = rownames(tmp)),
               show_rownames = T,show_colnames = F,legend = T)
invisible(dev.off())

# check DMP localisation
dmp.tmp <- myDMP.FGFR3_curated$Wild_to_Mutated
dmp.tmp <- cbind.data.frame(dmp.tmp,feature[rownames(dmp.tmp),])
write.table(dmp.tmp,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.FGFR3_curated_addFeatures.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

gainMethInFGFR3 <- rownames(dmp.tmp[which(dmp.tmp$Mutated_AVG >= 0.4 
                                    & dmp.tmp$Wild_AVG <= 0.2 
                                    & dmp.tmp$adj.P.Val < 0.05),])
tmp <- dmp.tmp[gainMethInFGFR3,21:31]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:31])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/length(gainMethInFGFR3)*100,0)
# FALSE. TRUE.           Names percentage
# 1      10     0         isShelf          0
# 2       6     4         isShore         40
# 3      10     0        isUTRCGI          0
# 4      10     0       isBodyCGI          0
# 5      10     0       isPromCGI          0
# 6      10     0           isCGI          0
# 7       9     1          isTFBS         10
# 8       9     1 isOpenChromatin         10
# 9       9     1           isDMR         10
# 10      9     1      isEnhancer         10
# 11      1     9         isDNase         90
tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage")
tmp3$Class <- rep(c("GainMethInFGFR3","EPIC"),each=11)
ggplot(tmp3,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_continuous(breaks=seq(0, 90, 10), limits=c(0, 90))
ggsave(file.path(resDir.NOFILTER.C12vsN,"pair-wise histogram for gainMethInFGFR3 and EPIC array.pdf"))

# enrich of enhancer
fisher.test(matrix(c(1,27485,9,838374),nrow = 2,byrow = T)) # 0.2757
# enrich of shore
fisher.test(matrix(c(4,154546,6,711313),nrow = 2,byrow = T)) # 0.08
# enrich of DNase
fisher.test(matrix(c(9,497025,1,368834),nrow = 2,byrow = T)) #0.05163
# enrich of DMR (yes)
fisher.test(matrix(c(1,39245,9,826614),nrow = 2,byrow = T)) #0.3711


lossMethInFGFR3 <- rownames(dmp.tmp[which(dmp.tmp$Mutated_AVG <= 0.2 
                                          & dmp.tmp$Wild_AVG >= 0.4 
                                          & dmp.tmp$adj.P.Val < 0.05),])
tmp <- dmp.tmp[lossMethInFGFR3,21:31]
freq <- table(col(tmp),as.matrix(tmp)) #get the counts of True and False level
Names <- c(colnames(dmp.res)[21:31])
tmp <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp$percentage <- round(as.numeric(tmp$TRUE.)/length(lossMethInFGFR3)*100,0)
# FALSE. TRUE.           Names percentage
# 1   16492  1195         isShelf          7
# 2   15834  1853         isShore         10
# 3   17657    30        isUTRCGI          0
# 4   17576   111       isBodyCGI          1
# 5   17643    44       isPromCGI          0
# 6   17424   263           isCGI          1
# 7   14892  2795          isTFBS         16
# 8   15821  1866 isOpenChromatin         11
# 9   17282   405           isDMR          2
# 10  17017   670      isEnhancer          4
# 11   6419 11268         isDNase         64
tmp2 <- feature
freq <- table(col(tmp2),as.matrix(tmp2)) #get the counts of True and False level
Names <- colnames(tmp2)
tmp2 <- data.frame(cbind(freq,Names),stringsAsFactors = F)
tmp2$percentage <- round(as.numeric(tmp2$TRUE.)/865859*100,0)

tmp3 <- rbind(tmp,tmp2); colnames(tmp3) <- c("False","True","Region","Percentage")
tmp3$Class <- rep(c("LossMethInFGFR3","EPIC"),each=11)
ggplot(tmp3,aes(x=Region,y=Percentage,fill=Class)) + geom_bar(stat="identity",width=1,position = position_dodge(0.7)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_continuous(breaks=seq(0, 70, 10), limits=c(0, 70))
ggsave(file.path(resDir.NOFILTER.C12vsN,"pair-wise histogram for lossMethInFGFR3 and EPIC array.pdf"))

# enrich of DNase (yes)
fisher.test(matrix(c(11268,497025,6419,368834),nrow = 2,byrow = T)) #<0.001
# enrich of Enhancer (yes)
fisher.test(matrix(c(670,27485,17017,838374),nrow = 2,byrow = T)) #<0.001
# enrich of TFBS
fisher.test(matrix(c(2795,134170,14892,731689),nrow = 2,byrow = T)) #0.266


# DMR GSEA pathways to plot
gaby.FGFR3.pathway <- read.table("F:/Project/UTUC_BLCA/InputData/Gaby_DMR_FGFR3_ssgsea_pathways_to_plot.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gaby.FGFR3.pathway.list <- list()
for (i in 1:nrow(gaby.FGFR3.pathway)) {
  lab <- rownames(gaby.FGFR3.pathway)[i]
  tmp <- gaby.FGFR3.pathway[i,"Genes"]
  tmp <- toupper(unlist(strsplit(tmp," ")))
  tmp <- rownames(annoPb[which(annoPb$Simple_UCSC_RefGene_Name %in% tmp),])
  gaby.FGFR3.pathway.list[[lab]] <- tmp
}

FGFR3.anno <- data.frame(FGFR3 = pheno.FGFR3_curated,SampleID = samples.FGFR3_curated,stringsAsFactors = F)
FGFR3.anno <- FGFR3.anno[order(FGFR3.anno$FGFR3),]
rownames(FGFR3.anno) <- FGFR3.anno$SampleID
gaby.FGFR3.enrichscore <- gsva(as.matrix(myCombat[,FGFR3.anno$SampleID]),gaby.FGFR3.pathway.list,method="gsva") # better than ssgsea
gaby.FGFR3.enrichscore.backup <- gaby.FGFR3.enrichscore
# gaby.FGFR3.enrichscore <- gaby.FGFR3.enrichscore.backup
p <- vector()
for (i in 1:nrow(gaby.FGFR3.enrichscore)) {
  tmp <- wilcox.test(as.numeric(gaby.FGFR3.enrichscore[i,1:15]),as.numeric(gaby.FGFR3.enrichscore[i,16:31]),alternative = "less")$p.value
  p <- append(p, tmp)
}
rownames(gaby.FGFR3.enrichscore) <- paste0(rownames(gaby.FGFR3.enrichscore),"_p_",p)

plotscore <- as.data.frame(na.omit(standarize.fun(indata = gaby.FGFR3.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotscore)), "euclidean"), "ward.D")
mycol <- colorpanel(64,low=cyan,mid = "black",high=peach)
annCol$FGFR3 <- FGFR3.anno[rownames(annCol),"FGFR3"]
annColors[["FGFR3"]] <- c("Wild" = lightgrey,"Mutated"=cherry)
pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_FGFR3_enrichmentscore_to_plot.pdf"), width = 15,height = 3.5)
hv <- pheatmap(as.matrix(plotscore),border_color = NA,
               color = mycol,cluster_rows = dendsort(hcg),treeheight_row = 0,cluster_cols = F,
               annotation_col = annCol[colnames(plotscore),c("FGFR3","MeTIL.score")], annotation_colors = annColors,
               show_rownames = T,show_colnames = F,legend = T)
invisible(dev.off())

############################################################
### DMP analysis between FGFR3 curated mutation and Normal ###
samples.FGFR3_curated_vs_Normal <- c(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "1" & UTUC.annotation$`DNA Methylation data` == "yes"),]),
                                     "06N","12N","02N","03N","07N","16N","09N","04N")
pheno.FGFR3_curated_vs_Normal <- rep(c("FGFR3_Mutated","Normal"),c(15,8))
myDMP.FGFR3_curated_vs_Normal <-   champ.DMP(beta=myCombat[,samples.FGFR3_curated_vs_Normal],
                                             pheno=pheno.FGFR3_curated_vs_Normal,
                                             arraytype = "EPIC",adjPVal = 1); 
write.table(myDMP.FGFR3_curated_vs_Normal$FGFR3_Mutated_to_Normal,file = file.path("H:/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C12vsN/ChAMP_RES_NOFILTER_C12vsN","myDMP.FGFR3_curated_vs_Normal.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
save(myDMP.FGFR3_curated_vs_Normal,file = file.path(resDir.NOFILTER.C12vsN,"myDMP.FGFR3_curated_vs_Normal.rda"))

###########################################################################
### supervised clustering of entire UTUC sample with pbForSpClusterTCGA ###
library(dendextend)
annCol.entire <- annCol[,1]
annCol.entire <- c(annCol.entire, rep("Normal",8))
annCol.entire <- data.frame(Group = annCol.entire,row.names = c(rownames(annCol),c("06N","12N","02N","03N","07N","16N","09N","04N")),stringsAsFactors = F)
annColors.entire <- list()
annColors.entire[["Group"]] = c("C1" = red,"C2" = blue, "Normal"=green)

tmp <- myCombat[gainMethInC1,rownames(annCol.entire)]
hcg <- hclust(distanceMatrix(as.matrix(t(tmp)), "euclidean"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(tmp), "euclidean"), "ward.D")

hcs.dend <- click_rotate(as.dendrogram(dendsort(hcs)),continue = TRUE)
# 
# group <- cutree(hcs,k = 7)
# annCol.entire$Cluster <- paste0("C",group[rownames(annCol.entire)])
# annColors.entire[["Cluster"]] = c("C1" = red,"C2" = blue, "C3"=green,"C4"=yellow,"C5"=lightred,"C6"=darkblue,"C7"="black")

pdf(file.path(resDir.NOFILTER.C12vsN, "UTUC_heatmap_SupervisedCluster_C1vsC2_DMPs_to_entirecohort.pdf"), width = 8,height = 8)
# pheatmap(as.matrix(tmp),
#          color = greenred(64),
#          cluster_cols = dendsort(hcs),
#          cluster_rows = hcs.dend,
#          annotation_col = annCol.entire,
#          annotation_colors = annColors.entire,
#          treeheight_col = 40,
#          treeheight_row = 0,
#          show_rownames = F)

hv <- aheatmap(as.matrix(tmp),
               color = c("#0074FE","#96EBF9","#FEE900","#F00003"),
               Rowv = dendsort(hcg),
               Colv = hcs.dend,
               annCol = annCol.entire, annColors = annColors.entire)
invisible(dev.off())

# MeTIL boxplot between EpiC-high and low
# MeTIL
tmp1 <- MeTIL.score[rownames(annCol.entire[which(annCol.entire$Group == "C1"),,drop = F])]
tmp2 <- MeTIL.score[rownames(annCol.entire[which(annCol.entire$Group == "C2"),,drop = F])]
p <- wilcox.test(tmp1,tmp2) # 0.002525
df <- data.frame("MeTIL"=c(tmp1,tmp2),"DMBcluster"=rep(c("EpiC-high","EpiC-low"),c(length(tmp1),length(tmp2))))
p2 <- ggplot(df,aes(x=DMBcluster,y=MeTIL,fill=DMBcluster)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c(red,blue)) +
  ggtitle("pValue = 0.003") +
  geom_jitter(aes(color = DMBcluster),alpha=1) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(resDir.NOFILTER.C12vsN,"pcaMeTIL between UTUC DMBcluster C1 and C2 .pdf"),width = 2,height = 4.5)

# save image
save.image(file.path(workdir,"chAMP.RData"))
