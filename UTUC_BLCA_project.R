tumorname <- "UTUC_BLCA"
workdir <- "H:/Jan2020/UTUC_2020/UTUC_BLCA"; setwd(workdir)
#workdir <- "D:/UTUC_BLCA"; setwd(workdir)
#workdir <- "H:/UTUC_BLCA"; setwd(workdir)
tumor.path <- workdir
res.path    <- file.path(tumor.path, "Results")
fig.path    <- file.path(tumor.path, "Figures")
data.path   <- file.path(tumor.path, "InputData")
if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }

GinfoFile <- "overlapTable.txt"  ## Gene annotation file created by Jack
CMGFile <- "CMG_Genes_list.txt"
comRFun.path <- file.path(tumor.path,"commonFun") ## universal functions
shareFun.path <- file.path(tumor.path,"SharedScripts") ## functions shared by this type of analysis only
script.path <- file.path(tumor.path,"Scripts")
comAnn.path <- file.path(tumor.path,"Annotation") ## path containning annotation files created by Jack
Ginfo <- read.table(file.path(comAnn.path,GinfoFile),sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
UTUC.annotation <- read.table(file.path(data.path,"UTUC_annotation.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
createGinfoFileFlag <- F ## to create Gene annotation file? FALSE by default
options(warn =0)
tailrows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique") ## tail rows not related to gene names in read counts input file
lowcut.mRNA <- 1
lowcut.lnc <- 0.25 # FPKM cutoff for low expressed genes
lowpct <- 0.9 # if a gene is below lowcut in 90% of samples, then it will be not be included for the analysis
DMBIIscore <- read.table(file = file.path(comAnn.path,"DNA-methylation-based-immune-infiltration-scores_XLu.txt"), header=T, row.names=NULL, sep="\t", quote="", stringsAsFactors=F, check.names=F)
DMBIIscore <- DMBIIscore[which(DMBIIscore$Tumorname=="BLCA"),]

# create annotation
# tmp1 <- read.table(file.path(data.path,"UTUC_annotation.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# tmp2 <- read.table(file.path(data.path,"UTUC_survival.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# tmp2$Age <- 2019 - as.numeric(sapply(strsplit(tmp2$`Date of birth`,"/"),"[",1)) + 1
# tmp2 <- tmp2[,setdiff(colnames(tmp2),"Date of birth")]
# tmp <- merge(tmp1,tmp2,by="SimpleID",all.x=T); rownames(tmp) <- tmp$SimpleID
# tmp <- tmp[rownames(UTUC.annotation),]
# write.table(tmp,file.path(data.path,"UTUC_annotation.txt"),sep = "\t",row.names = F)

source(file.path(shareFun.path,"LoadCommonFunLib.R")) 
source(file.path(script.path,"createList.R"))
library(survminer)
library(estimate)
library(dendsort)
library(dendextend)
library(pgirmess)
library(ggplot2)
library(cowplot)
library(enrichplot)
library(grid)
library(gridExtra)
library(easyGgplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(CMScaller) #containing ntp (Nearest Templates Prediction)
library(sva)
library(edgeR)
library(ComplexHeatmap)
library(reshape)
library(RColorBrewer)

kmplotfun <- function(group=NULL, Sinfo, fig.path, outFigFile=NULL, p.txt=NULL, xlim=NULL, lgdpos="topright", col=NULL, width=4, height=4, cex=1, line=1, at=NA, xunit="Day",legend.title=""){
  
  # line: controls the vertical position of the mtext
  # at: controls the horizontal position of the mtext
  # xunit: c("Day", "Year"): Unit of X axis
  
  if(!is.element(xunit, c("Day", "Year"))) { stop("xunit must be Day or Year!") }
  
  NAcols <- c("NotApplicable", "Not Applicable", "Not Available", "NotAvailable", "[Not Available]", "Discrepancy", "[Discrepancy]", "Unknown", "Completed","-","N/A")  
  
  common = intersect(names(group), row.names(Sinfo))
  group <- group[common]
  clin <- Sinfo[common,]  
  
  clin$time=as.character(clin$OS.time)
  #clin[clin$time %in% NAcols,"time"] <- as.character(clin[clin$time %in% NAcols, "days_to_last_followup"]) #error?
  clin <- clin[!(clin$time %in% NAcols), ]
  if( !setequal( unique(clin$vital_status), c("LIVING", "DECEASED") ) && !setequal( unique(clin$vital_status), c("Dead", "Alive") ) ) {stop("vital_status error!")}
  clin$time=as.numeric(clin$OS.time)
  
  for (j in 1: nrow(clin)){
    if (clin$vital_status[j]=="LIVING" || clin$vital_status[j]=="Alive") clin$death[j]<-0
    else clin$death[j]<-1
  }
  
  grplevels <- levels(factor(group))
  #if (is.null(col)) { col=rainbow(length(grplevels)) }
  #if (is.null(col)) { col=  c(darkred, lightred, blue) }
  if (is.null(col)) { col=  c(red,blue) }
  data = cbind(group[common], clin[common, c("time", "death")])
  data$group <- factor(data$group)
  names(data)=c("group", "time", "death")
  xlab="Days"
  if (xunit=="Year") {
    data$time=data$time/365
    xlab="Years"
  }
  
  fitd=survdiff(Surv(time, death)~ group, data=data, na.action=na.exclude)
  p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  
  fit=survfit(Surv(time, death)~ group, data=data, type="kaplan-meier",error="greenwood", conf.type="plain", na.action=na.exclude)
  summary <- summary(fit)$table[,-c(2,3)]
  
  HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  if(is.na(at)) {at=max(data$time)/2}
  
  if (!is.null(outFigFile)) {
    legend.txt = NULL
    for (k in 1:length(unique(group))) {
      legend.txt =c(legend.txt, paste(grplevels[k], ", N=", as.character(fitd$n[k]), sep=""))
    }
    
    survp <- ggsurvplot(fit, conf.int=F,
                        risk.table=F, risk.table.col="strata",
                        palette=col,data=data,
                        size=1.5,font.legend = 11,
                        xlim=c(0,3),break.time.by=1,risk.table.y.text = FALSE,
                        legend.title=legend.title,
                        legend.labs=legend.txt,
                        pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = "")),
                                     HR, CI, sep = "\n")
    )
    ggsave(file = file.path(fig.path,outFigFile),width = width,height = height)
  }
  
  return(list(summary=summary, pval=p.val,legend=legend.txt))
}

kmplotfun2 <- function(group=NULL, Sinfo, showHR = F, fig.path, outFigFile=NULL, p.txt=NULL, xlim=NULL, lgdpos="topright", col=NULL, width=4, height=4, cex=1, line=1, at=NA, xunit="Day",cutsurv = 3,legend.title=""){
  
  if(!is.element(xunit, c("Day", "Month"))) { stop("xunit must be Day or Month!") }
  
  common = intersect(names(group), row.names(Sinfo))
  group <- group[common]
  clin <- Sinfo[common,]  
  clin$time=as.numeric(clin$OS.time)
  clin$death=as.numeric(clin$OS)

  grplevels <- levels(factor(group))
  #if (is.null(col)) { col=rainbow(length(grplevels)) }
  #if (is.null(col)) { col=  c(darkred, lightred, blue) }
  if (is.null(col)) { col=  c(red,blue) }
  data = cbind(group[common], clin[common, c("time", "death")])
  data$group <- factor(data$group)
  names(data)=c("group", "time", "death")
  xlab="Days"
  if (xunit=="Month") {
    data$time=data$time/30.5
    xlab="Months"
  }
  
  fitd=survdiff(Surv(time, death)~ group, data=data, na.action=na.exclude)
  p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  
  fit=survfit(Surv(time, death)~ group, data=data, type="kaplan-meier",error="greenwood", conf.type="plain", na.action=na.exclude)
  summary <- summary(fit)$table[,-c(2,3)]
  
  HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  if(is.na(at)) {at=max(data$time)/2}
  
  if (!is.null(outFigFile)) {
    legend.txt = NULL
    for (k in 1:length(unique(group))) {
      legend.txt =c(legend.txt, paste(grplevels[k], ", N=", as.character(fitd$n[k]), sep=""))
    }
    
    if (isTRUE(showHR)) {
      survp <- ggsurvplot(fit, conf.int=F,
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = "")),
                                     HR, CI, sep = "\n"))
    } else {
      survp <- ggsurvplot(fit, conf.int=F,
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = ""))))
    }
    
    ggsave(file = file.path(fig.path,outFigFile),width = width,height = height)  }
  
  return(list(summary=summary, pval=p.val,legend=legend.txt))
}

kmplotfun.event1 <- function(group=NULL, Sinfo, showHR = F, fig.path, outFigFile=NULL, p.txt=NULL, xlim=NULL, lgdpos="topright", col=NULL, width=4, height=4, cex=1, line=1, at=NA, xunit="Day",cutsurv = 3,legend.title=""){
  
  if(!is.element(xunit, c("Day", "Month"))) { stop("xunit must be Day or Month!") }
  
  common = intersect(names(group), row.names(Sinfo))
  group <- group[common]
  clin <- Sinfo[common,]  
  clin$time=as.numeric(clin$PFS.time)
  clin$death=as.numeric(clin$PFS)
  
  grplevels <- levels(factor(group))
  #if (is.null(col)) { col=rainbow(length(grplevels)) }
  #if (is.null(col)) { col=  c(darkred, lightred, blue) }
  if (is.null(col)) { col=  c(red,blue) }
  data = cbind(group[common], clin[common, c("time", "death")])
  data$group <- factor(data$group)
  names(data)=c("group", "time", "death")
  xlab="Days"
  if (xunit=="Month") {
    data$time=data$time/30.5
    xlab="Months"
  }
  
  fitd=survdiff(Surv(time, death)~ group, data=data, na.action=na.exclude)
  p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  
  fit=survfit(Surv(time, death)~ group, data=data, type="kaplan-meier",error="greenwood", conf.type="plain", na.action=na.exclude)
  summary <- summary(fit)$table[,-c(2,3)]
  
  HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  if(is.na(at)) {at=max(data$time)/2}
  
  if (!is.null(outFigFile)) {
    legend.txt = NULL
    for (k in 1:length(unique(group))) {
      legend.txt =c(legend.txt, paste(grplevels[k], ", N=", as.character(fitd$n[k]), sep=""))
    }
    
    if (isTRUE(showHR)) {
      survp <- ggsurvplot(fit, conf.int=F,#fun = "event",
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = "")),
                                     HR, CI, sep = "\n"))
    } else {
      survp <- ggsurvplot(fit, conf.int=F,#fun = "event",
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = ""))))
    }
    
    ggsave(file = file.path(fig.path,outFigFile),width = width,height = height)  }
  
  return(list(summary=summary, pval=p.val,legend=legend.txt))
}

kmplotfun.event2 <- function(group=NULL, Sinfo, showHR = F, fig.path, outFigFile=NULL, p.txt=NULL, xlim=NULL, lgdpos="topright", col=NULL, width=4, height=4, cex=1, line=1, at=NA, xunit="Day",cutsurv = 3,legend.title=""){
  
  if(!is.element(xunit, c("Day", "Month"))) { stop("xunit must be Day or Month!") }
  
  common = intersect(names(group), row.names(Sinfo))
  group <- group[common]
  clin <- Sinfo[common,]  
  clin$time=as.numeric(clin$PFS.time)
  clin$death=as.numeric(clin$PFS2)
  
  grplevels <- levels(factor(group))
  #if (is.null(col)) { col=rainbow(length(grplevels)) }
  #if (is.null(col)) { col=  c(darkred, lightred, blue) }
  if (is.null(col)) { col=  c(red,blue) }
  data = cbind(group[common], clin[common, c("time", "death")])
  data$group <- factor(data$group)
  names(data)=c("group", "time", "death")
  xlab="Days"
  if (xunit=="Month") {
    data$time=data$time/30.5
    xlab="Months"
  }
  
  fitd=survdiff(Surv(time, death)~ group, data=data, na.action=na.exclude)
  p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  
  fit=survfit(Surv(time, death)~ group, data=data, type="kaplan-meier",error="greenwood", conf.type="plain", na.action=na.exclude)
  summary <- summary(fit)$table[,-c(2,3)]
  
  HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  if(is.na(at)) {at=max(data$time)/2}
  
  if (!is.null(outFigFile)) {
    legend.txt = NULL
    for (k in 1:length(unique(group))) {
      legend.txt =c(legend.txt, paste(grplevels[k], ", N=", as.character(fitd$n[k]), sep=""))
    }
    
    if (isTRUE(showHR)) {
      survp <- ggsurvplot(fit, conf.int=F,#fun = "event",
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = "")),
                                     HR, CI, sep = "\n"))
    } else {
      survp <- ggsurvplot(fit, conf.int=F,#fun = "event",
                          risk.table=F, risk.table.col="strata",
                          palette=col,data=data,
                          size=0.8,font.legend = 10,
                          xlim=c(0,cutsurv),break.time.by=10,risk.table.y.text = FALSE,
                          legend.title=legend.title,
                          legend.labs=legend.txt,
                          pval=paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                   paste("p = ",round(p.val,3), sep = ""))))
    }
    
    ggsave(file = file.path(fig.path,outFigFile),width = width,height = height)  }
  
  return(list(summary=summary, pval=p.val,legend=legend.txt))
}

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
darkblue <- "#1d00ff"
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
#################################################################
## UTUC DNA Methylation analysis plan                 ###########
# comfused samples were removed : "01N","01T","07T","11N","11T" #
# All analysis should remove sex chromosome probes ##############
#################################################################
# save image : DNA2

#load and preprocess data
methFile <- file.path(data.path,"1493_Malouf-beta_values-noob.csv")
annoFile <- file.path(comAnn.path,"MethylationEPIC_v-1-0_B4.csv")
SinfoFile<- file.path(data.path,"Sinfo.txt")

UTUCsinfo <- read.table(SinfoFile,sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
annoRaw <- read.csv(annoFile,check.names = F,stringsAsFactors = F,header = T,row.names = 1)
anno <- annoRaw[,c(1,15:19,20,21,23,22,33,35,36,37,38)]

tmp <- sapply(strsplit(anno$UCSC_RefGene_Name,";"),"[",1)
anno$Simple_UCSC_RefGene_Name <- tmp
tmp <- sapply(strsplit(anno$UCSC_RefGene_Group,";"),"[",1)
anno$Simple_UCSC_RefGene_Group <- tmp
anno[which(anno$Simple_UCSC_RefGene_Name == "1-Dec"),"Simple_UCSC_RefGene_Name"] <- "DEC1"
anno[which(anno$Simple_UCSC_RefGene_Name == "1-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH1"
anno[which(anno$Simple_UCSC_RefGene_Name == "1-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT1"
anno[which(anno$Simple_UCSC_RefGene_Name == "10-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH10" #wrong with gene name, should be modified in next project
anno[which(anno$Simple_UCSC_RefGene_Name == "11-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH11"
anno[which(anno$Simple_UCSC_RefGene_Name == "11-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT11"
anno[which(anno$Simple_UCSC_RefGene_Name == "13-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT13"
anno[which(anno$Simple_UCSC_RefGene_Name == "14-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT14"
anno[which(anno$Simple_UCSC_RefGene_Name == "2-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH2"
anno[which(anno$Simple_UCSC_RefGene_Name == "3-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH3"
anno[which(anno$Simple_UCSC_RefGene_Name == "4-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH4"
anno[which(anno$Simple_UCSC_RefGene_Name == "5-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH5"
anno[which(anno$Simple_UCSC_RefGene_Name == "6-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH6"
anno[which(anno$Simple_UCSC_RefGene_Name == "7-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH7"
anno[which(anno$Simple_UCSC_RefGene_Name == "8-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH8"
anno[which(anno$Simple_UCSC_RefGene_Name == "9-Mar"),"Simple_UCSC_RefGene_Name"] <- "MARTCH9"
anno[which(anno$Simple_UCSC_RefGene_Name == "5-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT5"
anno[which(anno$Simple_UCSC_RefGene_Name == "8-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT8"
anno[which(anno$Simple_UCSC_RefGene_Name == "9-Sep"),"Simple_UCSC_RefGene_Name"] <- "SEPT9"
write.table(anno,file.path(comAnn.path,"Methylation_probes_annotation_simplified.txt"),sep = "\t",row.names = T,col.names = NA)

orgmeth <- read.csv(methFile,check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
colnames(orgmeth) <- orgmeth["Sample ID",]
orgmeth <- orgmeth[-c(1:3),]; colnames(orgmeth)[1] <- "01N" ;tmp <- rownames(orgmeth) #866091
orgmeth <- as.data.frame(sapply(orgmeth, as.numeric)); rownames(orgmeth) <- tmp
orgmeth <- orgmeth[,setdiff(colnames(orgmeth),c("01N","01T","07T","11N","11T"))]
orgmeth <- as.data.frame(na.omit(orgmeth)) #855447
##########################################
### New DNA methylation plan 6/25/2018 ###
##########################################
#1.	Please generate two files of DNA methylation data. One contains CpG island probes only, and another contains all of other probes (non-CpG island probes). For both files, please remove all probes in chromosome X and Y.
#2.	Please remove samples 01T 7T 11T and 01N 11N
#3.	For clustering on all samples (tumor plus normal), please don't perform filtering of probes
#4.	For clustering on tumor samples only, please perform filtering of probes as specified before.
probe.noXY <- intersect(rownames(orgmeth),rownames(annoRaw[which(!(annoRaw$CHR %in% c("X","Y"))),]))
orgmeth <- orgmeth[probe.noXY,] #836691
probe.CpG <- intersect(probe.noXY,rownames(anno[which(anno$Relation_to_UCSC_CpG_Island == "Island"),]))
probe.noCpG <- intersect(probe.noXY,rownames(anno[which(anno$Relation_to_UCSC_CpG_Island != "Island"),]))

meth.CpG <- orgmeth[probe.CpG,]
meth.noCpG <- orgmeth[probe.noCpG,]

write.table(meth.CpG,file.path(res.path,"meth.CpG.noXY.txt"),sep = "\t",row.names = T,col.names = NA)
write.table(meth.noCpG,file.path(res.path,"meth.noCpG.noXY.txt"),sep = "\t",row.names = T,col.names = NA)

ind_ord1 <- order(apply(meth.CpG, 1, sd),decreasing=T)
ind_ord2 <- order(apply(meth.noCpG, 1, sd),decreasing=T)

#Transform to numeric!
k=5000
sel_idx = ind_ord1[1:k]
tmp=meth.CpG[sel_idx,]
meth.CpG.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.CpG.sel) <- rownames(tmp)

meth.CpG.sel.withName <- meth.CpG.sel
meth.CpG.sel.withName$GeneSymbol <- anno[rownames(meth.CpG.sel.withName),"Simple_UCSC_RefGene_Name"]; meth.CpG.sel.withName <- meth.CpG.sel.withName[,c(44,1:43)]
write.table(meth.CpG.sel.withName,file.path(res.path,"meth.CpG.noXY.top5000.with.GeneSymbol.nofilter.txt"),sep = "\t",row.names = T,col.names = NA)
rm(meth.CpG.sel.withName)

k=5000
sel_idx = ind_ord2[1:k]
tmp=meth.noCpG[sel_idx,]
meth.noCpG.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.noCpG.sel) <- rownames(tmp)

meth.noCpG.sel.withName <- meth.noCpG.sel
meth.noCpG.sel.withName$GeneSymbol <- anno[rownames(meth.noCpG.sel.withName),"Simple_UCSC_RefGene_Name"]; meth.noCpG.sel.withName <- meth.noCpG.sel.withName[,c(44,1:43)]
write.table(meth.noCpG.sel.withName,file.path(res.path,"meth.noCpG.noXY.top5000.with.GeneSymbol.nofilter.txt"),sep = "\t",row.names = T,col.names = NA)
rm(meth.noCpG.sel.withName)

#	Perform clustering analysis on all samples (tumor plus normal) by using most variable Probes in CpG islands without probe filtering
hcs <- hclust(distanceMatrix(as.matrix(meth.CpG.sel), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.CpG.sel)), "pearson"), "ward.D")
group.CpG.TvsN <- cutree(hcs,k=4)
group.CpG.TvsN <- paste0("cluster",group.CpG.TvsN); names(group.CpG.TvsN) <- colnames(meth.CpG.sel)
pdf(file.path(fig.path,"meth.CpG.noXY.nofilter.top5000.pdf"))
hv <- aheatmap(as.matrix(meth.CpG.sel), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
               annCol=data.frame("cluster"=group.CpG.TvsN), annColors=list("cluster"=c("cluster1"=red,"cluster2"=blue,"cluster3"=yellow,"cluster4"=green)), 
               color=c("blue","green","yellow","red"), revC=TRUE, fontsize=5, labRow = NA, cexCol = 1,cexAnn = 0.6)
invisible(dev.off())

# Perform clustering analysis on all samples (tumor plus normal) by using most variable Probes outside CpG islands without probe filtering
hcs <- hclust(distanceMatrix(as.matrix(meth.noCpG.sel), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.noCpG.sel)), "pearson"), "ward.D")
group.noCpG.TvsN <- cutree(hcs,k=4)
group.noCpG.TvsN <- paste0("cluster",group.noCpG.TvsN); names(group.noCpG.TvsN) <- colnames(meth.noCpG.sel)
pdf(file.path(fig.path,"meth.noCpG.noXY.nofilter.top5000.pdf"))
hv <- aheatmap(as.matrix(meth.noCpG.sel), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
               annCol=data.frame("cluster"=group.noCpG.TvsN), annColors=list("cluster"=c("cluster1"=red,"cluster2"=blue,"cluster3"=yellow,"cluster4"=green)), 
               color=c("blue","green","yellow","red"), revC=TRUE, fontsize=5, labRow = NA, cexCol = 1,cexAnn = 0.6)
invisible(dev.off())

meth.CpG.sel.scale <- as.data.frame(t(scale(t(meth.CpG.sel))))
meth.noCpG.sel.scale <- as.data.frame(t(scale(t(meth.noCpG.sel))))

batchPCA(indata = meth.CpG.sel.scale,batch = colnames(meth.CpG.sel.scale),fig.dir = fig.path,PCA.fig.title = "meth.CpG.noXY.nofilter.top5000.AllSamples.PCA",cols = rep(blue,43),showID = T,cex = 0.7,showLegend = F)
batchPCA(indata = meth.noCpG.sel.scale,batch = colnames(meth.noCpG.sel.scale),fig.dir = fig.path,PCA.fig.title = "meth.noCpG.noXY.nofilter.top5000.AllSamples.PCA",cols = rep(blue,43),showID = T,cex = 0.7,showLegend = F)

#For clustering on tumor samples only, please perform filtering of probes as specified before.
#filter probes
UTUC.meth.tumor.samples <- rownames(UTUCsinfo[which(UTUCsinfo$Type=="Tumor"),]); UTUC.meth.tumor.samples <- setdiff(UTUC.meth.tumor.samples,c("01T","07T","11T"))
UTUC.meth.normal.samples <- rownames(UTUCsinfo[which(UTUCsinfo$Type=="Normal"),]); UTUC.meth.normal.samples <- setdiff(UTUC.meth.normal.samples,c("01N","11N"))

meth <- rmHighMethInNorm(indata = orgmeth,samInfo = UTUCsinfo,medcut = 0.2,betacut = 0.2,highpct = 0.5)
meth.CpG.filtered <- as.data.frame(na.omit(meth[intersect(probe.CpG,rownames(meth)),]))
meth.noCpG.filtered <- as.data.frame(na.omit(meth[intersect(probe.noCpG,rownames(meth)),]))

#Attention: use tumor samples only here before calculate high variable probes!!!
ind_ord3 <- order(apply(meth.CpG.filtered[,UTUC.meth.tumor.samples], 1, sd),decreasing=T)
ind_ord4 <- order(apply(meth.noCpG.filtered[,UTUC.meth.tumor.samples], 1, sd),decreasing=T)

k=5000
sel_idx = ind_ord3[1:k]
tmp=meth.CpG.filtered[sel_idx,]
meth.CpG.filtered.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.CpG.filtered.sel) <- rownames(tmp)

k=5000
sel_idx = ind_ord4[1:k]
tmp=meth.noCpG.filtered[sel_idx,]
meth.noCpG.filtered.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.noCpG.filtered.sel) <- rownames(tmp)

hcs <- hclust(distanceMatrix(as.matrix(meth.CpG.filtered.sel[,UTUC.meth.tumor.samples]), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.CpG.filtered.sel[,UTUC.meth.tumor.samples])), "pearson"), "ward.D")
group.CpG.Tumor <- cutree(hcs,k=2)
group.CpG.Tumor <- paste0("cluster",group.CpG.Tumor); names(group.CpG.Tumor) <- UTUC.meth.tumor.samples
pdf(file.path(fig.path,"meth.CpG.noXY.filter.top5000.pdf"))
hv <- aheatmap(as.matrix(meth.CpG.filtered.sel[,UTUC.meth.tumor.samples]), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
               annCol=data.frame("cluster"=group.CpG.Tumor), annColors=list("cluster"=c("cluster1"=red,"cluster2"=blue)), 
               color=c("blue","green","yellow","red"), revC=TRUE, fontsize=5, labRow = NA, cexCol = 1,cexAnn = 0.6)
invisible(dev.off())

hcs <- hclust(distanceMatrix(as.matrix(meth.noCpG.filtered.sel[,UTUC.meth.tumor.samples]), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.noCpG.filtered.sel[,UTUC.meth.tumor.samples])), "pearson"), "ward.D")
group.noCpG.Tumor <- cutree(hcs,k=2)
group.noCpG.Tumor <- paste0("cluster",group.noCpG.Tumor); names(group.noCpG.Tumor) <- UTUC.meth.tumor.samples
pdf(file.path(fig.path,"meth.noCpG.noXY.filter.top5000.pdf"))
hv <- aheatmap(as.matrix(meth.noCpG.filtered.sel[,UTUC.meth.tumor.samples]), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
               annCol=data.frame("cluster"=group.noCpG.Tumor), annColors=list("cluster"=c("cluster1"=red,"cluster2"=blue)), 
               color=c("blue","green","yellow","red"), revC=TRUE, fontsize=5, labRow = NA, cexCol = 1,cexAnn = 0.6)
invisible(dev.off())

meth.CpG.sel.scale <- as.data.frame(t(scale(t(meth.CpG.filtered.sel[,UTUC.meth.tumor.samples]))))
meth.noCpG.sel.scale <- as.data.frame(t(scale(t(meth.noCpG.filtered.sel[,UTUC.meth.tumor.samples]))))

batchPCA(indata = meth.CpG.sel.scale,batch = colnames(meth.CpG.sel.scale),fig.dir = fig.path,PCA.fig.title = "meth.CpG.noXY.filter.top5000.TumorSamples.PCA",cols = rep(blue,35),showID = T,cex = 0.7,showLegend = F)
batchPCA(indata = meth.noCpG.sel.scale,batch = colnames(meth.noCpG.sel.scale),fig.dir = fig.path,PCA.fig.title = "meth.noCpG.noXY.filter.top5000.TumorSamples.PCA",cols = rep(blue,35),showID = T,cex = 0.7,showLegend = F)

##############################################
### Analysis plan 6/28/2018 ##################
##############################################
probe.cpg.promoter <- intersect(probe.CpG,intersect(rownames(anno[which(anno$Simple_UCSC_RefGene_Group %in% c("TSS1500","TSS200")),]),rownames(orgmeth))) #60483
probe.cpg.genebody <- intersect(probe.CpG,intersect(rownames(anno[which(anno$Simple_UCSC_RefGene_Group == "Body"),]),rownames(orgmeth))) #37410

#1.	Identify hypermethylated genes with at least 10% (at least 3 samples) samples with beta value>=0.3 after filtering out the probes methylated in normal samples
Gmeth <- apply(orgmeth, 2, function(x) tapply(x, INDEX=factor(anno[rownames(orgmeth), "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
Gmeth <- as.data.frame(na.omit(Gmeth))
write.table(Gmeth,file = file.path(data.path,"UTUC_medianAvg_Gmeth_noXY.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

methSinfo <- read.table(file.path(data.path,"methSinfo.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
UTUC.FPKM <- read.table(file.path(res.path,"UTUC.FPKM.all.features.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)

combat.UTUC.FPKM = read.table(file.path(res.path,"UTUC.lowcut1.lowpct0.9.combatFPKM.mRNA.LowExpfiltered.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
colnames(combat.UTUC.FPKM) <- c("11T","13T","19T","01T","21T","22T","25T","26T","29T","02T","31T","04T","08T","09T","10T","16T","37T","03T","05T","07T")
# in tumor
samples <- intersect( intersect( rownames(UTUCsinfo[which(UTUCsinfo$Type=="Tumor"),]), colnames(combat.UTUC.FPKM) ), colnames(Gmeth) )
MethPct <- CalMethPct(gdata=Gmeth, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)

# in normal
normsamples <- intersect(colnames(Gmeth), rownames(UTUCsinfo[which(UTUCsinfo$Type=="Normal"),]))
NormMethPct <- CalMethPct(gdata=Gmeth, samples=normsamples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)

if(all(names(MethPct)==names(NormMethPct))) {
  write.table(data.frame(HighMethPctInTumor=MethPct, HighMethPctInNormal=NormMethPct), file=file.path(res.path, "high_methylation_pct_of_genes_noXY_UTUC.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
}else{
  stop("MethPct and NormMethPct names mismatch")
}

#2.	Identify epigenetically silenced genes for promoter CpG probes only after filtering out the probes methylated in normal samples.

samples <- intersect( intersect( rownames(UTUCsinfo[which(UTUCsinfo$Type=="Tumor"),]), colnames(combat.UTUC.FPKM) ), colnames(Gmeth) )
tmp <- orgmeth[probe.cpg.promoter,]
Gmeth.cpg.promoter <- apply(tmp, 2, function(x) tapply(x, INDEX=factor(anno[rownames(tmp), "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
Gmeth.cpg.promoter <- as.data.frame(na.omit(Gmeth.cpg.promoter))
write.table(Gmeth.cpg.promoter,file = file.path(data.path,"UTUC_medianAvg_noXY_Gmeth.cpg.promoter.txt"),row.names=T, col.names=NA, sep="\t", quote=F)
test <- CalExpDownPct(gdata=Gmeth.cpg.promoter, expdata=combat.UTUC.FPKM, outfile=file.path(res.path, "Frequency_methylated_and_underexpressed_of_genes_UTUC.txt"), samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=methSinfo, highMethCut=0.3, Ginfo=Ginfo, pcut=0.05)

#3.	Identify differentially methylated probes with beta value mean difference (not folder change) for CpG island-based two clustering on tumor samples only with the following column:
#a.	Promoter status
#b.	Enhancer status
#c.	DMR status
#d.	OpenChromatin
#e.	TFBS status
tmp <- as.data.frame(na.omit(meth))
tmp$GeneSymbol <- toupper(anno[rownames(tmp),"Simple_UCSC_RefGene_Name"])
tmp <- tmp[,c(ncol(tmp),1:(ncol(tmp)-1))]
group.CpG.Tumor2 <- group.CpG.Tumor
group.CpG.Tumor2[which(group.CpG.Tumor2=="cluster1")] <- "CpGTumorHigh"  
group.CpG.Tumor2[which(group.CpG.Tumor2=="cluster2")] <- "CpGTumorLow"  
methcomplist <- createList.CpG.Tumor(group=group.CpG.Tumor2)
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

#4.	Identify differentially methylated probes with beta value mean difference (not folder change) for non-CpG island-based two clustering on tumor samples only with the following column:
group.noCpG.Tumor2 <- group.noCpG.Tumor
group.noCpG.Tumor2[which(group.noCpG.Tumor2=="cluster1")] <- "noCpGTumorHigh"  
group.noCpG.Tumor2[which(group.noCpG.Tumor2=="cluster2")] <- "noCpGTumorLow"  
methcomplist <- createList.noCpG.Tumor(group=group.noCpG.Tumor2)
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

#5.	For CpG-based clustering with all samples, Identify differentially methylated probes with beta value mean difference (not folder change) between normal sample clustering and one of three tumor clustering only with the following column:
group.CpG.TvsN2 <- group.CpG.TvsN
group.CpG.TvsN2[which(group.CpG.TvsN2=="cluster1")] <- "CpGNormal"
group.CpG.TvsN2[which(group.CpG.TvsN2=="cluster2")] <- "CpGTumorHigh"
group.CpG.TvsN2[which(group.CpG.TvsN2=="cluster3")] <- "CpGTumorMid"
group.CpG.TvsN2[which(group.CpG.TvsN2=="cluster4")] <- "CpGTumorLow"
methcomplist <- createList.CpG.THvsN(group=group.CpG.TvsN2[which(group.CpG.TvsN2 %in% c("CpGTumorHigh","CpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

methcomplist <- createList.CpG.TMvsN(group=group.CpG.TvsN2[which(group.CpG.TvsN2 %in% c("CpGTumorMid","CpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

methcomplist <- createList.CpG.TLvsN(group=group.CpG.TvsN2[which(group.CpG.TvsN2 %in% c("CpGTumorLow","CpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

#6.	For non-CpG-based clustering with all samples, Identify differentially methylated probes with beta value mean difference (not folder change) between normal sample clustering and one of three tumor clustering only with the following column:
group.noCpG.TvsN2 <- group.noCpG.TvsN
group.noCpG.TvsN2[which(group.noCpG.TvsN2=="cluster1")] <- "noCpGNormal"
group.noCpG.TvsN2[which(group.noCpG.TvsN2=="cluster2")] <- "noCpGTumorMid"
group.noCpG.TvsN2[which(group.noCpG.TvsN2=="cluster3")] <- "noCpGTumorHigh"
group.noCpG.TvsN2[which(group.noCpG.TvsN2=="cluster4")] <- "noCpGTumorLow"
methcomplist <- createList.noCpG.THvsN(group=group.noCpG.TvsN2[which(group.noCpG.TvsN2 %in% c("noCpGTumorHigh","noCpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

methcomplist <- createList.noCpG.TMvsN(group=group.noCpG.TvsN2[which(group.noCpG.TvsN2 %in% c("noCpGTumorMid","noCpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)

methcomplist <- createList.noCpG.TLvsN(group=group.noCpG.TvsN2[which(group.noCpG.TvsN2 %in% c("noCpGTumorLow","noCpGNormal"))])
methDEtest(indata=tmp, res.path=res.path, complist=methcomplist, featType="Diff methylation", overwt=TRUE)


a <- read.table("F:/Project/UTUC_BLCA/Annotation/Methylation_probes_annotation_simplified.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
b <- read.table("C:/Users/XLu6/Desktop/UTUC_BLCA/Results/Diff methylation_wilcox_test_result.CpGTumorHigh_vs_CpGTumorLow.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
c<- data.frame("isShelf"=rep("FALSE",nrow(a)),"isShore"=rep("FALSE",nrow(a)),"isUTRCGI"=rep("FALSE",nrow(a)),"isBodyCGI"=rep("FALSE",nrow(a)),"isPromCGI"=rep("FALSE",nrow(a)),"isCGI"=rep("FALSE",nrow(a)),"isTFBS"=rep("FALSE",nrow(a)),"isOpenChromatin"=rep("FALSE",nrow(a)),"isDMR"=rep("FALSE",nrow(a)),"isEnhancer"=rep("FALSE",nrow(a)),"isDNase"=rep("FALSE",nrow(a)),stringsAsFactors = F,row.names = rownames(a))

for (i in 1:nrow(a)) {
  # if(a[i,"Phantom4_Enhancers"]!="") {
  #   c[i,"isEnhancer"] <- "TRUE"
  # }
  if(a[i,"Phantom5_Enhancers"]!="") {
    c[i,"isEnhancer"] <- "TRUE"
  }
  # if(isTRUE(a[i,"450k_Enhancer"])) {
  #   c[i,"isEnhancer"] <- "TRUE"
  # }
  if((a[i,"Simple_UCSC_RefGene_Group"] %in% c("TSS1500","TSS200")) & (a[i,"Relation_to_UCSC_CpG_Island"]=="Island")) {
    c[i,"isPromCGI"] <- "TRUE"
  }
  if((a[i,"Simple_UCSC_RefGene_Group"] %in% c("3'UTR","5'UTR")) & (a[i,"Relation_to_UCSC_CpG_Island"]=="Island")) {
    c[i,"isUTRCGI"] <- "TRUE"
  }
  if((a[i,"Simple_UCSC_RefGene_Group"] %in% "Body") & (a[i,"Relation_to_UCSC_CpG_Island"]=="Island")) {
    c[i,"isBodyCGI"] <- "TRUE"
  }
  if(a[i,"Relation_to_UCSC_CpG_Island"] %in% c("N_Shelf","S_Shelf")) {
    c[i,"isShelf"] <- "TRUE"
  }
  if(a[i,"Relation_to_UCSC_CpG_Island"] %in% c("N_Shore","S_Shore")) {
    c[i,"isShore"] <- "TRUE"
  }
  if(is.na(a[i,"OpenChromatin_Evidence_Count"])==F) {
    c[i,"isOpenChromatin"] <- "TRUE"
  }
  if(is.na(a[i,"TFBS_Evidence_Count"])==F) {
    c[i,"isTFBS"] <- "TRUE"
  }
  if(a[i,"DMR"]!="") {
    c[i,"isDMR"] <- "TRUE"
  }
  if(a[i,"Relation_to_UCSC_CpG_Island"]=="Island") {
    c[i,"isCGI"] <- "TRUE"
  }
  if(a[i,"DNase_Hypersensitivity_NAME"]!="") {
    c[i,"isDNase"] <- "TRUE"
  }
  display.progress(index = i,totalN = nrow(a))
}
# remove rs and control
c <- c[1:865859,]
write.table(c,file.path(comAnn.path,"Methylation_probes_features_logical2.txt"),sep = "\t",row.names = T,col.names = NA)

feature <- read.table(file.path(comAnn.path,"Methylation_probes_features_logical2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
methDEfile <- "Diff methylation_wilcox_test_result.*.txt"
methDEfile <- dir(res.path,pattern = methDEfile)
for (i in 1:length(methDEfile)) {
  tmp <- read.table(file.path(res.path,methDEfile[i]),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
  tmp <- cbind(tmp,feature[rownames(tmp),])
  write.table(tmp,file.path(res.path,paste0(substr(methDEfile[i],start = 1,stop = nchar(methDEfile[i])-4),"_addFeature.txt")),sep = "\t",row.names = T,col.names = NA)
}

##############################################
### Analysis plan 7/16/2018 ##################
##############################################
# 1. Perform unsupervised clustering using all probes (cpg and nocpg) without XY of UTUC using allsamples with and without normal samples
# Tumor vs Normal: no filtering high probe in normal samples
ind_ord5 <- order(apply(orgmeth, 1, sd),decreasing=T)
k <- 0.01
sel_idx = ind_ord5[1:(length(ind_ord5)*k)]
tmp=orgmeth[sel_idx,]
meth.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.sel) <- rownames(tmp)

hcs <- hclust(distanceMatrix(as.matrix(meth.sel), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.sel)), "pearson"), "ward.D")
group.TvsN <- cutree(hcs,k=3)
group.TvsN <- paste0("Subgroup",group.TvsN); names(group.TvsN) <- colnames(meth.sel)

tmp <- data.frame("Cluster"=group.TvsN)
tmp <- cbind.data.frame(tmp,UTUC.annotation[rownames(tmp),c("Type","SEX","GRADE","pathological stage","FGFR3_mut")])
#tmp <- cbind.data.frame(tmp,UTUC.annotation[rownames(tmp),c("Type","FGFR3","SWI_SNF_noACTL6B")])
colnames(tmp) <- c("Cluster","Type","Gender","Grade","Pathological stage","FGFR3")
tmp$Tissue <- ifelse(is.na(tmp$Type),"Normal","Tumour")
tmp[is.na(tmp)] <- "N/A"
ClusterColor = c("#B31E22","#529B40","#020105","#383D8E")

#pdf(file.path(fig.path,paste0("test.meth.noXY.nofilter.top",k*100,"pct_TumorNormalSamples.pdf")))
pdf(file.path(fig.path,paste0("meth2.noXY.nofilter.top",k*100,"pct_TumorNormalSamples_FGFR3_curated.pdf")),height = 9,width = 9)
#pdf(file.path("H:/UTUC_BLCA/Figures",paste0("meth.noXY.nofilter.top",k*100,"pct_TumorNormalSamples_FGFR3_raw.pdf")),height = 5.5,width = 6)
# hv <- aheatmap(as.matrix(meth.sel), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
#                annCol=tmp[c(5,2,3,4,1)], annColors=list("Cluster"=c("Subgroup1"=ClusterColor[1],"Subgroup2"=ClusterColor[2],"Subgroup3"=nake),#,"Subgroup4"=ClusterColor[4]),
#                                                         "Type"=c("muscle-invasive"=soil,"Non-muscle invasive"=grey,"N/A"="white"),
#                                                         "FGFR3"=c("1"=purple,"0"=lightgrey,"N/A"="white"),
#                                                         "SWI_SNF_noACTL6B"=c("1"=purple,"0"=lightgrey,"N/A"="white"),
#                                                         "Tissue"=c("Tumour"="black","Normal"=nake)), 
#                color=c("#0074FE","#96EBF9","#FEE900","#F00003"),
#                revC=TRUE, fontsize=7, labRow = NA, cexCol = 0.9)

pheatmap(as.matrix(meth.sel),
         cluster_rows = hcg,
         cluster_cols = as.hclust(dendsort(as.dendrogram(hcs))),
         annotation_col = tmp[,c(1,6,5,4,3,2,7)],
         annotation_colors = list("Cluster"=c("Subgroup1"=ClusterColor[1],"Subgroup2"=ClusterColor[2],"Subgroup3"=nake),#,"Subgroup4"=ClusterColor[4]),
                                  "Type"=c("muscle-invasive"=soil,"Non-muscle invasive"=grey,"N/A"="white"),
                                  "FGFR3"=c("1"=purple,"0"=lightgrey,"N/A"="white"),                                 
                                  "Tissue"=c("Tumour"="black","Normal"=nake),
                                  "Gender" = c("M"=blue,"F"=sun,"N/A"="white"),
                                  "Grade" = c("Low"=nake,"High"="black","N/A"="white"),
                                  "Pathological stage" = c("pTa"=white,"pT1"=yellow,"pT2"=green,"pT3"=sun,"pT4"=cherry,"N/A"="white")),
         color = NMF:::ccRamp(c("#0074FE","#96EBF9","#FEE900","#F00003"),n=64),
         show_rownames = F,show_colnames = T,
         treeheight_row = 0,
)
invisible(dev.off())

# test FGFR3 independency
table(tmp$Cluster,tmp$FGFR3)
fisher.test(matrix(c(5,14,9,0),byrow = T,ncol = 2)) # 0.0006
# test FGFR3 curated independency
table(tmp$Cluster,tmp$FGFR3)
fisher.test(matrix(c(7,15,9,0),byrow = T,ncol = 2)) # 0.0008
# test MI independency
table(tmp$Cluster,tmp$Type)
fisher.test(matrix(c(10,14,7,4),byrow = T,ncol = 2)) # 0.289

# group.TvsN2 <- group.TvsN
# group.TvsN2[which(group.TvsN2=="Subgroup1")] <- "Normal"
# group.TvsN2[which(group.TvsN2 %in% c("Subgroup2","Subgroup3"))] <- "Tumor"
# write.csv(data.frame("Sample_Name"=names(group.TvsN2),"Sample_Group"=group.TvsN2),file.path(res.path,"Sinfo of Tumor Normal for RnBeads.csv"),row.names = F)

group.TvsN3 <- group.TvsN
group.TvsN3[which(group.TvsN3=="Subgroup1")] <- "Normal"
group.TvsN3[which(group.TvsN3=="Subgroup2")] <- "FGFR3_Wild"
group.TvsN3[which(group.TvsN3=="Subgroup3")] <- "FGFR3_Mutated"

write.csv(data.frame("Sample_Name"=names(group.TvsN3),"Sample_Group"=group.TvsN3),file.path(res.path,"Sinfo of detailed Tumor Normal for RnBeads.csv"),row.names = F)

# Tumor C1 VS C2: filter high probe in normal samples
#Attention: use tumor samples only here before calculate high variable probes!!!
# ind_ord6 <- order(apply(meth[,UTUC.meth.tumor.samples], 1, sd),decreasing=T)
# k <- 0.05
# sel_idx = ind_ord6[1:(length(ind_ord6)*k)]
# tmp=meth[sel_idx,]
# meth.sel <- as.data.frame(sapply(tmp, as.numeric))
# rownames(meth.sel) <- rownames(tmp)
# 
# hcs <- hclust(distanceMatrix(as.matrix(meth.sel[,UTUC.meth.tumor.samples]), "euclidean"), "ward.D")
# hcg <- hclust(distanceMatrix(as.matrix(t(meth.sel[,UTUC.meth.tumor.samples])), "pearson"), "ward.D")
# group.Tumor <- cutree(hcs,k=2)
# group.Tumor <- paste0("cluster",group.Tumor); names(group.Tumor) <- UTUC.meth.tumor.samples
# pdf(file.path(fig.path,paste0("meth.noXY.filter.top",k*100,"pct_TumorSamples.pdf")))
# hv <- aheatmap(as.matrix(meth.sel[,UTUC.meth.tumor.samples]), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
#                annCol=data.frame("cluster"=group.Tumor), annColors=list("cluster"=c("cluster1"=red,"cluster2"=blue)), 
#                color=c("blue","green","yellow","red"), revC=TRUE, fontsize=5, labRow = NA, cexCol = 1,cexAnn = 0.6)
# invisible(dev.off())
# 
# group.Tumor2 <- group.Tumor
# group.Tumor2[which(group.Tumor2=="cluster1")] <- "C1"
# group.Tumor2[which(group.Tumor2=="cluster2")] <- "C2"
# write.csv(data.frame("Sample_Name"=names(group.Tumor2),"Sample_Group"=group.Tumor2),file.path(res.path,"Sinfo of Tumor C1 C2 for RnBeads.csv"),row.names = F)

# Tumor C1 VS C2: nofilter high probe in normal samples
ind_ord7 <- order(apply(orgmeth[,UTUC.meth.tumor.samples], 1, sd),decreasing=T)
k <- 0.01
sel_idx = ind_ord7[1:(length(ind_ord7)*k)]
tmp=orgmeth[sel_idx,]
meth.sel <- as.data.frame(sapply(tmp, as.numeric))
rownames(meth.sel) <- rownames(tmp)

hcs <- hclust(distanceMatrix(as.matrix(meth.sel[,UTUC.meth.tumor.samples]), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(meth.sel[,UTUC.meth.tumor.samples])), "pearson"), "ward.D")
group.Tumor <- cutree(hcs,k=2)
group.Tumor <- paste0("cluster",group.Tumor); names(group.Tumor) <- UTUC.meth.tumor.samples
hcs.DMBcluster <- hcs
#save(hcs.DMBcluster,file = file.path(res.path,"hcs.DMBcluster.rda"))

annCol.C1vsC2 <- data.frame("cluster"=group.Tumor)
annCol.C1vsC2 <- cbind(annCol.C1vsC2,UTUC.annotation[rownames(annCol.C1vsC2),c("SEX","Type","GRADE","pathological stage","Lymph nodes","Localisation","FGFR3","SWI_SNF_noACTL6B")])
MeTIL.UTUC <- read.table("F:/Project/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C12vsN/ChAMP_RES_NOFILTER_C12vsN/MeTIL in DMBcluster C1 and C2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
annCol.C1vsC2$MeTIL <- MeTIL.UTUC[rownames(annCol.C1vsC2),"MeTIL.score"]
colnames(annCol.C1vsC2)[4] <- "Grade"
annCol.C1vsC2$cluster <- ifelse(annCol.C1vsC2$cluster == "cluster1","C1","C2")

annColors.C1vsC2 <- list()
annColors.C1vsC2[["cluster"]] <- c("C1"=red,"C2"=blue)
annColors.C1vsC2[["Type"]] <- c("muscle-invasive"=soil,"Non-muscle invasive"=grey)
annColors.C1vsC2[["SEX"]] <- c("M"=blue,"F"=sun)
annColors.C1vsC2[["Grade"]] <- c("Low"=nake,"High"="black")
annColors.C1vsC2[["pathological stage"]] <- c("pTa"=white,"pT1"=yellow,"pT2"=green,"pT3"=sun,"pT4"=cherry)
annColors.C1vsC2[["Lymph  nodes"]] <- c("N0"=lightred,"N1"=darkred,"N2"=cherry)
annColors.C1vsC2[["Localisation"]] <- c("renal pelvis"=red,"ureter"=blue,"N/A"="white")
annColors.C1vsC2[["FGFR3"]] <- c("1"=purple,"0"=lightgrey,"N/A"="white")
annColors.C1vsC2[["SWI_SNF_noACTL6B"]] <- c("1"=purple,"0"=lightgrey,"N/A"="white")
annColors.C1vsC2[["MeTIL"]] <- bluered(64)

pdf(file.path(fig.path,paste0("meth.noXY.nofilter.top",k*100,"pct_TumorSamples.pdf")),height = 5.5,width = 6)
#pdf(file.path("H:/UTUC_BLCA/Figures",paste0("meth.noXY.nofilter.top",k*100,"pct_TumorSamples.pdf")),height = 5.5,width = 6)
hv <- aheatmap(as.matrix(meth.sel[,UTUC.meth.tumor.samples]), 
               Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), 
               annCol=annCol.C1vsC2[,c(3,4,5,8,9,10,1)], annColors=annColors.C1vsC2, 
               color=greenred(64), revC=TRUE, fontsize=7, labRow = NA, cexCol = 1,cexAnn = 1)
invisible(dev.off())
save(meth.sel,file = file.path(res.path,"probes_For_UTUC_DMBcluster.rda"))
# MeTIL
tmp1 <- annCol.C1vsC2[which(annCol.C1vsC2$cluster == "C1"),"MeTIL"]
tmp2 <- annCol.C1vsC2[which(annCol.C1vsC2$cluster == "C2"),"MeTIL"]
p <- wilcox.test(tmp1,tmp2) #0.002525

# Type
fisher.test(table(annCol.C1vsC2$Type,annCol.C1vsC2$cluster)) #0.000945
# cluster1 cluster2
# muscle-invasive           16        1
# Non-muscle invasive        7       11

# Grade
fisher.test(table(annCol.C1vsC2$GRADE,annCol.C1vsC2$cluster)) #0.1061
# cluster1 cluster2
# High       22        9
# Low         1        3

# pathological stage
fisher.test(table(annCol.C1vsC2$`pathological stage`,annCol.C1vsC2$cluster)) #0.005439
# cluster1 cluster2
# pT1        2        4
# pT2        3        1
# pT3       10        0
# pT4        3        0
# pTa        5        7

fisher.test(matrix(c(10,12,13,0),byrow = T,ncol = 2)) #0.0008562
# cluster1 cluster2
# pTa12       10        12
# pT34        13        0

# lymph nodes
fisher.test(table(annCol.C1vsC2$`Lymph nodes`,annCol.C1vsC2$cluster)) #0.6908
# cluster1 cluster2
# N0       20       12
# N1        1        0
# N2        2        0

# localistation
fisher.test(table(annCol.C1vsC2$Localisation,annCol.C1vsC2$cluster)) #0.09532
# cluster1 cluster2
# N/A                 3        0
# renal pelvis       17        7
# ureter              3        5

group.Tumor3 <- group.Tumor
group.Tumor3[which(group.Tumor3=="cluster1")] <- "C1"
group.Tumor3[which(group.Tumor3=="cluster2")] <- "C2"
write.csv(data.frame("Sample_Name"=names(group.Tumor3),"Sample_Group"=group.Tumor3),file.path(res.path,"Sinfo of Tumor C1 C2 derived from nofiltering for ChAMP.csv"),row.names = F)


###########################################################################################################
# supervised clustering for TCGA 450 methylation data based on strict C1 vs C2 DMP 
# succeed
orgmeth.tcga <- read.table(file.path(data.path,"TCGA_BLCA_HumanMethylation450.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- colnames(orgmeth.tcga)
tmp <- paste0("BLCA",substr(tmp,start = 8,stop = 15))
colnames(orgmeth.tcga) <- tmp
normsamples.tcga <- tmp[which(sapply(strsplit(tmp,"-"),"[",3) == "11")]
tumorsamples.tcga <- tmp[which(sapply(strsplit(tmp,"-"),"[",3) == "01")]

load("F:/Project/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C12vsN/ChAMP_RES_NOFILTER_C12vsN/pbForSpClusterTCGA_14209.rda")
#mapping to probes for supervised clustering
pb <- intersect(rownames(orgmeth.tcga),pbForSpClusterTCGA)
methForSpCluster.tcga <- as.data.frame(na.omit(orgmeth.tcga[pb,tumorsamples.tcga])) #3790 412

indata <- methForSpCluster.tcga
N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))
map.res.path <- file.path(res.path, paste("pb_strict_diff0.2_fdr0.05_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
featType <- "UTUC_C1vsC2_DMP"
options(warn=-1)
sp.cluster2.tcga.ans <- plot.common.cluster(indata, 
                                            tumorname="BLCA", 
                                            N.cluster=N.cluster, 
                                            N.bootstrap=N.bootstrap, 
                                            N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                            N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                            map.res.path=map.res.path, fig.path=fig.path, 
                                            featType=featType,annCol = annCol.C1vsC2.tcga[colnames(indata),],annColors = annColors.C1vsC2.tcga,
                                            seed=123456, dist0="manhattan", dist="pearson",link0 = "ward.D2",link = "ward.D2",
                                            clstCol=c(red,blue), 
                                            namecvt=NULL, height = 25, fontsize=5, labRow = F, labCol = F,cexAnn = 0.5)

save(sp.cluster2.tcga.ans,file = file.path(res.path,"sp.cluster2.tcga.ans.rda"))
#use OS.time
tmp <- sp.cluster2.tcga.ans$group
Sinfo.tcga$DMPsCluster <- tmp[rownames(Sinfo.tcga)]
#write.table(Sinfo.tcga,file.path(data.path,"BLCA_ClinicalInfo.final.txt"),sep = "\t",row.names = T,col.names = NA)
#names(tmp) <- substr(names(tmp),start = 1,stop = 9)
#Sinfo.tcga <- read.table(file.path(data.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,fill = T)
#Sinfo.tcga <- Sinfo.tcga[which(Sinfo.tcga$type == "BLCA"),]
#rownames(Sinfo.tcga) <- Sinfo.tcga$SampleName
tmp <- ifelse(tmp == "cluster1","C1","C2"); names(tmp) <- names(sp.cluster2.tcga.ans$group)
tmp <- factor(tmp,levels = c("C2","C1"))
kmplotfun(group = tmp,col = c("#5bc0eb","#f25f5c"),Sinfo = Sinfo.tcga,fig.path = fig.path,outFigFile="KM2 plot supervised cluster of TCGA BLCA.pdf", xunit="Year",legend.title = "DMBcluster")

# $`summary`
# records events   *rmean *se(rmean)   median  0.95LCL  0.95UCL
# factor(group)=cluster1     219    112 4.870167  0.3980611 2.353425 1.876712 3.745205
# factor(group)=cluster2     192     68 5.756249  0.5404646 5.400000 2.476712 8.720548
# 
# $pval
# [1] "0.035"

# HR = 0.7234818

#heatmap for methForSpCluster.tcga
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")

#TIDE need work?
TIDE.res <- read.csv(file.path(res.path,"BLCA_TIDE_output2.csv"),header = T,row.names = 1,check.names = F,stringsAsFactors = F)
rownames(TIDE.res) <- substr(rownames(TIDE.res),1,12); TIDE.res$samples <- rownames(TIDE.res)
TIDE.res <- TIDE.res[,c("samples","Responder")]
tmp <- annCol.C1vsC2.tcga; tmp$samples <- rownames(tmp)
tmp <- merge(tmp,TIDE.res,by="samples",all.x=T)
rownames(tmp) <- tmp$samples; tmp <- tmp[,-1]
tmp$Responder <- ifelse(is.na(tmp$Responder),"N/A",ifelse(tmp$Responder == "True","True","False"))
tmp$SupervisedCluster <- sp.cluster2.tcga.ans$group[rownames(tmp)]
tmp$SupervisedCluster <- ifelse(tmp$SupervisedCluster == "cluster1","C1","C2")
annColors.C1vsC2.tcga[["Responder"]] <- c("True"=purple,"False"="white","N/A"=lightgrey)
annColors.C1vsC2.tcga[["SupervisedCluster"]] <- c("C1"=red,"C2"=blue)
annColors.C1vsC2.tcga[["FGFR3"]] <- c("1"="black","0"="white")
annColors.C1vsC2.tcga[["swisnf"]] <- c("1"="black","0"="white")
annColors.C1vsC2.tcga[["FGFR3_expr"]] <- bluered(64)
annColors.C1vsC2.tcga[["T cells"]] <- annColors.C1vsC2.tcga[["CD8 T cells"]] <- greenred(16)
annColors.C1vsC2.tcga[["Cytotoxic lymphocytes"]] <- annColors.C1vsC2.tcga[["NK cells"]] <-annColors.C1vsC2.tcga[["B lineage"]] <- greenred(16)
annColors.C1vsC2.tcga[["Monocytic lineage"]] <- annColors.C1vsC2.tcga[["Myeloid dendritic cells"]] <- greenred(16)
annColors.C1vsC2.tcga[["Neutrophils"]] <- annColors.C1vsC2.tcga[["Endothelial cells"]] <- annColors.C1vsC2.tcga[["Fibroblasts"]] <- greenred(16)

spheat.anno <- tmp
sig.mut <- c("FGFR3","KDM6A","KMT2D","KMT2C","ZFP36L1","ARID1A","STAG2","TP53","CRIPAK","GANAB")
tmp <- read.table(file.path(data.path,"BLCA_UTUC_sig_mut10.tsv"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
colnames(tmp) <- paste0("BLCA",substr(colnames(tmp),8,15))
tmp[is.na(tmp)] <- ""
tmp[tmp != ""] <- 1
tmp[tmp == ""] <- 0
tmp <- cbind.data.frame(t(tmp[names(sp.cluster2.tcga.ans$group)]),group = sp.cluster2.tcga.ans$group)
blca.mut10 <- tmp[,1:10]
spheat.anno$FGFR3 <- blca.mut10[rownames(spheat.anno),"FGFR3"]

tcga.mut.swisnf <- read.table(file.path(data.path,"BLCA_SWI_SNF20.tsv"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
colnames(tcga.mut.swisnf) <- paste0("BLCA",substr(colnames(tcga.mut.swisnf),8,15))
tcga.mut.swisnf[is.na(tcga.mut.swisnf)] <- ""
tcga.mut.swisnf[tcga.mut.swisnf != ""] <- 1
tcga.mut.swisnf[tcga.mut.swisnf == ""] <- 0
tmp <- rownames(tcga.mut.swisnf)
tcga.mut.swisnf <- as.data.frame(sapply(tcga.mut.swisnf,as.numeric)); rownames(tcga.mut.swisnf) <- tmp
tcga.mut.swisnf <- as.data.frame(t(tcga.mut.swisnf))
tcga.mut.swisnf$swisnf <- ifelse(rowSums(tcga.mut.swisnf) >= 1,1,0)
tcga.mut.swisnf <- cbind.data.frame(tcga.mut.swisnf[names(sp.cluster2.tcga.ans$group),],group = sp.cluster2.tcga.ans$group)
spheat.anno$swisnf <- tcga.mut.swisnf[rownames(spheat.anno),"swisnf"]

FGFR3.expr <- read.table(file.path(res.path,"FPKM.BLCA.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(FGFR3.expr) <- substr(colnames(FGFR3.expr),1,12)
FGFR3.expr <- data.frame(FGFR3 = scale(as.numeric(log2(FGFR3.expr["FGFR3",]+1))),
                         row.names = colnames(FGFR3.expr))
tmp.add <- data.frame(FGFR3 = c(0,0,0,0,0),
                      row.names = setdiff(rownames(spheat.anno),rownames(FGFR3.expr)))
FGFR3.expr <- rbind.data.frame(FGFR3.expr,tmp.add)
FGFR3.expr <- FGFR3.expr[rownames(spheat.anno),]

BLCA.mcpcounter <- read.table(file.path(res.path,"MCPscore.BLCA.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
colnames(BLCA.mcpcounter) <- substr(colnames(BLCA.mcpcounter),1,12)
BLCA.mcpcounter <- standarize.fun(log2(BLCA.mcpcounter+1),halfwidth = 2)
tmp.add <- matrix(0,ncol = 5,nrow = 10,dimnames = list(rownames(BLCA.mcpcounter),setdiff(rownames(spheat.anno),colnames(BLCA.mcpcounter))))
BLCA.mcpcounter <- cbind.data.frame(BLCA.mcpcounter,tmp.add)
BLCA.mcpcounter <- BLCA.mcpcounter[,rownames(spheat.anno)]

spheat.anno$FGFR3_expr <- FGFR3.expr
spheat.anno <- cbind.data.frame(spheat.anno,as.data.frame(t(BLCA.mcpcounter)))

mycol <- colorRampPalette(brewer.pal(11,"Spectral"))(11)[11:1]
pdf(file = file.path(fig.path,"BLCA2_heatmap_SupervisedCluster_UTUC_C1vsC2_DMPs_annotation_legend.pdf"),height = 35,width = 8)
# pdf(file.path("H:/UTUC_BLCA/Figures","BLCA2_heatmap_SupervisedCluster_UTUC_C1vsC2_DMPs.pdf"),height = 32,width = 24)
# aheatmap(as.matrix(indata), 
#          Rowv=dendsort(as.dendrogram(hcg)), 
#          Colv=sp.cluster2.tcga.ans$dendro, 
#          annCol=tmp[colnames(indata),c(3,5,8,9,10,11, # clinic
#                                        21,22,27,30,23,28,26,24,25,29,31, # mut
#                                        32,34,33,35,37,36,38,39,47,45,40,41,44,43,42,46,
#                                        49,14:19,52)], # cnv
#          annColors = annColors.C1vsC2.tcga, 
#          #color=c("blue","blue","green","yellow","red","red"), 
#          color=c("#0074FE","#96EBF9","#FEE900","#F00003"),
#          revC=TRUE, fontsize=7.2,labRow = NA,labCol = NA,cexAnn = 1)

indata2 <- indata[1:2,]
pheatmap(as.matrix(indata2), 
         #cluster_rows = hcg, 
         cluster_cols = as.hclust(sp.cluster2.tcga.ans$dendro), 
         annotation_col = spheat.anno[colnames(indata),c(3,5,8,9,10,11, # clinic
                                       38,53,55,54,56:65,
                                       49,16:19,52)], # cnv
         annotation_colors = annColors.C1vsC2.tcga, 
         treeheight_col = 4,
         #color=c("blue","blue","green","yellow","red","red"), 
         color=NMF:::ccRamp(c("#0074FE","#96EBF9","#FEE900","#F00003"),n=64),
         show_colnames = F,show_rownames = F,
         fontsize_col = 6)
invisible(dev.off())

pdf(file = file.path(fig.path,"BLCA2_heatmap_SupervisedCluster_UTUC_C1vsC2_DMPs_body.pdf"),height = 6,width = 8)
pheatmap(as.matrix(indata), 
         #cluster_rows = hcg, 
         cluster_cols = as.hclust(sp.cluster2.tcga.ans$dendro), 
         # annotation_col = spheat.anno[colnames(indata),c(3,5,8,9,10,11, # clinic
         #                                                 38,53,55,54,56:65,
         #                                                 49,16:19,52)], # cnv
         annotation_colors = annColors.C1vsC2.tcga, 
         treeheight_col = 0,
         treeheight_row = 0,
         #color=c("blue","blue","green","yellow","red","red"), 
         color=NMF:::ccRamp(c("#0074FE","#96EBF9","#FEE900","#F00003"),n=64),
         show_colnames = F,show_rownames = F,
         fontsize_col = 6)
invisible(dev.off())
# DMI DNA methylation based immune infiltration score (MBII score) 
tmp1 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C1"),"DMI scores"]
tmp2 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C2"),"DMI scores"]
p <- wilcox.test(tmp1,tmp2) # < 2.2e-16
df <- data.frame("MBII"=c(tmp1,tmp2),"DMBcluster"=rep(c("C1","C2"),c(length(tmp1),length(tmp2))))
p1 <- ggplot(df,aes(x=DMBcluster,y=MBII,fill=DMBcluster)) + geom_boxplot(alpha=0.6) + scale_fill_manual(values = c(red,blue)) + ggtitle("pValue < 2.2e-16") + geom_jitter(aes(color = DMBcluster),alpha=1) + stat_summary(fun.y=mean, geom="point", shape=18, size=4) + theme(legend.position="none")

# MeTIL
tmp1 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C1"),"pca.MeTIL"]
tmp2 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C2"),"pca.MeTIL"]
p <- wilcox.test(tmp1,tmp2) # < 2.2e-16
df <- data.frame("MeTIL"=c(tmp1,tmp2),"DMBcluster"=rep(c("C1","C2"),c(length(tmp1),length(tmp2))))
# p2 <- ggplot(df,aes(x=DMBcluster,y=MeTIL,fill=DMBcluster)) + 
#   geom_boxplot(alpha=0.6) + 
#   scale_fill_manual(values = c(red,blue)) + 
#   ggtitle("pValue < 2.2e-16") + 
#   geom_jitter(aes(color = DMBcluster),alpha=1) + 
#   stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")

ggplot(data = df,aes(x = DMBcluster, y = MeTIL, fill = DMBcluster))+
  scale_fill_manual(values = c(red,blue)) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8,color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+ # 边框线黑色
  theme_pubr() +
  ylab("MeTIL") +
  xlab("Supervised DMBcluster")  +
  annotate(geom="text", cex=6,
           x=1.5, y=0.11, 
           label="P < 2.2e-16", # 添加P值
           color="black") + 
  rremove("legend.title")+ # 移除图例
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = "none")+
  font("xylab",size=15)+  
  font("xy",size=15)+ 
  font("xy.text", size = 15) +  
  font("legend.text",size = 15)
ggsave(file.path(fig.path,"pcaMeTIL between TCGA supervised DMBcluster C1 and C2 .pdf"),width = 3,height = 4.5)

# Neutrophils
tmp1 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C1"),"Neutrophils"]
tmp2 <- spheat.anno[which(spheat.anno$SupervisedCluster == "C2"),"Neutrophils"]
p <- wilcox.test(tmp1,tmp2) # 1.436e-08
df <- data.frame("Neutrophils"=c(tmp1,tmp2),"DMBcluster"=rep(c("C1","C2"),c(length(tmp1),length(tmp2))))
ggplot(data = df,aes(x = DMBcluster, y = Neutrophils, fill = DMBcluster))+
  scale_fill_manual(values = c("#2874C5","#EABF00")) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8,color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+ # 边框线黑色
  theme_pubr() +
  ylab("Neutrophils") +
  xlab("Supervised DMBcluster")  +
  annotate(geom="text", cex=6,
           x=1.5, y=2, 
           label="P = 1.4e-8", # 添加P值
           color="black") + 
  rremove("legend.title")+ # 移除图例
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = "none")+
  font("xylab",size=15)+  
  font("xy",size=15)+ 
  font("xy.text", size = 15) +  
  font("legend.text",size = 15)
ggsave(file.path(fig.path,"Neutrophils between TCGA supervised DMBcluster C1 and C2 .pdf"),width = 3,height = 4.5)

# p <- plot_grid(p1, p2, labels = c("a", "b"))
# ggsave(file.path(fig.path,"DMI and MeTIL between TCGA supervised DMBcluster C1 and C2 .pdf"),width = 6,height = 5)
# ggsave("H:/UTUC_BLCA/Figures/DMI and MeTIL between TCGA supervised DMBcluster C1 and C2 .pdf",width = 6,height = 5)
# ggsave("H:/UTUC_BLCA/Figures/DMI and pcaMeTIL between TCGA supervised DMBcluster C1 and C2 .pdf",width = 6,height = 5)

# check if TCGA FGFR3 and MLL3 mutation is associated with METIL score
tmp1 <- tmp[which(tmp$`mutation in FGFR3` == "yes"),"pca.MeTIL"]
tmp2 <- tmp[which(tmp$`mutation in FGFR3` == "no"),"pca.MeTIL"]
wilcox.test(tmp1,tmp2) # 9.039e-6
df <- data.frame("MeTIL"=c(tmp1,tmp2),"FGFR3_Mutation"=rep(c("Mutated","Wild"),c(length(tmp1),length(tmp2))))
p2 <- ggplot(df,aes(x=FGFR3_Mutation,y=MeTIL,fill=FGFR3_Mutation)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(purple,grey)) + 
  ggtitle("pValue = 9.039e-6") + 
  geom_jitter(aes(color = FGFR3_Mutation),alpha=1) +   scale_color_manual(values = c(purple,grey)) + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave("H:/UTUC_BLCA/Figures/pcaMeTIL between TCGA FGFR3 Mutation .pdf",width = 4,height = 4.5)

tmp1 <- tmp[which(tmp$`mutation in FGFR3` == "yes"),"DMI scores"]
tmp2 <- tmp[which(tmp$`mutation in FGFR3` == "no"),"DMI scores"]
wilcox.test(tmp1,tmp2) # 2.828e-08

tcga_MLL3 <- read.table("H:/UTUC_BLCA/InputData/cBioPortal_mutation_FGFR3_KMT2C.txt",row.names = 1,sep = "\t",check.names = F,header = T,stringsAsFactors = F)
colnames(tcga_MLL3) <- paste0("BLCA",substr(colnames(tcga_MLL3),8,15))
tcga_MLL3[tcga_MLL3 == ""] <- "no"
tcga_MLL3[tcga_MLL3 != "no"] <- "yes"
tcga_MLL3[is.na(tcga_MLL3)] <- "no"

tmp$`mutation in MLL3` <- NA
tmp[intersect(rownames(tmp),colnames(tcga_MLL3)),"mutation in MLL3"] <- as.character(tcga_MLL3["KMT2C",intersect(rownames(tmp),colnames(tcga_MLL3))])

tmp1 <- tmp[which(tmp$`mutation in MLL3` == "yes"),"pca.MeTIL"]
tmp2 <- tmp[which(tmp$`mutation in MLL3` == "no"),"pca.MeTIL"]
wilcox.test(tmp1,tmp2) # 0.8146

tmp1 <- tmp[which(tmp$`mutation in MLL3` == "yes"),"DMI scores"]
tmp2 <- tmp[which(tmp$`mutation in MLL3` == "no"),"DMI scores"]
wilcox.test(tmp1,tmp2) # 0.441

# fisher test for mutation and cnv in supervised BLCA cluster
# order by high mutated in cluster1 (top)
# order by high mutated in cluster2 (bottom)

# TP53
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in TP53`)) #0.04817 top1
# no yes
# cluster1 104 116
# cluster2 110  82

# RB1
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in RB1`)) #0.002755 top2
# no yes
# cluster1 170  50
# cluster2 170  22

# KDM6A
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in KDM6A`)) #0.02455 bottom1
# no yes
# cluster1 173  47
# cluster2 132  60

# ELF3
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in ELF3`)) #0.003913 bottom4
# no yes
# cluster1 203  17
# cluster2 159  33

# CDKN1A
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in CDKN1A`)) #0.02423 bottom5
# no yes
# cluster1 207  13
# cluster2 168  24

# STAG2
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in STAG2`)) #0.02153 bottom3
# no yes
# cluster1 198  22
# cluster2 157  35

# KMT2A
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in KMT2A`)) #0.001494 top3
# no yes
# cluster1 186  34
# cluster2 181  11

# FGFR3
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in FGFR3`)) #1.019e-08 bottom2
# no yes
# cluster1 209  11
# cluster2 145  47

# PSIP1
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in PSIP1`)) #0.006362 bottom6
# no yes
# cluster1 215   5
# cluster2 176  16

# ERBB2
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in ERBB2`)) #0.03371 top4
# no yes
# cluster1 186  34
# cluster2 176  16

# TMCO4
fisher.test(table(tmp$SupervisedCluster,tmp$`mutation in TMCO4`)) #0.01533 bottom7
# no yes
# cluster1 218   2
# cluster2 182  10

tmp2 <- tmp[,c(34:47,51)]
tmp2 <- tmp[-which(tmp2$`focal amplification in PPARG` == "NaN"),]

# amplification
# YWHAZ
fisher.test(table(ifelse(tmp2$`focal amplification in YWHAZ` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.01102 No.1
# cluster1 cluster2
# 0      182      138
# 1       36       52

# CCND1
fisher.test(table(ifelse(tmp2$`focal amplification in CCND1` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.007735 No.3
# cluster1 cluster2
# 0      202      160
# 1       16       30

# PPARG
fisher.test(table(ifelse(tmp2$`focal amplification in PPARG` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.002942 No.2
# cluster1 cluster2
# 0      200      155
# 1       18       35

# BCL2L1
fisher.test(table(ifelse(tmp2$`focal amplification in BCL2L1` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.004557 No.4
# cluster1 cluster2
# 0      205      162
# 1       13       28

# ERBB2
fisher.test(table(ifelse(tmp2$`focal amplification in ERBB2` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.006321 No.6
# cluster1 cluster2
# 0      213      174
# 1        5       16

# AHR
fisher.test(table(ifelse(tmp2$`focal amplification in AHR` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.01907 No.5
# cluster1 cluster2
# 0      209      170
# 1        9       20

# FGFR3
fisher.test(table(ifelse(tmp2$`focal amplification in FGFR3` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.0004213 No.7
# cluster1 cluster2
# 0      218      180
# 1        0       10

# deletion
# CDKN2A
fisher.test(table(ifelse(tmp2$`focal deletion in CDKN2A` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.0006324 No.1
# cluster1 cluster2
# 0      142       91
# 1       76       99

# PDE4D
fisher.test(table(ifelse(tmp2$`focal deletion in PDE4D` == ".","0","1"),
                  tmp2$SupervisedCluster)) #3.184e-05 No.4
# cluster1 cluster2
# 0      186      129
# 1       32       61

# CCSER1
fisher.test(table(ifelse(tmp2$`focal deletion in CCSER1` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.01352 No.5
# cluster1 cluster2
# 0      189      146
# 1       29       44

# CREBBP
fisher.test(table(ifelse(tmp2$`focal deletion in CREBBP` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.00261 No.8
# cluster1 cluster2
# 0      198      152
# 1       20       38

# WWOX
fisher.test(table(ifelse(tmp2$`focal deletion in WWOX` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.0005644 No.7
# cluster1 cluster2
# 0      197      148
# 1       21       42

# PTEN
fisher.test(table(ifelse(tmp2$`focal deletion in PTEN` == ".","0","1"),
                  tmp2$SupervisedCluster)) #1.963e-06 No.6
# cluster1 cluster2
# 0      205      148
# 1       13       42

# NCOR1
fisher.test(table(ifelse(tmp2$`focal deletion in NCOR1` == ".","0","1"),
                   tmp2$SupervisedCluster)) #2.887e-07 No.3
# cluster1 cluster2
# 0      190      125
# 1       28       65

# RAD51B
fisher.test(table(ifelse(tmp2$`focal deletion in RAD51B` == ".","0","1"),
                  tmp2$SupervisedCluster)) #0.0005886 No.9
# cluster1 cluster2
# 0      204      157
# 1       14       33

# PTPRD
fisher.test(table(ifelse(tmp2$`focal deletion in PTPRD` == ".","0","1"),
                  tmp2$SupervisedCluster)) #1.32e-07 No.2
# cluster1 cluster2
# 0      182      114
# 1       36       76

# supervised clustering by using strict dmp.enhancer
#mapping to probes for supervised clustering
# pb <- intersect(rownames(orgmeth.tcga),pbForSpCluster_strict_diff0.2_fdr0.05)
# pb <- intersect(pb,rownames(feature[which(feature$isEnhancer == "TRUE"),]))
# enhancerForSpCluster.tcga <- as.data.frame(na.omit(orgmeth.tcga[pb,tumorsamples.tcga])) #1526 412
# save(enhancerForSpCluster.tcga,file = file.path(data.path,"enhancerForSpCluster.tcga.rda"))
# 
# indata <- enhancerForSpCluster.tcga
# N.cluster <- 2
# N.bootstrap <- 500
# N.gene.per.bootstrap <- round(0.8*nrow(indata))
# N.sample.per.bootstrap <- round(0.8*ncol(indata))
# map.res.path <- file.path(res.path, paste("enhancer_strict_diff0.2_fdr0.05_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
# featType <- "UTUC_C1vsC2_DMP.enhancer"
# 
# options(warn=-1)
# sp.cluster2.enhancer.tcga.ans <- plot.common.cluster(indata, 
#                                             tumorname="BLCA", 
#                                             N.cluster=N.cluster, 
#                                             N.bootstrap=N.bootstrap, 
#                                             N.gene.per.bootstrap=N.gene.per.bootstrap, 
#                                             N.sample.per.bootstrap=N.sample.per.bootstrap, 
#                                             map.res.path=map.res.path, fig.path=fig.path, 
#                                             featType=featType, 
#                                             annCol = annCol.C1vsC2.tcga[colnames(indata),],annColors = annColors.C1vsC2.tcga,
#                                             seed=123456, dist0="manhattan", dist="pearson",link0 = "ward.D2",link = "ward.D2",
#                                             clstCol=c(red,blue), 
#                                             namecvt=NULL, height = 22, fontsize=5, labRow = F, labCol = F,cexAnn = 0.5)

# tmp <- sp.cluster2.enhancer.tcga.ans$group
# tmp <- substr(names(tmp),start = 1,stop = 9)
# tmp1 <- rownames(Sinfo.tcga)
# tmp2 <- intersect(tmp,tmp1)
# index <- match(tmp2,tmp)
# Sinfo.tcga <- Sinfo.tcga[tmp2,]
# rownames(Sinfo.tcga) <- names(sp.cluster2.tcga.ans$group)[index] # final Sinfo with BLCA-XXXX-01
# 
# Sinfo.tcga.gaby <- read.table(file.path(data.path,"BLCA_ClinicalInfo_Gaby.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# rownames(Sinfo.tcga.gaby) <- paste0("BLCA-",substr(Sinfo.tcga.gaby$`Case ID`,9,12))
# tmp <- sp.cluster2.enhancer.tcga.ans$group
# tmp <- substr(names(tmp),start = 1,stop = 9)
# tmp1 <- rownames(Sinfo.tcga.gaby)
# tmp2 <- intersect(tmp,tmp1)
# index <- match(tmp2,tmp)
# Sinfo.tcga.gaby <- Sinfo.tcga.gaby[tmp2,]
# rownames(Sinfo.tcga.gaby) <- names(sp.cluster2.enhancer.tcga.ans$group)[index] # final Sinfo with BLCA-XXXX-01
# Sinfo.tcga.gaby$Enhancer_cluster <- sp.cluster2.enhancer.tcga.ans$group[rownames(Sinfo.tcga.gaby)]

#names(tmp) <- substr(names(tmp),start = 1,stop = 9)
#Sinfo.tcga <- read.table(file.path(data.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,fill = T)
#Sinfo.tcga <- Sinfo.tcga[which(Sinfo.tcga$type == "BLCA"),]
#rownames(Sinfo.tcga) <- Sinfo.tcga$SampleName
#kmplotfun(group = tmp,Sinfo = Sinfo.tcga,fig.path = fig.path,outFigFile="KM plot supervised cluster of TCGA BLCA using DMP enhancer.pdf", xunit="Year")
# $`summary`
# records events   *rmean *se(rmean)   median  0.95LCL  0.95UCL
# factor(group)=cluster1     247    126 4.877483  0.3788050 2.287671 1.876712 3.693151
# factor(group)=cluster2     164     54 6.001688  0.5803412 7.235616 2.600000 8.720548
# 
# $pval
# [1] "0.011"

#cox HR = 0.6598

# multivariateCox
# rt <- as.data.frame(t(enhancerForSpCluster.tcga))
# rownames(rt) <- substr(rownames(rt),start = 1,stop = 9)
# rt$futime <- Sinfo.tcga[rownames(rt),"OS.time"]/365
# rt$fustat <- Sinfo.tcga[rownames(rt),"OS"]
# rt <- rt[,c(1527:1528,1:1526)]
# outTab=data.frame()
# for(i in colnames(rt[,3:ncol(rt)])){
#   cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
#   coxSummary = summary(cox)
#   outTab=rbind.data.frame(outTab,cbind.data.frame(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
#                             z=coxSummary$coefficients[,"z"],
#                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]),stringsAsFactors = F)
# }
# outTab$GeneSymbol <- PbInfo[outTab$gene]
# write.table(outTab,file.path(res.path,"UnivariateCox for strict enhancers.txt"),sep = "\t",row.names = F)
# geneList <- as.character(outTab[which(outTab$pvalue<0.05),"gene"])
# rt2 <- rt[,c("futime","fustat",geneList)]
# cox <- coxph(Surv(futime, fustat) ~ ., data = rt2)
# cox.origin <- cox
# cox <- step(cox,direction = "both")
# cox.step.both <- cox
# coxSummary.origin <- summary(cox.origin)
# outTab.cox.origin <- cbind.data.frame(HR=coxSummary.origin$coefficients[,"exp(coef)"],
#                                         z=coxSummary.origin$coefficients[,"z"],
#                                         pvalue=coxSummary.origin$coefficients[,"Pr(>|z|)"])
#outTab.cox.origin <- outTab.cox.origin[which(outTab.cox.origin$pvalue<0.05),]
#anno <- read.table(file.path(comAnn.path,"Methylation_probes_annotation_simplified.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
#PbInfo <- anno$Simple_UCSC_RefGene_Name; names(PbInfo) <- rownames(anno)
#rm(anno)
# outTab.cox.origin$GeneSymbol <-PbInfo[rownames(outTab.cox.origin)]
# write.table(outTab.cox.origin,file.path(res.path,"multivariateCox for univariate significant strict enhancers.txt"),sep = "\t",row.names = T,col.names = NA)
# 
# supervised clustering by using CpGlist.Glist
load("F:/Project/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C1vsC2/ChAMP_RES_NOFILTER_C1vsC2/CpGlist.Glist.rda")
pval <- vector()
for (i in 1:length(CpGlist.Glist)) {
  tmp <- CpGlist.Glist[[i]]
  label <- names(CpGlist.Glist)[i]
  pb <- intersect(rownames(orgmeth.tcga),tmp)
  if(length(pb) < 300) {cat(paste0(label,"'s probe is less than 300!\n")); next()}
  indata <- as.data.frame(na.omit(orgmeth.tcga[pb,tumorsamples.tcga]))
  N.cluster <- 2
  N.bootstrap <- 500
  N.gene.per.bootstrap <- round(0.8*nrow(indata))
  N.sample.per.bootstrap <- round(0.8*ncol(indata))
  map.res.path <- file.path(res.path, paste(label,"_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
  featType <- label
  options(warn=-1)
  cluster2.CpGlist.tcga.ans <- plot.common.cluster(indata, 
                                                   tumorname="BLCA", 
                                                   N.cluster=N.cluster, 
                                                   N.bootstrap=N.bootstrap, 
                                                   N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                   N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                   map.res.path=map.res.path, fig.path=fig.path, 
                                                   featType=featType, 
                                                   seed=123456, dist0="manhattan", dist="pearson",link0 = "ward.D2",link = "ward.D2",
                                                   clstCol=c(white,purple), 
                                                   namecvt=NULL, height = 6, fontsize=6, labRow = F, labCol = F)
  
  tmp <- sp.cluster2.enhancer.tcga.ans$group
  #names(tmp) <- substr(names(tmp),start = 1,stop = 9)
  #Sinfo.tcga <- read.table(file.path(data.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,fill = T)
  #Sinfo.tcga <- Sinfo.tcga[which(Sinfo.tcga$type == "BLCA"),]
  #rownames(Sinfo.tcga) <- Sinfo.tcga$SampleName
  km <- kmplotfun(group = tmp,Sinfo = Sinfo.tcga,fig.path = fig.path,outFigFile=paste0("KM plot supervised cluster of TCGA BLCA using ",label, ".pdf"), xunit="Year")
  p <- km$pval; names(p) <- label; cat(paste0(label," KM pvalue = "),p,"\n")
  pval <- append(pval,p)
}

# > pval
# hyperC1_hyperN1.24         hyperC1_nodiff.T24vsN1          hyperC1_nodiff.T3vsN1 
# "0.018"                         "0.26"                         "0.12" 
# hyperC1_nodiff.T3vsN1.enhancer            hyperC1_nodiff.TvsN   hyperC1_nodiff.TvsN.enhancer 
# "0.14"                         "0.26"                         "0.23" 

###################################
#### create clinical table ########
#### using Sinfo.tcga.gaby ########
clinic.label <- c("gender",
                  "Age at diagnosis",
                  "AJCC Tumor category",
                  "Lymphovascular invasion",
                  "AJCC LN category",
                  "Number of LNs examined",
                  "AJCC metastasis category",
                  "AJCC pathologic tumor stage",
                  "Histologic grade",
                  "Histologic subtype",
                  "Squamous pathology",
                  "Neuroendocrine pathology",
                  "Plasmacytoid pathology")
cluster.label <- c("Hypermethylation cluster",
                   "Hypomethylation cluster",
                   "mRNA cluster",
                   "microRNA cluster",
                   "lncRNA cluster",
                   "RPPA cluster",
                   "Mutation process cluster")
mut.cnv.label <- colnames(Sinfo.tcga)[c(95:152,172:204)]

mytmp <- Sinfo.tcga.gaby[,c(clinic.label,cluster.label,"Enhancer_cluster")]
tmp <- mytmp
median.number.LNs <- median(as.numeric(tmp$`Number of LNs examined`[which(!(tmp$`Number of LNs examined` %in% "ND"))])) #18
for (i in 1:nrow(tmp)) {
  if(tmp$`Number of LNs examined`[i] == "ND") {
    next()}else {if(as.numeric(tmp$`Number of LNs examined`[i]) > median.number.LNs ){
      tmp$`Number of LNs examined`[i] <- ">medianLNs"
    } else  {
      tmp$`Number of LNs examined`[i] <- "<=medianLNs"}
  }
}
  
table(tmp$`Histologic grade`,tmp$Enhancer_cluster)      
fisher.test(matrix(c(120,22,9,17,61,17,12,127,13,7),ncol = 2,byrow = T),workspace = 1e7)

#############################################################################################################################
## perform an univariate and multivariate analysis for BLCA TCGA using the following variables:
## - age, sex, stage and our enhancer classification to see if our enhancer signature is independent from these parameters.

# Sinfo.tcga$Enhancer_cluster <- sp.cluster2.enhancer.tcga.ans$group[rownames(Sinfo.tcga)]
# Sinfo.tcga <- cbind(Sinfo.tcga,Sinfo.tcga.gaby[rownames(Sinfo.tcga),])
rt <- data.frame("futime"=Sinfo.tcga$OS.time/365,"fustat"=Sinfo.tcga$OS)
rt <- cbind(rt,Sinfo.tcga[,c(clinic.label,cluster.label,"DMPsCluster")])
rt$`Age at diagnosis` <- ifelse(as.numeric(rt$`Age at diagnosis`) > median(as.numeric(rt$`Age at diagnosis`)),"UPMedian","DNMedian")
for (i in 1:nrow(rt)) {
  if(rt$`Number of LNs examined`[i] == "ND") {
    next()}else {if(as.numeric(rt$`Number of LNs examined`[i]) > median.number.LNs ){
      rt$`Number of LNs examined`[i] <- "UPMedian"
    } else  {
      rt$`Number of LNs examined`[i] <- "DNMedian"}
    }
}
for (i in 1:nrow(rt)) {
  if(rt[i,"AJCC Tumor category"] %in% c("T1","T2")){
    rt[i,"AJCC Tumor category"] <- "T1_T2"
  } else if (rt[i,"AJCC Tumor category"] %in% c("T3","T3b")){
    rt[i,"AJCC Tumor category"] <- "T3"
  } else if (rt[i,"AJCC Tumor category"] %in% c("T4","T4a")){
    rt[i,"AJCC Tumor category"] <- "T4"
  }  else (rt[i,"AJCC Tumor category"] <- "ND")
}
for (i in 1:nrow(rt)) {
  if(rt[i,"AJCC LN category"] %in% c("N2","N3")){
    rt[i,"AJCC LN category"] <- "N2_N3"
  } else if (rt[i,"AJCC LN category"] %in% c("ND","NX")){
    rt[i,"AJCC LN category"] <- "ND"
  }
}

for (i in 1:nrow(rt)) {
  if(rt[i,"AJCC metastasis category"] %in% c("MD","MX")){
    rt[i,"AJCC metastasis category"] <- "ND"
  }
}
for (i in 1:nrow(rt)) {
  if(rt[i,"AJCC pathologic tumor stage"] %in% c("I","II")){
    rt[i,"AJCC pathologic tumor stage"] <- "I_II"
  }
}
#rt <- rt[,c(1,2,4,7,8,6)]
#rt <- rt[which(!(rt$Stage %in% NACols)),]
univariate.outTab=data.frame()
NACols <- c("[Not Applicable]","[Not Available]","Unknown","ND") #because all patients are alive in stage I
for(i in colnames(rt[,3:ncol(rt)])){
  tmp <- as.data.frame(na.omit(rt[,c("futime","fustat",i)]))
  tmp <- tmp[which(!(tmp[,i] %in% NACols)),]
  cox <- coxph(Surv(futime, fustat) ~ tmp[,i], data = tmp)
  #cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  
  coxSummary = summary(cox)
  univariate.outTab=rbind.data.frame(univariate.outTab,cbind.data.frame(feature=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                                                        z=coxSummary$coefficients[,"z"],
                                                                        pvalue=coxSummary$coefficients[,"Pr(>|z|)"]),stringsAsFactors = F)
}
#rownames(univariate.outTab) <- c("Gender","Age","StageIII","StageIV","Enhancer_cluster")
write.table(univariate.outTab,file.path(res.path,"univariateCox for clinical information and DMPsCluster.txt"),sep = "\t",row.names = T,col.names = NA)

#rt <- data.frame("futime"=Sinfo.tcga$OS.time/365,"fustat"=Sinfo.tcga$OS)
#rt <- cbind(rt,Sinfo.tcga[,c("age_at_initial_pathologic_diagnosis","gender","ajcc_pathologic_tumor_stage","Enhancer_cluster")])
#rt <- rt[which(!(rt$ajcc_pathologic_tumor_stage %in% NACols)),]
rt2 <- rt[,c(1,2,4,10,12,16)]
tmp <- c()
for (i in 3:ncol(rt2)) {
  tmp <- union(tmp,which(rt2[,i] %in% "ND"))
}
rt2 <- rt2[-tmp,]
cox <- coxph(Surv(futime, fustat) ~ ., data = rt2)
coxSummary = summary(cox)
multivariate.outTab=cbind.data.frame(HR=coxSummary$coefficients[,"exp(coef)"],
                                     z=coxSummary$coefficients[,"z"],
                                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
write.table(multivariate.outTab,file.path(res.path,"multivariateCox for clinical information and DMPsCluster.txt"),sep = "\t",row.names = T,col.names = NA)

#Indepent test for mut.cnv label
p <- c()
tmp <- Sinfo.tcga[,c(mut.cnv.label,"DMPsCluster")]
for (i in colnames(tmp)) {
  if(i == "DMPsCluster") {next()}
  tmp1 <- tmp[,c(i,"DMPsCluster")]
  tmp2 <- fisher.test(table(tmp1[,1],tmp1[,2]))$p.value
  names(tmp2) <- i
  p <- c(p,tmp2)
}
write.table(data.frame(p),file.path(res.path,"Independent test for MutCnv.ClinicInfo with DMPsCluster.txt"),sep = "\t",row.names = T,col.names = NA)
sig.mut.cnv <- names(p[p<0.05])

tmp3 <- data.frame()
tmp <- Sinfo.tcga[,c(sig.mut.cnv,"DMPsCluster")]
for (i in sig.mut.cnv) {
  genename <- strsplit(i," ")[[1]]
  genename <- genename[length(genename)]
  tmp1 <- tmp[,c(i,"DMPsCluster")]
  tmp2 <- as.data.frame.array(table(tmp1[,1],tmp1[,2]))
  if(grepl("mutation",i)){
    rownames(tmp2) <- paste0(genename,"_Mut_",rownames(tmp2))
  } else if (grepl("amplification",i)) {
    rownames(tmp2) <- paste0(genename,"_Amp_",rownames(tmp2))
  } else if (grepl("deletion",i)) {
    rownames(tmp2) <- paste0(genename,"_Del_",rownames(tmp2))
}
  tmp3 <- rbind.data.frame(tmp3,tmp2)
}
write.table(tmp3,file.path(res.path,"Details for significant MutCnv.ClinicInfo with DMPsCluster.txt"),sep = "\t",row.names = T,col.names = NA)
##################################################
## create annCol for BLCA C1 vs C2 ###############
#DMBII score
DMB <- DMBIIscore[which(substr(DMBIIscore$Name,1,12) %in% tumorsamples.tcga),]
DMB$Name <- substr(DMB$Name,1,12)
DMB <- DMB[!duplicated(DMB$Name),]
DMBtmp <- DMB[,-c(1:2)]; names(DMBtmp) <- DMB$Name
DMBtmp <- annTrackScale(indata = DMBtmp, halfwidth = 2, poolsd = F)

annCol.C1vsC2.tcga <- rt[,3:15]
annCol.C1vsC2.tcga <- cbind(annCol.C1vsC2.tcga,Sinfo.tcga[rownames(annCol.C1vsC2.tcga),c(cluster.label,sig.mut.cnv)])
annCol.C1vsC2.tcga$`Hypermethylation cluster` <- as.character(annCol.C1vsC2.tcga$`Hypermethylation cluster`)
annCol.C1vsC2.tcga$`Hypomethylation cluster` <- as.character(annCol.C1vsC2.tcga$`Hypomethylation cluster`)
annCol.C1vsC2.tcga$`Mutation process cluster` <- as.character(annCol.C1vsC2.tcga$`Mutation process cluster`)
annCol.C1vsC2.tcga$`DMI scores` <- DMBtmp[rownames(annCol.C1vsC2.tcga)]

annColors.C1vsC2.tcga <- list()
annColors.C1vsC2.tcga[["DMI scores"]] <- bluered(64)
annColors.C1vsC2.tcga[["gender"]] <- c("MALE"=blue,"FEMALE"=sun)
annColors.C1vsC2.tcga[["Age at diagnosis"]] <- c("UPMedian"=purple,"DNMedian"="white")
annColors.C1vsC2.tcga[["AJCC Tumor category"]] <- c("ND"=lightgrey,"T1_T2"=nake,"T3"=seagreen,"T4"=sun)
annColors.C1vsC2.tcga[["Lymphovascular invasion"]] <- c("ND"=lightgrey,"NO"="white","YES"=purple)
annColors.C1vsC2.tcga[["AJCC LN category"]] <- c("ND"=lightgrey,"N0"=nake,"N1"=seagreen,"N2_N3"=sun)
annColors.C1vsC2.tcga[["Number of LNs examined"]] <- c("UPMedian"=purple,"DNMedian"="white","ND"=lightgrey)
annColors.C1vsC2.tcga[["AJCC metastasis category"]] <- c("ND"=lightgrey,"M0"=blue,"M1"=sun)
annColors.C1vsC2.tcga[["AJCC pathologic tumor stage"]] <- c("ND"=lightgrey,"I_II"=nake,"III"=seagreen,"IV"=sun)
annColors.C1vsC2.tcga[["Histologic grade"]] <- c("ND"=lightgrey,"High Grade"=purple,"Low Grade"="white")
annColors.C1vsC2.tcga[["Histologic subtype"]] <- c("ND"=lightgrey,"Non-Papillary"="white","Papillary"=purple)
annColors.C1vsC2.tcga[["Squamous pathology"]] <- c("no"="white","yes"=purple)
annColors.C1vsC2.tcga[["Neuroendocrine pathology"]] <- c("no"="white","yes"=purple)
annColors.C1vsC2.tcga[["Plasmacytoid pathology"]] <- c("no"="white","yes"=purple)
annColors.C1vsC2.tcga[["Enhancer_cluster"]] <- c("cluster1"=blue,"cluster2"=red)
annColors.C1vsC2.tcga[["Hypermethylation cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold,"5"=brown)
annColors.C1vsC2.tcga[["Hypomethylation cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold,"5"=brown)
annColors.C1vsC2.tcga[["microRNA cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold,"ND"=lightgrey)
annColors.C1vsC2.tcga[["lncRNA cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold,"ND"=lightgrey)
annColors.C1vsC2.tcga[["RPPA cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold,"5"=brown,"ND"=lightgrey)
annColors.C1vsC2.tcga[["Mutation process cluster"]] <- c("1"=blue,"2"=sun,"3"=green,"4"=gold)
annColors.C1vsC2.tcga[["mRNA cluster"]] <- c("Basal_squamous"=blue,"Luminal"=sun,"Luminal_papillary"=green,"Neuronal"=brown,"Luminal_infiltrated"=gold,"ND"=lightgrey)

annColors.C1vsC2.tcga[["mutation in TP53"]] <- annColors.C1vsC2.tcga[["mutation in RB1"]] <- annColors.C1vsC2.tcga[["mutation in CDKN1A"]] <-
annColors.C1vsC2.tcga[["mutation in KDM6A"]] <- annColors.C1vsC2.tcga[["mutation in ELF3"]] <- annColors.C1vsC2.tcga[["mutation in PSIP1"]] <-
annColors.C1vsC2.tcga[["mutation in STAG2"]] <- annColors.C1vsC2.tcga[["mutation in KMT2A"]] <- annColors.C1vsC2.tcga[["mutation in FGFR3"]] <-
annColors.C1vsC2.tcga[["mutation in ERBB2"]] <- annColors.C1vsC2.tcga[["mutation in TMCO4"]] <- c("no"="white","yes"=purple)

annColors.C1vsC2.tcga[["focal amplification in YWHAZ"]] <- annColors.C1vsC2.tcga[["focal amplification in CCND1"]] <- annColors.C1vsC2.tcga[["focal amplification in PPARG"]] <-
annColors.C1vsC2.tcga[["focal amplification in ERBB2"]] <- annColors.C1vsC2.tcga[["focal amplification in BCL2L1"]] <- annColors.C1vsC2.tcga[["focal amplification in FGFR3"]] <-
annColors.C1vsC2.tcga[["focal amplification in AHR"]] <- c("."="white","Gain"=lightred,"Amp"=darkred,"NaN"=lightgrey)

annColors.C1vsC2.tcga[["focal deletion in CDKN2A"]] <- annColors.C1vsC2.tcga[["focal deletion in PDE4D"]] <- annColors.C1vsC2.tcga[["focal deletion in WWOX"]] <-
annColors.C1vsC2.tcga[["focal deletion in PTEN"]] <- annColors.C1vsC2.tcga[["focal deletion in NCOR1"]] <- annColors.C1vsC2.tcga[["focal deletion in RAD51B"]] <-
annColors.C1vsC2.tcga[["focal deletion in PTPRD"]] <- annColors.C1vsC2.tcga[["focal deletion in CCSER1"]] <- annColors.C1vsC2.tcga[["focal deletion in CREBBP"]] <-
c("."="white","Loss"=blue,"Del"=darkblue,"NaN"=lightgrey)

# for(i in 28:40){
#   annCol.C1vsC2.tcga[,i] <- gsub("\\.","ND",annCol.C1vsC2.tcga[,i])
# }

#calclulate MeTIL signature markers NPCA score

MeTIL.marker <- c("cg20792833","cg20425130","cg23642747","cg12069309","cg21554552")
MeTIL <- orgmeth.tcga[MeTIL.marker,tumorsamples.tcga]
MeTIL <- t(scale(t(MeTIL)))
pca.MeTIL <- prcomp(MeTIL,center = F,scale. = F)
MeTIL.score <- annTrackScale(indata = pca.MeTIL$rotation[,1], halfwidth = 2, poolsd = F)

annCol.C1vsC2.tcga$MeTIL <- MeTIL.score[rownames(annCol.C1vsC2.tcga)]
annCol.C1vsC2.tcga$pca.MeTIL <- pca.MeTIL$rotation[,1][rownames(annCol.C1vsC2.tcga)]
annColors.C1vsC2.tcga[["MeTIL"]] <- bluered(64)
write.table(annCol.C1vsC2.tcga,"H:/UTUC_BLCA/Results/annCol.C1vsC2.tcga.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# test if 9p21.3 enriched in C1 or C2 in UTUC data
tmp <- read.table(file.path(tumor.path,"CNA/GISTIC/1732580/CNA_segment_forGISTIC2.0.all_lesions.conf_90.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- tmp[which(tmp$Descriptor == "9p21.3  "),9:43][1,]
colnames(tmp) <- substr(colnames(tmp),start = 1,stop = 8)
tmp <- as.data.frame(t(tmp))
tmp$DMBcluster <- UTUC.annotation[substr(rownames(tmp),start = 6,stop = 8),"Cluster DNA methylation no filter without normal"]
cnv9p21.3 <- tmp
table(cnv9p21.3$`Deletion Peak 2`,cnv9p21.3$DMBcluster)
# C1 C2
# 0 18  8
# 1  5  4
fisher.test(table(cnv9p21.3$`Deletion Peak 2`,cnv9p21.3$DMBcluster))
# Fisher's Exact Test for Count Data
# 
# data:  table(cnv9p21.3$`Deletion Peak 2`, cnv9p21.3$DMBcluster)
# p-value = 0.6855
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.2737408 10.9489797
# sample estimates:
# odds ratio 
#    1.76859 

##########################################################################
######## link mutation to methylation by SWI/SNF and MAPK pathway ########

# tcga.mut <- read.table(file.path(data.path,"BLCA_Mutation_forJack.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# colnames(tcga.mut) <- substr(colnames(tcga.mut),start = 1,stop = 12)
# samples.tcga.mut <- intersect(colnames(tcga.mut),names(sp.cluster2.tcga.ans$group))
# 
# ###########
# # SWI/SNF #
# SWI_SNF <- c("ARID1A",	"PBRM1",	"SMARCA4",	"ARID1B",	"ARID2",
#              "SMARCA2",	"SMARCB1",	"SMARCC1",	"SMARCC2",	"DPF2",
#              "ACTL6A",	"SMARCD3",	"PHF10",	"BRD7",	"SMARCD1",
#              "SMARCD2",	"SMARCE1",	"ACTL6B",	"DPF1",	"DPF3")
# #SWI_SNF_6gene <- c("ARID1A","ARID1B","DPF3","PBRM1","SMARCA2","SMARCA4")
# 
# tcga.mut.swisnf <- tcga.mut[intersect(SWI_SNF,rownames(tcga.mut)),]
# #tcga.mut.swisnf_6gene <- tcga.mut[intersect(SWI_SNF_6gene,rownames(tcga.mut)),]
#   
# #DMBspcluster C1 and C2
# group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
# plotdata <- as.data.frame(t(tcga.mut.swisnf[,names(group)]))
# plotdata$Class <- group[rownames(plotdata)]
# plotdata <- apply(plotdata[,1:16], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
# plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
# plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/77,2),round(plotdata$Count/53,2))
# p2 <- ggplot() + 
#   geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(16)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
#   scale_y_continuous(breaks=seq(0, 0.9, 0.1), limits=c(0, 0.9)) + theme(legend.position="bottom")
# ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 16genes Pct between TCGA DMBspcluster C1 and C2.pdf"),width = 6)

# Using cBioportal mutation data

# #DMBspcluster C1 and C2
tcga.mut.swisnf <- read.table(file.path(data.path,"BLCA_SWI_SNF_mutation_cBioportal.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
tcga.mut.swisnf[is.na(tcga.mut.swisnf)] <- ""
tcga.mut.swisnf[grepl("passenger",tcga.mut.swisnf)] <- ""
tcga.mut.swisnf <- as.data.frame(t(apply(tcga.mut.swisnf,1,function(x) {ifelse(x!="",1,0)})))
tcga.mut.swisnf <- tcga.mut.swisnf[rowSums(tcga.mut.swisnf)>0,]
colnames(tcga.mut.swisnf) <- paste0("BLCA",substr(colnames(tcga.mut.swisnf),start = 8,stop = 12),"-01")
samples.tcga.mut <- intersect(colnames(tcga.mut.swisnf),names(sp.cluster2.tcga.ans$group))

group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
plotdata <- as.data.frame(t(tcga.mut.swisnf[,names(group)]))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:5], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/220,2),round(plotdata$Count/192,2))
p2 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 0.3, 0.1), limits=c(0, 0.3)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 5genes putative driver Pct between TCGA DMBspcluster C1 and C2.pdf"),width = 6)

tcga.mut.swisnf.status <- data.frame("Mut"=colSums(tcga.mut.swisnf),"Class"=group[colnames(tcga.mut.swisnf)],row.names = colnames(tcga.mut.swisnf))
tcga.mut.swisnf.status$Mut <- ifelse(tcga.mut.swisnf.status$Mut >= 1,1,0)
table(tcga.mut.swisnf.status$Mut,tcga.mut.swisnf.status$Class)
# cluster1 cluster2
# 0      177      165
# 1       43       27
fisher.test(table(tcga.mut.swisnf.status$Mut,tcga.mut.swisnf.status$Class)) #0.1498

# tcga.mut <- read.table("C:/Users/Sugus/Desktop/mutation_bcgsc_gene",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# colnames(tcga.mut) <- paste0("BLCA",substr(colnames(tcga.mut),start = 8,stop = 15))
# tcga.mut.swisnf <- tcga.mut[intersect(rownames(tcga.mut),SWI_SNF),]
# samples.tcga.mut <- intersect(colnames(tcga.mut.swisnf),names(sp.cluster2.tcga.ans$group))
# 
# group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
# plotdata <- as.data.frame(t(tcga.mut.swisnf[,names(group)]))
# plotdata$Class <- group[rownames(plotdata)]
# plotdata <- apply(plotdata[,1:20], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
# plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
# plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/214,2),round(plotdata$Count/190,2))
# p2 <- ggplot() + 
#   geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(20)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
#   scale_y_continuous(breaks=seq(0, 0.8, 0.1), limits=c(0, 0.8)) + theme(legend.position="bottom")
# ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 20genes Pct between XENA TCGA DMBspcluster C1 and C2.pdf"),width = 6)

################
# MAPK pathway #
# MAPK <- read.table(file.path(data.path,"MAPK.geneset.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# MAPK <- MAPK$`GENE LIST`
# 
# #DMBspcluster C1 and C2
# 
# tcga.mut.mapk <- tcga.mut[intersect(MAPK,rownames(tcga.mut)),]
# 
# group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
# plotdata <- as.data.frame(t(tcga.mut.mapk[,names(group)]))
# plotdata$Class <- group[rownames(plotdata)]
# plotdata <- apply(plotdata[,1:205], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
# plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
# plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/77,2),round(plotdata$Count/53,2))
# p4 <- ggplot() + 
#   geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(205)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
#   scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4)) + theme(legend.position="bottom")
# ggsave(file.path(fig.path,"Stacked bar of MAPK 205genes Pct between TCGA DMBspcluster C1 and C2.pdf"),width = 6,height = 15)

tcga.mut.mapk <- read.table(file.path(data.path,"BLCA_MAPK_mutation_cBioportal.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
tcga.mut.mapk[is.na(tcga.mut.mapk)] <- ""
tcga.mut.mapk[grepl("passenger",tcga.mut.mapk)] <- ""
tcga.mut.mapk <- as.data.frame(t(apply(tcga.mut.mapk,1,function(x) {ifelse(x!="",1,0)})))
tcga.mut.mapk <- tcga.mut.mapk[rowSums(tcga.mut.mapk)>0,]
colnames(tcga.mut.mapk) <- paste0("BLCA",substr(colnames(tcga.mut.mapk),start = 8,stop = 12),"-01")
samples.tcga.mut <- intersect(colnames(tcga.mut.mapk),names(sp.cluster2.tcga.ans$group))

group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
plotdata <- as.data.frame(t(tcga.mut.mapk[,names(group)]))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:6], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/220,2),round(plotdata$Count/192,2))
p4 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 0.2, 0.1), limits=c(0, 0.2)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of MAPK 6genes putative driver Pct between TCGA DMBspcluster C1 and C2.pdf"),width = 6,height = 15)

# tcga.mut.mapk <- tcga.mut[intersect(rownames(tcga.mut),MAPK),]
# group <- sp.cluster2.tcga.ans$group[samples.tcga.mut]
# plotdata <- as.data.frame(t(tcga.mut.mapk[,names(group)]))
# plotdata$Class <- group[rownames(plotdata)]
# plotdata <- apply(plotdata[,1:265], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
# plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
# plotdata$Percentage <- ifelse(plotdata$Class == "cluster1",round(plotdata$Count/214,2),round(plotdata$Count/190,2))
# p4 <- ggplot() + 
#   geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(267)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
#   scale_y_continuous(breaks=seq(0, 5, 0.1), limits=c(0, 5)) + theme(legend.position="bottom")
# ggsave(file.path(fig.path,"Stacked bar of MAPK 265genes Pct between XENA TCGA DMBspcluster C1 and C2.pdf"),width = 6,height = 15)

########################################################################
### Recalculate methylation percentage in all 35 methylation samples ###
# in tumor

Gmeth.noXY <- apply(orgmeth, 2, function(x) tapply(x, INDEX=factor(anno[rownames(orgmeth), "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
Gmeth.noXY <- as.data.frame(na.omit(Gmeth.noXY))
write.table(Gmeth.noXY,file = file.path(data.path,"UTUC_medianAvg_Gmeth_noXY.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

Gmeth.noXY.cpg.promoter <- apply(orgmeth[probe.cpg.promoter,], 2, function(x) tapply(x, INDEX=factor(anno[probe.cpg.promoter, "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
Gmeth.noXY.cpg.promoter <- as.data.frame(na.omit(Gmeth.noXY.cpg.promoter))
write.table(Gmeth.noXY.cpg.promoter,file = file.path(data.path,"UTUC_medianAvg_Gmeth_noXY_cpg_promoter.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

samples <- intersect( rownames(UTUCsinfo[which(UTUCsinfo$Type=="Tumor"),]), colnames(Gmeth.noXY.cpg.promoter) )
MethPct.noXY.cpg.promoter <- CalMethPct(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)

# in normal
normsamples <- intersect(colnames(Gmeth.noXY.cpg.promoter), rownames(UTUCsinfo[which(UTUCsinfo$Type=="Normal"),]))
NormMethPct.noXY.cpg.promoter <- CalMethPct(gdata=Gmeth.noXY.cpg.promoter, samples=normsamples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)

if(all(names(MethPct.noXY.cpg.promoter)==names(NormMethPct.noXY.cpg.promoter))) {
  write.table(data.frame(HighMethPctInTumor=MethPct.noXY.cpg.promoter, HighMethPctInNormal=NormMethPct.noXY.cpg.promoter), file=file.path(res.path, "high_methylation_pct_of_genes_noXY_cpg_promoter_UTUC.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
}else{
  stop("MethPct and NormMethPct names mismatch")
}

# calulate methylation percentage regarding FGFR3 mutation and muscle-invasiveness
# FGFR3
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$FGFR3=="1"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp1 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$FGFR3=="0"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp2 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
Gmeth.noXY.cpg.promoter.FGFR3 <- cbind.data.frame(tmp1,tmp2)
colnames(Gmeth.noXY.cpg.promoter.FGFR3) <- paste0(c("FGFR3_Mutated_","FGFR3_Mutated_","FGFR3_Wild_","FGFR3_Wild_"),colnames(Gmeth.noXY.cpg.promoter.FGFR3))

p <- c()
for (i in 1:nrow(Gmeth.noXY.cpg.promoter.FGFR3)) {
  tmp <- Gmeth.noXY.cpg.promoter.FGFR3[i,]
  tmp <- matrix(c(tmp[1,1],tmp[1,2],tmp[1,3],tmp[1,4]),byrow = T,ncol = 2)
  p <- c(p,fisher.test(tmp)$p.value)
}
Gmeth.noXY.cpg.promoter.FGFR3$nominal.p <- p
Gmeth.noXY.cpg.promoter.FGFR3$FDR <- p.adjust(p,method = "fdr")
write.table(Gmeth.noXY.cpg.promoter.FGFR3,file = file.path("H:/UTUC_BLCA/Results","differentially_methylated_genes_regarding_FGFR3_noXY_cpg_promoter_UTUC.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

# FGFR3 curated
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut=="1"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp1 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut=="0"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp2 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
Gmeth.noXY.cpg.promoter.FGFR32 <- cbind.data.frame(tmp1,tmp2)
colnames(Gmeth.noXY.cpg.promoter.FGFR32) <- paste0(c("FGFR3_curated_Mutated_","FGFR3_curated_Mutated_","FGFR3_curated_Wild_","FGFR3_curated_Wild_"),colnames(Gmeth.noXY.cpg.promoter.FGFR32))

p <- c()
for (i in 1:nrow(Gmeth.noXY.cpg.promoter.FGFR32)) {
  tmp <- Gmeth.noXY.cpg.promoter.FGFR32[i,]
  tmp <- matrix(c(tmp[1,1],tmp[1,2],tmp[1,3],tmp[1,4]),byrow = T,ncol = 2)
  p <- c(p,fisher.test(tmp)$p.value)
}
Gmeth.noXY.cpg.promoter.FGFR32$nominal.p <- p
Gmeth.noXY.cpg.promoter.FGFR32$FDR <- p.adjust(p,method = "fdr")
write.table(Gmeth.noXY.cpg.promoter.FGFR32,file = file.path("H:/UTUC_BLCA/Results","differentially_methylated_genes_regarding_FGFR3_curated_noXY_cpg_promoter_UTUC.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

# muscle-invasiveness
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$Type=="muscle-invasive"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp1 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$Type=="Non-muscle invasive"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp2 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
Gmeth.noXY.cpg.promoter.MI <- cbind.data.frame(tmp1,tmp2)
colnames(Gmeth.noXY.cpg.promoter.MI) <- paste0(c("MI_","MI_","NMI_","NMI_"),colnames(Gmeth.noXY.cpg.promoter.MI))

p <- c()
for (i in 1:nrow(Gmeth.noXY.cpg.promoter.MI)) {
  tmp <- Gmeth.noXY.cpg.promoter.MI[i,]
  tmp <- matrix(c(tmp[1,1],tmp[1,2],tmp[1,3],tmp[1,4]),byrow = T,ncol = 2)
  p <- c(p,fisher.test(tmp)$p.value)
}
Gmeth.noXY.cpg.promoter.MI$nominal.p <- p
Gmeth.noXY.cpg.promoter.MI$FDR <- p.adjust(p,method = "fdr")
write.table(Gmeth.noXY.cpg.promoter.MI,file = file.path("H:/UTUC_BLCA/Results","differentially_methylated_genes_regarding_MI_noXY_cpg_promoter_UTUC.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

# swi/snf
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$SWI_SNF_noACTL6B=="1"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp1 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
samples <- intersect( rownames(UTUC.annotation[which(UTUC.annotation$SWI_SNF_noACTL6B=="0"),]), colnames(Gmeth.noXY.cpg.promoter) )
tmp2 <- CalMethCount(gdata=Gmeth.noXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
Gmeth.noXY.cpg.promoter.SWISNF <- cbind.data.frame(tmp1,tmp2)
colnames(Gmeth.noXY.cpg.promoter.SWISNF) <- paste0(c("SWISNF_","SWISNF_","SWISNF_","SWISNF_"),colnames(Gmeth.noXY.cpg.promoter.SWISNF))

p <- c()
for (i in 1:nrow(Gmeth.noXY.cpg.promoter.SWISNF)) {
  tmp <- Gmeth.noXY.cpg.promoter.SWISNF[i,]
  tmp <- matrix(c(tmp[1,1],tmp[1,2],tmp[1,3],tmp[1,4]),byrow = T,ncol = 2)
  p <- c(p,fisher.test(tmp)$p.value)
}
Gmeth.noXY.cpg.promoter.SWISNF$nominal.p <- p
Gmeth.noXY.cpg.promoter.SWISNF$FDR <- p.adjust(p,method = "fdr")
write.table(Gmeth.noXY.cpg.promoter.SWISNF,file = file.path("H:/UTUC_BLCA/Results","differentially_methylated_genes_regarding_SWISNF_noXY_cpg_promoter_UTUC.txt"),row.names=T, col.names=NA, sep="\t", quote=F)

# probe.withXY.CpG <- rownames(anno[which(anno$Relation_to_UCSC_CpG_Island == "Island"),])
# probe.withXY.CpG.promoter <- intersect(probe.withXY.CpG,intersect(rownames(anno[which(anno$Simple_UCSC_RefGene_Group %in% c("TSS1500","TSS200")),]),rownames(orgmeth.withXY))) #62227
#   
# orgmeth.withXY <- read.csv(methFile,check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
# colnames(orgmeth.withXY) <- orgmeth.withXY["Sample ID",]
# orgmeth.withXY <- orgmeth.withXY[-c(1:3),]; colnames(orgmeth.withXY)[1] <- "01N" ;tmp <- rownames(orgmeth.withXY) #866091
# tmp <- rownames(orgmeth.withXY)
# orgmeth.withXY <- as.data.frame(sapply(orgmeth.withXY, as.numeric)); rownames(orgmeth.withXY) <- tmp
# orgmeth.withXY <- orgmeth.withXY[,setdiff(colnames(orgmeth.withXY),c("01N","01T","07T","11N","11T"))]
# orgmeth.withXY <- as.data.frame(na.omit(orgmeth.withXY)) #855447
# 
# Gmeth.withXY <- apply(orgmeth.withXY, 2, function(x) tapply(x, INDEX=factor(anno[rownames(orgmeth.withXY), "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
# Gmeth.withXY <- as.data.frame(na.omit(Gmeth.withXY))
# write.table(Gmeth.withXY,file = file.path(data.path,"UTUC_medianAvg_Gmeth_withXY.txt"),row.names=T, col.names=NA, sep="\t", quote=F)
# 
# Gmeth.withXY.cpg.promoter <- apply(orgmeth.withXY[probe.withXY.CpG.promoter,], 2, function(x) tapply(x, INDEX=factor(anno[probe.withXY.CpG.promoter, "Simple_UCSC_RefGene_Name", drop=T]), FUN=median, na.rm=TRUE))
# Gmeth.withXY.cpg.promoter <- as.data.frame(na.omit(Gmeth.withXY.cpg.promoter))
# write.table(Gmeth.withXY.cpg.promoter,file = file.path(data.path,"UTUC_medianAvg_Gmeth_withXY_cpg_promoter.txt"),row.names=T, col.names=NA, sep="\t", quote=F)
# 
# samples <- intersect( rownames(UTUCsinfo[which(UTUCsinfo$Type=="Tumor"),]), colnames(Gmeth.withXY.cpg.promoter) )
# MethPct.withXY.cpg.promoter <- CalMethPct(gdata=Gmeth.withXY.cpg.promoter, samples=samples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
# 
# # in normal
# normsamples <- intersect(colnames(Gmeth.withXY.cpg.promoter), rownames(UTUCsinfo[which(UTUCsinfo$Type=="Normal"),]))
# NormMethPct.withXY.cpg.promoter <- CalMethPct(gdata=Gmeth.withXY.cpg.promoter, samples=normsamples, betacut=0.2, medcut=0.2, highpct=0.5, methSinfo=UTUCsinfo, highMethCut=0.3)
# 
# if(all(names(MethPct.withXY.cpg.promoter)==names(NormMethPct.withXY.cpg.promoter))) {
#   write.table(data.frame(HighMethPctInTumor=MethPct.withXY.cpg.promoter, HighMethPctInNormal=NormMethPct.withXY.cpg.promoter), file=file.path(res.path, "high_methylation_pct_of_genes_withXY_cpg_promoter_UTUC.txt"), row.names=T, col.names=NA, sep="\t", quote=F)
# }else{
#   stop("MethPct and NormMethPct names mismatch")
# }
# write.table(t(data.frame(Gmeth.withXY.cpg.promoter["STAG2",])),file.path(res.path,"Gmeth_STAG2.txt"),sep = "\t")

##############################################################################################################################
### Differentially methylation analysis to genes that have highmethylation in Tumor but no methylation in Normal in C1vsC2 ###
tmp <- read.table(file.path(res.path,"high_methylation_pct_of_genes_noXY_cpg_promoter_UTUC.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- rownames(tmp[which(tmp$HighMethPctInTumor > 0 & tmp$HighMethPctInNormal == 0),])

UTUC.annotation <- read.table(file.path(data.path,"UTUC_annotation.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
C1 <- rownames(UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C1"),])
C2 <- rownames(UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C2"),])
group <- UTUC.annotation$`Cluster DNA methylation no filter without normal`; names(group) <- rownames(UTUC.annotation); group <- group[-which(group %in% "N/A")]
complist <- createList.C1vsC2(group = group)

indata <- cbind.data.frame("GeneSymbol"=gene,Gmeth.noXY.cpg.promoter[gene,names(group)])
methDEtest(indata=indata, res.path=res.path, complist=complist, featType="methylation_noXY.cpg.promoter", overwt=T)
 
################################################################################
### Determine methylation status by using highmethylation genes in C1 and C2 ### 
tmp <- Gmeth.noXY.cpg.promoter[gene,c(C1,C2)]
tmp[tmp >= 0.3] <- "Methylated"
tmp[tmp < 0.3] <- "Unmethylated"

tmp1 <- tmp[,C1]
tmp2 <- tmp[,C2]

outTable <- data.frame()
for (i in gene) {
  a = length(which(tmp1[i,] == "Methylated")) # methylated in C1
  b = length(which(tmp1[i,] == "Unmethylated")) # unmethylated in C1
  
  c = length(which(tmp2[i,] == "Methylated")) # methylated in C2
  d = length(which(tmp2[i,] == "Unmethylated")) # unmethylated in C2
  
  f = fisher.test(matrix(c(a,c,b,d),byrow = T,ncol = 2))
  tmp <- data.frame(MethInC1=a,UnmethInC1=b,MethInC2=c,UnmethInC2=d,pvalue=f$p.value,stringsAsFactors = F,row.names = i)
  outTable <- rbind.data.frame(outTable,tmp)
}
outTable$padj <- p.adjust(outTable$pvalue,method = "BH")
write.table(outTable,file.path(res.path,"Methylation status counts of noXY_PromCGIHighPctInTumor in C1 and C2 with independent test.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

#################################################
### derive results for swi/snf complex ACTL6B ###
tmp1 <- Gmeth.noXY.cpg.promoter["ACTL6B",c(C1,C2)]
ACTL6B <- data.frame(DMBcluster=rep(c("C1","C2"),time=c(length(C1),length(C2))),BetaValue=as.numeric(tmp1),MethStatus=ifelse(as.numeric(tmp1) >= 0.3,"Methylated","Unmethylated"),row.names = names(tmp1),stringsAsFactors = F)
ACTL6B.probes <- intersect(rownames(anno[which(anno$UCSC_RefGene_Name == "ACTL6B"),]),probe.cpg.promoter)
#"cg16024950" "cg17592984"
ACTL6B$cg16024950 <- as.numeric(orgmeth["cg16024950",c(C1,C2)])
ACTL6B$cg17592984 <- as.numeric(orgmeth["cg17592984",c(C1,C2)])
write.table(ACTL6B,file.path(res.path,"Methylation details for ACTL6B in C1 and C2.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

##########################################
### check if ACTL6B methylated in BLCA ###
anno.450k <- read.table(file.path(comAnn.path,"Annotation_PromCGI_450k.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,quote = "")
ACTL6B.probes.450k <- rownames(anno.450k[which(anno.450k$mRNA == "ACTL6B"),])

##########################################################################################
### Supervised clustering using DMPs derived from UTUC FGFR3 curated mutation and wild ###
DMPs.FGFR3_curated <- read.delim("F:/Project/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C12vsN/ChAMP_RES_NOFILTER_C12vsN/myDMP.FGFR3_curated.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
pb.FGFR3_curated.dn <- rownames(DMPs.FGFR3_curated[which(DMPs.FGFR3_curated$Wild_AVG >= 0.4 & DMPs.FGFR3_curated$Mutated_AVG <= 0.2 & DMPs.FGFR3_curated$adj.P.Val < 0.05),]) # 17687
pb.FGFR3_curated.up <- rownames(DMPs.FGFR3_curated[which(DMPs.FGFR3_curated$Wild_AVG <= 0.2 & DMPs.FGFR3_curated$Mutated_AVG >= 0.4 & DMPs.FGFR3_curated$adj.P.Val < 0.05),]) # 10

pb.FGFR3_curated <- intersect(rownames(orgmeth.tcga),c(pb.FGFR3_curated.dn,pb.FGFR3_curated.up))
length(pb.FGFR3_curated) #5655

indata <- as.data.frame(na.omit(orgmeth.tcga[pb.FGFR3_curated,tumorsamples.tcga])) # 4305
N.cluster <- 2
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))
#map.res.path <- file.path(res.path, paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_mpww_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_epww_BLCA_ClusterNum",N.cluster,sep = ""))
map.res.path <- file.path(res.path, paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_eeww_BLCA_ClusterNum",N.cluster,sep = ""))

featType <- "UTUC_FGFR3_Mutated_vs_Wild_DMP"
options(warn=-1)
sp.cluster2.tcga.fgfr3.mutatedvswild <- plot.common.cluster(indata, 
                                                            tumorname="BLCA", 
                                                            N.cluster=N.cluster, 
                                                            N.bootstrap=N.bootstrap, 
                                                            N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                            N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                            map.res.path=map.res.path, fig.path="H:/UTUC_BLCA/Figures", 
                                                            featType=featType,annCol = annCol.C1vsC2.tcga[colnames(indata),c("MeTIL","mutation in FGFR3")],annColors = annColors.C1vsC2.tcga,
                                                            seed=123456, dist0="euclidean", dist="euclidean",link0 = "ward.D",link = "ward.D",
                                                            clstCol=c("black",darkred), 
                                                            namecvt=NULL, height = 6, fontsize=9, labRow = F, labCol = F,cexAnn = 1)

tmp <- data.frame(cluster=sp.cluster2.tcga.fgfr3.mutatedvswild$group,
                  FGFR3 = annCol.C1vsC2.tcga[names(sp.cluster2.tcga.fgfr3.mutatedvswild$group),"mutation in FGFR3"])
table(tmp$cluster,tmp$FGFR3)
# no yes
# cluster1 264  13
# cluster2  90  45
fisher.test(table(tmp$cluster,tmp$FGFR3)) #3.926e-14

hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")

tmp <- annCol.C1vsC2.tcga; tmp$samples <- rownames(tmp)
# tmp$SupervisedCluster <- sp.cluster2.tcga.ans$group[rownames(tmp)]
# tmp$SupervisedCluster <- ifelse(tmp$SupervisedCluster == "cluster1","C1","C2")
tmp$Cluster <- ifelse(sp.cluster2.tcga.fgfr3.mutatedvswild$group[rownames(tmp)] == "cluster1","FGFR3-scarce","FGFR3-enriched")

# table(tmp$Cluster,tmp$SupervisedCluster)
# C1  C2
# FGFR3-enriched   0 135
# FGFR3-scarce   220  57
# fisher.test(table(tmp$Cluster,tmp$SupervisedCluster))# p < 2.2e-16

annColors.C1vsC2.tcga[["SupervisedCluster"]] <- c("C1"=red,"C2"=blue)
annColors.C1vsC2.tcga[["Cluster"]] <- c("FGFR3-enriched"=darkred,"FGFR3-scarce"="black")

tmp <- tmp[colnames(indata),c(49,28,52)]
colnames(tmp)[2] <- "FGFR3"
tmp$FGFR3 <- ifelse(tmp$FGFR3 == "yes","Mutated","Wild")
pdf(file.path(fig.path,"BLCA_heatmap_SupervisedCluster_UTUC_FGFR3_MutatedvsWild_DMPs.pdf"),height = 8.5,width = 7)
aheatmap(as.matrix(indata), 
         Rowv=dendsort(as.dendrogram(hcg)), 
         Colv=sp.cluster2.tcga.fgfr3.mutatedvswild$dendro, 
         annCol=tmp[colnames(indata),], # cnv
         annColors = append(annColors.C1vsC2.tcga[c("SupervisedCluster","Cluster","MeTIL")],
                            list("FGFR3"=c("Mutated"=cherry,"Wild"=lightgrey))), 
         #color=c("blue","blue","green","yellow","red","red"), 
         color=c("#0074FE","#96EBF9","#FEE900","#F00003"),
         revC=TRUE, fontsize=8,labRow = NA,labCol = NA,cexAnn = 1)
invisible(dev.off())

# MeTIL
tmp1 <- tmp[which(tmp$Cluster == "FGFR3-enriched"),"pca.MeTIL"]
tmp2 <- tmp[which(tmp$Cluster == "FGFR3-scarce"),"pca.MeTIL"]
wilcox.test(tmp1,tmp2) # <2.2e-16
df <- data.frame("MeTIL"=c(tmp1,tmp2),"Cluster"=rep(c("FGFR3_Mutated","FGFR3_Wild"),c(length(tmp1),length(tmp2))))
p2 <- ggplot(df,aes(x=Cluster,y=MeTIL,fill=Cluster)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(cherry,grey)) + 
  ggtitle("pValue < 2.2e-16") + 
  geom_jitter(aes(color = Cluster),alpha=1) +   scale_color_manual(values = c(cherry,grey)) + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave("F:/Project/UTUC_BLCA/Figures/pcaMeTIL between TCGA FGFR3Cluster by supervised clustering of UTUC DMPs .pdf",width = 4,height = 3.5)

#use OS.time
tmp <- sp.cluster2.tcga.fgfr3.mutatedvswild$group
Sinfo.tcga$FGFR3cluster <- tmp[rownames(Sinfo.tcga)]
tmp <- ifelse(tmp == "cluster1","FGFR3_Wild","FGFR3_Mutated"); names(tmp) <- names(sp.cluster2.tcga.fgfr3.mutatedvswild$group)
tmp <- factor(tmp,levels = c("FGFR3_Wild","FGFR3_Mutated"))
kmplotfun(group = tmp,Sinfo = Sinfo.tcga,fig.path = fig.path,outFigFile="KM plot supervised cluster of TCGA BLCA regarding UTUC FGFR3 mutation.pdf", xunit="Year",legend.title = "FGFR3cluster",col = c("black",darkred))

##########################################################################################
### Supervised clustering using DMPs derived from UTUC FGFR3 curated mutation and wild ###
DMPs.FGFR3_curated_vs_Normal <- read.delim("F:/Project/UTUC_BLCA/ChAMP/IDAT_NOFILTER_C12vsN/ChAMP_RES_NOFILTER_C12vsN/myDMP.FGFR3_curated_vs_Normal.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
pb.FGFR3_curated_vs_Normal.dn <- rownames(DMPs.FGFR3_curated_vs_Normal[which(DMPs.FGFR3_curated_vs_Normal$Normal_AVG >= 0.4 & DMPs.FGFR3_curated_vs_Normal$FGFR3_Mutated_AVG <= 0.2 & DMPs.FGFR3_curated_vs_Normal$adj.P.Val < 0.05),]) #57549
pb.FGFR3_curated_vs_Normal.up <- rownames(DMPs.FGFR3_curated_vs_Normal[which(DMPs.FGFR3_curated_vs_Normal$Normal_AVG <= 0.2 & DMPs.FGFR3_curated_vs_Normal$FGFR3_Mutated_AVG >= 0.4 & DMPs.FGFR3_curated_vs_Normal$adj.P.Val < 0.05),]) #4809

pb.FGFR3_curated_vs_Normal <- intersect(rownames(orgmeth.tcga),c(pb.FGFR3_curated_vs_Normal.dn,pb.FGFR3_curated_vs_Normal.up))
length(pb.FGFR3_curated_vs_Normal) # 22717

indata <- as.data.frame(na.omit(orgmeth.tcga[pb.FGFR3_curated_vs_Normal,c(rownames(Sinfo.tcga[which(Sinfo.tcga$`mutation in FGFR3` == "yes"),]),
                                                                          normsamples.tcga)]))
N.cluster <- 3
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))
#map.res.path <- file.path(res.path, paste("pb_FGFR3_carated_mutated_vs_wild_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_normal_0.8_0.8_mpw2w2_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_normal_0.8_0.8_mpww_BLCA_ClusterNum",N.cluster,sep = ""))
#map.res.path <- file.path("H:/UTUC_BLCA/Results", paste("pb_FGFR3_carated_mutated_vs_normal_0.8_0.8_epww_BLCA_ClusterNum",N.cluster,sep = ""))
map.res.path <- file.path(res.path, paste("pb_FGFR3_carated_mutated_vs_normal_0.8_0.8_eeww_BLCA_ClusterNum",N.cluster,sep = ""))

featType <- "UTUC_FGFR3_Mutated_vs_Normal_DMP"
options(warn=-1)
sp.cluster2.tcga.fgfr3.mutatedvsnormal <- plot.common.cluster(indata, 
                                                              tumorname="BLCA", 
                                                              N.cluster=N.cluster, 
                                                              N.bootstrap=N.bootstrap, 
                                                              N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                              N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                              map.res.path=map.res.path, fig.path="H:/UTUC_BLCA/Figures", 
                                                              featType=featType,
                                                              annCol = data.frame(Tissue = rep(c("BLCA","Normal"),c(58,21)),
                                                                                  FGFR3 = rep(c("Mutated","N/A"),c(58,21)),
                                                                                  row.names = colnames(indata)),
                                                              annColors = list("Tissue"=c("BLCA"=purple,"Normal"=white),
                                                                               "FGFR3"=c("Mutated"=cherry,"N/A"="white")),
                                                              seed=123456, dist0="euclidean", dist="euclidean",link0 = "ward.D",link = "ward.D",
                                                              clstCol=c("black",nake,seagreen), 
                                                              namecvt=NULL, height = 6, fontsize=9, labRow = F, labCol = F,cexAnn = 1)

tmp <- data.frame(Cluster=sp.cluster2.tcga.fgfr3.mutatedvsnormal$group,
                  Tissue = rep(c("BLCA","Normal"),c(58,21)))
tmp$Cluster <- ifelse(tmp$Cluster %in% c("cluster1","cluster2"),"FGFR3-like","Normal-like")
table(tmp$Cluster,tmp$Tissue)
# BLCA Normal
# FGFR3-like    55      0
# Normal-like    3     21
fisher.test(table(tmp$Cluster,tmp$Tissue)) #2.717e-16

hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")

pdf(file.path(fig.path,"BLCA_heatmap_SupervisedCluster_UTUC_FGFR3_MutatedvsNormal_DMPs.pdf"),height = 5,width = 6)
aheatmap(as.matrix(indata), 
         Rowv=dendsort(as.dendrogram(hcg)), 
         Colv=sp.cluster2.tcga.fgfr3.mutatedvsnormal$dendro, 
         annCol = data.frame(Tissue = rep(c("BLCA","Normal"),c(58,21)),
                             FGFR3 = rep(c("Mutated","N/A"),c(58,21)),
                             Cluster = tmp$Cluster,
                             row.names = colnames(indata)),
         annColors = list("Tissue"=c("BLCA"=purple,"Normal"=white),
                          "FGFR3"=c("Mutated"=cherry,"N/A"="white"),
                          "Cluster"=c("FGFR3-like"=nake,"Normal-like"=seagreen)),
         #color=c("blue","blue","green","yellow","red","red"), 
         color=c("#0074FE","#96EBF9","#FEE900","#F00003"),
         revC=TRUE, fontsize=8,labRow = NA,labCol = NA,cexAnn = 1)
invisible(dev.off())

# save
save.image(file.path(tumor.path,"Malouf.DNA2.RData"))
#save.image(file.path("H:/UTUC_BLCA","Malouf.DNA2.RData"))


##############################################
####### UTUC RNASeq analysis plan ############
##############################################

###########################################
############ 1. DATA PROCESS  #############
###########################################
HERVKInFile <- NULL
exp.file <- "UTUC_BLCA_raw_count_gencode15.txt"
SinfoFile <- "UTUC_BLCA_clinial_Info.txt"
res <- loadData(data.path, shareFun.path, script.path, comAnn.path, createGinfoFileFlag, SinfoFile, GinfoFile, HERVKInFile, GENEInFile = exp.file, tailrows)
Sinfo=res$Sinfo; Ginfo<-res$Ginfo; Vids<-res$Vids; Lids<-res$Lids; Mids<-res$Mids; countsTable<-res$countsTable

#batch <- rep(c("BLCA","UTUC"),times=c(426,20))
#countsTable <- removeBatchEffect(as.matrix(countsTable),batch=batch)
tmp <- as.data.frame(t(scale(t(countsTable))))
batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUC"),times=c(407,19,20)),fig.dir = fig.path,PCA.fig.title = "PCA.rawCounts.batch",cols = c(blue,red,yellow),showID = F,cex = 0.7,showLegend = T)

res <-normData(countsTable, tailrows, res.path, tumorname)
countsNorm <- res$countsNorm; vsd <- res$vsd

tmp <- as.data.frame(t(scale(t(countsNorm))))
batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUC"),times=c(407,19,20)),fig.dir = fig.path,PCA.fig.title = "PCA.countsNorm.batch",cols = c(blue,red,yellow),showID = F,cex = 0.7,showLegend = T)

FPKM <- calFPKM( countsTable, Vids, tailrows, Ginfo, outfile=file.path(res.path, paste(tumorname, ".FPKM.all.features.txt", sep="")))
tmp <- as.data.frame(t(scale(t(FPKM))))
batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUC"),times=c(407,19,20)),fig.dir = fig.path,PCA.fig.title = "PCA.FPKM.batch",cols = c(blue,red,yellow),showID = F,cex = 0.7,showLegend = T)

FPKM.BLCA <- FPKM[,1:407]
FPKM.BLCA$genename <- Ginfo[rownames(FPKM.BLCA),"genename"] 
FPKM.BLCA <- FPKM.BLCA[!duplicated(FPKM.BLCA$genename),]; rownames(FPKM.BLCA) <- FPKM.BLCA$genename; FPKM.BLCA <- FPKM.BLCA[,-ncol(FPKM.BLCA)]

FPKM.BLCA.Nromal <- FPKM[,408:426]
FPKM.BLCA.Nromal$genename <- Ginfo[rownames(FPKM.BLCA.Nromal),"genename"] 
FPKM.BLCA.Nromal <- FPKM.BLCA.Nromal[!duplicated(FPKM.BLCA.Nromal$genename),]; rownames(FPKM.BLCA.Nromal) <- FPKM.BLCA.Nromal$genename; FPKM.BLCA.Nromal <- FPKM.BLCA.Nromal[,-ncol(FPKM.BLCA.Nromal)]

PASSFlag.lnc <- LowExpFilter(FPKM, countsNorm, lowcut.lnc, lowpct, Mids, Lids, Ginfo, res.path, paste(tumorname,".lowcut",lowcut.lnc,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste(tumorname,".lowcut",lowcut.lnc,".lowpct",lowpct,".Count.mRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste(tumorname,".lowcut",lowcut.lnc,".lowpct",lowpct,".FPKM.mRNA.LowExpfiltered.txt",sep = "")))

PASSFlag.mRNA <- LowExpFilter(FPKM, countsNorm, lowcut.mRNA, lowpct, Mids, Lids, Ginfo, res.path, paste(tumorname,".lowcut",lowcut.mRNA,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste(tumorname,".lowcut",lowcut.mRNA,".lowpct",lowpct,".Count.lncRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste(tumorname,".lowcut",lowcut.mRNA,".lowpct",lowpct,".FPKM.lncRNA.LowExpfiltered.txt",sep = "")))

#remove batch effect by combat
rowids <- intersect(Mids, names(PASSFlag.mRNA[PASSFlag.mRNA==TRUE]) )
indata <- FPKM[rowids,]
Sinfo$batch <- rep(c("BLCA","UTUCrna","UTUCr"),times=c(426,14,6))
modcombat = model.matrix(~1, data=Sinfo)
combat.FPKM = ComBat(dat=as.matrix(indata), batch=Sinfo$batch, mod=modcombat)
tmp <- as.data.frame(t(scale(t(combat.FPKM))))
batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUCrna","UTUCr"),times=c(407,19,14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.combatFPKM.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)

#rowids <- intersect(Mids, names(PASSFlag.mRNA[PASSFlag.mRNA==TRUE]) )
#indata <- log2(FPKM[rowids,] + 1)
#combat.logFPKM = ComBat(dat=as.matrix(indata), batch=Sinfo$batch, mod=modcombat)
#tmp <- as.data.frame(t(scale(t(combat.logFPKM))))
#batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUCrna","UTUCr"),times=c(407,19,14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.combatlogFPKM.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)

#rowids <- intersect(Mids, names(PASSFlag.mRNA[PASSFlag.mRNA==TRUE]) )
#indata <- countsNorm[rowids,]
#combat.countsNorm = ComBat(dat=as.matrix(indata), batch=Sinfo$batch, mod=modcombat)
#tmp <- as.data.frame(t(scale(t(combat.countsNorm))))
#batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUCrna","UTUCr"),times=c(407,19,14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.combatcountsNorm.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)

#rowids <- intersect(Mids, names(PASSFlag.mRNA[PASSFlag.mRNA==TRUE]) )
#indata <- log2(countsNorm[rowids,] + 1)
#combat.logcountsNorm = ComBat(dat=as.matrix(indata), batch=Sinfo$batch, mod=modcombat)
#tmp <- as.data.frame(t(scale(t(combat.logcountsNorm))))
#batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("BLCA","Normal","UTUCrna","UTUCr"),times=c(407,19,14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.combatlogcountsNorm.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)

#create non-negative matrix for NMF or log input (only nneg.combat.FPKM and nneg.combat.countsNorm)
nneg.combat.FPKM <- pmax(combat.FPKM,0)*2
#nneg.combat.logFPKM <- pmax(combat.logFPKM,0)*2
#nneg.combat.countsNorm <- pmax(combat.countsNorm,0)*2
#nneg.combat.logcountsNorm <- pmax(combat.logcountsNorm,0)*2

#remove batch effect by normlize.quantile
#rowids <- intersect(Mids, names(PASSFlag.mRNA[PASSFlag.mRNA==TRUE]) )
#indata <- countsNorm[rowids,]
#norqua.countsNorm <- normalize.quantiles(as.matrix(indata)); rownames(norqua.countsNorm) <- rownames(indata); colnames(norqua.countsNorm) <- colnames(indata)
#combat.norqua.countsNorm = ComBat(dat=as.matrix(norqua.countsNorm), batch=Sinfo$batch, mod=modcombat)

#tmp <- as.data.frame(t(scale(t(combat.norqua.countsNorm))))
#batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("Basal_squamous","Luminal","Normal","UTUC"),times=c(141,266,19,20)),fig.dir = fig.path,PCA.fig.title = "PCA.combatnorquacountsNorm.Basal.Luminal.UTUC.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)

#norqua.combat.countsNorm <- normalize.quantiles(as.matrix(combat.countsNorm));rownames(norqua.combat.countsNorm) <- rownames(combat.countsNorm); colnames(norqua.combat.countsNorm) <- colnames(combat.countsNorm)
#tmp <- as.data.frame(t(scale(t(norqua.combat.countsNorm))))
#batchPCA(indata = tmp[,rownames(Sinfo)],batch = rep(c("Basal_squamous","Luminal","Normal","UTUC"),times=c(141,266,19,20)),fig.dir = fig.path,PCA.fig.title = "PCA.norquacombatcountsNorm.Basal.Luminal.UTUC.batch",cols = c(blue,red,yellow,green),showID = F,cex = 0.7,showLegend = T)
 
#1.	Perform unsupervised consensus clustering (K=5) by mRNA with combined TCGA Bladder RNASeq data (tumor only) with UTUC UTUC RNAseq data (please normalize the raw count first) with bar of TCGA_Subtype
samples <- intersect(rownames(Sinfo[which(Sinfo$SampleType=="PrimaryTumor"),]), colnames(countsTable))
annCol <- data.frame("Source"=Sinfo[samples,"batch"],"TCGA_Subtype"=Sinfo[samples,"TCGA_Subtype"],row.names = samples,stringsAsFactors = F)
annCol$Source <- ifelse(annCol$Source == "BLCA","BLCA","UTUC")
annCol$TCGA_Subtype[422:427] <- "UTUC_MI"
annColors <- list()
annColors[["Source"]] <- c("BLCA"=purple,"UTUC"=white)
annColors[["TCGA_Subtype"]] <- c("Basal_squamous"=darkred,"Luminal"=yellow,"Luminal_infiltrated"=green,"Luminal_papillary"=blue,"Neuronal"=lightred,"UTUC_MI"=purple,"UTUC_NMI"=grey)

#1
#indata <- combat.logFPKM[,samples]
#map.res.path <- file.path(tumor.path, paste("commonCluster_CombatLogFPKM_ClusterNum",N.cluster,sep = ""))
#featType <- "CombatLogFPKM"
#2 choosed!
indata <- log2(nneg.combat.FPKM[,samples] + 1)
map.res.path <- file.path(tumor.path, paste("commonCluster_LognnegCombatFPKM_batch3_ClusterNum",N.cluster,sep = ""))
featType <- "LognnegCombatFPKM"
#3
#indata <- nneg.combat.logFPKM[,samples]
#map.res.path <- file.path(tumor.path, paste("commonCluster_nnegCombatLogFPKM_ClusterNum",N.cluster,sep = ""))
#featType <- "nnegCombatLogFPKM"
#4
#indata <- combat.logFPKM[,samples]
#indata <- sweep(indata,2, apply(indata,2,median,na.rm=T))
#indata <- sweep(indata,1, apply(indata,1,median,na.rm=T))
#map.res.path <- file.path(tumor.path, paste("commonCluster_mediancenteredCombatLogFPKM_ClusterNum",N.cluster,sep = ""))
#featType <- "mediancenteredCombatLogFPKM"
#5
#indata <- combat.logcountsNorm[,samples]
#map.res.path <- file.path(tumor.path, paste("commonCluster_CombatLogcountsNorm_ClusterNum",N.cluster,sep = ""))
#featType <- "CombatLogcountsNorm"
#6 2 batch
#indata <- log2(nneg.combat.countsNorm[,samples] + 1)
#map.res.path <- file.path(tumor.path, paste("commonCluster_LognnegCombatcountsNorm_batch2_ClusterNum",N.cluster,sep = ""))
#featType <- "LognnegCombatcountsNorm_batch2"

#7 test 3 batch
#indata <- log2(nneg.combat.countsNorm[,samples] + 1)
#map.res.path <- file.path(tumor.path, paste("commonCluster_LognnegCombatcountsNorm_batch3_ClusterNum",N.cluster,sep = ""))
#featType <- "LognnegCombatcountsNorm_batch3"

#8
#indata <- as.data.frame(t(scale(t(norqua.combat.countsNorm))))
#map.res.path <- file.path(tumor.path, paste("commonCluster_scalednorquacombatcountsNorm_batch3_ClusterNum",N.cluster,sep = ""))
#featType <- "scalednorquaCombatcountsNorm_batch3"

#test consensus HC
N.cluster <- 5
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))

cluster5.UTUCBLCA.ans <- plot.common.cluster(indata, 
                                             tumorname=tumorname, 
                                             N.cluster=N.cluster, 
                                             N.bootstrap=N.bootstrap, 
                                             N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                             N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                             map.res.path=map.res.path, fig.path=fig.path, 
                                             featType=featType,
                                             annCol=annCol, annColors=annColors, 
                                             seed=123456, dist0="euclidean", dist="pearson", 
                                             clstCol=c(yellow, green, blue, red, purple), 
                                             namecvt=NULL, height = 7, fontsize=6, labRow = F, labCol = T,dendsort = T, cexCol = 0.06)

annCol$cluster = cluster5.UTUCBLCA.ans$group[rownames(annCol)]
table(annCol$Source,annCol$cluster)
fisher.test(matrix(c(89+42+76+74,126,0+1+0+0,19),byrow = T,ncol = 2)) # 7.57e-09
# [,1] [,2]
# [1,]  281  126
# [2,]    1   19

table(annCol$cluster,annCol$TCGA_Subtype)
fisher.test(matrix(c(88+41,1+1,7+5,3+4+43+23+73+24+74),byrow = T,ncol = 2),workspace = 1e9)

#log2(nnen.combat.countsFPKM) can derive results making sense

#########################################
# choose nneg.combat.countsNorm for NMF #
#########################################
N.cluster <- 5
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))

indata <- log2(nneg.combat.countsNorm[,samples] + 1)
map.res.path <- file.path(tumor.path, paste("NMFCluster_LognnegCombatcountsNorm_batch3_ClusterNum",N.cluster,sep = ""))
featType <- "NMF_LognnegCombatcountsNorm_batch3"


NMF5.mRNA.tumor.ans <- plot.common.cluster.NMF.addbar(indata,
                                                      tumorname=tumorname,
                                                      N.cluster=N.cluster,
                                                      N.bootstrap=N.bootstrap,
                                                      N.gene.per.bootstrap=N.gene.per.bootstrap,
                                                      N.sample.per.bootstrap=N.sample.per.bootstrap,
                                                      dist0 = "pearson",link0 = "ward.D",
                                                      map.res.path,
                                                      fig.path,
                                                      featType=featType,
                                                      annCol=annCol[samples,],
                                                      annColors=annColors, clstCol = c(purple, red, yellow, green, blue),
                                                      height=7, fontsize=8, labRow = F, labCol = F,dendsort = T)

#########################################
# choose nneg.combat.FPKM for NMF #
#########################################
N.cluster <- 5
N.bootstrap <- 500
N.gene.per.bootstrap <- round(0.8*nrow(indata))
N.sample.per.bootstrap <- round(0.8*ncol(indata))

indata <- log2(nneg.combat.FPKM[,samples] + 1)
map.res.path <- file.path(tumor.path, paste("NMFCluster_LognnegCombatFPKM_batch3_ClusterNum",N.cluster,sep = ""))
featType <- "NMF_LognnegCombatFPKM_batch3"


NMF5.mRNA.tumor.ans <- plot.common.cluster.NMF.addbar(indata,
                                                      tumorname=tumorname,
                                                      N.cluster=N.cluster,
                                                      N.bootstrap=N.bootstrap,
                                                      N.gene.per.bootstrap=N.gene.per.bootstrap,
                                                      N.sample.per.bootstrap=N.sample.per.bootstrap,
                                                      dist0 = "pearson",link0 = "ward.D",
                                                      map.res.path,
                                                      fig.path,
                                                      featType=featType,
                                                      annCol=annCol[samples,],
                                                      annColors=annColors, clstCol = c(purple, red, yellow, green, blue),
                                                      height=7, fontsize=8, labRow = F, labCol = F,dendsort = T)


#2.	Perform supervised clustering (K=5) by 50 genes (ASCO Poster) with combined TCGA Bladder RNASeq data (tumor only) with UTUC UTUC RNAseq data with bar of TCGA_Subtype

ASCO.genes <- read.table(file.path(data.path,"ASCOgenes.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)

combat.FPKM.HUGO <- as.data.frame(combat.FPKM)
combat.FPKM.HUGO$genename <- Ginfo[rownames(combat.FPKM.HUGO),"genename"] 
combat.FPKM.HUGO <- combat.FPKM.HUGO[!duplicated(combat.FPKM.HUGO$genename),]; rownames(combat.FPKM.HUGO) <- combat.FPKM.HUGO$genename; combat.FPKM.HUGO <- combat.FPKM.HUGO[,-ncol(combat.FPKM.HUGO)]

combat.logFPKM.HUGO <- as.data.frame(combat.logFPKM)
combat.logFPKM.HUGO$genename <- Ginfo[rownames(combat.logFPKM.HUGO),"genename"] 
combat.logFPKM.HUGO <- combat.logFPKM.HUGO[!duplicated(combat.logFPKM.HUGO$genename),]; rownames(combat.logFPKM.HUGO) <- combat.logFPKM.HUGO$genename; combat.logFPKM.HUGO <- combat.logFPKM.HUGO[,-ncol(combat.logFPKM.HUGO)]

nneg.combat.FPKM.HUGO <- as.data.frame(nneg.combat.FPKM)
nneg.combat.FPKM.HUGO$genename <- Ginfo[rownames(nneg.combat.FPKM.HUGO),"genename"] 
nneg.combat.FPKM.HUGO <- nneg.combat.FPKM.HUGO[!duplicated(nneg.combat.FPKM.HUGO$genename),]; rownames(nneg.combat.FPKM.HUGO) <- nneg.combat.FPKM.HUGO$genename; nneg.combat.FPKM.HUGO <- nneg.combat.FPKM.HUGO[,-ncol(nneg.combat.FPKM.HUGO)]

nneg.combat.logFPKM.HUGO <- as.data.frame(nneg.combat.logFPKM)
nneg.combat.logFPKM.HUGO$genename <- Ginfo[rownames(nneg.combat.logFPKM.HUGO),"genename"] 
nneg.combat.logFPKM.HUGO <- nneg.combat.logFPKM.HUGO[!duplicated(nneg.combat.logFPKM.HUGO$genename),]; rownames(nneg.combat.logFPKM.HUGO) <- nneg.combat.logFPKM.HUGO$genename; nneg.combat.logFPKM.HUGO <- nneg.combat.logFPKM.HUGO[,-ncol(nneg.combat.logFPKM.HUGO)]

nneg.combat.countsNorm.HUGO <- as.data.frame(nneg.combat.countsNorm)
nneg.combat.countsNorm.HUGO$genename <- Ginfo[rownames(nneg.combat.countsNorm.HUGO),"genename"] 
nneg.combat.countsNorm.HUGO <- nneg.combat.countsNorm.HUGO[!duplicated(nneg.combat.countsNorm.HUGO$genename),]; rownames(nneg.combat.countsNorm.HUGO) <- nneg.combat.countsNorm.HUGO$genename; nneg.combat.countsNorm.HUGO <- nneg.combat.countsNorm.HUGO[,-ncol(nneg.combat.countsNorm.HUGO)]

nneg.combat.logcountsNorm.HUGO <- as.data.frame(nneg.combat.logcountsNorm)
nneg.combat.logcountsNorm.HUGO$genename <- Ginfo[rownames(nneg.combat.logcountsNorm.HUGO),"genename"] 
nneg.combat.logcountsNorm.HUGO <- nneg.combat.logcountsNorm.HUGO[!duplicated(nneg.combat.logcountsNorm.HUGO$genename),]; rownames(nneg.combat.logcountsNorm.HUGO) <- nneg.combat.logcountsNorm.HUGO$genename; nneg.combat.logcountsNorm.HUGO <- nneg.combat.logcountsNorm.HUGO[,-ncol(nneg.combat.logcountsNorm.HUGO)]


rowids <- intersect(rownames(nneg.combat.countsNorm.HUGO),ASCO.genes$`GENE LIST`)

tmp <- log2(nneg.combat.countsNorm.HUGO[rowids,samples] + 1)
hcg <- hclust(distanceMatrix(as.matrix(t(tmp)), "pearson"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(tmp), "pearson"), "ward.D")
plotdata <- as.data.frame(na.omit(standarize.fun(indata = tmp,halfwidth = 3)))

outFigFile <- "BLCA_UTUC_lognnegCombatCountsNorm_supervised_ASCOgenes.pdf"
pdf(file.path(fig.path, outFigFile), height=10)
hv <- aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol, annColors=annColors, color=greenred(128), revC=TRUE, fontsize=4, cexCol = 0.5, cexRow = 0.8,cexAnn = 0.5)
invisible(dev.off())

res <- nmf(nneg.combat.logFPKM.HUGO[rowids,samples], rank=5, nrun=50, seed=1248103)
consensusmap(res, annCol = annCol, annColors = annColors)

#3.	Perform unsupervised consensus clustering (K=2) by mRNA with UTUC RNAseq data alone with bar of TCGA_Subtype

samples.UTUC <- rownames(annCol[which(annCol$Source=="UTUC"),])

res <-normData(countsTable[,samples.UTUC], tailrows, res.path, "UTUC")
UTUC.countsNorm <- res$countsNorm; vsd <- res$vsd

UTUC.FPKM <- calFPKM( countsTable[,samples.UTUC], Vids, tailrows, Ginfo, outfile=file.path(res.path,  "UTUC.FPKM.all.features.txt"))

PASSFlag.lnc.UTUC <- LowExpFilter(UTUC.FPKM, UTUC.countsNorm, lowcut.lnc, lowpct, Mids, Lids, Ginfo, res.path, paste("UTUC.lowcut",lowcut.lnc,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste("UTUC.lowcut",lowcut.lnc,".lowpct",lowpct,".Count.mRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste("UTUC.lowcut",lowcut.lnc,".lowpct",lowpct,".FPKM.mRNA.LowExpfiltered.txt",sep = "")))

PASSFlag.mRNA.UTUC <- LowExpFilter(UTUC.FPKM, UTUC.countsNorm, lowcut.mRNA, lowpct, Mids, Lids, Ginfo, res.path, paste("UTUC.lowcut",lowcut.mRNA,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste("UTUC.lowcut",lowcut.mRNA,".lowpct",lowpct,".Count.lncRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste("UTUC.lowcut",lowcut.mRNA,".lowpct",lowpct,".FPKM.lncRNA.LowExpfiltered.txt",sep = "")))

#remove batch effect
rowids <- intersect(Mids, names(PASSFlag.mRNA.UTUC[PASSFlag.mRNA.UTUC==TRUE]))
Sinfo$batch <- rep(c("BLCA","UTUCrna","UTUCr"),times=c(426,14,6))
modcombat = model.matrix(~1, data=Sinfo[c(427:446),])

indata <- UTUC.FPKM[rowids,]
combat.UTUC.FPKM = ComBat(dat=as.matrix(indata), batch=Sinfo[c(427:446),"batch"], mod=modcombat)
nneg.combat.UTUC.FPKM <- pmax(combat.UTUC.FPKM,0)*2
tmp <- as.data.frame(t(scale(t(UTUC.FPKM[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.FPKM.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)
tmp <- as.data.frame(t(scale(t(combat.UTUC.FPKM[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.combatFPKM.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)
tmp <- as.data.frame(t(scale(t(nneg.combat.UTUC.FPKM[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.nnegcombatFPKM.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)

indata <- log2(UTUC.FPKM[rowids,] + 1)
combat.UTUC.logFPKM = ComBat(dat=as.matrix(indata), batch=Sinfo[c(427:446),"batch"], mod=modcombat)

indata <- UTUC.countsNorm[rowids,]
combat.UTUC.countsNorm = ComBat(dat=as.matrix(indata), batch=Sinfo[c(427:446),"batch"], mod=modcombat)
nneg.combat.UTUC.countsNorm <- pmax(combat.UTUC.countsNorm,0)*2
tmp <- as.data.frame(t(scale(t(UTUC.countsNorm[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.countsNorm.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)
tmp <- as.data.frame(t(scale(t(combat.UTUC.countsNorm[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.combatcountsNorm.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)
tmp <- as.data.frame(t(scale(t(nneg.combat.UTUC.countsNorm[rowids,]))))
batchPCA(indata = tmp,batch = rep(c("UTUCrna","UTUCr"),times=c(14,6)),fig.dir = fig.path,PCA.fig.title = "PCA.UTUC.nnegcombatcountsNorm.batch",cols = c(yellow,green),showID = F,cex = 0.7,showLegend = T)

indata <- log2(UTUC.countsNorm[rowids,] + 1)
combat.UTUC.logcountsNorm = ComBat(dat=as.matrix(indata), batch=Sinfo[c(427:446),"batch"], mod=modcombat)

#indata <- log2(nneg.combat.UTUC.FPKM[,samples.UTUC] + 1)
#indata <- log2(nneg.combat.countsNorm[,samples.UTUC] + 1)
#indata <- combat.UTUC.logcountsNorm
#indata <- log2(nneg.combat.UTUC.FPKM[rowids,] + 1)

indata <- log2(nneg.combat.UTUC.countsNorm[rowids,] + 1)

tmp <- apply(indata, 1, sd)
tmp <- names(tmp[tmp>quantile(tmp,probs=seq(0,1,0.25))[4]])

plotdata <- as.data.frame(na.omit(standarize.fun(indata = indata[tmp,],halfwidth = 3)))
#hcg <- hclust(distanceMatrix(as.matrix(t(indata[tmp,])), "pearson"), "ward.D")
#hcs <- hclust(distanceMatrix(as.matrix(indata[tmp,]), "pearson"), "ward.D")

hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "ward.D2")

outFigFile <- "UTUC_unsupervised_quantile0.9_test_lognnegcombatUTUCcountsNorm_euclidean_ward.D.pdf"
pdf(file.path(fig.path, outFigFile), height=6)
hv <- aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol[samples.UTUC,], annColors=annColors, color=greenred(128), revC=TRUE, fontsize=4, cexCol = 0.5, cexRow = 0.8,cexAnn = 0.5,labRow = NA)
invisible(dev.off())

indata <- plotdata[tmp,]
N.cluster <- 4
N.bootstrap <- 500
N.gene.per.bootstrap <- round(1*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))
map.res.path <- file.path(tumor.path, paste("commonCluster_quantile0.75_LognnegCombatUTUCcountsNorm_ClusterNum",N.cluster,sep = ""))
featType <- "UTUC_quantile0.75LognnegCombatcountsNorm"

cluster4.UTUC.ans <- plot.common.cluster(indata[tmp,], 
                                         tumorname="UTUC", 
                                         N.cluster=N.cluster, 
                                         N.bootstrap=N.bootstrap, 
                                         N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                         N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                         map.res.path=map.res.path, fig.path=fig.path, 
                                         featType=featType,
                                         annCol=annCol[samples.UTUC,], annColors=annColors, 
                                         seed=123456, dist0="pearson", dist="pearson", link0 = "ward.D2",link = "average",
                                         clstCol=c(yellow, green, blue, red), 
                                         namecvt=NULL, height = 6, fontsize=6, labRow = F, labCol = F,dendsort = T)

indata <- log2(nneg.combat.UTUC.countsNorm[tmp,] + 1)
N.cluster <- 4
N.bootstrap <- 500
N.gene.per.bootstrap <- round(1*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))
map.res.path <- file.path(tumor.path, paste("NMFCluster_quantile0.75_LognnegCombatUTUCcountsNorm_ClusterNum",N.cluster,sep = ""))
featType <- "NMF_UTUC_quantile0.75LognnegCombatcountsNorm"
cluster4.UTUC.nmf.ans <- plot.common.cluster.NMF.addbar(indata, 
                                                          tumorname="UTUC", 
                                                          N.cluster=N.cluster, 
                                                          N.bootstrap=N.bootstrap, 
                                                          N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                          N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                          dist0 = "pearson",link0 = "ward.D",
                                                          map.res.path, 
                                                          fig.path, 
                                                          featType=featType, 
                                                          annCol=annCol[samples.UTUC,], 
                                                          annColors=annColors, clstCol = c( red, yellow, green, blue),
                                                          height=6, fontsize=6, labRow = F, labCol = F,dendsort = T) 

#UTUCrna only
samples.UTUCrna <- rownames(annCol[which(annCol$TCGA_Subtype %in% c("UTUC_NMI","UTUC_MI")),])

res <-normData(countsTable[,samples.UTUCrna], tailrows, res.path, "UTUCrna")
UTUCrna.countsNorm <- res$countsNorm; vsd <- res$vsd

UTUCrna.FPKM <- calFPKM( countsTable[,samples.UTUCrna], Vids, tailrows, Ginfo, outfile=file.path(res.path,  "UTUCrna.FPKM.all.features.txt"))

PASSFlag.lnc.UTUCrna <- LowExpFilter(UTUCrna.FPKM, UTUCrna.countsNorm, lowcut.lnc, lowpct, Mids, Lids, Ginfo, res.path, paste("UTUCrna.lowcut",lowcut.lnc,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste("UTUCrna.lowcut",lowcut.lnc,".lowpct",lowpct,".Count.mRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste("UTUCrna.lowcut",lowcut.lnc,".lowpct",lowpct,".FPKM.mRNA.LowExpfiltered.txt",sep = "")))

PASSFlag.mRNA.UTUCrna <- LowExpFilter(UTUCrna.FPKM, UTUCrna.countsNorm, lowcut.mRNA, lowpct, Mids, Lids, Ginfo, res.path, paste("UTUCrna.lowcut",lowcut.mRNA,".lowpct",lowpct,sep = ""))
file.remove(file = file.path(res.path,paste("UTUCrna.lowcut",lowcut.mRNA,".lowpct",lowpct,".Count.lncRNA.LowExpfiltered.txt",sep = "")))
file.remove(file = file.path(res.path,paste("UTUCrna.lowcut",lowcut.mRNA,".lowpct",lowpct,".FPKM.lncRNA.LowExpfiltered.txt",sep = "")))

rowids <- intersect(Mids, names(PASSFlag.mRNA.UTUCrna[PASSFlag.mRNA.UTUCrna==TRUE]))
indata <- log2(UTUCrna.countsNorm[rowids,] + 1)
tmp <- apply(indata, 1, sd)
tmp <- names(tmp[tmp>quantile(tmp,probs=seq(0,1,0.25))[4]])
indata <- indata[tmp,]

N.cluster <- 3
N.bootstrap <- 500
N.gene.per.bootstrap <- round(1*nrow(indata))
N.sample.per.bootstrap <- round(0.9*ncol(indata))

map.res.path <- file.path(tumor.path, paste("commonCluster_quantile0.75_logcountsNorm_UTUCrna_ClusterNum",N.cluster,sep = ""))
featType <- "UTUCrna_quantile0.75LogcountsNorm"
cluster3.UTUCrna.ans <- plot.common.cluster(indata, 
                                            tumorname="UTUCrna", 
                                            N.cluster=N.cluster, 
                                            N.bootstrap=N.bootstrap, 
                                            N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                            N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                            map.res.path=map.res.path, fig.path=fig.path, 
                                            featType=featType,
                                            annCol=annCol[samples.UTUCrna,], annColors=annColors, 
                                            seed=123456, dist0="euclidean", dist="pearson", link0 = "ward.D",link = "ward.D",
                                            clstCol=c(yellow, green, blue, red), 
                                            namecvt=NULL, height = 6, fontsize=6, labRow = F, labCol = F,dendsort = T)


##########################################################
## Differentially expression analysis between C1 and C2 ##
group <- read.csv(file.path(res.path,"Sinfo of Tumor C1 C2 derived from nofiltering for ChAMP.csv"),row.names = NULL,header = T,stringsAsFactors = F)
rownames(group) <- group$Sample_Name
samples.rna <- colnames(UTUC.countsNorm)
index <- sapply(strsplit(samples.rna,"_"), "[",2) #attention for "nT" should be "0nT"
index <- c("11T", "13T", "19T", "01T",  "21T", "22T", "25T", "26T", "29T", "02T",  "31T", "04T",  "08T",  "09T",  "10T", "16T", "37T", "03T",  "05T",  "07T")
samples.dna <- rownames(group)
samples <- intersect(index,samples.dna)
group <- group[samples,]
tmp <- group[samples,"Sample_Group"]
samples <- c("UTUCrna_13T","UTUCrna_19T","UTUCrna_21T","UTUCrna_22T","UTUCrna_25T","UTUCrna_26T","UTUCrna_29T","UTUCrna_2T",
             "UTUCrna_31T","UTUCrna_4T","UTUCrna_8T","UTUCrna_9T","UTUCr_16T","UTUCr_37T","UTUCr_3T","UTUCr_5T")
names(tmp) <- samples 
group$Sample_LName <- samples
#group <- group[order(group$Sample_Group),]
group$batch <- "UTUCr"
group$batch <- ifelse(grepl("rna",group$Sample_LName)>0,"UTUCrna","UTUCr")
rownames(group) <- samples

complist <- createList.C1vsC2(group=tmp)
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
#twoclasscomp(res.path, Ginfo, countsTable=countsTable[, samples], tailrows, features=Mids, featType="mRNA", complist, PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, samples], tailrows, Groupinfo=group, features=Mids, featType="mRNA_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

tmp <- as.data.frame(t(scale(t(UTUC.countsNorm[rowids,samples]))))
colnames(tmp) <- intersect(index,samples.dna)
batchPCA(indata = tmp[,rownames(group)],batch = rep(c("C1","C2"),times=c(11,5)),fig.dir = fig.path,PCA.fig.title = "PCA.C1vsC2.batch",cols = c(blue,red),showID = F,cex = 0.7,showLegend = T)

#volcano plot
xceiling=10
yceiling=5
fccut=2
pcut=0.25
DEres <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
x <- DEres$logFC
y <- -log10(DEres$FDR)
y[(!is.na(y)) & y>=yceiling]  <-  yceiling

x[(!is.na(x)) & x>=xceiling]  <-  xceiling
x[(!is.na(x)) & x<=-xceiling] <- -xceiling

sigupindex <- which( (!is.na(y)) & y>=(-log10(pcut)) & x >= fccut)
sigdnindex <- which( (!is.na(y)) & y>=(-log10(pcut)) & x < -fccut)

naindex <- which(is.na(y))
unsigindex <- setdiff( 1:length(y), c(naindex, sigupindex,sigdnindex) )
pdf(file.path(fig.path,"volcano plot for mRNA_rmBatch_edgeR_test_result.C1_vs_C2.pdf"))
plot(x[unsigindex], y[unsigindex], xlim=c(-max(x,na.rm=T), max(x, na.rm=T)), ylim=c(0, max(y,na.rm=TRUE)), type="p", col="black", pch=16, main="", xlab="Gene Expression log2(Fold Change)", ylab="-log10(FDR)", lwd=1, cex=0.2)
points(x[sigupindex], y[sigupindex], type="p", col=red, pch=20, lwd=1, cex=0.5)
points(x[sigdnindex], y[sigdnindex], type="p", col=blue, pch=20, lwd=1, cex=0.5)
points(x=c(-100, 100), y=c(-log10(pcut), -log10(pcut)), type='l', col="grey", lty="dotted", lwd=1.3)
points(x=c(0, 0), y=c(-100, 100), type='l', col="grey", lty="dotted", lwd=1.3)
invisible(dev.off())

rm(x);rm(y);rm(xceiling);rm(yceiling);rm(fccut);rm(pcut)
#Using edgeR

#attention batch=c(rep("C2",5),rep("C1",11))  means C1 - C2
# group <- group[order(group$Sample_Group,decreasing = T),]
# batch=factor(c(rep("C2",5),rep("C1",11)))    
# batch2=factor(group$batch)
# 
# y <- DGEList(counts=countsTable[rowids, group$Sample_LName],group=batch)
# y <- calcNormFactors(y)
# 
# design <- model.matrix(~batch2+batch)
# rownames(design) <- colnames(y)
# 
# y <- estimateCommonDisp(y)
# y <- estimateTagwiseDisp(y)
# et <- exactTest(y,pair = c("C2","C1"))
# topTags(et)
# ordered_tags <- topTags(et, n=100000)
# 
# allDiff=ordered_tags$table
# allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
# diff=allDiff
# diff$GeneSymbol <- Ginfo[rownames(diff),"genename"]
# newData=y$pseudo.counts
# write.table(diff,file.path(res.path,"mRNA_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",row.names = T,col.names = NA)

##########################################################
## Differentially expression analysis between T3 and T24 ##
group <- read.csv(file.path(res.path,"Sinfo of detailed Tumor Normal.csv"),row.names = NULL,header = T,stringsAsFactors = F)
rownames(group) <- group$Sample_Name
samples.rna <- colnames(UTUC.countsNorm)
index <- sapply(strsplit(samples.rna,"_"), "[",2) #attention for "nT" should be "0nT"
index <- c("11T", "13T", "19T", "01T",  "21T", "22T", "25T", "26T", "29T", "02T",  "31T", "04T",  "08T",  "09T",  "10T", "16T", "37T", "03T",  "05T",  "07T")
samples.dna <- rownames(group)
samples <- intersect(index,samples.dna)
group <- group[samples,]
tmp <- group[samples,"Sample_Group"]
samples <- c("UTUCrna_13T","UTUCrna_19T","UTUCrna_21T","UTUCrna_22T","UTUCrna_25T","UTUCrna_26T","UTUCrna_29T","UTUCrna_2T",
             "UTUCrna_31T","UTUCrna_4T","UTUCrna_8T","UTUCrna_9T","UTUCr_16T","UTUCr_37T","UTUCr_3T","UTUCr_5T")
names(tmp) <- samples 
group$Sample_LName <- samples
group$Sample_Group2 <- group$Sample_Group
group[which(group$Sample_Group2 %in% c("Tumor2","Tumor4")),"Sample_Group2"] <- "T24"
group[which(group$Sample_Group2 %in% c("Tumor3")),"Sample_Group2"] <- "T3"
#group <- group[order(group$Sample_Group2),]
group$batch <- "UTUCr"
group$batch <- ifelse(grepl("rna",group$Sample_LName)>0,"UTUCrna","UTUCr")

complist <- createList.T3vsT24(group=tmp)
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
#twoclasscomp(res.path, Ginfo, countsTable=countsTable[, samples], tailrows, features=Mids, featType="mRNA", complist, PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, samples], tailrows, Groupinfo=group, features=Mids, featType="mRNA_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

#Using edgeR
#attention batch=c(rep("T24",5),rep("T3",11))  means T3 - T24
# group <- group[order(group$Sample_Group2),]
# batch=c(rep("T24",11),rep("T3",5)) 
# batch2=factor(group$batch)
# 
# y <- DGEList(counts=countsTable[rowids, group$Sample_LName],group=batch)
# y <- calcNormFactors(y)
# 
# design <- model.matrix(~batch2+batch)
# rownames(design) <- colnames(y)
# 
# y <- estimateCommonDisp(y)
# y <- estimateTagwiseDisp(y)
# et <- exactTest(y,pair = c("T24","T3"))
# topTags(et)
# ordered_tags <- topTags(et, n=100000)
# 
# allDiff=ordered_tags$table
# allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
# diff=allDiff
# diff$GeneSymbol <- Ginfo[rownames(diff),"genename"]
# newData=y$pseudo.counts
# write.table(diff,file.path(res.path,"mRNA_edgeR_test_result.T3_vs_T24.txt"),sep = "\t",row.names = T,col.names = NA)

###############################
# perform GSEA for edgeR DEGs #
tmp <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
#tmp$GeneSymbol <- toupper(tmp$GeneSymbol)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.DEGs.C1vsC2 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res <- data.frame(GSEA.DEGs.C1vsC2)
write.table(as.data.frame(GSEA.DEGs.C1vsC2),file.path(res.path,"GSEA results for C1vsC2 DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

###########################
# inverse GSEA for C2vsC1 #
tmp <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- -geneList # inverse logFC value
geneList <- sort(geneList,decreasing = T)
GSEA.DEGs.C2vsC1 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res2 <- data.frame(GSEA.DEGs.C2vsC1)
write.table(as.data.frame(GSEA.DEGs.C2vsC1),file.path(res.path,"Inverse GSEA results for C1vsC2 DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

eg <- bitr(names(geneList), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneList <- geneList[eg$SYMBOL] 
names(geneList) <- eg$ENTREZID
GSEAKEGG.DEGs.C2vsC1 <- gseKEGG(geneList = geneList,organism = 'hsa',seed = T,pvalueCutoff = 0.25)
res3 <- as.data.frame(GSEAKEGG.DEGs.C2vsC1)
GSEAMKEGG.DEGs.C2vsC1 <- gseMKEGG(geneList = geneList,organism = 'hsa',seed = T,pvalueCutoff = 0.25)
res4 <- as.data.frame(GSEAMKEGG.DEGs.C2vsC1)

library("pathview")
#detach dplyr and succeed!
hsa04060 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04060",
                     species    = "hsa",
                     #gene.idtype = "ENTREZID",
                     #gene.annotpkg = "org.Hs.eg.db",
                     kegg.dir = fig.path,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
################################
# Heatmap for H3K27ME3 pathway #

Gaby.pathway <- res[grep("H3K27ME3",res$Description),]
rownames(Gaby.pathway)
#"BENPORATH_ES_WITH_H3K27ME3"                 "MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"
Gaby.pathway.gene <- list()
# Gaby.pathway.gene[["BENPORATH_ES_WITH_H3K27ME3"]] <- strsplit(Gaby.pathway[1,"core_enrichment"],"/")[[1]]
# Gaby.pathway.gene[["MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"]] <- strsplit(Gaby.pathway[2,"core_enrichment"],"/")[[1]]
Gaby.pathway.gene[["BENPORATH_ES_WITH_H3K27ME3"]] <- toupper(MSigDB[which(MSigDB$ont %in% "BENPORATH_ES_WITH_H3K27ME3"),"gene"])
Gaby.pathway.gene[["MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"]] <- toupper(MSigDB[which(MSigDB$ont %in% "MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3"),"gene"])

indata <- combat.UTUC.FPKM[intersect(rownames(Ginfo[which(Ginfo$genename %in% Gaby.pathway.gene[[1]]),]),rownames(combat.UTUC.FPKM)),rownames(annCol.C1vsC2.rna)]
rownames(indata) <- Ginfo[rownames(indata),"genename"]
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward.D")
plotdata <- standarize.fun(indata=indata, halfwidth=1, centerFlag=T, scaleFlag=T)

pdf(file.path(fig.path,"BENPORATH_ES_WITH_H3K27ME3_withoutColv.pdf"),height=7)
hv = aheatmap(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=NA, annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=blueyellow(64), revC=TRUE, fontsize=5,cexRow = 0.3,cexCol = 1.5)
invisible(dev.off())

pdf(file.path(fig.path,"BENPORATH_ES_WITH_H3K27ME3_withColv.pdf"),height=7)
hv = aheatmap(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=as.dendrogram(hcs), annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=blueyellow(64), revC=TRUE, fontsize=5,cexRow = 0.3,cexCol = 1.5)
invisible(dev.off())

indata <- combat.UTUC.FPKM[intersect(rownames(Ginfo[which(Ginfo$genename %in% Gaby.pathway.gene[[2]]),]),rownames(combat.UTUC.FPKM)),rownames(annCol.C1vsC2.rna)]
rownames(indata) <- Ginfo[rownames(indata),"genename"]
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward.D")
plotdata <- standarize.fun(indata=indata, halfwidth=1, centerFlag=T, scaleFlag=T)

pdf(file.path(fig.path,"MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3_withoutColv.pdf"),height=7)
hv = aheatmap(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=NA, annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=blueyellow(64), revC=TRUE, fontsize=5,cexRow = 1,cexCol = 1.5)
invisible(dev.off())

pdf(file.path(fig.path,"MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3_withColv.pdf"),height=7)
hv = aheatmap(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=as.dendrogram(hcs), annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=blueyellow(64), revC=TRUE, fontsize=5,cexRow = 1,cexCol = 1.5)
invisible(dev.off())

g1 <- gseaplot2(GSEA.DEGs.C2vsC1,"BENPORATH_ES_WITH_H3K27ME3",color = darkgreen,pvalue_table = T)
ggsave(filename = file.path(fig.path,"BENPORATH_ES_WITH_H3K27ME3_gseaplot.pdf"),width = 10,height = 7)
g2 <- gseaplot2(GSEA.DEGs.C2vsC1,"MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3",color = darkgreen,pvalue_table = T)
ggsave(filename = file.path(fig.path,"MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3_gseaplot.pdf"),width = 12,height = 8)

#
label <- c("MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3",
           "MIKKELSEN_MCV6_HCP_WITH_H3K27ME3",
           "HATADA_METHYLATED_IN_LUNG_CANCER_UP",
           "MODULE_220",
           "YCATTAA_UNKNOWN",
           "HOQUE_METHYLATED_IN_CANCER",
           "MEISSNER_NPC_HCP_WITH_H3K4ME2",
           "MODULE_176",
           "WONG_ENDMETRIUM_CANCER_DN",
           "REACTOME_GPCR_LIGAND_BINDING",
           "MODULE_112",
           "SCHLESINGER_METHYLATED_DE_NOVO_IN_CANCER")
for (i in label) {
  gseaplot(GSEA.DEGs.C1vsC2,geneSetID = i,title = i)
  ggsave(filename = file.path(fig.path,paste0("GSEA_",i,".pdf")))
}

####################################################
# calculate ssGSEA foe 293 Gaby interested pathway #
tmp <- as.data.frame(combat.UTUC.FPKM)
tmp$genename <- Ginfo[rownames(tmp),"genename"] 
tmp <- tmp[!duplicated(tmp$genename),]; rownames(tmp) <- tmp$genename; tmp <- tmp[,-ncol(tmp)]
combat.UTUC.FPKM.HUGO <- tmp

gaby.interested.pathway <- read.table(file.path(data.path,"Gaby_interested_pathway.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
gaby.interested.pathway.list <- list()
for (i in 1:nrow(gaby.interested.pathway)) {
  tmp <- gaby.interested.pathway$PathwayID[i]
  gaby.interested.pathway.list[[tmp]] <- intersect(rownames(combat.UTUC.FPKM.HUGO),toupper(MSigDB[which(MSigDB$ont %in% tmp),"gene"]))
}
gaby.enrichscore <- gsva(as.matrix(combat.UTUC.FPKM.HUGO),gaby.interested.pathway.list,method="ssgsea")
gaby.enrichscore <- as.data.frame(gaby.enrichscore[,rownames(annCol.C1vsC2.rna)])

p <- vector()
for (i in 1:nrow(gaby.enrichscore)) {
  tmp <- round(wilcox.test(as.numeric(gaby.enrichscore[i,1:11]),as.numeric(gaby.enrichscore[i,12:15]),"greater")$p.value,4)
  p <- append(p, tmp)
}
rownames(gaby.enrichscore) <- paste0(rownames(gaby.enrichscore),"_p_",p)

plotscore <- as.data.frame(na.omit(standarize.fun(indata = gaby.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotscore)), "euclidean"), "ward.D")
pdf(file.path(fig.path, "UTUC_heatmap_mRNA_enrichmentscore293_for annotation.pdf"), height=12)
hv <- aheatmap(as.matrix(plotscore), Rowv=dendsort(as.dendrogram(hcg)), Colv=NA, annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=5, cexCol = 0.8, cexRow = 0.4,cexAnn = 0.6)
invisible(dev.off())
pdf(file.path(fig.path, "UTUC_heatmap_mRNA_enrichmentscore293.pdf"), height=6)
hv <- aheatmap(as.matrix(plotscore), Rowv=dendsort(as.dendrogram(hcg)), Colv=NA, annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=5, cexCol = 0.8, cexRow = 0.4,cexAnn = 0.6,labRow = NA)
invisible(dev.off())

########################################################
## Annotation MCPcounter to RNAseq analysis of C1vsC2 ##

# calculate MCPcounter (use orginal FPKM then remove combat)
tmp <- as.data.frame(UTUC.FPKM)
tmp$genename <- Ginfo[rownames(tmp),"genename"] 
tmp <- tmp[!duplicated(tmp$genename),]; rownames(tmp) <- tmp$genename; tmp <- tmp[,-ncol(tmp)]
write.table(tmp,file = file.path(res.path,paste("UTUC","_FPKM_all_features_HUGO_symbols.txt",sep = "")),row.names=T, col.names=NA, sep="\t", quote=F)
MCPscore <- MCPcounter.estimate(expression = tmp,featuresType = "HUGO_symbols")

write.table(MCPscore,file = file.path(res.path,paste("UTUC","_FPKM_MCPscore_HUGO_symbols.txt",sep = "")),row.names=T, col.names=NA, sep="\t", quote=F)

filterCommonGenes(input.f=file.path(res.path, paste("UTUC","_FPKM_all_features_HUGO_symbols.txt", sep="")) , output.f=file.path(res.path,paste("UTUC","_FPKM_all_features_HUGO_symbol_ESTIMATE.txt",sep = "")), id="GeneSymbol")
estimateScore(file.path(res.path,paste("UTUC","_FPKM_all_features_HUGO_symbol_ESTIMATE.txt",sep = "")), file.path(res.path,paste("UTUC","_FPKM_all_features_estimate_score.txt",sep = "")), platform="affymetrix")
est.score <- read.table(file = file.path(res.path,paste("UTUC","_FPKM_all_features_estimate_score.txt",sep = "")),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.score) <- est.score[,2]; colnames(est.score) <- est.score[1,]; est.score <- est.score[-1,c(-1,-2)];

# create annCol and annColors
MCPscore.combat <- ComBat(dat=as.matrix(MCPscore), batch=Sinfo[c(427:446),"batch"], mod=modcombat)
MCPscore.combat <- annTrackScale(indata = MCPscore.combat, halfwidth = 2, poolsd = F); MCPscore.combat <- as.data.frame(t(MCPscore.combat))

est.score.combat <- ComBat(dat=as.matrix(sapply(est.score, as.numeric)), batch=Sinfo[c(427:446),"batch"], mod=modcombat)
rownames(est.score.combat) <- rownames(est.score)

tmp <- as.data.frame(t(est.score.combat)); tmp <- tmp[,setdiff(colnames(tmp),c("TumorPurity"))]; tmp <- as.data.frame(t(tmp),stringsAsFactors=F)
tmp <- as.data.frame(sapply(tmp, as.numeric)); rownames(tmp) <- setdiff(rownames(est.score.combat),"TumorPurity")
tmp <- annTrackScale(indata = tmp, halfwidth = 2, poolsd = F); tmp <- as.data.frame(t(tmp))

annCol.C1vsC2.rna <- UTUC.annotation
annCol.C1vsC2.rna <- annCol.C1vsC2.rna[which(annCol.C1vsC2.rna$DetailID %in% samples),]
annCol.C1vsC2.rna <- annCol.C1vsC2.rna[order(annCol.C1vsC2.rna$`Cluster DNA methylation no filter without normal`,decreasing = F),]
annCol.C1vsC2.rna <- cbind.data.frame(MCPscore.combat[annCol.C1vsC2.rna$DetailID,],tmp[annCol.C1vsC2.rna$DetailID,],annCol.C1vsC2.rna)

annColors.C1vsC2.rna <- list()
annColors.C1vsC2.rna[["T cells"]] <- annColors.C1vsC2.rna[["CD8 T cells"]] <- greenred(64)
annColors.C1vsC2.rna[["Cytotoxic lymphocytes"]] <- annColors.C1vsC2.rna[["NK cells"]] <-annColors.C1vsC2.rna[["B lineage"]] <- greenred(64)
annColors.C1vsC2.rna[["Monocytic lineage"]] <- annColors.C1vsC2.rna[["Myeloid dendritic cells"]] <- greenred(64)
annColors.C1vsC2.rna[["Neutrophils"]] <- annColors.C1vsC2.rna[["Endothelial cells"]] <- annColors.C1vsC2.rna[["Fibroblasts"]] <- greenred(64)
annColors.C1vsC2.rna[["Cluster DNA methylation no filter without normal"]] <- c("C1"=red,"C2"=blue)
annColors.C1vsC2.rna[["ImmuneScore"]] <- bluered(64)
annColors.C1vsC2.rna[["StromalScore"]] <- bluered(64)
# using DEGs between C1vsC2 to plot heatmap with MCPcounter annotation
tmp <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.C1_vs_C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
DEGs <- tmp[which(tmp$FDR < 0.25 & (abs(tmp$logFC) >2)),]
indata <- combat.UTUC.FPKM[intersect(rownames(Ginfo[which(Ginfo$genename %in% rownames(DEGs)),]),rownames(combat.UTUC.FPKM)),rownames(annCol.C1vsC2.rna)]
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D")
plotdata <- standarize.fun(indata=indata, halfwidth=1, centerFlag=T, scaleFlag=T)

pdf(file.path(fig.path,"UTUC_C1vsC2_DEG_heatmap_withMCPcounterandEstimate.pdf"),height=7)
hv = aheatmap(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=NA, annCol=annCol.C1vsC2.rna[,c(1:12,21)], annColors=annColors.C1vsC2.rna, color=blueyellow(64), revC=TRUE, fontsize=5,labRow = NA)
#hv = heatmap.2(as.matrix(plotdata), Rowv=as.dendrogram(hcg), Colv=NULL,dendrogram = "row", col=greenred(64), trace="none",srtRow = 30)
invisible(dev.off())

###################################################################################################################
##### Indepent test for recurrent exome mutation between Non-muscle invasive and muscle-invasive UTUC samples #####
#sig.mut.gene <- c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","REPIN1","ARID1A","TP53","GANAB","RARG","ARID1B") #REPIN1 and RARG had duplications and equals to 3 cases mutation
sig.mut.gene <- c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB","ARID1B")
exome.mut <- read.table("F:/Project/UTUC_BLCA/exome UTUC/Results/UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_wideformat.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
exome.mut <- exome.mut[,c(14,17:46)]; colnames(exome.mut)[1] <- "GeneName"
# ARID1A.mut <- read.table(file.path(data.path,"recurrent mutations.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# ARID1A.mut <- ARID1A.mut[which(ARID1A.mut$GeneName == "ARID1A"),c(2,8:37)]
# if(all(colnames(exome.mut)==colnames(ARID1A.mut))) {exome.mut <- rbind.data.frame(exome.mut,ARID1A.mut)}
exome.mut <-exome.mut[which(exome.mut$GeneName %in% sig.mut.gene),]

exome.mut <- aggregate(exome.mut[2:31],by=list(exome.mut$GeneName),FUN=sum)
rownames(exome.mut) <- exome.mut$Group.1
exome.mut <- as.matrix(exome.mut[,-1])
exome.mut <- as.data.frame(floor(pmin(exome.mut,1.1)))

index <- sapply(strsplit(colnames(exome.mut),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
colnames(exome.mut) <- index
exome.mut <- rbind.data.frame(exome.mut,ARID1AB=exome.mut["ARID1A",]+exome.mut["ARID1B",])

write.table(exome.mut,file.path(res.path,"significant exome mutation trimed.txt"),sep = "\t",row.names = T,col.names = NA)

# group for muscle-invasive and non-musle-invasive
group <- UTUC.annotation[colnames(exome.mut),"Type"]; names(group) <- colnames(exome.mut)

genelist <- rownames(exome.mut[setdiff(rownames(exome.mut),c("ARID1B","ARID1AB")),])
binarymut <- as.data.frame(matrix(0, nrow=length(group), ncol=length(genelist)))
rownames(binarymut) <- names(group)
colnames(binarymut) <- genelist
for (k in 1:length(genelist)) {
  res <- create.anntrack(samples=names(group), subtype=createMutSubtype(exome.mut[setdiff(rownames(exome.mut),c("ARID1B","ARID1AB")),], names(group), genelist[k]))
  binarymut[, genelist[k]] <- res$subtype
}

out <- matrix(0, nrow=length(genelist), ncol=length(unique(group))+1)
colnames(out) <- c(levels(factor(group)), "pvalue")
rownames(out) <- paste(genelist, "Mutated", sep="_")
for (k in 1:length(genelist)) {
  genek <- genelist[k]
  x <- group
  y <- binarymut[names(x), genek, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["Mutated", ]
  
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out <- as.data.frame(out)
out$FDR <- p.adjust(out[,"pvalue"],method = "BH")
write.table(out, file.path(res.path, "Independence test between significant exome mutation and muscle-invasive type.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

# group for C1 and C2
group <- UTUC.annotation[colnames(exome.mut),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut)
group <- group[-c(9,28)] #remove 23T and 07 T which no methylation cluster
genelist <- rownames(exome.mut[setdiff(rownames(exome.mut),c("ARID1A","ARID1B")),])
binarymut <- as.data.frame(matrix(0, nrow=length(group), ncol=length(genelist)))
rownames(binarymut) <- names(group)
colnames(binarymut) <- genelist
for (k in 1:length(genelist)) {
  res <- create.anntrack(samples=names(group), subtype=createMutSubtype(exome.mut[setdiff(rownames(exome.mut),c("ARID1A","ARID1B")),setdiff(colnames(exome.mut),c("23T","07T"))], names(group), genelist[k]))
  binarymut[, genelist[k]] <- res$subtype
}

out <- matrix(0, nrow=length(genelist), ncol=length(unique(group))+1)
colnames(out) <- c(levels(factor(group)), "pvalue")
rownames(out) <- paste(genelist, "Mutated", sep="_")
for (k in 1:length(genelist)) {
  genek <- genelist[k]
  x <- group
  y <- binarymut[names(x), genek, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["Mutated", ]
  
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out <- as.data.frame(out)
out$FDR <- p.adjust(out[,"pvalue"],method = "BH")
write.table(out, file.path(res.path, "Independence test between significant exome mutation and methylation C1 C2 type.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

####################################################################
############### oncoprint for trimed mutation data #################

#use fill to set color

# mygene <- exome.mut; mygene <- mygene[-which(rownames(mygene) %in% c("TTN","HYDIN")),] # remove bad genes
# mut.order <- rownames(mygene)[order(rowSums(mygene),decreasing = T)]
# mygene[mygene == 1] <- "mut"
# mygene[mygene == 0] <- ""
# 
# alter_fun = list(
#   background = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = white, col = NA)) #
#   },
#   mut = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = darkblue, col = "black")) # mut color
#   })
# 
# #mutation bar plot color
# col = c("mut" = blue)
# 
# tmp <- UTUC.annotation[colnames(mygene),c(9:11,14,15,17)]
# colnames(tmp) <- c("Gender","Type","Grade","Pathological_Stage","Lymph_Nodes","Localisation")
# my_annotation = HeatmapAnnotation(df = tmp,col = list(Gender = c("M" = blue, "F" = darkred),
#                                                      Type = c("Non-muscle invasive" = nake, "muscle-invasive" = soil),
#                                                      Grade = c("High" = sun, "Low" = lightblue),
#                                                      Pathological_Stage = c("pTa" = yellow, "pT1" = orange, "pT2" =lightred, "pT3" =peach, "pT4" =darkred),
#                                                      Lymph_Nodes = c("N0"=white,"N1"=cyan,"N2"=violet),
#                                                      Localisation = c("renal pelvis" = gold,"ureter"=cherry,"N/A"=grey) ))
# 
# p <- oncoPrint_plus(mygene, get_type = function(x) x,
#           alter_fun = alter_fun, col = col,
#           #row_order = NULL, # no order by mutation frequency
#           #column_order = NULL, #no order by mutation frequency of sample
#           show_pct = T, #show pct in left
#           show_row_barplot = T, #show bar plot in right
#           #remove_empty_columns = TRUE,
#           column_title = "OncoPrint for UTUC exome mutation",
#           bottom_annotation = my_annotation,
#           show_heatmap_legend=F) # hide heatmap legend
# 
# oncoprint.sample.order <- colnames(mygene)[p$sampleOrder]
# pdf(file.path(fig.path,"oncoprint_UTUC_exome_mutation2.pdf"),width = 10,height = 6)
# draw(p$oncoP, heatmap_legend_side = "bottom")
# dev.off()

# oncoprint with detailed mutation type
mygene <- read.table(file.path(res.path,"significant exome mutation trimed with mutation type.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
mygene <- mygene[setdiff(rownames(mygene),c("RARG","REPIN1","ARID1AB","ARID1B")),]
mygene <- mygene[rowSums(mygene!="0") > 1,] # just keep mutation that has >= 2 cases mutation
# index <- sapply(strsplit(colnames(mygene),"_"), "[",2)
# index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
# colnames(mygene) <- index
mygene[mygene == "nonsynonymous SNV"] <- "Nonsynonymous"
mygene[mygene == "synonymous SNV"] <- "Synonymous"
mygene[mygene == "stopgain"] <- "Stopgain"
mygene[mygene == "frameshift deletion"] <- "Frameshift_Deletion"
mygene[mygene == "frameshift insertion"] <- "Frameshift_Insertion"
mygene[mygene == "nonframeshift deletion"] <- "Nonframeshift_Deletion"
mygene[mygene == "nonframeshift insertion"] <- "Nonframeshift_Insertion"
#mygene[mygene == "multi-hit"] <- "Multi_Hit"
mygene[mygene == "splicing"] <- "Splicing"
mygene[mygene == "0"] <- ""

#mycol <- brewer.pal(8,"RdYlBu") # color with gradually change
mycol <- colorRampPalette(brewer.pal(8,'Spectral'))(8) #more like rainbow color
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = lightgrey, col = NA)) #不要背景色
  },
  Nonsynonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[1], col = "black")) 
  },
  Stopgain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[2], col = "black")) 
  },
  Frameshift_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[3], col = "black")) 
  },
  Frameshift_Insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[4], col = "black")) 
  },
  Nonframeshift_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[5], col = "black")) 
  },
  Nonframeshift_Insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[6], col = "black")) 
  },
  Splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[7], col = "black")) 
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[8], col = "black"))
  }
)

col = c("Nonsynonymous" = mycol[1], 
        "Stopgain" = mycol[2], 
        "Frameshift_Deletion" = mycol[3], 
        "Frameshift_Insertion" = mycol[4], 
        "Nonframeshift_Deletion" = mycol[5], 
        "Nonframeshift_Insertion" = mycol[6],
        "Splicing" = mycol[7],
        "Synonymous" = mycol[8])

tmp <- UTUC.annotation[colnames(mygene),c(9:11,14,15,17,19)]
ACTL6B <- read.table(file.path(res.path,"Methylation details for ACTL6B in C1 and C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,quote = "")
tmp$ACTL6B_Meth <- "N/A"
tmp[intersect(rownames(tmp),rownames(ACTL6B)),"ACTL6B_Meth"] <- ACTL6B[intersect(rownames(tmp),rownames(ACTL6B)),"MethStatus"]
# tmp$Del9p21.3 <- "N/A"
# tmp[intersect(rownames(tmp),rownames(del9p)),"Del9p21.3"] <- del9p[intersect(rownames(tmp),rownames(del9p)),"del9p"]
# tmp$Del9p21.3 <- ifelse(tmp$Del9p21.3 == "1","Del",ifelse(tmp$Del9p21.3 == "0",".","N/A"))

colnames(tmp) <- c("Gender","Type","Grade","Pathological_Stage","Lymph_Nodes","Localisation","Del9p21.3","ACTL6B_Meth")
my_annotation = HeatmapAnnotation(df = tmp,col = list(Gender = c("M" = blue, "F" = red),
                                                      Type = c("Non-muscle invasive" = grey, "muscle-invasive" = soil),
                                                      Grade = c("High" = sun, "Low" = lightblue),
                                                      Pathological_Stage = c("pTa"=white,"pT1"=yellow,"pT2"=green,"pT3"=sun,"pT4"=cherry),
                                                      Lymph_Nodes = c("N0"=white,"N1"=cyan,"N2"=violet),
                                                      Localisation = c("renal pelvis" = gold,"ureter"=cherry,"N/A"=lightgrey),
                                                      ACTL6B_Meth = c("Methylated"=purple,"Unmethylated"=white,"N/A"=lightgrey),
                                                      Del9p21.3 = c("Presence"=darkblue,"Absence"=white,"N/A"=lightgrey)))

p <- oncoPrint_plus(mygene, get_type = function(x) x,
                    alter_fun = alter_fun, col = col,
                    #row_order = NULL, # no order by mutation frequency
                    #column_order = c(meta.mut,local.mut), #no order by mutation frequency of sample
                    show_pct = T, #show pct in left
                    show_row_barplot = T, #show bar plot in right
                    #remove_empty_columns = TRUE,
                    column_title = "",
                    bottom_annotation = my_annotation,
                    show_heatmap_legend=T,
                    axis_gp = gpar(fontsize = 9),
                    column_title_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 7),
                    column_names_gp = gpar(fontsize = 6)) 
pdf(file.path(fig.path,"oncoprint_significant_UTUC_exome_mutation_with_detailedType_addACTL6B_adddel9p21.3.pdf"),width = 10,height = 5.5)
draw(p$oncoP)
dev.off()

oncoprint.sample.order <- colnames(mygene)[p$sampleOrder]
tmp <- mutload
index <- sapply(strsplit(names(tmp),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
names(tmp) <- index

tmp <- tmp[oncoprint.sample.order]
pdf(file = file.path(fig.path,"oncoprint top add mutation freqency bar.pdf"),width = 6,height = 3.5)
barplot(tmp,ylim = c(0,820),col = red,xaxt="n")
#text(x = a, y = tmp, label = as.numeric(tmp), pos = 3, cex = 0.4, col = "black")
invisible(dev.off())
##########################################################################
#### Make mutation comparison Table 1 with Guo, Hurst and TCGA oohorts ### 

#create mutation rate for genes that mutated >= 10% (3samples) or significant gene that pvalue < 0.05 in MutSigCV output
exome.mut.all <- read.table(file.path(data.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_wideformat.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
exome.mut.all <- exome.mut.all[,c(14,17:46)]
freq <- as.data.frame(table(exome.mut.all$Gene_Location))
exome.mut.all <- aggregate(exome.mut.all[2:31],by=list(exome.mut.all$Gene_Location),FUN=sum)
rownames(exome.mut.all) <- exome.mut.all$Group.1
exome.mut.all <- as.matrix(exome.mut.all[,-1])
exome.mut.all <- as.data.frame(floor(pmin(exome.mut.all,1.1)))
index <- sapply(strsplit(colnames(exome.mut.all),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
colnames(exome.mut.all) <- index
if(all(rownames(exome.mut.all)==freq$Var1)) {exome.mut.all.withFreq <- cbind.data.frame(exome.mut.all,cases.mut=freq$Freq,sam.mut=rowSums(exome.mut.all))}

write.table(exome.mut.all,file.path(res.path,"all exome mutation trimed.txt"),sep = "\t",row.names = T,col.names = NA)
write.table(exome.mut.all.withFreq,file.path(res.path,"all exome mutation trimed with case and sam frequency.txt"),sep = "\t",row.names = T,col.names = NA)


# check if utuc 12T,14T,25T,03T,07T have no mutation with others ######

# tmp <- exome.mut.all 
# UTUC.nomut.list <- list()
# 
# for (i in c("12T","14T","25T","03T","07T")) {
#   tmp1 <- rowSums(tmp[which(tmp[,i]==1),])
#   UTUC.nomut.list[[i]] <- tmp1[tmp1 >= 2]
# }
# only FBXW11 in 14T is one of the sig mut genes.

exome.mut.cut10 <- rownames(exome.mut.all[rowSums(exome.mut.all)>=3,])
exome.mut.p0.05 <- read.table(file.path(data.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_longformat_forMutSigCV_simplified.sig_genes.txt"),sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
exome.mut.p0.05 <- rownames(exome.mut.p0.05[which(exome.mut.p0.05$p < 0.05),])
exome.mut.sel <- exome.mut.all[union(exome.mut.cut10,exome.mut.p0.05),]

exome.mut.sel <- exome.mut.sel[order(rowSums(exome.mut.sel),decreasing = T),]
exome.mut.sel$Freq <- rowSums(exome.mut.sel)
exome.mut.sel$res <- paste0(exome.mut.sel$Freq," (",round(exome.mut.sel$Freq/30*100,0),")")
write.table(exome.mut.sel,file.path(res.path,"selected 198 cut10 and p0.05 ordered exome mutation trimed for Table 1.txt"),sep = "\t",row.names = T,col.names = NA)

# group for muscle-invasive and non-musle-invasive
exome.mut.sel <- exome.mut.sel[,setdiff(colnames(exome.mut.sel),c("Freq","res"))]
group <- UTUC.annotation[colnames(exome.mut.sel),"Type"]; names(group) <- colnames(exome.mut.sel)
exome.mut.sel <- as.data.frame(t(exome.mut.sel))
exome.mut.sel$group <- group[rownames(exome.mut.sel)]
outTab <- data.frame()
for (i in 1:198) {
  genename <- colnames(exome.mut.sel)[i]
  tmp1 <- table(exome.mut.sel$group,exome.mut.sel[,i])[2,2] #mutcount in non-muscle-invasive
  tmp2 <- table(exome.mut.sel$group,exome.mut.sel[,i])[1,2] #mutcount in muscle-invasive
  tmp <- data.frame("Nonagressive"=paste0(tmp1," (",round(tmp1/15*100,0),")"), "Agressive"=paste0(tmp2," (",round(tmp2/15*100,0),")"),stringsAsFactors = F,row.names = paste0(genename," (",round((tmp1+tmp2)/30*100,0),")"))
  outTab <- rbind.data.frame(outTab,tmp)
}
write.table(outTab,file.path(res.path,"Mutation rate of 198 cut10 and p0.05 ordered exome mutation in non muscle and muscle invasive.txt"),sep = "\t",row.names = T,col.names = NA)

mut.sel <- colnames(exome.mut.sel)[1:198]
mut.sel[mut.sel == "MLL3"] <- "KMT2C"; mut.sel[mut.sel == "MLL2"] <- "KMT2D"
Guo_nonaggressive_37 <- read.table(file.path(res.path,"BLCA_BGI_Nonaggressive_37_mutationRate.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
Guo_nonaggressive_37 <- Guo_nonaggressive_37[mut.sel,]
Guo_nonaggressive_37$rate <- round(as.numeric(gsub("%","",Guo_nonaggressive_37$Freq)),0)
Guo_nonaggressive_37$res <- paste0(Guo_nonaggressive_37$MutTrue," (",Guo_nonaggressive_37$rate,")")
write.table(Guo_nonaggressive_37,file.path(res.path,"Mutation rate of 198 genes in Guo_nonaggressive_37.txt"),sep = "\t",row.names = T,col.names = NA)

Guo_aggressive_62 <- read.table(file.path(res.path,"BLCA_BGI_Aggressive_62_mutationRate.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T) 
Guo_aggressive_62 <- Guo_aggressive_62[mut.sel,]
Guo_aggressive_62$rate <- round(as.numeric(gsub("%","",Guo_aggressive_62$Freq)),0)
Guo_aggressive_62$res <- paste0(Guo_aggressive_62$MutTrue," (",Guo_aggressive_62$rate,")")
write.table(Guo_aggressive_62,file.path(res.path,"Mutation rate of 198 genes in Guo_aggressive_62.txt"),sep = "\t",row.names = T,col.names = NA)

TCGA_aggressive_412 <- read.table(file.path(res.path,"BLCA_TCGA_Aggressive_412_mutationRate.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T) 
TCGA_aggressive_412 <- TCGA_aggressive_412[mut.sel,]
TCGA_aggressive_412$rate <- round(as.numeric(gsub("%","",TCGA_aggressive_412$Freq)),0)
TCGA_aggressive_412$res <- paste0(TCGA_aggressive_412$MutTrue," (",TCGA_aggressive_412$rate,")")
write.table(TCGA_aggressive_412,file.path(res.path,"Mutation rate of 198 genes in TCGA_aggressive_412.txt"),sep = "\t",row.names = T,col.names = NA)

Hurst_nonaggresive_raw <- read.table(file.path(res.path,"BLCA_Hurst_Nonaggressive_rawData.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,fill = T) 
Hurst_nonaggresive_raw$`Sample ID` <- as.character(Hurst_nonaggresive_raw$`Sample ID`)
Hurst_nonaggresive_raw$sampleMut <- paste0(Hurst_nonaggresive_raw$`Gene symbol`,"_",Hurst_nonaggresive_raw$`Sample ID`)
Hurst_nonaggresive_raw <- Hurst_nonaggresive_raw[!duplicated(Hurst_nonaggresive_raw$sampleMut),]
Hurst_gene <- unique(Hurst_nonaggresive_raw$`Gene symbol`)
Hurst_sample <- unique(Hurst_nonaggresive_raw$`Sample ID`)
Hurst_nonaggresive_matrix <- matrix(0,nrow=length(Hurst_gene),ncol=length(Hurst_sample),dimnames = list(Hurst_gene,Hurst_sample))
for (i in 1:nrow(Hurst_nonaggresive_raw)) {
 tmp <- Hurst_nonaggresive_raw[i,]
 Hurst_nonaggresive_matrix[tmp$`Gene symbol`,tmp$`Sample ID`] <- 1
}
Hurst_nonaggresive_24 <- as.data.frame(Hurst_nonaggresive_matrix)
Hurst_nonaggresive_24 <- Hurst_nonaggresive_24[mut.sel,]
Hurst_nonaggresive_24$Freq <- rowSums(Hurst_nonaggresive_24)
Hurst_nonaggresive_24$res <- paste0(Hurst_nonaggresive_24$Freq, " (",round(Hurst_nonaggresive_24$Freq/24*100,0),")")
write.table(Hurst_nonaggresive_24,file.path(res.path,"Mutation rate of 198 genes in Hurst_nonaggresive_24.txt"),sep = "\t",row.names = T,col.names = NA)

exome.mut.sel <- as.data.frame(exome.mut.sel[,-199])
genelist <- colnames(exome.mut.sel)
binarymut <- as.data.frame(matrix(0, nrow=length(group), ncol=length(genelist)))
rownames(binarymut) <- names(group)
colnames(binarymut) <- genelist
for (k in 1:length(genelist)) {
  res <- create.anntrack(samples=names(group), subtype=createMutSubtype(t(exome.mut.sel), names(group), genelist[k]))
  binarymut[, genelist[k]] <- res$subtype
}

out <- matrix(0, nrow=length(genelist), ncol=length(unique(group))+1)
colnames(out) <- c(levels(factor(group)), "pvalue")
rownames(out) <- paste(genelist, "Mutated", sep="_")
for (k in 1:length(genelist)) {
  genek <- genelist[k]
  x <- group
  y <- binarymut[names(x), genek, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["Mutated", ]
  
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
write.table(out, file.path(res.path, "Independence test between 198 cut10 and p0.05 exome mutation and muscle-invasive type.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

##########################################################################
######## link mutation to methylation by SWI/SNF and MAPK pathway ########

###########
# SWI/SNF #
SWI_SNF <- c("ARID1A",	"PBRM1",	"SMARCA4",	"ARID1B",	"ARID2",
             "SMARCA2",	"SMARCB1",	"SMARCC1",	"SMARCC2",	"DPF2",
             "ACTL6A",	"SMARCD3",	"PHF10",	"BRD7",	"SMARCD1",
             "SMARCD2",	"SMARCE1",	"ACTL6B",	"DPF1",	"DPF3")

exome.mut.swisnf <- exome.mut.all[intersect(SWI_SNF,rownames(exome.mut.all)),]

#non-muscle and muscle invasive
group <- UTUC.annotation[colnames(exome.mut.swisnf),"Type"]; names(group) <- colnames(exome.mut.swisnf)
plotdata <- as.data.frame(t(exome.mut.swisnf))

plotdata$ACTL6B <- 0
plotdata[intersect(rownames(plotdata),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(plotdata),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)

plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:7], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Type","Gene","Count")
plotdata$Percentage <- round(plotdata$Count/15,2)
p1 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = Type, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 6genes Pct between non-muscle and muscle invasieve addACTL6B.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.swisnf),"Type"]; names(group) <- colnames(exome.mut.swisnf)
plotdata <- as.data.frame(t(exome.mut.swisnf))
plotdata <- as.data.frame(t(exome.mut.swisnf[,names(group)]))
plotdata$ACTL6B <- 0
plotdata[intersect(rownames(plotdata),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(plotdata),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)

tmp <- data.frame("SWI_SNF_Mut"=rowSums(plotdata),"Type"=UTUC.annotation[rownames(plotdata),"Type"],row.names = rownames(plotdata),"DMBcluster"=UTUC.annotation[rownames(plotdata),"Cluster DNA methylation no filter without normal"])
tmp[tmp$SWI_SNF_Mut>=1,"SWI_SNF_Mut"] <- 1
# muscle-invasive Non-muscle invasive
# 0              8                  9
# 1              7                  6
fisher.test(table(tmp$SWI_SNF_Mut,tmp$Type)) #p=1
table(tmp$SWI_SNF_Mut,tmp$DMBcluster)
# C1 C2 N/A
# 0  7  9   1
# 1 11  1   1
write.table(tmp,file.path(res.path,"SWI_SNF_Mut status in UTUC combined ACTL6B meth.txt"),sep = "\t",row.names = T,col.names = NA)

group <- UTUC.annotation[colnames(exome.mut.swisnf),"Type"]; names(group) <- colnames(exome.mut.swisnf)
plotdata <- as.data.frame(t(exome.mut.swisnf))
plotdata$ACTL6B <- 0
plotdata[intersect(rownames(plotdata),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(plotdata),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)
SWISNF_mut_addACTL6B <- rowSums(plotdata); SWISNF_mut_addACTL6B <- ifelse(SWISNF_mut_addACTL6B >= 1,1,0)

plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:7], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Type","Gene","Count")
p1 <- ggplot() + 
  geom_bar(aes(y = Count, x = Type, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0,13, 1), limits=c(0, 13)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 6genes Count between non-muscle and muscle invasieve addACTL6B.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.swisnf),"Type"]; names(group) <- colnames(exome.mut.swisnf)
plotdata <- as.data.frame(t(exome.mut.swisnf))
SWISNF_mut_noACTL6B <- rowSums(plotdata); SWISNF_mut_noACTL6B <- ifelse(SWISNF_mut_noACTL6B >= 1,1,0)
save(SWISNF_mut_noACTL6B,file=file.path(res.path,"SWISNF_mut_noACTL6B.rda"))
#DMBcluster C1 and C2
# with ACTL6B
group <- UTUC.annotation[colnames(exome.mut.swisnf),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.swisnf[,names(group)]))

plotdata$ACTL6B <- 0
plotdata[intersect(rownames(plotdata),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(plotdata),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)

plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:7], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("DMBcluster","Gene","Count")
plotdata$Percentage <- ifelse(plotdata$DMBcluster == "C1",round(plotdata$Count/18,2),round(plotdata$Count/10,2))
p1 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = DMBcluster, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(7)) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 1.1, 0.1), limits=c(0,1.1)) + guides(fill=guide_legend(title="SWI/SNF Pathway")) #+ theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 6genes Pct between DMBcluster C1 and C2 addACTL6B.pdf"),width = 6)

p2 <- ggplot() + 
  geom_bar(aes(y = Count, x = DMBcluster, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(7)) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 18, 1), limits=c(0,18)) + guides(fill=guide_legend(title="SWI/SNF Pathway")) #+ theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 6genes Count between DMBcluster C1 and C2 addACTL6B.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.swisnf),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.swisnf[,names(group)]))
plotdata$ACTL6B <- 0
plotdata[intersect(rownames(plotdata),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(plotdata),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)
tmp <- data.frame("SWI_SNF_Mut"=rowSums(plotdata),"DMBcluster"=UTUC.annotation[rownames(plotdata),"Cluster DNA methylation no filter without normal"],row.names = rownames(plotdata))
tmp[tmp$SWI_SNF_Mut>=1,"SWI_SNF_Mut"] <- 1

#   C1 C2
# 0  7  9
# 1 11  1
tmp$SWI_SNF_Mut <- ifelse(tmp$SWI_SNF_Mut == "1","Presence","Absence")
fisher.test(table(tmp$SWI_SNF_Mut,tmp$DMBcluster)) #p=0.01587
tmp <- as.data.frame(table(tmp$SWI_SNF_Mut,tmp$DMBcluster))
colnames(tmp) <- c("Mutation","DMBcluster","Frequency")
tmp <- ddply(tmp,"DMBcluster",transform,label_y=cumsum(Frequency)); tmp$Mutation <- factor(tmp$Mutation,levels = c("Presence","Absence"))
p3 <- ggplot(tmp,aes(x=DMBcluster,y=Frequency,fill=Mutation)) + geom_bar(stat="identity") + scale_fill_manual(values = c(cherry,lightgrey)) +
      scale_y_continuous(breaks=seq(0, 18, 3), limits=c(0,18) ) + guides(fill=guide_legend(title="SWI/SNF Mutation"))
ggsave(filename = file.path(fig.path,"SWI_SNF_mutation addACTL6B frequency between DMBcluster C1 and C2 dodged plot with pvalue = 0.01587.pdf"),width = 4,height = 5)

p <- plot_grid(p2, p3, labels = c("a", "b"))
ggsave(file.path(fig.path,"SWI_SNF_mutation addACTL6B between DMBcluster C1 and C2 .pdf"),width = 8,height = 5)

# without ACTL6B
group <- UTUC.annotation[colnames(exome.mut.swisnf),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.swisnf[,names(group)]))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:6], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
p2 <- ggplot() + 
  geom_bar(aes(y = Count, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 13, 1), limits=c(0, 13))
ggsave(file.path(fig.path,"Stacked bar of SWI_SNF 6genes Count between DMBcluster C1 and C2.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.swisnf),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.swisnf[,names(group)]))
tmp <- data.frame("SWI_SNF_Mut"=rowSums(plotdata),"DMBcluster"=UTUC.annotation[rownames(plotdata),"Cluster DNA methylation no filter without normal"],row.names = rownames(plotdata))
tmp[tmp$SWI_SNF_Mut>=1,"SWI_SNF_Mut"] <- 1
tmp$SWI_SNF_Mut <- ifelse(tmp$SWI_SNF_Mut == "1","Presence","Absence")
# C1 C2
# Absence   9 10
# Presence  9  0
fisher.test(table(tmp$SWI_SNF_Mut,tmp$DMBcluster)) #p=0.009816
tmp <- as.data.frame(table(tmp$SWI_SNF_Mut,tmp$DMBcluster))
colnames(tmp) <- c("Mutation","DMBcluster","Frequency")
tmp <- ddply(tmp,"DMBcluster",transform,label_y=cumsum(Frequency)); tmp$Mutation <- factor(tmp$Mutation,levels = c("Presence","Absence"))
p3 <- ggplot(tmp,aes(x=DMBcluster,y=Frequency,fill=Mutation)) + geom_bar(stat="identity") + scale_fill_manual(values = c(cherry,lightgrey)) +
  scale_y_continuous(breaks=seq(0, 18, 3), limits=c(0,18) ) + guides(fill=guide_legend(title="SWI/SNF Mutation"))
ggsave(filename = file.path(fig.path,"SWI_SNF_mutation frequency between DMBcluster C1 and C2 dodged plot with pvalue = 0.009816.pdf"),width = 4,height = 5)

p <- plot_grid(p2, p3, labels = c("a", "b"))
ggsave(file.path(fig.path,"SWI_SNF_mutation between DMBcluster C1 and C2 .pdf"),width = 8,height = 5)

# mutual exclusivity for SWI/SNF and ACTL6B
exome.mut.swisnf <- exome.mut.all[intersect(SWI_SNF,rownames(exome.mut.all)),]
plotdata <- as.data.frame(t(exome.mut.swisnf))
tmp <- data.frame("SWI_SNF_Mut"=rowSums(plotdata),row.names = rownames(plotdata))
tmp[tmp$SWI_SNF_Mut>=1,"SWI_SNF_Mut"] <- 1
tmp$ACTL6B <- 0
tmp[intersect(rownames(tmp),rownames(ACTL6B)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(tmp),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)
fisher.test(tmp$SWI_SNF_Mut,tmp$ACTL6B,alternative = "greater")

# Fisher's Exact Test for Count Data
# 
# data:  tmp$SWI_SNF_Mut and tmp$ACTL6B
# p-value = 0.3064
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.3453276       Inf
# sample estimates:
# odds ratio 
#   2.351417 

################
# MAPK pathway #
MAPK <- read.table(file.path(data.path,"MAPK.geneset.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
MAPK <- MAPK$`GENE LIST`

exome.mut.mapk <- exome.mut.all[intersect(MAPK,rownames(exome.mut.all)),] #53 30

#non-muscle and muscle invasive
group <- UTUC.annotation[colnames(exome.mut.mapk),"Type"]; names(group) <- colnames(exome.mut.mapk)
plotdata <- as.data.frame(t(exome.mut.mapk))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:53], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Type","Gene","Count")
plotdata$Percentage <- round(plotdata$Count/15,2)
p3 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = Type, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(53)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())  + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of MAPK 53genes Pct between non-muscle and muscle invasive.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.mapk),"Type"]; names(group) <- colnames(exome.mut.mapk)
plotdata <- as.data.frame(t(exome.mut.mapk))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:53], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Type","Gene","Count")
p3 <- ggplot() + 
  geom_bar(aes(y = Count, x = Type, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(53)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) + theme(legend.position="bottom")+ scale_y_continuous(breaks=seq(0, 55, 5), limits=c(0, 55))
ggsave(file.path(fig.path,"Stacked bar of MAPK 53genes Count between non-muscle and muscle invasive.pdf"),width = 6)

#DMBcluster C1 and C2
group <- UTUC.annotation[colnames(exome.mut.mapk),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.mapk[,names(group)]))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:53], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
plotdata$Percentage <- ifelse(plotdata$Class == "C1",round(plotdata$Count/18,2),round(plotdata$Count/10,2))
p4 <- ggplot() + 
  geom_bar(aes(y = Percentage, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(53)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of MAPK 53genes Pct between DMBcluster C1 and C2.pdf"),width = 6)

group <- UTUC.annotation[colnames(exome.mut.mapk),"Cluster DNA methylation no filter without normal"]; names(group) <- colnames(exome.mut.swisnf)
group <- group[setdiff(names(group),c("23T","07T"))]
plotdata <- as.data.frame(t(exome.mut.mapk[,names(group)]))
plotdata$Class <- group[rownames(plotdata)]
plotdata <- apply(plotdata[,1:53], 2, FUN = function(x)tapply(x,plotdata$Class,FUN = sum))
plotdata <- melt(plotdata); colnames(plotdata) <- c("Class","Gene","Count")
p4 <- ggplot() + 
  geom_bar(aes(y = Count, x = Class, fill = Gene), data = plotdata,stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(53)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 60, 5), limits=c(0,60)) + theme(legend.position="bottom")
ggsave(file.path(fig.path,"Stacked bar of MAPK 53genes Count between DMBcluster C1 and C2.pdf"),width = 6)

#############################################################################
#### Mutual exclusivity analysis for pair-wise significant mutated genes ####
tmp <- t(exome.mut)
tmp <- as.data.frame(tmp[,c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB","ARID1B")])
tmp <- tmp[,-11]
# combine ARID1A and ARID1B
#tmp$'ARID1A/B' <- tmp$ARID1A + tmp$ARID1B
#tmp <- tmp[,-c(8,11)]
#addACTL6B
#tmp$ACTL6B <- 0
#tmp[intersect(rownames(ACTL6B),rownames(tmp)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(tmp),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)
#tmp <- tmp[,colSums(tmp!=0) > 1] # just keep mutation that has >= 2 cases mutation
#tmp <- tmp[rowSums(tmp!=0) >= 1,] # just keep samples that has >= 1 mutation
# perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
# "less" for mutual-exclusivity analysis, or "greater" for co-occurrence analysis 

res <- NULL
for(i in 1:ncol(tmp)){
  for(j in i:ncol(tmp)){
    tab <- table(tmp[,i], tmp[,j])
    if(i!=j){
      f <- fisher.test(tab) 
      # deal with zero in 2x2 cell
      if(f$estimate == 0) {
        f <- fisher.test(tab,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      } else {
        if(is.infinite(f$estimate)) {
          f <- fisher.test(tab + 1)
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab + 1,alternative = "greater")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
          } else{
            f <- fisher.test(tab + 1,alternative = "less")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
          }
          
        } else {
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab,alternative = "greater")
          } else{f <- fisher.test(tab,alternative = "less")}
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
        }   
      }
    }
  }
}
# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
#res$geneA <- factor(res$geneA,levels=c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","TP53","GANAB","ARID1A","ARID1B","ARID1A/B","ACTL6B"))
#res$geneB <- factor(res$geneB,levels=c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","TP53","GANAB","ARID1A","ARID1B","ARID1A/B","ACTL6B"))
res$geneA <- factor(res$geneA,levels=c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB"))
res$geneB <- factor(res$geneB,levels=c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB"))
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$log2OR <- log2(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
# plot with ggplot 2
#write.table(res,file.path(res.path,"Mutual_exclusivity_for_significant_mutated_genes_addACTL6B_rmNonMut_samples.txt"),sep = "\t",row.names = F)
write.table(res,file.path(res.path,"Mutual_exclusivity_for_significant_mutated_genes_keepNonMut_samples.txt"),sep = "\t",row.names = F)

p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = log2OR),colour = "white") + scale_fill_gradient2(low = "darkblue",mid = lightgrey,high = "darkred",midpoint=0) + 
  geom_text(aes(label=stars), color="white", size=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank())
p
ggsave(file.path(fig.path,"Mutual_exclusivity_for_significant_mutated_genes_keepNonMut_samples.pdf"),width = 5,height = 5)
#ggsave(file.path(fig.path,"Mutual_exclusivity_for_significant_mutated_genes_combinedARID1AB_addACTL6B_rmNonMut_samples.pdf"),width = 5,height = 5)
#ggsave(file.path(fig.path,"Mutual_exclusivity_for_significant_mutated_genes_combinedARID1AB_addACTL6B_keepNonMut_samples.pdf"),width = 5,height = 5)

# mutual exclusivity between del9p21.3 and significant mutation genes and ACTL6B meth
tmp$Del9p21.3 <- NA
#tmp[intersect(rownames(del9p),rownames(tmp)),"Del9p21.3"] <- del9p[intersect(rownames(del9p),rownames(tmp)),"del9p"]
tmp[intersect(rownames(del9p.exome),rownames(tmp)),"Del9p21.3"] <- del9p.exome[intersect(rownames(del9p.exome),rownames(tmp)),"del9p"]
tmp <- as.data.frame(na.omit(tmp))

res <- NULL
for (i in 1:(ncol(tmp)-1)) {
  tab <- table(tmp[,i], tmp[,ncol(tmp)])
    f <- fisher.test(tab) 
    # deal with zero in 2x2 cell
    if(f$estimate == 0) {
      f <- fisher.test(tab,alternative = "less")
      res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                   Neither=tab[1,1],
                                                   AnotB=tab[2,1],
                                                   BnotA=tab[1,2],
                                                   Both=tab[2,2],
                                                   oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
    } else {
      if(is.infinite(f$estimate)) {
        f <- fisher.test(tab + 1)
        if(log(f$estimate) > 0) {
          f <- fisher.test(tab + 1,alternative = "greater")
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
        } else{
          f <- fisher.test(tab + 1,alternative = "less")
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
        }
        
      } else {
        if(log(f$estimate) > 0) {
          f <- fisher.test(tab,alternative = "greater")
        } else{f <- fisher.test(tab,alternative = "less")}
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
      
    }   
  }
}
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$logOR <- log(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
write.table(res,file.path(res.path,"Mutual_exclusivity_for_Del9p23 and significant_mutated_genes_addACTL6B_rmNonMut_samples.txt"),sep = "\t",row.names = F)
#write.table(res,file.path(res.path,"Mutual_exclusivity_for_Del9p23 and significant_mutated_genes_addACTL6B_keepNonMut_samples.txt"),sep = "\t",row.names = F)


# mutual exclusivity between SWI/SNF and significant mutation genes and ACTL6B meth
tmp <- t(exome.mut)
tmp <- as.data.frame(tmp[,c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB","ARID1B")])
# combine ARID1A and ARID1B
tmp$'ARID1A/B' <- tmp$ARID1A + tmp$ARID1B
#tmp <- tmp[,-c(8,11)]
#addACTL6B
tmp$ACTL6B <- 0
tmp[intersect(rownames(ACTL6B),rownames(tmp)),"ACTL6B"] <- ifelse(ACTL6B[intersect(rownames(tmp),rownames(ACTL6B)),"MethStatus"] == "Methylated",1,0)
#tmp <- tmp[,colSums(tmp!=0) > 1] # just keep mutation that has >= 2 cases mutation
tmp$SWI_SNF <- SWISNF_mut_addACTL6B[rownames(tmp)]
tmp <- tmp[rowSums(tmp!=0) >= 1,] # just keep samples that has >= 1 mutation

res <- NULL
for (i in 1:(ncol(tmp)-1)) {
  tab <- table(tmp[,i], tmp[,ncol(tmp)])
  f <- fisher.test(tab) 
  # deal with zero in 2x2 cell
  if(f$estimate == 0) {
    f <- fisher.test(tab,alternative = "less")
    res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                 Neither=tab[1,1],
                                                 AnotB=tab[2,1],
                                                 BnotA=tab[1,2],
                                                 Both=tab[2,2],
                                                 oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
  } else {
    if(is.infinite(f$estimate)) {
      f <- fisher.test(tab + 1)
      if(log(f$estimate) > 0) {
        f <- fisher.test(tab + 1,alternative = "greater")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
      } else{
        f <- fisher.test(tab + 1,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      }
      
    } else {
      if(log(f$estimate) > 0) {
        f <- fisher.test(tab,alternative = "greater")
      } else{f <- fisher.test(tab,alternative = "less")}
      res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[ncol(tmp)],
                                                   Neither=tab[1,1],
                                                   AnotB=tab[2,1],
                                                   BnotA=tab[1,2],
                                                   Both=tab[2,2],
                                                   oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
      
    }   
  }
}
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$logOR <- log(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
write.table(res,file.path(res.path,"Mutual_exclusivity_for_SWISNF and significant_mutated_genes_addACTL6B_rmNonMut_samples.txt"),sep = "\t",row.names = F)
#write.table(res,file.path(res.path,"Mutual_exclusivity_for_SWISNF and significant_mutated_genes_addACTL6B_keepNonMut_samples.txt"),sep = "\t",row.names = F)

########################################################################################################################
#### Mutual exclusivity analysis for pair-wise significant mutated genes use non/muscle invasive separately in UTUC ####

# HAVENT UPDATE !!!

# muscle invasive
samples <- rownames(UTUC.annotation[which(UTUC.annotation$`Exome-data` == "yes" & UTUC.annotation$Type == "muscle-invasive"),])
tmp <- t(exome.mut[,samples])
tmp <- tmp[,c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B")]
tmp <- tmp[,colSums(tmp!=0) >= 1] # just keep mutation that has >= 1 cases mutation
tmp <- tmp[rowSums(tmp!=0) >= 1,] # just keep samples that has >= 1 mutation

# perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
# "less" for mutual-exclusivity analysis, or "greater" for co-occurrence analysis (need work?)
res <- NULL
for(i in 1:ncol(tmp)){
  for(j in i:ncol(tmp)){
    tab <- table(tmp[,i], tmp[,j])
    if(i!=j){
      f <- fisher.test(tab) 
      # deal with zero in 2x2 cell
      if(f$estimate == 0) {
        f <- fisher.test(tab,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      } else {
        if(is.infinite(f$estimate)) {
          f <- fisher.test(tab + 1)
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab + 1,alternative = "greater")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
          } else{
            f <- fisher.test(tab + 1,alternative = "less")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
          }
          
        } else {
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab,alternative = "greater")
          } else{f <- fisher.test(tab,alternative = "less")}
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
        }   
      }
    }
  }
}
# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$geneA <- factor(res$geneA,levels=c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B"))
res$geneB <- factor(res$geneB,levels=c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B"))
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$logOR <- log(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
# plot with ggplot 2
write.table(res,file.path(res.path,"Mutual_exclusivity_for_significant_mutated_genes_use_muscleinvasive_UTUC_samples.txt"),sep = "\t",row.names = F)
p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = logOR),colour = "white") + scale_fill_gradient2(low = "darkblue",mid = "white",high = "darkred",midpoint=0) + 
  geom_text(aes(label=stars), color="white", size=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank())
p
ggsave(file.path(fig.path,"Mutual_exclusivity_for_significant_mutated_genes_use_muscleinvasive_UTUC_samples.pdf"),width = 5,height = 5)

# nonmuscle invasive
samples <- rownames(UTUC.annotation[which(UTUC.annotation$`Exome-data` == "yes" & UTUC.annotation$Type == "Non-muscle invasive"),])
tmp <- t(exome.mut[,samples])
tmp <- tmp[,c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B")]
tmp <- tmp[,colSums(tmp!=0) >= 1] # just keep mutation that has >= 1 cases mutation
tmp <- tmp[rowSums(tmp!=0) >= 1,] # just keep samples that has >= 1 mutation

# perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
# "less" for mutual-exclusivity analysis, or "greater" for co-occurrence analysis (need work?)
res <- NULL
for(i in 1:ncol(tmp)){
  for(j in i:ncol(tmp)){
    tab <- table(tmp[,i], tmp[,j])
    if(i!=j){
      f <- fisher.test(tab) 
      # deal with zero in 2x2 cell
      if(f$estimate == 0) {
        f <- fisher.test(tab,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      } else {
        if(is.infinite(f$estimate)) {
          f <- fisher.test(tab + 1)
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab + 1,alternative = "greater")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
          } else{
            f <- fisher.test(tab + 1,alternative = "less")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
          }
          
        } else {
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab,alternative = "greater")
          } else{f <- fisher.test(tab,alternative = "less")}
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
        }   
      }
    }
  }
}
# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$geneA <- factor(res$geneA,levels=c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B"))
res$geneB <- factor(res$geneB,levels=c("FGFR3","KDM6A","MLL2","ZFP36L1","ARID1A","STAG2","GANAB","EPHA5","GABRA5","GLYR1","ARID1B"))
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$logOR <- log(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
# plot with ggplot 2
write.table(res,file.path(res.path,"Mutual_exclusivity_for_significant_mutated_genes_use_nonmuscleinvasive_UTUC_samples.txt"),sep = "\t",row.names = F)
p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = logOR),colour = "white") + scale_fill_gradient2(low = "darkblue",mid = "white",high = "darkred",midpoint=0) + 
  geom_text(aes(label=stars), color="white", size=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank())
p
ggsave(file.path(fig.path,"Mutual_exclusivity_for_significant_mutated_genes_use_nonmuscleinvasive_UTUC_samples.pdf"),width = 5,height = 5)

###########################################################################################
#### Mutual exclusivity analysis for Hurst cohort (others could be down by cBioportal) ####
tmp <- read.table(file.path(res.path,"Mutation rate of 198 genes in Hurst_nonaggresive_24.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- as.data.frame(na.omit(tmp))
rownames(tmp) <- tmp[,1]
tmp <- tmp[,2:25]
tmp <- tmp[intersect(c("KMT2D","FGFR3","KMT2C","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB","ARID1B"),rownames(tmp)),]
tmp <- t(tmp)
#tmp <- tmp[rowSums(tmp!=0) >= 1,] # just keep samples that has >= 1 mutation
# perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
# "less" for mutual-exclusivity analysis, or "greater" for co-occurrence analysis (need work?)
res <- NULL
for(i in 1:ncol(tmp)){
  for(j in i:ncol(tmp)){
    tab <- table(tmp[,i], tmp[,j])
    if(i!=j){
      f <- fisher.test(tab) 
      # deal with zero in 2x2 cell
      if(f$estimate == 0) {
        f <- fisher.test(tab,alternative = "less")
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
      } else {
        if(is.infinite(f$estimate)) {
          f <- fisher.test(tab + 1)
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab + 1,alternative = "greater")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
          } else{
            f <- fisher.test(tab + 1,alternative = "less")
            res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                         Neither=tab[1,1],
                                                         AnotB=tab[2,1],
                                                         BnotA=tab[1,2],
                                                         Both=tab[2,2],
                                                         oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
          }
          
        } else {
          if(log(f$estimate) > 0) {
            f <- fisher.test(tab,alternative = "greater")
          } else{f <- fisher.test(tab,alternative = "less")}
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp)[i],geneB=colnames(tmp)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
        }   
      }
    }
  }
}
# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$geneA <- factor(res$geneA,levels=c( "KMT2D",   "FGFR3",   "KMT2C",   "KDM6A",   "ZFP36L1", "STAG2",   "CRIPAK",  "ARID1A",  "TP53"))
res$geneB <- factor(res$geneB,levels=c( "KMT2D",   "FGFR3",   "KMT2C",   "KDM6A",   "ZFP36L1", "STAG2",   "CRIPAK",  "ARID1A",  "TP53"))
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$logOR <- log(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))

write.table(res,file.path(res.path,"BLCA_Hurst_Nonaggresive_24_Mutual_exclusivity.tsv"),sep = "\t",row.names = F,quote = F)

###################################################################
#### Mutual exclusivity analysis for FGFR3 and SWI/SNF 6 genes ####

# tmp <- t(rbind(exome.mut["FGFR3",],exome.mut.swisnf))
# res <- NULL
# for(i in 2:ncol(tmp)){
#      f <- fisher.test(tmp[,"FGFR3"], tmp[,i])
#      res <- rbind(res,cbind(geneA="FGFR3",geneB=colnames(tmp)[i],oddsRatio=f$estimate,pvalue=f$p.value))
# }

# nothing significant

######################################################################
#### Mutual exclusivity analysis for FGFR3 and mutation signature ####
# tmp <- read.table("F:/Project/UTUC_BLCA/exome UTUC/Results/Mut.sig.cluster3 results of DMBcluster NMFbasis MWBcluster.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# rownames(tmp) <- gsub("UTUC_","",rownames(tmp))
# tmp$FGFR3_Mut <- as.numeric(exome.mut["FGFR3",rownames(tmp)])
# tmp$Sig13 <- ifelse(tmp$NMFbasis == 2, 1,0)
# tmp$Sig16 <- ifelse(tmp$NMFbasis == 1, 1,0)
# tmp$Sig1 <- ifelse(tmp$NMFbasis == 3, 1,0)
# 
# fisher.test(tmp$FGFR3_Mut,tmp$Sig13)
# fisher.test(tmp$FGFR3_Mut,tmp$Sig16)
# fisher.test(tmp$FGFR3_Mut,tmp$Sig1)

# nothing significant

###########################################################################
#### Unsupervised clustering for UTUC RNAseq with batch effect removel ####

annCol.rna <- UTUC.annotation
annCol.rna <- annCol.rna[which(annCol.rna$`RNAseq (please add)` %in% "yes"),]; rownames(annCol.rna) <- annCol.rna$DetailID

tmp <- as.data.frame(t(est.score.combat)); tmp <- tmp[,setdiff(colnames(tmp),c("TumorPurity"))]; tmp <- as.data.frame(t(tmp),stringsAsFactors=F)
tmp <- as.data.frame(sapply(tmp, as.numeric)); rownames(tmp) <- setdiff(rownames(est.score.combat),"TumorPurity")
tmp <- annTrackScale(indata = tmp, halfwidth = 2, poolsd = F); tmp <- as.data.frame(t(tmp))

annCol.rna <- cbind.data.frame(MCPscore.combat[rownames(annCol.rna),],tmp[rownames(annCol.rna),1:2],annCol.rna[,c(8,9,10,11,14,15,17,19,20,21,22,23,30:45)])

index <- sapply(strsplit(rownames(annCol.rna),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
rownames(annCol.rna) <- index

# annCol.rna$sample <- rownames(annCol.rna)
# tmp <- as.data.frame(t(exome.mut[, intersect(colnames(exome.mut),rownames(annCol.rna))])); tmp$sample <- rownames(tmp)
# annCol.rna <- merge(annCol.rna,tmp,by="sample",all.x=T)
# rownames(annCol.rna) <- annCol.rna$sample; annCol.rna <- annCol.rna[,-1]
annCol.rna[is.na(annCol.rna)] <- "N/A"
colnames(annCol.rna)[13:24] <- c("DMBcluster","Gender","Type","Grade","Pathological_Stage","Lymph_Nodes","Localisation","Del9p21.3","RNAcluster","MWBcluster","iCluster","FGFR3_Mut_add")
annCol.rna$Age <- as.numeric(UTUC.annotation[rownames(annCol.rna),"Age"])
annCol.rna$Age <- ifelse(annCol.rna$Age > 82,">82","<=82")

annColors.rna <- annColors.C1vsC2.rna
annColors.rna[["Gender"]] <- c("M" = blue, "F" = darkred)
annColors.rna[["Type"]] <- c("Non-muscle invasive" = grey, "muscle-invasive" = soil)
annColors.rna[["Grade"]] <- c("High" = sun, "Low" = lightblue)
annColors.rna[["Pathological_Stage"]] <- c("pTa" = white,"pT1" = yellow,"pT2" = green,"pT3" = sun,"pT4" = cherry)
annColors.rna[["Lymph_Nodes"]] <- c("N0"=white,"N1"=cyan,"N2"=violet)
annColors.rna[["Localisation"]] <- c("renal pelvis" = gold,"ureter"=cherry)
annColors.rna[["FGFR3_Mut_add"]] <- c("1"=purple,"0"=lightgrey)
annColors.rna[["Del9p21.3"]] <- c("Presence"=purple,"Absence"=lightgrey,"N/A"="white")
#annColors.rna[["RNAcluster"]] <- c("C1"=purple,"C2"=lightgrey,"N/A"="white")
annColors.rna[["RNAcluster"]] <- c("Muscle_Enriched"=cherry,"Non-muscle_Enriched"=soil,"N/A"="white")
annColors.rna[["iCluster"]] <- c("C1"=jco[2],"C2"=jco[2],"N/A"="white")
annColors.rna[["WMBcluster"]] <- c("C1"="#F09300","C2"=lightgreen,"C3"="#0068B5")
annColors.rna[["DMBcluster"]] <- c("C1"=red,"C2"=blue,"N/A"="white")
annColors.rna[["MLL2"]] <- annColors.rna[["FGFR3"]] <- annColors.rna[["MLL3"]] <- 
  annColors.rna[["KDM6A"]] <- annColors.rna[["ZFP36L1"]] <- annColors.rna[["GABRA5"]] <- 
  annColors.rna[["CRIPAK"]] <- annColors.rna[["ARID1A"]] <- annColors.rna[["TP53"]] <-
  annColors.rna[["GANAB"]] <- annColors.rna[["ARID1B"]] <- annColors.rna[["ARID1A/B"]] <- 
  annColors.rna[["ZFP36L2"]] <- annColors.rna[["ZFP36"]] <- annColors.rna[["ZFP36_family"]] <- 
  annColors.rna[["SWI_SNF"]] <- c("1"=purple,"0"=lightgrey,"N/A"="white")
annColors.rna[["Age"]] <- c(">82"="black","<=82"=nake)
# k = 1000
# sel_idx <- order(rowVars(combat.UTUC.FPKM), decreasing=T)[1:k]
# vst_sel <- combat.UTUC.FPKM[sel_idx,]
# index <- sapply(strsplit(colnames(vst_sel),"_"), "[",2)
# index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
# colnames(vst_sel) <- index
# 
# hcg <- hclust(distanceMatrix(as.matrix(t(vst_sel)), "pearson"), "ward.D")
# hcs <- hclust(distanceMatrix(as.matrix(vst_sel), "euclidean"), "ward.D")
# plotdata <- standarize.fun(indata=vst_sel, halfwidth=3, centerFlag=T, scaleFlag=T)
# 
# pdf(file.path(fig.path,paste0("Unsupervised cluster on UTUC 20 RNAseq samples with removel of Batcheffect for top ",k," high variation genes.pdf")),height=12)
# hv = aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol.rna[colnames(plotdata),], annColors=annColors.rna, color=greenred(128), revC=TRUE, fontsize=5,cexCol = 0.2,cexRow = 0.2,cexAnn = 0.6, labRow = NA,labCol = NA)
# invisible(dev.off())

rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC[PASSFlag.mRNA.UTUC==TRUE]) )
indata <- combat.UTUC.FPKM[rowids,]
index <- sapply(strsplit(colnames(indata),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
colnames(indata) <- index
indata <- sweep(indata,2,apply(indata,2,median))
indata <- sweep(indata,1,apply(indata,1,median))

batchPCA(indata = indata,batch = colnames(indata),fig.dir = fig.path,PCA.fig.title = "UTUC.combat_FPKM.RNAseq.Samples.PCA",cols = rep(blue,20),showID = T,cex = 0.7,showLegend = F)

# dist.list <- c("euclidean","pearson","manhattan")
# link.list <- c("ward.D","ward.D2","average")
# 
# N.cluster <- 3
# N.bootstrap <- 500
# N.gene.per.bootstrap <- round(0.8*nrow(indata))
# N.sample.per.bootstrap <- round(0.8*ncol(indata))
# 
# res.list <- list()
# for (clusterNum in N.cluster) {
#   for (dist0 in dist.list) {
#     for(dist in dist.list) {
#       for (link0 in link.list) {
#         for (link in link.list) {
# 
#           a = ifelse(dist0 == "euclidean","e",ifelse(dist == "pearson","p","m"))
#           b = ifelse(dist == "euclidean","e",ifelse(dist0 == "pearson","p","m"))
#           c = ifelse(link0 == "ward.D","w",ifelse(link == "ward.D2","w2","a"))
#           d = ifelse(link == "ward.D","w",ifelse(link0 == "ward.D2","w2","a"))
# 
#           cat(paste(a,b,c,d,"\n",sep = "_"))
#           map.res.path <- file.path(res.path, paste0("CC_combat.FPKM_20RNAseq_",a,b,c,d,"_ClusterNum",clusterNum))
#           featType <- paste0("combat.FPKM_20RNAseq_",clusterNum,a,b,c,d)
# 
#           cluster.UTUC.20RNAseq.ans <- plot.common.cluster(indata,
#                                                         tumorname="UTUC",
#                                                         N.cluster=clusterNum,
#                                                         N.bootstrap=N.bootstrap,
#                                                         N.gene.per.bootstrap=N.gene.per.bootstrap,
#                                                         N.sample.per.bootstrap=N.sample.per.bootstrap,
#                                                         map.res.path=map.res.path, fig.path=fig.path,
#                                                         featType=featType,
#                                                         annCol=annCol.rna[colnames(indata),], annColors=annColors.rna,
#                                                         seed=123456, dist0=dist0, dist=dist,link0 = link0,link = link,
#                                                         clstCol=rainbow(clusterNum),
#                                                         namecvt=NULL, height = 7, fontsize=6, labRow = F, labCol = T,dendsort = T, cexCol = 0.1)
#           tmp <- data.frame("Type"=annCol.rna$Type,"Cluster"=cluster.UTUC.20RNAseq.ans$group[rownames(annCol.rna)],stringsAsFactors = F)
#           res.list[[featType]] <- table(tmp$Type,tmp$Cluster)
#         }
#       }
#     }
#   }
# }

map.res.path <- file.path(res.path, "CC_combat.FPKM_20RNAseq_mmww_ClusterNum3")
featType <- "combat.FPKM_20RNAseq"

cluster.UTUC.20RNAseq.ans <- plot.common.cluster(indata, 
                                                 tumorname="UTUC2", 
                                                 N.cluster=N.cluster, 
                                                 N.bootstrap=N.bootstrap, 
                                                 N.gene.per.bootstrap=N.gene.per.bootstrap, 
                                                 N.sample.per.bootstrap=N.sample.per.bootstrap, 
                                                 map.res.path=map.res.path, fig.path=fig.path, 
                                                 featType=featType,
                                                 annCol=annCol.rna[colnames(indata),], annColors=annColors.rna, 
                                                 seed=123456, dist0="manhattan", dist="manhattan",link0 = "ward.D",link = "ward.D", 
                                                 clstCol=c(blue,red,green), 
                                                 namecvt=NULL, height = 17, fontsize=5, labRow = F, labCol = T,dendsort = T, cexCol = 1.3,cexAnn = 0.8)

tmp <- data.frame(cluster = cluster.UTUC.20RNAseq.ans$group, feature = ifelse(cluster.UTUC.20RNAseq.ans$group == "cluster2","Muscle_Enriched","Non-muscle_Enriched"))
write.table(tmp,file.path(res.path,"cluster.assignment_by.20RNAseq.clustering_tumor.samples_UTUC.txt"),row.names = T,col.names = NA,sep = "\t")
tmp <- cluster.UTUC.20RNAseq.ans$group
group <- ifelse(tmp == "cluster2","Muscle_Enriched","Non-muscle_Enriched"); names(group) <- samples.UTUC
complist <- createList.muscle(group=group)
group <- data.frame("group"=group,"batch"=ifelse(grepl("rna",names(group)),"UTUCrna","UTUCr"),row.names = names(group))
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[rowids, samples.UTUC], tailrows, Groupinfo=group, features=Mids, featType="mRNA_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

DEFile <- file.path(res.path,"mRNA_rmBatch_edgeR_test_result.Muscle_Enriched_vs_Non-muscle_Enriched.txt")
figfile <- file.path(fig.path,"RNAcluster edgeR test volcano plot.pdf")
plotvolcano(DEFile = DEFile,"expression",figfile = figfile,pcut = 0.25,logfccut = 2)

index <- sapply(strsplit(rownames(group),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
rownames(group) <- index
annCol.rna$RNAcluster <- group[rownames(annCol.rna),"group"]
annColors.rna[["RNAcluster"]] <- c("Muscle_Enriched"=cherry,"Non-muscle_Enriched"=nake)

tmp <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.Muscle_Enriched_vs_Non-muscle_Enriched.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
#MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.DEGs.Muscle_Enriched_vs_Nonmuscle_Enriched <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res <- data.frame(GSEA.DEGs.Muscle_Enriched_vs_Nonmuscle_Enriched)
write.table(as.data.frame(GSEA.DEGs.Muscle_Enriched_vs_Nonmuscle_Enriched),file.path(res.path,"GSEA results for Muscle_Enriched_vs_Nonmuscle_Enriched DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

rowids <- rownames(tmp[which(tmp$FDR < 0.25 & abs(tmp$logFC) >= 2),])
rowids <- rownames(Ginfo[which(Ginfo$genename %in% rowids),])

indata <- combat.UTUC.FPKM[rowids,]
index <- sapply(strsplit(colnames(indata),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
colnames(indata) <- index

plotdata <- standarize.fun(indata,halfwidth = 1)
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D")
pdf(file = file.path(fig.path,"heatmap of CC_combat.FPKM_20RNAseq.pdf"),height = 5.5,width = 6)
mycol <- colorpanel(64,low=blue,mid = "black",high=gold)
aheatmap(as.matrix(plotdata), 
         Rowv=dendsort(as.dendrogram(hcg)), 
         Colv=cluster.UTUC.20RNAseq.ans$dendro, 
         annCol=annCol.rna[colnames(plotdata),c(12:11,15,17,24,21)],
         annColors = annColors.rna, 
         color = mycol,
         #color=blueyellow(64), 
         revC=TRUE, fontsize=7,cexAnn = 1,cexCol = 1,labRow = NA)
invisible(dev.off())

# pdf(file = file.path(fig.path,"heatmap of CC_combat.FPKM_20RNAseq_for annotation.pdf"),height = 15,width = 15)
# aheatmap(as.matrix(plotdata), 
#          Rowv=dendsort(as.dendrogram(hcg)), 
#          Colv=cluster.UTUC.20RNAseq.ans$dendro, 
#          annCol=annCol.rna[colnames(plotdata),c(1:12,15,16,24,28,25,21)],
#          annColors = annColors.rna, 
#          color=blueyellow(64), 
#          revC=TRUE, fontsize=8,cexAnn = 1,cexCol = 1.8,labRow = NA)
# invisible(dev.off())

tmp1 <- annCol.rna[which(annCol.rna$RNAcluster == "Muscle_Enriched"),"ImmuneScore"]
tmp2 <- annCol.rna[which(annCol.rna$RNAcluster == "Non-muscle_Enriched"),"ImmuneScore"]
p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.45
df <- data.frame("ImmuneScore"=c(tmp1,tmp2),"RNAcluster"=rep(c("MIeC","NMIeC"),c(length(tmp1),length(tmp2))))
p1 <- ggplot(df,aes(x=RNAcluster,y=ImmuneScore,fill=RNAcluster)) + geom_boxplot() + scale_fill_manual(values = c(cherry,nake)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")

tmp1 <- annCol.rna[which(annCol.rna$RNAcluster == "Muscle_Enriched"),"StromalScore"]
tmp2 <- annCol.rna[which(annCol.rna$RNAcluster == "Non-muscle_Enriched"),"StromalScore"]
p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.45
df <- data.frame("StromalScore"=c(tmp1,tmp2),"RNAcluster"=rep(c("MIeC","NMIeC"),c(length(tmp1),length(tmp2))))
p2 <- ggplot(df,aes(x=RNAcluster,y=StromalScore,fill=RNAcluster)) + geom_boxplot() + scale_fill_manual(values = c(cherry,nake)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")
p <- plot_grid(p1, p2, labels = c("a", "b"))
ggsave(file.path(fig.path,"boxplot of Immune and stromal score between RNAcluster C1 and C2 .pdf"),width = 5,height = 4)

#res <- as.data.frame(GSEA.DEGs.Muscle_Enriched_vs_Nonmuscle_Enriched)
#gsea.label <- res[grep("IMMU",res$Description),]
#gsea.label <- gsea.label[order(gsea.label$NES,decreasing = T),]
#gsea.label <- gsea.label[1:9,"Description"]
gsea.label <- c("NanoImmune_Cell_Functions",
                "NanoImmune_Regulation",
                "NanoImmune_B-Cell_Functions",
                "NanoImmune_T-Cell_Functions",
                "NanoImmune_Chemokines",
                "NanoImmune_Cytokines",
                "GO_ADAPTIVE_IMMUNE_RESPONSE",
                "KEGG_PRIMARY_IMMUNODEFICIENCY",
                "GO_LYMPHOCYTE_MEDIATED_IMMUNITY")
gsea.plot <- list()
for (i in gsea.label) {
  gsea.plot[[i]] <- gseaplot2(GSEA.DEGs.Muscle_Enriched_vs_Nonmuscle_Enriched,i,color = darkgreen,pvalue_table = T)
}
p <- plot_grid(gsea.plot[[1]], gsea.plot[[2]], gsea.plot[[3]], gsea.plot[[4]], gsea.plot[[5]], gsea.plot[[6]], gsea.plot[[7]], gsea.plot[[8]], gsea.plot[[9]], labels = c("a", "b", "c", "d", "e" ,"f", "g", "h", "i"),nrow = 3)
ggsave(filename = file.path(fig.path,"gseaplot for immune-realted pathway between RNAcluster C1 and C2_gseaplot.pdf"),width = 22,height = 20)

# Type
table(annCol.rna$Type,annCol.rna$RNAcluster)
fisher.test(matrix(c(10,5,1,4),byrow = T,ncol = 2)) #0.1273
# Muscle_Enriched Non-muscle_Enriched
# muscle-invasive                  10                   5
# Non-muscle invasive               1                   4

# grade
table(annCol.rna$Pathological_Stage,annCol.rna$RNAcluster)
fisher.test(matrix(c(1,4,10,5),byrow = T,ncol = 2)) #0.1273

# age
table(annCol.rna$Age,annCol.rna$RNAcluster)
fisher.test(matrix(c(6,2,5,7),byrow = T,ncol = 2)) #0.1968

# sex
table(annCol.rna$Gender,annCol.rna$RNAcluster)
fisher.test(matrix(c(1,1,10,8),byrow = T,ncol = 2)) #1

# ARID1A mut
table(annCol.rna$ARID1A_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(6,5,3,0),byrow = T,ncol = 2)) #0.2582

# EPHA5 ZFP36L1 mut
table(annCol.rna$EPHA5_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(9,3,0,2),byrow = T,ncol = 2)) #0.1099

# FGFR3 mut
table(annCol.rna$FGFR3_Mut_add,annCol.rna$RNAcluster)
fisher.test(matrix(c(6,6,5,3),byrow = T,ncol = 2)) #0.6699

# GABRA5 STAG2 mut
table(annCol.rna$GABRA5_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(8,5,1,0),byrow = T,ncol = 2)) #1

# GANAB KDM6A mut
table(annCol.rna$GANAB_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(8,4,1,1),byrow = T,ncol = 2)) #1

# GLYR1 mut
table(annCol.rna$GLYR1_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(7,5,2,0),byrow = T,ncol = 2)) #0.5055

# MLL2 mut
table(annCol.rna$MLL2_Mut,annCol.rna$RNAcluster)
fisher.test(matrix(c(7,3,2,2),byrow = T,ncol = 2)) #0.5804

###################################################################################
#### Differentially expression analysis for ZFP36L1 and CELSR1 mutation status ####
tmp <- UTUC.annotation[which(UTUC.annotation$`Exome-data`== "yes" &UTUC.annotation$`RNAseq (please add)` == "yes"),]
samples.exome.rna <- intersect(colnames(exome.mut.all),rownames(tmp))

#ZFP36L1
ZFP36L1 <- as.character(exome.mut.all["ZFP36L1",samples.exome.rna])
ZFP36L1 <- ifelse(ZFP36L1 == "1","Mutation","Normal"); names(ZFP36L1) <- tmp[samples.exome.rna,"DetailID"]
complist <- createList.mut(group=ZFP36L1)
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
#twoclasscomp(res.path, Ginfo, countsTable=countsTable[, names(ZFP36L1)], tailrows, features=Mids, featType="mRNA_ZFP36L1", complist, PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

ZFP36L1 <- data.frame("group"=ZFP36L1,"batch"=ifelse(grepl("rna",names(ZFP36L1)),"UTUCrna","UTUCr"),row.names = names(ZFP36L1))
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, rownames(ZFP36L1)], tailrows, Groupinfo=ZFP36L1, features=Mids, featType="mRNA_ZFP36L1_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

#ZFP36_family
ZFP36 <- exome.mut.all[c("ZFP36L1","ZFP36L2","ZFP36"),samples.exome.rna]
ZFP36 <- as.character(colSums(ZFP36))
ZFP36 <- ifelse(ZFP36 == "1","Mutation","Normal"); names(ZFP36) <- tmp[samples.exome.rna,"DetailID"]
complist <- createList.mut(group=ZFP36)
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
#twoclasscomp(res.path, Ginfo, countsTable=countsTable[, names(ZFP36)], tailrows, features=Mids, featType="mRNA_ZFP36", complist, PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

ZFP36 <- data.frame("group"=ZFP36,"batch"=ifelse(grepl("rna",names(ZFP36)),"UTUCrna","UTUCr"),row.names = names(ZFP36))
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, rownames(ZFP36)], tailrows, Groupinfo=ZFP36, features=Mids, featType="mRNA_ZFP36_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

# BatchInfo <- as.character(ZFP36$batch); names(BatchInfo) <- rownames(ZFP36)
# twoclassDESeq2(res.path, Ginfo, countsTable=countsTable[, rownames(ZFP36)], tailrows, Batchinfo=BatchInfo, features=Mids, featType="mRNA_ZFP36_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

DEFile <- file.path(res.path,"mRNA_ZFP36_rmBatch_edgeR_test_result.Mutation_vs_Normal.txt")
figfile <- file.path(fig.path,"ZFP36 edgeR test volcano plot.pdf")
plotvolcano(DEFile = DEFile,"expression",figfile = figfile,pcut = 0.25,logfccut = 1.5)
#CELSR1
# CELSR1 <- as.character(exome.mut.all["CELSR1",samples.exome.rna])
# CELSR1 <- ifelse(CELSR1 == "1","Mutation","Normal"); names(CELSR1) <- tmp[samples.exome.rna,"DetailID"]
# complist <- createList.mut(group=CELSR1)
# rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
# twoclasscomp(res.path, Ginfo, countsTable=countsTable[, names(CELSR1)], tailrows, features=Mids, featType="mRNA_CELSR1", complist, PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)
# 
# CELSR1 <- data.frame("group"=CELSR1,"batch"=ifelse(grepl("rna",names(CELSR1)),"UTUCrna","UTUCr"),row.names = names(CELSR1))
# twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, rownames(CELSR1)], tailrows, Groupinfo=CELSR1, features=Mids, featType="mRNA_CELSR1_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

############################################
# perform GSEA for ZFP36 family edgeR DEGs #
tmp <- read.table(file.path(res.path,"mRNA_ZFP36_rmBatch_edgeR_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.edgeR.ZFP36 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res <- data.frame(GSEA.edgeR.ZFP36)
write.table(as.data.frame(GSEA.edgeR.ZFP36),file.path(res.path,"GSEA results for ZFP63 family Mutation vs Normal DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

gsea.label <- c("NanoImmune_Cell_Functions",
                "NanoImmune_Regulation",
                "NanoImmune_B-Cell_Functions",
                "NanoImmune_T-Cell_Functions",
                "NanoImmune_Chemokines",
                "NanoImmune_Cytokines",
                "GO_ADAPTIVE_IMMUNE_RESPONSE",
                "KEGG_PRIMARY_IMMUNODEFICIENCY",
                "GO_LYMPHOCYTE_MEDIATED_IMMUNITY")
gsea.plot <- list()
for (i in gsea.label) {
  gsea.plot[[i]] <- gseaplot2(GSEA.edgeR.ZFP36,i,color = darkgreen,pvalue_table = T)
}
p <- plot_grid(gsea.plot[[1]], gsea.plot[[2]], gsea.plot[[3]], gsea.plot[[4]], gsea.plot[[5]], gsea.plot[[6]], gsea.plot[[7]], gsea.plot[[8]], gsea.plot[[9]], labels = c("a", "b", "c", "d", "e" ,"f", "g", "h", "i"),nrow = 3)
ggsave(filename = file.path(fig.path,"gseaplot for immune-realted pathway in ZFP36 mutation_gseaplot.pdf"),width = 22,height = 20)

#############################################
# perform GSEA for ZFP36 family DESeq2 DEGs #
tmp <- read.table(file.path(res.path,"mRNA_ZFP36_rmBatch_deseq2_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$log2FoldChange
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.deseq2.ZFP36 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res <- data.frame(GSEA.deseq2.ZFP36)
write.table(as.data.frame(GSEA.deseq2.ZFP36),file.path(res.path,"GSEA results for ZFP63 family Mutation vs Normal DEGs of DESeq2 by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

##################################################################
# test immune score difference bewteen ZFP36 mutation and normal #

# no significant but GSEA show enrichement in immunity
tmp <- annCol.rna[samples.exome.rna,"ImmuneScore"]; names(tmp) <- samples.exome.rna
wilcox.test(tmp[c("29T","31T","02T")],tmp[setdiff(names(tmp),c("29T","31T","02T"))])

p <- c()
for (i in 1:12) {
  tmp <- annCol.rna[samples.exome.rna,i]; names(tmp) <- samples.exome.rna
  tmp <- wilcox.test(tmp[c("29T","31T","02T")],tmp[setdiff(names(tmp),c("29T","31T","02T"))])
  p <- c(p,tmp$p.value)
}

#######################################################
### check BLCA ZFP36 family mutation and MCPcounter ###
blca.mut <- read.table(file.path(res.path,"BLCA_TCGA_ZFP36_Family_Mutation.tsv"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
colnames(blca.mut) <- paste0("BLCA",substr(colnames(blca.mut),start = 8,stop = 15))

tmp <- substr(colnames(FPKM.BLCA),1,12)
tmp1 <- colnames(blca.mut)
tmp2 <- intersect(tmp,tmp1)
index <- match(tmp2,tmp)
blca.mut <- blca.mut[,tmp2]
colnames(blca.mut) <- colnames(FPKM.BLCA)[index]

blca.mut[is.na(blca.mut)] <- "0"
blca.mut[blca.mut == ""] <- "0"
blca.mut[blca.mut != "0"] <- "1"
blca.mut <- as.data.frame(sapply(blca.mut, as.numeric)); rownames(blca.mut) <- c("ZFP36","ZFP36L1","ZFP36L2")

# ZFP36 family
blca.ZFP36 <- colSums(blca.mut[c("ZFP36","ZFP36L1","ZFP36L2"),])
blca.ZFP36 <- ifelse(blca.ZFP36 >= 1,1,0)
blca.ZFP36.mutsam <- intersect(names(which(blca.ZFP36 == 1)),colnames(FPKM.BLCA))
blca.ZFP36.norsam <- intersect(names(which(blca.ZFP36 == 0)),colnames(FPKM.BLCA))
MCPscore.BLCA <- MCPcounter.estimate(expression = FPKM.BLCA,featuresType = "HUGO_symbols")
#MCPscore.BLCA <- log2(MCPscore.BLCA + 1)
#MCPscore.BLCA <- annTrackScale(indata = MCPscore.BLCA, halfwidth = 2, poolsd = F)
write.table(MCPscore.BLCA,file.path(res.path,"MCPscore.BLCA.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

table(Sinfo[which(rownames(Sinfo) %in% blca.ZFP36.mutsam),"TCGA_Subtype"])
p <- c()
for (i in 1:10) {
  tmp <- MCPscore.BLCA[i,]
  tmp <- wilcox.test(tmp[blca.ZFP36.mutsam],tmp[blca.ZFP36.norsam])
  p <- c(p,tmp$p.value)
}
names(p) <- rownames(MCPscore.BLCA)

blca.ZFP36.mut <- blca.ZFP36[c(blca.ZFP36.mutsam,blca.ZFP36.norsam)]; blca.ZFP36.mut <- ifelse(blca.ZFP36.mut == 1,"Mutation","Normal")
complist <- createList.mut(group=blca.ZFP36.mut)
blca.ZFP36.mut <- data.frame("group"=blca.ZFP36.mut,samples=names(blca.ZFP36.mut),row.names = names(blca.ZFP36.mut))
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, rownames(blca.ZFP36.mut)], tailrows, Groupinfo=blca.ZFP36.mut, features=Mids, featType="mRNA_ZFP36_family_BLCA", complist,PASSFlag=PASSFlag.mRNA, overwt=TRUE)
# twoclassDESeq2(res.path, Ginfo, countsTable=countsTable[, rownames(blca.ZFP36.mut)], tailrows, Batchinfo=NULL, features=Mids, featType="mRNA_ZFP36_family_BLCA", complist,PASSFlag=PASSFlag.mRNA, overwt=TRUE)
# twoclasscomp(res.path, Ginfo, countsTable=countsTable[, rownames(blca.ZFP36.mut)], tailrows, features=Mids, featType="mRNA_ZFP36_BLCA", complist, PASSFlag=PASSFlag.mRNA, overwt=TRUE)
# 
# tmp <- read.table(file.path(res.path,"mRNA_ZFP36_BLCA_deseq_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
# geneList <- tmp$foldChange
# names(geneList) <- rownames(tmp)
# geneList <- sort(geneList,decreasing = T)
# #MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
# GSEA.DESeq.ZFP36.BLCA <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
# res <- data.frame(GSEA.DESeq.ZFP36.BLCA)
# write.table(as.data.frame(GSEA.DESeq.ZFP36.BLCA),file.path(res.path,"GSEA results for BLCA ZFP63 family Mutation vs Normal DEGs of DESeq by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

tmp <- read.table(file.path(res.path,"mRNA_ZFP36_family_BLCA_edgeR_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
#MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.edgeR.ZFP36.BLCA <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F,pvalueCutoff = 0.25)
res <- data.frame(GSEA.edgeR.ZFP36.BLCA)
write.table(as.data.frame(GSEA.edgeR.ZFP36.BLCA),file.path(res.path,"GSEA results for BLCA ZFP63 family Mutation vs Normal DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
immune.id <- c("NanoImmune_T-Cell_Functions",
               "NanoImmune_Cytokines",
               "NanoImmune_Chemokines",
               "NanoImmune_Interleukins",
               "GO_INNATE_IMMUNE_RESPONSE")
gseaplot2(GSEA.edgeR.ZFP36.BLCA, geneSetID = immune.id,color = mycol[1:length(immune.id)],pvalue_table = F)
ggsave(file.path(fig.path,"GSEA immune id for BLCA ZFP63 family Mutation vs Normal.pdf"),width = 5.5,height = 5)

# ZFP36L1
blca.ZFP36L1 <- colSums(blca.mut["ZFP36L1",])
blca.ZFP36L1.mutsam <- intersect(names(which(blca.ZFP36L1 == 1)),colnames(FPKM.BLCA))
blca.ZFP36L1.norsam <- intersect(names(which(blca.ZFP36L1 == 0)),colnames(FPKM.BLCA))
blca.ZFP36L1.mut <- blca.ZFP36L1[c(blca.ZFP36L1.mutsam,blca.ZFP36L1.norsam)]; blca.ZFP36L1.mut <- ifelse(blca.ZFP36L1.mut == 1,"Mutation","Normal")
complist <- createList.mut(group=blca.ZFP36L1.mut)
blca.ZFP36L1.mut <- data.frame("group"=blca.ZFP36L1.mut,samples=names(blca.ZFP36L1.mut),row.names = names(blca.ZFP36L1.mut))
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, rownames(blca.ZFP36L1.mut)], tailrows, Groupinfo=blca.ZFP36L1.mut, features=Mids, featType="mRNA_ZFP36L1_BLCA2", complist,PASSFlag=PASSFlag.mRNA, overwt=TRUE)

tmp <- read.table(file.path(res.path,"mRNA_ZFP36L1_BLCA_edgeR_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
#MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.edgeR.ZFP36L1.BLCA <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F,pvalueCutoff = 0.25)
res <- data.frame(GSEA.edgeR.ZFP36L1.BLCA)
write.table(as.data.frame(GSEA.edgeR.ZFP36L1.BLCA),file.path(res.path,"GSEA results for BLCA ZFP63L1 Mutation vs Normal DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

#APOBEC3A expression between BLCA ZFP36L1 mutation and Normal

tmp1 <- log2(countsNorm[rownames(Ginfo[which(Ginfo$genename == "APOBEC3A"),]),blca.ZFP36L1.mutsam] + 1)
tmp2 <- log2(countsNorm[rownames(Ginfo[which(Ginfo$genename == "APOBEC3A"),]),blca.ZFP36L1.norsam] + 1)

std <- function(x) sd(x)/sqrt(length(x))

avg1 <- mean(as.numeric(tmp1))
avg2 <- mean(as.numeric(tmp2))

std1 <- std(as.numeric(tmp1))
std2 <- std(as.numeric(tmp2))

avg <- c(avg1,avg2)
var <- c(std1,std2)

pdf(file.path(fig.path,"log2 normalized counts of APOBEC3A in TCGA_BLCA between ZFP36L1 mutation and normal.pdf"),width = 3,height = 6)
par(bty="o", mgp = c(2,0.5,0), mar = c(7.1,4.1,2.1,4.1),tcl=-.25)
par(xpd = T)
bar <- barplot(avg, # save object to get x-coordinate
               border = F,
               ylab = "log2 normalized counts of APOBEC3A in TCGA-BLCA",
               #border = c(T,F,T,T,T,F,T,T,T,T,T),
               ylim = c(0,6),
               yaxt = "n",
               xaxt = "n",
               col = c(cherry,lightgrey))
axis(side = 2,at = seq(0,6,1))

# get vertical error bar
segments(bar, avg - var * 2, bar,
         avg + var * 2, lwd = 1.5)

# add horizon line to error bar
arrows(bar, avg - var * 2, bar,
       avg + var * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

# add x label, x postion is bar object
text(x = bar, y = par("usr")[3] - 0.05, srt = 45,
     adj = 1, labels = c("ZFP36L1 mutated","ZFP36L1 wild"), xpd = TRUE, cex = 1)
invisible(dev.off())

#############################################################
######### Pie chart for mutation type distribution ##########
tmp <- read.table(file.path(data.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_wideformat.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T)
tmp <- tmp[-which(tmp$Variation_Classification %in% "unknown"),]
tmp <- as.data.frame(table(tmp$Variation_Classification))
mygene[mygene == "nonsynonymous SNV"] <- "Nonsynonymous"
mygene[mygene == "stopgain"] <- "Stopgain"
mygene[mygene == "frameshift deletion"] <- "Frameshift Deletion"
mygene[mygene == "frameshift insertion"] <- "Frameshift Insertion"
mygene[mygene == "nonframeshift deletion"] <- "Non-frameshift Deleletion"
mygene[mygene == "nonframeshift insertion"] <- "Non-frameshift Insertion"
#mygene[mygene == "multi-hit"] <- "Multi_Hit"
mygene[mygene == "splicing"] <- "Splicing"

df <- tmp$Freq; names(df) <- c("Frameshift Deletion","Frameshift Insertion","Non-frameshift Deletion","Non-frameshift Insertion","Nonsynonymous","Splicing","Stopgain","Stoploss","Synonymous")
df <- df[c("Nonsynonymous","Frameshift Deletion","Frameshift Insertion","Non-frameshift Deletion","Non-frameshift Insertion","Splicing","Stopgain","Stoploss","Synonymous")]
pdf(file.path(fig.path,"pie chart of mutation type distribution.pdf"),width = 9,height = 9)
pie(df, cex=1.3, col=colorRampPalette(brewer.pal(11,'Spectral'))(9), border="white",  radius=0.75, cex.main=0.7, labels=ifelse(df/sum(df) > 0.01, paste("",names(df), "(" ,round(df/sum(df),2)*100, ")",sep=""),""), init.angle=90)
symbols(0,0,circles=.25, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
invisible(dev.off())
#title(main = "Distribution of mutation type",cex.main=0.95)

######################################################
### mutation load and permegabase in UTUC and BLCA ###

# UTUC
# median mutatonal load per exome all mutation
tmp <- read.table(file.path(data.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_wideformat.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T)
mutload <- c()
for (i in 17:46) {
  ml <- tmp[,c(i,47)]
  #ml$Variation_Classification <- ifelse(ml$Variation_Classification == "synonymous SNV","silence","nonsilence")
  #ml <- as.data.frame(table(ml[,1],ml$Variation_Classification))
  #ml <- ml[which(ml$Var1 == "1" & ml$Var2 == "nonsilence"),"Freq"]
  ml <- ml[which(ml$Variation_Classification !="unknown"),]
  ml <- as.data.frame(table(ml[,1],ml$Variation_Classification))
  ml <- sum(ml[which(ml$Var1 == "1"),"Freq"])
  names(ml) <- colnames(tmp)[i]
  mutload <- c(mutload,ml)
}
range(mutload) #13-786
median.mutload <- median(mutload) #73
avg.mutload <- mean(mutload) # 141.8

median.mutload/50 #1.46
range(mutload)/50 #0.26-15.72

# average of all mutations per megabase
tmp <- read.table(file.path(data.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_wideformat.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T)
mutpermeg <- c()
for (i in 17:46) {
  ml <- round(sum(tmp[,i])/50,2) #megabase is 50
  names(ml) <- colnames(tmp)[i]
  mutpermeg <- c(mutpermeg,ml)
}
range(mutpermeg) # 0.26-15.86
avg.mutpermeg <- mean(mutpermeg) # 2.867
sd.mutpermeg <- sd(mutpermeg) #3.450

utuc.ZFP36.mutsam <- c("UTUC_23T","UTUC_28T","UTUC_29T","UTUC_2T","UTUC_30T","UTUC_31T","UTUC_34T","UTUC_6T")
utuc.ZFP36.norsam <- setdiff(names(mutpermeg),utuc.ZFP36.mutsam)

wilcox.test(mutpermeg[utuc.ZFP36.mutsam],mutpermeg[utuc.ZFP36.norsam],alternative = "greater") # 0.007857
wilcox.test(mutload[utuc.ZFP36.mutsam],mutload[utuc.ZFP36.norsam],alternative = "greater") # 0.008377

# df <- data.frame("MB"=c(mutpermeg[utuc.ZFP36.mutsam],mutpermeg[utuc.ZFP36.norsam]),"ZFP36_Family_Mutation"=rep(c("Presence","Absence"),c(8,22)))
# df$ZFP36_Family_Mutation <- factor(df$ZFP36_Family_Mutation,levels = c("Presence","Absence"))
# p1 <- ggplot(df,aes(x=ZFP36_Family_Mutation,y=MB,fill=ZFP36_Family_Mutation)) + geom_boxplot() + scale_fill_manual(values = c(cherry,lightgrey)) + ggtitle("pValue = 0.02") + theme(legend.position="none") + labs(y = "Mutation/MB",x = "ZFP36 Family Mutation")

df <- data.frame("Load"=c(mutload[utuc.ZFP36.mutsam],mutload[utuc.ZFP36.norsam]),"ZFP36_Family_Mutation"=rep(c("Presence","Absence"),c(8,22)))
df$ZFP36_Family_Mutation <- factor(df$ZFP36_Family_Mutation,levels = c("Presence","Absence"))
p1 <- ggplot(df,aes(x=ZFP36_Family_Mutation,y=Load,fill=ZFP36_Family_Mutation)) + geom_boxplot() + scale_fill_manual(values = c(cherry,lightgrey)) + ggtitle("UTUC") + theme(legend.position="none") + labs(y = "Mutation Load",x = "ZFP36 Family Mutation") + ylim(0,400) + annotate("text",x=1.5,y=380,label="p = 0.008",size=3)

# mutation load between NMI and MI
tmp1 <- as.numeric(na.omit(mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$Type == "muscle-invasive"),"ID"]),"T")]))
tmp2 <- as.numeric(na.omit(mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),"ID"]),"T")]))
wilcox.test(tmp1,tmp2,alternative = "greater",rm.na=T) # 0.06241
df <- data.frame("Load"=c(tmp1,tmp2),"Invasiveness"=rep(c("MI","NMI"),c(15,15)))
df$Invasiveness <- factor(df$Invasiveness,levels = c("MI","NMI"))
p <- ggplot(df,aes(x=Invasiveness,y=Load,fill=Invasiveness)) + geom_boxplot() + scale_fill_manual(values = c(soil,lightgrey)) + ggtitle("UTUC") + theme(legend.position="none") + labs(y = "Mutation Load",x = "Invasiveness") + ylim(0,400) + annotate("text",x=1.5,y=380,label="p = 0.062",size=3)
ggsave(file.path(fig.path,"Mutation load between invasiveness in UTUC.pdf"),width = 2,height = 3.5)


# BLCA
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)

blca.mut.maf <- read_tsv(file.path(data.path,"BLCA_Mutation_maf_cBioportal.txt"), comment = "#")
blca.mut.maf <- as.data.frame(blca.mut.maf)
blca.mut.maf$Tumor_Sample_Barcode <- paste0("BLCA",substr(blca.mut.maf$Tumor_Sample_Barcode,start = 8,stop = 15))
mutect.dataframe <- function(x){
  # delete rows of Silent
  #cut_id <- x$Variant_Classification == "Silent"
  #x <- x[!cut_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}
mutload.blca <- as.data.frame(mutect.dataframe(blca.mut.maf)); rownames(mutload.blca) <- mutload.blca$Tumor_Sample_Barcode
dim(mutload.blca)
head(mutload.blca)

wilcox.test(mutload.blca[substr(blca.ZFP36.mutsam,1,12),"TCGA_sum"],mutload.blca[substr(blca.ZFP36.norsam,1,12),"TCGA_sum"],alternative = "greater") # 0.0009173
df <- data.frame("Load"=c(mutload.blca[substr(blca.ZFP36.mutsam,1,12),"TCGA_sum"],mutload.blca[substr(blca.ZFP36.norsam,1,12),"TCGA_sum"]),"ZFP36_Family_Mutation"=rep(c("Presence","Absence"),c(39,368)))
df$ZFP36_Family_Mutation <- factor(df$ZFP36_Family_Mutation,levels = c("Presence","Absence"))
p2 <- ggplot(df,aes(x=ZFP36_Family_Mutation,y=Load,fill=ZFP36_Family_Mutation)) + geom_boxplot() + scale_fill_manual(values = c(cherry,lightgrey)) + ggtitle("BLCA") + theme(legend.position="none") + labs(y = "Mutation Load",x = "ZFP36 Family Mutation") + ylim(0,400) + annotate("text",x=1.5,y=380,label="p = 0.0009",size=3)
p3 <- plot_grid(p1, p2, labels = c("a", "b"))
ggsave(file.path(fig.path,"Mutation load between UTUC and BLCA ZFP36 Family mutation.pdf"),width = 5.5,height = 4)

#######################
# CNA 9p21.3 deletion #
cna.meth.path <- "F:/Project/UTUC_BLCA/CNA/GISTIC/GISTIC_del0.4_amp1.7_qvalue0.0001"
cna.meth <- read.table(file.path(cna.meth.path,"CNA_segment_forGISTIC2.0.all_lesions.conf_90.txt"),sep = "\t",check.names = F,header = T,row.names = NULL,stringsAsFactors = F)
del9p <- as.data.frame(t(cna.meth[which(cna.meth$`Unique Name` == "Deletion Peak 2"),10:44]))
rownames(del9p) <- substr(rownames(del9p),start = 6,stop = 8);colnames(del9p) <- "del9p"

cna.exome.path <- "F:/Project/UTUC_BLCA/CNA"
cna.exome <- read.table(file.path(cna.exome.path,"UTUC.rm.all_lesions.conf_75.txt"),sep = "\t",check.names = F,header = T,row.names = NULL,stringsAsFactors = F)
del9p.exome <- as.data.frame(t(cna.exome[9,10:38])); colnames(del9p.exome) <- "del9p"
rownames(del9p.exome) <- sapply(strsplit(rownames(del9p.exome),"_"),"[",2)
rownames(del9p.exome) <- paste0("0",rownames(del9p.exome))
rownames(del9p.exome) <- substr(rownames(del9p.exome),nchar(rownames(del9p.exome))-2,nchar(rownames(del9p.exome)))

#########################################################
### compare mutation signature in UTUC ZFP63 mutation ###
wt <- read.table("F:/Project/UTUC_BLCA/exome UTUC/Results/mutation.signature.weightMatrix.bydeconstructSigs.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
wt <- as.data.frame(t(wt[,-31]))

# no significance
p <- c()
for (i in rownames(wt)) {
  tmp <- wt[i,]
  tmp <- wilcox.test(as.numeric(tmp[utuc.ZFP36.mutsam]),as.numeric(tmp[utuc.ZFP36.norsam]))$p.value; names(tmp) <- i
  p <- c(p,tmp)
}

####################################
### coexpression of UTUC ZFP36L1 ###
rowids <- intersect(Mids,names(PASSFlag.mRNA[which(PASSFlag.mRNA == T)]))
rowids <- intersect(unique(Ginfo[rowids,"genename"]),rownames(combat.UTUC.FPKM.HUGO))

cor.mat <- NULL
indata <- combat.UTUC.FPKM.HUGO[rowids,]
a <- log2(as.numeric(indata["ZFP36L1",]) + 1)
for (i in 1:nrow(indata)) {
  if(rownames(indata)[i] != "ZFP36L1"){
    b <- log2(as.numeric(indata[i,]) + 1)
    #tmp <- cor.test(a,b,method = "spearman")
    tmp <- cor.test(a,b)
    p <- tmp$p.value
    r <- tmp$estimate
    cor.mat <- rbind.data.frame(cor.mat,data.frame(geneA="ZFP36L1",geneB=rownames(indata)[i],r=r,p=p))
  }
}
utuc.cor.mat <- cor.mat
utuc.cor.mat$fdr <- p.adjust(utuc.cor.mat$p,method = "BH")
write.table(utuc.cor.mat,file.path(res.path,"UTUC pearson correlation matrix between lowfiltered gene and ZFP36L1.txt"),sep = "\t",row.names = F)
utuc.poscor.gene <- as.character(utuc.cor.mat[which(utuc.cor.mat$r > 0.25 & utuc.cor.mat$p < 0.05),"geneB"])
utuc.negcor.gene <- as.character(utuc.cor.mat[which(utuc.cor.mat$r < -0.25 & utuc.cor.mat$p < 0.05),"geneB"])

cor.mat <- NULL
indata <- FPKM.BLCA[rowids,]
a <- log2(as.numeric(indata["ZFP36L1",]) + 1)
for (i in 1:nrow(indata)) {
  if(rownames(indata)[i] != "ZFP36L1"){
    b <- log2(as.numeric(indata[i,]) + 1)
    #tmp <- cor.test(a,b,method = "spearman")
    tmp <- cor.test(a,b)
    p <- tmp$p.value
    r <- tmp$estimate
    cor.mat <- rbind.data.frame(cor.mat,data.frame(geneA="ZFP36L1",geneB=rownames(indata)[i],r=r,p=p))
  }
}
blca.cor.mat <- cor.mat
blca.cor.mat$fdr <- p.adjust(blca.cor.mat$p,method = "BH")
write.table(blca.cor.mat,file.path(res.path,"BLCA pearson correlation matrix between lowfiltered gene and ZFP36L1.txt"),sep = "\t",row.names = F)
blca.poscor.gene <- as.character(blca.cor.mat[which(blca.cor.mat$r > 0.25 & blca.cor.mat$p < 0.05),"geneB"])
blca.negcor.gene <- as.character(blca.cor.mat[which(blca.cor.mat$r < -0.25 & blca.cor.mat$p < 0.05),"geneB"])

cor.gene <- list(UTUCp_BLCAp=intersect(utuc.poscor.gene,blca.poscor.gene),
                 UTUCn_BLCAn=intersect(utuc.negcor.gene,blca.negcor.gene),
                 UTUCp_BLCAn=intersect(utuc.poscor.gene,blca.negcor.gene),
                 UTUCn_BLCAp=intersect(utuc.negcor.gene,blca.poscor.gene))
library(VennDiagram)
venn.diagram(cor.gene, filename=file.path(fig.path,"venn of UTUC and BLCA pearson coexpression gene with ZFP36L1.tiff"),imagetype = "tiff", lty=1, col=c(sun,blue,gold,green),fill=c(sun,blue,gold,green), alpha=0.05, euler.d=TRUE, fontfamily="Helvetica", cat.fontfamily="Helvetica", cat.fontface="italic",cat.cex=0.7)

test <- bitr(as.character(cor.gene[[1]]), fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = 'org.Hs.eg.db')
ego.MF <- enrichGO(gene          = as.character(test$ENSEMBL),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.25,
                   minGSSize     = 10,
                   readable      = T)
#ego.MF2 <- simplify(ego.MF)
MF <- data.frame(summary(ego.MF))
write.table(MF,file.path(res.path,"UTUCp_BLCAp MF clusterprofiler.txt"),sep = "\t",row.names = F)

ego.BP <- enrichGO(gene          = as.character(test$ENSEMBL),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.25,
                   minGSSize     = 10,
                   readable      = T)
#ego.BP2 <- simplify(ego.BP)
BP <- data.frame(summary(ego.BP))
write.table(BP,file.path(res.path,"UTUCp_BLCAp BP clusterprofiler.txt"),sep = "\t",row.names = F)

ego.CC <- enrichGO(gene          = as.character(test$ENSEMBL),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.25,
                   minGSSize     = 10,
                   readable      = T)
#ego.CC2 <- simplify(ego.CC)
CC <- data.frame(summary(ego.CC))
write.table(CC,file.path(res.path,"UTUCp_BLCAp CC clusterprofiler.txt"),sep = "\t",row.names = F)

#######################################
### immunity of BLCA ZFP63 mutation ###
immune.checkpoint <- c("PDCD1","CD274","PDCD1LG2","CD80",
                       "CD86","CD28","CTLA4","ICOSLG",
                       "ICOS","CD276","VTCN1","TNFRSF14",
                       "BTLA","CD160","TNFSF14","LAG3",
                       "TNFSF9","TNFRSF9","TNFSF4","TNFRSF4",
                       "CD70","CD27","CD40","CD40LG",
                       "LGALS9","HAVCR2","ADORA2A","TNFRSF18",
                       "TNFSF18","IDO1","C10ORF54","CEACAM1",
                       "CD47","SIRPA","CD226","PVR",
                       "TIGIT","CD244","CD48","TMIGD2",
                       "HHLA2","BTN2A1","CD209","BTN2A2",
                       "BTN3A1","BTNL3","BTNL9","CD96",
                       "TDO2","CD200","CD200R1")

indata <- FPKM.BLCA[immune.checkpoint,c(blca.ZFP36.mutsam,blca.ZFP36.norsam)]
plotdata <- log2(indata + 1)

p <- c()
for (i in 1:length(immune.checkpoint)) {
  tmp <- plotdata[i,]
  tmp <- wilcox.test(as.numeric(tmp[blca.ZFP36.mutsam]),as.numeric(tmp[blca.ZFP36.norsam]),alternative = "greater")
  p <- c(p,tmp$p.value)
}
names(p) <- immune.checkpoint

#plotdata <- sweep(plotdata,2,apply(plotdata, 2, median))
#plotdata <- sweep(plotdata,1,apply(plotdata, 1, median))
plotdata <- as.data.frame(na.omit(standarize.fun(indata = plotdata,halfwidth = 2)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")

lm=rbind(c(1,1,2,2,2))
mycol <- colorpanel(256,low=blue,mid = "black",high=gold)
p1 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.mutsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = F,show_colnames = F,legend = F)
p2 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.norsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = T,show_colnames = F,legend = T)

outFigFile <- "BLCA_heatmap_immune.checkpoint_51genes_.pdf"
pdf(file.path(fig.path, outFigFile), height=6)
grid.arrange(grobs = list(p1[[4]],p2[[4]]), layout_matrix = lm)
invisible(dev.off())

immune.signature <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig <- list()
for (i in cell.type) {
  immune.sig[[i]] <- toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"])
}
immune.sig[["Angiogenesis"]] <- c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")

immune.enrichscore<- gsva(as.matrix(FPKM.BLCA[,c(blca.ZFP36.mutsam,blca.ZFP36.norsam)]),immune.sig,method="ssgsea")
p <- vector()
for (i in 1:nrow(immune.enrichscore)) {
  tmp <- immune.enrichscore[i,]
  tmp <- wilcox.test(tmp[blca.ZFP36.mutsam],tmp[blca.ZFP36.norsam],alternative = "greater")$p.value
  p <- append(p, tmp)
}
names(p) <- names(immune.sig)
plotdata <- as.data.frame(na.omit(standarize.fun(indata = immune.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
p1 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.mutsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = F,show_colnames = F,legend = F)
p2 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.norsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = T,show_colnames = F,legend = T)
grid.arrange(grobs = list(p1[[4]],p2[[4]]), layout_matrix = lm)

# immune.sig.curated <- list()
# for (i in names(immune.sig)) {
#   tmp <- immune.sig[[i]]
#   tmp3 <- c()
#   for (j in tmp) {
#     if(is.element(j,rownames(FPKM.BLCA))){
#       tmp1 <- as.numeric(FPKM.BLCA[j,c(blca.ZFP36.mutsam,blca.ZFP36.norsam)])
#       tmp2 <- wilcox.test(tmp1[1:39],tmp1[40:407],alternative = "greater")$p.value
#       if(tmp2 < 0.05) {
#         tmp3 <- c(tmp3,j)
#       }else {next()}
#     }else {next()}
#   }
#   if(length(tmp3) > 0) {
#     immune.sig.curated[[i]] <- tmp3
#   }else {next()}
# }
# 
# immune.enrichscore.curated<- gsva(as.matrix(FPKM.BLCA[,c(blca.ZFP36.mutsam,blca.ZFP36.norsam)]),immune.sig.curated,method="ssgsea")
# 
# p <- vector()
# for (i in 1:nrow(immune.enrichscore.curated)) {
#   tmp <- wilcox.test(immune.enrichscore.curated[i,1:39],immune.enrichscore.curated[i,40:407],alternative = "greater")$p.value
#   p <- append(p, tmp)
# }
# names(p) <- names(immune.sig.curated)
# 
# plotdata <- as.data.frame(na.omit(standarize.fun(indata = immune.enrichscore.curated,halfwidth = 1)))
# hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
# p1 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.mutsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = F,show_colnames = F,legend = F)
# p2 <- pheatmap(as.matrix(plotdata[,blca.ZFP36.norsam]),color = mycol,cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = T,show_colnames = F,legend = T)
# grid.arrange(grobs = list(p1[[4]],p2[[4]]), layout_matrix = lm)

#########################
### metabolic pathway ###
# tcga.DMPsCluster.C1 <- intersect(rownames(Sinfo.BLCA[which(Sinfo.BLCA$DMPsCluster == "cluster1"),]),substr(colnames(FPKM.BLCA),1,12))
# tcga.DMPsCluster.C2 <- intersect(rownames(Sinfo.BLCA[which(Sinfo.BLCA$DMPsCluster == "cluster2"),]),substr(colnames(FPKM.BLCA),1,12))
# 
# meta.signature <- read.table(file.path(comAnn.path,"Metabolism_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# meta.type <- unique(meta.signature$`Metabolic Function`)
# meta.sig <- list()
# for (i in meta.type) {
#   meta.sig[[i]] <- toupper(meta.signature[which(meta.signature$`Metabolic Function` == i),"Gene Name"])
# }
# 
# meta.enrichscore<- gsva(as.matrix(FPKM.BLCA),meta.sig,method="ssgsea")
# colnames(meta.enrichscore) <- substr(colnames(meta.enrichscore),1,12)
# p <- vector()
# for (i in 1:nrow(meta.enrichscore)) {
#   tmp <- meta.enrichscore[i,]
#   #tmp <- wilcox.test(as.numeric(tmp[blca.ZFP36.mutsam]),as.numeric(tmp[blca.ZFP36.norsam]))$p.value
#   tmp <- round(wilcox.test(tmp[tcga.DMPsCluster.C1],tmp[tcga.DMPsCluster.C2])$p.value,3)
#   p <- append(p, tmp)
# }
# names(p) <- names(meta.sig)
# 
# plotdata <- as.data.frame(na.omit(standarize.fun(indata = meta.enrichscore[,c(tcga.DMPsCluster.C1,tcga.DMPsCluster.C2)],halfwidth = 1)))
# hcg <- hclust(distanceMatrix(as.matrix(t(plotdata[names(p[p<0.1]),])), "pearson"), "ward.D")
# p1 <- pheatmap(as.matrix(plotdata[,tcga.DMPsCluster.C1]),color = bluered(128),cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = F,show_colnames = F,legend = F)
# p2 <- pheatmap(as.matrix(plotdata[,tcga.DMPsCluster.C2]),color = bluered(128),cluster_rows = hcg,treeheight_row = 0,cluster_cols = F,show_rownames = T,show_colnames = F,legend = T)
# lm=rbind(c(1,1,2,2,2,2))
# grid.arrange(grobs = list(p1[[4]],p2[[4]]), layout_matrix = lm)
# 
# meta.enrichscore.utuc<- gsva(as.matrix(combat.UTUC.FPKM.HUGO),meta.sig,method="ssgsea")
# p <- vector()
# for (i in 1:nrow(meta.enrichscore.utuc)) {
#   tmp <- meta.enrichscore.utuc[i,]
#   #tmp <- wilcox.test(tmp[UTUC.annotation[which(UTUC.annotation$RNAcluster == "Muscle_Enriched"),"DetailID"]],tmp[UTUC.annotation[which(UTUC.annotation$RNAcluster == "Non-muscle_Enriched"),"DetailID"]])$p.value
#   tmp <- wilcox.test(tmp[UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C1"),"DetailID"]],tmp[UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C2"),"DetailID"]])$p.value
#   #tmp <- wilcox.test(tmp[UTUC.annotation[which(UTUC.annotation$iCluster == "C1"),"DetailID"]],tmp[UTUC.annotation[which(UTUC.annotation$iCluster == "C2"),"DetailID"]])$p.value
#   p <- append(p, tmp)
# }
# names(p) <- names(meta.sig)
# cna.meth.path <- "F:/Project/UTUC_BLCA/CNA/GISTIC/1750662"
# cna.meth <- read.table(file.path(cna.meth.path,"CNA.all_lesions.conf_99.txt"),sep = "\t",check.names = F,header = T,row.names = NULL,stringsAsFactors = F)
# del9p <- as.data.frame(t(cna.meth[which(cna.meth$`Unique Name` == "Deletion Peak 12"),10:44]))
# rownames(del9p) <- substr(rownames(del9p),start = 6,stop = 8);colnames(del9p) <- "del9p"

# ####################################################################################
# ############################## GSEA for phenotype C1 and C2 ########################
# sam_info <- data.frame("DMBcluster"=annCol.C1vsC2.rna$`Cluster DNA methylation no filter without normal`)
# rownames(sam_info) <- rownames(annCol.C1vsC2.rna)
# 
# in_gct <- combat.UTUC.FPKM.HUGO[,rownames(sam_info)]
# gct_file <- file.path(res.path,"UTUC_DMBclusterC1vsC2.gct_file.gct")
# cls_file <- file.path(res.path,"UTUC_DMBclusterC1vsC2.cls_file.cls")
# generateInputFileForGSEA(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "DMBcluster")
# 
# # GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute 
# #
# # R script to run GSEA Analysis of the DowningAll33OldETVPAX5 vs C2 example (cut and paste into R console)
# 
# GSEA.program.location <- "E:/Project/myOV/GSEA_Analysis_XiaofanLu/GSEA.1.0.R"   #  R source program (change pathname to the rigth location in local machine)
# source(GSEA.program.location, verbose=T, max.deparse.length=9999)
# 
# GSEA( # Input/Output Files :-------------------------------------------
#      input.ds =  "D:/UTUC_BLCA/Results/UTUC_DMBclusterC1vsC2.gct_file.gct",               # Input gene expression Affy dataset file in RES or GCT format
#      input.cls = "D:/UTUC_BLCA/Results/UTUC_DMBclusterC1vsC2.cls_file.cls",               # Input class vector (phenotype) file in CLS format
#      gs.db =     "E:/Project/myOV/GSEA_Analysis_XiaofanLu/msigdb.v6.0.symbols.gmt",           # Gene set database in GMT format
#      output.directory      = "D:/UTUC_BLCA/Results/UTUC.C1vsC2.GSEA/",            # Directory where to store output and results (default: "")
#      #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
#      doc.string            = "UTUC.C1vsC2",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
#      non.interactive.run   = T,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
#      reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
#      nperm                 = 1000,            # Number of random permutations (default: 1000)
#      weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
#      nom.p.val.threshold   = 0.05,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
#      fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
#      fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
#      topgs                 = 25,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
#      adjust.FDR.q.val      = T,               # Adjust the FDR q-vals (default: F)
#      gs.size.threshold.min = 10,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
#      gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
#      reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
#      preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
#      random.seed           = 1248103,             # Random number generator seed. (default: 123456)
#      perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
#      fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
#      replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
#      save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
#      OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
#      use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
# )
# #--------------------------------------------------------------------------------------------------------------------------------------------------
# 
# # Overlap and leading gene subset assignment analysis of the GSEA results
# 
# GSEA.Analyze.Sets(
#   directory           = "D:/UTUC_BLCA/Results/UTUC.C1vsC2.GSEA/",        # Directory where to store output and results (default: "")
#   topgs = 20,                                                           # number of top scoring gene sets used for analysis
#   height = 16,
#   width = 16
# )

####################
#### GISTIC cnv ####
gistic.path <- "F:/Project/UTUC_BLCA/CNA/GISTIC"
gistic <- read.table(file.path(gistic.path,"1725620_0.4/CNA_segment_forGISTIC2.0.all_thresholded.by_genes.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
for (i in 3:37) {
  write.table(table(gistic[,2],gistic[,i]),file.path(gistic.path,paste0(colnames(gistic)[i],"_all_thresholded.by_genes_table.txt")),sep = "\t",row.names = T,col.names = NA,quote = F)
}


#####################
#### TIDE for BLCA ##
# TIDE <- log2(FPKM.BLCA + 1)
# tmp <- apply(log2(FPKM.BLCA.Nromal + 1),1,mean,na.rm=T)
# TIDE <- sweep(TIDE,1, tmp)
# #TIDE <- sweep(TIDE,2, apply(TIDE,2,mean,na.rm=T))
# write.table(TIDE,file.path(res.path,"BLCA_TIDE_input2.txt"),sep = "\t",row.names = T,col.names = NA)

#############################
### UTUC cluster survival ###

# DMBcluster
tmp <- as.character(UTUC.annotation$`Cluster DNA methylation no filter without normal`); names(tmp) <- rownames(UTUC.annotation)
tmp <- tmp[which(tmp != "N/A")]
tmp <- factor(tmp,levels = c("C2","C1"))

kmplotfun2(group = tmp,Sinfo = UTUC.annotation,showHR = T,cutsurv = 50,fig.path = fig.path,outFigFile="KM2 plot DMBcluster.OS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(blue,red))
#0.03840005
kmplotfun.event1(group = tmp,Sinfo = UTUC.annotation,showHR = T,cutsurv = 50,fig.path = fig.path,outFigFile="KM2 plot DMBcluster.PFS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(blue,red))
#0.049
kmplotfun.event2(group = tmp,Sinfo = UTUC.annotation,showHR = T,cutsurv = 50,fig.path = fig.path,outFigFile="KM2 plot DMBcluster.PFS from distance.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(blue,red))
#0.036

# DMBcluster's NMI samples
tmp <- as.character(UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),"Cluster DNA methylation no filter without normal"]); names(tmp) <- rownames(UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),])
tmp <- tmp[which(tmp != "N/A")]
tmp <- factor(tmp,levels = c("C1","C2"))

kmplotfun2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot DMBcluster NMI only.OS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(red,blue))
#0.06
kmplotfun.event1(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot DMBcluster NMI only.PFS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(red,blue))
#0.07
kmplotfun.event2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot DMBcluster NMI only.PFS from distance.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(red,blue))
#0.09

# RNAcluster redo
tmp <- as.character(UTUC.annotation$RNAcluster); names(tmp) <- rownames(UTUC.annotation)
tmp <- tmp[which(tmp != "N/A")]
tmp <- ifelse(tmp == "Muscle_Enriched","C1","C2")
tmp <- factor(tmp,levels = c("C1","C2"))

kmplotfun2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot RNAcluster.OS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,nake))
#0.022
kmplotfun.event1(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot RNAcluster.PFS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,nake))

kmplotfun.event2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot RNAcluster.PFS from distance.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,nake))
#0.058

# MWBcluster #redo
tmp <- as.character(UTUC.annotation$MWBcluster); names(tmp) <- rownames(UTUC.annotation)
tmp <- tmp[which(tmp != "N/A")]
tmp <- factor(tmp,levels = c("C1","C2","C3"))

kmplotfun2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot MWBcluster.OS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c("#F09300",lightgreen,"#0068B5"))
#0.753
kmplotfun.event1(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot MWBcluster.PFS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c("#F09300",lightgreen,"#0068B5"))
#0.289
kmplotfun.event2(group = tmp,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot MWBcluster.PFS from distance.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c("#F09300",lightgreen,"#0068B5"))
#0.625
#######################################
### determine median follow-up time ###

# reverse the OS status and do survival analysis
tmp <- UTUC.annotation
tmp$OS <- abs(tmp$OS-1)
survfit(Surv(OS.time, OS) ~ 1, data = tmp)
# n  events  median 0.95LCL 0.95UCL 
# 40      28    1620    1082    2603

range(tmp$OS.time)/30.5
#######################################################
### Fisher test between del9p and invasiveness type ###
table(UTUC.annotation$Del9p21.3,UTUC.annotation$Type)
fisher.test(matrix(c(7,9,7,5),byrow = T,ncol = 2))

##############################################################
### Association between clinical information and surivival ###
tmp <- UTUC.annotation
#sex
tmp1 <- tmp$SEX; names(tmp1) <- rownames(tmp)
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot SEX.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot SEX.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot SEX.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
coxph(Surv(OS.time, OS) ~ UTUC.annotation$SEX, data = UTUC.annotation)
coxph(Surv(PFS.time, PFS) ~ UTUC.annotation$SEX, data = UTUC.annotation)

#invasive
tmp1 <- tmp$Type; names(tmp1) <- rownames(tmp)
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot INVASIVENESS.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot INVASIVENESS.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot INVASIVENESS.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

#age
tmp1 <- tmp$Age; names(tmp1) <- rownames(tmp)
tmp1 <- ifelse(tmp1 > 82,"H","L")
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot AGE.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot AGE.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot AGE.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

#stage
tmp1 <- tmp$`pathological stage`; names(tmp1) <- rownames(tmp)
tmp1 <- ifelse(tmp1 %in% c("pTa","pT1"),"pTa1","pT234"); names(tmp1) <- rownames(tmp)
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot STAGE.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot STAGE.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot STAGE.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

#localisation
tmp1 <- tmp$Localisation; names(tmp1) <- rownames(tmp)
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

#grade
tmp1 <- tmp$GRADE; names(tmp1) <- rownames(tmp)
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot LOCAL.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

#Del9p
tmp1 <- tmp$Del9p21.3; names(tmp1) <- rownames(tmp)
tmp1 <- tmp1[tmp1 != "N/A"]
kmplotfun2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot DEL9p.OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event1(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot DEL9p.PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))
kmplotfun.event2(group = tmp1,Sinfo = UTUC.annotation,showHR = F,cutsurv = 5,fig.path = fig.path,outFigFile=paste0("KM plot DEL9p.PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = c(cherry,lightgrey))

########################################################################
### Association between mutation and clinical information (survival) ###
tmp <- t(exome.mut)
tmp <- as.data.frame(tmp[,c("MLL2","FGFR3","MLL3","KDM6A","ZFP36L1","STAG2","CRIPAK","ARID1A","TP53","GANAB","ARID1B")])
# combine ARID1A and ARID1B
tmp$'ARID1A/B' <- tmp$ARID1A + tmp$ARID1B
#tmp <- tmp[,-c(8,11)]
tmp <- cbind.data.frame(tmp,t(exome.mut.all[c("ZFP36L2","ZFP36"),]))
tmp$ZFP36_family <- tmp$ZFP36L1 + tmp$ZFP36L2 + tmp$ZFP36
#tmp <- tmp[,-c(11,12)]
#tmp$SWI_SNF <- colSums(exome.mut.swisnf)
#tmp$SWI_SNF <- ifelse(tmp$SWI_SNF >= 1,1,0)
tmp$SWI_SNF <- SWISNF_mut_addACTL6B[rownames(tmp)]
tmp$SimpleID <- rownames(tmp)

tmp1 <- UTUC.annotation[,1:29]
tmp1$SimpleID <- rownames(tmp1)
#tmp2 <- merge(tmp1,tmp,by="SimpleID",all.x=T); tmp2[is.na(tmp2)] <- "N/A"; rownames(tmp2) <- tmp2$SimpleID; tmp2 <- tmp2[,-1]; tmp2 <- tmp2[rownames(UTUC.annotation),]
#write.table(tmp2,file.path(data.path,"UTUC.annotation2.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
tmp <- merge(tmp,tmp1,by="SimpleID",all.x=T)
tmp$Age <- ifelse(tmp$Age > 82,"H","L")
tmp$stage <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1"),"pTa1","pT234")

mutlab <- colnames(tmp)[2:17]
clilab <- colnames(tmp)[c(26,47,27,46,28,34)]

outTab <- as.data.frame(matrix(c(0),ncol = 6,nrow = 16,dimnames = list(mutlab,clilab)))
for (mut in mutlab) {
  cat(paste0(mut)," start!\n")
  for (cli in clilab) {
    p <- fisher.test(table(tmp[,cli],tmp[,mut]))
    outTab[mut,cli] = -log10(p$p.value)
    #if(p$p.value<0.1) {  
      print(table(tmp[,cli],tmp[,mut]))
      cat(paste0(cli," p=",p$p.value),"\n")
      cat("\n")
    #}
  }
}
# MLL2  start!
#   
#   0  1
# F  0  2
# M 22  6
# SEX p=0.064367816091954 

# 0  1
# renal pelvis 20  4
# ureter        2  4
# Localisation p=0.0293545534924845 

# FGFR3  start!
#   
#   0  1
# muscle-invasive     11  4
# Non-muscle invasive  4 11
# Type p=0.0268377292262022 
# 
# 
# 0  1
# pT234 11  4
# pTa1   4 11
# stage p=0.0268377292262022 
# 
# 
# 0  1
# High 15 11
# Low   0  4
# GRADE p=0.0996168582375479 

# ARID1A  start!
# 0  1
# muscle-invasive     10  5
# Non-muscle invasive 15  0
# Type p=0.0421455938697318 
# 
# 0  1
# pT234 10  5
# pTa1  15  0
# stage p=0.0421455938697318

# curated FGFR3
table(UTUC.annotation$FGFR3_mut,UTUC.annotation$Type)
# muscle-invasive Non-muscle invasive
# 0                13                   5
# 1                 5                  13
# N/A               2                   2
fisher.test(matrix(c(13,5,5,13),byrow = T,nrow = 2)) #0.018

for (mut in mutlab) {
  tmp2 <- tmp[,mut]; names(tmp2) <- tmp$SimpleID
  tmp2 <- ifelse(tmp2 == 1,paste0(mut,"_Mut"),paste0(mut,"_Wild"))
  if(mut == "ARID1A/B") {
    tmp2 <- gsub("/","_",tmp2)
    mut <- "ARID1A_B"
  }
  tmp2 <- factor(tmp2,levels = c(paste0(mut,"_Mut"),paste0(mut,"_Wild")))
  kmplotfun2(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile=paste0("KM2 plot ",mut,".OS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])
  kmplotfun.event1(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile=paste0("KM2 plot ",mut,".PFS.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])
  kmplotfun.event2(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile=paste0("KM2 plot ",mut,".PFS from distance.pdf"), xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])
}

tmp2 <- UTUC.annotation$FGFR3_mut; names(tmp2) <- rownames(UTUC.annotation); tmp2 <- tmp2[tmp2!="N/A"]
tmp2 <- ifelse(tmp2 == 1,"FGFR3_Mut","FGFR3_Wild")
tmp2 <- factor(tmp2,levels = c("FGFR3_Mut","FGFR3_Wild"))
kmplotfun2(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot curated FGFR3.OS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])
kmplotfun.event1(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot curated FGFR3.PFS.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])
kmplotfun.event2(group = tmp2,Sinfo = UTUC.annotation,showHR = F,cutsurv = 50,fig.path = fig.path,outFigFile="KM plot curated FGFR3.PFS from distance.pdf", xunit="Month",legend.title = "",width = 4.5,height = 4.5,col = jco[2:1])

# just calculate independency between genes and muscle-invasiveness and adjust pvalue
p.res <- c()
for (mut in mutlab) {
  cat(mut,"\n")
  for (cli in "Type") {
    print(table(tmp[,cli],tmp[,mut]))
    cat("\n")
    p <- fisher.test(table(tmp[,cli],tmp[,mut]))
    p.res <- c(p.res,p$p.value)
  }
}
names(p.res) <- mutlab
p.res <- c(p.res,FGFR3_curated = 0.018)

write.table(data.frame(nominal.p = round(as.numeric(p.res[1:10]),3),FDR=round(p.adjust(as.numeric(p.res[1:10]),"fdr"),3),row.names = mutlab[1:10],stringsAsFactors = F),
            file.path(res.path,"independent test between top 10 mutation and musle invasiveness with raw FGFR3.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

write.table(data.frame(nominal.p = round(as.numeric(p.res[c(1,3:10,17)]),3),FDR=round(p.adjust(as.numeric(p.res[c(1,3:10,17)]),"fdr"),3),row.names = c(mutlab[c(1,3:10)],"FGFR3"),stringsAsFactors = F),
            file.path(res.path,"independent test between top 10 mutation and musle invasiveness with curated FGFR3.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

#########################################################################
### association between clinical information and RNAcluster C1 and C2 ###
tmp <- UTUC.annotation[which(UTUC.annotation$RNAcluster %in% c("Non-muscle_Enriched","Muscle_Enriched")),]
tmp$Age <- ifelse(tmp$Age > 82,">82","<=82")
#tmp$stage <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1","pT2"),"pTa12","pT34")
tmp$stage <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1"),"pTa1","pT234")

tmp$OS.time <- tmp$OS.time/365
tmp$PFS.time <- tmp$PFS.time/365
#sex
fisher.test(table(tmp$SEX,tmp$RNAcluster)) #1
#age
fisher.test(table(tmp$Age,tmp$RNAcluster)) #0.1968
#type
fisher.test(table(tmp$Type,tmp$RNAcluster)) #0.1273
#pathological
fisher.test(table(tmp$`pathological stage`,tmp$RNAcluster)) 
fisher.test(table(tmp$stage,tmp$RNAcluster)) #0.1273

#grade
fisher.test(table(tmp$GRADE,tmp$RNAcluster)) #
#localisation
fisher.test(table(tmp$Localisation,tmp$RNAcluster)) #0.3618

#########################################################################
### association between clinical information and DMBcluster C1 and C2 ###
tmp <- UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` %in% c("C1","C2")),]
tmp$DMBcluster <- tmp$`Cluster DNA methylation no filter without normal`
tmp$Age <- ifelse(tmp$Age > 82,">82","<=82")
#tmp$stage <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1","pT2"),"pTa12","pT34")
tmp$stage <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1"),"pTa1","pT234")

tmp$OS.time <- tmp$OS.time/365
tmp$PFS.time <- tmp$PFS.time/365
#sex
fisher.test(table(tmp$SEX,tmp$DMBcluster))
#age
fisher.test(table(tmp$Age,tmp$DMBcluster))
#type
fisher.test(table(tmp$Type,tmp$DMBcluster)) #0.000945
#pathological
fisher.test(table(tmp$stage,tmp$DMBcluster)) #0.000945
fisher.test(table(tmp$`pathological stage`,tmp$DMBcluster)) #0.005439

#grade
fisher.test(table(tmp$GRADE,tmp$DMBcluster)) #0.1061
#localisation
fisher.test(table(tmp$Localisation,tmp$DMBcluster)) #0.4506

#unvariateCox
outTab <- NULL
clilab <- c("SEX","Type","GRADE","stage","Localisation","Age","DMBcluster")
for (i in clilab) {
  cox <- coxph(Surv(OS.time, OS) ~ tmp[,i], data = tmp)
  coxSummary = summary(cox)
  outTab=rbind.data.frame(outTab,cbind.data.frame(Cha=i,HR=round(coxSummary$coefficients[,"exp(coef)"],2),
                                                  lower.95=coxSummary$conf.int[3],upper.95=coxSummary$conf.int[4],
                                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"]),stringsAsFactors = F)
}

#multivariateCox
rt <- as.data.frame(na.omit(tmp[,c(clilab,"PFS2","PFS.time")]))
cox <- coxph(Surv(PFS.time, PFS2) ~ ., data = rt)


############################################
### association mutation with DMBcluster ###
UTUC.annotation$SWI_SNF_noACTL6B <- SWISNF_mut_noACTL6B[rownames(UTUC.annotation)]
UTUC.annotation[is.na(UTUC.annotation$SWI_SNF_noACTL6B),"SWI_SNF_noACTL6B"] <- "N/A"
tmp <- UTUC.annotation[which(UTUC.annotation$`DNA Methylation data` == "yes"),c(8,19,23,30:46)]
for (i in 2:20) {
  tmp1 <- tmp[,c(1,i)]
  if(is.element("N/A",tmp1[,2])) {
    tmp1 <- tmp1[which(tmp1[,2] != "N/A"),]
    if(length(unique(tmp1[,2])) == 1) {next()}
  }
  cat(paste0(colnames(tmp1)[2]," starts! p=",fisher.test(table(tmp1[,1],tmp1[,2]))$p.value,"\n"))
  print(table(tmp1[,1],tmp1[,2]))
}

# FGFR3_mut starts! p=0.0031887282165369
# 
# 0  1
# C1 14  5
# C2  2 10

# FGFR3 starts! p=0.00442410373760489
# 
# 0  1
# C1 13  5
# C2  1  9

# ARID1A/B starts! p=0.030171277997365
# 
# 0  1
# C1 11  7
# C2 10  0

# SWI_SNF starts! p=0.015870780630506
# 
# 0  1
# C1  7 11
# C2  9  1

# SWI_SNF_noACTL6B starts! p=0.0098162706858359
# 
# 0  1
# C1  9  9
# C2 10  0

p.adjust(c(0.4,0.004,0.207,0.626,1,1,1,0.128,0.626,1,0.0098),method = "fdr")
############################################
### association mutation with RNAcluster ###
tmp <- annCol.rna[,c(21,24:40)]

for (i in 2:18) {
  tmp1 <- tmp[,c(1,i)]
  if(is.element("N/A",tmp1[,2])) {
    tmp1 <- tmp1[which(tmp1[,2] != "N/A"),]
    if(length(unique(tmp1[,2])) == 1) {next()}
  }
  cat(paste0(colnames(tmp1)[2]," starts! p=",fisher.test(table(tmp1[,1],tmp1[,2]))$p.value,"\n"))
  print(table(tmp1[,1],tmp1[,2]))
}

# FGFR3 starts! p=0.0859
# 
# 0 1
# Muscle_Enriched     4 5
# Non-muscle_Enriched 5 0

# TP53 starts! p=0.0949
# 
# 0 1
# Muscle_Enriched     8 1
# Non-muscle_Enriched 2 3

############################################
### association mutation with MWBcluster ###
tmp <- UTUC.annotation[which(UTUC.annotation$`Exome-data` == "yes"),c(21,19,23,30:45)]
#tmp$MWBcluster <- ifelse(tmp$MWBcluster == "C3","C3","Others") # no significance
for (i in 2:19) {
  tmp1 <- tmp[,c(1,i)]
  if(is.element("N/A",tmp1[,2])) {
    tmp1 <- tmp1[which(tmp1[,2] != "N/A"),]
    if(length(unique(tmp1[,2])) == 1) {next()}
  }
  cat(paste0(colnames(tmp1)[2]," starts! p=",fisher.test(table(tmp1[,1],tmp1[,2]))$p.value,"\n"))
  print(table(tmp1[,1],tmp1[,2]))
} # no significance
###################################
### supplement dated 02/23/2019 ###

#box plot for immune/stromal score between invasiveness
tmp1 <- annCol.rna[which(annCol.rna$Type == "muscle-invasive"),"ImmuneScore"]
tmp2 <- annCol.rna[which(annCol.rna$Type == "Non-muscle invasive"),"ImmuneScore"]
p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.03
df <- data.frame("ImmuneScore"=c(tmp1,tmp2),"Invasiveness"=rep(c("MI","NMI"),c(length(tmp1),length(tmp2))))
p1 <- ggplot(df,aes(x=Invasiveness,y=ImmuneScore,fill=Invasiveness)) + geom_boxplot() + scale_fill_manual(values = c(soil,grey)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")

tmp1 <- annCol.rna[which(annCol.rna$Type == "muscle-invasive"),"StromalScore"]
tmp2 <- annCol.rna[which(annCol.rna$Type == "Non-muscle invasive"),"StromalScore"]
p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.03
df <- data.frame("StromalScore"=c(tmp1,tmp2),"Invasiveness"=rep(c("MI","NMI"),c(length(tmp1),length(tmp2))))
p2 <- ggplot(df,aes(x=Invasiveness,y=StromalScore,fill=Invasiveness)) + geom_boxplot() + scale_fill_manual(values = c(soil,grey)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")
p <- plot_grid(p1, p2, labels = c("a", "b"))
ggsave(file.path(fig.path,"boxplot of Immune and stromal score between Invasiveness MI and NMI .pdf"),width = 5,height = 4)

#association between immune/stromal score and pathological stage (see if advanced stage increase immunity)

# use original est.score.combat
tmp <- as.data.frame(t(est.score.combat)); tmp <- tmp[,setdiff(colnames(tmp),c("ESTIMATEScore","TumorPurity"))]; tmp <- as.data.frame(t(tmp),stringsAsFactors=F)
tmp <- as.data.frame(sapply(tmp, as.numeric)); rownames(tmp) <- setdiff(rownames(est.score.combat),c("ESTIMATEScore","TumorPurity"))
colnames(tmp) <- sapply(strsplit(colnames(tmp),"_"),"[",2)
colnames(tmp) <- paste0("0",colnames(tmp))
colnames(tmp) <- substr(colnames(tmp),nchar(colnames(tmp))-2,nchar(colnames(tmp)))
tmp <- as.data.frame(t(tmp[,rownames(annCol.rna)]))
tmp <- scale(tmp)
tmp <- cbind.data.frame(tmp,annCol.rna[,c("Pathological_Stage","Type")])

# tmp1 <- tmp[which(tmp$Type == "muscle-invasive"),"ImmuneScore"]
# tmp2 <- tmp[which(tmp$Type == "Non-muscle invasive"),"ImmuneScore"]
# p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.03
# df <- data.frame("ImmuneScore"=c(tmp1,tmp2),"Invasiveness"=rep(c("MI","NMI"),c(length(tmp1),length(tmp2))))
# p1 <- ggplot(df,aes(x=Invasiveness,y=ImmuneScore,fill=Invasiveness)) + geom_boxplot() + scale_fill_manual(values = c(soil,grey)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")
# 
# tmp1 <- tmp[which(tmp$Type == "muscle-invasive"),"StromalScore"]
# tmp2 <- tmp[which(tmp$Type == "Non-muscle invasive"),"StromalScore"]
# p <- round(wilcox.test(tmp1,tmp2,alternative = "greater")$p.value,2) # 0.03
# df <- data.frame("StromalScore"=c(tmp1,tmp2),"Invasiveness"=rep(c("MI","NMI"),c(length(tmp1),length(tmp2))))
# p2 <- ggplot(df,aes(x=Invasiveness,y=StromalScore,fill=Invasiveness)) + geom_boxplot() + scale_fill_manual(values = c(soil,grey)) + ggtitle(paste0("pValue = ",p)) + theme(legend.position="none")
# p <- plot_grid(p1, p2, labels = c("a", "b"))
# ggsave(file.path(fig.path,"boxplot of Immune and stromal score between Invasiveness MI and NMI .pdf"),width = 5,height = 4)
tmp$Stage <- ifelse(tmp$Pathological_Stage %in% c("pTa","pT1"),"pTa1","pT234")
tmp$Stage <- factor(tmp$Stage,levels = c("pTa1","pT234"))
tmp$stage <- sapply(tmp$Stage,function(x) {switch(as.character(x),
                                                  "pTa1"=1,
                                                  "pT234"=2)})
cor.test(tmp$ImmuneScore,tmp$stage,method = "spearman")
cor.test(tmp$StromalScore,tmp$stage,method = "spearman")

# use twoAxis.R plot

# use standard.fun tailed score

# tmp <- annCol.rna[,c("ImmuneScore","StromalScore","Pathological_Stage","Type")]
# tmp$Stage <- ifelse(tmp$Pathological_Stage %in% c("pTa","pT1"),"pTa1","pT234")
#                    
# tmp$Stage <- factor(tmp$Stage,levels = c("pTa1","pT234"))
# tmp$stage <- sapply(tmp$Stage,function(x) {switch(as.character(x),
#                                                                "pTa1"=1,
#                                                                "pT234"=2)})
# cor.test(tmp$ImmuneScore,tmp$stage,alternative = "greater",method = "spearman")
# cor.test(tmp$StromalScore,tmp$stage,alternative = "greater",method = "spearman")

# differentially expression analysis between invasiveness and perform GSEA
tmp <- UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes"),"Type"]; names(tmp) <- UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"]
complist <- createList.invasiveness(group=tmp)
group <- data.frame("group"=as.character(tmp),"batch"=ifelse(grepl("rna",names(tmp)),"UTUCrna","UTUCr"),row.names = names(tmp))
rowids <- intersect( Mids, names(PASSFlag.mRNA.UTUC) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[rowids, rownames(group)], tailrows, Groupinfo=group, features=Mids, featType="mRNA_rmBatch", complist,PASSFlag=PASSFlag.mRNA.UTUC, overwt=TRUE)

tmp <- read.table(file.path(res.path,"mRNA_rmBatch_edgeR_test_result.muscle-invasive_vs_Non-muscle invasive.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$logFC
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
#MSigDB=read.gmt("F:/Project/gsea.xlu/GeneSetDataBases/msigdb.v6.0.symbols.gmt")
GSEA.edgeR.invasiveness <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F)
res <- data.frame(GSEA.edgeR.invasiveness)
write.table(as.data.frame(GSEA.edgeR.invasiveness),file.path(res.path,"GSEA results for muscle-invasive vs non-muscle invasive DEGs of edgeR by clusterprofiler.txt"),row.names = T,col.names = NA,sep = "\t")

#################################################
### supervised clustering by immune signature ###
tmp <- read.table(file = file.path(res.path,paste("UTUC","_FPKM_all_features_HUGO_symbols.txt",sep = "")),row.names=1, check.names = F,stringsAsFactors = F,header = T, sep="\t")
immune.signature <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig <- list()
for (i in cell.type) {
  immune.sig[[i]] <- toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"])
}
immune.sig[["Angiogenesis"]] <- c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")
immune.enrichscore<- gsva(as.matrix(tmp),immune.sig,method="gsva")
immune.enrichscore.combat <- ComBat(dat=as.matrix(immune.enrichscore), batch=Sinfo[c(427:446),"batch"], mod=modcombat)

indata <- immune.enrichscore.combat
index <- sapply(strsplit(colnames(indata),"_"), "[",2)
index <- paste0("0",index); index <- substr(index,start = nchar(index)-2,stop = nchar(index))
colnames(indata) <- index

hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward.D")
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward.D")

plotdata <- standarize.fun(indata,halfwidth = 1)
pdf(file = file.path(fig.path,"supervised clustering by immune signature.pdf"),height = 11.5,width = 10)
aheatmap(as.matrix(plotdata), 
         Rowv=dendsort(as.dendrogram(hcg)), 
         Colv=dendsort(as.dendrogram(hcs)), 
         annCol=annCol.rna[colnames(plotdata),c(1:12,14,16,19,28,29,32)],
         annColors = annColors.rna, 
         color=blueyellow(64), 
         revC=TRUE, fontsize=7,cexAnn = 1,cexCol = 2.5,labRow = NA)
invisible(dev.off())

############################################################################################
### test if any MCPcounter cells differentially expression in TCGA ZFP36 family mutation ###

for (i in rownames(MCPscore.BLCA)) {
  tmp1 <- log2(as.numeric(MCPscore.BLCA[i,blca.ZFP36.mutsam])+1)
  tmp2 <- log2(as.numeric(MCPscore.BLCA[i,blca.ZFP36.norsam])+1)
  p1 <- t.test(tmp1,tmp2)$p.value
  
  tmp1 <- as.numeric(MCPscore.BLCA[i,blca.ZFP36.mutsam])
  tmp2 <- as.numeric(MCPscore.BLCA[i,blca.ZFP36.norsam])
  p2 <- wilcox.test(tmp1,tmp2)$p.value
  
  cat(paste0(i," t.test.p=",p1,"; wilcox.p=",p2,"\n"))
}
# T cells t.test.p=0.180047359285845; wilcox.p=0.113356110284598
# CD8 T cells t.test.p=0.514374773378745; wilcox.p=0.288326727053068
# Cytotoxic lymphocytes t.test.p=0.84099608908021; wilcox.p=0.640213629984086
# NK cells t.test.p=0.480280875871112; wilcox.p=0.550988283011256
# B lineage t.test.p=0.735672973487219; wilcox.p=0.713476028917447
# Monocytic lineage t.test.p=0.820666719590444; wilcox.p=0.976017276238565
# Myeloid dendritic cells t.test.p=0.99057566479307; wilcox.p=0.890126475737925
# Neutrophils t.test.p=0.284450357772327; wilcox.p=0.24915978111778
# Endothelial cells t.test.p=0.290379661415647; wilcox.p=0.700174385901698
# Fibroblasts t.test.p=0.481830953816042; wilcox.p=0.464025677182883

immune.signature <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig <- list()
for (i in cell.type) {
  immune.sig[[i]] <- toupper(immune.signature[which(immune.signature$CellType == i),"Symbol"])
}
immune.sig[["Angiogenesis"]] <- c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")
immune.enrichscore.blca <- gsva(as.matrix(FPKM.BLCA),immune.sig,method="gsva")

for (i in rownames(immune.enrichscore.blca)) {
  tmp1 <- log2(as.numeric(immune.enrichscore.blca[i,blca.ZFP36.mutsam])+1)
  tmp2 <- log2(as.numeric(immune.enrichscore.blca[i,blca.ZFP36.norsam])+1)
  p1 <- t.test(tmp1,tmp2)$p.value
  
  tmp1 <- as.numeric(immune.enrichscore.blca[i,blca.ZFP36.mutsam])
  tmp2 <- as.numeric(immune.enrichscore.blca[i,blca.ZFP36.norsam])
  p2 <- wilcox.test(tmp1,tmp2)$p.value
  
  cat(paste0(i," t.test.p=",p1,"; wilcox.p=",p2,"\n"))
}
# B cells t.test.p=0.674126785466186; wilcox.p=0.582025713805639
# T cells t.test.p=0.916357364904426; wilcox.p=0.927569693388747
# T helper cells t.test.p=0.169130089538286; wilcox.p=0.181904614535725
# Tcm t.test.p=0.324313890458089; wilcox.p=0.382141746229681
# Tem t.test.p=0.742857860113229; wilcox.p=0.790584885290089
# Th1 cells t.test.p=0.858932561491971; wilcox.p=0.904852446437714
# Th2 cells t.test.p=0.919764084421855; wilcox.p=0.821615227170397
# TFH t.test.p=0.971355321169718; wilcox.p=0.910525198663887
# Th17 cells t.test.p=0.925276386459663; wilcox.p=0.930982838257677
# TReg t.test.p=0.463153829264126; wilcox.p=0.457499320335952
# CD8 T cells t.test.p=0.178978136527293; wilcox.p=0.222316985522686
# Tgd t.test.p=0.118312468564883; wilcox.p=0.106824279131796
# Cytotoxic cells t.test.p=0.900390033001902; wilcox.p=0.883342664377539
# NK cells t.test.p=0.965013778623453; wilcox.p=0.986865225728077
# NK CD56dim cells t.test.p=0.678775416253008; wilcox.p=0.587929463834527
# NK CD56bright cells t.test.p=0.0829236068521403; wilcox.p=0.0620391264675752
# DC t.test.p=0.265898676008766; wilcox.p=0.218543231423487
# iDC t.test.p=0.890307723403679; wilcox.p=0.70176596372275
# aDC t.test.p=0.679310005501274; wilcox.p=0.677510127182292
# Eosinophils t.test.p=0.331051132834577; wilcox.p=0.427317423133529
# Macrophages t.test.p=0.696014114877893; wilcox.p=0.754438620949819
# Mast cells t.test.p=0.875617978063437; wilcox.p=0.876567164513363
# Neutrophils t.test.p=0.711071602819494; wilcox.p=0.695407773419942
# SW480 cancer cells t.test.p=0.631391554760953; wilcox.p=0.714544783219586
# Normal mucosa t.test.p=0.30643257034447; wilcox.p=0.379807254441795
# Blood vessels t.test.p=0.430324648085804; wilcox.p=0.934966429495586
# Lymph vessels t.test.p=0.368020523689352; wilcox.p=0.43654337107995
# Angiogenesis t.test.p=0.239880843871616; wilcox.p=0.190970223555811

#######################################################################################
### test if any MCPcounter cells differentially expression in TCGA ZFP36L1 mutation ###

for (i in rownames(MCPscore.BLCA)) {
  tmp1 <- log2(as.numeric(MCPscore.BLCA[i,blca.ZFP36L1.mutsam])+1)
  tmp2 <- log2(as.numeric(MCPscore.BLCA[i,blca.ZFP36L1.norsam])+1)
  p1 <- t.test(tmp1,tmp2)$p.value
  
  tmp1 <- as.numeric(MCPscore.BLCA[i,blca.ZFP36L1.mutsam])
  tmp2 <- as.numeric(MCPscore.BLCA[i,blca.ZFP36L1.norsam])
  p2 <- wilcox.test(tmp1,tmp2)$p.value
  
  cat(paste0(i," t.test.p=",p1,"; wilcox.p=",p2,"\n"))
}

# T cells t.test.p=0.413643810935741; wilcox.p=0.594904143127274
# CD8 T cells t.test.p=0.124339718638968; wilcox.p=0.0739314338364606
# Cytotoxic lymphocytes t.test.p=0.0776933079685501; wilcox.p=0.511038279734592
# NK cells t.test.p=0.345391024653083; wilcox.p=0.44154330265968
# B lineage t.test.p=0.724036034843636; wilcox.p=0.567843744551942
# Monocytic lineage t.test.p=0.788936610258508; wilcox.p=0.967103937575407
# Myeloid dendritic cells t.test.p=0.260129080735773; wilcox.p=0.434319608292085
# Neutrophils t.test.p=0.486799967878724; wilcox.p=0.387904286538476
# Endothelial cells t.test.p=0.619703441535744; wilcox.p=0.983898462258915
# Fibroblasts t.test.p=0.96716361799388; wilcox.p=0.96710395513117

# create high methylation pct file for Gaby
tmp1 <- read.table(file.path(res.path,"high_methylation_pct_of_genes_noXY_cpg_promoter_UTUC.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp1 <- tmp1[which(tmp1$HighMethPctInTumor > 0.1 & tmp1$HighMethPctInNormal == 0),]

tmp2 <- read.table(file.path("F:/Project/PanCan/BLCA/Results/high_methylation_pct_of_genes_BLCA.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- cbind.data.frame(tmp1,tmp2[rownames(tmp1),])
colnames(tmp) <- c("HighMethPctInTumor.UTUC","HighMethPctInNormal.UTUC","HighMethPctInTumor.BLCA","HighMethPctInNormal.BLCA")
write.table(tmp,file.path(res.path,"905 highmethpct in UTUC mapping TCGA BLCA.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
highmethpct.genes <- rownames(tmp)

# differentially methylation analysis on 905 highmethpct genes regarding invasiveness and DMBcluster C1 and C2
tmp <- read.table(file.path(data.path,"UTUC_medianAvg_Gmeth_noXY.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- tmp[highmethpct.genes,]
indata <- tmp

tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),]),rownames(UTUC.annotation[which(UTUC.annotation$`DNA Methylation data` == "yes"),]))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "muscle-invasive"),]),rownames(UTUC.annotation[which(UTUC.annotation$`DNA Methylation data` == "yes"),]))

am <- bm <- pvalues <- padj <- fold <- diff <- log2fold <- c()
Ngene <- nrow(indata)
for (m in 1:Ngene) {
  x <- as.numeric(indata[m, tmp1])
  y <- as.numeric(indata[m, tmp2])
  res <- wilcox.test(x=x, y=y, alternative="two.sided", paired=FALSE)
  am <- c(am, mean(x,na.rm=TRUE))
  bm <- c(bm, mean(y,na.rm=TRUE))      
  pvalues <- c(pvalues, res$p.value)
  display.progress(index=m, totalN=Ngene)
}
fold <- bm/am
diff <- bm-am
log2fold <- log2(fold)
padj = p.adjust( pvalues, method="BH" )

out <- data.frame(am, bm, fold, log2fold, diff, pvalues, padj,row.names = rownames(indata))
colnames(out) <- c(paste("meanbeta.", "NMI", sep=""), paste("meanbeta.", "MI", sep=""), "folcChange", "log2FoldChange", "diff", "pval", "padj")
write.table(out,file.path(res.path,"Differentially methylation results between MI and NMI by 905 highmethpct genes.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

#
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C2"),]),rownames(UTUC.annotation[which(UTUC.annotation$`DNA Methylation data` == "yes"),]))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C1"),]),rownames(UTUC.annotation[which(UTUC.annotation$`DNA Methylation data` == "yes"),]))

am <- bm <- pvalues <- padj <- fold <- diff <- log2fold <- c()
Ngene <- nrow(indata)
for (m in 1:Ngene) {
  x <- as.numeric(indata[m, tmp1]) 
  y <- as.numeric(indata[m, tmp2]) 
  res <- wilcox.test(x=x, y=y, alternative="two.sided", paired=FALSE)
  am <- c(am, mean(x,na.rm=TRUE)) 
  bm <- c(bm, mean(y,na.rm=TRUE)) 
  pvalues <- c(pvalues, res$p.value)
  display.progress(index=m, totalN=Ngene)
}
fold <- bm/am
diff <- bm-am
log2fold <- log2(fold)
padj = p.adjust( pvalues, method="BH" )

out <- data.frame(am, bm, fold, log2fold, diff, pvalues, padj,row.names = rownames(indata))
colnames(out) <- c(paste("meanbeta.", "DMBcluster_C2", sep=""), paste("meanbeta.", "DMBcluster_C1", sep=""), "folcChange", "log2FoldChange", "diff", "pval", "padj")
write.table(out,file.path(res.path,"Differentially methylation results between DMBcluster C1 and C2 by 905 highmethpct genes.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# create ZFP36 family mutation and expression file for Gaby
tmp <- FPKM.BLCA[,blca.ZFP36.mut$samples]
tmp1 <- read.table(file.path(res.path,"mRNA_ZFP36_family_BLCA_edgeR_test_result.Mutation_vs_Normal.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
tmp <- tmp[rownames(tmp1),]
tmp <- cbind.data.frame(tmp,tmp1)
tmp <- rbind.data.frame(tmp,as.character(blca.ZFP36.mut$group))
tmp <- tmp[c(nrow(tmp),1:(nrow(tmp)-1)),]
rownames(tmp)[1] <- "ZNF36_Family_Mutation"
write.table(tmp,file.path(res.path,"FPKM expression matrix with ZFP36 family mutation DE results and mutation status annotation.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

#######################################################################################
### test difference of pcaMeTil score in over 20% mutations and muscle-invasiveness ###
pca.MeTIL.utuc <- read.table(file.path(res.path,"pca.MeTIL.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

# FGFR3
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 7.822e-05
df <- data.frame("MeTIL"=c(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]),"FGFR3_Mutation"=rep(c("Mutated","Wild"),c(14,14)))
p2 <- ggplot(df,aes(x=FGFR3_Mutation,y=MeTIL,fill=FGFR3_Mutation)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(purple,grey)) + scale_color_manual(values = c(purple,grey)) + 
  ggtitle("pValue = 7.822e-5") + 
  geom_jitter(aes(color = FGFR3_Mutation),alpha=1) + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"pcaMeTIL between UTUC FGFR3 Mutation .pdf"),width = 4,height = 4.5)

# FGFR3_curated
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 9.118e-05

# KDM6A
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$KDM6A == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$KDM6A == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.6405

# MLL2
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL2 == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL2 == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.9009

# MLL3
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL3 == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL3 == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.02737

# ZFP36L1
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.4474

# ZFP36_family
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36_family == "1"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36_family == "0"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.6405

# muscle-invasiveness
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "muscle-invasive"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.08933
df <- data.frame("MeTIL"=c(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]),"Type"=rep(c("MI","NMI"),c(17,18)))
p2 <- ggplot(df,aes(x=Type,y=MeTIL,fill=Type)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(soil,grey)) + scale_color_manual(values = c(soil,grey)) + 
  ggtitle("pValue = 0.08933") + 
  geom_jitter(aes(color = Type),alpha=1) + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"pcaMeTIL between UTUC MI .pdf"),width = 4,height = 4.5)

# swi/snf
tmp1 <- intersect(names(SWISNF_mut_noACTL6B[SWISNF_mut_noACTL6B == 1]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(names(SWISNF_mut_noACTL6B[SWISNF_mut_noACTL6B == 0]),rownames(pca.MeTIL.utuc))
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 0.9615

# test METIL in MI and NMI regarding FGFR mutation
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "1" & UTUC.annotation$Type == "muscle-invasive"),]),rownames(pca.MeTIL.utuc))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "1" & UTUC.annotation$Type == "Non-muscle invasive"),]),rownames(pca.MeTIL.utuc))
tmp3 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "0" & UTUC.annotation$Type == "muscle-invasive"),]),rownames(pca.MeTIL.utuc))
tmp4 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "0" & UTUC.annotation$Type == "Non-muscle invasive"),]),rownames(pca.MeTIL.utuc))

wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp2,1]) # 
wilcox.test(pca.MeTIL.utuc[tmp3,1],pca.MeTIL.utuc[tmp4,1]) # 
wilcox.test(pca.MeTIL.utuc[tmp1,1],pca.MeTIL.utuc[tmp3,1]) # 0.007882
wilcox.test(pca.MeTIL.utuc[tmp2,1],pca.MeTIL.utuc[tmp4,1]) # 0.02398

#############################################################################################
### test difference of immune/stromal score in over 20% mutations and muscle-invasiveness ###

# FGFR3
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3 == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.01593
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.003216
df <- data.frame("ImmuneScore"=c(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]),"FGFR3_Mutation"=rep(c("Mutated","Wild"),c(5,9)))
p2 <- ggplot(df,aes(x=FGFR3_Mutation,y=ImmuneScore,fill=FGFR3_Mutation)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(purple,grey)) + scale_color_manual(values = c(purple,grey)) + 
  ggtitle("pValue = 0.016") + 
  #geom_jitter(aes(color = FGFR3_Mutation),alpha=1) + 
  #stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"Immunescore between UTUC FGFR3 Mutation .pdf"),width = 4.5,height = 4)

df <- data.frame("StromalScore"=c(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]),"FGFR3_Mutation"=rep(c("Mutated","Wild"),c(5,9)))
p2 <- ggplot(df,aes(x=FGFR3_Mutation,y=StromalScore,fill=FGFR3_Mutation)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(purple,grey)) + scale_color_manual(values = c(purple,grey)) + 
  ggtitle("pValue = 0.003") + 
  #geom_jitter(aes(color = FGFR3_Mutation),alpha=1) + 
  #stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"StromalScore between UTUC FGFR3 Mutation .pdf"),width = 4.5,height = 4)

# FGFR3_curated
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.05277
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.01632

# KDM6A
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$KDM6A == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$KDM6A == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.7833
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 1

# MLL2
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL2 == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL2 == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.4347
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.2273

# MLL3
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL3 == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$MLL3 == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.4092
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.4092

# ZFP36L1
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.521
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.6466

# ZFP36_family
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36_family == "1"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$ZFP36_family == "0"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.1593
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.7544

# swisnf
tmp1 <- intersect(names(SWISNF_mut_noACTL6B[SWISNF_mut_noACTL6B == 1]),rownames(annCol.rna))
tmp2 <- intersect(names(SWISNF_mut_noACTL6B[SWISNF_mut_noACTL6B == 0]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.7544
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.7544

# muscle-invasiveness
tmp1 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "muscle-invasive"),]),rownames(annCol.rna))
tmp2 <- intersect(rownames(UTUC.annotation[which(UTUC.annotation$Type == "Non-muscle invasive"),]),rownames(annCol.rna))
wilcox.test(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]) # 0.06567
wilcox.test(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]) # 0.05382

df <- data.frame("ImmuneScore"=c(annCol.rna[tmp1,"ImmuneScore"],annCol.rna[tmp2,"ImmuneScore"]),"Type"=rep(c("MI","NMI"),c(15,5)))
p2 <- ggplot(df,aes(x=Type,y=ImmuneScore,fill=Type)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(soil,grey)) + scale_color_manual(values = c(soil,grey)) + 
  ggtitle("pValue = 0.066") + 
  #geom_jitter(aes(color = FGFR3_Mutation),alpha=1) + 
  #stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"ImmuneScore between UTUC muscle-invasiveness .pdf"),width = 4.5,height = 4)

df <- data.frame("StromalScore"=c(annCol.rna[tmp1,"StromalScore"],annCol.rna[tmp2,"StromalScore"]),"Type"=rep(c("MI","NMI"),c(15,5)))
p2 <- ggplot(df,aes(x=Type,y=StromalScore,fill=Type)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = c(soil,grey)) + scale_color_manual(values = c(soil,grey)) + 
  ggtitle("pValue = 0.054") + 
  #geom_jitter(aes(color = FGFR3_Mutation),alpha=1) + 
  #stat_summary(fun.y=mean, geom="point", shape=18, size=4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
ggsave(file.path(fig.path,"StromalScore between UTUC muscle-invasiveness .pdf"),width = 4.5,height = 4)

##################################################################################################
### make heatmap for differentia expression genes' DAVID pathways between RNAcluster C1 and C2 ###
david <- read.table(file.path(data.path,"DAVID pathways related to genes up and down in RNAcluster C1 and C2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
david.list <- list()
for (i in 1:nrow(david)) {
  label <- sapply(strsplit(david[i,"Term"],"~"),"[",1)
  david.list[[label]] <- intersect(rownames(combat.UTUC.FPKM.HUGO),unlist(strsplit(david[i,"Genes"],", ")))
}

david.enrichscore <- gsva(as.matrix(combat.UTUC.FPKM.HUGO),david.list,method="ssgsea")
colnames(david.enrichscore) <- paste0(0,sapply(strsplit(colnames(david.enrichscore),"_"),"[",2))
colnames(david.enrichscore) <- substr(colnames(david.enrichscore),start = nchar(colnames(david.enrichscore))-2,stop = nchar(colnames(david.enrichscore)))

p <- vector()
for (i in 1:nrow(david.enrichscore)) {
  tmp <- round(wilcox.test(as.numeric(david.enrichscore[i,rownames(UTUC.annotation[which(UTUC.annotation$RNAcluster == "Muscle_Enriched"),])]),
                           as.numeric(david.enrichscore[i,rownames(UTUC.annotation[which(UTUC.annotation$RNAcluster == "Non-muscle_Enriched"),])]))$p.value,4)
  p <- append(p, tmp)
}
rownames(david.enrichscore) <- paste0(rownames(david.enrichscore),"_p_",p)

plotscore <- as.data.frame(na.omit(standarize.fun(indata = david.enrichscore,halfwidth = 1)))
hcg <- hclust(distanceMatrix(as.matrix(t(plotscore)), "euclidean"), "ward.D")
mycol <- colorpanel(64,low=cyan,mid = "black",high=peach)
pdf(file.path(fig.path, "UTUC_heatmap_mRNA_DAVID_pathways_for_DEGs_of_RNACluster.pdf"), width=7,height = 4)
hv <- aheatmap(as.matrix(plotscore), 
               Rowv=NA, Colv=cluster.UTUC.20RNAseq.ans$dendro, 
               annCol=annCol.rna[colnames(plotdata),21,drop = F],
               annColors = annColors.rna,
               color=mycol, 
               revC=TRUE, fontsize=9, cexCol = 1, cexRow = 1,cexAnn = 1)
invisible(dev.off())


###################################################################
### use NTP to test if UTUC samples belongs to luminal or basal ###
tmp <- Sinfo$TCGA_Subtype; names(tmp) <- rownames(Sinfo)
tmp <- tmp[tmp %in% c("Luminal_papillary","Luminal_infiltrated","Basal_squamous","Luminal")]; tmp1 <- names(tmp)
tmp <- ifelse(grepl("Luminal",tmp),"Luminal","Basal"); names(tmp) <- tmp1
complist <- createList.luminal.basal(group=tmp)
tmp1 <- data.frame(tmp,sample=names(tmp)); colnames(tmp1)[1] <- "group"
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, names(tmp)], tailrows, Groupinfo=tmp1, features=Mids, featType="mRNA_BLCA", complist,PASSFlag=PASSFlag.mRNA, overwt=F)

tmp <- read.table(file.path(res.path,"mRNA_BLCA_edgeR_test_result.Luminal_vs_Basal.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
luminal.marker <- tmp[order(tmp$logFC,decreasing = T),][1:200,]
basal.marker <- tmp[order(tmp$logFC,decreasing = F),][1:200,]
templates.luminal.basal <- data.frame(probe=c(rownames(luminal.marker),rownames(basal.marker)),Description="na",class=rep(c("Luminal","Basal"),c(200,200)))

indata <- as.data.frame(combat.UTUC.logcountsNorm)
rownames(indata) <- Ginfo[rownames(indata),"genename"]
indata <- t(scale(t(indata),center = T,scale = T))
ntp.utuc <- ntp(indata, templates.luminal.basal, nPerm=1000,distance = "cosine")
table(ntp.utuc$prediction)

###################################################################################
### use NTP to test if UTUC samples belongs to luminal papillary or infiltrated ###
tmp <- Sinfo$TCGA_Subtype; names(tmp) <- rownames(Sinfo)
tmp <- tmp[tmp %in% c("Luminal_papillary","Luminal_infiltrated")]
complist <- createList.luminal(group=tmp)
tmp1 <- data.frame(tmp,sample=names(tmp)); colnames(tmp1)[1] <- "group"
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[, names(tmp)], tailrows, Groupinfo=tmp1, features=Mids, featType="mRNA_BLCA", complist,PASSFlag=PASSFlag.mRNA, overwt=F)

tmp <- read.table(file.path(res.path,"mRNA_BLCA_edgeR_test_result.Luminal_papillary_vs_Luminal_infiltrated.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
papillary.marker <- tmp[order(tmp$logFC,decreasing = T),][1:200,]
infiltrated.marker <- tmp[order(tmp$logFC,decreasing = F),][1:200,]

papillary.marker2 <- c("KRT20","PPARG","FOXA1","GATA3","SNX31","UPK1A","UPK2","FGFR3")
infiltrated.marker2 <- c("PGM5","DES","C7","SFRP4","COMP","SGCD","ZEB1","ZEB2","TWIST1","CDH2","CLDN3","CLDN4","CLDN7")

templates.papillary.infiltrated <- data.frame(probe=c(rownames(papillary.marker),rownames(infiltrated.marker)),Description="na",class=rep(c("Luminal_papillary","Luminal_infiltrated"),c(200,200)))
#templates2 <- data.frame(probe=c(papillary.marker2,infiltrated.marker2),Description="na",class=rep(c("Luminal_papillary","Luminal_infiltrated"),c(8,13)))

# load different data
# indata <- read.table(file.path(res.path,"UTUC_BLCA.lowcut1.lowpct0.9.FPKM.mRNA.LowExpfiltered.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# indata <- indata[,427:446]

# indata <- as.data.frame(combat.UTUC.logFPKM)
# indata <- as.data.frame(combat.UTUC.FPKM.HUGO)
indata <- as.data.frame(combat.UTUC.logcountsNorm)
# indata <- as.data.frame(combat.UTUC.countsNorm)

indata[1:3,1:3]
range(indata)

# postive
# indata <- pmax(as.matrix(indata),0)*2
# range(indata)

# transfer ensembl to hugo
rownames(indata) <- Ginfo[rownames(indata),"genename"]

# log transformation
# indata <- log2(indata + 1)
# range(indata)

# scale data
# indata <- standarize.fun(indata,scaleFlag = T,centerFlag = T,halfwidth = 3)

# ntp
indata <- t(scale(t(indata),center = T,scale = T))
ntp.utuc <- ntp(indata, templates.papillary.infiltrated, nPerm=1000,distance = "cosine")
tmp <- UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes"),]; rownames(tmp) <- tmp$DetailID
ntp.utuc$MI <- tmp[rownames(ntp.utuc),"Type"]
ntp.utuc$FGFR3 <- tmp[rownames(ntp.utuc),c("FGFR3")]
ntp.utuc$FGFR3_mut <- tmp[rownames(ntp.utuc),c("FGFR3_mut")]

write.table(ntp.utuc,file.path(res.path,"NTP for luminal papillary or infiltrated.txt"),sep = "\t",row.names = T,col.names = NA)

table(ntp.utuc$prediction)
table(ntp.utuc$prediction,ntp.utuc$MI)
fisher.test(table(ntp.utuc$prediction,ntp.utuc$MI)) # 0.008127
fisher.test(matrix(c(9,0,4,5),nrow = 2,byrow = T)) # 0.02941

table(ntp.utuc$prediction,ntp.utuc$FGFR3)
table(ntp.utuc$prediction,ntp.utuc$FGFR3_mut)
fisher.test(table(ntp.utuc$prediction,ntp.utuc$FGFR3_mut)) # 0.004525


# indata <- combat.logFPKM[,rownames(annCol[which(annCol$Source == "UTUC" | annCol$TCGA_Subtype %in% c("Luminal_papillary","Luminal_infiltrated")),])]
# rownames(indata) <- Ginfo[rownames(indata),"genename"]
# 
# plotdata <- indata[intersect(rownames(indata),c(rownames(papillary.marker),rownames(infiltrated.marker))),]
# plotdata <- standarize.fun(plotdata,halfwidth = 3)
# hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "ward.D")
# pheatmap(plotdata,
#          cluster_cols = hcs,
#          color = greenred(64),
#          annotation_col = annCol[colnames(plotdata),],
#          annotation_colors = annColors)
# 
# indata <- log2(nneg.combat.FPKM[,rownames(annCol[which(annCol$Source == "UTUC"),])] + 1)
# rownames(indata) <- Ginfo[rownames(indata),"genename"]
# indata <- t(scale(t(indata),center = T,scale = T))
# templates <- data.frame(probe=c(rownames(papillary.marker),rownames(infiltrated.marker)),Description="na",class=rep(c("Luminal_papillary","Luminal_infiltrated"),c(200,200)))
# templates3 <- data.frame(probe=c("RNF128","GPD1L","RAB15","FAM174B","CYP4B1","ADIRF",
#                                 "PPARG","FBP1","GATA3","PPFIBP2","TOX3","SPINK1",
#                                 "CAPN5","UPK2","UPK1A","SCNN1G","SCNN1B","TBX2",
#                                 "HMGCS2","VGLL1","PLEKHG6","GDPD3","SEMA5A","TMEM97",
#                                 "SLC9A2","BHMT","CYP2J2","GAREM","SLC27A2","TUBB6",
#                                 #
#                                 "CHST15","MSN","CD14","ALOX5AP","GLIPR1","EMP3",
#                                 "PRRX1","FAP","PALLD","PDGFC","PRKCDBP",
#                                 "MT2A","MT1X","CDK6","AYNAK2"),Description="na",class=rep(c("Luminal","Basal"),c(30,15)))
# templates4 <- data.frame(probe = c("KRT20","PPARG","FOXA1","GATA3","SNX31","UPK1A","UPK2","FGFR3",
#                                     #
#                                     "CD44","KRT6A","KRT5","KRT14","COL17A1"),
#                          Description="na",class=rep(c("Luminal","Basal"),c(8,5)))

# use consensusMIBC
library(consensusMIBC)
indata <- as.data.frame(combat.UTUC.logcountsNorm)
rownames(indata) <- Ginfo[rownames(indata),"genename"]
range(indata)
# scale
# indata <- t(scale(t(indata),center = T,scale = F))
# log
# indata <- log2(indata + 1)

tmp <- UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes"),]; rownames(tmp) <- tmp$DetailID
MIBC <- getConsensusClass(indata, minCor = .2, gene_id = "hgnc_symbol")
MIBC$MI <- tmp[rownames(MIBC),"Type"]
MIBC$FGFR3 <- tmp[rownames(MIBC),c("FGFR3")]
MIBC$FGFR3_mut <- tmp[rownames(MIBC),c("FGFR3_mut")]
table(MIBC$consensusClass,MIBC$MI)
plotCorrelations(MIBC)
write.table(MIBC,file.path(res.path,"consensusMIBC.outTab.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

rownames(MIBC) <- sapply(strsplit(rownames(MIBC),"_"),"[",2)
rownames(MIBC) <- paste0("0",rownames(MIBC))
rownames(MIBC) <- substr(rownames(MIBC),nchar(rownames(MIBC))-2,nchar(rownames(MIBC)))

consensus.mut <- mygene
consensus.mut <- consensus.mut[,intersect(colnames(consensus.mut),rownames(annCol.rna))]

type <- c()
for (i in 1:nrow(consensus.mut)) {
  tmp <- as.character(consensus.mut[i,])
  type <- unique(c(type,tmp))
}
cat(paste0("You have a total of ",length(setdiff(type,""))," types of mutation!\nPlease check as the following:\n"))
print(setdiff(type,""))# this is all the mutation type now you have and you have to set the color one by one

mycol <- RColorBrewer::brewer.pal(n = 10, name = 'Paired')[c(2,1,4,3,6,5)]
lightgrey <- "#dcddde" # background color was set to grey
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = lightgrey, col = NA))
  },
  Nonsynonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[1], col = "black")) # e.g., the first color in mycol vector is for Missense_Mutation
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[2], col = "black")) 
  },
  Frameshift_Insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[3], col = "black")) 
  },
  Frameshift_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[4], col = "black")) 
  },
  Stopgain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[5], col = "black")) 
  },
  Splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[6], col = "black"))
  }
)

# set color again for mutation annotation which should exactely the same with above
col = c("Nonsynonymous" = mycol[1], 
        "Synonymous" = mycol[2], 
        "Frameshift_Insertion" = mycol[3], 
        "Frameshift_Deletion" = mycol[4], 
        "Stopgain" = mycol[5],
        "Splicing" = mycol[6])
mycol2 <- c("#B24745","#DF8F44","#374E55","#00A1D5")

p1 <- oncoPrint(consensus.mut[,rownames(my_ann[which(my_ann$consensusClass == "LumP"),,drop =F])], # cluster1 oncoprint
                alter_fun = alter_fun,
                col = col)
p2 <- oncoPrint(consensus.mut[,rownames(my_ann[which(my_ann$consensusClass == "LumU"),,drop =F])], # cluster1 oncoprint
                alter_fun = alter_fun,
                col = col)

my_ann <- data.frame(MIBC[colnames(consensus.mut),"consensusClass",drop = F])
my_ann <- my_ann[c(rownames(my_ann[which(my_ann$consensusClass == "LumP"),,drop =F])[p1@column_order],
                   rownames(my_ann[which(my_ann$consensusClass == "LumU"),,drop =F])[p2@column_order],
                   rownames(my_ann[which(my_ann$consensusClass == "Ba/Sq"),,drop =F]),
                   rownames(my_ann[which(my_ann$consensusClass == "Stroma-rich"),,drop =F])),,drop = F] # sort sample by the suboncoprint order

#my_ann$consensusClass <- factor(my_ann$consensusClass,levels = c("LumP","LumU","Ba/Sq","Stroma-rich"))
#my_ann <- my_ann[order(my_ann$consensusClass),,drop = F]
my_annotation = HeatmapAnnotation(df = data.frame(my_ann), # annotation data frame
                                  col = list(consensusClass = c("LumP" = red, "Stroma-rich" = blue, "Ba/Sq" = seagreen,"LumU" = yellow)))

p <- oncoPrint(consensus.mut[,rownames(my_ann)], # this is the detailed mutation matrix
               alter_fun = alter_fun,  # this is the alteration function we already set
               col = col, # this is the color for mutation
               bottom_annotation = my_annotation, # this is the annotation bar which was put bottom
               column_order = rownames(my_ann), # sort your sample as your subtype order
               #row_order = mut.order, # sort mutation as mutation frequency
               show_pct = T, #show percentage in left
               column_title = "", # no title shown
               show_heatmap_legend=T, # show legend in the oncoprint
               #column_split = my_ann$Subtype,
               # some detailed size below and you may not have to change it
               top_annotation = NULL,
               #right_annotation = NULL,
               column_title_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8)) 

# output this oncoprint to pdf in the local working path
pdf(file.path(fig.path,"consensusMIBC oncoprint.pdf"),width = 7.5,height = 3)
p
invisible(dev.off())

##############################################################
### test mutation load regarding significant mutated genes ###

# FGFR3
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$FGFR3 == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$FGFR3 == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.5897

# KDM6A
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$KDM6A == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$KDM6A == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.439

# MLL2
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$MLL2 == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$MLL2 == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.04132

# MLL3
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$MLL3 == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$MLL3 == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.8063

# ZFP36L1
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$ZFP36L1 == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.05504

# ARID1A
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$ARID1A == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$ARID1A == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.01433

# STAG2
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$STAG2 == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$STAG2 == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.8238

# CRIPAK
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$CRIPAK == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$CRIPAK == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.8072

# GANAB
tmp1 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$GANAB == "1"),"ID"]),"T")]
tmp2 <- mutload[paste0(gsub("-","_",UTUC.annotation[which(UTUC.annotation$GANAB == "0"),"ID"]),"T")]
wilcox.test(tmp1,tmp2) # 0.02213

#########################################################################
### check if TET1 and TET2 is overexpressed in FGFR3 mutation samples ###
tmp1 <- as.numeric(combat.UTUC.FPKM.HUGO["TET1",UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "1" & UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"]])
tmp2 <- as.numeric(combat.UTUC.FPKM.HUGO["TET1",UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "0" & UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"]])
wilcox.test(tmp1,tmp2)

tmp <- FPKM.BLCA; colnames(tmp) <- substr(colnames(tmp),1,12)
tmp1 <- as.numeric(tmp["TET1",intersect(colnames(tmp),rownames(Sinfo.BLCA[which(Sinfo.BLCA$`mutation in FGFR3` == "yes"),]))])
tmp2 <- as.numeric(tmp["TET1",intersect(colnames(tmp),rownames(Sinfo.BLCA[which(Sinfo.BLCA$`mutation in FGFR3` == "no"),]))])
wilcox.test(tmp1,tmp2,alternative = "less")

tmp1 <- as.numeric(combat.UTUC.FPKM.HUGO["TET2",UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "1" & UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"]])
tmp2 <- as.numeric(combat.UTUC.FPKM.HUGO["TET2",UTUC.annotation[which(UTUC.annotation$FGFR3_mut == "0" & UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"]])
wilcox.test(tmp1,tmp2)

tmp <- FPKM.BLCA; colnames(tmp) <- substr(colnames(tmp),1,12)
tmp1 <- as.numeric(tmp["TET2",intersect(colnames(tmp),rownames(Sinfo.BLCA[which(Sinfo.BLCA$`mutation in FGFR3` == "yes"),]))])
tmp2 <- as.numeric(tmp["TET2",intersect(colnames(tmp),rownames(Sinfo.BLCA[which(Sinfo.BLCA$`mutation in FGFR3` == "no"),]))])
wilcox.test(tmp1,tmp2,alternative = "greater")

blca.fgfr3.mut <- Sinfo.BLCA$`mutation in FGFR3`; names(blca.fgfr3.mut) <- rownames(Sinfo.BLCA)
blca.fgfr3.mut <- ifelse(blca.fgfr3.mut == "yes","Mutation","Normal")
tmp <- countsTable; colnames(tmp) <- substr(colnames(tmp),1,12); tmp <- tmp[,intersect(colnames(tmp),names(blca.fgfr3.mut))]
blca.fgfr3.mut <- blca.fgfr3.mut[colnames(tmp)]
complist <- createList.mut(group=blca.fgfr3.mut)
blca.fgfr3.mut <- data.frame("group"=blca.fgfr3.mut,samples=names(blca.fgfr3.mut),row.names = names(blca.fgfr3.mut))
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=tmp, tailrows, Groupinfo=blca.fgfr3.mut, features=Mids, featType="mRNA_FGFR3_BLCA", complist,PASSFlag=PASSFlag.mRNA, overwt=TRUE)

utuc.fgfr3.mut <- UTUC.annotation$FGFR3_mut; names(utuc.fgfr3.mut) <- UTUC.annotation$DetailID
utuc.fgfr3.mut <- utuc.fgfr3.mut[utuc.fgfr3.mut %in% c("0","1")]
utuc.fgfr3.mut <- ifelse(utuc.fgfr3.mut == "1","Mutation","Normal")
utuc.fgfr3.mut <- utuc.fgfr3.mut[intersect(names(utuc.fgfr3.mut),UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes"),"DetailID"])]
complist <- createList.mut(group=utuc.fgfr3.mut)
utuc.fgfr3.mut <- data.frame("group"=utuc.fgfr3.mut,samples=names(utuc.fgfr3.mut),row.names = names(utuc.fgfr3.mut))
rowids <- intersect( Mids, names(PASSFlag.mRNA) )
twoclassedgeR(res.path, Ginfo, countsTable=countsTable[,rownames(utuc.fgfr3.mut)], tailrows, Groupinfo=utuc.fgfr3.mut, features=Mids, featType="mRNA_FGFR3_UTUC", complist,PASSFlag=PASSFlag.mRNA, overwt=TRUE)

##############################################################################
# test if CNA associated with clinical features and frequently mutated genes #
cna.exome.path <- "H:/Jan2020/UTUC_2020/UTUC_BLCA/CNA"
cna.exome <- read.table(file.path(cna.exome.path,"UTUC.rm.all_lesions.conf_75.txt"),sep = "\t",check.names = F,header = T,row.names = NULL,stringsAsFactors = F)
cna.binary <- cna.exome[1:10,10:38]; rownames(cna.binary) <- paste0(cna.exome$Descriptor[1:10],"_",1:10)
cna.binary[cna.binary>1] <- 1
cna.binary <- cna.binary[-6,]

amp <- rownames(cna.binary)[1:6]
del <- rownames(cna.binary)[7:9]

cna_cli_mut <- as.data.frame(t(cna.binary))
rownames(cna_cli_mut) <- sapply(strsplit(rownames(cna_cli_mut),"_"),"[",2)
rownames(cna_cli_mut) <- paste0("0",rownames(cna_cli_mut))
rownames(cna_cli_mut) <- substr(rownames(cna_cli_mut),nchar(rownames(cna_cli_mut))-2,nchar(rownames(cna_cli_mut)))

tmp <- cbind.data.frame(cna_cli_mut,UTUC.annotation[rownames(cna_cli_mut),])
tmp <- tmp[,c(1:9,39:48,19,20,23,33:37)]
tmp$`pathological stage` <- ifelse(tmp$`pathological stage` %in% c("pTa","pT1"),"a1","234")
# check amplification
for (a in amp) {
  for (g in setdiff(sig.mut.gene,"ARID1B")) {
    cat(paste0(a," and ",g,"\n"))
    tmp1 <- tmp[,c(a,g)]
    tmp1 <- tmp1[which(tmp[,g] != "N/A"),]
    tmp1[,a] <- ifelse(tmp1[,a] == 1,"amp","wild")
    tmp1[,g] <- ifelse(tmp1[,g] == 1,"mut","wild")
    #if(fisher.test(table(tmp1[,a],tmp1[,g]))$p.value < 0.1) {
      print(table(tmp1[,g],tmp1[,a]))
      print(fisher.test(table(tmp1[,g],tmp1[,a]))$p.value)
    #}
    cat("\n")
  }
  
  for (c in c("Type","GRADE","pathological stage")) {
    cat(paste0(a," and ",c,"\n"))
    tmp1 <- tmp[,c(a,c)]
    tmp1[,a] <- ifelse(tmp1[,a] == 1,"amp","wild")
    #if(fisher.test(table(tmp1[,a],tmp1[,c]))$p.value < 0.1) {
      print(table(tmp1[,c],tmp1[,a]))
      print(fisher.test(table(tmp1[,c],tmp1[,a]))$p.value)
    #}
    cat("\n")
  }
  
  # tmp1 <- tmp[,c(a,"OS","OS.time","PFS","PFS2","PFS.time")]
  # fitd=survdiff(Surv(OS.time, OS)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.os <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # fitd=survdiff(Surv(PFS.time, PFS)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.pfs <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # fitd=survdiff(Surv(PFS.time, PFS2)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.pfs2 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # cat(paste0("OS: ",p.os,"; PFS: ",p.pfs,"; PFS2: ",p.pfs2,"\n"))
}

# check deletion
for (a in del) {
  for (g in setdiff(sig.mut.gene,"ARID1B")) {
    cat(paste0(a," and ",g,"\n"))
    tmp1 <- tmp[,c(a,g)]
    tmp1 <- tmp1[which(tmp[,g] != "N/A"),]
    tmp1[,a] <- ifelse(tmp1[,a] == 1,"del","wild")
    tmp1[,g] <- ifelse(tmp1[,g] == 1,"mut","wild")
    #if(fisher.test(table(tmp1[,a],tmp1[,g]))$p.value < 0.1) {
      print(table(tmp1[,g],tmp1[,a]))
      print(fisher.test(table(tmp1[,g],tmp1[,a]))$p.value)
    #}
    cat("\n")
  }
  
  for (c in c("Type","GRADE","pathological stage")) {
    cat(paste0(a," and ",c,"\n"))
    tmp1 <- tmp[,c(a,c)]
    tmp1[,a] <- ifelse(tmp1[,a] == 1,"del","wild")
    #if(fisher.test(table(tmp1[,a],tmp1[,c]))$p.value < 0.1) {
      print(table(tmp1[,c],tmp1[,a]))
      print(fisher.test(table(tmp1[,c],tmp1[,a]))$p.value)
    #}
    cat("\n")
  }
  
  # tmp1 <- tmp[,c(a,"OS","OS.time","PFS","PFS2","PFS.time")]
  # fitd=survdiff(Surv(OS.time, OS)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.os <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # fitd=survdiff(Surv(PFS.time, PFS)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.pfs <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # fitd=survdiff(Surv(PFS.time, PFS2)~ tmp1[,a], data=tmp1, na.action=na.exclude)
  # p.pfs2 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
  # cat(paste0("OS: ",p.os,"; PFS: ",p.pfs,"; PFS2: ",p.pfs2,"\n"))
}

tmp1 <- tmp[,1:9]
tmp2 <- tmp[,c("FGFR3","TP53","ARID1A","MLL3")]
# perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
# "less" for mutual-exclusivity analysis, or "greater" for co-occurrence analysis 

res <- NULL
for(i in 1:ncol(tmp1)){
  for(j in 1:ncol(tmp2)){
    tab <- data.frame("cna"=as.character(tmp1[,i]),"mut"=as.character(tmp2[,j]),stringsAsFactors = F)
    tab <- tab[which(tab$mut != "N/A"),]
    tab <- table(tab$cna, tab$mut)
    f <- fisher.test(tab) 
    # deal with zero in 2x2 cell
    if(f$estimate == 0) {
      f <- fisher.test(tab,alternative = "less")
      res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp1)[i],geneB=colnames(tmp2)[j],
                                                   Neither=tab[1,1],
                                                   AnotB=tab[2,1],
                                                   BnotA=tab[1,2],
                                                   Both=tab[2,2],
                                                   oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
    } else {
      if(is.infinite(f$estimate)) {
        f <- fisher.test(tab + 1)
        if(log(f$estimate) > 0) {
          f <- fisher.test(tab + 1,alternative = "greater")
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp1)[i],geneB=colnames(tmp2)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "greater")$p.value),stringsAsFactors = F)
        } else{
          f <- fisher.test(tab + 1,alternative = "less")
          res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp1)[i],geneB=colnames(tmp2)[j],
                                                       Neither=tab[1,1],
                                                       AnotB=tab[2,1],
                                                       BnotA=tab[1,2],
                                                       Both=tab[2,2],
                                                       oddsRatio=f$estimate,pvalue=fisher.test(tab,alternative = "less")$p.value),stringsAsFactors = F)
        }
        
      } else {
        if(log(f$estimate) > 0) {
          f <- fisher.test(tab,alternative = "greater")
        } else{f <- fisher.test(tab,alternative = "less")}
        res <- rbind.data.frame(res,cbind.data.frame(geneA=colnames(tmp1)[i],geneB=colnames(tmp2)[j],
                                                     Neither=tab[1,1],
                                                     AnotB=tab[2,1],
                                                     BnotA=tab[1,2],
                                                     Both=tab[2,2],
                                                     oddsRatio=f$estimate,pvalue=f$p.value),stringsAsFactors = F)
      }   
    }
  }
}

# some formatting
res <- as.data.frame(res)
res$Tendency <- ifelse(as.numeric(res$oddsRatio) > 1,"Co-occurrence","Mutual-exclusivity")
res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
res$log2OR <- log2(as.numeric(as.character(res$oddsRatio))+0.1)
res$pvalue <- as.numeric(as.character(res$pvalue))
# use p.adjust to correct for multi testing using a FDR
res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
# change the FDR in labels for plotting
res$stars <- cut(res$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ".",""))
# plot with ggplot2
write.table(res,file.path(res.path,"Mutual_exclusivity_of cna and mut.txt"),sep = "\t",row.names = F)

p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = log2OR),colour = "white") + scale_fill_gradient2(low = "darkblue",mid = lightgrey,high = "darkred",midpoint=0) + 
  geom_text(aes(label=stars), color="white", size=5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank())
ggsave(file.path(fig.path,"Mutual_exclusivity_of cna and mut.pdf"),width = 5,height = 2.5)

dat <- tmp[,c(20,4,5,9)]
dat <- dat[order(dat$Type),]
dat$Type <- ifelse(dat$Type == "muscle-invasive","Muscle-invasive","Non-muscle invasive")
p.chisq <- c()
for (i in 2:ncol(dat)) { # 第一列是亚型，第二列例文没有做独立性检验，所以从第三列开始
  tmp <- table(dat[,1],dat[,i])
  p.chisq <- c(p.chisq,round(fisher.test(tmp)$p.value,3)) # 使用卡方检验
}
head(p.chisq)
ann_labels <- paste0(c("8q22.3 gain", "12q15 gain","11p15.1 loss"),"\nP",ifelse(p.chisq < 0.001," < 0.001",paste0(" = ",p.chisq)))
ann_labels <- c(NA,ann_labels) 

### 设置颜色 ###
darkblue1  <- "#292C7C"
lightblue1 <- "#97B9D7"
blue1      <- "#4582B5"
white1     <- "#FFFFFF"
dwhite1    <- "#B6D1E8"

pdf(file.path(fig.path,"AssociationHeatmap_cna_type.pdf"),width = 7.5,height = 3.5)

par(bty="n", mgp = c(2,0.15,0), mar = c(3.1,2.1,3.1,8.1),tcl=-.25, font.main=3) # 张开画布，设置基本参数
par(xpd=NA)
# 生成一个空白图
plot(c(0,nrow(dat)), c(0,ncol(dat)), # 产生左下和左上两个点，张开画布
     col = "white", # 点设置为白色
     xlab = "", xaxt = "n", # 不显示x坐标轴
     ylab = "", yaxt = "n") # 不显示y坐标轴)
title(bquote(underline("UTUC cohort")), adj = 0, line = 0) # 左对齐title，并添加下划线

# 生成纵坐标
axis(4, at = 0.475:(ncol(dat)-0.475), # 对应每行0.95的高度，因此中点在0.475的位置
     labels = ann_labels[4:1], las = 1, # 生成y坐标轴刻度，并位于每个间隔的中心位置
     col = "white") # 隐藏坐标轴

# 生成横坐标
#axis(1, at = c(seq(0, nrow(dat), 20), nrow(dat)), labels = T) # 生成x坐标轴刻度

# 产生对应颜色
input_matrix <- as.matrix(dat) 
col.mat <- matrix(NA, byrow = T, ncol = nrow(input_matrix), nrow = ncol(input_matrix)) # 生成空的颜色矩阵
rownames(col.mat) <- colnames(input_matrix)

for (i in 2:nrow(col.mat)) { # 由于例文第一列数据的颜色不同，所以单独分配，这里从第二列开始
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- ifelse(input_matrix[j,i] == 1,blue1,white1) # 根据输入数据（二分类数据）给对应的颜色矩阵赋值
  }
}
col.mat[1,] <- ifelse(input_matrix[,1] == "Muscle-invasive",darkblue1,lightblue1) # 分配第一列的颜色

# 通过矩形块产生热图
x_size <- nrow(input_matrix)
y_size <- ncol(input_matrix)

my_xleft = rep(0:(x_size-1), each = y_size) # 产生各个矩形的左x点，注意是从上到下从左往右的顺序
my_xright = my_xleft + 1 # 产生各个矩形的右x点
my_ybottom = rep((y_size-1):0, x_size) # 产生各个矩形的下y点
my_ytop = my_ybottom + 0.95 # 产生各个矩形的上y点, 每行余留0.05的空隙
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col= col.mat, # 填充颜色
     border = dwhite1) # 矩形的边界定义为深白色

# 根据原文图样，修饰第一行，隐去间隔
my_xleft = rep(0:(x_size-1)) # 产生各个矩形的左x点
my_xright = my_xleft + 1 # 产生各个矩形的右x点
my_ybottom = rep(y_size-1, x_size) # 产生各个矩形的下y点
my_ytop = my_ybottom + 0.95 # 产生各个矩形的上y点, 每行余留0.05的空隙

rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col= col.mat[1,], # 填充颜色
     border = NA) # 取消矩形的边界

# 第一行的文字
text(nrow(input_matrix)/4, ncol(input_matrix) - 0.55,
     substitute(italic("Muscle invasive")), cex = 1.3, col = "white") # 填充标签，注意位置，斜体和字体颜色
text(nrow(input_matrix)/4*3, ncol(input_matrix) - 0.55,
     substitute(italic("Non-mucle invasive")), cex = 1.3, col = "black") # 填充标签，注意位置，斜体和字体颜色

# 图例
legend("topright", 
       inset = c(-0.18,0.07), # 微调图例的位置，让图例在画布外侧
       xjust = 1, # 左对齐
       cex = 1,
       y.intersp = 1.3, # 两个图例的间距大一点
       fill = c(blue1, white1), # 图例颜色
       legend = c("Presence", "Absence"), 
       bty = "n") # 不绘制外部矩形

invisible(dev.off())
################################################################
### check association between mutations and BLCA DMB cluster ###
load(file.path(res.path,"sp.cluster2.tcga.ans.rda"))
sig.mut <- c("FGFR3","KDM6A","KMT2D","KMT2C","ZFP36L1","ARID1A","STAG2","TP53","CRIPAK","GANAB")
tmp <- read.table(file.path(data.path,"BLCA_UTUC_sig_mut10.tsv"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
colnames(tmp) <- paste0("BLCA",substr(colnames(tmp),8,15))
tmp[is.na(tmp)] <- ""
tmp[tmp != ""] <- 1
tmp[tmp == ""] <- 0
library(tableone)
tmp <- cbind.data.frame(t(tmp[names(sp.cluster2.tcga.ans$group)]),group = sp.cluster2.tcga.ans$group)
blca.mut10 <- tmp[,1:10]

stabl <- CreateTableOne(vars=sig.mut,strata="group",
                        data=tmp,factorVars=sig.mut)
print(stabl,showAllLevels = TRUE)
write.table(print(stabl,showAllLevels = TRUE),
            file.path(res.path,"BLCA_UTUC_sig_mut10_statistic_table.txt"),sep = "\t",quote = F)

tcga.mut.swisnf <- read.table(file.path(data.path,"BLCA_SWI_SNF20.tsv"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1,fill = T)
colnames(tcga.mut.swisnf) <- paste0("BLCA",substr(colnames(tcga.mut.swisnf),8,15))
tcga.mut.swisnf[is.na(tcga.mut.swisnf)] <- ""
tcga.mut.swisnf[tcga.mut.swisnf != ""] <- 1
tcga.mut.swisnf[tcga.mut.swisnf == ""] <- 0
tmp <- rownames(tcga.mut.swisnf)
tcga.mut.swisnf <- as.data.frame(sapply(tcga.mut.swisnf,as.numeric)); rownames(tcga.mut.swisnf) <- tmp
tcga.mut.swisnf <- as.data.frame(t(tcga.mut.swisnf))
tcga.mut.swisnf$swisnf <- ifelse(rowSums(tcga.mut.swisnf) >= 1,1,0)
tcga.mut.swisnf <- cbind.data.frame(tcga.mut.swisnf[names(sp.cluster2.tcga.ans$group),],group = sp.cluster2.tcga.ans$group)
fisher.test(table(tcga.mut.swisnf$swisnf,tcga.mut.swisnf$group)) # 0.3744
# cluster1 cluster2
# 0      101       97
# 1      119       95
blca.swi.snf <- tcga.mut.swisnf

# blca heatmap
load(file = file.path(res.path,"indata_BLCA_heatmap_SupervisedCluster.rda"))
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")

annCol.blca <- data.frame(SpCluster = sp.cluster2.tcga.ans$group,row.names = names(sp.cluster2.tcga.ans$group),stringsAsFactors = F)
annCol.blca <- cbind.data.frame(annCol.blca,blca.mut10[rownames(annCol.blca),],blca.swi.snf[rownames(annCol.blca),"swisnf",drop = F])
annCol.blca$SpCluster <- ifelse(annCol.blca$SpCluster == "cluster1","C1","C2")
annCol.blca <- annCol.blca[,c("SpCluster",names(sort(colSums(sapply(blca.mut10,as.numeric)))),"swisnf")]
annCol.blca$swisnf <- as.character(annCol.blca$swisnf)
annColors.blca <- list()
annColors.blca[["SpCluster"]] <- c("C1"=red,"C2"=blue)
annColors.blca[["FGFR3"]] <- c("1"="black","0"="white")
annColors.blca[["KDM6A"]] <- c("1"="black","0"="white")
annColors.blca[["KMT2D"]] <- c("1"="black","0"="white")
annColors.blca[["KMT2C"]] <- c("1"="black","0"="white")
annColors.blca[["ZFP36L1"]] <- c("1"="black","0"="white")
annColors.blca[["ARID1A"]] <- c("1"="black","0"="white")
annColors.blca[["STAG2"]] <- c("1"="black","0"="white")
annColors.blca[["TP53"]] <- c("1"="black","0"="white")
annColors.blca[["CRIPAK"]] <- c("1"="black","0"="white")
annColors.blca[["GANAB"]] <- c("1"="black","0"="white")
annColors.blca[["swisnf"]] <- c("1"="black","0"="white")

mycol <- colorRampPalette(brewer.pal(11,"Spectral"))(11)[11:1]
pdf(file = file.path(fig.path,"BLCA3_heatmap_SupervisedCluster_UTUC_C1vsC2_DMPs.pdf"),height = 8,width = 8)
pheatmap(as.matrix(indata), 
         cluster_rows=hcg, 
         cluster_cols=as.hclust(sp.cluster2.tcga.ans$dendro), 
         annotation_col = annCol.blca[colnames(indata),],
         annotation_colors = annColors.blca, 
         #color=c("blue","blue","green","yellow","red","red"), 
         color=NMF:::ccRamp(c("#0074FE","#96EBF9","#FEE900","#F00003"),n = 64),
         show_colnames = F,show_rownames = F)
invisible(dev.off())

#################################################
### MCPcounter between UTUC EpiC-high and low ###
tmp <- as.data.frame(UTUC.FPKM)
tmp$genename <- Ginfo[rownames(tmp),"genename"] 
tmp <- tmp[!duplicated(tmp$genename),]; rownames(tmp) <- tmp$genename; tmp <- tmp[,-ncol(tmp)]
colnames(tmp) <- sapply(strsplit(colnames(tmp),"_"),"[",2)
colnames(tmp) <- paste0("0",colnames(tmp))
colnames(tmp) <- substr(colnames(tmp),nchar(colnames(tmp))-2,nchar(colnames(tmp)))
epich.exp <- rownames(UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes" & UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C1"),])
epicl.exp <- rownames(UTUC.annotation[which(UTUC.annotation$`RNAseq (please add)` == "yes" & UTUC.annotation$`Cluster DNA methylation no filter without normal` == "C2"),])

tmp <- log2(tmp[,c(epich.exp,epicl.exp)] + 1)
MCPscore.epicexp <- MCPcounter.estimate(expression = tmp,featuresType = "HUGO_symbols")
MCPscore.epicexp.std <- standarize.fun(MCPscore.epicexp,halfwidth = 2)
pheatmap(MCPscore.epicexp.std,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = 11,
         color = greenred(64))
wilcox.test(MCPscore.epicexp[8,epich.exp],MCPscore.epicexp[1,epicl.exp])

#save
save.image(file.path(tumor.path,"Malouf.RNA2.RData"))
