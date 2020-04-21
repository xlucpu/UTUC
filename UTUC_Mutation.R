workdir <- "F:/Project/UTUC_BLCA/exome UTUC"; setwd(workdir)
fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(deconstructSigs)
library(ClassDiscovery)
library(NMF)
library(dendsort)
library(gplots)
library(RColorBrewer)
source("./getNearbyBases_XS.R")

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
#https://stackoverflow.com/questions/52310844/how-to-reshape-data-from-wide-to-long-with-multipe-variables
#################################
#### 1. Load MAF DATA ##########
# maf <- read.table("F:/Project/UTUC_BLCA/exome UTUC/utuc_EXOME_mutation_data.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
# maf$Chr <- paste0("chr",maf$Chr)
# sel_mut <- maf %>%
#   filter(MutStatus %in% c("nonsynonymous SNV","stopgain SNV","stoploss SNV"),
#          AlleleB %in% c("A","T","C","G"))
# # reshape data
# sel_mut <- sel_mut %>%
#   gather( sample, sample_status, colnames(sel_mut)[8:37]) %>%
#   filter( sample_status == 1 )
# 
# ################################################
# ## 2. Calculate weights of mutation signature ##
# 
# # 2.1 Build the base substitution column.
# sel_mut <- sel_mut %>%
#   mutate(baseSubstitution = factor(str_c(AlleleA,">", AlleleB)),
#          baseSubstitution = fct_recode(baseSubstitution,
#                                        "T>G" = "A>C",
#                                        "T>C" = "A>G",
#                                        "T>A" = "A>T",
#                                        "C>T" = "G>A",
#                                        "C>G" = "G>C",
#                                        "C>A" = "G>T"),
#          baseSubstitution = factor(baseSubstitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G") )
# 
#   )
# # 2.2 Get the strand of genes
# g2strand <- getStrandByHugoSymbol(unique(sel_mut$GeneName),use.hg38 = F)
# table(g2strand, useNA = "ifany")
# # g2strand
# #         -    + <NA> 
# #   92 1226 1205    7
# 
# g2strand[is.na(g2strand)] <- ""
# sel_mut$strand <- g2strand[match(sel_mut$GeneName, names(g2strand))]
# write.table(sel_mut,"sel_mut.txt",sep = "\t",row.names = F)
# # 2.3 Get the neighboring bases
# input_dat <- data.frame(chr = sel_mut$Chr,
#                         position = sel_mut$Start,
#                         strand = sel_mut$strand,
#                         stringsAsFactors=F )
# 
# sel_mut$flankingSeq =  getFlankingSeq(input_dat,use.hg38 = F)
# table(sel_mut$flankingSeq)
# # AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT 
# # 14  14  17  18  25  25 100  13 273  34  26  20   8  14  18  13  13  22  20  28  26  38  98  35 169 126  99 103  10   9  13  16   8  12  10  20  31  39 119  32 
# # GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT 
# # 159  45  38  25   8   2  13   9   5  12   4  35 255 102 111 157 246  38  54  28   3   8   4   4  
# 
# snv_nb_sam_v <- sel_mut %>% mutate(
#   mut_nb = paste0(str_sub(flankingSeq, 1, 1), "(", 
#                   baseSubstitution, ")",
#                   str_sub(flankingSeq, 3, 3))
# ) %>%
#   group_by(sample, mut_nb) %>%
#   count()
# 
# snv_nb_sam_h <- snv_nb_sam_v %>% spread(key = mut_nb, value = n, fill = 0)
# dim(snv_nb_sam_h)
# #[1] 30 94
# 
# # 2.4 Calculate weights of gene signatures per sample
# sanger <- read.delim("signatures_probabilities.txt", header = TRUE, stringsAsFactors = FALSE)
# sanger <- sanger[1:96, 1:33]
# rownames(sanger) <- as.character(sanger$Somatic.Mutation.Type)
# sanger[1:3, 1:5]
# #         Substitution.Type Trinucleotide Somatic.Mutation.Type Signature.1  Signature.2
# # A[C>A]A               C>A           ACA               A[C>A]A 0.011098326 6.827082e-04
# # A[C>A]C               C>A           ACC               A[C>A]C 0.009149341 6.191072e-04
# # A[C>A]G               C>A           ACG               A[C>A]G 0.001490070 9.927896e-05
# 
# sanger_mat <- as.matrix(sanger[, 4:33])
# rownames(sanger_mat) <- as.character(sanger$Somatic.Mutation.Type)
# sanger_mat[1:3, 1:3]
# #         Signature.1  Signature.2 Signature.3
# # A[C>A]A 0.011098326 6.827082e-04  0.02217231
# # A[C>A]C 0.009149341 6.191072e-04  0.01787168
# # A[C>A]G 0.001490070 9.927896e-05  0.00213834
# 
# snv_nb_mat <- t(as.matrix(snv_nb_sam_h[, 2:94]))
# colnames(snv_nb_mat) <- snv_nb_sam_h$sample
# rownames(snv_nb_mat) <- gsub("\\(", "\\[", rownames(snv_nb_mat))
# rownames(snv_nb_mat) <- gsub("\\)", "\\]", rownames(snv_nb_mat))
# snv_nb_mat[1:4, 1:4]
# #         UTUC_12T UTUC_13T UTUC_14T UTUC_15T
# # A[C>A]A        0        0        3        0
# # A[C>A]C        0        0        0        0
# # A[C>A]G        0        0        0        0
# # A[C>A]T        0        0        1        0
# 
# sanger_mat <- sanger_mat[rownames(snv_nb_mat), ]
# sanger <- sanger[rownames(snv_nb_mat), ]
# 
# sanger_wt <- calculateGeneSignatureWeight(snv_nb_mat, sanger_mat)
# range(sanger_wt)
# #[1] -1.670344e-15  6.376708e-01
# 
# sanger_wt[1:5, 1:5]
# write.table(sanger_wt,"UTUC_exome_mutation_signature.txt",sep = "\t",row.names = T,col.names = NA)
# 
# nneg.sanger_wt <- as.data.frame(pmax(sanger_wt,0))
# range(nneg.sanger_wt)
# nneg.sanger_wt[1:5,1:5]
# 
# #estim.r <- nmf(nneg.sanger_wt, 2:3, nrun=10, seed=123456)
# mut.nmf <- nmf(t(nneg.sanger_wt),2)
# group <- predict(mut.nmf)
# 
# featureScore(mut.nmf)
# extractFeatures(mut.nmf,method="max")
# 
# plotdata <- as.data.frame(na.omit(standarize.fun(indata = t(sanger_wt),halfwidth = 1)))
# hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
# hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "ward.D")
# 
# outFigFile <- "heatmap_UTUC_exome_mutation_signature.pdf"
# pdf(outFigFile, height=6)
# par(mfrow=c(2,2))
# basismap(mut.nmf, cexCol = 0.6, cexRow = 0.6)
# coefmap(mut.nmf,cexCol = 0.6, cexRow = 0.3)
# consensusmap(mut.nmf,cexCol = 0.6, cexRow = 1.5)
# hv <- aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=NULL, annColors=NULL, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=6, cexCol = 0.6, cexRow = 0.6)
# invisible(dev.off())
# 
# #
# save.image("MutSig.RData")

#########################
#### deconstructSigs ####

maf <- read.table(file.path(res.path,"UTUC_exome_mutation_with_silence_curated_addARID1A_del01T_2019315_modified_annovar_longformat_forMutSigCV_simplified.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T)
maf$Chromosome <- paste0("chr",maf$Chromosome)
unique(maf$Variant_Classification)
# [1] "Missense_Mutation" "Synonymous"        "Nonsense_Mutation" "Frame_Shift_Ins"   "Frame_Shift_Del"  
# [6] "Splice"            "In_Frame_Del"      "In_Frame_Ins"      "Nonstop_Mutation"
#maf <- maf[-which(maf$Variant_Classification %in% c("Synonymous")),]
#maf <- maf[which(maf$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation")),]
maf <- maf[which(maf$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Synonymous")),]

# sample.mut.ref <- maf %>%
#   filter(MutStatus %in% c("nonsynonymous SNV","stopgain SNV","stoploss SNV"),
#          AlleleB %in% c("A","T","C","G"))
# reshape data
# sample.mut.ref <- sample.mut.ref %>% 
#   gather( sample, sample_status, colnames(sample.mut.ref)[8:37]) %>%
#   filter( sample_status == 1 )

sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele1")
write.table(sigs.input,file.path(res.path,"mutation.sig.input.snp.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

cut.off <- 0.06
mut.wt <- data.frame()
sigs.out.list <- list()
for (sample in rownames(sigs.input)) {
  tmp <- whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  #Plot output
  pdf(file.path(fig.path,paste0(sample,"_plotSignatures.pdf")))
  plotSignatures(tmp)
  invisible(dev.off())

  pdf(file.path(fig.path,paste0(sample,"_weightPie.pdf")))
  makePie(tmp)
  invisible(dev.off())
  # 
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt,tmp)
}
write.table(mut.wt,file.path(res.path,"mutation.snp.signature.weightMatrix.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

#trim mut.wt by k (10,20,30% etc..)
rownames(mut.wt) <- paste0(0,sapply(strsplit(rownames(mut.wt),"_"),"[",2))
rownames(mut.wt) <- substr(rownames(mut.wt),start = nchar(rownames(mut.wt))-2,stop = nchar(rownames(mut.wt)))
rownames(mut.wt) <- paste0("UTUC_",rownames(mut.wt))
k = 0.2
MB <- 1
mut.wt.trim <- mut.wt[,1:30]
mut.wt.trim[mut.wt.trim < k] <- 0
mut.wt.trim <- mut.wt.trim[,colSums(mut.wt.trim) > 0] #10signature for snp

mut.wt.trim.backup <- mut.wt.trim

mut.wt.trim <- mut.wt.trim * MB
mut.wt.trim$Subject <- rownames(mut.wt.trim)
mut.wt.trim <- mut.wt.trim[order(mut.wt.trim$Signature.13,mut.wt.trim$Signature.16,mut.wt.trim$Signature.1,decreasing = T),]
sample.level2 <- rownames(mut.wt.trim)
sig.level2 <- colnames(mut.wt.trim)[order(colSums(mut.wt.trim[,1:10]),decreasing = T)]

mut.wt.trim.long <- reshape(mut.wt.trim,
                            idvar = "Subject",
                            varying = list(1:10),
                            v.names = "Mutations_MB",
                            direction = "long")
mut.wt.trim.long$SigName <- rep(paste0("Signature.",c(1,2,3,5,6,7,13,16,20,26)),each=30)


mut.wt.trim.long$Subject <- factor(mut.wt.trim.long$Subject,levels = sample.level2)
mut.wt.trim.long$SigName <- factor(mut.wt.trim.long$SigName,levels = sig.level2)
mut.wt.trim.long$SigName <- gsub("Signature","Sig",mut.wt.trim.long$SigName)

mycol <- c(cyan,darkred,lightgreen,gold,purple,soil,nake,brown,grey,darkblue,cherry,seagreen,purple,white)
ggplot(mut.wt.trim.long, aes(x = Subject, y = Mutations_MB, fill = SigName, width=0.6)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=colorRampPalette(brewer.pal(11,"Spectral"))(12)) +
  scale_fill_manual(values=mycol) +
  theme(axis.text.x=element_text(angle=45, hjust=1,size = 5.5),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1.1)) +
  theme(axis.title.x = element_blank())+ theme(legend.text=element_text(size=6),legend.title = element_blank()) + scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) + ylab("Cumulativ Weight of Mutation Signature")
ggsave(filename = file.path(fig.path,"stacked.snp histogram for 0.2 trimed mutation signature derived by deconstructSigs.pdf"),height = 3.5)


# using trimed mut.wt.trim.backup to do nmf again

nmf.input <- t(mut.wt.trim.backup)
nmf.input <- nmf.input[,colSums(nmf.input) > 0] #10 signatures 30samples
#nmf.input <- nmf.input + 0.001
#determine optimal rank : 3
#estim.r <- nmf(nmf.input, 2:5, nrun=10, seed=123456)
ranks <- 2:5
estim <- lapply(ranks, function(r){
  fit <- nmf(nmf.input, r, nrun = 30, seed = 4, method = "lee")
  list(fit = fit, consensus = consensus(fit), .opt = "vp",coph = cophcor(fit))
})
names(estim) <- paste('rank', ranks)

pdf(file.path(fig.path,"consensusmap for rank discovert 10sig.pdf"))
consensusmap(sapply(estim, '[[', 'consensus', simplify = FALSE))
invisible(dev.off())

pdf(file.path(fig.path, "Cophenetic coefficient for nmf rank 10sig.pdf"))
par(cex.axis=1.5)
plot(ranks, sapply(estim, '[[', 'coph'), xlab="", ylab="", type="b", col="red", lwd=4,xaxt="n")
axis(side = 1, at=1:5)
title(xlab="number of clusters", ylab="Cophenetic coefficient", cex.lab=1.5)
invisible(dev.off())
silhouette.plot(nmf.input,2:5,fig.path,"euclidean","ward.D","silhouette for10sig")

plotdata <- t(mut.wt.trim.backup)
annCol <- read.table("F:/Project/UTUC_BLCA/InputData/UTUC_annotation.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
annCol <- annCol[which(annCol$`Exome-data` == "yes"),]
rownames(annCol) = paste0("UTUC_",rownames(annCol))
annCol <- annCol[colnames(plotdata),]
annCol <- data.frame("DMBcluster"=annCol$`Cluster DNA methylation no filter without normal`,row.names = rownames(annCol),stringsAsFactors = F)

annColors <- list()
annColors[["DMBcluster"]] <- c("C1"=red,"C2"=blue,"N/A"=lightgrey)

###########################################################################
########################   cluster for 2   ################################
#percentage weight matrix cluster
# cluster <-cutree(hcs,2)
# cluster <- ifelse(cluster==1,"C2","C1")
# 
# #percentage weight matrix consensus nmf
# mut.nmf <- nmf(nmf.input,2,seed = 4) #seed=4 is the best
# group <- predict(mut.nmf)
# 
# annCol$NMFbasis <- group[rownames(annCol)]
# annColors[["NMFbasis"]] <- c("1"="#9FA0FF","2"="#EE7BFF")
# annCol$MWBcluster <- cluster[rownames(annCol)]
# annColors[["MWBcluster"]] <- c("C1"="#0068B5","C2"="#F09300")
# 
# outFigFile <- file.path(fig.path,"heatmap_trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_13sig_cluster2.pdf")
# pdf(outFigFile, height=6,width=8)
# par(mfrow=c(2,2))
# basismap(mut.nmf, cexCol = 0.8, cexRow = 0.3)
# coefmap(mut.nmf,cexCol = 0.6, cexRow = 0.3)
# consensusmap(mut.nmf,cexCol = 0.5, cexRow = 1.0,annCol=data.frame("basis"=group))
# aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol, annColors=annColors, color=c("#EAF0FA","#6081C3","#3454A7"), revC=TRUE, fontsize=6, cexCol = 1, cexRow = 1)
# invisible(dev.off())
# 
# cluster <-cutree(hcs,2)
# cluster <- ifelse(cluster==1,"C2","C1")

###########################################################################
########################   cluster for 3   ################################
#UTUC annotation
UTUC.annotation <- read.table("F:/Project/UTUC_BLCA/InputData/UTUC_annotation.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
anno <- UTUC.annotation[,c(18:19,24)]
anno$sample <- paste0("UTUC_",rownames(anno))
#RBcluster 
RBcluster <- read.table("F:/Project/UTUC_BLCA/Results/cluster.assignment_by.20RNAseq.clustering_tumor.samples_UTUC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
RBcluster$sample <- paste0("UTUC_",rownames(RBcluster))
colnames(RBcluster) <- c("RBcluster","RBCfeature","sample")
RBcluster$RNAcluster <- RBcluster$RBCfeature

#percentage weight matrix consensus nmf
rownames(nmf.input) <- gsub("Signature","Sig",rownames(nmf.input))
mut.nmf <- nmf(nmf.input,3,seed = 8,method = "lee") #seed=8 is the best
group <- predict(mut.nmf)

index <- extractFeatures(mut.nmf,"max")
sig.order <- unlist(index)
sample.order <- names(group[order(group)])

annCol$MWBcluster <- group[rownames(annCol)]; annCol$MWBcluster <- paste0("C",annCol$MWBcluster)
annCol$sample <- rownames(annCol)
annCol <- merge(annCol,RBcluster,by="sample",all.x=T)
annCol <- merge(annCol,anno,by="sample",all.x=T)

annCol[is.na(annCol)] <- "N/A"
rownames(annCol) <- annCol$sample; annCol <- annCol[,-1]
annColors[["RNAcluster"]] <- c("Muscle_Enriched"=cherry,"Non-muscle_Enriched"=soil,"N/A"=lightgrey)
annColors[["MWBcluster"]] <- c("C1"="#F09300","C2"=lightgreen,"C3"="#0068B5")
annColors[["ZFP36_mut"]] <- c("Presence"=cherry,"Absence"=white)
annColors[["Del9p21.3"]] <- c("Presence"=cherry,"Absence"=white,"N/A"=lightgrey)
annColors[["FGFR3_mut"]] <- c("1"=cherry,"0"=white,"N/A"=lightgrey)

annCol <- annCol[,c("DMBcluster","RNAcluster","FGFR3_mut","ZFP36_mut","Del9p21.3","MWBcluster")]
group <- annCol$MWBcluster; names(group) <- rownames(annCol)

coef.matrix <- coef(mut.nmf)
basis.matrix <- basis(mut.nmf)

pdf(file = file.path(fig.path,"basismap.snp_0.2trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_11sig_cluster3.pdf"))
basismap(mut.nmf, cexCol = 0.8, cexRow = 0.16,annColors=list(c("1"="#F09300","2"=lightgreen,"3"="#0068B5")))
invisible(dev.off())

pdf(file = file.path(fig.path,"coefmap.snp_0.2trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_11sig_cluster3.pdf"))
coefmap(mut.nmf,cexCol = 0.8, cexRow = 0.1,annColors=list(c("1"="#F09300","2"=lightgreen,"3"="#0068B5")))
invisible(dev.off())

pdf(file = file.path(fig.path,"consensusmap.snp_0.2trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_11sig_cluster3.pdf"))
consensusmap(mut.nmf,cexCol = 0.8, cexRow = 0.8,annCol=data.frame("MWBcluster"=group[colnames(nmf.input)]),annColors=list(MWBcluster=c("C1"="#F09300","C2"=lightgreen,"C3"="#0068B5")))
invisible(dev.off())

plotdata <- t(mut.wt.trim.backup)
rownames(plotdata) <- gsub("Signature","Sig",rownames(plotdata))
pdf(file = file.path(fig.path,"heatmap.snp_0.2trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_11sig_cluster3_without_annobar.pdf"),width = 5,height = 3)
aheatmap(as.matrix(plotdata[sig.order,sample.order]), 
         Rowv=NA, 
         Colv=NA, 
         annCol=data.frame(MWBcluster=annCol[sample.order,"MWBcluster"],row.names = rownames(annCol)), annColors=annColors, 
         color=c("#EAF0FA","#6081C3","#3454A7"), 
         revC=TRUE, fontsize=6, cexCol = 1, cexRow = 1.2)
invisible(dev.off())

write.table(annCol,file.path(res.path,"Mut.sig.snp cluster3 results of DMBcluster MWBcluster RNAcluster.txt"),sep = "\t",row.names = T,col.names = NA)

# outFigFile <- file.path(fig.path,"heatmap_trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_12sig_cluster3.pdf")
# pdf(outFigFile, height=7,width=9)
# par(mfrow=c(2,2))
# 
# basismap(mut.nmf, cexCol = 0.8, cexRow = 0.7,annColors=list(c("1"="#F09300","2"=lightgreen,"3"="#0068B5")),fontsize=6,cexAnn = 0.8)
# 
# coefmap(mut.nmf,cexCol = 0.7, cexRow = 0.8,annColors=list(c("1"="#F09300","2"=lightgreen,"3"="#0068B5")),fontsize=6,cexAnn = 0.8)
# 
# consensusmap(mut.nmf,cexCol = 0.4, cexRow = 1,labCol=NA,labRow=NA,annCol=data.frame("MWBcluster"=group[colnames(nmf.input)]),annColors=list(MWBcluster=c("C1"="#F09300","C2"=lightgreen,"C3"="#0068B5")),fontsize=6,cexAnn = 0.8)
# 
# plotdata <- t(mut.wt.trim.backup)
# rownames(plotdata) <- gsub("Signature","Sig",rownames(plotdata))
# aheatmap(as.matrix(plotdata[sig.order,sample.order]), 
#          Rowv=NA, 
#          Colv=NA, 
#          annCol=annCol[sample.order,], annColors=annColors, 
#          color=c("#EAF0FA","#6081C3","#3454A7"), 
#          revC=TRUE, fontsize=6.5, cexCol = 0.45, cexRow = 0.6,cexAnn = 0.8)
# invisible(dev.off())

tmp <- basis(mut.nmf); colnames(tmp) <- c("basis1","basis2","basis3")
write.table(tmp,file.path(res.path,"basis matrix for 11sig of Mut.sig.cluster3.txt"),sep = "\t",row.names = T,col.names = NA)

tmp <- coef(mut.nmf); rownames(tmp) <- c("basis1","basis2","basis3")
write.table(tmp,file.path(res.path,"coef matrix for 11sig of Mut.sig.cluster3.txt"),sep = "\t",row.names = T,col.names = NA)

# ###################################################
# ### clustering Using significant mutation genes ###
# library(vegan)
# exome.mut.all <- read.table("F:/Project/UTUC_BLCA/Results/all exome mutation trimed.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
# exome.mut.cut10 <- rownames(exome.mut.all[rowSums(exome.mut.all)>=3,])
# exome.mut.p0.05 <- read.table("F:/Project/UTUC_BLCA/InputData/UTUC_exome_mutation_data_with_silence_curated_modified_longformat_forMutSigCV_simplified_addARID1A.sig_genes.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
# exome.mut.p0.05 <- rownames(exome.mut.p0.05[which(exome.mut.p0.05$p < 0.05),])
# #exome.mut.sel <- exome.mut.all[union(exome.mut.cut10,exome.mut.p0.05),]
# 
# tmp1 <- exome.mut.all[setdiff(union(exome.mut.cut10,exome.mut.p0.05),c("TTN","HYDIN")),]
# 
# #tmp1 <- exome.mut.all[exome.mut.p0.05,]
# #tmp1 <- exome.mut.all[setdiff(exome.mut.cut10,c("TTN","HYDIN")),]
# #tmp1 <- exome.mut.all[exome.mut.cut10,]
# 
# #tmp1 <- exome.mut.all[intersect(exome.mut.cut10,exome.mut.p0.05),]
# 
# tmp1 <- tmp1[rowSums(tmp1)>0,colSums(tmp1)>0]
# hcg <- hclust(vegdist(as.matrix(tmp1), method = "jaccard"), "ward.D")
# hcs <- hclust(vegdist(as.matrix(t(tmp1)), method = "jaccard"), "ward.D")
# 
# # hcg <- hclust(distanceMatrix(as.matrix(t(tmp1)), "euclidean"), "ward.D")
# # hcs <- hclust(distanceMatrix(as.matrix(tmp1), "euclidean"), "ward.D2")
# 
# tmp <- cutree(hcs,2); tmp <- ifelse(tmp == "1","C2","C1"); names(tmp) <- names(cutree(hcs,2))
# annCol$SMBcluster <- "N/A"
# annCol[paste0("UTUC_",names(tmp)),"SMBcluster"] <- tmp
# annColors[["SMBcluster"]] <- c("C1"=darkred,"C2"=darkblue)
# 
# outFigFile <- file.path(fig.path,"heatmap for 119 p0.05cut0.1 mutated genes removed TTN and HYDIN.pdf")
# pdf(outFigFile, height=6,width=6)
# aheatmap(as.matrix(tmp1), 
#          Rowv=dendsort(as.dendrogram(hcg)), 
#          Colv=dendsort(as.dendrogram(hcs)), 
#          annCol=annCol[paste0("UTUC_",colnames(tmp1)),], annColors=annColors, 
#          color=c(lightgrey,cherry), 
#          revC=TRUE, fontsize=6.5, cexCol = 1, cexRow = 0.3,cexAnn = 1)
# invisible(dev.off())
# 
# table(annCol$DMBcluster,annCol$SMBcluster)
# fisher.test(matrix(c(10,6,2,8),byrow = T,ncol = 2)) # 0.0511
#################################################
# test if signature is enriched in DMBcluster

# tmp <- annCol
# tmp$Sig13 <- ifelse(tmp$MWBcluster=="C2","Sig13","Others") 
# tmp$Sig16 <- ifelse(tmp$MWBcluster=="C1","Sig16","Others")
# tmp$Sig1  <- ifelse(tmp$MWBcluster=="C3","Sig1","Others")
# fisher.test(table(tmp$Sig13,tmp$DMBcluster)) #0.6159
# fisher.test(matrix(c(10,6,8,4),byrow = T,ncol = 2)) #1
# 
# fisher.test(table(tmp$Sig16,tmp$DMBcluster)) #0.4134
# fisher.test(matrix(c(13,5,5,5),byrow = T,ncol = 2)) #0.4119
# 
# fisher.test(table(tmp$Sig1,tmp$DMBcluster)) #0.3076
# fisher.test(matrix(c(13,9,5,1),byrow = T,ncol = 2)) #0.3746
# 
# fisher.test(table(tmp$Sig13,tmp$RBcluster)) #0.8133
# fisher.test(matrix(c(1,6,1,2,3,1),byrow = T,ncol = 3)) #0.7483
# 
# fisher.test(table(tmp$Sig16,tmp$RBcluster)) #0.7575
# fisher.test(matrix(c(2,7,1,1,2,1),byrow = T,ncol = 3)) #0.7483
# 
# fisher.test(table(tmp$Sig1,tmp$RBcluster))#0.3316
# fisher.test(matrix(c(3,5,2,0,4,0),byrow = T,ncol = 3)) #0.3646

#test
# mut.nmf <- nmf(nmf.input,2,seed = 4) #seed=4 is the best
# group <- predict(mut.nmf)
# 
# C1 <- paste0("UTUC_",c(12,13,14,15,17,18,19,23,29,32,34,37,39,03,06,07,08),"T") #17 sig16
# C2 <- paste0("UTUC_",c(16,25,28,02,30,31,33,35,36,38,04,05,09),"T") # sig13
# 
# intersect(names(group[group=="1"]),C1)
# intersect(names(group[group=="2"]),C2)
# 
# intersect(names(group[group=="1"]),C2)
# intersect(names(group[group=="2"]),C1)
# 
# basismap(mut.nmf, cexCol = 0.6, cexRow = 0.3)
# coefmap(mut.nmf,cexCol = 0.6, cexRow = 0.3)
# 
# sanger_wt.trim <- sanger_wt[colnames(nmf.input),paste0("Signature.",c(1,3,5,6,13,16,26,30))]
# plotdata <- as.data.frame(na.omit(standarize.fun(indata = t(sanger_wt.trim),halfwidth = 1)))
# hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
# hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "ward.D")
# annCol <- data.frame("Cluster"=paste0("C",group),row.names = names(group))
# hv <- aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol, annColors=NULL, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=6, cexCol = 1, cexRow = 1)
# 
# #########################
# 
# mut.nmf <- nmf(nmf.input,2,seed = 4) #seed=4 is the best
# group <- predict(mut.nmf)
# 
# sanger_wt.trim <- sanger_wt[colnames(nmf.input),paste0("Signature.",c(1,3,5,6,13,16,26,30))]
# plotdata <- as.data.frame(na.omit(standarize.fun(indata = t(sanger_wt.trim),halfwidth = 1)))
# hcg <- hclust(distanceMatrix(as.matrix(t(plotdata)), "pearson"), "ward.D")
# hcs <- hclust(distanceMatrix(as.matrix(plotdata), "pearson"), "ward.D")
# 
# cluster <-cutree(hcs,2)
# cluster <- ifelse(cluster==1,"C2","C1")
# annCol <- data.frame("basis"=paste0("C",group),"Cluster"=cluster,row.names = names(group))
# write.table(annCol,file.path(res.path,"cluster result for trimed deconstructSigs percentage 8sig.txt"),sep = "\t",row.names = T,col.names = NA)
# 
# annColors <- list()
# annColors[["Cluster"]] <- c("C1"="#0068B5","C2"="#F09300")
# outFigFile <- file.path(fig.path,"heatmap_trimed_UTUC_exome_mutation_signature_fromdeconstructSigs_8sig.pdf")
# pdf(outFigFile, height=6,width=8)
# par(mfrow=c(2,2))
# basismap(mut.nmf, cexCol = 0.8, cexRow = 0.3)
# coefmap(mut.nmf,cexCol = 0.6, cexRow = 0.3)
# consensusmap(mut.nmf,cexCol = 0.5, cexRow = 1.0,annCol=data.frame("basis"=group))
# hv <- aheatmap(as.matrix(plotdata), Rowv=dendsort(as.dendrogram(hcg)), Colv=dendsort(as.dendrogram(hcs)), annCol=annCol, annColors=annColors, color=c("#6699CC","white","#FF3C38"), revC=TRUE, fontsize=6, cexCol = 0.5, cexRow = 1)
# invisible(dev.off())

#save image
save.image(file.path(workdir,"MutSig.deconstructSigs.RData"))
