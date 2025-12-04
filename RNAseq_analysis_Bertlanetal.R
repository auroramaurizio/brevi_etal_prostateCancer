#Beltran et al. 2016 article.
#There are 47 samples, 13 CRPC-NE and 34 CRPC-Adeno. 
################ load libraries ################ 
suppressMessages(library("edgeR"))
suppressMessages(library(data.table))
suppressMessages(library("DESeq2"))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("RColorBrewer"))
#suppressMessages(library(enrichR))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library('VennDiagram'))
suppressMessages(library(venn))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library(patchwork))
suppressMessages(library(magick))
suppressMessages(library('philentropy'))
#suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library(assertr))
suppressMessages(library("remotes"))
suppressMessages(library(GeneOverlap))
suppressMessages(library('stringr'))
suppressMessages(library(ggrepel))

#setwd("/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/ric.bellone/BelloneM_1435_Neuroendocrine_TCGA/")

setwd("/Users/tascini.annasofia/OneDrive - Ospedale San Raffaele/ric.cosr/ric.bellone/BelloneM_1435_Neuroendocrine_TCGA/")

#-------- load counts ---------
load("../Neuroendocrine_cancer_analysis/Beltran2016_RNAseq_counts_toShare.RData")
counts = share_obj$counts
tpm = share_obj$tpm
# create txi object
txi <- setClass("txi", slots=list(counts="matrix"))
txi = list(counts = as.matrix(counts), countsFromAbundance = "scaledTPM")
class(txi) <- "txi"
str(txi$counts)
head(txi$counts)

metadata = share_obj$metadata
write.xlsx(metadata, "./metadata.xlsx")
#metadata <- metadata[,c("bind_ID", "Purity", "NEPC.vs..CRPC")]
#colnames(metadata)<- c("samplesID", "Purity", "condition")
row.names(metadata) <- metadata$samplesID
head(metadata)


annotation_SI = read.xlsx("Neuroendocrine_cancer_analysis/metadata_modified.xlsx", sheet = "annotation")
annotation_SI_RNAsample = annotation_SI[annotation_SI$Tumor.Sample.ID %in% metadata$CBIO_Common_ID,]
row.names(metadata) <- metadata$CBIO_Common_ID
row.names(annotation_SI_RNAsample) <- annotation_SI_RNAsample$Tumor.Sample.ID
metadata$Tumor.Sample.ID = metadata$CBIO_Common_ID
full_metadata = base::merge(metadata, annotation_SI_RNAsample, by = "Tumor.Sample.ID")
colnames(full_metadata)
write.xlsx(full_metadata, "./combined_metadata.xlsx")

order_meta =c("Tumor.Sample.ID" ,           "bind_ID"    ,                "Patient.ID.x"     ,         
"Ethnicity"   ,               "Treatment"    ,              "Pathology.Classification.x",
 "Purity"       ,              "Ploidy"        ,             "Genomic_Burden.x"          ,
"CBIO_Common_ID" ,            "NEPC.vs..CRPC"   ,           "Integrated NEPC Score"     ,
 "AR Signaling"   ,            "Patient.ID.y"    ,           "Pathology.Classification.y",
 "Body.Site"       ,           "Purity.(CLONET)"  ,          "Ploidy.(CLONET)"           ,
 "Genomic_Burden.y" ,          "#.non-silent.SNV"  ,         "Control.Site"              ,
"Exome.Capture.Kit" )


full_metadata$NEPC.vs..CRPC = as.factor(full_metadata$NEPC.vs..CRPC)
p <- ggplot(full_metadata, aes(x = NEPC.vs..CRPC, fill = Body.Site)) + geom_bar(position = 'dodge') +
  theme_minimal() + 
  geom_label(aes(label = as.character(c(NCRPC, NNEPC))))
p

NCRPC = as.numeric(table(full_metadata[full_metadata$NEPC.vs..CRPC %in% "CRPC",]$Body.Site))
NNEPC = as.numeric(table(full_metadata[full_metadata$NEPC.vs..CRPC %in% "NEPC",]$Body.Site))
NNEPC
table(full_metadata$Body.Site)

library(dplyr)
full_metadata %>%
  group_by(NEPC.vs..CRPC, Body.Site) %>%
  tally() %>%
  ggplot(., aes(x = NEPC.vs..CRPC, y = n, fill=Body.Site)) +
  geom_bar(stat = "identity", position ="dodge") +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9)) +
  theme_minimal() + xlab("Tumor Type") + ylab( "#N")
ggsave(filename = "barplot_BodySite.png")


full_metadata %>%
  group_by(NEPC.vs..CRPC, Pathology.Classification.x) %>%
  tally() %>%
  ggplot(., aes(x = NEPC.vs..CRPC, y = n, fill=Pathology.Classification.x)) +
  geom_bar(stat = "identity", position ="dodge") +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9)) +
  theme_minimal() + xlab("Tumor Type") + ylab( "#N")

ggsave(filename = "barplot_Pathology.Classification.png")

table(full_metadata[full_metadata$NEPC.vs..CRPC %in% "CRPC",]$Pathology.Classification.x)

full_metadata %>%
  group_by(NEPC.vs..CRPC) %>%
  tally() %>%
  ggplot(., aes(x = NEPC.vs..CRPC, y = n)) +
  geom_bar(stat = "identity", position ="dodge", fill = c("#CC0000", "#3366CC")) +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9)) +
  theme_minimal() + xlab("Tumor Type") + ylab( "#N")

ggsave(filename = "barplot_TT.png")


full_metadata$condition = full_metadata$NEPC.vs..CRPC

library("DESeq2")

# ---- Deseq2 --------
row.names(full_metadata) <- full_metadata$bind_ID
full_metadata_ro = full_metadata[colnames(txi$counts),]
ddsTxi <- DESeqDataSetFromTximport(txi$counts,
                                   colData = full_metadata_ro,
                                   design = ~condition)

Nreplica = min(table(full_metadata$condition))

filter <- rowSums(cpm(counts(ddsTxi)) >= 1) >= Nreplica
table(filter)
ddsFiltered <- ddsTxi[filter,]

ddsFiltered[['condition']] <- relevel(
  ddsFiltered[['condition']] , ref = 'CRPC')

dga <- DESeq(
  object = ddsFiltered,
  test = "Wald",
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf)

plotDispEsts(dga)


# ---------- pca ---------
vsd <- vst(dga, blind=FALSE)
rld <- rlog(dga, blind=FALSE)
#BiocManager::install("vsn")
library("vsn") 
meanSdPlot(assay(vsd))

select <- order(rowMeans(counts(dga,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dga)[,c("condition", "Purity")])
ntd <- normTransform(dga)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE,  annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$Body.Site, vsd$Patient.ID.x ,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cellwidth = 15, cellheight = 15,
         filename = "SampleDist.pdf")
dev.off()
plotPCA(vsd, intgroup=c("condition"))
vsd$Pathology.Classification.x
colnames(colData(dga))
pcaData <- plotPCA(vsd, intgroup=c("condition", "Patient.ID.x",
                                   "Integrated NEPC Score",
                                   "Body.Site", "Pathology.Classification.x",
                                   "Purity", "Genomic_Burden.x"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Body.Site, shape=condition, label = Patient.ID.x)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave("PCA_BodySite.png")

ggplot(pcaData, aes(PC1, PC2, color=Integrated.NEPC.Score, shape=condition, label = Patient.ID.x)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave("PCA_Integrated.NEPC.Score.png")

ggplot(pcaData, aes(PC1, PC2, color=Body.Site, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave("PCA_BodySite.png")

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()

################ DGE - comparisons ########################
resultsNames(dga)
contrasts = resultsNames(dga)[- which(resultsNames(dga) %in% 'Intercept')]
alpha = 0.05

dgeResults <- list()
for (contrast in contrasts) {
  dgeResults[[contrast]] <- results(
    dga,
    name                 = contrast,
    cooksCutoff          = Inf,
    independentFiltering = TRUE,
    alpha                = alpha,
    pAdjustMethod        = "BH")
  print(summary(dgeResults[[contrast]]))
  # sorting gene list according to significance
  dgeResults[[contrast]] <- dgeResults[[contrast]][order(dgeResults[[contrast]]$pvalue, decreasing = F),]
}
contrast

################ DGE - Save results ######################## 
# Save results
f="dgeResults"
dir.create(f, showWarnings=TRUE, recursive=TRUE)
geneidColname <- 'Geneid'
lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))


dgeResults_table = list()    
dgeResults_table = lapply(
  names(dgeResults),
  function(x) 
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = paste(f,'/DGE_results.xlsx', sep=''), 
           row.names = F,
           asTable = T, 
           sheetName = str_sub(names(dgeResults),1,31)) 



################ DGE - MAplot and Vulcano Plots ######################## 
f="dgeResults"
n.label = 10
FDR = T
pvalue = 0.01
for (condition in names(dgeResults)) {
  results = as.data.frame(dgeResults[[condition]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }        
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }        
  
  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))
  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(1,.5); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)
  
  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]
  
  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", condition)) +
    scale_color_manual(values = mycolors) +
    #scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +  
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")  
  
  #print(MAplot)
  #pdf(paste(f,'/','MAplot_',condition,'.pdf',sep=''),width=8, height=6.5)
  #plot(MAplot)
  #dev.off()
  
  # Vulcano plot 
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    ggtitle(paste("Vulcano Plot,", condition)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[c(up_label, down_label),], 
                     aes(x = log2FoldChange, y = -log10(padj), color = DE), 
                     label = row.names(results[c(up_label, down_label),]), size = 2, max.overlaps = 100) + #,
                     #box.padding = unit(0.35, "lines"), point.padding = unit(0.6, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value")   
  
  print(VP)
  pdf(paste(f,'/','VulcanoPlot_',condition,'.pdf',sep=''),width=8, height=8)
  plot(VP)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;   
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  
  pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  dev.off()
  
  A1 <- image_read_pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''), density = 140)
  image_write(A1, path = paste(f,'/','MA_VP_',condition,'.tiff',sep=''), format = "tiff")
  
}

####### FDR ########
# save FDR genes in a list
fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)       

fdrDW = list()
fdrDW = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)

####### SEQC ########
# save SEQC genes in a list
seqcUP = list()

pvalue = 0.01
seqcUP = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange > 1])
names(seqcUP)= names(dgeResults)               

seqcDW = list()
seqcDW = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange < -1])
names(seqcDW)= names(dgeResults)
names(dgeResults)

print('FDRup');lengths(fdrUP); print('FDRdw');lengths(fdrDW)
print('SEQCup');lengths(seqcUP); print('SEQCdw');lengths(seqcDW)

print('FDR_tot');lengths(fdrUP) + lengths(fdrDW)
print('SEQC_tot'); lengths(seqcUP) + lengths(seqcDW)

# -------- heatmap DGE fdr ----------
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors_hm  = crp(255)
heatmap_dir_as = "./Paper/"
dir.create(heatmap_dir_as)
cpm = cpm(counts(dga), log = T)

annotation_column <- as.data.frame(full_metadata_ro[,c("condition")])
colnames(annotation_column) <- c("condition")
annotation_column$pz <- row.names(full_metadata_ro)
annotation_column = annotation_column[order(annotation_column$condition),]

annotation_column <- as.data.frame(annotation_column$condition, row.names = annotation_column$pz)
colnames(annotation_column) <- c("condition")
#annotation_column$Body.Site = factor(annotation_column$Body.Site)
mycolors_c <- c("#9900FF", "#339900");     
names(mycolors_c) = levels(annotation_column$condition)
#mycolors_bs <- rainbow(10);
#names(mycolors_bs) = levels(annotation_column$Body.Site)

ann_colors = list(
  condition = mycolors_c #,
  #Body.Site = mycolors_bs
)


myb = seq(-3.5,3.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
HPv <- pheatmap::pheatmap(cpm[c(fdrUP$condition_NEPC_vs_CRPC, 
                                fdrDW$condition_NEPC_vs_CRPC),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = T, 
                          cluster_cols = T, 
                          #cutree_cols  = 5,
                          #cutree_row  = 2,
                          show_rownames = F,
                          show_colnames = F,
                          border_color = NA,
                          cellwidth=10, cellheight=0.1,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_allFDR.pdf',sep=''))

dev.off()


file = "./enrichR/condition_NEPC_vs_CRPC_fdr_both_.xlsx"
df = read.xlsx(file, sheet = "GO_Biological_Process_2018")
g = unlist(str_split(df[grep(pattern = "GO:0030198",
        x = df$Term),"Genes"], pattern = ";"))

HPv <- pheatmap::pheatmap(cpm[g,
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = T, 
                          cluster_cols = F, 
                          #cutree_cols  = 5,
                          #cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = F,
                          border_color = NA,
                          cellwidth=10, cellheight=5,
                          fontsize = 10, fontsize_row = 5, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_ExMO.pdf',sep=''))


#### ------
table(metadata$condition)

g = c('ITGA2', 'YAP1', 'ILK', 'MST1', 'TAZ')

df = metadata
my_comparisons <- list(c("CRPC", "NEPC") )
frpkm <- cpm(round(counts)); ylabel = "CPM";dir= 'gene_CPM/'
#frpkm <- tpm; ylabel = "TPM"; dir= 'gene_TPM/'
dir.create(dir, recursive = TRUE)
for (gene in g) {
  df$rpkm = as.numeric(frpkm[gene,])
  row.names(df) <- df$samplesID
  df$ID = as.character(1:length(df$samplesID))   
  p <- ggboxplot(df, x = "condition", y = "rpkm", 
                 #label = str_sub(metadata$SampleID, 1, 2), repel = T, label.rectangle = T,
                 fill = "condition", color = "condition", width = 0.8,
                 add = c('jitter'), add.params = list(size =3, shape = 18),
                 alpha = 0.5, palette = viridis(2),
                 xlab = "tumor type", ylab = ylabel,
                 short.panel.labs = FALSE, title = gene, 
                 font = list(size = 16, face = 'bold', color = "black")) + 
    theme(text = element_text(size=18))

  assign(paste('p',gene,sep='_'),
         p) 
  pdf(paste(dir,'boxplot_gene_',gene,'.pdf',sep=''), width = 6, height = 5)
  print(p)
  dev.off()
  pdf(paste(dir,'boxplot_gene_',gene,'_pvalue.pdf',sep=''), width = 6, height = 5)
  print(p + stat_compare_means(label = "p.value", 
                               comparisons = my_comparisons))
  dev.off()
}

###### enrichment #######

# Prepare lists of genes to run GSEA
# ------- GSEA ----------
outdir = 'GSEA/'
dir.create(outdir)
for (c in names(dgeResults)) {
  l_ranked = data.frame(GeneID= row.names(dgeResults[[c]]), LogFC = dgeResults[[c]][,'log2FoldChange'])
  l_ranked = l_ranked[order(l_ranked$LogFC, decreasing = T),]
  write.table(l_ranked, file = paste(outdir, c,'_ranked_list.rnk', sep =''), 
              quote = F, row.names= F, col.names = F, sep ='\t')
}

# ------- enrichR ----------
databases <- listEnrichrDbs()
enrichf = 'enrichR/' 
dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
# enrichment Parameters
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016")
# alpha used in DGE
padj.cutoff = alpha; 

# Perform Enrichment
enrichr.list <- list()
for (i in 1:length(dgeResults)){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUP[[i]]
  down.genes <- fdrDW[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste(enrichf,'/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  write.table(down.genes, paste(enrichf,'/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''), 
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste(enrichf,'/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  
  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)

keggup = enrichr.list[["condition_NEPC_vs_CRPC"]][["fdr_up"]][["KEGG_2016"]]
keggdw = enrichr.list[["condition_NEPC_vs_CRPC"]][["fdr_down"]][["KEGG_2016"]]
keggup_f = keggup[keggup$Adjusted.P.value < 0.05,]
keggup_f$condition = "up"
keggup_f$cn = 1
keggdw_f = keggdw[keggdw$Adjusted.P.value < 0.05,]
keggdw_f$condition = "down"
keggdw_f$cn = -1
kegg_f = merge(keggup_f, keggdw_f, all.y = TRUE, all.x = TRUE)
kegg_f$Enrichment = -log10(kegg_f$P.value)*kegg_f$cn
kegg_f = kegg_f[order(kegg_f$Enrichment, decreasing = T),]
kegg_f$Pathway.num = dim(kegg_f)[1]:1 
kegg_f$Pathway.num = as.factor(kegg_f$Pathway.num)

fN <- function(x){as.numeric(unlist(str_split(string = x, pattern = "/"))[1])}
fp <- function(x){
  n = as.numeric(unlist(str_split(string = x, pattern = "/"))[1])
  t = as.numeric(unlist(str_split(string = x, pattern = "/"))[2])
  return(n/t)
}
kegg_f$GeneCount = sapply(kegg_f$Overlap,fN)
kegg_f$GenePercent = sapply(kegg_f$Overlap,fp)

patplot = ggplot(data=kegg_f, aes(x=Pathway.num, y=Enrichment, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge(), color = 'black') +
  scale_fill_manual(values=c('#3366CC','#CC0000'), name = "Modulation") +
  coord_flip() + 
  scale_x_discrete(breaks=kegg_f$Pathway.num, 
                   labels=str_remove(string = as.character(kegg_f$Term), 
                                     pattern = "Homo sapiens"),
                   position = "bottom") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "italic", color = 'black', size=8, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(angle = 0, face = "italic", color = 'black', size=12),#, family = 'mono'),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        legend.position = 'top',
        legend.text = element_text(face = "bold", color = "black", size = 12),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))+
  labs(x = "Pathways", y = "Enrichment") 
ggsave(filename = paste(heatmap_dir_as, "keggPathway.pdf", sep = ''),
       patplot, dev = "pdf", width = 20, height = 12)


df = kegg_f[kegg_f$condition == "down",]
df = df[order(df$P.value, decreasing = F), ]
df$Pathway.num = dim(df)[1]:1 
df$Pathway.num = as.factor(df$Pathway.num)
GCplotPV = ggplot(data=df, aes(x=Pathway.num, y=GeneCount, fill=-Enrichment)) +
  geom_bar(stat="identity", position=position_dodge(), color = 'black') +
  scale_fill_gradient(low ="white", high = "red") +
  #scale_fill_manual(values=c('#3366CC','#CC0000'), name = "Modulation") +
  coord_flip() + 
  scale_x_discrete(breaks=df$Pathway.num, 
                   labels=str_remove(string = as.character(df$Term), 
                                     pattern = "Homo sapiens"),
                   position = "bottom") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "italic", color = 'black', size=8, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(angle = 0, face = "italic", color = 'black', size=12),#, family = 'mono'),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        legend.position = 'right',
        legend.text = element_text(face = "bold", color = "black", size = 12),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))+
  labs(x = "Pathways", y = "Enrichment") 
ggsave(filename = paste(heatmap_dir_as, "keggPathway_GC_pvalue.pdf", sep = ''),
       GCplotPV, dev = "pdf", width = 20, height = 15)

# Write excels files
for (i in 1:length(dgeResults)){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path('enrichR', 
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
  }
}




 #----------- Jaccard distances -----------
#' Jaccard distances
dir = paste('enrichR_padj_',padj.cutoff,'_LgFC_',LgFC.cutoff,'/JaccardPlots/', 
            sep = '')
dir.create(dir)
p.value.thr = 0.005

breaksList = seq(0, 1, by = 0.001)

cutree_rows_values = c(5,10,5)
lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))

lf
for (file in lf) {    
  enrichR.file = paste("./enrichr/",
                       file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  Pathways.Table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways 
    colnames(M) = gene.all 
    
    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1 
        }            
      }    
    }
    
    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }
      
      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      
      w=5; h=5; fs = 4; cutree_rows_N = 5
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc, 
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               cellwidth = w, cellheight = h,
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               #main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 12,
               fontsize_row = fs,
               filename = paste(dir,file,'_JaccardDist.pdf', sep=''))
    }
  }
}

dev.off()

# 
############## Plot top N pathways in a plot ################
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))
N=25
plotdir = paste(enrichf,'','/TOP_',N,'_pathways/',sep='')
dir.create(plotdir, recursive = T)

path.sigs = list()

lf = list.files(paste(enrichf,'','/', sep=''), pattern=glob2rx("*.xlsx"))
s 
for (file in lf) { 
  path.sign_dataframe = data.frame()
  enrichR.file = paste(enrichf,'','/',file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  enrichR.table = data.frame()
  #s[1:length(s)]
  for (dat in c("KEGG_2016")) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    enrichR.table = rbind(enrichR.table, Table)
  }
  p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.1,])
  if (length(p)>0) {
    path.sign_dataframe.entry = data.frame(Pathway = p,
                                          gene.ratio = sapply(p, fx),
                                          p.value = enrichR.table[p,]$P.value,
                                          p.value.adj = enrichR.table[p,]$Adjusted.P.value,
                                          genes = enrichR.table[p,]$Genes)
    path.sign_dataframe <- rbind(path.sign_dataframe, path.sign_dataframe.entry) 
    
    path.sigs[[file]] = path.sign_dataframe
    
    #up
    
    #path.sign_dataframe = path.sigs[[file]]
    path.sign_dataframe = path.sign_dataframe[order(path.sign_dataframe$p.value.adj),]
    path.sign_dataframe$Pathway.num = dim(path.sign_dataframe)[1]:1
    path.sign_dataframe$Pathway.num = as.factor(path.sign_dataframe$Pathway.num )
    
    if (str_detect(file, pattern = "up")) {
      up_down = '#CC0033'
    } else if (str_detect(file, pattern = "down")) {
      up_down = '#0066CC'
    } else {
      up_down = '#FF99FF'
    }
    
    pd = ggplot(path.sign_dataframe[1:N,], aes(Pathway.num,-log10(p.value.adj))) + 
      geom_point(aes(size = gene.ratio), color = up_down) +
      scale_size_continuous(range = c(2,8), name = "Gene ratio") +
      #scale_colour_viridis_c(begin = 0, end = 40) + 
      #ylim(c(10,40)) +
      coord_flip() +
      scale_x_discrete(breaks=path.sign_dataframe$Pathway.num, 
                       labels=stringr::str_wrap(path.sign_dataframe$Pathway,width = 85),
                       position = "top") +
      theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
            plot.subtitle = element_text(color="black", size=16, face="italic"),
            axis.text.x = element_text(angle = 90, color = 'black', size=14, hjust =1, family = "mono"), 
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(angle = 0, color = 'black', size=14, family = "mono"),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            legend.text = element_text(color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14),
            panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
      labs(title = paste("Top", N, 'pathways -', str_sub(file,1,-7)), 
           x = "Pathways", 
           y = "-Log10(p.adj.value)") 
    print(pd)
    pdf(paste(plotdir, str_sub(file,1,-7),'_TOP',N,'.pdf',sep=''), 
        width = 14, height = 8)
    print(pd)
    dev.off()
  }
  
}    

path.sigs$condition_NEPC_vs_CRPC_fdr_up_.xlsx
path.sigs_m = path.sigs[c(2,3)]
path.sigs$condition_NEPC_vs_CRPC_fdr_up_.xlsx$condition = 'up'
path.sigs$condition_NEPC_vs_CRPC_fdr_down_.xlsx$condition = 'down'
path.sign_dataframe = rbind(path.sigs$condition_NEPC_vs_CRPC_fdr_up_.xlsx,
                            path.sigs$condition_NEPC_vs_CRPC_fdr_down_.xlsx)
path.sign_dataframe = path.sign_dataframe[path.sign_dataframe$p.value.adj<0.005,]
dim(path.sign_dataframe)
path.sign_dataframe$Enrichment = -log10(path.sign_dataframe$p.value)
path.sign_dataframe[path.sign_dataframe$condition == 'down',]$Enrichment = path.sign_dataframe[path.sign_dataframe$condition == 'down',]$Enrichment * -1
path.sign_dataframe = path.sign_dataframe[order(path.sign_dataframe$Enrichment, decreasing = TRUE),]
path.sign_dataframe$Pathway.num = dim(path.sign_dataframe)[1]:1
path.sign_dataframe$Pathway.num = as.factor(path.sign_dataframe$Pathway.num)


patplot = ggplot(data=path.sign_dataframe, aes(x=Pathway.num, y=Enrichment, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge(), color = 'black') +
  scale_fill_manual(values=c('#3366CC','#CC0000'), name = "Modulation") +
  coord_flip() + 
  scale_x_discrete(breaks=path.sign_dataframe$Pathway.num, 
                   labels=str_remove(string = as.character(path.sign_dataframe$Pathway), 
                                     pattern = "Homo sapiens"),
                   position = "bottom") +
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        plot.subtitle = element_text(color="black", size=16, face="italic"),
        axis.text.x = element_text(angle = 0, face = "italic", color = 'black', size=14, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "italic", color = 'black', size=14),#, family = 'mono'),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.position = 'top',
        legend.text = element_text(face = "bold", color = "black", size = 12),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))+
  labs(x = "Pathways", y = "Enrichment") 
patplot

pdf(paste(enrichf,'/kegg_pathway_padj_0.005.pdf',sep=''), width = 18, height = 10)
print(patplot)
dev.off()
A1 <- image_read_pdf(paste(enrichf,'/kegg_pathway_padj_0.005.pdf',sep=''), density = 140)
image_write(A1, path = paste(enrichf,'/kegg_pathway_padj_0.005.tiff',sep=''), format = "tiff")

lfp = list.files(paste("./enrichR/JaccardPlots/", sep=''), pattern=glob2rx("*.pdf"))
lfp

for (file in lfp) {    
  plot.file = paste("./enrichr/JaccardPlots/",
                       file,sep='')
  A1 <- image_read_pdf(plot.file, density = 140)
  image_write(A1, path = paste(plot.file,'.tiff',sep=''), format = "tiff")
  
}

# ---------- heatmap Hippo --------- 
gene_Hippo = sort(unlist(str_split(path.sign_dataframe["Hippo signaling pathway Homo sapiens hsa04390",]$genes, 
                       pattern = ";")))

gene_PI3KAkt = sort(unlist(str_split(path.sign_dataframe["PI3K-Akt signaling pathway Homo sapiens hsa04151",]$genes, 
                              pattern = ";")))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors_hm  = crp(255)
heatmap_dir_as = "./CustomHM/"
dir.create(heatmap_dir_as)
cpm = cpm(counts(dga), log = T)

annotation_column <- as.data.frame(metadata[,2:(dim(metadata)[2])])

mycolors_c <- c("#9900FF", "#339900");     names(mycolors_c) = levels(annotation_column$condition)

ann_colors = list(
  condition = mycolors_c,
  Purity = scale_color_gradient(low = "yellow", high = "darkblue")
)

row.names(annotation_column) <- metadata[,1]

HPv <- pheatmap::pheatmap(cpm[gene_Hippo,],
                          scale = 'row',
                          annotation_col = annotation_column,
                          #annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 5,
                          #cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          cellwidth=10, cellheight=10,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_Hippo_pathway.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[gene_PI3KAkt,],
                          scale = 'row',
                          annotation_col = annotation_column,
                          #annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 5,
                          #cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          cellwidth=10, cellheight=10,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_PI3KAkt_pathway.pdf',sep=''))
dev.off()

lfp = list.files(paste("./CustomHM/", sep=''), pattern=glob2rx("*.pdf"))
lfp

for (file in lfp) {    
  plot.file = paste("./CustomHM/",
                    file,sep='')
  A1 <- image_read_pdf(plot.file, density = 140)
  image_write(A1, path = paste(plot.file,'.tiff',sep=''), format = "tiff")
  
}

a = as.list(read_tsv(file = "GSEA/costum_hsa.gmx"))
for (col in names(a)) {
  gene = as.character(na.omit(a[[col]][2:length(a[[col]])]))
  g_exp = intersect(gene,row.names(dgeResults[[1]]))
  write.xlsx(x = as.data.frame(dgeResults[[1]][g_exp,]), 
             file = paste(col, "DGEstatistic_NEPC_vs_CRPC.xlsx", sep="_"),
             asTable = T, row.names = TRUE)
  HPv <- pheatmap::pheatmap(cpm[g_exp,],
                            scale = 'row',
                            annotation_col = annotation_column,
                            #annotation_colors = ann_colors, 
                            cluster_rows = T, 
                            cluster_cols = T, 
                            cutree_cols  = 5,
                            cutree_row  = 4,
                            show_rownames = T,
                            show_colnames = T,
                            cellwidth=10, cellheight=10,
                            #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                            display_numbers = F,
                            col=colors_hm,
                            filename = paste(heatmap_dir_as,col,'_allGenes.pdf',sep=''))
}

table(metadata$condition)
dgeResults[[1]][gene,]

# ---------- Purity -------
# filter out NA sample
metadata_f = metadata[!is.na(metadata$Purity),]
txi_f = txi
txi_f$counts = txi$counts[, metadata_f$samplesID]
ddsTxi <- DESeqDataSetFromTximport(txi_f,
                                   colData = metadata_f,
                                   design = ~condition+Purity)

Nreplica = min(table(metadata_f$condition))
filter <- rowSums(cpm(counts(ddsTxi)) >= 1) >= Nreplica
table(filter)
ddsFiltered <- ddsTxi[filter,]

ddsFiltered[['condition']] <- relevel(
  ddsFiltered[['condition']] , ref = 'CRPC')

dga <- DESeq(
  object = ddsFiltered,
  test = "Wald",
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf)

plotDispEsts(dga)

################ DGE - comparisons ########################
resultsNames(dga)
contrasts = resultsNames(dga)[- which(resultsNames(dga) %in% 'Intercept')]
alpha = 0.05

dgeResults <- list()
for (contrast in contrasts) {
  dgeResults[[contrast]] <- results(
    dga,
    name                 = contrast,
    cooksCutoff          = Inf,
    independentFiltering = TRUE,
    alpha                = alpha,
    pAdjustMethod        = "BH")
  print(summary(dgeResults[[contrast]]))
  # sorting gene list according to significance
  dgeResults[[contrast]] <- dgeResults[[contrast]][order(dgeResults[[contrast]]$pvalue, decreasing = F),]
}

# ---------- pca ---------
vsd <- vst(dga, blind=FALSE)
rld <- rlog(dga, blind=FALSE)
#BiocManager::install("vsn")
library("vsn") 
meanSdPlot(assay(vsd))

select <- order(rowMeans(counts(dga,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dga)[,c("condition", "Purity")])
ntd <- normTransform(dga)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE,  annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$samplesID, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))

pcaData <- plotPCA(vsd, intgroup=c("condition", "Purity"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Purity, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()

pcaData <- plotPCA(vsd, intgroup=c("condition", "Purity"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()


####### FDR ########
# save FDR genes in a list
fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)       

fdrDW = list()
fdrDW = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)
lengths(fdrDW)
lengths(fdrUP)
###### enrichment #######

# Prepare lists of genes to run GSEA

# ------- enrichR ----------
databases <- listEnrichrDbs()
dir.create('enrichR_purity/', showWarnings=FALSE, recursive=TRUE)
# enrichment Parameters
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016")
# alpha used in DGE
padj.cutoff = alpha; 

# Perform Enrichment
enrichr.list <- list()
#for (i in 1:length(dgeResults)){
for (i in 1){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUP[[i]]
  down.genes <- fdrDW[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste('./enrichR_purity/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  write.table(down.genes, paste('./enrichR_purity/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''), 
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste('./enrichR_purity/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  
  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)[1]

# Write excels files
#for (i in 1:length(dgeResults)){
for (i in 1){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path('enrichR_purity', 
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
  }
}



#----------- Jaccard distances -----------
#' Jaccard distances
dir = "./enrichR_purity/Jaccard_dist/"
dir.create(dir)
p.value.thr = 0.005

breaksList = seq(0, 1, by = 0.001)

cutree_rows_values = c(5,10,5)
lf = list.files(paste("./enrichR_purity/", sep=''), pattern=glob2rx("*.xlsx"))
 
lf
for (file in lf) {    
  enrichR.file = paste("./enrichr_purity/",
                       file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  Pathways.Table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways 
    colnames(M) = gene.all 
    
    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1 
        }            
      }    
    }
    
    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }
      
      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      
      w=5; h=5; fs = 4; cutree_rows_N = 5
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc, 
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               cellwidth = w, cellheight = h,
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               #main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 12,
               fontsize_row = fs,
               filename = paste(dir,file,'_JaccardDist.pdf', sep=''))
    }
  }
}

dev.off()

# ------
row.names(full_metadata) <- full_metadata$bind_ID
full_metadata_ro = full_metadata[colnames(txi$counts),]

full_metadata_ro$Body.Site = str_replace_all(string = full_metadata_ro$Body.Site, pattern = "-", replacement = "")
full_metadata_ro$Body.Site = str_replace_all(string = full_metadata_ro$Body.Site, pattern = " ", replacement = ".")

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = full_metadata_ro,
                                   design = ~condition+Body.Site)

Nreplica = min(table(full_metadata_ro$condition))

filter <- rowSums(cpm(counts(ddsTxi)) >= 1) >= Nreplica
table(filter)
ddsFiltered <- ddsTxi[filter,]

ddsFiltered[['condition']] <- relevel(
  ddsFiltered[['condition']] , ref = 'CRPC')

dga <- DESeq(
  object = ddsFiltered,
  test = "Wald",
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf)

plotDispEsts(dga)

resultsNames(dga)
contrasts = resultsNames(dga)[- which(resultsNames(dga) %in% 'Intercept')]
alpha = 0.05

dgeResults <- list()
for (contrast in contrasts) {
  print(contrast)
  dgeResults[[contrast]] <- results(
    dga,
    name                 = contrast,
    cooksCutoff          = Inf,
    independentFiltering = TRUE,
    alpha                = alpha,
    pAdjustMethod        = "BH")
  print(summary(dgeResults[[contrast]]))
  # sorting gene list according to significance
  dgeResults[[contrast]] <- dgeResults[[contrast]][order(dgeResults[[contrast]]$pvalue, decreasing = F),]
}

contrast



# ------ select lymph node samples -------
folder = "LymphNode/"
dir.create(folder)
row.names(full_metadata) <- full_metadata$bind_ID
full_metadata_ro = full_metadata[colnames(txi$counts),]
# samples selection
full_metadata_sel = full_metadata_ro[full_metadata_ro$Body.Site == "Lymph node",]
txi_sel = txi; txi_sel$counts =  txi_sel$counts[,row.names(full_metadata_sel)]
head(txi_sel)
head(full_metadata_sel)
dim(txi_sel$counts)
# ---- Deseq2 --------
ddsTxi <- DESeqDataSetFromTximport(txi_sel,
                                   colData = full_metadata_sel,
                                   design = ~condition)

Nreplica = min(table(full_metadata_sel$condition))

filter <- rowSums(cpm(counts(ddsTxi)) >= 1) >= Nreplica
table(filter)
ddsFiltered <- ddsTxi[filter,]

ddsFiltered[['condition']] <- relevel(
  ddsFiltered[['condition']] , ref = 'CRPC')

dga <- DESeq(
  object = ddsFiltered,
  test = "Wald",
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf)

plotDispEsts(dga)


# ---------- pca ---------
vsd <- vst(dga, blind=FALSE)

select <- order(rowMeans(counts(dga,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dga)[,c("condition", "Purity")])
ntd <- normTransform(dga)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE,  annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition,vsd$Patient.ID.x ,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cellwidth = 15, cellheight = 15,
         filename = paste(folder,"SampleDist.pdf", sep=''))
dev.off()
plotPCA(vsd, intgroup=c("condition"))
vsd$Pathology.Classification.x
colnames(colData(dga))
pcaData <- plotPCA(vsd, intgroup=c("condition", "Patient.ID.x",
                                   "Integrated NEPC Score",
                                   "Treatment",
                                   "Body.Site", "Pathology.Classification.x",
                                   "Purity", "Genomic_Burden.y"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition, label = Patient.ID.x)) +
  geom_point(size=3) +
  scale_color_manual(values = c("#CC0000", "#3366CC")) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave(paste(folder,"PCA.png", sep =''))

ggplot(pcaData, aes(PC1, PC2, color=Purity, shape=condition, label = Patient.ID.x)) +
  geom_point(size=3) +
  #scale_color_distiller(palette = "YlOrBr") +
  scale_color_viridis_c(option = "plasma", begin = 0, end = 0.9) + 
  #scale_color_manual(values = c("#CC0000", "#3366CC")) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave(paste(folder,"PCA_purity.png", sep =''))


ggplot(pcaData, aes(PC1, PC2, color=Integrated.NEPC.Score, shape=condition, label = Patient.ID.x)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_linedraw()
ggsave(paste(folder,"PCA_Integrated.NEPC.Score.png", sep =''))


################ DGE - comparisons ########################
resultsNames(dga)
contrasts = resultsNames(dga)[- which(resultsNames(dga) %in% 'Intercept')]
alpha = 0.05

dgeResults <- list()
for (contrast in contrasts) {
  dgeResults[[contrast]] <- results(
    dga,
    name                 = contrast,
    cooksCutoff          = Inf,
    independentFiltering = TRUE,
    alpha                = alpha,
    pAdjustMethod        = "BH")
  print(summary(dgeResults[[contrast]]))
  # sorting gene list according to significance
  dgeResults[[contrast]] <- dgeResults[[contrast]][order(dgeResults[[contrast]]$pvalue, decreasing = F),]
}
contrast

################ DGE - Save results ######################## 
# Save results
f=paste(folder,"dgeResults", sep='')
dir.create(f, showWarnings=TRUE, recursive=TRUE)
geneidColname <- 'Geneid'
lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))


dgeResults_table = list()    
dgeResults_table = lapply(
  names(dgeResults),
  function(x) 
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = paste(f,'/DGE_results_LN.xlsx', sep=''), 
           row.names = F,
           asTable = T, 
           sheetName = str_sub(names(dgeResults),1,31)) 



################ DGE - MAplot and Vulcano Plots ######################## 
n.label = 10
FDR = T
pvalue = 0.01
for (condition in names(dgeResults)) {
  results = as.data.frame(dgeResults[[condition]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }        
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }        
  
  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))
  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(1,.5); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)
  
  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]
  
  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", condition)) +
    scale_color_manual(values = mycolors) +
    #scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +  
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")  
  
  #print(MAplot)
  #pdf(paste(f,'/','MAplot_',condition,'.pdf',sep=''),width=8, height=6.5)
  #plot(MAplot)
  #dev.off()
  
  # Vulcano plot 
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    ggtitle(paste("Vulcano Plot,", condition)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[c(up_label, down_label),], 
                     aes(x = log2FoldChange, y = -log10(padj), color = DE), 
                     label = row.names(results[c(up_label, down_label),]), size = 2, max.overlaps = 100) + #,
    #box.padding = unit(0.35, "lines"), point.padding = unit(0.6, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value")   
  
  print(VP)
  pdf(paste(f,'/','VulcanoPlot_',condition,'.pdf',sep=''),width=8, height=8)
  plot(VP)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;   
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  
  pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  dev.off()
  
  A1 <- image_read_pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''), density = 140)
  image_write(A1, path = paste(f,'/','MA_VP_',condition,'.tiff',sep=''), format = "tiff")
  
}

####### FDR ########
# save FDR genes in a list
fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)       

fdrDW = list()
fdrDW = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)

####### SEQC ########
# save SEQC genes in a list
seqcUP = list()

pvalue = 0.01
seqcUP = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange > 1])
names(seqcUP)= names(dgeResults)               

seqcDW = list()
seqcDW = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange < -1])
names(seqcDW)= names(dgeResults)
names(dgeResults)

print('FDRup');lengths(fdrUP); print('FDRdw');lengths(fdrDW)
print('SEQCup');lengths(seqcUP); print('SEQCdw');lengths(seqcDW)

print('FDR_tot');lengths(fdrUP) + lengths(fdrDW)
print('SEQC_tot'); lengths(seqcUP) + lengths(seqcDW)

table(metadata$condition)

g = c('ITGA2', 'YAP1', 'ILK', 'MST1', 'TAZ')

df = full_metadata_sel
my_comparisons <- list(c("CRPC", "NEPC") )
frpkm <- cpm(round(txi_sel$counts)); ylabel = "CPM";dir= paste(folder,'gene_CPM/',sep='')
#frpkm <- tpm; ylabel = "TPM"; dir= 'gene_TPM/'
dir.create(dir, recursive = TRUE)

for (gene in g) {
  df$rpkm = as.numeric(frpkm[gene,])
  row.names(df) <- df$samplesID
  df$ID = as.character(1:length(df$samplesID))   
  p <- ggboxplot(df, x = "condition", y = "rpkm", 
                 #label = str_sub(metadata$SampleID, 1, 2), repel = T, label.rectangle = T,
                 fill = "condition", color = "condition", width = 0.8,
                 add = c('jitter'), add.params = list(size =3, shape = 18),
                 alpha = 0.5, palette = viridis(2),
                 xlab = "tumor type", ylab = ylabel,
                 short.panel.labs = FALSE, title = gene, 
                 font = list(size = 16, face = 'bold', color = "black")) + 
    theme(text = element_text(size=18))
  
  assign(paste('p',gene,sep='_'),
         p) 
  pdf(paste(dir,'boxplot_gene_',gene,'.pdf',sep=''), width = 6, height = 5)
  print(p)
  dev.off()
  pdf(paste(dir,'boxplot_gene_',gene,'_pvalue.pdf',sep=''), width = 6, height = 5)
  print(p + stat_compare_means(label = "p.value", 
                               comparisons = my_comparisons))
  dev.off()
}

###### enrichment #######

# Prepare lists of genes to run GSEA
# ------- GSEA ----------
outdir = paste(folder, 'GSEA/', sep='')
dir.create(outdir)
for (c in names(dgeResults)) {
  l_ranked = data.frame(GeneID= row.names(dgeResults[[c]]), LogFC = dgeResults[[c]][,'log2FoldChange'])
  l_ranked = l_ranked[order(l_ranked$LogFC, decreasing = T),]
  write.table(l_ranked, file = paste(outdir, c,'_ranked_list.rnk', sep =''), 
              quote = F, row.names= F, col.names = F, sep ='\t')
}

# ------- enrichR ----------
databases <- listEnrichrDbs()
alpha = 0.05
LogFC_T = 3
fdrUPe = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > LogFC_T])
fdrDWe = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < -LogFC_T])
names(fdrUPe)= names(dgeResults); names(fdrDWe)= names(dgeResults)
lengths(fdrUPe)
lengths(fdrDWe)
enrichf = paste(folder, 'enrichR_LogFC',LogFC_T,'/', sep='')
dir.create(enrichf, showWarnings=FALSE, recursive=TRUE)
# enrichment Parameters
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016")
# alpha used in DGE
padj.cutoff = alpha; 

# Perform Enrichment
enrichr.list <- list()

for (i in 1:length(dgeResults)){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUPe[[i]]
  down.genes <- fdrDWe[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste(enrichf,'/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  write.table(down.genes, paste(enrichf,'/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''), 
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste(enrichf,'/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  
  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)

# Write excels files
for (i in 1:length(dgeResults)){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path(enrichf, 
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
    
    # let's plot
    p = "GO_Biological_Process_2018"; st = 20
    pp = plotEnrich(enrichr.list[[i]][[j]][[p]],
                    showTerms = st, numChar = 40, y = "Count", 
                    orderBy = "P.value", title= p)
    ggsave(paste(enrichf, p, "_TOP", st, "_pathways_",j,".png", sep =''),pp)
    p = "KEGG_2016"; st = 20
    pp = plotEnrich(enrichr.list[[i]][[j]][[p]],
                    showTerms = st, numChar = 40, y = "Count", 
                    orderBy = "P.value", title= p)
    pp
    ggsave(paste(enrichf, p, "_TOP", st, "_pathways_",j,".png", sep =''),pp)
    
  }
}

save.image("Results_BelloneM1435.RData")
load("Neuroendocrine_cancer_analysis/Results_BelloneM1435.RData")

# ---------- heatmap FDR --------- 
lengths(fdrUPe)
lengths(fdrDWe)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors_hm  = crp(255)
heatmap_dir_as = paste(enrichf,"CustomHM/", sep='')
dir.create(heatmap_dir_as)
cpm = cpm(counts(dga), log = T)
colnames(full_metadata_sel)
annotation_column <- as.data.frame(full_metadata_sel[,c("condition", "Body.Site", 
                                                        "Purity", "Pathology.Classification.y")])
annotation_column$condition = as.factor(annotation_column$condition)
mycolors_c <- c("#CC0000", "#3366CC");     names(mycolors_c) = levels(annotation_column$condition)
ann_colors
ann_colors = list(
  condition = mycolors_c
  #Purity = scale_color_gradient(low = "yellow", high = "darkblue")
)

row.names(annotation_column) <- full_metadata_sel[,"bind_ID"]
annotation_column
head(cpm)
minH = -2; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
HPv <- pheatmap::pheatmap(cpm[c(fdrUPe[[i]], fdrDWe[[i]]),],
                          scale = 'row',
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 2,
                          cutree_row  = 2,
                          show_rownames = F,
                          show_colnames = T,
                          cellwidth=20, cellheight=0.1,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          breaks = myb,
                          col=myc,
                          filename = paste(heatmap_dir_as,'HM_FDR_LogFC',LogFC_T,'.pdf',sep=''))
dev.off()

# ---------- heatmap Hippo --------- 
names(enrichr.list)

j= "fdr_down"; p = "KEGG_2016"

hippo_n = grep(x =enrichr.list[[i]][[j]][[p]]$Term, pattern = "Hippo")
PI3KAkt_n = grep(x =enrichr.list[[i]][[j]][[p]]$Term, pattern = "PI3K-Akt")
PI3KAkt_n

gene_Hippo = sort(unlist(str_split(enrichr.list[[i]][[j]][[p]][hippo_n,]$Genes, 
                                   pattern = ";")))
gene_PI3KAkt = sort(unlist(str_split(enrichr.list[[i]][[j]][[p]][PI3KAkt_n,]$Genes, 
                                   pattern = ";")))


HPv <- pheatmap::pheatmap(cpm[gene_Hippo,],
                          scale = 'row',
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 4,
                          #cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          cellwidth=10, cellheight=10,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_Hippo_pathway.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[gene_PI3KAkt,],
                          scale = 'row',
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 4,
                          #cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          cellwidth=10, cellheight=10,
                          #fontsize = 12, fontsize_row = 12, fontsize_col = 12, 
                          display_numbers = F,
                          col=colors_hm,
                          filename = paste(heatmap_dir_as,'HM_PI3KAkt_pathway.pdf',sep=''))
dev.off()


#----------- Jaccard distances -----------
#' Jaccard distances
dir = paste(enrichf,'/JaccardPlots/', 
            sep = '')
dir.create(dir)
p.value.thr = 0.005

breaksList = seq(0, 1, by = 0.001)

cutree_rows_values = c(5,10,5)
lf = list.files(paste(enrichf, sep=''), pattern=glob2rx("*.xlsx"))

lf
for (file in lf) {    
  enrichR.file = paste(enrichf,
                       file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  #s = c("GO_Biological_Process_2018","KEGG_2016")
  Pathways.Table = data.frame()
  for (dat in s) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways 
    colnames(M) = gene.all 
    
    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1 
        }            
      }    
    }
    
    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }
      
      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      
      w=5; h=5; fs = 4; cutree_rows_N = 10
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc, 
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               cellwidth = w, cellheight = h,
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               #main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 12,
               fontsize_row = fs,
               filename = paste(dir,file,'_JaccardDist.pdf', sep=''))
    }
  }
}

dev.off()



# --------- CRPCsig51 -------
file_CRPCsig51 = ".../CRPCsig51.xlsx"
CRPCsig51_df = read.xlsx(file_CRPCsig51)
CRPCsig51 = CRPCsig51_df$Gene
CRPCsig51

annotation_column <- as.data.frame(full_metadata_ro[,c("condition")])
colnames(annotation_column) <- c("condition")
annotation_column$pz <- row.names(full_metadata_ro)
annotation_column = annotation_column[order(annotation_column$condition),]

full_metadata_ro = full_metadata_ro[order(full_metadata_ro$condition),]
annotation_column <- full_metadata_ro[,c("condition", "Body.Site","Integrated NEPC Score")]
annotation_column$condition = as.factor(annotation_column$condition)
annotation_column$Body.Site = as.factor(annotation_column$Body.Site)

row.names(annotation_column) = row.names(full_metadata_ro)


myb = seq(-3.5,3.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))

library("MetBrewer")

mycolors_c <- c("#9900FF", "#339900");     
names(mycolors_c) = levels(annotation_column$condition)
mycolors_bs <- met.brewer("Signac", 10);
names(mycolors_bs) = levels(annotation_column$Body.Site)

ann_colors = list(
  condition = mycolors_c ,
  Body.Site = mycolors_bs
)


HPv <- pheatmap::pheatmap(cpm[intersect(CRPCsig51_df$Gene,row.names(cpm)),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPCsig51_cluster.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[intersect(CRPCsig51_df$Gene,row.names(cpm)),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPCsig51.pdf',sep=''))


HPv



load("Neuroendocrine_cancer_analysis/Results_BelloneM1435.RData")

save.image("BelloneM1435.RData")
# ------ prediction ---------
#setwd("~/Documents/KhuranaLab/EX_0208_prostateRnaCluster")

library(devtools)
# devtools::install_github("Lothelab/CMScaller")
library(CMScaller)
library(edgeR)
# BiocManager::install("sva")
library(sva)
library(biomaRt)


################################################################################################################################

# NTP IPM

tmp_org <- read.xlsx("/Users/tascini.annasofia/Downloads/science.abe1505_tables_s1_to_s17/science.abe1505_table_s5.xlsx", 
                     sheet = "Signature gene")
tmp_org_m <- tmp_org[,c("class", "gene_name")]
tmp_org_m$class <- as.factor(tmp_org_m$class)
colnames(tmp_org_m) <- c("class", "probe")

tmp_org_m <- tmp_org_m[,c("probe", "class")]

matrixData <- txi$counts

emat <- ematAdjust(matrixData, normMethod = 'RLE')
res <- ntp(emat, tmp_org_m, doPlot=TRUE, nPerm=1000, seed = 42)
View(res)

write.csv(res,'prediction_Tang_Beltranetal.csv', quote = F)

res <- read.csv('prediction_Tang_Beltranetal.csv')
View(res)
res$bind_ID = res$X
table(res$prediction)

meta = read.xlsx("metadata.xlsx")
View(meta)
prova = merge(meta, res, by = "bind_ID")
prova = prova[,c("bind_ID", "prediction", "NEPC.vs..CRPC", "FDR")]
write.xlsx(prova[prova$FDR > 0.05,], "prediction_NS.samples_2.xlsx")

res$condition <- full_metadata_ro$condition
res_NS <- res[res$FDR > 0.05,]
write.xlsx(res_NS, "prediction_NS.samples.xlsx")
metadata$prediction_Tang <- res$prediction

library(ggplot2)

full_metadata_ro %>%
  group_by(prediction_Tang, condition) %>%
  tally() %>%
 # ggplot(., aes(x = condition, y = n, group = prediction_Tang)) +
#  geom_bar(stat = "identity", position ="dodge", 
 #          fill = c("darkgoldenrod2", "blue3","#FF99FF","chartreuse3"), color = "black") +
#  theme_classic2(16) + xlab("prediction") + ylab( "#N")
  ggplot(., aes(x = condition, y = n, fill = prediction_Tang)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("CRPC-AR" = "#FAD02E", "CRPC-NE" = "#4575B4", "CRPC-SCL" = "#EA4C89", "CRPC-WNT" = "#4DAF4A")) +
  labs(x = "condition", y = "#N", fill = "prediction_Tang") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

full_metadata_ro %>%
  group_by(prediction_Tang, condition) %>%
  tally() %>%
  ggplot(., aes(x = condition, y = n, fill = prediction_Tang)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 4) +
  scale_fill_manual(values = mycolors_p) +
  labs(x = "Condition", y = "#N", fill = "Prediction Tang") +
  theme_classic(14) 

venn(list(row.names(full_metadata_ro), res$X), ilabels = 'counts')
head(res)
res$bind_ID <- res$X
head(full_metadata_ro)
full_metadata_ro
#full_metadata_ro$prediction_Tang <- res$prediction
metal <- merge(full_metadata_ro, res, by = "bind_ID")
row.names(metal) <- metal$bind_ID
metal <- metal[row.names(full_metadata_ro),]
full_metadata_ro$prediction_Tang <- metal$prediction

library(dplyr)
library(ggplot2)
full_metadata_ro %>%
  group_by(condition, prediction_Tang) %>%
  tally() %>%
  ggplot(., aes(x = condition, y = n, fill=prediction_Tang)) +
  geom_bar(stat = "identity", position ="dodge") +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = mycolors_p) +
  theme_classic(16) + xlab("condition") + ylab( "#N")



## ------- CRPC_SCL signature ----------
#library(openxlsx)
file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_20240924_ML.xlsx"
CRPC_SCL_df = read.xlsx(file_Signatures, sheet = "CRPC_SCL", startRow = 4, colNames = FALSE)
CRPC_SCL_df
CRPC_SCL = CRPC_SCL_df$X1
CRPC_SCL

# annotation_column <- as.data.frame(full_metadata_ro[,c("condition", "prediction_Tang")])
# colnames(annotation_column) <- c("condition", "prediction_Tang")
# annotation_column$pz <- row.names(full_metadata_ro)
# annotation_column = annotation_column[order(annotation_column$condition),]

full_metadata_ro = full_metadata_ro[order(full_metadata_ro$condition),]
annotation_column <- full_metadata_ro[,c("condition", "Body.Site","Integrated NEPC Score", "prediction_Tang")]
annotation_column$condition = as.factor(annotation_column$condition)
annotation_column$Body.Site = as.factor(annotation_column$Body.Site)
annotation_column$prediction_Tang = as.factor(annotation_column$prediction_Tang)
row.names(annotation_column) = row.names(full_metadata_ro)


myb = seq(-3.5,3.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))

library("MetBrewer")
colors()
mycolors_c <- c("#9900FF", "chartreuse3");     
names(mycolors_c) = levels(annotation_column$condition)
mycolors_bs <- met.brewer("Signac", 10);
names(mycolors_bs) = levels(annotation_column$Body.Site)
mycolors_p <- c("#eb8700", "#005fa7","#cc5897","#00925c")   
names(mycolors_p) = levels(annotation_column$prediction_Tang)


ann_colors = list(
  condition = mycolors_c ,
  Body.Site = mycolors_bs,
  prediction_Tang = mycolors_p
)


HPv <- pheatmap::pheatmap(cpm[c("YAP1", "ITGA2", intersect(CRPC_SCL,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 3,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_SCL_cluster.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[c("YAP1","ITGA2", intersect(CRPC_SCL,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_SCL.pdf',sep=''))


HPv

## ------- CRPC_NE signature ----------

file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
sheet_names <- getSheetNames(file_Signatures)
print(sheet_names)
CRPC_NE_df = read.xlsx(file_Signatures, sheet = "CRPC_NE", startRow = 4, colNames = FALSE)
CRPC_NE_df
CRPC_NE = CRPC_NE_df$X1
CRPC_NE


HPv <- pheatmap::pheatmap(cpm[c("YAP1", "ITGA2",intersect(CRPC_NE,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 3,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_NE_cluster.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[c("YAP1","ITGA2", intersect(CRPC_NE,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_NE.pdf',sep=''))


HPv

## ------- CRPC_AR signature ----------

file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
sheet_names <- getSheetNames(file_Signatures)
print(sheet_names)
CRPC_AR_df = read.xlsx(file_Signatures, sheet = "CRPC_AR", startRow = 4, colNames = FALSE)
CRPC_AR_df
CRPC_AR = CRPC_AR_df$X1
CRPC_AR



HPv <- pheatmap::pheatmap(cpm[c("YAP1", "ITGA2", intersect(CRPC_AR,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 3,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_AR_cluster.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[c("YAP1",  "ITGA2",intersect(CRPC_AR,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_AR.pdf',sep=''))


HPv

## ------- CRPC_AR signature ----------

file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
sheet_names <- getSheetNames(file_Signatures)
print(sheet_names)
CRPC_WNT_df = read.xlsx(file_Signatures, sheet = "CRPC_WNT", startRow = 4, colNames = FALSE)
CRPC_WNT_df
CRPC_WNT = CRPC_WNT_df$X1
CRPC_WNT


HPv <- pheatmap::pheatmap(cpm[c("YAP1", "ITGA2", intersect(CRPC_WNT,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 3,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_WNT_cluster.pdf',sep=''))
dev.off()

HPv <- pheatmap::pheatmap(cpm[c("YAP1", "ITGA2",intersect(CRPC_WNT,row.names(cpm))),
                              row.names(annotation_column)],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          annotation_col = annotation_column,
                          annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = T,
                          border_color = NA,
                          cellwidth=10, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPC_WNT.pdf',sep=''))


HPv

