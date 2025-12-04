# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE65502", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "111000"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("NE","PC"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE65502", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE65502", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE65502")

# ------ heatmap -----------
count.matrix <- exprs(gset)
head(count.matrix)

featuresD = fData(gset)
table(is.na(featuresD$UniGene.title))

file_CRPCsig51 = "CRPCsig51.xlsx"
CRPCsig51_df = read.xlsx(file_CRPCsig51)
CRPCsig51 = CRPCsig51_df$Gene

CRPCsig51

musGenes <- c("Hmmr", "Tlx3", "Cpeb4")
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

convertMouseGeneList(musGenes)
CRPCsig51_mouse = convertHumanGeneList(CRPCsig51)

info51 = featuresD[featuresD$Gene.symbol %in% CRPCsig51_mouse,]
info51_g = info51$Gene.symbol
info51_s = as.character(info51$ID)
setdiff(CRPCsig51_mouse,info51_g)


myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
annotation_column = data.frame(condition = c("TPIN","TPIN","TPIN",
                                                "TNE","TNE","TNE"),
                                  row.names = colnames(count.matrix))

mycolors_c <- c("#339900","#9900FF");     
names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c)

pheatmap::pheatmap(count.matrix[info51_s,],
                   scale = 'row',
                   color = myc, 
                   breaks = myb,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors, 
                   cluster_rows = F, 
                   cluster_cols = T, 
                   cutree_cols  = 2,
                   cutree_row  = 2,
                   show_rownames = T,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth=40, cellheight=10,
                   fontsize = 10, 
                   #gaps_col = 34,
                   fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F,
                   labels_row = info51_g,
                   #col=colors_hm,
                   filename = paste('TRAMP_CRPCsig51.pdf',sep=''))


# Install and load homologene
if (!requireNamespace("homologene", quietly = TRUE)) {
  install.packages("homologene")
}
library(homologene)

# Convert human gene symbols to mouse
human_genes <- c("TP53", "BRCA1", "MYC")
mouse_genes <- homologene(human_genes, inTax = 9606, outTax = 10090)
print(mouse_genes)

# --------- CRPC_SCL -------
file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
CRPC_SCL_df = read.xlsx(file_Signatures, sheet = "CRPC_SCL", startRow = 4, colNames = FALSE)
CRPC_SCL_df
CRPC_SCL = CRPC_SCL_df$X1
CRPC_SCL
CRPC_SCL_mouse_df <- homologene(CRPC_SCL, inTax = 9606, outTax = 10090)
CRPC_SCL_mouse <- CRPC_SCL_mouse_df$`10090`


info = featuresD[featuresD$Gene.symbol %in% c("Yap1",CRPC_SCL_mouse),]
info_g = info$Gene.symbol
info_s = as.character(info$ID)
setdiff(CRPC_SCL_mouse,info_g)

myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
annotation_column = data.frame(condition = c("TPIN","TPIN","TPIN",
                                             "TNE","TNE","TNE"),
                               row.names = colnames(count.matrix))

annotation_column = data.frame(condition = c("PAC-SC","PAC-SC","PAC-SC",
                                             "PNE-SC","PNE-SC","PNE-SC"),
                               row.names = colnames(count.matrix))
annotation_column$condition <- as.factor(annotation_column$condition)

mycolors_c <- c("#9900FF","#339900");     
names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c)

pheatmap::pheatmap(count.matrix[info_s,],
                   scale = 'row',
                   color = myc, 
                   breaks = myb,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   cutree_cols  = 3,
                   cutree_row  = 2,
                   show_rownames = T,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth=40, cellheight=10,
                   fontsize = 10, 
                   #gaps_col = 34,
                   fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F,
                   labels_row = info_g,
                   #col=colors_hm,
                   filename = paste('TRAMP_CRPC_SCL.pdf',sep=''))

# --------- CRPC_NE -------
file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
CRPC_NE_df = read.xlsx(file_Signatures, sheet = "CRPC_NE", startRow = 4, colNames = FALSE)
CRPC_NE_df
CRPC_NE = CRPC_NE_df$X1
CRPC_NE
CRPC_NE_mouse_df <- homologene(CRPC_NE, inTax = 9606, outTax = 10090)
CRPC_NE_mouse <- CRPC_NE_mouse_df$`10090`

info = featuresD[featuresD$Gene.symbol %in% c("Yap1",CRPC_NE_mouse),]
info_g = info$Gene.symbol
info_s = as.character(info$ID)
setdiff(CRPC_NE_mouse,info_g)

myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
annotation_column = data.frame(condition = c("PAC-SC","PAC-SC","PAC-SC",
                                             "PNE-SC","PNE-SC","PNE-SC"),
                               row.names = colnames(count.matrix))
annotation_column$condition <- as.factor(annotation_column$condition)

mycolors_c <- c("#9900FF","#339900");       
names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c)

pheatmap::pheatmap(count.matrix[info_s,],
                   scale = 'row',
                   color = myc, 
                   breaks = myb,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   cutree_cols  = 2,
                   cutree_row  = 2,
                   show_rownames = T,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth=40, cellheight=10,
                   fontsize = 10, 
                   #gaps_col = 34,
                   fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F,
                   labels_row = info_g,
                   #col=colors_hm,
                   filename = paste('TRAMP_CRPC_NE.pdf',sep=''))


# --------- CRPC_AR -------
file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
CRPC_AR_df = read.xlsx(file_Signatures, sheet = "CRPC_AR", startRow = 4, colNames = FALSE)
CRPC_AR_df
CRPC_AR = CRPC_AR_df$X1
CRPC_AR
CRPC_AR_mouse_df <- homologene(CRPC_AR, inTax = 9606, outTax = 10090)
CRPC_AR_mouse <- CRPC_AR_mouse_df$`10090`

info = featuresD[featuresD$Gene.symbol %in% c("Yap1",CRPC_AR_mouse),]
info_g = info$Gene.symbol
info_s = as.character(info$ID)
setdiff(CRPC_AR_mouse,info_g)

myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
annotation_column = data.frame(condition = c("PAC-SC","PAC-SC","PAC-SC",
                                             "PNE-SC","PNE-SC","PNE-SC"),
                               row.names = colnames(count.matrix))
annotation_column$condition <- as.factor(annotation_column$condition)

mycolors_c <- c("#9900FF","#339900");     
names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c)

pheatmap::pheatmap(count.matrix[info_s,],
                   scale = 'row',
                   color = myc, 
                   breaks = myb,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   cutree_cols  = 2,
                   cutree_row  = 2,
                   show_rownames = T,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth=40, cellheight=10,
                   fontsize = 10, 
                   #gaps_col = 34,
                   fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F,
                   labels_row = info_g,
                   #col=colors_hm,
                   filename = paste('TRAMP_CRPC_AR.pdf',sep=''))

# --------- CRPC_WNT -------
file_Signatures = "AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx"
sheet_names <- getSheetNames(file_Signatures)
print(sheet_names)
CRPC_WNT_df = read.xlsx(file_Signatures, sheet = "CRPC_WNT", startRow = 4, colNames = FALSE)
CRPC_WNT_df
CRPC_WNT = CRPC_WNT_df$X1
CRPC_WNT
CRPC_WNT_mouse_df <- homologene(CRPC_WNT, inTax = 9606, outTax = 10090)
CRPC_WNT_mouse <- CRPC_WNT_mouse_df$`10090`

info = featuresD[featuresD$Gene.symbol %in% c("Yap1",CRPC_WNT_mouse),]
info_g = info$Gene.symbol
info_s = as.character(info$ID)
setdiff(CRPC_WNT_mouse,info_g)

myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
annotation_column = data.frame(condition = c("PAC-SC","PAC-SC","PAC-SC",
                                             "PNE-SC","PNE-SC","PNE-SC"),
                               row.names = colnames(count.matrix))
annotation_column$condition <- as.factor(annotation_column$condition)

mycolors_c <- c("#9900FF","#339900");       
names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c)

pheatmap::pheatmap(count.matrix[info_s,],
                   scale = 'row',
                   color = myc, 
                   breaks = myb,
                   annotation_col = annotation_column,
                   annotation_colors = ann_colors, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   cutree_cols  = 2,
                   cutree_row  = 2,
                   show_rownames = T,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth=40, cellheight=10,
                   fontsize = 10, 
                   #gaps_col = 34,
                   fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F,
                   labels_row = info_g,
                   #col=colors_hm,
                   filename = paste('TRAMP_CRPC_WNT.pdf',sep=''))
