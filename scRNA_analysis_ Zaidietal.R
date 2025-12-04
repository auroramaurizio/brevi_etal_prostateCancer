# Zaidi
# https://www.pnas.org/doi/epub/10.1073/pnas.2322203121 
# Data at GSE210358 and GSE264573

setwd("/beegfs/scratch/ric.cosr/tascini.annasofia/BelloneM_1435_Neuroendocrine_TCGA")
library(plyr)
library(dplyr)

library(Seurat)
library(patchwork)
library(stringr)
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(RColorBrewer))
# suppressMessages(library(repr))
suppressMessages(library(cowplot))
# suppressMessages(library("viridis"))
# suppressMessages(library('pals'))
suppressMessages(library("stringr"))
# suppressMessages(library(enrichR))
# suppressMessages(library(ggpubr))
# suppressMessages(library("svglite"))
# suppressMessages(library('SeuratWrappers'))
# suppressMessages(library(SeuratWrappers))
suppressMessages(library(SeuratObject))
suppressMessages(library('SeuratObject'))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
# suppressMessages(library("ggvenn"))
# suppressMessages(library("MetBrewer"))
#install.packages("scCustomize")
#suppressMessages(library("scCustomize"))
suppressPackageStartupMessages({
  library(rlang)
})
# library(pals)


allcells <- readRDS("data/GEO/GSE210358_GSE264573/GSE264573_msk.integrated.remove.cellcycle.allcells.rds")
head(allcells) 

tumor_cells <- readRDS("data/GEO/GSE210358_GSE264573/GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds")

u1 <- DimPlot(allcells, group.by = "coarse_ano") + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))
colorType = c("red2", "blue", "gold1")
names(colorType) <- c("NEPC", "CRPC", "CSPC")
allcells$MSH_ID = str_sub(allcells$patient, start = 5, end = -1)
u2 <- DimPlot(allcells, group.by = "subtype", cols = colorType) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

u3 <- DimPlot(allcells, group.by = "MSH_ID") + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))


head(tumor_cells)
tumor_cells$MSH_ID = str_sub(tumor_cells$patient, start = 5, end = -1)
colorType = c("red2", "blue", "gold1")
names(colorType) <- c("NEPC", "CRPC", "CSPC")
u5 <- DimPlot(tumor_cells, group.by = "subtype", cols = colorType) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

u6 <- DimPlot(tumor_cells, group.by = "MSH_ID") + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

head(allcells)
metaA <- allcells@meta.data
metaT <- tumor_cells@meta.data

table(metaA[row.names(metaT),"coarse_ano"])
head(tumor_cells)
u7 <- DimPlot(tumor_cells, group.by = "coarse_ano") + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))

u5 | u6

(u1 | u2 | u3) / (u5 | u6)
u1 | (u2 / u5) | (u3 / u6)
ggsave("data/Zaidietal_data_description.pdf", device = 'pdf', 
       plot =u1 | (u2 / u5) | (u3 / u6),
       width = 21, height = 10)
ggsave("data/Zaidietal_data_description.png", device = 'png', 
       plot = u1 | (u2 / u5) | (u3 / u6),
       width = 21, height = 10)
ggsave("data/Zaidietal_data_description.svg", device = 'svg', 
       plot = u1 | (u2 / u5) | (u3 / u6),
       width = 21, height = 10)


##### Signatures #####
CRPCsig51_t = read.xlsx("custom/input/CRPCsig51.xlsx")

CRPC_SCL <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", 
                      sheet = "CRPC_SCL", startRow = 3, colNames = FALSE)
YAP_TAZ_pathway <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", 
                             sheet = "YAP_TAZ_pathway", startRow = 3, colNames = FALSE)
YAP_TAZ_Targets <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", 
                             sheet = "YAP_TAZ_Targets", startRow = 3, colNames = FALSE)
CRPC_NE <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                      sheet = "CRPC_NE", startRow = 4, colNames = FALSE)
NEPC_dge <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                      sheet = "NEPC_dge", startRow = 4, colNames = FALSE)

Signatures = list(CRPCsig51 = CRPCsig51_t$Gene,
                  CRPC_SCL = CRPC_SCL$X1, 
                  YAP_TAZ_pathway = YAP_TAZ_pathway$X1,
                  YAP_TAZ_Targets = YAP_TAZ_Targets$X1,
                  CRPC_NE = CRPC_NE$X1,
                  NEPC_dge = NEPC_dge$X1)


allcells  <- AddModuleScore(object =  allcells, 
                            features = Signatures,
                            ctrl = 5, 
                            name = names(Signatures))

tumor_cells  <- AddModuleScore(object =  tumor_cells, 
                            features = Signatures,
                            ctrl = 5, 
                            name = names(Signatures))
  
pf = FeaturePlot(tumor_cells, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                 order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) #& 
#scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
pf
pf = FeaturePlot(tumor_cells, features = c("YAP1"), 
                 min.cutoff = 'q5', max.cutoff = 'q95', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf  

pf_YAP1 = FeaturePlot(tumor_cells, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                      order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
pf_EZH1 = FeaturePlot(tumor_cells, features = c("EZH1"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_EZH2 = FeaturePlot(tumor_cells, features = c("EZH2"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_SYP = FeaturePlot(tumor_cells, features = c("SYP"), 
                     min.cutoff = 'q10', max.cutoff = 'q90', 
                     order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_ITGA2 = FeaturePlot(tumor_cells, features = c("ITGA2"), 
                       min.cutoff = 'q10', max.cutoff = 'q90', 
                       order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_NDRG1 = FeaturePlot(tumor_cells, features = c("NDRG1"), 
                       min.cutoff = 'q10', max.cutoff = 'q90', 
                       order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

(pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1)
ggsave("data/Zaidietal_expression_of_genes_of_interest.pdf", device = 'pdf', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Zaidietal_expression_of_genes_of_interest.png", device = 'png', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Zaidietal_expression_of_genes_of_interest.svg", device = 'svg', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)


pf_YAP1_f = FeaturePlot(tumor_cells, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                        order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_SYP_f = FeaturePlot(tumor_cells, features = c("SYP"), min.cutoff = "q5", max.cutoff = 'q95', 
                       order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_YAP1_f | pf_SYP_f
dir.create("data/Paper_Figure/", recursive = TRUE)
ggsave("data/Paper_Figure/Figure_B_Zaidietal_YAP1_SYP_expression.pdf", device = 'pdf', 
       plot = pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)
ggsave("data/Paper_Figure/Figure_B_Zaidietal_YAP1_SYP_expression.png", device = 'png', 
       plot =  pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)
ggsave("data/Paper_Figure/Figure_B_Zaidietal_YAP1_SYP_expression.svg", device = 'svg', 
       plot = pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)



Idents(tumor_cells) <- "subtype"
vp = VlnPlot(tumor_cells, features = "YAP1", cols = colorType)
vp = vp + theme(axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=15, hjust =1, vjust = .5), 
                axis.title.x = element_text(face = "bold", color = "black", size = 15),
                axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
                axis.title.y = element_text(face = "bold", color = "black", size = 15),
                legend.text = element_text(face = "bold", color = "black", size = 15),
                legend.position="top") + NoLegend()

vp

AveSeu = AverageExpression(object = tumor_cells, return.seurat = T)
#levels(AveSeu) <- c("Basal", "Luminal (PSA-low)","Luminal (PSA-high)", "SCNC", "CRPC (AR-high)", "mCRPC (PSA-high)", "Non epithelial")
hm = DoHeatmap(AveSeu, 
               size = 10,
               features = "YAP1",#slot = "data",
               #group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
hm
avehm = DoHeatmap(AveSeu, 
                  size = 5,
                  features = "YAP1",#slot = "data",
                  #group.colors = col,
                  #disp.min = -1.5,
                  #disp.max = 1.5,
                  draw.lines = F,
                  group.bar.height = 0.05,
                  angle = 90) + NoLegend() & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

avehm 
vp2 = vp + theme(axis.text.x = element_text(size=0)) 
vp2 / avehm

head(tumor_cells)
CRPCsig51 <- FeaturePlot(tumor_cells, 
            features = "CRPCsig511", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

CRPC_SCL <- FeaturePlot(tumor_cells, 
                          features = "CRPC_SCL2", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

YAP_TAZ_pathway <- FeaturePlot(tumor_cells, 
                        features = "YAP_TAZ_pathway3", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

YAP_TAZ_Targets <- FeaturePlot(tumor_cells, 
                        features = "YAP_TAZ_Targets4", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

CRPC_NE <- FeaturePlot(tumor_cells, 
                        features = "CRPC_NE5", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

NEPC_dge <- FeaturePlot(tumor_cells, 
                       features = "NEPC_dge6", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

(CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
  (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | u5 + ggtitle("UMAP plot"))


ggsave("data/Zaidietal_expression_of_signatures_of_interest.pdf", device = 'pdf', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | u5 + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Zaidietal_expression_of_signatures_of_interest.png", device = 'png', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | u5 + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Zaidietal_expression_of_signatures_of_interest.svg", device = 'svg', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | u5 + ggtitle("UMAP plot")),
       width = 18, height = 10)


Signature_Zaidi_vertical = (CRPCsig51 + ggtitle("CRPCsig51")) / 
  (CRPC_SCL + ggtitle("CRPC_SCL")) / 
  (CRPC_NE + ggtitle("CRPC_NE")) / 
  (YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))

Signature_Zaidi_horizontal = (CRPCsig51 + ggtitle("CRPCsig51")) | 
  (CRPC_SCL + ggtitle("CRPC_SCL")) |
  (CRPC_NE + ggtitle("CRPC_NE")) |
  (YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))


ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest_v.pdf", 
       device = 'pdf', 
       plot = Signature_Zaidi_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest_v.png", 
       device = 'png', 
       plot = Signature_Zaidi_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest_v.svg", 
       device = 'svg', 
       plot = Signature_Zaidi_vertical,
       width = 5, height = 18)


ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest.pdf", 
       device = 'pdf', 
       plot = Signature_Zaidi_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest.png", 
       device = 'png', 
       plot = Signature_Zaidi_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Zaidietal_expression_of_signatures_of_interest.svg", 
       device = 'svg', 
       plot = Signature_Zaidi_horizontal,
       width = 20, height = 5)


table(tumor_cells$subtype)

sub_CRPC <- subset(tumor_cells, subset = `subtype` == "CRPC")
DimPlot(sub_CRPC)
FeaturePlot(sub_CRPC, features = "CRPC_NE5",order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(sub_CRPC, features = "NEPC_dge6",order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(allcells, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE)
ggsave("data/Zaidietal_Blend_YAP1_ITGA2.pdf", device = 'pdf', 
       plot =FeaturePlot(tumor_cells, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)
ggsave("data/Zaidietal_Blend_YAP1_ITGA2.png", device = 'png', 
       plot = FeaturePlot(tumor_cells, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)
ggsave("data/Zaidietal_Blend_YAP1_ITGA2.svg", device = 'svg', 
       plot = FeaturePlot(tumor_cells, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)

FeaturePlot(tumor_cells, features = c("YAP1", "TNC"), order = TRUE, pt.size = 2, blend = TRUE)



FeaturePlot(subset(tumor_cells, subset = YAP1 > 0), 
            features = c("ITGA2"), order = TRUE, pt.size = 2)
ggsave("data/Zaidietal_YAP1p_ITGA2_expression.pdf", device = 'pdf', 
       plot = FeaturePlot(subset(tumor_cells, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2),
       width = 6, height = 5)
ggsave("data/Zaidietal_YAP1p_ITGA2_expression.png", device = 'png', 
       plot = FeaturePlot(subset(tumor_cells, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2), 
       width = 6, height = 5)
ggsave("data/Zaidietal_YAP1p_ITGA2_expression.svg", device = 'svg', 
       plot = FeaturePlot(subset(tumor_cells, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2),
       width = 6, height = 5)

FeaturePlot(tumor_cells, features = c("YAP1", "HCFC1"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "ADAM10"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "PLXNB2"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "COL4A4"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "LAMC1"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(tumor_cells, features = c("YAP1", "CDH1"), order = TRUE, pt.size = 2, blend = TRUE)


pf = FeaturePlot(sub_CRPC, features = c("YAP1"), 
                 min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

Idents(tumor_cells) <- "subtype"
hm_tumor_genes = DoHeatmap(tumor_cells, 
               slot = "data",
               size = 6, group.colors = colorType,
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 2,
               draw.lines = T,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

ggsave("data/Zaidietal_heatmap_of_gene_of_interest.png", device = 'png', 
       plot =hm_tumor_genes,
       width = 10, height = 5)
ggsave("data/Zaidietal_heatmap_of_gene_of_interest.pdf", device = 'pdf', 
       plot =hm_tumor_genes,
       width = 10, height = 5)
ggsave("data/Zaidietal_heatmap_of_gene_of_interest.svg", device = 'svg', 
       plot =hm_tumor_genes,
       width = 10, height = 5)

# ggsave("data/Zaidietal_hm_interestingGenes.png", hm, width = 20, height = 8)

hm_CRPCsig51 = DoHeatmap(tumor_cells, 
               slot = "data",
               size = 6, group.colors = colorType,
               features = c("YAP1",Signatures$CRPCsig51),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_CRPC_SCL = DoHeatmap(tumor_cells, 
                         slot = "data",
                         size = 6, group.colors = colorType,
                         features = c("YAP1",Signatures$CRPC_SCL),
                         #group.colors = col,
                         #disp.min = -1.5,
                         disp.max = 1.5,
                         draw.lines = T,
                         group.bar.height = 0.05,
                        raster = FALSE,
                         angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_CRPC_NE = DoHeatmap(tumor_cells, 
                        slot = "data",
                       raster = FALSE,
                        size = 6, group.colors = colorType,
                        features = c("YAP1",Signatures$CRPC_NE),
                        #group.colors = col,
                        #disp.min = -1.5,
                        disp.max = 1.5,
                        draw.lines = T,
                        group.bar.height = 0.05,
                        angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)


ggsave("data/Zaidietal_heatmap_of_signature_of_interest_1.pdf", device = 'pdf', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)

ggsave("data/Zaidietal_heatmap_of_signature_of_interest_1.png", device = 'png', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)

ggsave("data/Zaidietal_heatmap_of_signature_of_interest_1.svg", device = 'svg', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)

# ggsave("data/Zaidietal_hm_interestingSign_1.png", hm_CRPCsig51 | hm_CRPC_SCL | hm_CRPC_NE , width = 18, height = 16)


hm_YAP_TAZ_pathway = DoHeatmap(tumor_cells, 
                        slot = "data",
                        size = 6, group.colors = colorType,
                        features = c("YAP1",Signatures$YAP_TAZ_pathway),
                        #group.colors = col,
                        #disp.min = -1.5,
                        disp.max = 1.5,
                        draw.lines = T,
                        raster = FALSE,
                        group.bar.height = 0.05,
                        angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_YAP_TAZ_Targets = DoHeatmap(tumor_cells, 
                       slot = "data",
                       size = 6, group.colors = colorType,
                       features = c("YAP1",Signatures$YAP_TAZ_Targets),
                       #group.colors = col,
                       #disp.min = -1.5,
                       disp.max = 1.5,
                       raster = FALSE,
                       draw.lines = T,
                       group.bar.height = 0.05,
                       angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

ggsave("data/Zaidietal_heatmap_of_signature_of_interest_2.pdf", device = 'pdf', 
       plot = hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 9)
ggsave("data/Zaidietal_heatmap_of_signature_of_interest_2.png", device = 'png', 
       plot = hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 9)

# ggsave("data/Zaidietal_hm_interestingSign_2.png", hm_YAP_TAZ_pathway | hm_YAP_TAZ_Targets , width = 18, height = 8)


sub_CRPC_YAP = subset(sub_CRPC, subset = YAP1 > 0)
FeaturePlot(sub_CRPC_YAP, features = "CRPC_NE5",order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

head()
p <- DimPlot(sub_CRPC, group.by = "MSH_ID")
df = p[[1]]$data
df_f <- df[df$UMAP_1 > 6 &df$UMAP_1 < 12 & df$UMAP_2 > 9,]
cellID <- row.names(df_f)
length(cellID)
magic67 <- subset(tumor_cells, cells = cellID)
magic67
DimPlot(magic67, group.by = "MSH_ID")
DimPlot(tumor_cells, cells.highlight = cellID)

ggsave("data/Zaidietal_67cells_highlights.pdf", device = 'pdf', 
       plot =DimPlot(tumor_cells, cells.highlight = cellID),
       width = 6, height = 4)
ggsave("data/Zaidietal_67cells_highlights.png", device = 'png', 
       plot = DimPlot(tumor_cells, cells.highlight = cellID),
       width = 6, height = 4)
ggsave("data/Zaidietal_67cells_highlights.svg", device = 'svg', 
       plot = DimPlot(tumor_cells, cells.highlight = cellID),
       width = 6, height = 4)


(CRPCsig51 | CRPC_SCL) / (YAP_TAZ_pathway | YAP_TAZ_Targets)
head(magic67.scaled)
Idents(magic67) <- "sample"
magic67.scaled <- ScaleData(magic67, features = c("EZH1","EZH2","YAP1", as.character(unlist(Signatures))))
hm = DoHeatmap(magic67.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", Signatures$CRPCsig51),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(magic67, slot = "data",
               size = 8,
               features = c("EZH1","EZH2","YAP1", Signatures$CRPC_NE),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

magic67_YAP1 <- subset(magic67, subset = YAP1 > 0)

magic67_YAP1.scaled <- ScaleData(magic67_YAP1, features = c("EZH1","EZH2","YAP1", as.character(unlist(Signatures))))
hm = DoHeatmap(magic67_YAP1, slot = "data",
               size = 8,
               features = c("EZH1","EZH2","YAP1",Signatures$CRPC_NE),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


hm = DoHeatmap(magic67, slot = "data",
               size = 8,
               features = c("EZH1","EZH2","YAP1",intersect(Signatures$CRPC_NE, Signatures$NEPC_dge)),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

row.names(tumor_cells)[grep("ITG", x = row.names(tumor_cells))]

FeaturePlot(tumor_cells, features = "ITGA2",  
            min.cutoff = 'q10', max.cutoff = 'q90', 
            order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf


Idents(magic67) <- "sample"

hm_CRPCsig51 = DoHeatmap(magic67, slot = "data",
               size = 4,
               features = c("YAP1", Signatures$CRPCsig51),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = F,
               raster = FALSE,
               group.bar.height = 0.05,
               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_CRPC_SCL = DoHeatmap(magic67, slot = "data",
                         size = 4,
                         features = c("YAP1", Signatures$CRPC_SCL),
                         #group.colors = col,
                         #disp.min = -1.5,
                         disp.max = 3,
                         draw.lines = F,
                        raster = FALSE,
                         group.bar.height = 0.05,
                         angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_CRPC_NE = DoHeatmap(magic67, slot = "data",
                        size = 4,
                       raster = FALSE,
                        features = c("YAP1", Signatures$CRPC_NE),
                        #group.colors = col,
                        #disp.min = -1.5,
                        disp.max = 3,
                        draw.lines = F,
                        group.bar.height = 0.05,
                        angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_YAP_TAZ_pathway = DoHeatmap(magic67, slot = "data",
                       size = 4,
                       raster = FALSE,
                       features = c("YAP1", Signatures$YAP_TAZ_pathway),
                       #group.colors = col,
                       #disp.min = -1.5,
                       disp.max = 3,
                       draw.lines = F,
                       group.bar.height = 0.05,
                       angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_YAP_TAZ_Targets = DoHeatmap(magic67, slot = "data",
                               size = 4,
                               raster = FALSE,
                               features = c("YAP1", Signatures$YAP_TAZ_Targets),
                               #group.colors = col,
                               #disp.min = -1.5,
                               disp.max = 3,
                               draw.lines = F,
                               group.bar.height = 0.05,
                               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)


ggsave("data/Zaidietal_67cells_heatmap_of_signatures_of_interest.pdf", device = 'pdf', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL")  | hm_CRPC_NE + ggtitle ("CRPC_NE") | hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 18, height = 14)
ggsave("data/Zaidietal_67cells_heatmap_of_signatures_of_interest.png", device = 'png', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL")  | hm_CRPC_NE + ggtitle ("CRPC_NE") | hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),  
       width = 18, height = 14)
ggsave("data/Zaidietal_67cells_heatmap_of_signatures_of_interest.svg", device = 'svg', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL")  | hm_CRPC_NE + ggtitle ("CRPC_NE") | hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 18, height = 14)

# ggsave("data/Zaidietal_magic67_hm_interestingSign_1.png", hm_CRPCsig51 | hm_CRPC_SCL | hm_CRPC_NE |  hm_YAP_TAZ_pathway | hm_YAP_TAZ_Targets, 
#       width = 24, height = 14)
Idents(magic67) <- "MSH_ID"
hm_tumor_genes_67 = DoHeatmap(magic67, 
               slot = "data",
               raster = FALSE,
               size = 6, group.colors = colorType,
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 2,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)


#ggsave("data/Zaidietal_magic67_hm_interestingGenes.png", hm, width = 18, height = 8)
ggsave("data/Zaidietal_67cells_heatmap_of_gene_of_interest.png", device = 'png', 
       plot =hm_tumor_genes_67,
       width = 12, height = 6)
ggsave("data/Zaidietal_67cells_heatmap_of_gene_of_interest.pdf", device = 'pdf', 
       plot =hm_tumor_genes_67,
       width = 12, height = 6)
ggsave("data/Zaidietal_67cells_heatmap_of_gene_of_interest.svg", device = 'svg', 
       plot =hm_tumor_genes_67,
       width = 12, height = 6)

hm_tumor_genes_67_f = DoHeatmap(magic67, 
                     slot = "data",
                     size = 6, label = TRUE, 
                     features = c("YAP1","EZH1","EZH2","SYP", "NDRG1"),
                     group.colors = c("aquamarine1", "darkorchid1"),
                     disp.max = 3,
                     draw.lines = T,
                     group.bar.height = 0.1,
                     raster = FALSE,
                     angle = 45) & scale_fill_gradient2(low = 'white', mid = "white", high = "firebrick2", midpoint = 0)

hm_tumor_genes_67_f + theme(text = element_text(size = 20))
hm_tumor_genes_67_f + theme(axis.text.y = element_text(size = 24))
ggsave("data/Paper_Figure/Supplementaty_Zaidietal_67_heatmap_of_genes_of_interest.pdf", device = 'pdf', 
       plot = hm_tumor_genes_67_f + theme(axis.text.y = element_text(size = 24)),
       width = 6, height = 5)
ggsave("data/Paper_Figure/Supplementaty_Zaidietal_67_heatmap_of_genes_of_interest.png", device = 'png', 
       plot = hm_tumor_genes_67_f + theme(axis.text.y = element_text(size = 24)),
       width = 6,  height = 5)
ggsave("data/Paper_Figure/Supplementaty_Zaidietal_67_heatmap_of_genes_of_interest.svg", device = 'svg', 
       plot = hm_tumor_genes_67_f + theme(text = element_text(size = 24)),
       width = 6,  height = 5)
