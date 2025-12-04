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
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library("stringr"))
suppressMessages(library(enrichR))
suppressMessages(library(ggpubr))
#suppressMessages(library("svglite"))
#suppressMessages(library('SeuratWrappers'))
#suppressMessages(library(SeuratWrappers))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(magrittr))
suppressMessages(library("ggvenn"))
suppressMessages(library("MetBrewer"))
suppressPackageStartupMessages({
  library(rlang)
})
library(pals)

dge_E <- readRDS(file = "data/GEO/GSE176031/dge_E.rds")
dge_pca <- readRDS(file = "data/GEO/GSE176031/dge_pca.rds")

head(dge_pca)
table(str_sub(dge_pca$sample, start = 1, end = 8))
DimPlot(dge_pca) | DimPlot(dge_E)

ggsave("data/Songetal_data_description.pdf", device = 'pdf', 
       plot = DimPlot(dge_pca) | DimPlot(dge_E),
       width = 16, height = 6)
ggsave("data/Songetal_data_description.png", device = 'png', 
       plot = DimPlot(dge_pca) | DimPlot(dge_E),
       width = 16, height = 6)
ggsave("data/Songetal_data_description.svg", device = 'svg', 
       plot = DimPlot(dge_pca) | DimPlot(dge_E),
       width = 16, height = 6)


#pf = FeaturePlot(samples.aggr, features = c("YAP1"), min.cutoff = 1, order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf_YAP1 = FeaturePlot(dge_E, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                      order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
pf_EZH1 = FeaturePlot(dge_E, features = c("EZH1"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_EZH2 = FeaturePlot(dge_E, features = c("EZH2"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_SYP = FeaturePlot(dge_E, features = c("SYP"), 
                     min.cutoff = 'q10', max.cutoff = 'q90', 
                     order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_ITGA2 = FeaturePlot(dge_E, features = c("ITGA2"), 
                       min.cutoff = 'q10', max.cutoff = 'q90', 
                       order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_NDRG1 = FeaturePlot(dge_E, features = c("NDRG1"), 
                       min.cutoff = 'q10', max.cutoff = 'q90', 
                       order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

(pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1)
(pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1)
ggsave("data/Songetal_expression_of_genes_of_interest.pdf", device = 'pdf', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Songetal_expression_of_genes_of_interest.png", device = 'png', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Songetal_expression_of_genes_of_interest.svg", device = 'svg', 
       plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)


pf_YAP1_f = FeaturePlot(dge_E, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                        order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_SYP_f = FeaturePlot(dge_E, features = c("SYP"), min.cutoff = "q5", max.cutoff = 'q95', 
                       order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_YAP1_f | pf_SYP_f
dir.create("data/Paper_Figure/", recursive = TRUE)
ggsave("data/Paper_Figure/Supplementary_Songetal_YAP1_SYP_expression.pdf", device = 'pdf', 
       plot = pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4.5)
ggsave("data/Paper_Figure/Supplementary_Songetal_YAP1_SYP_expression.png", device = 'png', 
       plot =  pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4.5)
ggsave("data/Paper_Figure/Supplementary_Songetal_YAP1_SYP_expression.svg", device = 'svg', 
       plot =  pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4.5)



pf = FeaturePlot(dge_pca, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                 order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) #& 
#scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
pf
pf = FeaturePlot(dge_pca, features = c("YAP1"), 
                 min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf 

pf = FeaturePlot(dge_pca, features = c("ITGA2"), 
                 min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf 
FeaturePlot(dge_pca, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE)


ggsave("data/Songetal_Blend_YAP1_ITGA2.pdf", device = 'pdf', 
       plot =FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)
ggsave("data/Songetal_Blend_YAP1_ITGA2.png", device = 'png', 
       plot = FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)
ggsave("data/Songetal_Blend_YAP1_ITGA2.svg", device = 'svg', 
       plot = FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 19, height = 5)

ggsave("data/Paper_Figure/Figure_C_Songetal_Blend_YAP1_ITGA2.pdf", device = 'pdf', 
       plot =FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 18, height = 5)
ggsave("data/Paper_Figure/Figure_C_Songetal_Blend_YAP1_ITGA2.png", device = 'png', 
       plot = FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 18, height = 5)
ggsave("data/Paper_Figure/Figure_C_Songetal_Blend_YAP1_ITGA2.svg", device = 'svg', 
       plot = FeaturePlot(dge_E, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE),
       width = 18, height = 5)


FeaturePlot(subset(dge_E, subset = YAP1 > 0), 
            features = c("ITGA2"), order = TRUE, pt.size = 2)
ggsave("data/Songetal_YAP1p_ITGA2_expression.pdf", device = 'pdf', 
       plot = FeaturePlot(subset(dge_E, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2),
       width = 6, height = 5)
ggsave("data/Songetal_YAP1p_ITGA2_expression.png", device = 'png', 
       plot = FeaturePlot(subset(dge_E, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2), 
       width = 6, height = 5)
ggsave("data/Songetal_YAP1p_ITGA2_expression.svg", device = 'svg', 
       plot = FeaturePlot(subset(dge_E, subset = YAP1 > 0), 
                          features = c("ITGA2"), order = TRUE, pt.size = 2),
       width = 6, height = 5)

DimPlot(dge_E)

vp = VlnPlot(dge_pca, features = "YAP1") #, cols = col)
vp = vp + theme(axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=15, hjust =1, vjust = .5), 
                axis.title.x = element_text(face = "bold", color = "black", size = 15),
                axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
                axis.title.y = element_text(face = "bold", color = "black", size = 15),
                legend.text = element_text(face = "bold", color = "black", size = 15),
                legend.position="top") + NoLegend()

vp
AveSeu = AverageExpression(object = dge_pca, return.seurat = T)
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

##### CRPCsig51.xlsx #####
CRPCsig51_t = read.xlsx("custom/input/CRPCsig51.xlsx")
CRPCsig51 = list(sign = CRPCsig51_t$Gene)
CRPCsig51
dge_pca  <- AddModuleScore(object = dge_pca, features = CRPCsig51,
                                ctrl = 5, name = "CRPCsig5")
head(dge_pca@meta.data)
summary(dge_pca@meta.data$CRPCsig51)
pf = FeaturePlot(dge_pca, features = "CRPCsig51", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
pf_CRPCsig51 = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_CRPCsig51

table(Idents(dge_pca))
dge_pca.scaled = ScaleData(dge_pca, features = row.names(dge_pca))
head(dge_pca.scaled@meta.data)
#Idents(dge_pca.scaled) <- "cluster"
hm = DoHeatmap(dge_pca.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPCsig51$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPCsig51$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


##### CRPC_SCL #####
CRPC_SCL <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", 
                      sheet = "CRPC_SCL", startRow = 3, colNames = FALSE)
CRPC_SCL = list(sign = CRPC_SCL$X1)
CRPC_SCL
dge_pca  <- AddModuleScore(object = dge_pca, features = CRPC_SCL,
                           ctrl = 5, name = "CRPC_SCL")
head(dge_pca@meta.data)
summary(dge_pca@meta.data$CRPC_SCL1)
pf = FeaturePlot(dge_pca, features = "CRPC_SCL1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

table(Idents(dge_pca))
dge_pca.scaled = ScaleData(dge_pca, features = row.names(dge_pca))
head(dge_pca.scaled@meta.data)
#Idents(dge_pca.scaled) <- "cluster"
hm = DoHeatmap(dge_pca.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_SCL$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_SCL$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

##### CRPC_NE #####
CRPC_NE <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                      sheet = "CRPC_NE", startRow = 4, colNames = FALSE)
CRPC_NE = list(sign = CRPC_NE$X1)
CRPC_NE
dge_pca  <- AddModuleScore(object = dge_pca, features = CRPC_NE,
                           ctrl = 5, name = "CRPC_NE")
head(dge_pca@meta.data)
summary(dge_pca@meta.data$CRPC_SCL1)
pf = FeaturePlot(dge_pca, features = "CRPC_NE1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

table(Idents(dge_pca))
dge_pca.scaled = ScaleData(dge_pca, features = row.names(dge_pca))
head(dge_pca.scaled@meta.data)
#Idents(dge_pca.scaled) <- "cluster"
hm = DoHeatmap(dge_pca.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_NE$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_NE$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


##### YAP_TAZ_pathway #####
YAP_TAZ_pathway <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", sheet = "YAP_TAZ_pathway", startRow = 3, colNames = FALSE)
YAP_TAZ_pathway = list(sign = YAP_TAZ_pathway$X1)
YAP_TAZ_pathway
dge_pca  <- AddModuleScore(object = dge_pca, features = YAP_TAZ_pathway,
                           ctrl = 5, name = "YAP_TAZ_pathway")
head(dge_pca@meta.data)
summary(dge_pca@meta.data$YAP_TAZ_pathway1)
pf = FeaturePlot(dge_pca, features = "YAP_TAZ_pathway1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

table(Idents(dge_pca))
dge_pca.scaled = ScaleData(dge_pca, features = row.names(dge_pca))
head(dge_pca.scaled@meta.data)
#Idents(dge_pca.scaled) <- "cluster"
hm = DoHeatmap(dge_pca.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_pathway$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_pathway$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


#YAP_TAZ_Targets
##### YAP_TAZ_Targets #####
YAP_TAZ_Targets <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", sheet = "YAP_TAZ_Targets", startRow = 3, colNames = FALSE)
YAP_TAZ_Targets = list(sign = YAP_TAZ_Targets$X1)
YAP_TAZ_Targets
dge_pca  <- AddModuleScore(object = dge_pca, features = YAP_TAZ_Targets,
                           ctrl = 5, name = "YAP_TAZ_Targets")
head(dge_pca@meta.data)
summary(dge_pca@meta.data$YAP_TAZ_Targets1)
pf = FeaturePlot(dge_pca, features = "YAP_TAZ_Targets1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

table(Idents(dge_pca))
#dge_pca.scaled = ScaleData(dge_pca, features = row.names(dge_pca))
head(dge_pca.scaled@meta.data)
#Idents(dge_pca.scaled) <- "cluster"
hm = DoHeatmap(dge_pca.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_Targets$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_Targets$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


# --------- ****Ephitelial cells only**** ----------
AveSeu = AverageExpression(object = dge_E, return.seurat = T)
##### CRPCsig51.xlsx #####
CRPCsig51_t = read.xlsx("custom/input/CRPCsig51.xlsx")
CRPCsig51 = list(sign = CRPCsig51_t$Gene)
CRPCsig51
dge_E  <- AddModuleScore(object = dge_E, features = CRPCsig51,
                           ctrl = 5, name = "CRPCsig5")
head(dge_E@meta.data)
summary(dge_E@meta.data$CRPCsig51)
pf_CRPCsig51 = FeaturePlot(dge_E, features = "CRPCsig51", order = T, pt.size = 2, cols = c("grey", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#pf_CRPCsig51 = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
pf_CRPCsig51

table(Idents(dge_E))
dge_E.scaled = ScaleData(dge_E, features = row.names(dge_E))
head(dge_E.scaled@meta.data)
#Idents(dge_E.scaled) <- "cluster"
hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPCsig51$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm_CRPCsig51 = DoHeatmap(dge_E, slot = "data",
               size = 8, label = FALSE,
               features = c("YAP1", CRPCsig51$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_CRPCsig51
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPCsig51$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


##### CRPC_SCL #####
CRPC_SCL <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", 
                      sheet = "CRPC_SCL", startRow = 3, colNames = FALSE)
CRPC_SCL = list(sign = CRPC_SCL$X1)
CRPC_SCL
dge_E  <- AddModuleScore(object = dge_E, features = CRPC_SCL,
                           ctrl = 5, name = "CRPC_SCL")
head(dge_E@meta.data)
summary(dge_E@meta.data$CRPC_SCL1)
pf_CRPC_SCL = FeaturePlot(dge_E, features = "CRPC_SCL1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#pfCRPC_SCL = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
pfCRPC_SCL

table(Idents(dge_E))
dge_E.scaled = ScaleData(dge_E, features = row.names(dge_E))
head(dge_E.scaled@meta.data)
#Idents(dge_E.scaled) <- "cluster"
hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_SCL$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm_CRPC_SCL = DoHeatmap(dge_E, slot = "data",
                         size = 8,label = FALSE,
                         features = c("YAP1", CRPC_SCL$sign),
                         #group.colors = col,
                         disp.min = -1.5,
                         disp.max = 1.5,
                         draw.lines = T,
                         group.bar.height = 0.05,
                        raster = FALSE,
                         angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_CRPC_SCL
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_SCL$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

##### CRPC_NE #####
CRPC_NE <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                     sheet = "CRPC_NE", startRow = 4, colNames = FALSE)
CRPC_NE = list(sign = CRPC_NE$X1)
CRPC_NE
dge_E  <- AddModuleScore(object = dge_E, features = CRPC_NE,
                           ctrl = 5, name = "CRPC_NE")
head(dge_E@meta.data)
summary(dge_E@meta.data$CRPC_SCL1)
pf_CRPC_NE = FeaturePlot(dge_E, features = "CRPC_NE1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#pf_CRPC_NE = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
pf_CRPC_NE

table(Idents(dge_E))
dge_E.scaled = ScaleData(dge_E, features = row.names(dge_E))
head(dge_E.scaled@meta.data)
#Idents(dge_E.scaled) <- "cluster"
hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_NE$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm_CRPC_NE = DoHeatmap(dge_E, slot = "data",
                        size = 8,label = FALSE,
                        features = c("YAP1", CRPC_NE$sign),
                        #group.colors = col,
                        disp.min = -1.5,
                        disp.max = 1.5,
                        draw.lines = T,
                        group.bar.height = 0.05,
                       raster = FALSE,
                        angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_CRPC_NE


# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPC_NE$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

##### NE_dge  #####
NEPC_dge <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                     sheet = "NEPC_dge", startRow = 4, colNames = FALSE)

venn(list(NEPC_dge = NEPC_dge$X1, CRPC_NE = CRPC_NE$X1))
subset_NEmarkers <- intersect( NEPC_dge$X1, CRPC_NE$X1)

NEPC_dge = list(sign = NEPC_dge$X1)
dge_E  <- AddModuleScore(object = dge_E, features = NEPC_dge,
                         ctrl = 5, name = "NEPC_dge")
head(dge_E@meta.data)
summary(dge_E@meta.data$CRPC_SCL1)
pf = FeaturePlot(dge_E, features = "NEPC_dge1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf

hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", subset_NEmarkers),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", subset_NEmarkers),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

##### YAP_TAZ_pathway #####
YAP_TAZ_pathway <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", sheet = "YAP_TAZ_pathway", startRow = 3, colNames = FALSE)
YAP_TAZ_pathway = list(sign = YAP_TAZ_pathway$X1)
YAP_TAZ_pathway
dge_E  <- AddModuleScore(object = dge_E, features = YAP_TAZ_pathway,
                           ctrl = 5, name = "YAP_TAZ_pathway")
head(dge_E@meta.data)
summary(dge_E@meta.data$YAP_TAZ_pathway1)
pf = FeaturePlot(dge_E, features = "YAP_TAZ_pathway1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf_YAP_TAZ_pathway = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_YAP_TAZ_pathway

table(Idents(dge_E))
dge_E.scaled = ScaleData(dge_E, features = row.names(dge_E))
head(dge_E.scaled@meta.data)
#Idents(dge_E.scaled) <- "cluster"
hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_pathway$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


hm_YAP_TAZ_pathway = DoHeatmap(dge_E, slot = "data",
                       size = 8,label = FALSE,
                       features = c("YAP1", YAP_TAZ_pathway$sign),
                       #group.colors = col,
                       disp.min = -1.5,
                       disp.max = 1.5,
                       draw.lines = T, 
                       group.bar.height = 0.05, 
                       raster = FALSE,
                       angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_YAP_TAZ_pathway
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_pathway$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


#YAP_TAZ_Targets
##### YAP_TAZ_Targets #####
YAP_TAZ_Targets <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA.xlsx", sheet = "YAP_TAZ_Targets", startRow = 3, colNames = FALSE)
YAP_TAZ_Targets = list(sign = YAP_TAZ_Targets$X1)
YAP_TAZ_Targets
dge_E  <- AddModuleScore(object = dge_E, features = YAP_TAZ_Targets,
                           ctrl = 5, name = "YAP_TAZ_Targets")
head(dge_E@meta.data)
summary(dge_E@meta.data$YAP_TAZ_Targets1)
pf_YAP_TAZ_Targets = FeaturePlot(dge_E, features = "YAP_TAZ_Targets1", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#pf_YAP_TAZ_Targets = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
pf_YAP_TAZ_Targets

table(Idents(dge_E))
#dge_E.scaled = ScaleData(dge_E, features = row.names(dge_E))
head(dge_E.scaled@meta.data)
#Idents(dge_E.scaled) <- "cluster"
hm = DoHeatmap(dge_E.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_Targets$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


hm_YAP_TAZ_Targets = DoHeatmap(dge_E, slot = "data",
                               size = 8,label = FALSE,
                               features = c("YAP1", YAP_TAZ_Targets$sign),
                               #group.colors = col,
                               disp.min = -1.5,
                               disp.max = 1.5,
                               draw.lines = T,
                               group.bar.height = 0.05,
                               raster = FALSE,
                               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_YAP_TAZ_Targets
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)

hm = DoHeatmap(AveSeu, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", YAP_TAZ_Targets$sign),
               #group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
head(dge_E)
dge_subset <- subset(dge_E, subset = YAP1 > 0)
#FeaturePlot(dge_E, features = "CRPC_SCL1", split.by = "malignancy")
FeaturePlot(dge_subset, features = c("YAP1","EZH1","EZH2", subset_NEmarkers), order = T, pt.size = 2)
DimPlot(dge_subset)
summary(dge_E$NEPC_dge1)

dge_subset
head(dge_subset)


(pf_CRPCsig51 + ggtitle("CRPCsig51") | pfCRPC_SCL + ggtitle("CRPC_SCL") | pf_CRPC_NE + ggtitle("CRPC_NE")) / 
  (pf_YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | DimPlot(dge_E) + ggtitle("UMAP plot"))

ggsave("data/Songetal_expression_of_signatures_of_interest.pdf", device = 'pdf', 
       plot = (pf_CRPCsig51 + ggtitle("CRPCsig51") | pfCRPC_SCL + ggtitle("CRPC_SCL") | pf_CRPC_NE + ggtitle("CRPC_NE")) / 
         (pf_YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | DimPlot(dge_E) + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Songetal_expression_of_signatures_of_interest.png", device = 'png', 
       plot = (pf_CRPCsig51 + ggtitle("CRPCsig51") | pfCRPC_SCL + ggtitle("CRPC_SCL") | pf_CRPC_NE + ggtitle("CRPC_NE")) / 
         (pf_YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | DimPlot(dge_E) + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Songetal_expression_of_signatures_of_interest.svg", device = 'svg', 
       plot = (pf_CRPCsig51 + ggtitle("CRPCsig51") | pfCRPC_SCL + ggtitle("CRPC_SCL") | pf_CRPC_NE + ggtitle("CRPC_NE")) / 
         (pf_YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | DimPlot(dge_E) + ggtitle("UMAP plot")),
       width = 18, height = 10)


Signature_Song_vertical = (pf_CRPCsig51 + ggtitle("CRPCsig51")) / 
  (pf_CRPC_SCL + ggtitle("CRPC_SCL")) / 
  (pf_CRPC_NE + ggtitle("CRPC_NE")) / 
  (pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))

Signature_Song_horizontal = (pf_CRPCsig51 + ggtitle("CRPCsig51")) | 
  (pf_CRPC_SCL + ggtitle("CRPC_SCL")) |
  (pf_CRPC_NE + ggtitle("CRPC_NE")) |
  (pf_YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))

ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest_v.pdf", 
       device = 'pdf', 
       plot = Signature_Song_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest_v.png", 
       device = 'png', 
       plot = Signature_Song_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest_v.svg", 
       device = 'svg', 
       plot = Signature_Song_vertical,
       width = 5, height = 18)


ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest.pdf", 
       device = 'pdf', 
       plot = Signature_Song_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest.png", 
       device = 'png', 
       plot = Signature_Song_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Songetal_expression_of_signatures_of_interest.svg", 
       device = 'svg', 
       plot = Signature_Song_horizontal,
       width = 20, height = 5)

hm_dgeE_genes = DoHeatmap(dge_E, 
               slot = "data",
               size = 6, label = FALSE,
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = T,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_dgeE_genes
ggsave("data/Songetal_heatmap_of_gene_of_interest.png", device = 'png', 
       plot =hm_dgeE_genes,
       width = 10, height = 5)
ggsave("data/Songetal_heatmap_of_gene_of_interest.pdf", device = 'pdf', 
       plot =hm_dgeE_genes,
       width = 10, height = 5)
ggsave("data/Songetal_heatmap_of_gene_of_interest.svg", device = 'svg', 
       plot =hm_dgeE_genes,
       width = 10, height = 5)

hm_dgeE_genes_f = DoHeatmap(dge_E, 
                          slot = "data",
                          size = 6, label = FALSE,
                          features = c("YAP1","EZH1","EZH2","SYP", "NDRG1"),
                          #group.colors = col,
                          #disp.min = -1.5,
                          disp.max = 3,
                          draw.lines = T,
                          group.bar.height = 0.1,
                          raster = FALSE,
                          angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_dgeE_genes_f + theme(axis.text.y = element_text(size = 24))

ggsave("data/Paper_Figure/Supplementaty_Songetal_heatmap_of_genes_of_interest.pdf", device = 'pdf', 
       plot = hm_dgeE_genes_f + theme(axis.text.y = element_text(size = 24)),
       width = 6, height = 5)
ggsave("data/Paper_Figure/Supplementaty_Songetal_heatmap_of_genes_of_interest.png", device = 'png', 
       plot = hm_dgeE_genes_f + theme(axis.text.y = element_text(size = 24)),
       width = 6,  height = 5)
ggsave("data/Paper_Figure/Supplementaty_Songetal_heatmap_of_genes_of_interest.svg", device = 'svg', 
       plot = hm_dgeE_genes_f + theme(text = element_text(size = 24)),
       width = 6,  height = 5)





hm_dgeE_genes_YAPp = DoHeatmap(subset(dge_E, subset = YAP1 > 0), 
               slot = "data",
               size = 0,
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = T, label = FALSE,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
ggsave("data/Songetal_YAPp_heatmap_of_gene_of_interest.png", device = 'png', 
       plot =hm_dgeE_genes_YAPp,
       width = 10, height = 5)
ggsave("data/Songetal_YAPp_heatmap_of_gene_of_interest.pdf", device = 'pdf', 
       plot =hm_dgeE_genes_YAPp,
       width = 10, height = 5)
ggsave("data/Songetal_YAPp_heatmap_of_gene_of_interest.svg", device = 'svg', 
       plot =hm_dgeE_genes_YAPp,
       width = 10, height = 5)

# hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE")

hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE")

ggsave("data/Songetal_heatmap_of_signature_of_interest_1.pdf", device = 'pdf', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)
ggsave("data/Songetal_heatmap_of_signature_of_interest_1.png", device = 'png', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)
ggsave("data/Songetal_heatmap_of_signature_of_interest_1.svg", device = 'svg', 
       plot = hm_CRPCsig51 + ggtitle ("CRPCsig51") | hm_CRPC_SCL + ggtitle ("CRPC_SCL") | hm_CRPC_NE + ggtitle ("CRPC_NE"),
       width = 18, height = 14)


hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets")
ggsave("data/Songetal_heatmap_of_signature_of_interest_2.pdf", device = 'pdf', 
       plot = hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)
ggsave("data/Songetal_heatmap_of_signature_of_interest_2.png", device = 'png', 
       plot = hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)
ggsave("data/Songetal_heatmap_of_signature_of_interest_2.svg", device = 'svg', 
       plot = hm_YAP_TAZ_pathway + ggtitle ("YAP_TAZ_Pathway") | hm_YAP_TAZ_Targets + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)

