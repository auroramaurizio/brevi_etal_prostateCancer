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

min_nFeature_RNA = 500
max_nFeature_RNA = 8000
max_percent_MT = 10

samples.aggr.data <- Read10X(data.dir = "data/ProstateCancer/outs/count/filtered_feature_bc_matrix/")
samples.aggr <- CreateSeuratObject(counts = samples.aggr.data, project = "Prostate", min.cells = 3, min.features = 200)
rm(samples.aggr.data)
metacsv = read.csv("custom/input/aggr.csv")

head(samples.aggr@meta.data)

samples.aggr$sample <- as.data.frame(str_split_fixed(rownames(samples.aggr@meta.data), pattern = "-", n = 2))[,2]
for (i in 1:length(metacsv$sample_id)) {
  print(i)
  samples.aggr@meta.data["sample"][samples.aggr@meta.data["sample"] == as.character(i)] <- metacsv$sample_id[i]
}


table(samples.aggr$sample)
samples.aggr$sample = as.factor(samples.aggr$sample)
lev = levels(samples.aggr$sample)
str_split(lev, pattern = "_")

samples.aggr$condition = samples.aggr$sample
meta = samples.aggr@meta.data
meta$sample = as.character(meta$sample)
meta$ID = str_remove(string = meta$sample, pattern = "_CD49L")
meta$ID = str_remove(string = meta$ID, pattern = "_CD49H")
meta$ID = str_remove(string = meta$ID, pattern = "_CXCR2")
table(meta$ID)
samples.aggr$ID = as.factor(meta$ID)
samples.aggr$ID = revalue(samples.aggr$ID, 
                          c("PCa2_benign"="PCa2_tumor_adjacentTissue", 
                            "PCa3_benign"="PCa3_tumor_adjacentTissue"))
levels(samples.aggr$ID)
samples.aggr$ID = factor(samples.aggr$ID, levels = c("PCa1_tumor",  
                                                     "PCa2_tumor_adjacentTissue", "PCa2_tumor", 
                                                     "PCa3_tumor_adjacentTissue", "PCa3_tumor",
                                                     "CRPC1", "CRPC2","CRPC3"))

meta$condition = meta$ID
meta$condition[grep(x = meta$ID, pattern = "tumor")] = "Tumor"
meta$condition[grep(x = meta$ID, pattern = "CRPC")] = "Tumor"
meta$condition[grep(x = meta$ID, pattern = "benign")] = "Tumor_adjacentTissue"
table(meta$condition)
samples.aggr$condition = as.factor(meta$condition)


# sample processing
samples.aggr[["percent.mt"]] <- PercentageFeatureSet(samples.aggr, pattern = "^MT-")
VlnPlot(samples.aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = F, group.by = "orig.ident")
#ggsave("plots/vio_plot.png", dpi = 330, scale = 0.8, width = 12, height = 5)


samples.aggr <- subset(samples.aggr, subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)  
samples.aggr

samples.aggr <- NormalizeData(samples.aggr, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = 'vst',
                       nfeatures = 2000,
                       verbose = FALSE)  %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 30, verbose = FALSE)%>%
  RunUMAP(reduction = "pca", dims = 1:20)%>%
  FindNeighbors(reduction = "pca", dims = 1:20)%>%
  FindClusters(resolution = 0.6)

ElbowPlot(object = samples.aggr, ndims = 40)
samples.aggr <- RunUMAP(samples.aggr, reduction = "pca", dims = 1:10) %>%
  FindNeighbors(reduction = "pca", dims = 1:10)%>%
  FindClusters(resolution = 0.6)

# ---------- plots -----------------
dir.create("data/Rplot", recursive = T, showWarnings = F)
colori = c("#9999FF","#33FFFF", "#0000CC","#33FF33", "#339900", "#CC3399", "#FF6666", "#990033")
#polychrome(length(levels(samples.aggr$ID)))
names(colori) = levels(samples.aggr$ID)
p1 = DimPlot(samples.aggr, reduction = "umap", label = F, pt.size = 0.5, group.by = "ID", cols = colori) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=15, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 15),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
        axis.title.y = element_text(face = "bold", color = "black", size = 15),
        legend.text = element_text(face = "bold", color = "black", size = 15),
        legend.position="right") +
  labs(x = "UMAP 1", y = "UMAP 2") 
p1 = p1 + theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
           axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
           axis.title.x = element_text(face = "bold", color = "black", size = 12),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
           axis.title.y = element_text(face = "bold", color = "black", size = 12),
           panel.background = element_rect(fill='transparent'),
           plot.background = element_rect(fill='transparent', color=NA),
           legend.text = element_text(face = "bold", color = "black", size = 12),
           legend.position="right") 
p1
ggsave("data/Rplot/umap_sampleID.pdf", plot = p1, device = "pdf", 
       width = 8, height = 5, bg = "transparent")
ggsave("data/Rplot/umap_sampleID.png", plot = p1, device = "png", 
       width = 8, height = 5, bg = "transparent")


p2 = DimPlot(samples.aggr, reduction = "umap", label = T, pt.size = 0.5) + 
  theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=12, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right") +
  labs(x = "UMAP 1", y = "UMAP 2") 

p2 
head(samples.aggr)

p3 = DimPlot(samples.aggr, reduction = "umap", label = F, pt.size = .2, group.by = "condition",cols = c("#CCCCCC", "#0000FF")) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=15, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 15),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
        axis.title.y = element_text(face = "bold", color = "black", size = 15),
        legend.text = element_text(face = "bold", color = "black", size = 15),
        legend.position="top") +
  labs(x = "UMAP 1", y = "UMAP 2") 
p3 = p3 + theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
           axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
           axis.title.x = element_text(face = "bold", color = "black", size = 12),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
           axis.title.y = element_text(face = "bold", color = "black", size = 12),
           panel.background = element_rect(fill='transparent'),
           plot.background = element_rect(fill='transparent', color=NA),
           legend.text = element_text(face = "bold", color = "black", size = 12),
           legend.position="top") 
p3
ggsave("data/Rplot/umap_Tumor.pdf", plot = p3, device = "pdf", 
       width = 6, height = 5, bg = "transparent")
ggsave("data/Rplot/umap_Tumor.png", plot = p3, device = "png", 
       width = 6, height = 5, bg = "transparent")



#pf = FeaturePlot(samples.aggr, features = c("YAP1"), min.cutoff = 1, order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf_YAP1 = FeaturePlot(samples.aggr, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                 order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

pf_EZH1 = FeaturePlot(samples.aggr, features = c("EZH1"), 
                 min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_EZH2 = FeaturePlot(samples.aggr, features = c("EZH2"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_SYP = FeaturePlot(samples.aggr, features = c("SYP"), 
                      min.cutoff = 'q10', max.cutoff = 'q90', 
                      order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_ITGA2 = FeaturePlot(samples.aggr, features = c("ITGA2"), 
                 min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf_NDRG1 = FeaturePlot(samples.aggr, features = c("NDRG1"), 
                       min.cutoff = 'q10', max.cutoff = 'q90', 
                       order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

(pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1)
ggsave("data/Chengetal_expression_of_genes_of_interest.pdf", device = 'pdf', plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Chengetal_expression_of_genes_of_interest.png", device = 'png', plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)
ggsave("data/Chengetal_expression_of_genes_of_interest.svg", device = 'svg', plot = (pf_YAP1 |  pf_ITGA2  | pf_EZH1) / (pf_EZH2 | pf_SYP | pf_NDRG1),
       width = 18, height = 10)


pf_YAP1_f = FeaturePlot(samples.aggr, features = c("YAP1"), min.cutoff = "q5", max.cutoff = 'q95', 
                      order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_SYP_f = FeaturePlot(samples.aggr, features = c("SYP"), min.cutoff = "q5", max.cutoff = 'q95', 
                        order = T, pt.size = 2, cols = c("lightgrey", "#FF0033")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu")))

pf_YAP1_f | pf_SYP_f
dir.create("data/Paper_Figure/", recursive = TRUE)
ggsave("data/Paper_Figure/Figure_A_Chengetal_YAP1_SYP_expression.pdf", device = 'pdf', 
       plot = pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)
ggsave("data/Paper_Figure/Figure_A_Chengetal_YAP1_SYP_expression.png", device = 'png', 
       plot =  pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)
ggsave("data/Paper_Figure/Figure_A_Chengetal_YAP1_SYP_expression.svg", device = 'svg', 
       plot = pf_YAP1_f | pf_SYP_f,
       width = 10, height = 4)


FeaturePlot(samples.aggr, features = c("YAP1", "ITGA2"), order = TRUE, pt.size = 2, blend = TRUE)
FeaturePlot(subset(samples.aggr, subset = YAP1 > 0), 
            features = c("ITGA2"), order = TRUE, pt.size = 2)


row.names(samples.aggr)[grep("EZH",row.names(samples.aggr))]
p1 + pf

saveRDS(samples.aggr, "data/samples.aggr.RDS")

# --------- load obj for plots and other analysis -------
samples.aggr = readRDS("data/samples.aggr.RDS")

cellCLU= paste("C", 0:17, sep = "")
cellCLU = c(cellCLU[1:3], cellCLU[6], cellCLU[4], cellCLU[7], cellCLU[5],cellCLU[8:18])
cellCLU
names(cellCLU) <- levels(samples.aggr)
samples.aggr <- RenameIdents(samples.aggr, cellCLU)
samples.aggr$cluster = Idents(samples.aggr)
levels(samples.aggr) <- paste("C", 0:17, sep = "")
levels(samples.aggr$cluster) <- paste("C", 0:17, sep = "")
p2 = DimPlot(samples.aggr, reduction = "umap", label = T, pt.size = 0.5, group.by = "cluster") + 
  theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="right") +
  #scale_color_manual(values = met.brewer("Egypt", n = length(levels(samples.aggr$cluster)))) + 
  labs(x = "UMAP 1", y = "UMAP 2") + NoLegend()

p2
ggsave("data/Rplot/umap_cluster.pdf", plot = p2, device = "pdf", 
       width = 6, height = 5, bg = "transparent")
ggsave("data/Rplot/umap_cluster.png", plot = p2, device = "png", 
       width = 6, height = 5, bg = "transparent")

cellID = c("Basal", #0
           "Luminal (PSA-high)", #1
           "mCRPC (PSA-high)", #2
           "Luminal (PSA-high)", #3
           "Basal", #4
           "Luminal (PSA-high)", #5
           "Luminal (PSA-low)", #6
           "Luminal (PSA-high)", #7
           "SCNC", #8
           "Basal", #9
           "Basal", #10
           "Luminal (PSA-low)", #11
           "CRPC (AR-high)", #12
           "Non epithelial", #13
           "Non epithelial", #14
           "Luminal (PSA-high)", #15
           "Non epithelial", #16
           "Non epithelial") #17


names(cellID) <- levels(samples.aggr)
samples.aggr <- RenameIdents(samples.aggr, cellID)
samples.aggr$lineage_subtypes = Idents(samples.aggr)
samples.aggr$celltype = Idents(samples.aggr)

levels(samples.aggr$celltype)
samples.aggr$celltype = factor(samples.aggr$celltype, 
                               c("Basal", "Luminal (PSA-low)","Luminal (PSA-high)", "SCNC", "CRPC (AR-high)", "mCRPC (PSA-high)", "Non epithelial"))

col = c("#339999","#99CCFF", "#0066CC","#FF6300", "#FF3399", "#990000","#CCCCCC")
names(col) = levels(samples.aggr$celltype)
p4 = DimPlot(samples.aggr, reduction = "umap",
             label = F, pt.size = 0.5, group.by = "celltype", cols = col) +
  theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA)
) + labs(x = "UMAP 1", y = "UMAP 2") 
p4 = p4 + theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
                axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =1), 
                axis.title.x = element_text(face = "bold", color = "black", size = 12),
                axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
                axis.title.y = element_text(face = "bold", color = "black", size = 12),
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent', color=NA),
                legend.text = element_text(face = "bold", color = "black", size = 12),
                legend.position="right") 


ggsave("data/Rplot/umap_combined.pdf", plot = p2 | p4 , 
       device = "pdf", width = 11, height = 5,
       bg = "transparent")
ggsave("data/Rplot/umap_combined.svg", plot = p2 | p4 , 
       device = "svg", width = 11, height = 5,
       bg = "transparent")

ggsave("data/Rplot/umap_lineage_subtypes.pdf", plot = p4, 
       device = "pdf", width = 7, height = 5,
       bg = "transparent")
ggsave("data/Rplot/umap_lineage_subtypes.png", plot = p4, 
       device = "png", width = 7, height = 5,
       bg = "transparent")
#p4 + pf

(p2 + ggtitle ("Seurat clusters") | p4 + ggtitle("Lineage subtypes")) / 
  (p3 + ggtitle ("Condition")| p1 + ggtitle ("Sample ID"))


ggsave("data/Chengetal_data_description.pdf", device = 'pdf', 
       plot = (p2 + ggtitle ("Seurat clusters") | p4 + ggtitle("Lineage subtypes")) / 
         (p3 + ggtitle ("Condition")| p1 + ggtitle ("Sample ID")),
       width = 12, height = 8)
ggsave("data/Chengetal_data_description.png", device = 'png', 
       plot = (p2 + ggtitle ("Seurat clusters") | p4 + ggtitle("Lineage subtypes")) / 
         (p3 + ggtitle ("Condition")| p1 + ggtitle ("Sample ID")),
       width = 12, height = 8)
ggsave("data/Chengetal_data_description.svg", device = 'svg', 
       plot = (p2 + ggtitle ("Seurat clusters") | p4 + ggtitle("Lineage subtypes")) / 
         (p3 + ggtitle ("Condition")| p1 + ggtitle ("Sample ID")),
       width = 12, height = 8)

Idents(samples.aggr) <- "celltype"
vp = VlnPlot(samples.aggr, features = "YAP1", cols = col)
vp = vp + theme(axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=15, hjust =1, vjust = .5), 
                axis.title.x = element_text(face = "bold", color = "black", size = 15),
                axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
                axis.title.y = element_text(face = "bold", color = "black", size = 15),
                legend.text = element_text(face = "bold", color = "black", size = 15),
                legend.position="top") + NoLegend()

vp
ggsave("data/Rplot/vp_YAP1_lineage_subtypes.pdf", vp, device = "pdf", width = 12, height = 8)
AveSeu = AverageExpression(object = samples.aggr, return.seurat = T)
levels(AveSeu) <- c("Basal", "Luminal (PSA-low)","Luminal (PSA-high)", "SCNC", "CRPC (AR-high)", "mCRPC (PSA-high)", "Non epithelial")
hm = DoHeatmap(AveSeu, 
               size = 10,
               features = "YAP1",#slot = "data",
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
hm
avehm = DoHeatmap(AveSeu, 
                  size = 5,
                  features = "YAP1",#slot = "data",
                  group.colors = col,
                  #disp.min = -1.5,
                  #disp.max = 1.5,
                  draw.lines = F,
                  group.bar.height = 0.05,
                  angle = 90) + NoLegend() & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

avehm 
ggsave("data/Rplot/hm_YAP1_lineage_subtypes.pdf", avehm , device = "pdf", width = 12, height = 7)

vp2 = vp + theme(axis.text.x = element_text(size=0)) 
vp2 / avehm
ggsave("data/Rplot/hmvp_YAP1_lineage_subtypes.pdf", vp2 / avehm, device = "pdf", width = 12, height =9)
ggsave("data/Rplot/hmvp_YAP1_lineage_subtypes.svg", vp2 / avehm, device = "svg", width = 12, height =9)


Idents(samples.aggr) <- "cluster"
levels(samples.aggr) <- paste("C", 0:17, sep = "")
vp = VlnPlot(samples.aggr, features = "YAP1", pt.size = 0.1)
vp = vp + theme(axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=15, hjust =0.5), 
           axis.title.x = element_text(face = "bold", color = "black", size = 15),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=15),
           axis.title.y = element_text(face = "bold", color = "black", size = 15),
           legend.text = element_text(face = "bold", color = "black", size = 15),
           legend.position="top") + NoLegend()

vp
ggsave("data/Rplot/vp_YAP1.pdf", vp, device = "pdf", width = 12, height = 4)
AveSeuclu = AverageExpression(object = samples.aggr, return.seurat = T)
#samples.aggr.scaled = ScaleData(samples.aggr, features = row.names(samples.aggr))
hm = DoHeatmap(AveSeuclu, 
               size = 10,
               features = "YAP1",#slot = "data",
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) + NoLegend() & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
hm
ggsave("data/Rplot/hm_YAP1.pdf", hm, device = "pdf", width = 12, height = 2)

vp /hm
ggsave("data/Rplot/hmvp_YAP1.pdf", vp /hm, device = "pdf", width = 12, height =8)
ggsave("data/Rplot/hmvp_YAP1.svg", vp /hm, device = "svg", width = 12, height =8)

hm = DoHeatmap(AveSeuclu, 
               size = 10,
               features = "EZH2",#slot = "data",
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) + NoLegend() & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
hm

Idents(samples.aggr) <- "condition"
vp = VlnPlot(samples.aggr, features = "YAP1")
vp
AveSeuco = AverageExpression(object = samples.aggr, return.seurat = T)
#samples.aggr.scaled = ScaleData(samples.aggr, features = row.names(samples.aggr))
hm = DoHeatmap(AveSeuco, 
               size = 5,
               hjust = .5,
               features = "YAP1",#slot = "data",
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
hm

Idents(samples.aggr) <- "ID"
vp = VlnPlot(samples.aggr, features = "YAP1", cols = colori)
AveSeuId = AverageExpression(object = samples.aggr, return.seurat = T)
#samples.aggr.scaled = ScaleData(samples.aggr, features = row.names(samples.aggr))
hm = DoHeatmap(AveSeuId, 
               size = 0,
               #hjust = .5,
               features = "YAP1",#slot = "data",
               group.colors = colori,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

vp / hm


cols.use <- list(celltype= col, 
                 #ID = colori,
                 condition = c("#CCCCCC", "#0000FF"))

hm = DoMultiBarHeatmap(samples.aggr.scaled, 
                       features="YAP1", 
                       group.by= 'celltype', 
                       additional.group.by = c('condition'), 
                       disp.min = -1.5,
                       disp.max = 1.5,
                       label = TRUE,
                       size = 10,
                       draw.lines = T,
                       lines.width = 20,
                       cols.use=cols.use) +
  scale_fill_gradientn(colours = coolwarm(200)) 
hm

# -------- CRPCsig51 --------
CRPCsig51_t = read.xlsx("custom/input/CRPCsig51.xlsx")
CRPCsig51 = list(sign = CRPCsig51_t$Gene)
CRPCsig51
samples.aggr  <- AddModuleScore(object = samples.aggr, features = CRPCsig51,
                                ctrl = 5, name = "CRPCsig5")
head(samples.aggr@meta.data)
pf = FeaturePlot(samples.aggr, features = "CRPCsig51", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300"))
pf_CRPCsig51 = pf  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


table(Idents(samples.aggr))
samples.aggr.scaled = ScaleData(samples.aggr, features = row.names(samples.aggr))
head(samples.aggr.scaled@meta.data)
Idents(samples.aggr.scaled) <- "cluster"
hm = DoHeatmap(samples.aggr.scaled, 
               size = 8,
               features = c("EZH1","EZH2","YAP1", CRPCsig51$sign),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = T,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
# scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# scale_fill_gradientn(colours = myc)
myb = seq(-3.5,3.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
myc = scale_fill_gradient2(low = muted('blue'), mid = "white", high = "red", midpoint = 0)

hm = DoHeatmap(AveSeuclu, 
               size = 8,
               features = c("YAP1", CRPCsig51$sign),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

#---------C12----------



Idents(samples.aggr) <- "cluster"
C12 = subset(samples.aggr, idents = "C12")
FeaturePlot(subset(C12, subset = YAP1 > 0), 
            features = c("ITGA2"), order = TRUE, pt.size = 2)


meta12 = C12@meta.data
meta12$ct = as.character(meta12$sample)
table(meta12$sample)
meta12$ct[grep("CRPC", meta12$sample)] = "CRPC"
meta12$ct[grep("PCa", meta12$sample)] = "PCa"
meta12$ct = as.factor(meta12$ct)
meta12$ct <- factor(meta12$ct, c("PCa", "CRPC"))
meta12 %>%
  group_by(ct) %>%
  tally() %>%
  ggplot(., aes(x="", y=n, fill=ct)) +
  geom_bar(stat="identity",  color="black") + #, position="fill") + 
  scale_fill_manual(values = c("#CC0033", "#3366FF")) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=12, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "C12", y = "Number of cells")


#VlnPlot(C12, features = c("YAP1", "EZH1", "EZH2"), split.by = "ct")
pfy = FeaturePlot(C12, features = c("YAP1"), 
                 min.cutoff = 'q5', max.cutoff = 'q95', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
#  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez1 = FeaturePlot(C12, features = c("EZH1"), 
                  min.cutoff = 'q5', max.cutoff = 'q95', 
                  order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez2 = FeaturePlot(C12, features = c("EZH2"), 
                  min.cutoff = 'q5', max.cutoff = 'q95', 
                  order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 
  

ggsave("data/Rplot/gene_expression.svg", plot = pfy | pfez1 | pfez2,
         device = "svg", 
       width = 12, height = 5, bg = "transparent")


pfy = FeaturePlot(samples.aggr, features = c("YAP1"), 
                  min.cutoff = 'q5', max.cutoff = 'q95', 
                  order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez1 = FeaturePlot(samples.aggr, features = c("EZH1"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez2 = FeaturePlot(samples.aggr, features = c("EZH2"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 
pfy | pfez1 | pfez2

ggsave("data/Rplot/gene_expression_app.svg", plot = pfy | pfez1 | pfez2,
       device = "svg", 
       width = 20, height = 8, bg = "transparent")


C12$ct = meta12$ct
DimPlot(C12, group.by = "ct", pt.size = 3) +
  xlim(c(-7,-4)) + ylim(c(-10,-3)) +
  scale_color_manual(values = c("#CC0033", "#3366FF"))
C12_2 <- RunUMAP(C12, dims = 1:20)

u = DimPlot(C12_2, group.by = "sample", pt.size = 2)

head(C12_2@meta.data)
pf = FeaturePlot(C12_2, features = c("YAP1"), 
                 min.cutoff = 'q10', max.cutoff = 'q50', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf + u

VlnPlot(C12_2, "YAP1", group.by = "sample", pt.size = 2)
VlnPlot(C12_2, "EZH2", group.by = "sample", pt.size = 2)

pf = FeaturePlot(C12_2, features = c("EZH2"), 
                 #min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf + u

C12.scaled = ScaleData(C12, c("EZH1", "EZH2","YAP1", CRPCsig51$sign))
hm = DoHeatmap(C12.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1","YAP1", CRPCsig51$sign),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

C12_data = GetAssay(C12)
metaC12 = C12@meta.data
annotation_column = as.data.frame(metaC12[,c("lineage_subtypes", "ID", "cluster")])
head(annotation_column)
myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
HPv <- pheatmap::pheatmap(C12_data[c("YAP1",CRPCsig51$sign),],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          #annotation_col = annotation_column,
                          #annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = F,
                          border_color = NA,
                          cellwidth=0.5, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPCsig51.pdf',sep=''))
HPv
dev.off()

# -------- prova --------



samples.aggr.scaled@meta.data

cd_features <- list(c(
  'CD79B',
  'CD79A',
  'CD19',
  'CD180',
  'CD200',
  'CD3D',
  'CD2',
  'CD3E',
  'CD7',
  'CD8A',
  'CD14',
  'CD1C',
  'CD68',
  'CD9',
  'CD247'
))
samples.aggr  <- AddModuleScore(
  object = samples.aggr,
  features = cd_features,
  ctrl = 5,
  name = 'CD_Features'
)
head(x = samples.aggr[])




FeaturePlot(samples.aggr, features = c("CD_Features1"), order = T, pt.size = 2)


# ------- function --------
suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}
suppressPackageStartupMessages({
  library(rlang)
})



# new signatures

##### Signatures #####
CRPCsig51_t = read.xlsx("custom/input/CRPCsig51.xlsx")

CRPC_SCL <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                      sheet = "CRPC_SCL", startRow = 4, colNames = FALSE)
YAP_TAZ_pathway <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                             sheet = "YAP_TAZ_pathway", startRow = 4, colNames = FALSE)
YAP_TAZ_Targets <- read.xlsx("custom/input/AdditionalSignatures_BelloneM_1435_Neuroendocrine_TCGA_202412.xlsx", 
                             sheet = "YAP_TAZ_Targets", startRow = 4, colNames = FALSE)
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


samples.aggr  <- AddModuleScore(object =  samples.aggr, 
                            features = Signatures,
                            ctrl = 5, 
                            name = names(Signatures))
head(samples.aggr)
CRPCsig51 <- FeaturePlot(samples.aggr, 
                         features = "CRPCsig511", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

CRPC_SCL <- FeaturePlot(samples.aggr, 
                        features = "CRPC_SCL2", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

YAP_TAZ_pathway <- FeaturePlot(samples.aggr, 
                               features = "YAP_TAZ_pathway3", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

YAP_TAZ_Targets <- FeaturePlot(samples.aggr, 
                               features = "YAP_TAZ_Targets4", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

CRPC_NE <- FeaturePlot(samples.aggr, 
                        features = "CRPC_NE5", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

NEPC_dge <- FeaturePlot(samples.aggr, 
                       features = "NEPC_dge6", order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) + 
  theme(panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



(CRPCsig51 | CRPC_NE | CRPC_SCL) / (YAP_TAZ_pathway | YAP_TAZ_Targets)
Signature_Cheng_vertical = (CRPCsig51 + ggtitle("CRPCsig51")) / 
  (CRPC_SCL + ggtitle("CRPC_SCL")) / 
  (CRPC_NE + ggtitle("CRPC_NE")) / 
  (YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))

Signature_Cheng_horizontal = (CRPCsig51 + ggtitle("CRPCsig51")) | 
  (CRPC_SCL + ggtitle("CRPC_SCL")) |
  (CRPC_NE + ggtitle("CRPC_NE")) |
  (YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets"))


ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest_v.pdf", 
       device = 'pdf', 
       plot = Signature_Cheng_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest_v.png", 
       device = 'png', 
       plot = Signature_Cheng_vertical,
       width = 5, height = 18)
ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest_v.svg", 
       device = 'svg', 
       plot = Signature_Cheng_vertical,
       width = 5, height = 18)


ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest.pdf", 
       device = 'pdf', 
       plot = Signature_Cheng_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest.png", 
       device = 'png', 
       plot = Signature_Cheng_horizontal,
       width = 20, height = 5)
ggsave("data/Paper_Figure/Supplementary_Chengetal_expression_of_signatures_of_interest.svg", 
       device = 'svg', 
       plot = Signature_Cheng_horizontal,
       width = 20, height = 5)

(CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
  (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | p2 + ggtitle("UMAP plot"))

ggsave("data/Chengetal_expression_of_signatures_of_interest.pdf", device = 'pdf', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | p2 + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Chengetal_expression_of_signatures_of_interest.png", device = 'png', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | p2 + ggtitle("UMAP plot")),
       width = 18, height = 10)
ggsave("data/Chengetal_expression_of_signatures_of_interest.svg", device = 'svg', 
       plot = (CRPCsig51 + ggtitle("CRPCsig51") | CRPC_SCL + ggtitle("CRPC_SCL") | CRPC_NE + ggtitle("CRPC_NE")) / 
         (YAP_TAZ_pathway + ggtitle("YAP_TAZ_pathway") | YAP_TAZ_Targets + ggtitle("YAP_TAZ_Targets") | p2 + ggtitle("UMAP plot")),
       width = 18, height = 10)


Idents(samples.aggr) <- "cluster"
C12 = subset(samples.aggr, idents = "C12")
meta12 = C12@meta.data
meta12$ct = as.character(meta12$sample)
table(meta12$sample)
meta12$ct[grep("CRPC", meta12$sample)] = "CRPC"
meta12$ct[grep("PCa", meta12$sample)] = "PCa"
meta12$ct = as.factor(meta12$ct)
meta12$ct <- factor(meta12$ct, c("PCa", "CRPC"))
meta12 %>%
  group_by(ct) %>%
  tally() %>%
  ggplot(., aes(x="", y=n, fill=ct)) +
  geom_bar(stat="identity",  color="black") + #, position="fill") + 
  scale_fill_manual(values = c("#CC0033", "#3366FF")) +
  theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=12, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
        axis.title.y = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(size = 0),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position = "right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "C12", y = "Number of cells")


#VlnPlot(C12, features = c("YAP1", "EZH1", "EZH2"), split.by = "ct")
pfy = FeaturePlot(C12, features = c("YAP1"), 
                  min.cutoff = 'q5', max.cutoff = 'q95', 
                  order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez1 = FeaturePlot(C12, features = c("EZH1"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez2 = FeaturePlot(C12, features = c("EZH2"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 


ggsave("data/Rplot/gene_expression.svg", plot = pfy | pfez1 | pfez2,
       device = "svg", 
       width = 12, height = 5, bg = "transparent")


pfy = FeaturePlot(samples.aggr, features = c("YAP1"), 
                  min.cutoff = 'q5', max.cutoff = 'q95', 
                  order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez1 = FeaturePlot(samples.aggr, features = c("EZH1"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 

pfez2 = FeaturePlot(samples.aggr, features = c("EZH2"), 
                    min.cutoff = 'q5', max.cutoff = 'q95', 
                    order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) +
  #xlim(c(-7,-4)) + ylim(c(-10,-3))&
  #  xlim(c(-7.5,3.5)) + ylim(c(-10,0))  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) 
pfy | pfez1 | pfez2

ggsave("data/Rplot/gene_expression_app.svg", plot = pfy | pfez1 | pfez2,
       device = "svg", 
       width = 20, height = 8, bg = "transparent")


C12$ct = meta12$ct
DimPlot(C12, group.by = "ct", pt.size = 3) +
  xlim(c(-7,-4)) + ylim(c(-10,-3)) +
  scale_color_manual(values = c("#CC0033", "#3366FF"))
C12_2 <- RunUMAP(C12, dims = 1:20)

u = DimPlot(C12_2, group.by = "sample", pt.size = 2)

head(C12_2@meta.data)
pf = FeaturePlot(C12_2, features = c("YAP1"), 
                 min.cutoff = 'q10', max.cutoff = 'q50', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf + u

VlnPlot(C12_2, "YAP1", group.by = "sample", pt.size = 2)
VlnPlot(C12_2, "EZH2", group.by = "sample", pt.size = 2)

pf = FeaturePlot(C12_2, features = c("EZH2"), 
                 #min.cutoff = 'q10', max.cutoff = 'q90', 
                 order = T, pt.size = 2, cols = c("#FFFFCC", "#FF3300")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
pf + u

C12.scaled = ScaleData(C12, c("EZH1", "EZH2","YAP1", as.character(unlist(Signatures))))
hm = DoHeatmap(C12.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1","YAP1", as.character(unlist(Signatures))),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm
library('pheatmap')
C12_data = GetAssay(C12)
metaC12 = C12@meta.data
annotation_column = as.data.frame(metaC12[,c("lineage_subtypes", "ID", "cluster")])
head(annotation_column)
myb = seq(-1.5,1.5,by = 0.01)
myc = colorRampPalette(c('blue','white','red'))(length(myb))
HPv <- pheatmap::pheatmap(C12_data[c("YAP1",as.character(unlist(Signatures)))],
                          scale = 'row',
                          color = myc, 
                          breaks = myb,
                          #annotation_col = annotation_column,
                          #annotation_colors = ann_colors, 
                          cluster_rows = F, 
                          cluster_cols = T, 
                          cutree_cols  = 4,
                          cutree_row  = 2,
                          show_rownames = T,
                          show_colnames = F,
                          border_color = NA,
                          cellwidth=0.5, cellheight=10,
                          fontsize = 10, 
                          gaps_col = 34,
                          fontsize_row = 10, fontsize_col = 10, 
                          display_numbers = F,
                          #col=colors_hm,
                          filename = paste('HM_CRPCsig51.pdf',sep=''))
HPv
dev.off()

Idents(C12) <- "ct"
C12_magic13 <- subset(C12, idents =  "PCa")
C12_magic13.scaled = ScaleData(C12_magic13, c("EZH1", "EZH2","YAP1", as.character(unlist(Signatures))))
hm = DoHeatmap(C12_magic13.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1","YAP1", as.character(unlist(Signatures))),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm = DoHeatmap(C12_magic13, slot = "data",
               size = 8,
               features = c("YAP1","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


hm_C12_magic13_genes = DoHeatmap(C12_magic13, slot = "data",
               size = 8, label = F, 
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = F,
               raster = FALSE,
               group.bar.height = 0.05,
               angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm_C12_magic13_genes
ggsave("data/Chengetal_C12_13PCa_heatmap_of_gene_of_interest.png", device = 'png', 
       plot =hm_C12_magic13_genes,
       width = 10, height = 5)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_gene_of_interest.pdf", device = 'pdf', 
       plot =hm_C12_magic13_genes,
       width = 10, height = 5)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_gene_of_interest.svg", device = 'svg', 
       plot =hm_C12_magic13_genes,
       width = 10, height = 5)

hm1 = DoHeatmap(C12_magic13.scaled, 
               size = 8, label = F, 
               slot = 'data',
               #slot = 'scale.data', 
               features = c("YAP1", Signatures$CRPC_SCL),
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm2 = DoHeatmap(C12_magic13.scaled, 
               size = 8, label = F, 
               slot = 'data',
               #slot = 'scale.data', 
               features = c("YAP1", Signatures$CRPC_NE),
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm3 = DoHeatmap(C12_magic13.scaled, 
               size = 8, label = F, 
               slot = 'data',
               #slot = 'scale.data', 
               features = c("YAP1", Signatures$CRPCsig51),
               group.colors = col,
               #disp.min = -1.5,
               #disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm1 + hm2 + hm3

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8, slot = 'scale.data', 
               features = c("EZH1", "EZH2", "YAP1","YAP1", Signatures$CRPC_NE),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8, slot = 'scale.data', 
               features = c("EZH1", "EZH2", "YAP1","YAP1", Signatures$CRPC_NE),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8, slot = 'scale.data', 
               features = c("EZH1", "EZH2", "YAP1","YAP1", Signatures$NEPC_dge),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8, slot = 'scale.data', 
               features = c("EZH1", "EZH2", "YAP1","YAP1", subset_NEmarkers),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm


hm = DoHeatmap(C12_magic13.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1","YAP1", Signatures$YAP_TAZ_Targets),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1","YAP1", Signatures$YAP_TAZ_pathway),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1",subset_NEmarkers),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

hm = DoHeatmap(C12_magic13.scaled, 
               size = 8,
               features = c("EZH1", "EZH2", "YAP1",NEPC_dge$X1),
               group.colors = col,
               disp.min = -1.5,
               disp.max = 1.5,
               draw.lines = F,
               group.bar.height = 0.05,
               angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)
hm

head(C12)
table(C12$sample)
Idents(C12) <- "ct"
hm_C12 = DoHeatmap(C12, 
               slot = "data",
               size = 8, label = TRUE, 
               features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
               #group.colors = col,
               #disp.min = -1.5,
               disp.max = 3,
               draw.lines = F,
               group.bar.height = 0.05,
               raster = FALSE,
               angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_C12
ggsave("data/Chengetal_C12_heatmap_of_genes_of_interest.pdf", device = 'pdf', 
       plot = hm_C12,
       width = 10, height = 6)
ggsave("data/Chengetal_C12_heatmap_of_genes_of_interest.png", device = 'png', 
       plot = hm_C12,
       width = 10,  height = 6)
ggsave("data/Chengetal_C12_heatmap_of_genes_of_interest.svg", device = 'svg', 
       plot = hm_C12,
       width = 18,  height = 6)


hm_C12_f = DoHeatmap(C12, 
                   slot = "data",
                   size = 6, label = TRUE, 
                   features = c("YAP1","EZH1","EZH2","SYP", "NDRG1"),
                   group.colors = c("darkgoldenrod1", "cornflowerblue"),
                   disp.max = 3,
                   draw.lines = T,
                   group.bar.height = 0.1,
                   raster = FALSE,
                   angle = 45) & scale_fill_gradient2(low = 'white', mid = "white", high = "firebrick2", midpoint = 0)

hm_C12_f + theme(text = element_text(size = 20))
hm_C12_f + theme(axis.text.y = element_text(size = 24))
ggsave("data/Paper_Figure/Supplementaty_Chengetal_C12_heatmap_of_genes_of_interest.pdf", device = 'pdf', 
       plot = hm_C12_f + theme(axis.text.y = element_text(size = 24)),
       width = 6, height = 5)
ggsave("data/Paper_Figure/Supplementaty_Chengetal_C12_heatmap_of_genes_of_interest.png", device = 'png', 
       plot = hm_C12_f + theme(axis.text.y = element_text(size = 24)),
       width = 6,  height = 5)
ggsave("data/Paper_Figure/Supplementaty_Chengetal_C12_heatmap_of_genes_of_interest.svg", device = 'svg', 
       plot = hm_C12_f + theme(text = element_text(size = 24)),
       width = 6,  height = 5)


hm1 = DoHeatmap(C12, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPC_SCL),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                raster = FALSE,
                angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm2 = DoHeatmap(C12, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPC_NE),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                raster = FALSE,
                angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm3 = DoHeatmap(C12, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPCsig51),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                raster = FALSE,
                angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm4 = DoHeatmap(C12, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$YAP_TAZ_pathway),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                raster = FALSE,
                angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm5 = DoHeatmap(C12, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$YAP_TAZ_Targets),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                raster = FALSE,
                group.bar.height = 0.05,
                angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)



hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE")

ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_1.pdf", device = 'pdf', 
       plot = hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE"),
       width = 18, height = 14)
ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_1.png", device = 'png', 
       plot = hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE"),
       width = 18, height = 14)
ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_1.svg", device = 'svg', 
       plot = hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE"),
       width = 18, height = 14)

hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets")

ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_2.pdf", device = 'pdf', 
       plot = hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)
ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_2.png", device = 'png', 
       plot = hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)
ggsave("data/Chengetal_C12_heatmap_of_signature_of_interest_2.svg", device = 'svg', 
       plot = hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 15, height = 5)

YAPp <- subset(C12, subset = YAP1 > 0) 
DimPlot(YAPp)

hm1 = DoHeatmap(subset(C12, subset = YAP1 > 0), 
                size = 8, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
                group.colors = col,
                #disp.min = -1.5,
                disp.max = 1.5,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm2 = DoHeatmap(C12, 
                size = 8, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPC_NE),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm3 = DoHeatmap(C12, 
                size = 8, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPCsig51),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 90) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)



hm_C12_m13 = DoHeatmap(C12_magic13, 
                   slot = "data",
                   size = 8, label = FALSE, 
                   features = c("YAP1","ITGA2","TNC","EZH1","EZH2","SYP", "NDRG1", "CHGA", "NPTX1"),
                   #group.colors = col,
                   #disp.min = -1.5,
                   disp.max = 3,
                   draw.lines = F,
                   group.bar.height = 0.05,
                   angle = 45) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm_C12_m13
ggsave("data/Chengetal_C12_13PCa_heatmap_of_genes_of_interest.pdf", device = 'pdf', 
       plot = hm_C12_m13,
       width = 10, height = 5)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_genes_of_interest.png", device = 'png', 
       plot = hm_C12_m13,
       width = 10,  height = 5)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_genes_of_interest.svg", device = 'svg', 
       plot = hm_C12_m13,
       width = 10,  height = 5)


hm1 = DoHeatmap(C12_magic13, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPC_SCL),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm2 = DoHeatmap(C12_magic13, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPC_NE),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm3 = DoHeatmap(C12_magic13, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$CRPCsig51),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm4 = DoHeatmap(C12_magic13, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$YAP_TAZ_pathway),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)

hm5 = DoHeatmap(C12_magic13, 
                size = 4, label = T, 
                slot = 'data',
                #slot = 'scale.data', 
                features = c("YAP1", Signatures$YAP_TAZ_Targets),
                group.colors = col,
                #disp.min = -1.5,
                #disp.max = 1.5,
                disp.max = 3,
                draw.lines = F,
                group.bar.height = 0.05,
                angle = 0) & scale_fill_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0)



hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE") | hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets")

ggsave("data/Chengetal_C12_13PCa_heatmap_of_signatures_of_interest.pdf", device = 'pdf', 
       plot = hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE") | hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 18, height = 14)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_signature_of_interest_1.png", device = 'png', 
       plot =  hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE") | hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 18, height = 14)
ggsave("data/Chengetal_C12_13PCa_heatmap_of_signature_of_interest_1.svg", device = 'svg', 
       plot =  hm3 + ggtitle ("CRPCsig51") | hm1 + ggtitle ("CRPC_SCL")  | hm2 + ggtitle ("CRPC_NE") | hm4 + ggtitle ("YAP_TAZ_Pathway") | hm5 + ggtitle ("YAP_TAZ_Targets"),
       width = 18, height = 14)

