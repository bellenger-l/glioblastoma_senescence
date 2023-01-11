#Computing senescence score (with ssGSEA) for Badhuri et al dataset
# Salam, Saliou et al, 2022 Nature Comm

#Import library
library(Matrix)
library(Seurat)
library(lava)
library(ggplot2)
library(readxl)
library(stringr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(scalop)
library(ggExtra)
library(gghighlight)
library(dplyr)
library(psych)
library(gridExtra)
library(GSVA)
library(biomaRt)

#Custom function
MetaModule_Assignment <- function(SeuratObject, SignatureLists, assay){
  
  ## Neftel methods and quadrant representation
  #sigScores
  df <- as.matrix(SeuratObject@assays[[assay]]@data)
  module_scores <- sigScores(m = df,
                             sigs = SignatureLists)
  module_scores$state <- colnames(module_scores)[max.col(module_scores[,1:6])]
  
  #hierarchy plot
  quadrants = list("AClike",
                   c("MESlike1", "MESlike2"), 
                   "OPClike", 
                   c("NPClike1", "NPClike2"))
  
  dat = as.data.frame(sapply(quadrants, function(col) { apply(subset(module_scores, select = col), 
                                                              1, 
                                                              max) 
  }))
  
  rows = rownames(module_scores)
  colnames(dat) = c('bl', 'br', 'tl', 'tr')
  
  dat = dplyr::mutate(dat,
                      bottom = pmax(bl, br),
                      top = pmax(tl, tr),
                      b.center = br - bl,
                      t.center = tr - tl,
                      x = ifelse(bottom > top, b.center, t.center), # dependent var
                      x.scaled = (sign(x) * log2(abs(x) + 1)),
                      y = top - bottom, # independent var
                      y.scaled = (sign(y) * log2(abs(y) + 1)))
  #log scale
  dat = dplyr::transmute(dat, X = x.scaled, Y = y.scaled)
  rownames(dat) = rows
  class(dat) = append(class(dat), 'hierarchy')
  
  #hierarchy plot
  hierarchy_dat <- merge(dat, module_scores, by = "row.names")
  rownames(hierarchy_dat) <- hierarchy_dat$Row.names
  hierarchy_dat <- hierarchy_dat[,-1]
  state <- rownames(hierarchy_dat)
  names(state) <- hierarchy_dat$state
  #plot_hierarchy
  print(plot_hierarchy(hierarchy_dat, 
                       quadrant.names = c("AC", "MES", "OPC", "NPC"), 
                       groups = state, legend.horiz = F, legend.pos = "left", 
                       main = "Hierarchy plot with Neftel signature method"))
  
  #add neftel signature info
  SeuratObject@meta.data <- merge(SeuratObject@meta.data, hierarchy_dat, by = "row.names")
  rownames(SeuratObject@meta.data) <- SeuratObject@meta.data$Row.names
  SeuratObject@meta.data <- SeuratObject@meta.data[,-1]
  
  return(SeuratObject)
}

##import modules info
MetaModule_gene_sets <- read.delim("./Inputs/Neftel_SingleCellPortal/scAdultPediatricGlio/IDHwt.GBM.MetaModules.tsv", stringsAsFactors = FALSE)


biomart_murin_humain <- read.table("./Inputs/201017_p16Hi_genelist_murine_GBM.tsv", header = T, sep = "\t")

# load Matrix
mat <- data.table::fread("Inputs/Badhuri/exprMatrix.tsv.gz")
meta <- read.table("Inputs/Badhuri/meta(1).tsv", header=T, sep="\t", as.is=T, row.names=1)
genes <- mat[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[,-1], row.names=genes)
badhuri <- CreateSeuratObject(counts = mat, project = "Badhuri", meta.data=meta)

# Add tsne Coordinates
tsne_coord <- data.table::fread("Inputs/Badhuri/Seurat_tMinusSNE.coords.tsv.gz")
colnames(tsne_coord) = c("barcodes", "tsne_1", "tsne_2")
tsne_coord <- data.frame(tsne_coord)
rownames(tsne_coord) = gsub("-",".", tsne_coord$barcodes)
tsne_coord <- tsne_coord[,-1]

tsne_coord$tsne_1 <- as.numeric(tsne_coord$tsne_1, length = 4)
tsne_coord$tsne_2 <- as.numeric(tsne_coord$tsne_2, length = 4)

tsne_coord_mat <- as(tsne_coord, "matrix")

badhuri[['tsne']] <- CreateDimReducObject(embeddings = tsne_coord_mat, key = "tSNE_", global = T, assay = "RNA")

# Add metadata
badhuri@meta.data = badhuri@meta.data[,1:2]
rownames(badhuri@meta.data) = gsub("-",".",rownames(badhuri@meta.data))

metadata = openxlsx::read.xlsx("Inputs/Badhuri/SupplementalTable3.xlsx")
metadata$Cell = gsub("-",".",metadata$Cell)
rownames(metadata) = metadata$Cell

badhuri = AddMetaData(badhuri,metadata = metadata)
Idents(badhuri) = badhuri@meta.data$Cell.Type.Assignment

#remove cell transcriptomes that expressed PTPRC.
badhuri <- subset(badhuri, PTPRC == 0)
badhuri <- subset(badhuri, idents =  c("Tumor Associated Macrophage", "Microglia", "B Cells", "Dividing B Cells"), invert = TRUE)

#Keep 10X samples
sample_10X <- c("SF11215", "SF11209", "SF11232", "SF11159", "SF11285", "SF11247")
badhuri <- subset(badhuri,  Tumor.ID %in% sample_10X)

DimPlot(badhuri,label = T)

#Compute Metamodule
badhuri <- MetaModule_Assignment(badhuri, MetaModule_gene_sets, assay = "RNA")
DimPlot(badhuri, split.by = c("Tumor.ID"), group.by = c("state"))

#Merge MES and NPC 
badhuri@meta.data$state <- gsub("[1-2]$",
                                "",
                                badhuri@meta.data$state)

## compute signature score

list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Human.gene.name))
badhuri_ssGsea = gsva(as.matrix(badhuri@assays$RNA@data), list_geneset, 
                      method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

badhuri@meta.data$score_ssGSEA <- badhuri_ssGsea[1,] * 100
badhuri@meta.data$rate_ssGSEA <- dplyr::ntile(badhuri@meta.data$score_ssGSEA, 10)
badhuri@meta.data$rate_ssGSEA <- ifelse(badhuri@meta.data$rate_ssGSEA == 1,
                                        "Low",
                                        ifelse(badhuri@meta.data$rate_ssGSEA == 10,
                                               "High",
                                               "Medium"))

#plots
plot1 <- FeaturePlot(badhuri, features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", brewer.pal(n = 9, name = "Reds")))
plot2 <- DimPlot(badhuri, group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

DimPlot(badhuri, group.by = c("Tumor.ID"), split.by = "rate_ssGSEA")

## export results
write.table(badhuri_ssGsea, "./figures_badhuri/badhuri_ssGsea.tsv", quote = FALSE, sep = "\t")

#violin plot
ggplot(badhuri@meta.data, 
       aes(x = Tumor.ID,
           y = score_ssGSEA,
           fill = Tumor.ID,
           colour = Tumor.ID)) +
  geom_violin(alpha=0.7, scale="width",adjust = .5) + 
  geom_point(size=0.3,alpha = 0.8, position = "jitter") +
  ggpubr::theme_pubr() +
  geom_hline(yintercept = max(badhuri@meta.data$score_ssGSEA[badhuri@meta.data$rate_ssGSEA == "Low"]),
             linetype="dashed") +
  geom_hline(yintercept = min(badhuri@meta.data$score_ssGSEA[badhuri@meta.data$rate_ssGSEA == "High"]),
             linetype="dashed") +
  ggtitle("Senescence score (ssGSEA) distribution between non-integrated patient datasets", "Badhuri dataset")

#pie charts
badhuri_summary <- subset(badhuri@meta.data, select = c("rate_ssGSEA", "state", "Tumor.ID"))
# by Module
badhuri_summary_byRate <- dcast(badhuri_summary, rate_ssGSEA ~ state)
badhuri_summary_byRate <- melt(badhuri_summary_byRate, id.vars = "rate_ssGSEA")
badhuri_summary_byRate <- badhuri_summary_byRate %>% group_by(rate_ssGSEA) %>% dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(badhuri_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ rate_ssGSEA") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

#return Seurat Object
saveRDS(badhuri, file = "./figures_badhuri/badhuri_SeuratObject.rds")