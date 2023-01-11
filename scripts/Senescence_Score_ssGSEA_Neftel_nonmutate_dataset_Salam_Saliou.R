#Computing senescence score (with ssGSEA) for Neftel et al dataset 
#Only non mutated tumors for CD4 and CDKN2A
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

biomart_murin_humain <- read.table("./Inputs/201017_p16Hi_genelist_murine_GBM.tsv", 
                                   header = T, sep = "\t")

## Read smartseq log10TPM
neftel_smartseq <- read.delim("./Inputs/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", row.names = 1)

#metadata SmartSeq2 SingleCellPortal
neftel_smarseq_metadata <- read.delim("./Inputs/Neftel_SingleCellPortal/scAdultPediatricGlio/IDHwt.GBM.Metadata.SS2.txt", check.names = FALSE, header = TRUE, row.names = 1)

#filter dataframe based on several variables
sample_PTPRC_neg <- colnames(neftel_smartseq)[neftel_smartseq["PTPRC",] == 0]
tumor_to_keep <- c("MGH66" , "MGH100", "MGH101", "MGH102", "MGH104", "MGH105",
                   "MGH106", "MGH110", "MGH113", "MGH115", "MGH121", "MGH122",
                   "MGH124", "MGH125", "MGH126", "MGH129") #not mutate for CD4 and/or CDKN2A

cat("Tumors kept : ", tumor_to_keep)

neftel_smarseq_metadata <- subset(neftel_smarseq_metadata, GBMType == "Adult" &
                                    Sample %in% tumor_to_keep & 
                                    rownames(neftel_smarseq_metadata) %in% sample_PTPRC_neg)

#Create Seurat objects
##smartseq2
neftel_smartseq_seurat <- CreateSeuratObject(neftel_smartseq[,rownames(neftel_smarseq_metadata)], project = "SmartSeq2", meta.data = neftel_smarseq_metadata, names.delim = "_")


#Basic Pipeline
neftel_smartseq_seurat <- FindVariableFeatures(neftel_smartseq_seurat)
neftel_smartseq_seurat <- ScaleData(neftel_smartseq_seurat)
neftel_smartseq_seurat <- RunPCA(neftel_smartseq_seurat, verbose = FALSE)
neftel_smartseq_seurat <- RunUMAP(neftel_smartseq_seurat, 
                                  reduction = "pca", 
                                  verbose = FALSE, 
                                  dims = 1:30)

## QC
VlnPlot(neftel_smartseq_seurat, features = c("nCount_RNA", "nFeature_RNA")) # QC plots

## Representation
DimPlot(neftel_smartseq_seurat, reduction = "umap", group.by = c("Sample"))


neftel_smartseq_seurat <- MetaModule_Assignment(neftel_smartseq_seurat, MetaModule_gene_sets, assay = "RNA")
DimPlot(neftel_smartseq_seurat, split.by = c("Sample"), group.by = c("state"))

#Merge MES and NPC 
neftel_smartseq_seurat@meta.data$state <- gsub("[1-2]$",
                                               "",
                                               neftel_smartseq_seurat@meta.data$state)

#tumor color palette
neftel_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(neftel_smartseq_seurat$Sample)))
names(neftel_colors) <- unique(neftel_smartseq_seurat$Sample)

## import 10X files
neftel_10X <- Read10X(data.dir = "./Inputs/Neftel_SingleCellPortal/scAdultPediatricGlio/10X_count1")

#metadata GEO
neftel_all_metadata <- read.delim("./Inputs/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.csv", header = TRUE, row.names = 1, sep = ";", check.names = F)

neftel_all_metadata <- subset(neftel_all_metadata, `adult/pediatric` == "adult" &
                                `tumour name` %in% tumor_to_keep)

neftel_10X_metadata <- subset(neftel_all_metadata, rownames(neftel_all_metadata) %in% colnames(neftel_10X)[neftel_10X["PTPRC",] == 0])

#Create Seurat object
##10X count_1
neftel_10X_seurat <- CreateSeuratObject(neftel_10X[,rownames(neftel_10X_metadata)], project = "10X", meta.data = neftel_10X_metadata, names.delim = "-")

## Basic pipeline
neftel_10X_seurat <- SCTransform(neftel_10X_seurat,
                                 verbose = FALSE, 
                                 bin_size = 250, 
                                 verbosity = 0)

neftel_10X_seurat <- RunPCA(neftel_10X_seurat, verbose = FALSE)
neftel_10X_seurat <- RunUMAP(neftel_10X_seurat, 
                             reduction = "pca", 
                             verbose = FALSE, 
                             dims = 1:30)

## Representation
DimPlot(neftel_10X_seurat, reduction = "umap", group.by = c("tumour.name"))


neftel_10X_seurat <- MetaModule_Assignment(neftel_10X_seurat, MetaModule_gene_sets, assay = "RNA")
DimPlot(neftel_10X_seurat, split.by = c("tumour.name"), group.by = c("state"))

#Merge MES and NPC 
neftel_10X_seurat@meta.data$state <- gsub("[1-2]$",
                                          "",
                                          neftel_10X_seurat@meta.data$state)

## colors
neftel_colors_10X <- neftel_colors[names(neftel_colors) %in% unique(neftel_10X_seurat$tumour.name)]
neftel_colors_10X <- c(neftel_colors_10X, MGH126 = "grey")


## compute signature score for SS2 dataset
list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Human.gene.name))
neftel_SS2_ssGsea = gsva(as.matrix(neftel_smartseq_seurat@assays$RNA@data), 
                         list_geneset, 
                         method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

neftel_smartseq_seurat@meta.data$score_ssGSEA <- neftel_SS2_ssGsea[1,] * 100
neftel_smartseq_seurat@meta.data$rate_ssGSEA <- dplyr::ntile(neftel_smartseq_seurat@meta.data$score_ssGSEA, 10)
neftel_smartseq_seurat@meta.data$rate_ssGSEA <- ifelse(neftel_smartseq_seurat@meta.data$rate_ssGSEA == 1,
                                                       "Low",
                                                       ifelse(neftel_smartseq_seurat@meta.data$rate_ssGSEA == 10,
                                                              "High",
                                                              "Medium"))

#plots
plot1 <- FeaturePlot(neftel_smartseq_seurat, 
                     features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", 
                                     brewer.pal(n = 9, name = "Reds")))

plot2 <- UMAPPlot(neftel_smartseq_seurat, 
                  group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

UMAPPlot(neftel_smartseq_seurat, 
         group.by = c("Sample"), 
         split.by = "rate_ssGSEA", 
         cols = neftel_colors)

## export results
write.table(neftel_SS2_ssGsea, 
            "./figures_neftel_10X_SS2/neftel_SS2_ssGsea_controle.tsv", 
            quote = FALSE, 
            sep = "\t")

#violin plot
ggplot(neftel_smartseq_seurat@meta.data, 
       aes(x = Sample,
           y = score_ssGSEA,
           fill = Sample,
           colour = Sample)) +
  geom_violin(alpha=0.7, scale="width",adjust = .5) + 
  geom_point(size=0.3,alpha = 0.8, position = "jitter")+
  scale_fill_manual(values = neftel_colors) +
  ggpubr::theme_pubr() +
  scale_color_manual(values = neftel_colors) +
  geom_hline(yintercept = max(neftel_smartseq_seurat@meta.data$score_ssGSEA[neftel_smartseq_seurat@meta.data$rate_ssGSEA == "Low"]),
             linetype="dashed") +
  geom_hline(yintercept = min(neftel_smartseq_seurat@meta.data$score_ssGSEA[neftel_smartseq_seurat@meta.data$rate_ssGSEA == "High"]),
             linetype="dashed") +
  ggtitle("Senescence score (ssGSEA) distribution between non-integrated patient datasets", 
          "SS2 dataset")

#pie charts
neftel_SS2_summary <- subset(neftel_smartseq_seurat@meta.data, 
                             select = c("rate_ssGSEA", "state", "Sample"))

# by Module
neftel_SS2_summary_byRate <- dcast(neftel_SS2_summary, rate_ssGSEA ~ state)
neftel_SS2_summary_byRate <- melt(neftel_SS2_summary_byRate, id.vars = "rate_ssGSEA")

neftel_SS2_summary_byRate <- neftel_SS2_summary_byRate %>% 
  group_by(rate_ssGSEA) %>% 
  dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(neftel_SS2_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ rate_ssGSEA") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), 
             position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

#return Seurat Object
saveRDS(neftel_smartseq_seurat, 
        file = "./figures_neftel_10X_SS2/neftel_smartseq_notmutate_SeuratObject.rds")

## compute signature score for 10X dataset
list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Human.gene.name))
neftel_10X_ssGsea = gsva(as.matrix(neftel_10X_seurat@assays$RNA@data), list_geneset, 
                         method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

neftel_10X_seurat@meta.data$score_ssGSEA <- neftel_10X_ssGsea[1,] * 100
neftel_10X_seurat@meta.data$rate_ssGSEA <- dplyr::ntile(neftel_10X_seurat@meta.data$score_ssGSEA, 10)
neftel_10X_seurat@meta.data$rate_ssGSEA <- ifelse(neftel_10X_seurat@meta.data$rate_ssGSEA == 1,
                                                  "Low",
                                                  ifelse(neftel_10X_seurat@meta.data$rate_ssGSEA == 10,
                                                         "High",
                                                         "Medium"))

#plots
plot1 <- FeaturePlot(neftel_10X_seurat, 
                     features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", 
                                     brewer.pal(n = 9, name = "Reds")))

plot2 <- UMAPPlot(neftel_10X_seurat, 
                  group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

UMAPPlot(neftel_10X_seurat, 
         group.by = c("tumour.name"), 
         split.by = "rate_ssGSEA", 
         cols = neftel_colors_10X)

## export results
write.table(neftel_10X_ssGsea, 
            "./figures_neftel_10X_SS2/neftel_10X_ssGsea_controle.tsv", 
            quote = FALSE, 
            sep = "\t")
