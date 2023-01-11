#Computing senescence score (with ssGSEA) for Johnson et al, 2021 dataset
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
MetaModule_gene_sets <- read.delim("./Inputs/Neftel_SingleCellPortal/scAdultPediatricGlio/IDHwt.GBM.MetaModules.tsv", 
                                   stringsAsFactors = FALSE)

#Import signature genes
biomart_murin_humain <- read.table("./Inputs/201017_p16Hi_genelist_murine_GBM.tsv", 
                                   header = T, sep = "\t")

# load Matrix
mat <- data.table::fread("Inputs/Johnson_nature_2021/analysis_scRNAseq_tumor_gene_expression.tsv")

#retrieve info at the end of the expression matrix (??)
mat_stat <- t(data.frame(tail(mat[,-1], 5), row.names = tail(mat[,1][[1]], 5), check.names = FALSE))
mat <- mat[1:(nrow(mat)-5),] #remove these lines from the expression matrix

#single cell metadata
meta <- read.table("Inputs/Johnson_nature_2021/analysis_scRNAseq_tumor_metadata.tsv", 
                   header = T, sep = "\t", as.is = T, row.names = 1)

metadata <- merge(meta, mat_stat, by = "row.names", sort = FALSE)

#tumour metadata
clinical_data <- read.delim("Inputs/Johnson_nature_2021/clinical_metadata.tsv", h = T)
meta2 <- merge(metadata[,-1], clinical_data[,-c(1:2)], by = "case_barcode", sort = FALSE)
rownames(meta2) <- metadata$Row.names

#load genes
genes <- read.delim("Inputs/Johnson_nature_2021/ref_genes.tsv", header = TRUE)

#some of the genes in the expression matrix are not in reference gene table
notinGenes <- mat[!mat[,1][[1]] %in% genes$gene_id,1][[1]]
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = "37") ## genes annotated in hg19
notinGenes_annotated <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                              filters = 'ensembl_gene_id',
                              values = notinGenes,
                              mart = ensembl)
colnames(notinGenes_annotated) <- c("gene_id", "gene_name")
#merge results
genes <- dplyr::bind_rows(genes, notinGenes_annotated)
rownames(genes) <- genes$gene_id #add rownames
genes <- genes[mat[,1][[1]],] #order lines based on expression matrix

genes$gene_name <- ifelse(is.na(genes$gene_name), genes$gene_id, genes$gene_name) #add ensemblGeneID as gene name if missing

mat <- data.frame(mat[,-1], row.names = make.unique(genes$gene_name), check.names = FALSE)
verhaak <- CreateSeuratObject(counts = mat, project = "Verhaak", meta.data = meta2)

# Add umap Coordinates
umap_coord_mat <- as(subset(meta, select = c("umap_1", "umap_2")), "matrix")

verhaak[['umap']] <- CreateDimReducObject(embeddings = umap_coord_mat, key = "UMAP_", global = T, assay = "RNA")


#retrieve only GBM
verhaak <- subset(verhaak, histological_classification == "Glioblastoma")

#set Idents as cell state identified by Johnson and colleagues
Idents(verhaak) = verhaak@meta.data$cell_state

#remove cell transcriptomes that expressed PTPRC.
verhaak <- subset(verhaak, PTPRC == 0)
verhaak <- subset(verhaak, idents =  c("Myeloid", "Granulocyte", "T cell", "Dendritic cell", "B cell"), invert = TRUE)


VlnPlot(verhaak, features = c("ENSGGENES", "ENSGUMI"))
UMAPPlot(verhaak,label = T)

#Compute Metamodule
verhaak[["percent.mt"]] <- PercentageFeatureSet(verhaak, pattern = "^MT-")
verhaak <- MetaModule_Assignment(verhaak, MetaModule_gene_sets, assay = "RNA")
UMAPPlot(verhaak, split.by = c("case_barcode"), group.by = c("state"))

#Merge MES and NPC 
verhaak@meta.data$state <- gsub("[1-2]$",
                                "",
                                verhaak@meta.data$state)

## compute signature score
list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Human.gene.name))
verhaak_ssGsea = gsva(as.matrix(verhaak@assays$RNA@data), list_geneset, 
                      method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

verhaak@meta.data$score_ssGSEA <- verhaak_ssGsea[1,] * 100
verhaak@meta.data$rate_ssGSEA <- dplyr::ntile(verhaak@meta.data$score_ssGSEA, 10)
verhaak@meta.data$rate_ssGSEA <- ifelse(verhaak@meta.data$rate_ssGSEA == 1,
                                        "Low",
                                        ifelse(verhaak@meta.data$rate_ssGSEA == 10,
                                               "High",
                                               "Medium"))

#plots
plot1 <- FeaturePlot(verhaak, features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", brewer.pal(n = 9, name = "Reds")))
plot2 <- UMAPPlot(verhaak, group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

UMAPPlot(verhaak, group.by = c("case_barcode"), split.by = "rate_ssGSEA")

#violin plot
ggplot(verhaak@meta.data, 
       aes(x = case_barcode,
           y = score_ssGSEA,
           fill = case_barcode,
           colour = case_barcode)) +
  geom_violin(alpha=0.7, scale="width",adjust = .5) + 
  geom_point(size=0.3,alpha = 0.8, position = "jitter") +
  ggpubr::theme_pubr() +
  geom_hline(yintercept = max(verhaak@meta.data$score_ssGSEA[verhaak@meta.data$rate_ssGSEA == "Low"]),
             linetype="dashed") +
  geom_hline(yintercept = min(verhaak@meta.data$score_ssGSEA[verhaak@meta.data$rate_ssGSEA == "High"]),
             linetype="dashed") +
  ggtitle("Senescence score (ssGSEA) distribution between non-integrated patient datasets", "Badhuri dataset")

#pie charts
verhaak_summary <- subset(verhaak@meta.data, select = c("rate_ssGSEA", "state", "case_barcode"))
# by Module
verhaak_summary_byRate <- dcast(verhaak_summary, rate_ssGSEA ~ state)
verhaak_summary_byRate <- melt(verhaak_summary_byRate, id.vars = "rate_ssGSEA")
verhaak_summary_byRate <- verhaak_summary_byRate %>% group_by(rate_ssGSEA) %>% dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(verhaak_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ rate_ssGSEA") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

#return Seurat Object
saveRDS(verhaak, file = "./figures_johnson_nature_2021/verhaak_SeuratObject.rds")
