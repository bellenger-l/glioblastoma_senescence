#Computing senescence score (with ssGSEA) for mouse dataset 
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

#color palette
cluster_colors <- c("NP-like"="#666666", "pri-OPC-like_1"="#AD7700",
                    "pri-OPC-like_2"= "#1C91D4", "pri-OPC-like_3"="#007756",
                    "pri-OPC-like_4"= "#D5C711", "G2/M_2"= "#005685",
                    "G1/S_1"="#A04700", "G1/S_2" ="#B14380",
                    "hypoxic cell"= "#4D4D4D", "astrocyte" ="#FFBE2D",
                    "COP"= "#80C7EF", "repressed cell" ="#00F6B3",
                    "mOL"= "#F4EB71","G2/M_1"= "#06A5FF",
                    "low quality cell"= "#FF8320", 
                    "ECM-remodelling cell"="#E69F00",
                    "neuron" ="#56B4E9")

##import modules info
MetaModule_gene_sets <- read.delim("./Inputs/Neftel_SingleCellPortal/scAdultPediatricGlio/IDHwt.GBM.MetaModules.tsv", stringsAsFactors = FALSE)

#Change gene case to fit mouse annotation
mouse_MetaModule_gene_sets <- data.frame(apply(MetaModule_gene_sets, 
                                               2, 
                                               str_to_sentence),
                                         stringsAsFactors = F)

biomart_murin_humain <- read.table("./Inputs/201017_p16Hi_genelist_murine_GBM.tsv", 
                                   header = T, sep = "\t")

#import files
load("./Inputs/Seurat_tumor_wocluster15.RData")

#upload cluster annotation
cluster_annotation_table <- t(read.delim("./Inputs/201015_annotations_clusters_tumoraux_res0.6_cluster_names.tsv", 
                                         sep = "\t", 
                                         header = FALSE, 
                                         row.names = 1))

cluster_annotation <- as.vector(cluster_annotation_table[,2])
names(cluster_annotation) <- as.vector(cluster_annotation_table[,1])

##Rename cluster based on manual annotation
Idents(obj3) <- "SCT_snn_res.0.6" # set the cell identifications to Seurat clusters
obj3 <- RenameIdents(obj3, cluster_annotation)

obj3@meta.data$cluster_names <- Idents(obj3)


#Score analysis for control data
mouse_ctl <- obj3[,rownames(subset(obj3@meta.data, condition == "ctl"))] #select only control

mouse_ctl <- MetaModule_Assignment(mouse_ctl, mouse_MetaModule_gene_sets, assay = "SCT")

#plot
p <- ggplot(mouse_ctl@meta.data, aes(X, Y)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

p_cluster <- p +
  geom_point(aes(colour = cluster_names)) +
  scale_color_manual(values = cluster_colors)

p_cluster <- ggMarginal(p_cluster, type = "violin")

p_G2M <- p +
  geom_point(aes(color = G2.M), alpha = 0.5) + scale_color_gradient2(low = "darkgrey", mid = "lightgrey", high = "red") 

gridExtra::grid.arrange(p_cluster, p_G2M, ncol=2)

p_priOPC <- ggplot(mouse_ctl@meta.data, aes(X, Y, colour = factor(cluster_names))) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
  geom_point() + 
  scale_color_manual(values = cluster_colors) + 
  gghighlight("pri-OPC-like_1" == cluster_names | 
                "pri-OPC-like_2" == cluster_names |
                "pri-OPC-like_3" == cluster_names | 
                "pri-OPC-like_4" == cluster_names, 
              keep_scales = TRUE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

p_cycling <- ggplot(mouse_ctl@meta.data, aes(X, Y, colour = factor(cluster_names))) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
  geom_point() + 
  scale_color_manual(values = cluster_colors) + 
  gghighlight("G1/S_1" == cluster_names | 
                "G1/S_2" == cluster_names |
                "G2/M_1" == cluster_names | 
                "G2/M_2" == cluster_names, 
              keep_scales = TRUE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

gridExtra::grid.arrange(p_priOPC, p_cycling, ncol=2)

#Merge MES and NPC 
mouse_ctl@meta.data$state <- gsub("[1-2]$",
                                  "",
                                  mouse_ctl@meta.data$state)

## compute signature score

list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Gene.name))
mouse_ctl_ssGsea = gsva(as.matrix(mouse_ctl@assays$SCT@data), list_geneset, 
                        method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

mouse_ctl@meta.data$score_ssGSEA <- mouse_ctl_ssGsea[1,] * 100
mouse_ctl@meta.data$rate_ssGSEA <- dplyr::ntile(mouse_ctl@meta.data$score_ssGSEA, 10)
mouse_ctl@meta.data$rate_ssGSEA <- ifelse(mouse_ctl@meta.data$rate_ssGSEA == 1,
                                          "Low",
                                          ifelse(mouse_ctl@meta.data$rate_ssGSEA == 10,
                                                 "High",
                                                 "Medium"))

#plots
plot1 <- FeaturePlot(mouse_ctl, features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", 
                                     brewer.pal(n = 9, name = "Reds")))
plot2 <- UMAPPlot(mouse_ctl, 
                  group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

UMAPPlot(mouse_ctl, group.by = c("cluster_names"), 
         split.by = "rate_ssGSEA", 
         cols = cluster_colors)


## export results
write.table(mouse_ctl_ssGsea, 
            "./figures_mouse_ctrl_p16/mouse_ctl_ssGSEA.tsv", 
            quote = FALSE, sep = "\t")

ggplot(mouse_ctl@meta.data, 
       aes(x = cluster_names,
           y = score_ssGSEA,
           fill = cluster_names,
           colour = cluster_names)) +
  geom_violin(alpha=0.7, scale="width",adjust = .5) + 
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  geom_point(size=0.3,alpha = 0.8, position = "jitter") +
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  geom_hline(yintercept = max(mouse_ctl@meta.data$score_ssGSEA[mouse_ctl@meta.data$rate_ssGSEA == "Low"]),
             linetype="dashed") +
  geom_hline(yintercept = min(mouse_ctl@meta.data$score_ssGSEA[mouse_ctl@meta.data$rate_ssGSEA == "High"]),
             linetype="dashed") +
  ggtitle("Senescence score (ssGSEA) distribution between non-integrated mouse dataset", 
          "Ctrl transcriptomes")


mouse_ctl_summary <- subset(mouse_ctl@meta.data, 
                            select = c("rate_ssGSEA", "state", "cluster_names"))
# by Module
mouse_ctl_summary_byRate <- dcast(mouse_ctl_summary, 
                                  rate_ssGSEA ~ state)
mouse_ctl_summary_byRate <- melt(mouse_ctl_summary_byRate, 
                                 id.vars = "rate_ssGSEA")
mouse_ctl_summary_byRate <- mouse_ctl_summary_byRate %>% 
  group_by(rate_ssGSEA) %>% 
  dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(mouse_ctl_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ rate_ssGSEA") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

#return Seurat Object
saveRDS(mouse_ctl, 
        file = "./figures_mouse_ctrl_p16/mouse_ctl_SeuratObject.rds")

#Score analysis for p16 data
mouse_p16 <- obj3[,rownames(subset(obj3@meta.data, condition == "p16"))] #select only p16

mouse_p16 <- MetaModule_Assignment(mouse_p16, mouse_MetaModule_gene_sets, assay = "SCT")

#plot
p <- ggplot(mouse_p16@meta.data, aes(X, Y)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

p_cluster <- p +
  geom_point(aes(colour = cluster_names)) +
  scale_color_manual(values = cluster_colors)

p_cluster <- ggMarginal(p_cluster, type = "violin")

p_G2M <- p +
  geom_point(aes(color = G2.M), alpha = 0.5) + 
  scale_color_gradient2(low = "darkgrey", mid = "lightgrey", high = "red") 

gridExtra::grid.arrange(p_cluster, p_G2M, ncol=2)

p_priOPC <- ggplot(mouse_p16@meta.data, aes(X, Y, 
                                            colour = factor(cluster_names))) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
  geom_point() + 
  scale_color_manual(values = cluster_colors) + 
  gghighlight("pri-OPC-like_1" == cluster_names | 
                "pri-OPC-like_2" == cluster_names |
                "pri-OPC-like_3" == cluster_names | 
                "pri-OPC-like_4" == cluster_names, 
              keep_scales = TRUE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

p_cycling <- ggplot(mouse_p16@meta.data, aes(X, Y, colour = factor(cluster_names))) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
  geom_point() + 
  scale_color_manual(values = cluster_colors) + 
  gghighlight("G1/S_1" == cluster_names | 
                "G1/S_2" == cluster_names |
                "G2/M_1" == cluster_names | 
                "G2/M_2" == cluster_names, 
              keep_scales = TRUE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") +
  annotate("text", x = 1.5, y = -1.5, label = "MES", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = -1.5, label = "AC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = 1.5, y = 1.5, label = "NPC", 
           vjust = "inward", hjust = "inward") +
  annotate("text", x = -1.5, y = 1.5, label = "OPC",
           vjust = "inward", hjust = "inward")

gridExtra::grid.arrange(p_priOPC, p_cycling, ncol=2)

#Merge MES and NPC 
mouse_p16@meta.data$state <- gsub("[1-2]$",
                                  "",
                                  mouse_p16@meta.data$state)
## compute signature score
list_geneset = list(senecsence_sig = na.omit(biomart_murin_humain$Gene.name))
mouse_p16_ssGsea = gsva(as.matrix(mouse_p16@assays$SCT@data), list_geneset, 
                        method="ssgsea", ssgsea.norm=TRUE, verbose = FALSE)

mouse_p16@meta.data$score_ssGSEA <- mouse_p16_ssGsea[1,] * 100
mouse_p16@meta.data$rate_ssGSEA <- dplyr::ntile(mouse_p16@meta.data$score_ssGSEA, 10)
mouse_p16@meta.data$rate_ssGSEA <- ifelse(mouse_p16@meta.data$rate_ssGSEA == 1,
                                          "Low",
                                          ifelse(mouse_p16@meta.data$rate_ssGSEA == 10,
                                                 "High",
                                                 "Medium"))

#plots
plot1 <- FeaturePlot(mouse_p16, features = c("score_ssGSEA")) + 
  scale_colour_gradientn(colours = c("gray90", brewer.pal(n = 9, name = "Reds")))
plot2 <- UMAPPlot(mouse_p16, group.by = c("rate_ssGSEA"))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

UMAPPlot(mouse_p16, group.by = c("cluster_names"), 
         split.by = "rate_ssGSEA", 
         cols = cluster_colors)

## export results
write.table(mouse_p16_ssGsea, "./figures_mouse_ctrl_p16/mouse_p16_ssGsea.tsv", 
            quote = FALSE, 
            sep = "\t")

ggplot(mouse_p16@meta.data, 
       aes(x = cluster_names,
           y = score_ssGSEA,
           fill = cluster_names,
           colour = cluster_names)) +
  geom_violin(alpha=0.7, scale="width",adjust = .5) + 
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  geom_point(size=0.3,alpha = 0.8, position = "jitter")+
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  geom_hline(yintercept = max(mouse_p16@meta.data$score_ssGSEA[mouse_p16@meta.data$rate_ssGSEA == "Low"]),
             linetype="dashed") +
  geom_hline(yintercept = min(mouse_p16@meta.data$score_ssGSEA[mouse_p16@meta.data$rate_ssGSEA == "High"]),
             linetype="dashed") +
  ggtitle("Senescence score (ssGSEA) distribution between non-integrated mouse dataset", 
          "p16 transcriptomes")


mouse_p16_summary <- subset(mouse_p16@meta.data, 
                            select = c("rate_ssGSEA", "state", "cluster_names"))
# by Module
mouse_p16_summary_byRate <- dcast(mouse_p16_summary, 
                                  rate_ssGSEA ~ state)
mouse_p16_summary_byRate <- melt(mouse_p16_summary_byRate, 
                                 id.vars = "rate_ssGSEA")
mouse_p16_summary_byRate <- mouse_p16_summary_byRate %>% 
  group_by(rate_ssGSEA) %>% 
  dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(mouse_p16_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ rate_ssGSEA") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

#return Seurat Object
saveRDS(mouse_p16, 
        file = "./figures_mouse_ctrl_p16/mouse_p16_SeuratObject.rds")

