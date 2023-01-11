## Zscores script
# Compute zscore of senescence scores of multiple datasets
# Salam, Saliou et al, 2022 Nature Comm


#Import Libraries
library(Seurat)
library(limma)
library(ggplot2)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(ggExtra)
library(gghighlight)
library(dplyr)
library(psych)
library(gridExtra)

# Custom functions

rate_decile <- function(score){
  #Create categories Low/Medium/High based on the deciles of a score
  rate = dplyr::ntile(score, 10)
  rate <- ifelse(rate == 1,
                 "Low",
                 ifelse(rate == 10,
                        "High",
                        "Medium"))
  rate <- factor(rate, levels = c("High", "Medium", "Low"), ordered = FALSE)
  return(rate)
}

plots <- function(df, tumour_colnames, dataset){
  # return several plots
  
  # chosen palette
  if(dataset == "mouse"){
    cols = cluster_colors
  } else{
    cols = col_vector
  }
  
  #violin plot
  print(ggplot(df,
               aes_string(x = tumour_colnames,
                          y = "zscore",
                          fill = tumour_colnames,
                          colour = tumour_colnames)) +
          geom_violin(alpha = 0.7,
                      scale = "width",
                      adjust = .5) +
          scale_fill_manual(values = cols) +
          geom_point(size = 0.3,
                     alpha = 0.8,
                     position = "jitter") +
          scale_colour_manual(values = cols) +
          ggpubr::theme_pubr() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Senescence z-score (ssGSEA) distribution between non-integrated patient datasets"))
  
  #violin plot
  print(ggplot(df,
               aes_string(
                 x = tumour_colnames,
                 y = "zscore",
                 fill = tumour_colnames,
                 colour = tumour_colnames
               )) +
          geom_violin(alpha = 0.7,
                      scale = "width",
                      adjust = .5) +
          scale_fill_manual(values = cols) +
          geom_point(size = 0.3,
                     alpha = 0.8,
                     position = "jitter") +
          scale_colour_manual(values = cols) +
          facet_grid( ~ orig.ident + status, scales = "free") +
          ggpubr::theme_pubr() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Senescence z-score (ssGSEA) distribution between non-integrated patient datasets"))
  
  #violin plot
  print(ggplot(df,
               aes_string(
                 x = tumour_colnames,
                 y = "score_ssGSEA",
                 fill = tumour_colnames,
                 colour = tumour_colnames
               )) +
          geom_violin(alpha = 0.7,
                      scale = "width",
                      adjust = .5) +
          scale_fill_manual(values = cols) +
          geom_point(size = 0.3,
                     alpha = 0.8,
                     position = "jitter") +
          scale_colour_manual(values = cols) +
          facet_grid( ~ orig.ident + status, scales = "free") +
          ggpubr::theme_pubr() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Senescence score ssGSEA distribution between non-integrated patient datasets"))  
  
  #violin plot
  print(ggplot(df,
               aes_string(
                 x = tumour_colnames,
                 y = "zscore",
                 fill = tumour_colnames,
                 colour = tumour_colnames
               )) +
          geom_violin(alpha = 0.7,
                      scale = "width",
                      adjust = .5) +
          scale_fill_manual(values = cols) +
          geom_point(size = 0.3,
                     alpha = 0.8,
                     position = "jitter") +
          scale_colour_manual(values = cols) +
          facet_grid( ~ orig.ident, scales = "free") +
          ggpubr::theme_pubr() + ylab("sen-z-score") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
          ggtitle("Senescence z-score (ssGSEA) distribution between non-integrated patient datasets"))
  
  print(ggplot(df, 
               aes_string(x = tumour_colnames, fill = "zscore_rate")) +
          geom_bar(position = "fill") +
          facet_grid( ~ orig.ident, scales = "free_x") +
          ggpubr::theme_pubr() + ylab("z-score rate distribution (%)") +
          scale_fill_manual(values = rate_colors) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
          ggtitle(
            "Distribution du z-score en fonction des tumeurs",
            "Deciles were computed based on all samples"
          ))
  

  hist_plot <- ggplot(df, 
                      aes_string(x = tumour_colnames, fill = "zscore_rate")) +
    geom_bar(position = "fill") +
    facet_grid( ~ orig.ident + status, scales = "free") +
    scale_fill_manual(values = rate_colors) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Distribution du z-score en fonction des tumeurs",
            "Deciles were computed based on all samples")
  print(hist_plot)
  
  #export frequency for each zscore rate
  write.table(as.data.frame(table(subset(df, 
                                         select = c("zscore_rate", 
                                                    tumour_colnames, 
                                                    ifelse(dataset == "human", 
                                                           "orig.ident", 
                                                           "status")))), 
                            stringsAsFactors = F), 
              file = paste0(dataset, "_zscores_rate_freq.tsv"),
              row.names = T, col.names = T,
              quote = F, sep = "\t")
  
  print(ggplot(df, 
               aes_string(x = "zscore_rate", fill = "state")) +
          geom_bar(position = "fill") +
          scale_fill_manual(values = module_colors) +
          ggpubr::theme_pubr() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Distribution des module en fonction du Rate des zscores",
                  "Deciles were computed based on all samples"))
  
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

rate_colors <- c("High" = "#0A2F51", "Medium" = "#DCDCDC", "Low" = "#BFE1B0")
module_colors <- c("AClike" = "#36989E", "OPClike" = "#9D8AF4", "MESlike" = "#F3C86F", "NPClike" = "#F3866F")

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_col_pals <- qual_col_pals[c("Set3", "Set2", "Pastel1", "Pastel2"),]
col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

#import Seurat object from RDS files
verhaak_object <- readRDS("./figures_johnson_nature_2021/verhaak_SeuratObject.rds")
mouse_ctl <- readRDS("./figures_mouse_ctrl_p16/mouse_ctl_SeuratObject.rds")
mouse_p16 <- readRDS("./figures_mouse_ctrl_p16/mouse_p16_SeuratObject.rds")
neftel_10X_notmutate <- readRDS("./figures_neftel_10X_SS2/neftel_10X_notmutate_SeuratObject.rds")
neftel_SS2_notmutate <- readRDS("./figures_neftel_10X_SS2/neftel_smartseq_notmutate_SeuratObject.rds")
neftel_SS2_mutate <- readRDS("./figures_neftel_10X_SS2_p16/neftel_smartseq_mutate_SeuratObject.rds")
badhuri <- readRDS("./figures_badhuri/badhuri_SeuratObject.rds")

#create a new column for all dataset where to find the tumour identity
verhaak_object$Sample <- verhaak_object$case_barcode
verhaak_object$orig.ident <- "Johnson"
verhaak_object$status <- "unknown"
DefaultAssay(verhaak_object) <- "RNA"

neftel_10X_notmutate$Sample <- neftel_10X_notmutate$tumour.name
neftel_10X_notmutate$orig.ident <- "Neftel (10X)"
neftel_10X_notmutate$status <- "not mutate"
DefaultAssay(neftel_10X_notmutate) <- "RNA"

neftel_SS2_notmutate$orig.ident <- "Neftel (SS2)"
neftel_SS2_notmutate$status <- "not mutate"
DefaultAssay(neftel_SS2_notmutate) <- "RNA"

neftel_SS2_mutate$orig.ident <- "Neftel (SS2)"
neftel_SS2_mutate$status <- "mutate"
DefaultAssay(neftel_SS2_mutate) <- "RNA"

badhuri$orig.ident <- "Badhuri"
badhuri$Sample <- badhuri$Tumor.ID
badhuri$status <- "unknown"
DefaultAssay(badhuri) <- "RNA"


#combine into one seurat object all human datasets
human_datasets <- merge(verhaak_object, 
                        list(neftel_10X_notmutate, neftel_SS2_notmutate, neftel_SS2_mutate, badhuri))

#combine into one seurat object all mouse datasets
mouse_dataset <- merge(mouse_ctl, mouse_p16)
mouse_dataset$status <- ifelse(grepl("ctl", mouse_dataset$orig.ident),
                               "WT+GCV",
                               "p16-3MR+GCV")
mouse_dataset$status <- factor(mouse_dataset$status, 
                               levels = c("WT+GCV", "p16-3MR+GCV"), 
                               ordered = TRUE)
mouse_dataset$Sample <- mouse_dataset$cluster_names
mouse_dataset$orig.ident <- "Mouse"


#combine into one seurat object all datasets
all_datasets <- merge(human_datasets, mouse_dataset)

#Compute zscore rate for human datasets
human_datasets$zscore <- scale(human_datasets$score_ssGSEA, center = TRUE, scale = TRUE)
human_datasets$zscore_rate <- rate_decile(human_datasets$zscore)
human_datasets$rate_ssGSEA <- factor(human_datasets$rate_ssGSEA, levels = c("High", "Medium", "Low"), ordered = FALSE)

plots(df = human_datasets@meta.data, tumour_colnames = "Sample", dataset = "human")

#pie charts
human_summary <- subset(human_datasets@meta.data, select = c("zscore_rate", "state", "Sample"))
# by Module
human_summary_byRate <- dcast(human_summary, zscore_rate ~ state)
human_summary_byRate <- melt(human_summary_byRate, id.vars = "zscore_rate")
human_summary_byRate <- human_summary_byRate %>% group_by(zscore_rate) %>% dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(human_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ zscore_rate") + 
  geom_col(position = 'fill') +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  scale_fill_manual(values = module_colors) +
  coord_polar(theta = 'y')

human_resume <- dcast(human_summary, state ~ zscore_rate)
human_resume$total <- apply(human_resume[,-1], 1, sum)
human_resume[nrow(human_resume)+1, ] <- c("Total", apply(human_resume[,-1], 2, sum))

knitr::kable(human_resume, caption = "Effectif des différents state en fonction des catégories de z-score (ssGSEA) tous datasets de patient confondus.")

#save data
write.table(human_datasets@meta.data, 
            file = "human_datasets_metadata.tsv", 
            row.names = T, col.names = T, 
            quote = F, sep = "\t")

human_datasets_list <- split(human_datasets@meta.data, human_datasets@meta.data$orig.ident)
lapply(human_datasets_list, function(x){
  print(ggplot(x, aes(X, Y, colour = Sample)) +
          theme_bw() + theme(legend.position = "bottom",
                             legend.title = element_text(size = 8),
                             legend.text = element_text(size = 6)) +
          labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
               y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
          geom_point() +
          gghighlight("High" == zscore_rate , 
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
                   vjust = "inward", hjust = "inward"))
})

#Compute zscore rate for mouse samples
mouse_dataset$zscore <- scale(mouse_dataset$score_ssGSEA, center = TRUE, scale = TRUE)
mouse_dataset$zscore_rate <- rate_decile(mouse_dataset$zscore)
mouse_dataset$rate_ssGSEA <- factor(mouse_dataset$rate_ssGSEA, levels = c("High", "Medium", "Low"), ordered = FALSE)

plots(df = mouse_dataset@meta.data, tumour_colnames = "Sample", dataset = "mouse")

#save data
write.table(mouse_dataset@meta.data, 
            file = "mouse_metadata.tsv", 
            row.names = T, col.names = T, 
            quote = F, sep = "\t")

#pie charts
mouse_summary <- subset(mouse_dataset@meta.data, select = c("zscore_rate", "state", "Sample"))
# by Module
mouse_summary_byRate <- dcast(mouse_summary, zscore_rate ~ state)
mouse_summary_byRate <- melt(mouse_summary_byRate, id.vars = "zscore_rate")
mouse_summary_byRate <- mouse_summary_byRate %>% group_by(zscore_rate) %>% dplyr::mutate(label = round(value / sum(value) * 100, 2))

ggplot(mouse_summary_byRate, aes('', value, fill = variable)) + 
  facet_wrap(". ~ zscore_rate") + 
  geom_col(position = 'fill') +
  scale_fill_manual(values = module_colors) +
  geom_label(aes(label = label), position = position_fill(vjust = 0.5), label.size = 0.1) +
  coord_polar(theta = 'y')

mouse_resume <- dcast(mouse_summary, state ~ zscore_rate)
mouse_resume$total <- apply(mouse_resume[,-1], 1, sum)
mouse_resume[nrow(mouse_resume)+1, ] <- c("Total", apply(mouse_resume[,-1], 2, sum))

ggplot(subset(mouse_dataset@meta.data, status == "WT+GCV"), aes(X, Y, colour = zscore_rate)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6)) +
  labs(x = "Relative meta-module score [log2(|SC1-SC2|+1)]",
       y = "Relative meta-module score [log2(|SC1-SC2|+1)]") + 
  geom_point() + 
  scale_color_manual(values = rate_colors) + 
  gghighlight("High" == zscore_rate , 
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
           vjust = "inward", hjust = "inward") + 
  ggtitle("WT+GCV", "High senescence score cells")