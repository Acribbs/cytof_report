# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)

library(dplyr)
library(reshape2)
library(viridisLite)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)

#### DATA PREPARATION ####
# load in data
md_ox_p13 <- read.delim("newdata/experiment_35143_annotations.tsv", 
                       stringsAsFactors = F)

md_ox_p13 <- md_ox_p13[,c(1,2,7)]
colnames(md_ox_p13) <- c("file_name", "condition", "sample_id")
md_ox_p13$patient_id <- substr(md_ox_p13$sample_id,1,7)
md_ox_p13$panel <- rep(c("P3", "P1"), each = 5)
md_ox_p13 <- md_ox_p13[order(md_ox_p13$panel, md_ox_p13$patient_id),]
md_ox_p13

# import fcs
setwd("newdata/")
fcs_raw_ox_p1 <- read.flowSet(md_ox_p13$file_name[md_ox_p13$panel == "P1"],
                              transformation = FALSE,
                              truncate_max_range = FALSE)
fcs_raw_ox_p1
setwd("../")

setwd("newdata/")
fcs_raw_ox_p3 <- read.flowSet(md_ox_p13$file_name[md_ox_p13$panel == "P3"],
                              transformation = FALSE,
                              truncate_max_range = FALSE)
fcs_raw_ox_p3
setwd("../")

panel1_filename <- "panel1_alt.csv"
panel1 <- read.csv(panel1_filename)
panel1 <- panel1[order(panel1$antigen),]
head(data.frame(panel1))

panel3_filename <- "panel3_alt.csv"
panel3 <- read.csv(panel3_filename)
panel3 <- panel3[order(panel3$antigen),]
head(data.frame(panel3))

# spot check that all panel columns are in the flowSet object
all(panel1$fcs_colname %in% colnames(fcs_raw_ox_p1))
all(panel3$fcs_colname %in% colnames(fcs_raw_ox_p3))

# create SingleCellExperiment
sce_ox_p1 <- prepData(
  fcs_raw_ox_p1, 
  panel = panel1, 
  md_ox_p13[md_ox_p13$panel == "P1",],
  md_cols = list(file = "file_name", id = "sample_id", 
                 factors =  c("patient_id", "condition")))

sce_ox_p3 <- prepData(
  fcs_raw_ox_p3, 
  panel = panel3, 
  md_ox_p13[md_ox_p13$panel == "P3",],
  md_cols = list(file = "file_name", id = "sample_id", 
                 factors =  c("patient_id", "condition")))

n_cells(sce_ox_p1)
n_cells(sce_ox_p3)

#### CLUSTERING ####
set.seed(42)
sce_ox_p1 <- cluster(sce_ox_p1, features = type_markers(sce_ox_p1), 
                     xdim = 10, ydim = 10, maxK = 20, 
                     verbose = TRUE, seed = 42)    

set.seed(42)
sce_ox_p3 <- cluster(sce_ox_p3, features = type_markers(sce_ox_p3), 
                     xdim = 10, ydim = 10, maxK = 20, 
                     verbose = TRUE, seed = 42)    

# metadata(sce_ox_p1)$delta_area
# metadata(sce_ox_p3)$delta_area

plotClusterHeatmap(sce_ox_p1, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
plotClusterHeatmap(sce_ox_p1, hm2 = c("CD45", "CD66b"), k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
# remove populations 1,2,3,4,7,17 because CD45-
# remove populations 2,7,15,18 because CD66b+
sce_ox_p1_qc <- filterSCE(sce_ox_p1, k = "meta20",
                          !(cluster_id %in%c(1,2,3,4,7,15,17,18)))

plotClusterHeatmap(sce_ox_p3, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
plotClusterHeatmap(sce_ox_p3, hm2 = c("CyIg_K", "CyIg_L"), k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
# remove populations 1,2,3,6,7,8,9,14 because CyIgK-CyIgL-
sce_ox_p3_qc <- filterSCE(sce_ox_p3, k = "meta20",
                          !(cluster_id %in%c(1,2,3,6,7,8,9,14)))

n_cells(sce_ox_p1_qc)
n_cells(sce_ox_p3_qc)

# Rerun the clustering
set.seed(42)
sce_ox_p1_qc <- cluster(sce_ox_p1_qc, features = type_markers(sce_ox_p1_qc), 
                        xdim = 10, ydim = 10, maxK = 20, 
                        verbose = TRUE, seed = 42)    

set.seed(42)
sce_ox_p3_qc <- cluster(sce_ox_p3_qc, features = type_markers(sce_ox_p3_qc), 
                        xdim = 10, ydim = 10, maxK = 20, 
                        verbose = TRUE, seed = 42)    

plotClusterHeatmap(sce_ox_p1_qc, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ox_p3_qc, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)

#### VISUALIZATION
set.seed(42)
sce_ox_p1_qc <- runDR(sce_ox_p1_qc, "UMAP", 
                      features = type_markers(sce_ox_p1_qc),
                      cells = min(50000, min(n_cells(sce_ox_p1_qc))))

# set.seed(42)
# sce_ox_p1_qc <- runDR(sce_ox_p1_qc, "TSNE",
#                       features = type_markers(sce_ox_p1),
#                       cells = min(50000, min(n_cells(sce_ox_p1_qc))), 
#                       perplexity = 250)

set.seed(42)
sce_ox_p3_qc <- runDR(sce_ox_p3_qc, "UMAP", 
                      features = type_markers(sce_ox_p3_qc),
                      cells = min(50000, min(n_cells(sce_ox_p3_qc))))

# set.seed(42)
# sce_ox_p3_qc <- runDR(sce_ox_p3_qc, "TSNE",
#                       features = type_markers(sce_ox_p3),
#                       cells = min(50000, min(n_cells(sce_ox_p3_qc))), 
#                       perplexity = 250)

plotDR(x = sce_ox_p1_qc, dr = "UMAP", color_by = "meta20")
plotDR(x = sce_ox_p3_qc, dr = "UMAP", color_by = "meta12") + theme_bw()

plotDR(x = sce_ox_p3_qc, dr = "UMAP", color_by = "meta12") + 
  facet_wrap("sample_id", ncol = 5) + theme_bw()

sce_ox_p3_qc_similar <- filterSCE(sce_ox_p3_qc, k = "meta12",
                                  !(cluster_id %in% c(2)))

p_ox <- plotMedExprs(sce_ox_p3_qc_similar, facet = "antigen", 
                     group_by = "sample_id", shape_by = "patient_id")
p_ox$facet$params$ncol <- 4
p_ox

plotClusterHeatmap(sce_ox_p3_qc, hm2 = "abundances", k = "meta12",
                   draw_freqs = T, cluster_anno = T, scale = F)

## Manual cluster annotation based on heatmaps
merging_table_ox <- "OX_P1_cluster_annot_200221.csv"
merging_table_ox <- read.csv(merging_table_ox)
dim(merging_table_ox)
head(data.frame(merging_table_ox))

# convert to factor with merged clusters in desired order
merging_table_ox$new_cluster_1 <- 
  factor(merging_table_ox$new_cluster_1)

sce_ox_p1_qc_annot <- mergeClusters(sce_ox_p1_qc, k = "meta20", 
                                    table = merging_table_ox[,c(1,2)],
                                    id = "merge01") 

plotDR(x = sce_ox_p1_qc_annot, dr = "UMAP", color_by = "merge01") +
  scale_colour_tableau("Classic 10") + theme_bw()

# plotDR(x = sce_ox_p1_qc_annot, dr = "UMAP", color_by = "merge01") +
#   scale_colour_tableau("Classic 10") +
#   facet_wrap("sample_id", ncol = 5)

# canonical marker heatmap
markers_hm <- c("CD11b", "CD14", "CD16", "CD123",
                "CD19", "CD3", "CD4", "CD8", "CD56")
sce_ox_p1_qc_annot_hm <- 
  sce_ox_p1_qc_annot[rownames(sce_ox_p1_qc_annot) %in% markers_hm,]

plotClusterHeatmap(sce_ox_p1_qc_annot_hm, hm2 = NULL, 
                   k = "merge01", draw_freqs = TRUE,
                   cluster_anno = FALSE, scale = T,
                   draw_dend = FALSE,
                   palette = brewer.pal(9, "Blues"))

# Get abundances
fq_ox <- prop.table(table(cluster_ids(sce_ox_p1_qc_annot, "merge01"), 
                          sample_ids(sce_ox_p1_qc_annot)), 2) * 100
df_ox <- melt(fq_ox, value.name = "freq", 
              varnames = c("cluster_id", "sample_id"))
m <- match(df_ox$sample_id, ei(sce_ox_p1_qc_annot)$sample_id)
cols <- setdiff(names(ei(sce_ox_p1_qc_annot)), names(df_ox))
df_ox <- data.frame(df_ox, ei(sce_ox_p1_qc_annot)[m, cols])
df_ox$site <- "OXFORD"

df_all <- rbind(df_ms, df_ox)

ggpaired(df_all, x = "site", y = "freq", id = "patient_id",
         color = "site", line.color = "gray", line.size = 0.4,
         palette = c("#D80B8C", "blue")) +
  stat_compare_means(paired = TRUE) + 
  facet_wrap("cluster_id", ncol = 5) + theme_bw()

df_all %>% ggplot(aes_string(y = "freq")) +
  geom_bar(aes_string(x = "sample_id", fill = "factor(cluster_id)"), 
           position = "fill", stat = "identity") + 
  facet_wrap("site") +
  scale_fill_tableau("Classic 10") + theme_bw() +
  labs(x = "frequency",
       y = "",fill = "")

## FILTER T CELLS
sce_ox_p1_t <- filterSCE(sce_ox_p1_qc_annot, k = "merge01",
                         cluster_id %in%c("T CELL CD4", "T CELL CD8"))

t_markers <- c("CCR7", "CD127", "CD137", "CD25", "CD28", "CD4",
               "CD45RA", "CD45RO", "CD56", "CD57", "CD8",
               "CD81", "CXCR3", "HLA_DR", "ICOS",
               "PD_1", "TIGIT", "TIM3_CD366")

# Rerun the clustering
set.seed(42)
sce_ox_p1_t <- cluster(sce_ox_p1_t, features = t_markers, 
                       xdim = 12, ydim = 12, maxK = 30, 
                       verbose = TRUE, seed = 42)

plotClusterHeatmap(sce_ox_p1_t, hm2 = "abundances", k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ox_p1_t, hm2 = c("CCR7", "CD45RA", "CD45RO"), 
                   k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ox_p1_t, hm2 = c("CD127", "CD25"), 
                   k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

set.seed(42)
sce_ox_p1_t <- runDR(sce_ox_p1_t, "UMAP", 
                     features = t_markers,
                     cells = min(25000, min(n_cells(sce_ox_p1_t))))

## Manual cluster annotation based on heatmaps
merging_table_ox_t <- "OX_P1_TCELL_cluster_annot_200221.csv"
merging_table_ox_t <- read.csv(merging_table_ox_t)
dim(merging_table_ox_t)
head(data.frame(merging_table_ox_t))

# convert to factor with merged clusters in desired order
merging_table_ox_t$new_cluster_2 <- 
  factor(merging_table_ox_t$new_cluster_2)

sce_ox_p1_t_annot <- mergeClusters(sce_ox_p1_t, k = "meta30", 
                                   table = merging_table_ox_t[,c(1,3)],
                                   id = "mergeTSUB") 

plotDR(x = sce_ox_p1_t_annot, dr = "UMAP", color_by = "mergeTSUB") +
  scale_colour_tableau("Classic 10") + theme_bw()

# Get abundances
fq_ox_t <- prop.table(table(cluster_ids(sce_ox_p1_t_annot, "mergeTSUB"), 
                            sample_ids(sce_ox_p1_t_annot)), 2) * 100
df_ox_t <- melt(fq_ox_t, value.name = "freq", 
                varnames = c("cluster_id", "sample_id"))
m <- match(df_ox_t$sample_id, ei(sce_ox_p1_t_annot)$sample_id)
cols <- setdiff(names(ei(sce_ox_p1_t_annot)), names(df_ox_t))
df_ox_t <- data.frame(df_ox_t, ei(sce_ox_p1_t_annot)[m, cols])
df_ox_t$site <- "OXFORD"

df_all_t <- rbind(df_ms_t, df_ox_t)

ggpaired(df_all_t, x = "site", y = "freq", id = "patient_id",
         color = "site", line.color = "gray", line.size = 0.4,
         palette = c("#D80B8C", "blue")) +
  facet_wrap("cluster_id", ncol = 5) + theme_bw()

df_all_t %>% ggplot(aes_string(y = "freq")) +
  geom_bar(aes_string(x = "sample_id", fill = "factor(cluster_id)"), 
           position = "fill", stat = "identity") + 
  facet_wrap("site", ncol = 1) +
  scale_fill_tableau("Classic 10") + theme_bw() +
  labs(x = "frequency",
       y = "",fill = "")
