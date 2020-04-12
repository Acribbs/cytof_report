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
md_ms_p1 <- read.delim("newdata/experiment_35142_annotations.tsv", 
                       stringsAsFactors = F)
md_ms_p3 <- read.delim("newdata/experiment_35141_annotations.tsv", 
                       stringsAsFactors = F)
md_ms_p13 <- rbind(md_ms_p1, md_ms_p3)

md_ms_p13 <- md_ms_p13[,c(1,2,7)]
md_ms_p13$Conditions[6:10] <- md_ms_p13$Conditions[1:5]
colnames(md_ms_p13) <- c("file_name", "condition", "sample_id")
md_ms_p13$patient_id <- substr(md_ms_p13$sample_id,1,7)
md_ms_p13$panel <- rep(c("P1", "P3"), each = 5)
md_ms_p13

# import fcs
setwd("newdata/")
fcs_raw_ms_p1 <- read.flowSet(md_ms_p13$file_name[md_ms_p13$panel == "P1"],
                              transformation = FALSE,
                              truncate_max_range = FALSE)
fcs_raw_ms_p1
setwd("../")

setwd("newdata/")
fcs_raw_ms_p3 <- read.flowSet(md_ms_p13$file_name[md_ms_p13$panel == "P3"],
                              transformation = FALSE,
                              truncate_max_range = FALSE)
fcs_raw_ms_p3
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
all(panel1$fcs_colname %in% colnames(fcs_raw_ms_p1))
all(panel3$fcs_colname %in% colnames(fcs_raw_ms_p3))

# create SingleCellExperiment
sce_ms_p1 <- prepData(
  fcs_raw_ms_p1, 
  panel = panel1, 
  md_ms_p13[md_ms_p13$panel == "P1",],
  md_cols = list(file = "file_name", id = "sample_id", 
                 factors =  c("patient_id", "condition")))

sce_ms_p3 <- prepData(
  fcs_raw_ms_p3, 
  panel = panel3, 
  md_ms_p13[md_ms_p13$panel == "P3",],
  md_cols = list(file = "file_name", id = "sample_id", 
                 factors =  c("patient_id", "condition")))

n_cells(sce_ms_p1)
n_cells(sce_ms_p3)

#### CLUSTERING ####
set.seed(42)
sce_ms_p1 <- cluster(sce_ms_p1, features = type_markers(sce_ms_p1), 
                     xdim = 10, ydim = 10, maxK = 20, 
                     verbose = TRUE, seed = 42)    

set.seed(42)
sce_ms_p3 <- cluster(sce_ms_p3, features = type_markers(sce_ms_p3), 
                     xdim = 10, ydim = 10, maxK = 20, 
                     verbose = TRUE, seed = 42)    

# metadata(sce_ms_p1)$delta_area
# metadata(sce_ms_p3)$delta_area

plotClusterHeatmap(sce_ms_p1, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
plotClusterHeatmap(sce_ms_p1, hm2 = c("CD45", "CD66b"), k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
# remove populations 1,3,6,12 because CD45-
# remove populations 11,18 because CD66b+
sce_ms_p1_qc <- filterSCE(sce_ms_p1, k = "meta20",
                          !(cluster_id %in%c(1,3,6,12,11,18)))

plotClusterHeatmap(sce_ms_p3, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
plotClusterHeatmap(sce_ms_p3, hm2 = c("CyIg_K", "CyIg_L"), k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)
# remove populations 1,7,8,12,13,16,17,18,19,20 because CyIgK-CyIgL-
sce_ms_p3_qc <- filterSCE(sce_ms_p3, k = "meta20",
                          !(cluster_id %in%c(1,7,8,12,13,16,17,18,19,20)))

n_cells(sce_ms_p1_qc)
n_cells(sce_ms_p3_qc)

# Rerun the clustering
set.seed(42)
sce_ms_p1_qc <- cluster(sce_ms_p1_qc, features = type_markers(sce_ms_p1_qc), 
                        xdim = 10, ydim = 10, maxK = 20, 
                        verbose = TRUE, seed = 42)    

set.seed(42)
sce_ms_p3_qc <- cluster(sce_ms_p3_qc, features = type_markers(sce_ms_p3_qc), 
                        xdim = 10, ydim = 10, maxK = 20, 
                        verbose = TRUE, seed = 42)    

plotClusterHeatmap(sce_ms_p1_qc, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ms_p3_qc, hm2 = "abundances", k = "meta20",
                   draw_freqs = T, cluster_anno = F, scale = T)

#### VISUALIZATION
set.seed(42)
sce_ms_p1_qc <- runDR(sce_ms_p1_qc, "UMAP", 
                      features = type_markers(sce_ms_p1_qc),
                      cells = min(50000, min(n_cells(sce_ms_p1_qc))))

# set.seed(42)
# sce_ms_p1_qc <- runDR(sce_ms_p1_qc, "TSNE",
#                       features = type_markers(sce_ms_p1),
#                       cells = min(50000, min(n_cells(sce_ms_p1_qc))), 
#                       perplexity = 250)

set.seed(42)
sce_ms_p3_qc <- runDR(sce_ms_p3_qc, "UMAP", 
                      features = type_markers(sce_ms_p3_qc),
                      cells = min(50000, min(n_cells(sce_ms_p3_qc))))

# set.seed(42)
# sce_ms_p3_qc <- runDR(sce_ms_p3_qc, "TSNE",
#                       features = type_markers(sce_ms_p3),
#                       cells = min(50000, min(n_cells(sce_ms_p3_qc))), 
#                       perplexity = 250)

plotDR(x = sce_ms_p1_qc, dr = "UMAP", color_by = "meta20")
plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "meta12") + theme_bw()

plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "meta12") + 
  facet_wrap("sample_id", ncol = 5) + theme_bw()

plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "CD319_SLAM7") + 
  scale_color_viridis_c(limits = c(0,4), na.value = "#FDE725FF") + 
  theme_bw()

plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "Aiolos") + 
  scale_color_viridis_c(limits = c(0,5), na.value = "#FDE725FF") + 
  theme_bw()

plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "BCL_xL") + 
  scale_color_viridis_c(limits = c(0,5), na.value = "#FDE725FF") + 
  theme_bw()

plotDR(x = sce_ms_p3_qc, dr = "UMAP", color_by = "CD56") + 
  # scale_color_viridis_c(limits = c(0,5), na.value = "#FDE725FF") + 
  theme_bw()

p_ms <- plotMedExprs(sce_ms_p3_qc, facet = "antigen", 
                     group_by = "sample_id", shape_by = "patient_id")
p_ms$facet$params$ncol <- 4
p_ms

## Manual cluster annotation based on heatmaps
merging_table_ms <- "MS_P1_cluster_annot_200221.csv"
merging_table_ms <- read.csv(merging_table_ms)
dim(merging_table_ms)
head(data.frame(merging_table_ms))

# convert to factor with merged clusters in desired order
merging_table_ms$new_cluster_1 <- 
  factor(merging_table_ms$new_cluster_1)

sce_ms_p1_qc_annot <- mergeClusters(sce_ms_p1_qc, k = "meta20", 
                                    table = merging_table_ms[,c(1,2)],
                                    id = "merge01") 

plotDR(x = sce_ms_p1_qc_annot, dr = "UMAP", color_by = "merge01") +
  scale_colour_tableau("Classic 10") + theme_bw()

# plotDR(x = sce_ms_p1_qc_annot, dr = "UMAP", color_by = "merge01") +
#   scale_colour_tableau("Classic 10") +
#   facet_wrap("sample_id", ncol = 5)

# canonical marker heatmap
markers_hm <- c("CD11b", "CD14", "CD16", "CD123",
                "CD19", "CD3", "CD4", "CD8", "CD56")
sce_ms_p1_qc_annot_hm <- 
  sce_ms_p1_qc_annot[rownames(sce_ms_p1_qc_annot) %in% markers_hm,]

plotClusterHeatmap(sce_ms_p1_qc_annot_hm, hm2 = NULL, 
                   k = "merge01", draw_freqs = TRUE,
                   cluster_anno = FALSE, scale = T,
                   draw_dend = FALSE,
                   palette = brewer.pal(9, "RdPu"))

# Get abundances
fq_ms <- prop.table(table(cluster_ids(sce_ms_p1_qc_annot, "merge01"), 
                          sample_ids(sce_ms_p1_qc_annot)), 2) * 100
df_ms <- melt(fq_ms, value.name = "freq", 
              varnames = c("cluster_id", "sample_id"))
m <- match(df_ms$sample_id, ei(sce_ms_p1_qc_annot)$sample_id)
cols <- setdiff(names(ei(sce_ms_p1_qc_annot)), names(df_ms))
df_ms <- data.frame(df_ms, ei(sce_ms_p1_qc_annot)[m, cols])
df_ms$site <- "MSSM"

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
sce_ms_p1_t <- filterSCE(sce_ms_p1_qc_annot, k = "merge01",
                         cluster_id %in%c("T CELL CD4", "T CELL CD8"))

t_markers <- c("CCR7", "CD127", "CD137", "CD25", "CD28", "CD4",
               "CD45RA", "CD45RO", "CD56", "CD57", "CD8",
               "CD81", "CXCR3", "HLA_DR", "ICOS",
               "PD_1", "TIGIT", "TIM3_CD366")

# Rerun the clustering
set.seed(42)
sce_ms_p1_t <- cluster(sce_ms_p1_t, features = t_markers, 
                        xdim = 12, ydim = 12, maxK = 30, 
                        verbose = TRUE, seed = 42)

plotClusterHeatmap(sce_ms_p1_t, hm2 = "abundances", k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ms_p1_t, hm2 = c("CCR7", "CD45RA", "CD45RO"), 
                   k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

plotClusterHeatmap(sce_ms_p1_t, hm2 = c("CD127", "CD25"), 
                   k = "meta30",
                   draw_freqs = T, cluster_anno = F, scale = T)

set.seed(42)
sce_ms_p1_t <- runDR(sce_ms_p1_t, "UMAP", 
                     features = t_markers,
                     cells = min(25000, min(n_cells(sce_ms_p1_t))))

## Manual cluster annotation based on heatmaps
merging_table_ms_t <- "MS_P1_TCELL_cluster_annot_200221.csv"
merging_table_ms_t <- read.csv(merging_table_ms_t)
dim(merging_table_ms_t)
head(data.frame(merging_table_ms_t))

# convert to factor with merged clusters in desired order
merging_table_ms_t$new_cluster_2 <- 
  factor(merging_table_ms_t$new_cluster_2)

sce_ms_p1_t_annot <- mergeClusters(sce_ms_p1_t, k = "meta30", 
                                   table = merging_table_ms_t[,c(1,3)],
                                   id = "mergeTSUB") 

plotDR(x = sce_ms_p1_t_annot, dr = "UMAP", color_by = "mergeTSUB") +
  scale_colour_tableau("Classic 10") + theme_bw()

# Get abundances
fq_ms_t <- prop.table(table(cluster_ids(sce_ms_p1_t_annot, "mergeTSUB"), 
                            sample_ids(sce_ms_p1_t_annot)), 2) * 100
df_ms_t <- melt(fq_ms_t, value.name = "freq", 
                varnames = c("cluster_id", "sample_id"))
m <- match(df_ms_t$sample_id, ei(sce_ms_p1_t_annot)$sample_id)
cols <- setdiff(names(ei(sce_ms_p1_t_annot)), names(df_ms_t))
df_ms_t <- data.frame(df_ms_t, ei(sce_ms_p1_t_annot)[m, cols])
df_ms_t$site <- "MSSM"

df_all_t <- rbind(df_ms_t, df_ox_t)

ggpaired(df_all_t, x = "site", y = "freq", id = "patient_id",
         color = "site", line.color = "gray", line.size = 0.4,
         palette = c("#D80B8C", "blue")) +
  facet_wrap("cluster_id", ncol = 5) + theme_bw()

df_all_t %>% ggplot(aes_string(y = "freq")) +
  geom_bar(aes_string(x = "sample_id", fill = "factor(cluster_id)"), 
           position = "fill", stat = "identity") + 
  facet_wrap("site") +
  scale_fill_tableau("Classic 10") + theme_bw() +
  labs(x = "frequency",
       y = "",fill = "")


