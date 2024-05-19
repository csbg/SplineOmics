library(here)
# library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(scales)
library(dplyr)


# load hits data ---------------------------------------------------------------

load("data/clust_hits_final_27_03_2024-15_49_10_PTX_exp_.RData")
clusters_exp <- clustered_hits %>%
  arrange(cluster)
load("data/clust_hits_final_27_03_2024-17_21_02_PTX_stat_.RData")
clusters_sta <- clustered_hits %>%
  arrange(cluster)

# load data matrix and meta -----------------------------------------------

load(here("data", "data_meta_tp.RData"))
load(here("data", "data_batch_corrected.RData"))

# limma results -----------------------------------------------

load(here("data", "results_dof2_splines_logfc.RData"))

# subset dataset and meta -------------------------------------------------
data.matrix.batch.exp <- data.matrix.batch[,1:18]
data.matrix.batch.exp.filt <- data.matrix.batch.exp[, -which(colnames(data.matrix.batch.exp) == "E12_TP05_Exponential")]
data.matrix.batch.sta <- data.matrix.batch[,19:36]
data.matrix.batch.sta.filt <- data.matrix.batch.sta[, !(colnames(data.matrix.batch.sta) %in% c("E10_TP10_Stationary"))]
data.matrix.batch.filt <- data.matrix.batch[, !(colnames(data.matrix.batch) %in% c("E12_TP05_Exponential", "E10_TP10_Stationary"))]

meta_exp <- meta %>% filter(phase_of_fermentation == "Exponential")
meta_exp_filt <- meta_exp %>% filter(!(sample_name %in% c("E12_TP05_Exponential")))
meta_sta <- meta %>% filter(phase_of_fermentation == "Stationary")
meta_sta_filt <- meta_sta %>% filter(!(sample_name %in% c("E10_TP10_Stationary")))
meta_filt <- meta %>% filter(!(sample_name %in% c("E12_TP05_Exponential", "E10_TP10_Stationary")))

# subset data needed for plotting -----------------------------------------
#Exponential phase
data.matrix.batch.filt.sig <- data.matrix.batch.exp.filt[as.numeric(clusters_exp$feature),]
z_score <- t(scale(t(data.matrix.batch.filt.sig)))
rownames(z_score) <- NULL

# plot heatmaps -----------------------------------------------------------

BASE_TEXT_SIZE_PT <- 5

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)

width = 20
height = 60

# png(filename = paste0("./figures/heatmap_filt_exp_noNames.png"),
#     width = width,
#     height = height,
#     units = "cm",
#     res = 600)

# knitr::plot_crop("./figures/heatmap_filt_exp_noNames.png")
# 
# pdf(file = paste0("./figures/heatmap_filt_exp.pdf"),
#     width = width / 2.54,
#     height = height / 2.54)

(ht <- Heatmap(z_score,
               column_split = meta_exp_filt$time_point,
               cluster_columns = FALSE,
               row_split = clusters_exp$cluster,
               cluster_rows = FALSE,
               heatmap_legend_param = list(title = "z-score of log2 intensity",title_position = "lefttop-rot"),
               row_gap = unit(2, "pt"),
               column_gap = unit(2, "pt"),
               width = unit(2, "mm") * ncol(z_score) + 5 * unit(2, "pt"), # to make each cell a square
               height = unit(2, "mm") * nrow(z_score) + 5 * unit(2, "pt"), # to make each cell a square
               show_row_names = TRUE
))


dev.off()

# plot normalized log2 intensities ----------------------------------------
# prepare data for plotting
data_to_plot <- data.frame(first_protein_description = rownames(data.matrix.batch.filt.sig), 
                           data.matrix.batch.filt.sig) %>%  # add row with the first protein descr
  mutate(cluster_number = clusters_exp$cluster) %>%
  pivot_longer(cols = colnames(data.matrix.batch.filt.sig),
               names_to = "sample_name",
               values_to = "log2_intensity") %>%
  separate(sample_name, into = c("reactor", "time_point", "phase_of_fermentation"),sep = "_") %>%
  mutate(time_to_feed = rep(meta_exp_filt$time_to_feed,length(clusters_exp$feature))) %>%
  group_by(first_protein_description) %>%
  mutate(log2_intensity = rescale(log2_intensity))# Normalize log2_intensity values to [0, 1]


data_to_plot_mean_protein  <- data_to_plot %>%
  group_by(first_protein_description, time_to_feed) %>%
  mutate(mean_intensity_protein = mean(log2_intensity)) %>%
  ungroup()

data_to_plot_mean_tp <- data_to_plot %>%
  group_by(cluster_number, time_to_feed) %>%
  summarise(mean_intensity_tp = mean(log2_intensity)) 

#count number of genes in each cluster
clusters_exp %>%
  count(cluster)


ggplot(data = data_to_plot_mean_protein) +
  geom_line(aes(x = time_to_feed, y = mean_intensity_protein, color = first_protein_description), alpha = 0.5) +
  geom_line(data = data_to_plot_mean_tp,aes(x = time_to_feed, y = mean_intensity_tp), linewidth = 0.8) +
  facet_wrap(~cluster_number, ncol = 2, labeller = labeller(cluster_number =  c("1" = "Cluster 1: 32 proteins",
                                                                                "2" = "Cluster 2: 48 proteins",
                                                                                "3" = "Cluster 3: 68 proteins",
                                                                                "4" = "Cluster 4: 29 proteins",
                                                                                "5" = "Cluster 5: 17 proteins",
                                                                                "6" = "Cluster 6: 13 proteins")))  +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red', linewidth = 0.5) +  
  ylab("normalized log2 intensity") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = data_to_plot$time_to_feed)

plot_file_path <- paste0("./figures/dof_2/curve_kmeans_clustered/filt_exp_composite_plot.png")

ggsave(plot_file_path,
       bg = 'white')





















#Stationaryphase
data.matrix.batch.filt.sig <- data.matrix.batch.sta.filt[as.numeric(clusters_sta$feature),]
z_score <- t(scale(t(data.matrix.batch.filt.sig)))
# rownames(z_score) <- NULL

# plot heatmaps -----------------------------------------------------------
png(filename = paste0("./figures/dof_2/curve_kmeans_clustered/heatmap_filt_sta.png"),
    width = 210,
    height = 180,
    units = "mm",
    res = 600)

pdf(file = paste0("./figures/dof_2/curve_kmeans_clustered/heatmap_filt_exp.pdf"),
    width = 8,
    height = 11,
    paper = "A4")

ht <- Heatmap(z_score,
              column_split = meta_sta_filt$time_point,
              cluster_columns = FALSE,
              row_split = clusters_sta$cluster,
              cluster_rows = FALSE,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              row_names_gp = grid::gpar(fontsize = 8),
              heatmap_legend_param = list(title = "z-score of log2 intensity",title_position = "lefttop-rot"),
              width = unit(10, "cm"),
              height = unit(10, "cm")
              # heatmap_width = unit(15, "cm"), 
              # heatmap_height = unit(27, "cm")
)

ht = draw(ht)
dev.off()
# plot normalized log2 intensities ----------------------------------------
# prepare data for plotting
data_to_plot <- data.frame(first_protein_description = rownames(data.matrix.batch.filt.sig), 
                           data.matrix.batch.filt.sig) %>%
  mutate(cluster_number = clusters_sta$cluster) %>%
  pivot_longer(cols = colnames(data.matrix.batch.filt.sig),
               names_to = "sample_name",
               values_to = "log2_intensity") %>%
  separate(sample_name, into = c("reactor", "time_point", "phase_of_fermentation"),sep = "_") %>%
  mutate(time_to_feed = rep(meta_sta_filt$time_to_feed,length(clusters_sta$feature))) %>%
  group_by(first_protein_description) %>%
  mutate(log2_intensity = rescale(log2_intensity))# Normalize log2_intensity values to [0, 1]


data_to_plot_mean_protein  <- data_to_plot %>%
  group_by(first_protein_description, time_to_feed) %>%
  mutate(mean_intensity_protein = mean(log2_intensity)) %>%
  ungroup()

data_to_plot_mean_tp <- data_to_plot %>%
  group_by(cluster_number, time_to_feed) %>%
  summarise(mean_intensity_tp = mean(log2_intensity)) 

#count number of genes in each cluster
clusters_sta %>%
  count(cluster)


ggplot(data = data_to_plot_mean_protein) +
  geom_line(aes(x = time_to_feed, y = mean_intensity_protein, color = first_protein_description), alpha = 0.5) +
  geom_line(data = data_to_plot_mean_tp,aes(x = time_to_feed, y = mean_intensity_tp), linewidth = 0.8) +
  facet_wrap(~cluster_number, ncol = 1, labeller = labeller(cluster_number =  c("1" = "Cluster 1: 11 proteins",
                                                                                "2" = "Cluster 2: 2 proteins",
                                                                                "3" = "Cluster 3: 1 protein"
                                                                                # "4" = "Cluster 4: 29 proteins",
                                                                                # "5" = "Cluster 5: 17 proteins",
                                                                                # "6" = "Cluster 6: 13 proteins"
  )))  +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red', linewidth = 0.5) +  
  ylab("normalized log2 intensity") +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = data_to_plot$time_to_feed)

plot_file_path <- paste0("./figures/dof_2/curve_kmeans_clustered/filt_sta_composite_plot.png")

ggsave(plot_file_path,
       bg = 'white')