library(devtools)
library(ggplot2)
library(SPACO)
library(Seurat)
library(SPARK)
devtools::load_all("/home/rstudio/data_dir/R_scripts/SpaCo_R_kia")

distances_bootstrap <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/distances_bootstrap.RDS")
spacs_bootstrap <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/spacs_bootstrap.RDS")
brain <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/brain.RDS")
orig_data <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/twin_data.RDS")
orig_names <- rownames(orig_data)
orig_data<- orig_data[rowSums(orig_data) != 0, ]
SpaCoObject <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/newSCAMoransIbrainnsim1000.RDS")
neighbours1to50 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/neighbours1to50.RDS")
coverage_adjusted_bootstrap <- function(orig_data, neighborhood_matrix){
  
  resample_rows_fixed_neighborhood <- function(data_matrix, neighborhood_matrix) {
    num_cols <- ncol(data_matrix)
    num_rows <- nrow(data_matrix)
    new_data <- matrix(data = NA, nrow = num_rows,ncol = num_cols)
    colnames(new_data) <- colnames(data_matrix)
    rownames(new_data) <- rownames(data_matrix)
    for (j in 1:num_rows) {
      neighbors <- which(neighborhood_matrix[,j ] == 1)
      for (i in 1:num_cols) {
        if (length(neighbors) > 0) {
          random_neighbor <- sample(neighbors, 1, replace = F)
          
          new_data[j, i ] <- data_matrix[random_neighbor, i ]
        }
      }
    }
    colnames(new_data) <- colnames(data_matrix)
    rownames(new_data) <- rownames(data_matrix)
    return( new_data)
  }
  resample_coverage_fixed_neighborhood <- function(coverage_vector, neighborhood_matrix) {
    num_cov <- length(coverage_vector)
    new_data <- vector(length=length(coverage_vector))
    names(new_data) <- names(coverage_vector)
    for (i in 1:num_cov) {
      neighbors <- which(neighborhood_matrix[,i ] == 1)
      if (length(neighbors) > 0) {
        random_neighbor <- sample(neighbors, 1, replace = F)
        
        new_data[i] <- coverage_vector[random_neighbor]
      }
    }
    names(new_data) <- names(coverage_vector)
    return( new_data)
  }
  
  coverage_vector <- colSums(orig_data)
  coverage_vector_rand <- resample_coverage_fixed_neighborhood(coverage_vector,neighborhood_matrix)
  
  orig_data_abund <- sweep(orig_data, 2, coverage_vector, "/")
  sample_matrix_p <- t(resample_rows_fixed_neighborhood(t(orig_data_abund),neighborhood_matrix))
  counts_sample <- t(t(sample_matrix_p)*(coverage_vector_rand/colSums(sample_matrix_p)))
  counts_sample <- round(counts_sample)
  return(counts_sample)
}
c(1,2,3,4,5,6,7,8,9,10,15,20,25)
set.seed(1000 + 5)
#bootstrapped_data <- coverage_adjusted_bootstrap(orig_data = orig_data,neighbours1to50[[5]])
#bootstrapped_data<- bootstrapped_data[rowSums(bootstrapped_data) != 0, ]

#adt_assay <- CreateAssayObject(counts = bootstrapped_data)
#brain[["shuffle"]] <- adt_assay
#brain <- PercentageFeatureSet(brain, pattern = "^mt-" ,col.name = "percent.mt")
#brain <- PercentageFeatureSet(brain, pattern = "^Hbb-" ,col.name = "percent.hbb")
#brain <- SCTransform(brain, assay = "shuffle", variable.features.n = 3000, verbose = F ,seed.use = 1000 + 5 )
#SpaCoObject_twin <- seurat_to_spaco(Seurat = brain, assay = "SCT", n_image= 1, slot = "scale.data")
#SpaCoObject_twin <- RunSCAI(SpaCoObject_twin, PC_criterion = "percent",
#                            PC_value = .8, compute_nSpacs = T,
 #                          compute_projections = TRUE, nSpacQuantile = 0.05, nSim = 1000 )
SpaCoObject_twin <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/Spaco_bootstrap_radius5.RDS")
slot(SpaCoObject_twin,"meta.data") <- new("data.frame")
DE_rad5 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/DE_rad5.RDS")
#DE_rad5 <- SVGTest(SpaCoObject_twin)

print(SpaCoObject_twin@nSpacs)
SpaCoObject_twin <- denoise_profiles(SpaCoObject_twin)

feature_plot(SpaCoObject,"Ttr")+ggtitle("original")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_ttr_orig.pdf", width = 297, height = 210, units = "mm")

feature_plot(SpaCoObject_twin,"Ttr")+ggtitle("bootstrapped")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_ttr_bootstrapped.pdf", width = 297, height = 210, units = "mm")

denoised_projection_plot(SpaCoObject_twin,"Ttr")+ggtitle("bootstrapped denoised")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_ttr_denoised.pdf", width = 297, height = 210, units = "mm")


feature_plot(SpaCoObject,"Cnn2")+ggtitle("original")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_Cnn2_orig.pdf", width = 297, height = 210, units = "mm")

feature_plot(SpaCoObject_twin,"Cnn2")+ggtitle("bootstrapped")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_Cnn2_bootstrapped.pdf", width = 297, height = 210, units = "mm")

denoised_projection_plot(SpaCoObject_twin,"Cnn2")+ggtitle("bootstrapped denoised")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_radius5_Cnn2_denoised.pdf", width = 297, height = 210, units = "mm")

orig <- feature_plot(SpaCoObject,c("Ttr","Cnn2"),ncol = 1)

boots <- feature_plot(SpaCoObject_twin,c("Ttr","Cnn2"),ncol = 1)

deno <- denoised_projection_plot(SpaCoObject_twin,c("Ttr","Cnn2"),ncol = 1)

gene <- boots|deno|orig

#bars <- box|line
#gene | bars 
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_full.pdf", width = 297*3, height = 210, units = "mm")

# Modify the layout to define the ratio
#combined_plot + plot_layout(widths = c(2, 1))

deltas <- list()
bootstraps <- c(1,2,3,4,5,6,7,8,9,10,15,20,25)
for (i in 1:length(bootstraps)){
deltas[[i]] <- distances_bootstrap[[bootstraps[i]]][,1]-distances_bootstrap[[bootstraps[i]]][,2]
}

SVG_brain0 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/DE_genes_spaco_orig.RDS")
SVG_brain0 <- SVG_brain0[order(SVG_brain0$score,decreasing = T),]
#SpaCoObject <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/newSCAMoransIbrainnsim1000.RDS")
top100_SPACO <- rownames(SVG_brain0[SVG_brain0$p.adjust < 0.05,])


names1 <- top100_SPACO[top100_SPACO %in% names(deltas[[1]])]
names2 <- top100_SPACO[top100_SPACO %in% names(deltas[[2]])]
names3 <- top100_SPACO[top100_SPACO %in% names(deltas[[3]])]
names4 <- top100_SPACO[top100_SPACO %in% names(deltas[[4]])]
names5 <- top100_SPACO[top100_SPACO %in% names(deltas[[5]])] 
names6 <- top100_SPACO[top100_SPACO %in% names(deltas[[6]])]
names7 <- top100_SPACO[top100_SPACO %in% names(deltas[[7]])]
names8 <- top100_SPACO[top100_SPACO %in% names(deltas[[8]])]
names9 <- top100_SPACO[top100_SPACO %in% names(deltas[[9]])]
names10 <- top100_SPACO[top100_SPACO %in% names(deltas[[10]])]
names15 <- top100_SPACO[top100_SPACO %in% names(deltas[[11]])] 
names20 <- top100_SPACO[top100_SPACO %in% names(deltas[[12]])]
names25 <- top100_SPACO[top100_SPACO %in% names(deltas[[13]])]


deltas1 <- deltas[[1]][names1,drop=F]
deltas2 <- deltas[[2]][c(names2,"Cnn2"),drop=F]
deltas3 <- deltas[[3]][c(names3,"Cnn2"),drop=F]
deltas4 <- deltas[[4]][names4,drop=F]
deltas5 <- deltas[[5]][names5,drop=F]
deltas6 <- deltas[[6]][names6,drop=F]
deltas7 <- deltas[[7]][names7,drop=F]
deltas8 <- deltas[[8]][names8,drop=F]
deltas9 <- deltas[[9]][names9,drop=F]
deltas10 <- deltas[[10]][names10,drop=F]
deltas15 <- deltas[[11]][names15,drop=F]
deltas20 <- deltas[[12]][names20,drop=F]
deltas25 <- deltas[[13]][names25,drop=F]


df <- stack(list(deltas1 = deltas1, deltas2 = deltas2, deltas3 = deltas3, deltas4 = deltas4, deltas5 = deltas5, deltas6 = deltas6,
                 deltas7 = deltas7, deltas8 = deltas8, deltas9 = deltas9, deltas10 = deltas10, deltas15 = deltas15, deltas20 = deltas20, deltas25 = deltas25))
label_order <- c("deltas1", "deltas2", "deltas3", "deltas4", "deltas5", "deltas6", "deltas7", "deltas8", "deltas9", "deltas10", "deltas15", "deltas20", "deltas25")
#library(ggplot2)
#box <- ggplot(df, aes(x = ind, y = values)) +
#  geom_boxplot() +
  #geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5) +
#  labs(x = "Bootstrap radius", y = "Delta", color = "Group") +
#  scale_x_discrete(
 #   breaks = label_order,
  #  labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25)
 # ) +
 # theme_classic(base_size = 14)
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/denoising_delta_boxplots.pdf", width = 297, height = 210, units = "mm")

# Create data frame function
create_df <- function(data, name) {
  df <- data.frame(values = data, names = names(data), ind = name)
  return(df)
}

# Apply function to each of your vectors
df_list <- lapply(seq_along(deltas), function(i) {
  create_df(deltas[[i]], paste0("deltas", c(1:10,15,20,25)[i]))
})

# Combine all data frames
df <- do.call(rbind, df_list)

# Convert ind to factor for ordered plotting
df$ind <- factor(df$ind, levels = label_order)

# Create boxplot
box <- ggplot(df, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(x = "Bootstrap radius", y = "Delta") +
  scale_x_discrete(breaks = label_order, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25)) +
  theme_classic(base_size = 14)+theme(axis.text = element_text(colour = "black"))

# Highlight "Ttr" values
box <- box +
  geom_jitter(data = subset(df, names == "Ttr"), color = "red", width = 0, size=3)

# Highlight "Cnn2" values
box <- box +
  geom_jitter(data = subset(df, names == "Cnn2"), color = "blue", width = 0,size=3)

print(box)





# Create boxplot
box <- ggplot(df, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(x = "Bootstrap radius", y = "Delta") +
  scale_x_discrete(breaks = label_order, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25)) +
  theme_classic(base_size = 11)

# Highlight "Ttr" values
box <- box +
  geom_jitter(data = subset(df, names == "Ttr"), aes(color = "Ttr"), width = 0.2,size=2)

# Highlight "Cnn2" values
box <- box +
  geom_jitter(data = subset(df, names == "Cnn2"), aes(color = "Cnn2"), width = 0.2,size=2)

# Manual color and legend
box <- box +
  scale_color_manual(values = c("Ttr" = "red", "Cnn2" = "blue"), name = "Gene") +
  theme(legend.position = "right")

print(box)








#####figure2
DE_rad5 <- DE_rad5[order(DE_rad5$score,decreasing = T),]
deltas5_rank <- distances_bootstrap[[5]]
table(rownames(DE_rad5) %in% rownames(deltas5_rank))

SpaCoObject_dist <- SpaCoObject
common_names <- colnames(SpaCoObject_dist@data)[colnames(SpaCoObject_dist@data) %in% colnames(SpaCoObject_twin@data)]

print(length(common_names))

distance_matrix_orig_denoised <- as.matrix(proxy::dist(scale(SpaCoObject_dist@data[,common_names]), SpaCoObject_twin@denoised[,common_names], method = "Euclidean",by_rows = F,pairwise = T))
rownames(distance_matrix_orig_denoised) <- common_names

distance_matrix_orig_bootstrap <- as.matrix(proxy::dist(scale(SpaCoObject_dist@data[,common_names]), scale(SpaCoObject_twin@data[,common_names]), method = "Euclidean",by_rows = F,pairwise = T))
rownames(distance_matrix_orig_bootstrap) <- common_names

distance <- cbind(distance_matrix_orig_denoised,distance_matrix_orig_bootstrap)
colnames(distance) <- c("orig_denoised", "orig_bootstrap")


DE_rad5 <- DE_rad5[common_names,]
DE_rad5 <- DE_rad5[order(DE_rad5$score,decreasing = T),]

library(dplyr)

distance <- as.data.frame(distance)
distance_squared <- distance %>% mutate_all(~(.^2) / 2696)
#View(distance_squared)

distance_squared_ord <- distance_squared[rownames(DE_rad5),]

distance_squared_ord$rank <- c(1:nrow(distance_squared_ord))
library(tidyverse)
library(rcartocolor)

table(DE_rad5$p.adjust < 0.05)
2004-286
df <- distance_squared_ord

# Convert rownames to a column
df <- df %>% rownames_to_column(var = "Gene")

# Gather data to long format
df_long <- df %>%
  gather(key = "variable", value = "value", -Gene, -rank)

# Get the 'Pastel' palette
colors <- carto_pal(name = "Pastel", n = 12)

# Match the colors to your variables
col_mapping <- c("orig_denoised" = colors[1], "orig_bootstrap" = colors[2])

# Plot with loess regression
ggplot(df_long, aes(x = rank, y = value)) +
  geom_point(aes(color = variable)) +
  scale_color_manual(values = col_mapping) +
  geom_smooth(data = df_long %>% filter(variable == "orig_denoised"), 
              method = "loess", se = TRUE, color = "black", fill = "red1") +
  geom_smooth(data = df_long %>% filter(variable == "orig_bootstrap"), 
              method = "loess", se = TRUE, color = "black", fill = "red1") +
  geom_vline(xintercept = 1718, color = "black", linetype = "solid") +
  theme_classic(base_size = 14)
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/radius5_squared_norm_dist.pdf", width = 297, height = 210, units = "mm")


# Identify points to be labeled
label_data <- df_long[df_long$Gene %in% c("Ttr", "Cnn2"), ]

# Add labels to the ggplot
library(ggrepel)

# Create and save the plot
line <- ggplot(df_long, aes(x = rank, y = value)) +
  geom_point(aes(color = variable)) +
  scale_color_manual(values = col_mapping) +
  geom_smooth(data = df_long %>% filter(variable == "orig_denoised"), 
              method = "loess", se = TRUE, color = "black", fill = "red1") +
  geom_smooth(data = df_long %>% filter(variable == "orig_bootstrap"), 
              method = "loess", se = TRUE, color = "black", fill = "red1") +
  geom_point(data = label_data, color = "black", size = 2) +
  geom_text_repel(data = label_data, aes(label = Gene),
                  nudge_x = 540.5, nudge_y = -.15,
                  direction = "both",
                  force = 1,
                  max.iter = 1000,
                  segment.color = "grey") +
  geom_vline(xintercept = 1718, color = "black", linetype = "solid") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",axis.text = element_text(colour = "black"))

box+line

ggsave("/home/rstudio/data_dir/R_scripts/SPACO_generated_figs/fig2ef.pdf", width = 297, height = 210, units = "mm")

