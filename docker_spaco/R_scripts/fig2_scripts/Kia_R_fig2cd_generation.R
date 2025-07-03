library(ggplot2)
library(devtools)
library(Seurat)
library(SPARK)
devtools::load_all("/home/rstudio/data_dir/R_scripts/SpaCo_R_kia")
bootstrapdemo_radius1 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/bootstrapdemo_radius1.RDS")
bootstrapdemo_radius5 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/bootstrapdemo_radius5.RDS")
bootstrapdemo_radius10 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/bootstrapdemo_radius10.RDS")
bootstrapdemo_radius20 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/bootstrapdemo_radius20.RDS")
bootstrapdemo_radius50 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/bootstrapdemo_radius50.RDS")

moran_i <- function(x, w) {
  # Calculate required components
  n <- length(x)
  S0 <- sum(w)
  xmean <- mean(x)
  deviations <- x - xmean
  S1 <- sum(deviations^2)
  
  # Calculate spatial lag
  xlag <- w %*% deviations
  
  # Calculate Moran's I
  I <- (n / S0) * (t(deviations) %*% xlag / S1)
  
  return(I)
}

result_mat1 <- apply(bootstrapdemo_radius1@data, 2, function(column) moran_i(column, bootstrapdemo_radius1@neighbours))



result_mat5 <- apply(bootstrapdemo_radius5@data, 2, function(column) moran_i(column, bootstrapdemo_radius5@neighbours))
# result_mat10 <- apply(bootstrapdemo_radius10@data, 2, function(column) moran_i(column, bootstrapdemo_radius10@neighbours))
# result_mat20 <- apply(bootstrapdemo_radius20@data, 2, function(column) moran_i(column, bootstrapdemo_radius20@neighbours))
# result_mat50 <- apply(bootstrapdemo_radius50@data, 2, function(column) moran_i(column, bootstrapdemo_radius50@neighbours))
SpaCoObject <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/newSCAMoransIbrainnsim1000.RDS")
# result_mat0 <- apply(SpaCoObject@data, 2, function(column) moran_i(column, SpaCoObject@neighbours))

#saveRDS(result_mat0,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat0.RDS")
#saveRDS(result_mat1,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat1.RDS")
#saveRDS(result_mat5,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat5.RDS")
#saveRDS(result_mat10,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat10.RDS")
#saveRDS(result_mat20,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat20.RDS")
#saveRDS(result_mat50,"/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat50.RDS")

result_mat0 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat0.RDS")
result_mat1 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat1.RDS")
result_mat50 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat50.RDS")
result_mat20 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat20.RDS")
result_mat10 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat10.RDS")
result_mat5 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/result_mat5.RDS") 



result_mat0 <- as.data.frame(result_mat0)
result_mat0<- result_mat0[order(result_mat0$result_mat0,decreasing = T),,drop=F]

result_mat1 <- as.data.frame(result_mat1)
result_mat1<- result_mat1[order(result_mat1$result_mat1,decreasing = T),,drop=F]

result_mat5 <- as.data.frame(result_mat5)
result_mat5<- result_mat5[order(result_mat5$result_mat5,decreasing = T),,drop=F]

result_mat10 <- as.data.frame(result_mat10)
result_mat10<- result_mat10[order(result_mat10$result_mat10,decreasing = T),,drop=F]

result_mat20 <- as.data.frame(result_mat20)
result_mat20<- result_mat20[order(result_mat20$result_mat20,decreasing = T),,drop=F]

result_mat50 <- as.data.frame(result_mat50)
result_mat50<- result_mat50[order(result_mat50$result_mat50,decreasing = T),,drop=F]


SpaCoObject <- denoise_profiles(SpaCoObject)

bootstrapdemo_radius5 <- RunSCA(bootstrapdemo_radius5,compute_nSpacs = T,nSim = 110) # switch back to 500
bootstrapdemo_radius5 <- denoise_profiles(bootstrapdemo_radius5)

bootstrapdemo_radius20 <- RunSCA(bootstrapdemo_radius20,compute_nSpacs = T,nSim = 110) # switch back to 500
bootstrapdemo_radius20 <- denoise_profiles(bootstrapdemo_radius20)

bootstrapdemo_radius50 <- RunSCA(bootstrapdemo_radius50,compute_nSpacs = T,nSim = 110) # switch back to 500
bootstrapdemo_radius50 <- denoise_profiles(bootstrapdemo_radius50)

g1 <- "Ttr"
g10 <-  "Cnn2"


Ttr_0 <- feature_plot(SpaCoObject,g1) + 
  theme(legend.position = "left")+ coord_cartesian()
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/g1_r0.pdf", width = 297, height = 210, units = "mm")
Siglech_0 <- feature_plot(SpaCoObject,g10)+ 
  theme(legend.position = "left")+ coord_cartesian()

#-----Generating SpaCoObject_twin (maybe by running bootstrapdemo151520.R)-------#
#brain <- readRDS('/home/rstudio/data_dir/R_scripts/SPACO_paper_data/brain.RDS')
SpaCoObject_twin <- RunSCA(bootstrapdemo_radius1, compute_nSpacs = TRUE, nSim = 110) # change this to 500
SpaCoObject_twin <- denoise_profiles(SpaCoObject_twin)
#NOTE: SpaCoObject_twin is literally bootstrapdemo_radius1 without RunSCA(, compute_nSpacs=T, nSims=500)
#     and the denoising profile
#---------------------EXPERIMENTAL WORKAROUND------------------------------------#

Ttr_0_de <- denoised_projection_plot(SpaCoObject_twin,g1)+ 
  theme(legend.position = "left")+ coord_cartesian()

Siglech_0_de <- denoised_projection_plot(SpaCoObject_twin,g10)+ 
  theme(legend.position = "left")+ coord_cartesian()

Ttr_5 <- feature_plot(bootstrapdemo_radius5,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/g1_r5.pdf", width = 297, height = 210, units = "mm")
Siglech_5 <- feature_plot(bootstrapdemo_radius5,g10)+ 
  theme(legend.position = "none")+ coord_cartesian()

Ttr_5_de <- denoised_projection_plot(bootstrapdemo_radius5,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()

Siglech_5_de <- denoised_projection_plot(bootstrapdemo_radius5,g10)+ 
  theme(legend.position = "none")+ coord_cartesian()



Ttr_20 <- feature_plot(bootstrapdemo_radius20,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/g1_r20.pdf", width = 297, height = 210, units = "mm")
Siglech_20 <- feature_plot(bootstrapdemo_radius20,"Cnn2")+ 
  theme(legend.position = "none")+ coord_cartesian()

Ttr_20_de <- denoised_projection_plot(bootstrapdemo_radius20,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()

Siglech_20_de <- denoised_projection_plot(bootstrapdemo_radius20,g10)+ 
  theme(legend.position = "none")+ coord_cartesian()



Ttr_50 <- feature_plot(bootstrapdemo_radius50,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/g1_r40.pdf", width = 297, height = 210, units = "mm")
Siglech_50 <- feature_plot(bootstrapdemo_radius50,"Cnn2")+ 
  theme(legend.position = "none")+ coord_cartesian()


Ttr_50_de <- denoised_projection_plot(bootstrapdemo_radius50,g1)+ 
  theme(legend.position = "none")+ coord_cartesian()

Siglech_50_de <- denoised_projection_plot(bootstrapdemo_radius50,g10)+ 
  theme(legend.position = "none")+ coord_cartesian()


row1 <- Ttr_0|Ttr_5|Ttr_20|Ttr_50
row2 <- Ttr_0_de|Ttr_5_de|Ttr_20_de|Ttr_50_de

row3 <- Siglech_0 |Siglech_5|Siglech_20|Siglech_50
row4 <- Siglech_0_de |Siglech_5_de|Siglech_20_de|Siglech_50_de

bootsis <- row1/row2/row3/row4

spatial_markers_brain <- read.csv2("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/cell_dart_marker_genes_10xvisium_brain.csv")

ltwo <- spatial_markers_brain[spatial_markers_brain$Cluster== "L2/3 IT","Genes"]
lfour <- spatial_markers_brain[spatial_markers_brain$Cluster== "L4","Genes"]
lfive <- spatial_markers_brain[spatial_markers_brain$Cluster=="L5 IT" | spatial_markers_brain$Cluster== "L5 PT","Genes"]
lsix <- spatial_markers_brain[spatial_markers_brain$Cluster== "L6 CT"| spatial_markers_brain$Cluster== "L6 IT"| spatial_markers_brain$Cluster== "L6b","Genes"]
astro <- spatial_markers_brain[spatial_markers_brain$Cluster== "Astro","Genes"]
meis <- spatial_markers_brain[spatial_markers_brain$Cluster== "Meis2","Genes"]
SMC <- spatial_markers_brain[spatial_markers_brain$Cluster== "SMC","Genes"]

spatialgenes <- unique(spatialgenes <- c(ltwo, lfour, lfive, lsix, astro, meis,SMC))
spatialgenes <- spatialgenes[spatialgenes %in% colnames(SpaCoObject@data)]







g36 <- "Nrgn"
g42 <- "Serpine2"
g69 <- "Foxp1"
g40 <- "Dgkb"
g39 <- "Syt6"
g12 <- "Rxrg"
g500 <- "Spp1"


Ttr <- c(result_mat0[g1,],result_mat1[g1,],result_mat5[g1,],result_mat10[g1,],result_mat20[g1,],result_mat50[g1,])
Rxrg <- c(result_mat0[g12,],result_mat1[g12,],result_mat5[g12,],result_mat10[g12,],result_mat20[g12,],result_mat50[g12,])
Syt6 <- c(result_mat0[g39,],result_mat1[g39,],result_mat5[g39,],result_mat10[g39,],result_mat20[g39,],result_mat50[g39,])
Dgkb <- c(result_mat0[g40,],result_mat1[g40,],result_mat5[g40,],result_mat10[g40,],result_mat20[g40,],result_mat50[g40,])
Foxp1 <- c(result_mat0[g69,],result_mat1[g69,],result_mat5[g69,],result_mat10[g69,],result_mat20[g69,],result_mat50[g69,])
Serpine2 <- c(result_mat0[g42,],result_mat1[g42,],result_mat5[g42,],result_mat10[g42,],result_mat20[g42,],result_mat50[g42,])
Nrgn <- c(result_mat0[g36,],result_mat1[g36,],result_mat5[g36,],result_mat10[g36,],result_mat20[g36,],result_mat50[g36,])
Cnn2 <- c(result_mat0[g10,],result_mat1[g10,],result_mat5[g10,],result_mat10[g10,],result_mat20[g10,],result_mat50[g10,])
Spp1 <- c(result_mat0[g500,],result_mat1[g500,],result_mat5[g500,],result_mat10[g500,],result_mat20[g500,],result_mat50[g500,])
# code breaks##
morans_plot <- data.frame(Ttr=Ttr, Cnn2 = Cnn2,Spp1=Spp1,Syt6=Syt6,Dgkb=Dgkb,Foxp1=Foxp1,Serpine2=Serpine2,Nrgn=Nrgn,Rxrg=Rxrg)
rownames(morans_plot) <- c("0","1","5","10","20","40")
morans_plot

median_moran <- median(result_mat50$result_mat50)

library("rcartocolor")
pastel_palette <- carto_pal(9, "Pastel")
pastel_palette[1] <- "blue"
pastel_palette[2] <- "red"

df <- morans_plot
# Make sure that your data frame is in a suitable format
df$radius <- row.names(df) # adds row names as a column named 'index'

# Convert the data from wide to long format
df_long <- df %>% gather(key = "gene", value = "value", -radius)

# Convert 'index' to numeric type for plotting
df_long$radius <- as.numeric(df_long$radius)

# Plot data
noise <- ggplot(df_long, aes(x = radius, y = value, color = gene)) +
  geom_line(linewidth = 1) +
  labs(x = "Bootstrap radius", y = "Moran's I", color = "Gene")+
  geom_hline(yintercept = 0.07615556, linetype = "dashed", color = "black",linewidth = 1) +
  scale_color_manual(values = pastel_palette) +
  theme_classic(base_size = 14)#+theme(legend.position = "none")
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/Morans1_radii.pdf", width = 297, height = 210, units = "mm")

SVG0 <- SVGTest(SpaCoObject)
SVG_vec <- SVG0[,"score",drop=F]
#result_mat0

score <- merge(SVG_vec, result_mat0, by = "row.names", all = TRUE)
colnames(score) <- c("gene", "SPACO_score","Moran's I")
score <- score[order(score$`Moran's I`),]
correlation <- cor(score$`Moran's I`, score$SPACO_score, use = "complete.obs")


ggplot(score, aes(x = SPACO_score, y =`Moran's I` )) +
  geom_point(size = 1) +
  theme_classic(base_size = 14)+
  labs(x = "SPACO score", y = "Moran's I")+
  theme(axis.text = element_text(colour = "black"))
#ggsave("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/corr_morans_spaco_score.pdf", width = 297, height = 210, units = "mm")




# Extract the unique genes
unique_genes <- unique(df_long$gene)

# Generate shades of grey for the other genes
grey_colors <- colorRampPalette(c("black", "darkgrey"))(length(unique_genes) - 2)

# Create the color palette
custom_palette <- c("Ttr" = "blue", "Cnn2" = "red", setNames(grey_colors, setdiff(unique_genes, c("Ttr", "Cnn2"))))

# Plot data
noise <- ggplot(df_long, aes(x = radius, y = value, color = gene)) +
  geom_line(linewidth = 1) +
  labs(x = "Bootstrap radius", y = "Moran's I", color = "Gene") +
  geom_hline(yintercept = 0.07615556, linetype = "dashed", color = "black",linewidth = 1) +
  scale_color_manual(values = custom_palette) +
  theme_classic(base_size = 14)+theme(axis.text = element_text(colour = "black"),legend.position = "right")



library(patchwork)

plot <- bootsis|noise
plot + plot_layout(widths = c(3, 1),heights = c(1/250,1/4))

#saving image
ggsave("/home/rstudio/data_dir/R_scripts/SPACO_generated_figs/fig2cd.pdf", width = 297, height = 210, units = "mm")







