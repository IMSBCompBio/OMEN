library(ggplot2)
library(devtools)
library(Seurat)
library(SPARK)
devtools::load_all("/home/rstudio/data_dir/R_scripts/SpaCo_R_kia")

DE_genes_SPACO1 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/DE_genes_SPACO_bootstrap1.RDS")
DE_genes_SPARK1 <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/DE_genes_SPARK_bootstrap1.RDS")
false_positive_names <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/twin_names.RDS")
true_positive_names <- readRDS("/home/rstudio/data_dir/R_scripts/SPACO_paper_data/totwin.RDS")

#false_positive_names <-false_positive_names[1:100]
#true_positive_names <- true_positive_names[1:100]

#radius <- 1
library(ROCR)
# Data Preparation
DE_rad1 <- na.omit(DE_genes_SPACO1[[1]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad1 <- DE_rad1[order(DE_rad1$score,decreasing = T),]

DE_rad5 <- na.omit(DE_genes_SPACO1[[5]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad5 <- DE_rad5[order(DE_rad5$score,decreasing = T),]

DE_rad20 <- na.omit(DE_genes_SPACO1[[20]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad20 <- DE_rad20[order(DE_rad20$score,decreasing = T),]

DE_rad40 <- na.omit(DE_genes_SPACO1[[40]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad40 <- DE_rad40[order(DE_rad40$score,decreasing = T),]


#DE_rad1 <- na.omit(DE_rad1[,"p.adjust",drop=F])
DE_rad1$names <- rownames(DE_rad1)
# Data Preparation
DE_rad_spark1 <- na.omit(DE_genes_SPARK1[[1]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad_spark1 <- DE_rad_spark1[order(DE_rad_spark1$adjustedPval,decreasing = F),]

DE_rad_spark5 <- na.omit(DE_genes_SPARK1[[5]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad_spark5 <- DE_rad_spark5[order(DE_rad_spark5$adjustedPval,decreasing = T),]

DE_rad_spark20 <- na.omit(DE_genes_SPARK1[[20]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad_spark20 <- DE_rad_spark20[order(DE_rad_spark20$adjustedPval,decreasing = T),]

DE_rad_spark40 <- na.omit(DE_genes_SPARK1[[40]][c(true_positive_names,false_positive_names),,drop=F])
DE_rad_spark40 <- DE_rad_spark40[order(DE_rad_spark40$adjustedPval,decreasing = T),]

get_perf_SPACO <- function(df, true_positive_names, false_positive_names) {
  true_labels <- ifelse(row.names(df) %in% true_positive_names, 1, 0)
  true_labels[row.names(df) %in% false_positive_names] <- 0
  pred <- prediction(df$score, true_labels)
  perf <- performance(pred, "tpr", "fpr")
  return(perf)
}
get_perf_SPARK <- function(df, true_positive_names, false_positive_names) {
  true_labels <- ifelse(row.names(df) %in% true_positive_names, 1, 0)
  true_labels[row.names(df) %in% false_positive_names] <- 0
  pred <- prediction(1 - df$adjustedPval, true_labels)
  perf <- performance(pred, "tpr", "fpr")
  return(perf)
}

# Get the performance objects for each data frame
perf_spac1 <- get_perf_SPACO(DE_rad1, true_positive_names, false_positive_names)
perf_spark1 <- get_perf_SPARK(DE_rad_spark1, true_positive_names, false_positive_names)

perf_spac5 <- get_perf_SPACO(DE_rad5, true_positive_names, false_positive_names)
perf_spark5 <- get_perf_SPARK(DE_rad_spark5, true_positive_names, false_positive_names)

perf_spac20 <- get_perf_SPACO(DE_rad20, true_positive_names, false_positive_names)
perf_spark20 <- get_perf_SPARK(DE_rad_spark20, true_positive_names, false_positive_names)

perf_spac40 <- get_perf_SPACO(DE_rad40, true_positive_names, false_positive_names)
perf_spark40 <- get_perf_SPARK(DE_rad_spark40, true_positive_names, false_positive_names)


colors <- c("darkgreen", "lightgreen", "black", "grey", "saddlebrown", "burlywood", "darkblue", "lightblue")
# Plot the ROC curves
plot(perf_spac1, col = colors[1]) # plot first ROC curve in blue
plot(perf_spac5, col = colors[3], add = TRUE)
plot(perf_spac20, col =colors[5], add = TRUE)
plot(perf_spac40, col = colors[7], add = TRUE)


plot(perf_spark1, add = TRUE, col = colors[2]) # add second ROC curve to the same plot in red
plot(perf_spark5, col = colors[4], add = TRUE)
plot(perf_spark20, col = colors[6], add = TRUE)
plot(perf_spark40, col = colors[8], add = TRUE)
# Add a legend
legend("bottomright", legend = c("SPACO r=1","SPACO r=5","SPACO r=20","SPACO r=40", "SPARK r=1","SPARK r=5","SPARK r=20","SPARK r=40"), col =c(colors[1],colors[3],colors[5],colors[7],colors[2],colors[4],colors[6],colors[6]) , lty = 1)


# Convert the performance objects to data frames
perf_spac1_df <- data.frame(FPR = unlist(slot(perf_spac1, "x.values")), TPR = unlist(slot(perf_spac1, "y.values")), Model = "SPACO r = 1")
perf_spac5_df <- data.frame(FPR = unlist(slot(perf_spac5, "x.values")), TPR = unlist(slot(perf_spac5, "y.values")), Model = "SPACO r = 5")
perf_spac20_df <- data.frame(FPR = unlist(slot(perf_spac20, "x.values")), TPR = unlist(slot(perf_spac20, "y.values")), Model = "SPACO r = 20")
perf_spac40_df <- data.frame(FPR = unlist(slot(perf_spac40, "x.values")), TPR = unlist(slot(perf_spac40, "y.values")), Model = "SPACO r = 40")


# Convert the performance objects to data frames
perf_spark1_df <- data.frame(FPR = unlist(slot(perf_spark1, "x.values")), TPR = unlist(slot(perf_spark1, "y.values")), Model = "SPARKX r = 1")
perf_spark5_df <- data.frame(FPR = unlist(slot(perf_spark5, "x.values")), TPR = unlist(slot(perf_spark5, "y.values")), Model = "SPARKX r = 5")
perf_spark20_df <- data.frame(FPR = unlist(slot(perf_spark20, "x.values")), TPR = unlist(slot(perf_spark20, "y.values")), Model = "SPARKX r = 20")
perf_spark40_df <- data.frame(FPR = unlist(slot(perf_spark40, "x.values")), TPR = unlist(slot(perf_spark40, "y.values")), Model = "SPARKX r = 40")


# Combine all data frames
roc_df <- bind_rows(perf_spac1_df, perf_spac5_df, perf_spac20_df, perf_spac40_df,
                    perf_spark1_df, perf_spark5_df, perf_spark20_df, perf_spark40_df)

# Add a new column to indicate whether each row is from a "SPACO" model or a "SPARK" model
roc_df$Method <- ifelse(grepl("SPACO", roc_df$Model), "SPACO", "SPARKX")

library(RColorBrewer)
pal <- brewer.pal(n = 8, name = "Paired")

# Define the color vector
colors <- c("SPACO r = 1" = pal[1], "SPACO r = 5" = pal[3], 
            "SPACO r = 20" = pal[5], "SPACO r = 40" = pal[7], 
            "SPARKX r = 1" = pal[2], "SPARKX r = 5" = pal[4], 
            "SPARKX r = 20" = pal[6], "SPARKX r = 40" = pal[8])

# Define the line types
linetypes <- c("SPACO" = "solid", "SPARKX" = "dashed")

order <- c("SPACO r = 1", "SPACO r = 5", "SPACO r = 20", "SPACO r = 40", 
           "SPARKX r = 1", "SPARKX r = 5", "SPARKX r = 20", "SPARKX r = 40")

# Plot the ROC curves using ggplot2
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model, linetype = Method)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = colors,limits = order) +
  scale_linetype_manual(values = linetypes) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curves") +
  theme_classic(base_size = 10) +
  ylim(0.6,1)

ggsave("/home/rstudio/data_dir/R_scripts/SPACO_generated_figs/fig2b.pdf", width = 297, height = 210, units = "mm")




