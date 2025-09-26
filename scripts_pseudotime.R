# Load required libraries
library(slingshot)
library(xgboost)
library(matrixStats)
library(clusterProfiler)
library(org.Hs.eg.db)  # For human gene annotations
library(ggplot2)
library(scater)
library(nnet)
library(ggeffects)


#first lets infere the pseudotime

# slinghsot
# Convert to SCE
sce <- as.SingleCellExperiment(PD)
# Add clustering if needed (e.g., microglia_type or seurat_clusters)
sce$cluster <- PD$microglia_type  
# run slingshot on the UMAP embedding space guided by cluster labels
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP', start.clus = 'Homeostatic Microglia')
# View pseudotime
head(slingPseudotime(sce))
plotUMAP(sce, colour_by = 'cluster')  # original clusters
plotUMAP(sce, colour_by = 'slingPseudotime_1')  # pseudotime trajectory
plotUMAP(sce, colour_by = 'slingPseudotime_2')  # pseudotime trajectory
plot(reducedDims(sce)$UMAP, col = viridis::viridis(100)[cut(slingPseudotime(sce)[,1], breaks = 100)], pch = 16)
lines(SlingshotDataSet(sce), lwd = 2)

#assign the lineages to the seurat object 
lineage_assignment <- apply(slingPseudotime(sce), 1, function(x) {
  if (is.na(x[1]) & !is.na(x[2])) {
    return("Lineage2")
  } else if (!is.na(x[1]) & is.na(x[2])) {
    return("Lineage1")
  } else if (!is.na(x[1]) & !is.na(x[2])) {
    return("Shared")
  } else {
    return(NA)
  }
})
PD$lineage <- lineage_assignment

# % of cells from each condition in each lineage
table(PD$lineage, PD$Genotype) / rowSums(table(PD$lineage, PD$Genotype))
table(PD$lineage, PD$Treatment) / rowSums(table(PD$lineage, PD$Treatment))
table(PD$lineage, PD$Timepoint) / rowSums(table(PD$lineage, PD$Timepoint))

#some statistics 
# Contingency table
geno_table <- table(PD$lineage, PD$Genotype)
# Chi-squared test
chisq.test(geno_table)

treat_table <- table(PD$lineage, PD$Treatment)
# Chi-squared test
chisq.test(treat_table)

#posthoc for the chisq
install.packages("rcompanion")
library(rcompanion)
pairwiseNominalIndependence(treat_table, method = "fdr")


time_table <- table(PD$lineage, PD$Timepoint)
# Chi-squared test
chisq.test(time_table)

# Set up dataframe
df <- data.frame(
  lineage = PD$lineage,
  genotype = PD$Genotype,
  treatment = PD$Treatment,
  timepoint = PD$Timepoint
)


# Fit multinomial model
model <- multinom(lineage ~ genotype + treatment + timepoint, data = df)
summary(model)
# to get pvalues for the multinomial model through z-scores
z <- summary(model)$coefficients / summary(model)$standard.errors
# p-values
p <- 2 * (1 - pnorm(abs(z)))
round(p, 5)

#to visulaize shifts trends
# Get predicted probabilities for each lineage, by treatment
preds <- ggpredict(model, terms = c("treatment", "genotype", "timepoint"))
plot(preds) + ggtitle("Predicted microglial lineage usage by condition")

preds_time <- ggpredict(model, terms = c("treatment", "timepoint", "genotype"))
plot(preds_time) + 
  ggtitle("Predicted microglial lineage usage by condition") +
  scale_color_manual(values = c("24h" = "#1b9e77", "7d" = "#d95f02"))



#now time for extracting genes associated with pseudotime by ML 
# Step 1: Extract pseudotime and gene expression data from SingleCellExperiment
# Assuming 'sce' is the SingleCellExperiment object from Slingshot
pseudotime_matrix <- slingPseudotime(sce)  # Matrix of pseudotime for all lineages
n_lineages <- ncol(pseudotime_matrix)     # Number of lineages
exprs <- assay(sce, "logcounts")          # Extract log-normalized gene expression
exprs <- t(exprs)                         # Transpose so genes are columns, cells are rows

# Step 2: Preprocess data
# Filter genes with low variance to reduce noise
gene_vars <- colVars(exprs)
exprs_filtered <- exprs[, gene_vars > quantile(gene_vars, 0.8)]  # Keep top 20% variable genes

# Scale the expression data
exprs_scaled <- scale(exprs_filtered)

# Optional: Downsample cells if dataset is too large (e.g., keep 50% of cells)
set.seed(123)  # For reproducibility
if (nrow(exprs_scaled) > 5000) {  # Adjust threshold based on your system's capacity
  sample_cells <- sample(nrow(exprs_scaled), 0.5 * nrow(exprs_scaled))
  exprs_scaled <- exprs_scaled[sample_cells, ]
  pseudotime_matrix <- pseudotime_matrix[sample_cells, ]
}

# Step 3: Process each lineage separately with memory management
for (lineage_idx in 1:n_lineages) {
  cat("\nProcessing Lineage", lineage_idx, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  
  # Extract pseudotime for the current lineage
  pseudotime <- pseudotime_matrix[, lineage_idx]
  
  # Handle NA/NaN/Inf in pseudotime
  valid_cells <- !is.na(pseudotime) & !is.infinite(pseudotime)
  pseudotime <- pseudotime[valid_cells]
  exprs_scaled_lineage <- exprs_scaled[valid_cells, ]
  
  # Verify no NA/NaN/Inf remain
  if (any(is.na(pseudotime)) || any(is.infinite(pseudotime))) {
    stop(paste("Pseudotime for lineage", lineage_idx, "still contains NA or Inf values after filtering."))
  }
  
  # Verify dimensions match
  if (nrow(exprs_scaled_lineage) != length(pseudotime)) {
    stop(paste("Mismatch between number of cells in exprs_scaled and pseudotime for lineage", lineage_idx))
  }
  
  # Step 4: Prepare data for XGBoost
  dtrain <- xgb.DMatrix(data = exprs_scaled_lineage, label = pseudotime)
  
  # Step 5: Set XGBoost parameters (optimized for lower memory usage)
  params <- list(
    objective = "reg:squarederror",  # Regression for pseudotime
    eta = 0.1,                      # Learning rate
    max_depth = 4,                  # Reduced from 6
    subsample = 0.6,                # Reduced from 0.8
    colsample_bytree = 0.6         # Reduced from 0.8
  )
  
  # Step 6: Train XGBoost model
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 50)  # Reduced from 100
  
  # Step 7: Extract feature (gene) importance
  importance <- xgb.importance(model = xgb_model)
  cat("Top 20 Genes by Importance for Lineage", lineage_idx, ":\n")
  print(head(importance, 20))
  
  # Save importance plot to file to avoid memory overload
  plot_title <- paste("Top_20_Genes_Importance_Lineage_", lineage_idx)
  png(filename = paste0(plot_title, ".png"), width = 800, height = 600, res = 100)  # High-res
  xgb.plot.importance(importance, top_n = 20, main = plot_title)
  dev.off()
  
  # Step 8: Optional cross-validation (comment out if too heavy)
  # cv_results <- xgb.cv(params = params, data = dtrain, nrounds = 50, nfold = 3, metrics = "rmse")
  # cat("Cross-validation RMSE for Lineage", lineage_idx, ":\n")
  # print(cv_results)
  
  # Step 9: Optional GO enrichment (comment out if too heavy)
  # top_genes <- importance$Feature[1:50]
  # go_enrich <- enrichGO(
  #   gene = top_genes,
  #   OrgDb = "org.Hs.eg.db",
  #   keyType = "SYMBOL",
  #   ont = "BP",
  #   pAdjustMethod = "BH",
  #   qvalueCutoff = 0.05
  # )
  # if (!is.null(go_enrich)) {
  #   go_plot_title <- paste("GO_Enrichment_Lineage_", lineage_idx)
  #   png(filename = paste0(go_plot_title, ".png"), width = 800, height = 600)
  #   barplot(go_enrich, showCategory = 10, title = go_plot_title)
  #   dev.off()
  # } else {
  #   cat("No significant GO terms found for Lineage", lineage_idx, ".\n")
  # }
  
  # Step 10: Save results and clear memory
  write.csv(importance, paste0("gene_importance_lineage_", lineage_idx, ".csv"))
  saveRDS(xgb_model, paste0("xgboost_model_lineage_", lineage_idx, ".rds"))
  rm(xgb_model, importance, dtrain, exprs_scaled_lineage)
  gc()  # Garbage collection to free memory
}

cat("Analysis completed at", format(Sys.time(), "%H:%M:%S"), "\n")

#evaluate xgboost

n_lineages <- ncol(pseudotime_matrix)

results <- data.frame(Lineage = integer(), RMSE = double(), MAE = double(), R2 = double())

for (lineage_idx in 1:n_lineages) {
  # Load model
  model_path <- paste0("xgboost_model_lineage_", lineage_idx, ".rds")
  xgb_model <- readRDS(model_path)
  
  # Get pseudotime and filter
  pseudotime <- pseudotime_matrix[, lineage_idx]
  valid_cells <- !is.na(pseudotime) & !is.infinite(pseudotime)
  pseudotime <- pseudotime[valid_cells]
  exprs_scaled_lineage <- exprs_scaled[valid_cells, ]
  
  # Recreate DMatrix
  dtrain <- xgb.DMatrix(data = exprs_scaled_lineage, label = pseudotime)
  
  # Predict
  predicted <- predict(xgb_model, newdata = dtrain)
  
  # Compute metrics
  rmse <- sqrt(mean((predicted - pseudotime)^2))
  mae <- mean(abs(predicted - pseudotime))
  r2 <- 1 - sum((pseudotime - predicted)^2) / sum((pseudotime - mean(pseudotime))^2)
  
  results <- rbind(results, data.frame(Lineage = lineage_idx, RMSE = rmse, MAE = mae, R2 = r2))
}

print(results)


ggplot(results, aes(x = factor(Lineage), y = R2)) +
  geom_col(fill = "darkgreen") +
  theme_minimal() +
  labs(title = "XGBoost Performance by Lineage", x = "Lineage", y = expression(R^2))


# Simple K-Fold Cross Validation for XGBoost
library(xgboost)
library(caret)

# Set parameters
k_folds <- 5
n_lineages <- ncol(pseudotime_matrix)
cv_results <- data.frame()

for (lineage_idx in 1:n_lineages) {
  # Get data for this lineage
  pseudotime <- pseudotime_matrix[, lineage_idx]
  valid_cells <- !is.na(pseudotime) & !is.infinite(pseudotime)
  pseudotime <- pseudotime[valid_cells]
  exprs_scaled_lineage <- exprs_scaled[valid_cells, ]
  
  # Create folds
  folds <- createFolds(pseudotime, k = k_folds, list = TRUE)
  
  # Store CV results for this lineage
  fold_results <- data.frame()
  
  for (i in 1:k_folds) {
    # Split data
    test_idx <- folds[[i]]
    train_idx <- setdiff(1:length(pseudotime), test_idx)
    
    # Create DMatrix
    dtrain <- xgb.DMatrix(data = exprs_scaled_lineage[train_idx, ], 
                          label = pseudotime[train_idx])
    dtest <- xgb.DMatrix(data = exprs_scaled_lineage[test_idx, ], 
                         label = pseudotime[test_idx])
    
    # Train model
    xgb_model <- xgb.train(
      data = dtrain,
      nrounds = 100,
      objective = "reg:squarederror",
      verbose = 0
    )
    
    # Predict
    predicted <- predict(xgb_model, dtest)
    actual <- pseudotime[test_idx]
    
    # Calculate metrics
    rmse <- sqrt(mean((predicted - actual)^2))
    mae <- mean(abs(predicted - actual))
    r2 <- 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
    
    # Store fold result
    fold_results <- rbind(fold_results, 
                          data.frame(Fold = i, RMSE = rmse, MAE = mae, R2 = r2))
  }
  
  # Average across folds
  avg_rmse <- mean(fold_results$RMSE)
  avg_mae <- mean(fold_results$MAE)
  avg_r2 <- mean(fold_results$R2)
  
  cv_results <- rbind(cv_results, 
                      data.frame(Lineage = lineage_idx, 
                                 CV_RMSE = avg_rmse, 
                                 CV_MAE = avg_mae, 
                                 CV_R2 = avg_r2))
}

print(cv_results)


# Visualization
library(ggplot2)
library(reshape2)

# 1. Bar plot of CV R² by lineage
ggplot(cv_results, aes(x = factor(Lineage), y = CV_R2)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = round(CV_R2, 3)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Cross-Validation R² by Lineage", 
       x = "Lineage", y = "CV R²") +
  ylim(0, max(cv_results$CV_R2) * 1.1)

# 2. Multiple metrics comparison
cv_long <- melt(cv_results, id.vars = "Lineage", 
                measure.vars = c("CV_RMSE", "CV_MAE", "CV_R2"))
ggplot(cv_long, aes(x = factor(Lineage), y = value, fill = variable)) +
  geom_col(position = "dodge") +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Cross-Validation Metrics by Lineage", 
       x = "Lineage", y = "Value") +
  theme(legend.position = "none")

