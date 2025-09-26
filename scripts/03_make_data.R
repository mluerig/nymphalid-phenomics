# 00 - header ---------------------------------------------------------------

## remove all vars
rm(list = ls())

## set wd
setwd("D:\\git-repos\\mluerig\\nymphalid-phenomics")
source("scripts\\utils\\utils_r.R")

## "pacman" does automatic package loading
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  # Data manipulation
  "data.table", "lubridate", "stringi", "stringr",
  "googlesheets4","fs",
  
  # Visualization
  "ggplot2", "cowplot", "ggridges", "ggnewscale", "pals",
  "plotly", "ggimage","gridExtra", "egg", "corrplot",
  "scales", "ks", 
  
  # Statistical modeling and multivariate analysis
  "mgcv", "gratia", "FactoMineR", "factoextra",
  "AICcmodavg", "caret", "sda", "fields",
  "umap", "Rtsne", "mclust", "MLmetrics",

  # Phylogenetics and tree visualization
  "phytools"
  
)

pacman::p_load(char = pkgs, install = F)
rm(pkgs)

theme_set(theme_cowplot())
theme_update(
  panel.grid.major = element_line(colour = "grey90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

# 01 - load primary data ---------------------------------------------------------------

## load files
data_emb = fread("data/data_primary/embeddings.csv", stringsAsFactors = T)
data_feat = fread("data/data_primary/features.csv", stringsAsFactors = T)
data_meta = fread("data/data_primary/meta.csv", stringsAsFactors = T)
data_labels = fread("data/data_primary/labels.csv", stringsAsFactors = T)
data_tree = read.nexus("data/data_primary/tree.nex")
feature_key = fread("data/data_primary/feature_key.csv")

## prep rename tree
data_tree_check = data.table(tip_labels = data_tree$tip.label)
data_tree_check[, species := sapply(strsplit(tip_labels, "_"), function(x) {paste(tail(x, 2), collapse = "_")})]
data_tree_check <- data_tree_check[!grepl("_(sp|nsp[0-9]*|nsp)$", species)]
data_tree_check[, species := tolower(species)]
data_tree_check = data_tree_check[!duplicated(species)]
data_tree_check = data_tree_check[order(species)]

## format per-img metadata
data_meta = merge(data_meta, data_labels[, c("species","tribe", "unpalatability")], all.x=T)
data_meta[is.na(unpalatability), unpalatability := "unknown"]
data_meta[is.na(sex) | sex=="", sex := "unknown"]
data_meta[, unpalatability := factor(unpalatability, levels=c("unknown","low","medium","high"))]
data_meta[, species_label := ifelse(species %in% data_tree_check$species & unpalatability!="unknown", "label", "no-label")]

## format per-uuid metadata
data_meta_uuid = unique(data_meta[, c("uuid", "species", "subfamily", "tribe", "unpalatability", "sex", "species_label")])
data_meta_uuid = data_meta_uuid[order(species_label, species, uuid)]
data_meta_uuid[, row := 1:.N,]

## format aggregated metadata
data_meta_agg = unique(data_meta_uuid[, c("species", "species_label","subfamily", "tribe", "unpalatability")])
data_meta_agg = merge(data_meta_agg, data_meta_uuid[, .N, by=.(species)])
setnames(data_meta_agg, "N", "n_specimens")

# 02 - variables / functions  --------------------------------------------------------------

## get colnames
cols_feat = names(data_feat)[substr(names(data_feat),1,4) == "feat"]
cols_emb = names(data_emb)[substr(names(data_emb),1,3) == "emb"]

## colors
colors_unpala = c("unknown"="gray80", "low"="#1b9e77","medium"="#7570b3","high"="#d95f02")
labels_unpala = c("unknown"="NA", "low"="Low","medium"="Medium","high"="High")
colors_moths = c("high"="firebrick1", "low"="darkolivegreen")
labels_moths = c("high"="Aposematism", "low"="Camouflage")

# fam_cols =  c("black", glasbey(15)[5:15])
fam_cols = c(`0` = "black", heliconiinae = "#FF00B6", satyrinae = "#005300", 
             charaxinae = "#FFD300", nymphalinae = "#009FFF", danainae = "#9A4D42", 
             pseudergolinae = "#00FFBE", apaturinae = "#783FC1", biblidinae = "#1F9698", 
             limenitidinae = "#FFACFD", cyrestinae = "#B1CC71", libytheinae = "#F1085C"
)
fam_labs = c("Root", tools::toTitleCase(names(fam_cols)[2:12])) 
names(fam_labs) = names(fam_cols)

# 03 - moths -------------------------------------------------------------------

## load moth dataset
data_moths = fread("data/data_primary/embeddings_moths.csv")
data_moths[, unpalatability := factor(ifelse(strategy=="aposematic", "high", "low"), levels=c("low","high"))]
data_meta_moths = data_moths[, -cols_emb, with=F]
data_meta_moths = data_meta_moths[, c("species",  "strategy",  "mask_name", "unpalatability")]

data_meta_moths_agg = unique(data_meta_moths[, .(species, strategy, unpalatability)])
fwrite(data_meta_moths_agg, "data/data_primary/meta_moths.csv")

# 04 - PCA -----------------------------------------------------------------

## prep
data_pca = merge(data_meta[,c("mask_name", "uuid","species", "unpalatability", "class_dv", "species_label", "sex")], data_emb, by="mask_name")
data_pca = data_pca[order(species_label, species, uuid, class_dv)]

# run
pca1 <- PCA(data_pca[, cols_emb, with=F], ncp=Inf, graph=FALSE, scale.unit=T)
# save(pca1, file="data/analyses_primary/pca1.RData")

## select components that expl. x% variance
pc_cutoff = which(pca1$eig[,3] >= 95)[1]
cols_pc_dims = str_remove(colnames(pca1$ind$coord[,1:pc_cutoff]), "\\.")
pca1_eigenvalues1 = data.table(pca1$eig)
setnames(pca1_eigenvalues1, c("eigenvalue","var_explained", "cum_var_explained"))

## dim reduction butterflies
pca_coords_butterflies = data.table(pca1$ind$coord[,1:pc_cutoff])
setnames(pca_coords_butterflies, cols_pc_dims)
pca_coords_butterflies = cbind(
  data_pca[, -cols_emb, with=F], # class_dv=="dorsal"
  pca_coords_butterflies)

## dim reduction moths
pca_coords_moths = predict(pca1, data_moths[, cols_emb, with=FALSE])
pca_coords_moths = data.table(pca_coords_moths$coord[,1:pc_cutoff])
setnames(pca_coords_moths, cols_pc_dims)
pca_coords_moths = cbind(
  data_moths[, -cols_emb, with=FALSE], 
  pca_coords_moths)

## aggregation
pca_coords_butterflies_agg = pca_coords_butterflies[
  ,lapply(.SD, mean), 
  .SDcols = cols_pc_dims, 
  by = c("species", "unpalatability", "species_label")]
pca_coords_moths_agg = pca_coords_moths[
  ,as.list(unlist(lapply(.SD, mean))),
  .SDcols=cols_pc_dims,
  by=c("species", "unpalatability")]



# 05 - outlier removal ----------------------------------------------------

result_list <- list()
unpala_classes <- c("low", "high")

for (up_class in unpala_classes) {
  
  # Subset data
  data_subset <- pca_coords_butterflies_agg[unpalatability == up_class]

  # Compute cutoff threshold (95% chi-square)
  df_dims <- length(cols_pc_dims)
  mahal_cutoff <- qchisq(0.95, df = df_dims)
  data_subset[, mahal_cutoff := mahal_cutoff]
  
  # Fit GMM with 1 component
  gmm_model <- Mclust(data_subset[, ..cols_pc_dims], G = 1)
  
  # Extract mean and covariance matrix
  mean_vec <- gmm_model$parameters$mean
  cov_mat <- gmm_model$parameters$variance$sigma[,,1]
  
  # Compute Mahalanobis distances
  data_subset[, mahal_dist := apply(.SD, 1, function(x) {
    mahalanobis(x, center = mean_vec, cov = cov_mat)
  }), .SDcols = cols_pc_dims]
  
  # Append to result list
  result_list[[up_class]] <- data_subset
}

# Combine all results
result_gmm_outlier <- rbindlist(result_list)
result_gmm_outlier = result_gmm_outlier[, -cols_pc_dims, with=F]
result_gmm_outlier = result_gmm_outlier[order(mahal_dist)]

## in-distribution vs. out-of-distribution
species_incl <- result_gmm_outlier[mahal_dist < mahal_cutoff]$species 
species_excl <- result_gmm_outlier[mahal_dist > mahal_cutoff]$species

# 06 - LDA --------------------------------------------------------

## subset 
pca_coords_butterflies_agg_incl <- pca_coords_butterflies_agg[species %in% species_incl]
pca_coords_butterflies_agg_incl[, unpalatability := droplevels(unpalatability)]

## split by dv
pca_coords_butterflies_agg_dv = pca_coords_butterflies[
  ,lapply(.SD, mean), 
  .SDcols = cols_pc_dims, 
  by = c("species", "unpalatability", "species_label", "class_dv")]
pca_coords_butterflies_agg_incl_dv <- pca_coords_butterflies_agg_dv[species %in% species_incl]
pca_coords_butterflies_agg_incl_dv[, unpalatability := droplevels(unpalatability)]
pca_coords_butterflies_agg_excl_dv <- pca_coords_butterflies_agg_dv[species %in% species_excl]
pca_coords_butterflies_agg_excl_dv[, unpalatability := droplevels(unpalatability)]

# Stratified train-test split
set.seed(42)
train_species <- pca_coords_butterflies_agg_incl[createDataPartition(
  pca_coords_butterflies_agg_incl[species %in% species_incl]$unpalatability, p = 0.8, list = FALSE),]$species

train_data <- pca_coords_butterflies_agg_incl[species %in% train_species]
train_data_dv <- pca_coords_butterflies_agg_incl_dv[species %in% train_species]
test_data_dv <- pca_coords_butterflies_agg_incl_dv[!species %in% train_species]

# Evaluate on all datasets
pca_coords_moths_agg$class_dv = "dorsal"
datasets <- list(
  "train" = train_data_dv,
  "test1_id" = test_data_dv, 
  "test2_ood" = pca_coords_butterflies_agg_excl_dv, 
  "test3_moths" = pca_coords_moths_agg
)

# Define range of PCs to test
pc_range <- 2:50

## run
results_metrics <- data.table()
results_ld_scores <- data.table()
for (pc in pc_range) {
  print(pc)
  # Fit LDA model
  lda_mod <- lda(as.matrix(train_data[, cols_pc_dims[1:pc], with=FALSE]), 
                    as.factor(train_data$unpalatability)) 
  for (dataset_name in names(datasets)) {
    dataset <- datasets[[dataset_name]]

    # Loop over unique values in `class_dv`
    for (dv_class in unique(dataset$class_dv)) {
      
      
      # Subset dataset for current class_dv
      subset_data <- dataset[class_dv == dv_class]
      
      # Get LDA predictions
      lda_pred <- predict(lda_mod, as.matrix(subset_data[, cols_pc_dims[1:pc], with=FALSE]))
      
      # Compute Precision, Recall, and F1-score (macro-averaged)
      y_pred <- lda_pred$class
      y_true <- subset_data$unpalatability
      conf_matrix <- table(Predicted = y_pred, Actual = y_true)
      acc <- sum(diag(conf_matrix)) / sum(conf_matrix) 
      precision_vals <- diag(conf_matrix) / rowSums(conf_matrix)
      recall_vals <- diag(conf_matrix) / colSums(conf_matrix)
      precision <- mean(precision_vals, na.rm = TRUE) 
      recall <- mean(recall_vals, na.rm = TRUE) 
      f1_score <- 2 * (precision * recall) / (precision + recall)  
      
      # Store performance metrics
      metrics <- data.table(
        pcs = pc, 
        dataset = dataset_name,
        class_dv = dv_class,
        accuracy = acc,
        precision = precision,
        recall = recall,
        F1_score = f1_score
      )
      results_metrics <- rbind(results_metrics, metrics, fill=TRUE)
      
      # Store LD scores
      ld_scores <- data.table(
        pcs = pc,
        dataset = dataset_name,
        class_dv = dv_class,
        species = subset_data$species,
        unpalatability = subset_data$unpalatability,
        lda_pred$x
      )
      results_ld_scores <- rbind(results_ld_scores, ld_scores, fill=TRUE)
    }
  }
}


fwrite(results_metrics, "data/data_secondary/lda_results.csv")

# 07 - LDA predict ---------------------------------------------------------

## use model with 7 components
final_components = 7
lda_mod <- lda(as.matrix(train_data[, cols_pc_dims[1:final_components], with=FALSE]), as.factor(train_data$unpalatability))
# save(lda_mod, file="data/analyses_primary/lda1.RData")

pred_ld_butterflies <- predict(lda_mod, as.matrix(pca_coords_butterflies[, cols_pc_dims[1:final_components], with=FALSE]))
data_ld_scores_butterflies = data.table(
  pca_coords_butterflies[,-cols_pc_dims, with=FALSE],
  LD1 = pred_ld_butterflies$x[,1],
  unpalatability_pred = pred_ld_butterflies$class
)
pred_ld_moths <- predict(lda_mod, as.matrix(pca_coords_moths[, cols_pc_dims[1:final_components], with=FALSE]))
data_ld_scores_moths = data.table(
  pca_coords_moths[,c("species", "activity", "diet", "strategy",  "unpalatability", "image_name", "id")],
  LD1 = pred_ld_moths$x[,1],
  unpalatability_pred = pred_ld_moths$class
)

# pred_ld_moths_agg <- predict(lda_mod, as.matrix(pca_coords_moths_agg[, cols_pc_dims[1:final_components], with=FALSE]))
# data_ld_scores_moths_agg = data.table(
#   pca_coords_moths_agg[,c("species", "unpalatability",  "Dim1", "Dim2")],
#   LD1 = pred_ld_moths_agg$x[,1],
#   unpalatability_pred = pred_ld_moths_agg$class
# )
# data_ld_scores_moths_agg[, fail := unpalatability != as.character(unpalatability_pred)]

# table S2 ----------------------------------------------------------------

for (dataset_name in names(datasets)) {
  cat(dataset_name)
  print(nrow(datasets[[dataset_name]]))
}

data_lda_results = results_metrics[pcs==final_components]
data_lda_results[dataset=="train", N:=215]
data_lda_results[dataset=="test1_id", N:=52]
data_lda_results[dataset=="test2_ood", N:=81]
data_lda_results[dataset=="test3_moths", N:=58]

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Loop over unique values in `class_dv`
  for (dv_class in unique(dataset$class_dv)) {
    
    # Subset dataset for current class_dv
    subset_data <- dataset[class_dv == dv_class]
    
    # Get LDA predictions
    lda_pred <- predict(lda_mod, as.matrix(subset_data[, cols_pc_dims[1:final_components], with=FALSE]))
    
    # Compute Precision, Recall, and F1-score (macro-averaged)
    y_pred <- lda_pred$class
    y_true <- subset_data$unpalatability
    conf_matrix <- table(Predicted = y_pred, Actual = y_true)
    print(paste(dataset_name, "-",dv_class))
    print(conf_matrix)
    
  }
}

# 08 - PCA + supp ---------------------------------------------------------

## identical to PCA1, but with supplementary variables (LD1 + features)

# prep
data_pca2 = merge(data_emb, data_feat)
data_pca2 = merge(data_pca2, data_ld_scores_butterflies[,c("mask_name", "LD1")])

## run
pca2 <- PCA(data_pca2[, c(cols_emb, cols_feat, "LD1"), with=F], ncp=Inf, graph=FALSE, scale.unit=T, quanti.sup=769:925)
# save(pca2, file="data/analyses_primary/pca2.RData")

## pca cols
cols_pc_dims_all = paste0("Dim", 1:176)
cols_pc_dims_apo = paste0("Dim", 1:7)

## eigenvalues
pca_eigenvalues = data.table(pca2$eig)
setnames(pca_eigenvalues, c("eigenvalue","var_explained", "cum_var_explained"))

## coordinates
pca_coords_butterflies = data.table(pca2$ind$coord[,1:176])
setnames(pca_coords_butterflies, cols_pc_dims_all)
pca_coords_butterflies = cbind(mask_name=data_pca2$mask_name, pca_coords_butterflies)
pca_coords_butterflies = merge(data_meta, pca_coords_butterflies, by="mask_name")

## suporting variables
data_pca_quanti_supp = data.table(
  feature_string=rownames(pca2$quanti.sup$coord),
  Dim1=pca2$quanti.sup$coord[,1],
  Dim2=pca2$quanti.sup$coord[,2])
data_pca_quanti_supp = merge(feature_key, data_pca_quanti_supp, all=T)

# 09 - similarity metrics  ------------------------------------------------

# pca_coords_butterflies = fread("data/data_secondary/pca2_coords_butterflies.csv")
# data_ld_scores_butterflies = fread("data/data_secondary/ld_scores_butterflies.csv")

## pca
pca_coords_butterflies_agg = pca_coords_butterflies[
  ,lapply(.SD, mean), 
  .SDcols = cols_pc_dims_all, 
  by = c("species")]
pca_coords_butterflies_agg_dv = pca_coords_butterflies[
  ,lapply(.SD, mean), 
  .SDcols = cols_pc_dims_all, 
  by = c("species", "class_dv")]
pca_coords_butterflies_agg_sex <- pca_coords_butterflies[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = cols_pc_dims_all, 
  by = c("species", "sex")]
pca_coords_butterflies_agg_dv_sex <- pca_coords_butterflies[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = cols_pc_dims_all, 
  by = c("species", "class_dv", "sex")]
pca_coords_butterflies_agg_dv_sex = pca_coords_butterflies_agg_dv_sex[
  species %in% pca_coords_butterflies_agg_dv_sex[, .N, by="species"][N==4]$species]

## LD1
data_ld1_agg <- data_ld_scores_butterflies[
  , lapply(.SD, mean),
  .SDcols = c("LD1"),  
  by = c("species")]
data_ld1_agg_dv <- data_ld_scores_butterflies[
  , lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "class_dv")]
data_ld1_agg_dv = data_ld1_agg_dv[
  species %in% data_ld1_agg_dv[, .N, by="species"][N==2]$species]
data_ld1_agg_sex <- data_ld_scores_butterflies[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "sex")]
data_ld1_agg_sex = data_ld1_agg_sex[
  species %in% data_ld1_agg_sex[, .N, by="species"][N==2]$species]
data_ld1_agg_dv_sex <- data_ld_scores_butterflies[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "class_dv", "sex")]
data_ld1_agg_dv_sex = data_ld1_agg_dv_sex[
  species %in% data_ld1_agg_dv_sex[, .N, by="species"][N==4]$species]


## dv sim
pca_coords_butterflies_agg_dv_wide = dcast(pca_coords_butterflies_agg_dv, species~class_dv, value.var = cols_pc_dims_all)
pca_coords_butterflies_agg_dv_wide[, dv_sim_apo := cosine_sim( 
  as.numeric(.SD[1, paste0("Dim", 1:7, "_dorsal"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:7, "_ventral"), with = FALSE])
), by = c("species")]
pca_coords_butterflies_agg_dv_wide[, dv_sim_all := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:176, "_dorsal"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:176, "_ventral"), with = FALSE])
), by = c("species")]
data_dv_sim_agg = pca_coords_butterflies_agg_dv_wide[,c("species", "dv_sim_apo", "dv_sim_all")]

pca_coords_butterflies_agg_dv_sex_wide = dcast(pca_coords_butterflies_agg_dv_sex, species+sex~class_dv, value.var = cols_pc_dims_all)
pca_coords_butterflies_agg_dv_sex_wide[, dv_sim_apo := cosine_sim( 
  as.numeric(.SD[1, paste0("Dim", 1:7, "_dorsal"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:7, "_ventral"), with = FALSE])
), by = c("species", "sex")]
pca_coords_butterflies_agg_dv_sex_wide[, dv_sim_all := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:176, "_dorsal"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:176, "_ventral"), with = FALSE])
), by = c("species", "sex")]
data_dv_sim_agg_sex = pca_coords_butterflies_agg_dv_sex_wide[,c("species", "sex", "dv_sim_apo", "dv_sim_all")]

## sex sim
pca_coords_butterflies_agg_sex_wide = dcast(pca_coords_butterflies_agg_sex, species~sex, value.var = cols_pc_dims_all)
pca_coords_butterflies_agg_sex_wide[, sex_sim_apo := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:7, "_male"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:7, "_female"), with = FALSE])
), by = c("species")]
pca_coords_butterflies_agg_sex_wide[, sex_sim_all := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:176, "_male"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:176, "_female"), with = FALSE])
), by = c("species")]
data_sex_sim_agg = pca_coords_butterflies_agg_sex_wide[!is.na(sex_sim_apo),c("species", "sex_sim_apo", "sex_sim_all")]

pca_coords_butterflies_agg_sex_dv_wide = dcast(pca_coords_butterflies_agg_dv_sex, species+class_dv~sex, value.var = cols_pc_dims_all)
pca_coords_butterflies_agg_sex_dv_wide[, sex_sim_apo := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:7, "_male"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:7, "_female"), with = FALSE])
), by = c("species", "class_dv")]
pca_coords_butterflies_agg_sex_dv_wide[, sex_sim_all := cosine_sim(
  as.numeric(.SD[1, paste0("Dim", 1:176, "_male"), with = FALSE]),
  as.numeric(.SD[1, paste0("Dim", 1:176, "_female"), with = FALSE])
), by = c("species", "class_dv")]
data_sex_sim_agg_dv = pca_coords_butterflies_agg_sex_dv_wide[,c("species", "class_dv", "sex_sim_apo", "sex_sim_all")]

pca_coords_butterflies_sub = pca_coords_butterflies[!sex=="unknown" ]
pca_coords_butterflies_sub_n = pca_coords_butterflies_sub[, .N, by=c("species", "sex", "class_dv")]
pca_coords_butterflies_sub_n = pca_coords_butterflies_sub_n[N>=3]
pca_coords_butterflies_sub = pca_coords_butterflies_sub[
  species %in% pca_coords_butterflies_sub_n[, .N, by="species"][N==4]$species]
data_mean_sim_agg_dv_sex <- pca_coords_butterflies_sub[, {
  mat_7   <- as.matrix(.SD[, cols_pc_dims_apo, with=F])
  mat_177 <- as.matrix(.SD[, cols_pc_dims_all, with=F])
  .(mean_sim_apo = mean_cosine_sim(mat_7), mean_sim_all = mean_cosine_sim(mat_177))
}, by = c("species", "class_dv", "sex")]


# EXPORT ------------------------------------------------------------------

fwrite(pca_coords_butterflies, "data/data_secondary/pca2_coords_butterflies.csv")
fwrite(pca_eigenvalues, "data/data_secondary/pca2_eigenvalues.csv")
fwrite(data_pca_quanti_supp, "data/data_secondary/pca2_quanti_supp.csv")
fwrite(pca_coords_moths_agg, "data/data_secondary/pca1_coords_moths_agg.csv")

fwrite(data_ld_scores_butterflies, "data/data_secondary/ld_scores_butterflies.csv")
fwrite(data_ld_scores_moths, "data/data_secondary/ld_scores_moths.csv")

fwrite(data_dv_sim_agg, "data/data_secondary/dv_sim_agg.csv")
fwrite(data_dv_sim_agg_sex, "data/data_secondary/dv_sim_agg_sex.csv")
fwrite(data_sex_sim_agg, "data/data_secondary/sex_sim_agg.csv")
fwrite(data_sex_sim_agg_dv, "data/data_secondary/sex_sim_agg_dv.csv")
fwrite(data_mean_sim_agg_dv_sex, "data/data_secondary/mean_sim_agg_dv_sex.csv")


# ROIs --------------------------------------------------------------------

pca_coords_butterflies[, roi_path := paste0( "data_raw//segmentation_masks_clean//all_masks", "//" , species, "//", mask_name)]

## class dv
pca_coords_butterflies_summ = pca_coords_butterflies[, lapply(.SD, mean),.SDcols = cols_pc_dims_all, by=c("species", "class_dv")]
setnames(pca_coords_butterflies_summ, cols_pc_dims_all, paste0(cols_pc_dims_all, "_c"))
pca_coords_butterflies_all = merge(pca_coords_butterflies, pca_coords_butterflies_summ, by=c("species", "class_dv"), all.x=T)
pca_coords_butterflies_all[, dist := sqrt(rowSums((.SD - pca_coords_butterflies_all[
  , paste0(cols_pc_dims_all, "_c"), with=FALSE])^2)), .SDcols = cols_pc_dims_all]
data_sel <- pca_coords_butterflies_all[pca_coords_butterflies_all[, .I[which.min(dist)], by=c("species", "class_dv")]$V1]
data_sel = data_sel[,-c(cols_pc_dims_all, paste0(cols_pc_dims_all, "_c")), with=F]
data_sel[, new_path := file.path(paste0("data_raw//segmentation_masks_clean//centroids_dv//", class_dv, "//", species, ".png"))]
data_sel[, roi_path_exists := file.exists(roi_path)]
data_sel[roi_path_exists == TRUE, {
  mapply(function(from, to) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    fs::file_copy(from, to, overwrite = TRUE)
  }, roi_path, new_path)
}]

## sex
pca_coords_butterflies_summ = pca_coords_butterflies[
  !sex=="unknown" & class_dv=="dorsal", 
  lapply(.SD, mean),.SDcols = cols_pc_dims_all, by=c("species", "sex")]
setnames(pca_coords_butterflies_summ, cols_pc_dims_all, paste0(cols_pc_dims_all, "_c"))
pca_coords_butterflies_all = merge(pca_coords_butterflies, pca_coords_butterflies_summ, by=c("species", "sex"), all.x=T)
pca_coords_butterflies_all[, dist := sqrt(rowSums((.SD - pca_coords_butterflies_all[
  , paste0(cols_pc_dims_all, "_c"), with=FALSE])^2)), .SDcols = cols_pc_dims_all]
data_sel <- pca_coords_butterflies_all[pca_coords_butterflies_all[, .I[which.min(dist)], by=c("species", "sex")]$V1]
data_sel = data_sel[,-c(cols_pc_dims_all, paste0(cols_pc_dims_all, "_c")), with=F]
data_sel[, new_path := file.path(paste0("data_raw//segmentation_masks_clean//centroids_sex//", sex, "//", species, ".png"))]
data_sel[, roi_path_exists := file.exists(roi_path)]
data_sel[roi_path_exists == TRUE, {
  mapply(function(from, to) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    fs::file_copy(from, to, overwrite = TRUE)
  }, roi_path, new_path)
}]

## class dv x sex
pca_coords_butterflies_summ = pca_coords_butterflies[!sex=="unknown", lapply(.SD, mean),.SDcols = cols_pc_dims_all, by=c("species", "class_dv", "sex")]
setnames(pca_coords_butterflies_summ, cols_pc_dims_all, paste0(cols_pc_dims_all, "_c"))
pca_coords_butterflies_all = merge(pca_coords_butterflies, pca_coords_butterflies_summ, by=c("species", "class_dv", "sex"), all.x=T)
pca_coords_butterflies_all[, dist := sqrt(rowSums((.SD - pca_coords_butterflies_all[
  , paste0(cols_pc_dims_all, "_c"), with=FALSE])^2)), .SDcols = cols_pc_dims_all]
data_sel <- pca_coords_butterflies_all[pca_coords_butterflies_all[, .I[which.min(dist)], by=c("species", "class_dv", "sex")]$V1]
data_sel = data_sel[,-c(cols_emb, paste0(cols_emb, "_c")), with=F]
data_sel[, new_path := file.path(paste0("data_raw//segmentation_masks_clean//centroids_dv_sex//",class_dv, "_", sex, "//", species, ".png"))]
data_sel[, roi_path_exists := file.exists(roi_path)]
data_sel[roi_path_exists == TRUE, {
  mapply(function(from, to) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    fs::file_copy(from, to, overwrite = TRUE)
  }, roi_path, new_path)
}]


