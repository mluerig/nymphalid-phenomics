# 00 - header ---------------------------------------------------------------

## remove all vars
rm(list = ls())

## set wd
setwd("D:\\git-repos\\mluerig\\nymphalid-phenomics")
source("scripts\\utils\\utils_r.R")

## non-cran packages
remotes::install_github("eric-pedersen/MRFtools")
remotes::install_github("Sebastien-Le/YesSiR")

## pacman
if (!require("pacman")) install.packages("pacman")
pkgs <- c(
  # Data manipulation
  "data.table", "lubridate", "stringi", "stringr",
  "googlesheets4","flextable","YesSiR", "matrixStats",
  
  # Visualization
  "ggplot2", "cowplot", "ggridges", "ggnewscale", "pals",
  "plotly", "ggimage","gridExtra", "egg", "corrplot",
  "scales", "grid", "ggrepel","ggtext","ggstar","ggtreeExtra",
  
  # Statistical modeling and multivariate analysis
  "mgcv", "gratia", "FactoMineR", "factoextra",
  "AICcmodavg", "caret", "MRFtools", "lsa",
  "sjPlot", "lmerTest","mclust","stringdist",
  
  # Phylogenetics and tree visualization
  "phytools", "ggtree", "mvMORPH"
)

pacman::p_load(char = pkgs, install = F)
rm(pkgs)

theme_set(theme_cowplot())
theme_update(
  panel.grid.major = element_line(colour = "grey90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

options(scipen = 5)

# 01 - primary data ---------------------------------------------------------------

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
data_meta[, roi_path := paste0( "data_raw//segmentation_masks_clean//all_masks//", species, "//", mask_name)]

## format per-uuid metadata
data_meta_uuid = unique(data_meta[, c("uuid", "species", "subfamily", "tribe", "unpalatability", "sex", "species_label")])
data_meta_uuid = data_meta_uuid[order(species_label, species, uuid)]
data_meta_uuid[, row := 1:.N,]

## format aggregated metadata
data_meta_agg = unique(data_meta_uuid[, c("species", "species_label","subfamily", "tribe", "unpalatability")])
data_meta_agg = merge(data_meta_agg, data_meta_uuid[, .N, by=.(species)])
setnames(data_meta_agg, "N", "n_specimens")

## rois
data_rois_dv = unique(data_meta[, c("species", "species_label","subfamily", "tribe", "unpalatability", "class_dv")])
data_rois_dv[, roi_path := paste0( "data_raw//segmentation_masks_clean//centroids_dv", "//" , class_dv, "//", species, ".png")]
data_rois_dv_sex = unique(data_meta[, c("species", "species_label","subfamily", "tribe", "unpalatability", "class_dv", "sex")])
data_rois_dv_sex = data_rois_dv_sex[!sex == "unknown"]
data_rois_dv_sex[, roi_path := paste0( "data_raw//segmentation_masks_clean//centroids_dv_sex", "//" , class_dv, "_", sex, "//", species, ".png")]

# 02 - secondary data ------------------------------------------------

## pca
pca_coords_butterflies = fread("data/data_secondary/pca2_coords_butterflies.csv", stringsAsFactors = T)
pca_eigenvalues = fread("data/data_secondary/pca2_eigenvalues.csv", stringsAsFactors = T)
pca_quanti_supp = fread("data/data_secondary/pca2_quanti_supp.csv", stringsAsFactors = T)
pca_coords_moths_agg = fread("data/data_secondary/pca1_coords_moths_agg.csv")
pca_coords_moths_agg[, unpalatability := factor(unpalatability, levels=c("low","high"))]
data_meta_moths_agg = fread("data/data_primary/meta_moths.csv")

# LDA
data_lda = fread("data/data_secondary/ld_scores_butterflies.csv", stringsAsFactors = T)
data_lda_moths = fread("data/data_secondary/ld_scores_moths.csv", stringsAsFactors = T)
data_lda_results = fread("data/data_secondary/lda_results.csv", stringsAsFactors = T)

# similarity metrics
data_dv_sim_agg = fread("data/data_secondary/dv_sim_agg.csv", stringsAsFactors = T)
data_dv_sim_agg_sex = fread("data/data_secondary/dv_sim_agg_sex.csv", stringsAsFactors = T)
data_sex_sim_agg = fread("data/data_secondary/sex_sim_agg.csv", stringsAsFactors = T)
data_sex_sim_agg_dv = fread("data/data_secondary/sex_sim_agg_dv.csv", stringsAsFactors = T)
data_mean_sim_agg_dv_sex = fread("data/data_secondary/mean_sim_agg_dv_sex.csv", stringsAsFactors = T)


# 03 - variables / definitions --------------------------------------------

## column names
cols_feat = names(data_feat)[substr(names(data_feat),1,4) == "feat"]
cols_emb = names(data_emb)[substr(names(data_emb),1,3) == "emb"]

## pca cols
cols_pc_dims_all = paste0("Dim", 1:176)
cols_pc_dims_apo = paste0("Dim", 1:7)

## colors
colors_unpala = c(
  "unknown"="gray80", 
  "low"="#1b9e77",
  "medium"="#7570b3",
  "high"="#d95f02")
labels_unpala = c(
  "unknown"="NA", 
  "low"="Low",
  "medium"="Medium",
  "high"="High")
labels_unpala_tree <- c(
  a = "Low",
  b = "Medium",
  c = "High")
colors_unpala_tree <- c(
  a = "#1b9e77",
  b = "#7570b3",
  c = "#d95f02")
colors_moths = c(
  "high"="firebrick1", 
  "low"="darkolivegreen")
labels_moths = c(
  "high"="Aposematism", 
  "low"="Camouflage")

# fam_cols =  c("black", glasbey(15)[5:15])
fam_cols = c(`0` = "black", heliconiinae = "#FF00B6", satyrinae = "#005300", 
  charaxinae = "#FFD300", nymphalinae = "#009FFF", danainae = "#9A4D42", 
  pseudergolinae = "#00FFBE", apaturinae = "#783FC1", biblidinae = "#1F9698", 
  limenitidinae = "#FFACFD", cyrestinae = "#B1CC71", libytheinae = "#F1085C"
)
fam_labs = c("Root", tools::toTitleCase(names(fam_cols)[2:12])) 
names(fam_labs) = names(fam_cols)

# 04 - aggregation / LD1 -------------------------------------------

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
data_ld1_agg <- data_lda[
  , lapply(.SD, mean),
  .SDcols = c("LD1"),  
  by = c("species")]
data_ld1_agg_dv <- data_lda[
  , lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "class_dv")]

data_ld1_agg_sex <- data_lda[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "sex")]
data_ld1_agg_sex = data_ld1_agg_sex[
  species %in% data_ld1_agg_sex[, .N, by="species"][N==2]$species]

data_ld1_agg_dv_sex <- data_lda[
  !sex=="unknown",
  lapply(.SD, mean),
  .SDcols = c("LD1"),
  by = c("species", "class_dv", "sex")]
data_ld1_agg_dv_sex = data_ld1_agg_dv_sex[
  species %in% data_ld1_agg_dv_sex[, .N, by="species"][N==4]$species]

data_lda_all = merge(
  data_meta, 
  data_lda[, .(mask_name, LD1, unpalatability_pred)],
  by="mask_name")

data_lda_summ1 <- data_lda_all[
  species_label=="label",
  {most_common <- .SD[, .N, by = unpalatability_pred][order(-N)][1, unpalatability_pred]
    .(LD1 = mean(LD1, na.rm = TRUE), unpalatability_pred = most_common)},
  by = .(species, subfamily, tribe, species_label, unpalatability)]
data_lda_summ1[!unpalatability=="medium", pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]

data_lda_summ2 <- data_lda_all[,
  {most_common <- .SD[, .N, by = unpalatability_pred][order(-N)][1, unpalatability_pred]
    .(LD1 = mean(LD1, na.rm = TRUE), unpalatability_pred = most_common)},
  by = .(species, subfamily, tribe, species_label, unpalatability)]
data_lda_summ2[!unpalatability=="medium", pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]

data_lda_summ3 <- data_lda_all[,
  {most_common <- .SD[, .N, by = unpalatability_pred][order(-N)][1, unpalatability_pred]
    .(LD1 = mean(LD1, na.rm = TRUE), unpalatability_pred = most_common)},
  by = .(species, subfamily, tribe, species_label, unpalatability, class_dv)]
data_lda_summ3[!unpalatability=="medium", pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]

data_lda_moths_summ1 <- data_lda_moths[,
  {most_common <- .SD[, .N, by = unpalatability_pred][order(-N)][1, unpalatability_pred]
  .(LD1 = mean(LD1, na.rm = TRUE), unpalatability_pred = most_common)},
  by = .(species, unpalatability)]
data_lda_moths_summ1[, pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]


data_lda[species %like% "helic", unique(species)]

# 05 - tree prep ----------------------------------------------------------

## rename tree 
data_tree_all = copy(data_tree)
data_tree_all$tip.label <- data_tree_check$species[
  match(data_tree_all$tip.label, 
        data_tree_check$tip_labels)]
data_tree_all <- drop.tip(data_tree_all, setdiff(
  data_tree_all$tip.label, 
  data_tree_check$species))

## sub 1 = all species in dataset
data_tree_sub1 <- drop.tip(data_tree_all, setdiff(
  data_tree_all$tip.label, 
  data_meta_agg$species))

## sub 2 = sex 
data_tree_sub2 <- drop.tip(data_tree_all, setdiff(
  data_tree_all$tip.label, 
  data_sex_sim_agg$species))

## sub 3 = chemical defense labels
data_tree_sub3 <- drop.tip(data_tree_all, setdiff(
  data_tree_all$tip.label, 
  data_meta_agg[species_label=="label"]$species))

## sub 4 = chemical defense labels / sex
data_tree_sub4 <- drop.tip(data_tree_all, setdiff(
  data_tree_all$tip.label, 
  data_meta_agg[species_label=="label" & species %in% data_sex_sim_agg$species]$species))

# figure 1 - stats --------------------------------------------------------

## this takes a very long time to run (multiple days on a HPC)

# ## across all contexts
# data_mod1 = merge(
#   data_meta_agg[!unpalatability=="unknown"],
#   pca_coords_butterflies_agg)
# data_mod1 = data_mod1[species %in% data_tree_sub1$tip.label]
# tree1 = drop.tip(data_tree_sub1, setdiff(
#   data_tree_sub1$tip.label,
#   data_mod1$species))
# 
# # Subset and construct trait matrix
# trait_mat <- as.matrix(data_mod1[, cols_pc_dims_all, with = FALSE])
# rownames(trait_mat) <- unique(data_mod1[, .(species)])$species
# 
# # Extract unpalatability by species in the same subset
# unpala <- factor(as.character(unique(data_mod1[, .(unpalatability, species)])$unpalatability))
# mv_dataset <- list(y = trait_mat, unpala = unpala)
# 
# # Fit model
# model <- mvgls(y ~ unpala, data = mv_dataset, tree = tree1, model = "lambda", method = "LOO")
# anova <- manova.gls(model, nperm = 1000, test = "Pillai", verbose = TRUE)
# 
# # Save output
# fname_prefix <- paste0("data/analyses_secondary/mvmorph_models/mvgls")
# save(model, file = paste0(fname_prefix, "_model.RData"))
# save(anova, file = paste0(fname_prefix, "_anova.RData"))

# ## across D/V
# data_mod2 = merge(
#   data_meta_agg[!unpalatability=="unknown"],
#   pca_coords_butterflies_agg_dv)
# data_mod2 = data_mod2[species %in% data_tree_sub1$tip.label]
# tree1 = drop.tip(data_tree_sub1, setdiff(
#   data_tree_sub1$tip.label,
#   data_mod2$species))
# v = "ventral"
# for (v in c("dorsal", "ventral")) {
#   
#   message("Fitting model for: ", v)
#   
#   # Subset and construct trait matrix
#   trait_mat <- as.matrix(data_mod2[class_dv == v , cols_pc_dims_all, with = FALSE])
#   rownames(trait_mat) <- unique(data_mod2[class_dv == v, .(species)])$species
#   
#   # Extract unpalatability by species in the same subset
#   unpala <- factor(as.character(
#     unique(data_mod2[class_dv == v, .(unpalatability, species)])$unpalatability
#   ))
#   mv_dataset <- list(y = trait_mat, unpala = unpala)
#   
#   # Fit model
#   model <- mvgls(y ~ unpala, data = mv_dataset, tree = tree1, model = "lambda", method = "LOO")
#   anova <- manova.gls(model, nperm = 1000, test = "Pillai", verbose = TRUE)
#   
#   # Save output
#   fname_prefix <- paste0("data/analyses_secondary/mvmorph_models/mvgls_", v)
#   save(model, file = paste0(fname_prefix, "_model.RData"))
#   save(anova, file = paste0(fname_prefix, "_anova.RData"))
# }
# 

# ## across D/V and sex
# data_mod3 = merge(
#   data_meta_agg[!unpalatability=="unknown"],
#   pca_coords_butterflies_agg_dv_sex)
# data_mod3 = data_mod3[species %in% data_tree_sub1$tip.label]
# tree2 = drop.tip(data_tree_sub1, setdiff(
#   data_tree_sub1$tip.label,
#   data_mod3$species))
# 
# for (s in  c("female", "male")) {
#   for (v in c("dorsal", "ventral")) {
#     
#     message("Fitting model for: ", v, " - ", s)
#     
#     # Subset and construct trait matrix
#     trait_mat <- as.matrix(data_mod3[class_dv == v & sex == s, cols_pc_dims_all, with = FALSE])
#     rownames(trait_mat) <- unique(data_mod3[class_dv == v & sex == s, .(species)])$species
#     
#     # Extract unpalatability by species in the same subset
#     unpala <- factor(as.character(
#       unique(data_mod3[class_dv == v & sex == s, .(unpalatability, species)])$unpalatability
#     ))
#     
#     mv_dataset <- list(y = trait_mat, unpala = unpala)
#     
#     # Fit model
#     model <- mvgls(y ~ unpala, data = mv_dataset, tree = tree2, model = "lambda", method = "LOO")
#     anova <- manova.gls(model, nperm = 1000, test = "Pillai", verbose = TRUE)
#     
#     # Save output
#     fname_prefix <- paste0("data/analyses_secondary/mvgls_", v, "_", s)
#     save(model, file = paste0(fname_prefix, "_model.RData"))
#     save(anova, file = paste0(fname_prefix, "_anova.RData"))
#   }
# }

## load previously run analyses 
results_list_model = list() 
results_list_anova = list() 
for(file_name in list.files("data/analyses_secondary/mvmorph_models")){
  model_name = str_split_fixed(str_split(file_name,"\\.")[[1]][1], "_", n=2)[[2]]
  load(file.path("data/analyses_secondary/mvmorph_models", file_name))
  if(str_detect(model_name, "model")){
    results_list_model[[model_name]] = model
  } else {
    results_list_anova[[model_name]] = anova
  }
}

## across all contexts
mod1 = results_list_model["model"]
res1 = results_list_anova["anova"]
effectsize(res1[[1]])

# ## across D/V
# mod1 = results_list_model["dorsal_model"]
# mod2 = results_list_model["ventral_model"]
# res1 = results_list_anova["dorsal_anova"]
# res2 = results_list_anova["ventral_anova"]
# effectsize(res1[[1]])
# effectsize(res2[[1]])

# ## across D/V and sex
# mod1 = results_list_model["dorsal_female_model"]
# mod2 = results_list_model["dorsal_male_model"]
# mod3 = results_list_model["ventral_female_model"]
# mod4 = results_list_model["ventral_male_model"]
# res1 = results_list_anova["dorsal_female_anova"]
# res2 = results_list_anova["dorsal_male_anova"]
# res3 = results_list_anova["ventral_female_anova"]
# res4 = results_list_anova["ventral_male_anova"]
# effectsize(res1[[1]])
# effectsize(res2[[1]])
# effectsize(res3[[1]])
# effectsize(res4[[1]])

# figure 1/S1 - prep ---------------------------------------------------------

## KDE on aggregated coords
data_plot = merge(data_meta_agg, pca_coords_butterflies_agg_dv)

## contours
data_KDE = data_plot[species_label=="label", {
  contour_dt <- kde_2d(Dim1,Dim2,levels = c(0.25, 0.5, 0.75),
                       n_grid = 100, bw_adjust = 1)
}, by = .(unpalatability, class_dv)]

## grid
data_KDE_grid <- data_plot[species_label == "label", {
  kde_grid_2d(Dim1, Dim2, n_grid = 100, bw_adjust = 1, expand=0.2)
}, by = .(unpalatability, class_dv)]
probs <- c(0.05, 0.1, 0.25, 0.5, 0.75)
thr_long <- data_KDE_grid[
  , {
    o   <- order(density)
    cum <- cumsum(pmass[o])
    tot <- sum(pmass)
    cut <- sapply(probs, function(p) density[o][ which.max(cum >= p*tot) ])
    data.table(prob = paste0("p",substr(as.character(probs), 3, 20)), cutoff = as.numeric(cut))
  },
  by = .(unpalatability, class_dv)
]
thr_wide <- dcast(thr_long, unpalatability + class_dv ~ prob, value.var = "cutoff")
data_KDE_grid <- thr_wide[data_KDE_grid, on=.(unpalatability, class_dv)]

## rois
roi_manual_selectíon = c(
  "actinote_surima", "anaeomorpha_splendida", "antillea_pelops", 
  "antirrhea_philoctetes", "archaeoprepona_demophon", "bematistes_alcinoe", 
  "brenthis_ino", "callicore_tolima", "chlosyne_janais", "danaus_plexippus", 
  "higginsius_fasciata", "hyposcada_anchialia", "libythea_celtis", 
  "microtia_elva", "morpho_helenor", "parthenos_sylvia", "phyciodes_cocyta", 
  "poladryas_arachne", "pseudacraea_poggei", "pteronymia_teresita", 
  "sephisa_dichroa", "thaumantis_klugius", "yramea_cytheris")
data_roi_subset = data_rois_dv[species %in% roi_manual_selectíon]
data_roi_subset = merge(
  data_plot[, c("species","class_dv", "Dim1", "Dim2")], 
  data_roi_subset, by=c("species", "class_dv"))

## scaling factor LD1
loadings_factor = 25

## supplementary vars
data_quanti_supp = pca_quanti_supp
data_quanti_supp[!feature_string=="LD1", Dim1_s := Dim1 * (loadings_factor + 5)]
data_quanti_supp[!feature_string=="LD1", Dim2_s := Dim2 * (loadings_factor + 5)]
data_quanti_supp[feature_string=="LD1", Dim1_s := Dim1 * loadings_factor]
data_quanti_supp[feature_string=="LD1", Dim2_s := Dim2 * loadings_factor]
data_quanti_supp_sel = data_quanti_supp[feature_string %in% c(
  "feat_color_moments_hue_var",
  "feat_color_moments_hue_uniformity",
  "feat_color_moments_red_var",
  "feat_color_moments_red_mean",
  "feat_color_moments_lum_var",
  "feat_color_hist_lab_binvar",
  "feat_color_dft_lum_meanmag"
)]
data_quanti_supp_sel[, feature_formatted2 := sub("[(]", "\n(", feature_formatted)]
data_quanti_supp_sel[feature_string == "feat_color_moments_hue_var", `:=` (Dim1_s=4, Dim2_s=5)]
data_quanti_supp_sel[feature_string == "feat_color_moments_hue_uniformity",`:=` (Dim1_s=-7, Dim2_s=4)]
data_quanti_supp_sel[feature_string == "feat_color_moments_red_mean",`:=` (Dim1_s=-10,Dim2_s=-17)]
data_quanti_supp_sel[feature_string == "feat_color_dft_lum_meanmag",`:=` (Dim1_s=-7,Dim2_s=-13)]
data_quanti_supp_sel[feature_string == "feat_color_moments_lum_var",`:=` (Dim1_s=15)]
data_quanti_supp_sel[, c("feature_string", "Dim1","Dim2","Dim1_s","Dim2_s")]

## insets
data_plot2 = merge(data_meta_agg, data_ld1_agg_dv)
inset_settings = list(
  coord_cartesian(ylim=c(0,0.5), xlim=c(-2.5, 6)),
  theme_cowplot(), 
  geom_density(aes(x = LD1, fill=unpalatability), alpha=0.5), 
  scale_fill_manual(values=colors_unpala, labels=labels_unpala, guide="none"),
  scale_y_continuous(expand=c(0,0)),
  theme(
    axis.text = element_text(size=10),
    axis.title.x = element_text(size=10, margin = margin(b=-2,t=2,unit="pt")),
    axis.title.y = element_text(size=10, margin = margin(l=-2,r=2,unit="pt")),
    plot.background = element_rect(colour = "black", fill = "white", size=1),
    plot.margin = margin(t=10,b=6,l=6,r=6, unit="pt")
  ) 
)
inset_d <-
  ggplot(data_plot2[class_dv=="dorsal" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings
inset_v <-
  ggplot(data_plot2[class_dv=="ventral" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings

# figure 1 - plot ---------------------------------------------------------


## figure
p1 =
  ggplot() + 
  facet_wrap(~class_dv, labeller = as_labeller(
    c("dorsal"="Dorsum", "ventral"="Ventrum")), scales="free_x") + 
  xlab(paste("Dim1 - Variance explained: ", round(pca_eigenvalues[1,]$var_explained,2), "%")) +
  ylab(paste("Dim2 - Variance explained: ", round(pca_eigenvalues[2,]$var_explained,2), "%")) +
  scale_x_continuous(breaks=seq(-100,100,10)) + 
  scale_y_continuous(breaks=seq(-100,100,10)) + 
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  
  ## points without labels
  geom_point(
    data=data_plot[!species_label=="label"],
    aes(x=Dim1, y=Dim2, color=unpalatability), size=1, inherit.aes = FALSE) +
  scale_color_manual("Chemical defense:", values=colors_unpala, labels=labels_unpala, drop=FALSE,
                     guide = guide_legend(keywidth = unit(1.75 , "cm"), order=1,
                                          override.aes = list(size=c(2,3,3,3)))) +
  
  ## kde grid
  geom_raster(
    data=data_KDE_grid[unpalatability == "low" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "medium" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "high" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  scale_alpha_continuous(
    range = c(0, 1),
    trans = "sqrt",
  ) +  
  
  ## kde contours
  geom_path(
    data=data_KDE,
    aes(x=x, y=y, color=unpalatability, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dashed")) +
  
  ## points with labels
  geom_point(
    data=data_plot[species_label=="label"],
    aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), size=2, pch=21) +
  scale_fill_manual("Chemical defense:", values=colors_unpala, labels=labels_unpala, guide="none") +
  
  ## pictograms
  geom_image(
    data=data_roi_subset,
    aes(image=roi_path, x=Dim1, y=Dim2), size=0.1, by="height",
    image_fun=function(img) {img = magick::image_background(img, "none")}) +

  ## points for examples
  geom_point(
    data=data_roi_subset,
    aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), size=3, pch=21) +
  
  ## loading arrow
  new_scale_fill() +
  geom_segment(data=pca_quanti_supp[feature_string=="LD1"],
               aes(x = 0, y = 0, xend =  Dim1*loadings_factor,
                   yend =  Dim2*loadings_factor), size=1,
               arrow = arrow(angle = 30, length = unit(2, "mm"),
                             ends = "last", type = "closed")) +
  geom_label(data=pca_quanti_supp[feature_string=="LD1"],
             aes(x =  Dim1*loadings_factor, y =  Dim2*loadings_factor,
                 label=feature_string, fill = feature_string), size=5, show.legend = FALSE) +
  scale_fill_manual(values=unname(colors_unpala["high"])) +
  
  ## select features
  geom_label(
    data=data_quanti_supp_sel,
    aes(x = Dim1_s, y = Dim2_s, label = feature_formatted2), 
    size=3.5, label.padding = unit(3.5, "pt"), lineheight = 0.8, alpha=0.85) +

  ## model results
  geom_label(data=data.table(
    x=-19, 
    y=-29, 
    class_dv="dorsal", 
    label=paste("Pillai effect size:", 0.144, "\nP-value: <0.001")),
   aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
   label.size=NA, alpha=0.85) +
  geom_label(data=data.table(
    x=-25, 
    y=-29, 
    class_dv="ventral", 
    label=paste("Pillai effect size:", 0.206, "\nP-value: <0.001")),
    aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
    label.size=NA, alpha=0.85) +
  
  ## theme settings
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    aspect.ratio = 1,
  )

p1 = ggdraw() +
  draw_plot(p1) +
  draw_plot(inset_d, x = 0.3, y = 0.1125, width = 0.21, height = 0.225) +
  draw_plot(inset_v, x = 0.775, y = 0.1125, width = 0.21, height = 0.225)

p_leg =
  ggplot(data_plot) + 
  
  geom_point(
    data=data_plot,
    aes(x=Dim1, y=Dim2, color=unpalatability), size=1) +
  scale_color_manual(
    "Chemical defense", values=colors_unpala, labels=labels_unpala, 
    guide = guide_legend(
      order=1, keywidth = unit(40 , "pt"),
      override.aes = list(size=c(3,3,3,1)), title.position="top")) +
  ## KD lines
  geom_path(
    data=data_KDE, 
    aes(x=x, y=y, color=unpalatability, linetype=factor(level), 
        group=interaction(level_idx,unpalatability, group)), linewidth=1.25) +
  
  scale_linetype_manual(
    "Spatial density quantiles (KDE Dim1&2)", values=c("solid","longdash","dashed"),
    labels=c("25%","50%","75%"), 
    guide = guide_legend(order=2, keywidth = unit(65 , "pt"), title.position="top")) +
  
  theme(
    legend.margin = margin(l = 50, unit = "pt"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "left")

leg <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]
pg_final = cowplot::plot_grid(p1, leg, nrow=2, rel_heights = c(0.9,0.1))

width = 4000
height = 2500
ggsave(paste0("figures/figure1.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")

# figure 2/S5 - prep ----------------------------------------------------------------

## prep tree 
tree = copy(data_tree_sub3)

## family labels
species_by_family <- data_meta_agg[, .(species_list = list(species)), by = "subfamily"]
species_by_family <- setNames(species_by_family$species_list, species_by_family$subfamily)
tree = groupOTU(tree, species_by_family, group_name = "subfamily")

## plot data
data_plot1 = dcast(data_ld1_agg_dv, species~class_dv, value.var = "LD1")
setnames(data_plot1, c("dorsal", "ventral"), c("LD1_d", "LD1_v"))
data_plot1 = merge(data_plot1, data_meta_agg[species_label == "label"])

## traits
trait_values1 <- data.frame(data_plot1[, c("LD1_v", "LD1_d")])
rownames(trait_values1) <- data_plot1$species

## labels unpala
unpalatability_labs <- data.frame(unpalatability=data_plot1$unpalatability)
unpalatability_labs$unpalatability = factor(
  unpalatability_labs$unpalatability, levels=c("low","medium","high"), labels=c("a","b","c"))
rownames(unpalatability_labs) <- data_plot1$species

## labels pred
pred_success_labs <- data.frame(pred_success=data_lda_summ1$pred_success)
pred_success_labs$pred_success = factor(
  pred_success_labs$pred_success, levels=c(T,F), labels=c("a","b"))
rownames(pred_success_labs) <- data_lda_summ1$species

## images
data_tree1 = data.table(ggtree(tree, layout = "circular")$data)
data_tree1[, species := label]
data_tree1 = data_tree1[isTip == TRUE,]
data_tree1 = data_tree1[order(y)]

## select species / images
y_sel = c(6, 17, 25, 37, 48, 63, 78, 93, 107, 123, 137, 154, 168, 182, 191, 200, 
          213, 228, 243, 258, 275, 293, 308, 324, 338, 349)
data_tree1 = data_tree1[y %in% y_sel]
data_tree1[, species_f := str_to_sentence(str_replace(as.character(species), "_", " "))]

data_tree1[, N := factor(1:.N)]
data_tree1[, y_img := y]
data_tree1[N==2, y_img := y_img-1]
data_tree1[N==4, y_img := y_img-1]
data_tree1[N==11, y_img := y_img+1]
data_tree1[N==14, y_img := y_img-2]
data_tree1[N==15, y_img := y_img-2]
data_tree1[N==16, y_img := y_img-1]
data_tree1[N==17, y_img := y_img-1]
data_tree1[N==19, y_img := y_img+1]

data_tree1[, roi_path_d := paste0("data_raw/segmentation_masks_clean/centroids_dv/", "dorsal", "/", species, ".png")]
data_tree1[, roi_path_v := paste0("data_raw/segmentation_masks_clean/centroids_dv/", "ventral", "/", species, ".png")]

data_tree1[species=="stichophthalma_howqua", species_f := "S. howqua"]
data_tree1[species=="archaeoprepona_demophon", species_f := "A. demophon"]
data_tree1[species=="pardopsis_punctatissima", species_f := "P. punctatissima"]
data_tree1[species=="etcheverrius_chiliensis", species_f := "E. chiliensis"]
data_tree1[species=="stygionympha_vigilans", species_f := "S. vigilans"]

## difference
data_ld1_agg_dv_sex_wide = dcast(data_ld1_agg_dv_sex, species+class_dv~sex, value.var = "LD1")
data_ld1_agg_dv_sex_wide[, diff := male-female]
data_ld1_agg_dv_sex_wide[, sign := ifelse(diff > 0, "female_higher", "male_higher")]
data_ld1_agg_dv_sex_wide[, sign := factor(sign, levels=c("male_higher", "female_higher"))]
data_ld1_agg_dv_sex_wide[, diff_abs := abs(diff)]
data_ld1_agg_dv_sex_wide = merge(data_meta_agg, data_ld1_agg_dv_sex_wide)
data_ld1_agg_dv_sex_wide[, subfamily_f := str_to_sentence(subfamily)]

## counts / props
data_ld1_diff_n = data_ld1_agg_dv_sex_wide[, .(N_total = .N), by = .(subfamily,class_dv)]
data_ld1_diff_n = merge(
  data_ld1_diff_n, 
  data_ld1_agg_dv_sex_wide[, .(N_sign = .N), by = .(subfamily, sign, class_dv)])
data_ld1_diff_n[, prop := N_sign / N_total]
data_ld1_diff_n[, subfamily_f := str_to_sentence(subfamily)]

## subfamilies to exclude 
subfamily_excl = data_ld1_diff_n[N_total<10, as.character(unique(subfamily))]

## null dist
set.seed(42)
n_sim <- 1000
obs_diffs <- data_ld1_agg_dv_sex_wide$diff
null_mat <- replicate(n_sim, {
  # for each row, randomly flip sign with 50% chance
  flip <- sample(c(-1, 1), length(obs_diffs), replace = TRUE)
  obs_diffs * flip
})
null_dist <- data.table(diff=as.numeric(null_mat))
q = quantile(null_dist$diff, c(0.1, 0.9))

## outside null
data_ld1_agg_dv_sex_wide[species %in% data_ld1_agg_dv_sex_wide[diff <= q[1] | diff >= q[2],species], out := T]
data_ld1_agg_dv_sex_wide[!species %in% data_ld1_agg_dv_sex_wide[diff <= q[1] | diff >= q[2],species], out := F]
data_ld1_agg_dv_sex_wide_sub = data_ld1_agg_dv_sex_wide[out==T]
data_ld1_diff_sub_n = data_ld1_agg_dv_sex_wide_sub[, .(N_total = .N), by = .(subfamily,class_dv)]
data_ld1_diff_sub_n = merge(
  data_ld1_diff_sub_n, 
  data_ld1_agg_dv_sex_wide_sub[, .(N_sign = .N), by = .(subfamily, sign, class_dv)])
data_ld1_diff_sub_n[, prop := N_sign / N_total]
data_ld1_diff_sub_n[, subfamily_f := str_to_sentence(subfamily)]

data_plot2 = data_ld1_agg_dv_sex_wide[class_dv=="dorsal" & species %in% tree$tip.label]
data_plot2[, sign_plot := sign]
data_plot2[out==F, sign_plot := NA]

# figure 2 - plot  --------------------------------------------------------

## tree
p1 =
  ggtree(tree, aes(colour = subfamily), 
         size=0.5, layout = "circular")  + 
  geom_segment(data = data_tree1,
               aes(x=x, xend=x+75, y=y, yend=y, color=subfamily),size=0.5) +
  geom_tiplab(data=data_tree1, aes(x=x+35, y=y_img, image=roi_path_v), geom="image",
              size=0.08, inherit.aes = FALSE, angle=0) +
  geom_tiplab(data=data_tree1, aes(x=x+75, y=y_img, image=roi_path_d), geom="image",
              size=0.1, inherit.aes = FALSE, angle=0) +
  scale_color_manual(values=fam_cols, "Subfamily") +
  geom_point(data=data_tree1, aes(x=x+55, y=y),
             size=6, pch=21, fill="white") +
  geom_text(data=data_tree1, aes(x=x+55, y=y, label=N),
            size=3.5, color="black") 

p1 <- p1 + new_scale_fill()
p1 = gheatmap(p1, unpalatability_labs, offset=-3,  width=0.05, colnames = F) + 
  scale_fill_manual(values=colors_unpala_tree) 

p1 <- p1 + new_scale_fill()
p1 = gheatmap(p1, pred_success_labs, offset=2,  width=0.05, colnames = F) + 
  scale_fill_manual(values=c("green","red", NA)) 

p1 <- p1 + new_scale_fill()
p1 =  gheatmap(p1, trait_values1, width=0.1, colnames = F, offset=7) +
  scale_fill_viridis_c(option="inferno", oob = scales::oob_squish,limits = c(-2.5,6.5)) +
  
  theme(legend.position = "none",
        plot.margin = margin(t=-35,b=0,l=-10,r=0)
        )

p2 =
  ggtree(tree, aes(colour = subfamily), 
         size=0.25, layout = "circular")  + 
  scale_color_manual(values=fam_cols, "Subfamily", guide="none") + 
  geom_fruit(
    data=data_plot2, 
    geom=geom_bar, mapping=aes(x=diff, y=species, fill=sign_plot), orientation="y",stat="identity",
    offset=0.4, pwidth=0.4)  +
  scale_fill_manual(values=c("sienna3", "deepskyblue")) 

p2 <- p2 + new_scale_fill()
p2 = gheatmap(p2, unpalatability_labs, offset=-3,  width=0.05, colnames = F) + 
  scale_fill_manual(values=colors_unpala_tree) +
  
  theme(legend.position = "none",
        plot.margin = margin(t=-100,b=-50,l=-70,r=-50),
        panel.background = element_rect(fill = "transparent", color = NA)
        )

p3 =
  ggplot(data_ld1_diff_sub_n[!subfamily %in% subfamily_excl & class_dv=="dorsal"]) + 
  facet_wrap(~class_dv, labeller = as_labeller(c("dorsal"="Dorsum", "ventral"="Ventrum"))) +
  ylab("Proportion") + 
  geom_bar(aes(x=subfamily_f, y=prop, fill=sign), position = "stack", stat = "identity") +
  geom_hline(yintercept = 0.5) +
  scale_fill_manual(values=c("sienna3", "deepskyblue")) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=-20,b=0,l=-20,r=70)
  )



p_leg <-
  ggplot(data_tree1, aes(x=1, y=N, fill=N)) +
  geom_tile() +  
  scale_fill_manual(
    name = NULL,
    values = rep("white", nrow(data_tree1)),  # invisible boxes
    labels = paste0(data_tree1$N,". ", data_tree1$species_f)
  ) +
  guides(fill = guide_legend(ncol=2, byrow = TRUE)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.1, 'cm'),
        legend.text=element_text(size=8),
        legend.margin = margin(t = -75),
        legend.box.margin = margin(t = 10, r = 5, b = 10, l = 5),
        # legend.box.background = element_rect(color = "black", linetype = "solid", size = 0.5)
  )


## legend
p_leg1 =
  ggplot(data_plot1) +
  
  geom_line(
    aes(x=1,y=1,color=subfamily)) +
  scale_color_manual(values=fam_cols, "Subfamily", guide=guide_legend(
    order=1, override.aes = list(linewidth=1), keywidth = unit(15, "pt"),
    title.position = "top", ncol=2, byrow=T),
    labels=fam_labs) +
  
  new_scale_fill() + 
  geom_bar(
    data=unpalatability_labs,
    aes(x=1,y=1,fill=unpalatability), stat="identity") +
  scale_fill_manual(
    "Chemical defense\n(labels)", 
    values=colors_unpala_tree, 
    labels=labels_unpala_tree, 
    guide=guide_legend(order=2, title.position = "top")) +

  new_scale_fill() + 
  geom_bar(
    data=pred_success_labs,
    aes(x=1,y=1,fill=pred_success), stat="identity") +
  scale_fill_manual(
    "Prediction of chemical defense\n(LD-classifier)", 
    values=c("green","red", NA), 
    labels=c("Correct", "Incorrect", "(NA [Medium])"), 
    guide=guide_legend(order=3, title.position = "top")) +
  
  new_scale_color() + 
  geom_point(
    data=data_ld1_agg_dv,
    aes(x=1,y=1,color=LD1)) +
  scale_color_viridis_c(
    "Aposematic color patterns (LD1) \n(outer=Dorsum, inner=Ventrum)",,
    option="inferno", 
    oob = scales::oob_squish,
    limits = c(-2.5,6.5),
    guide=guide_colorbar(
      order=4, keyheight=unit(0.5, "cm"), keywidth=unit(6,"cm"), 
      title.position = "top", direction = "horizontal")) +
  
  new_scale_fill() + 
  geom_tile(data=data_tree1,
            aes(x=1, y=1, fill=N)) +  
  scale_fill_manual(
    "  Selected species",
    values = rep("white", nrow(data_tree1)),  # invisible boxes
    labels = paste0("(",data_tree1$N,") ", data_tree1$species_f), 
    guide=guide_legend(
      order=5, title.position = "top",
      keywidth=unit(0,"pt"),
      theme = theme(
        legend.margin = margin(l=-8, unit="pt"),
        legend.text = element_text(size=10),
        # legend.key.size = unit(20,"pt"),
        legend.spacing.x = unit(0, "pt")
      ),
    )
  ) +
  
  new_scale_fill() + 
  geom_bar(
    data=data_ld1_diff_n,
    aes(x=1,y=1,fill=sign ), stat="identity") +
  scale_fill_manual(
    "\n\nSexual differences in\naposematic color patterns",
    values=c("deepskyblue", "sienna3"),
    labels=c("Higher male LD1", "Higher female LD1"),
    guide=guide_legend(order=6)) + 
  
  theme(
    legend.position = "right",
    legend.margin = margin(b = 20, unit="pt"),
    )

leg <- get_plot_component(p_leg1, 'guide-box', return_all = TRUE)[[1]]

pg1 = cowplot::plot_grid(p2, p3, ncol=2, rel_widths = c(0.6,0.4), labels=c("B","C"), vjust=-3)
pg2 = cowplot::plot_grid(p1, pg1, nrow=2, rel_heights = c(0.75,0.25), labels=c("A",NA))
pg_final = cowplot::plot_grid(pg2, leg, ncol=2, rel_widths = c(0.75,0.25))

width=4000
height=4000
ggsave(paste0("figures/figure2.png"), pg_final,bg="white",
       width=width, height=height, units = c("px"))


# figure 3/S5 - stats ---------------------------------------------------------------

data_mod1 = merge(data_meta_agg, data_ld1_agg)
data_mod1 = merge(data_mod1, data_dv_sim_agg, by=c("species"))
data_mod1_phylo = data_mod1[species %in% data_tree_sub1$tip.label]
data_tree_wahl2009_mod1 <- drop.tip(data_tree_sub1, setdiff(data_tree_sub1$tip.label, data_mod1_phylo$species))
phylo_penalty_mod1 <- mrf_penalty(data_tree_wahl2009_mod1)

data_mod2 = merge(data_meta_agg, data_ld1_agg_dv)
data_mod2 = merge(data_mod2, data_sex_sim_agg_dv, by=c("species", "class_dv"))
data_mod2_phylo = data_mod2[species %in% data_tree_sub1$tip.label]
data_tree_wahl2009_mod2 <- drop.tip(data_tree_sub1, setdiff(data_tree_sub1$tip.label, data_mod2_phylo$species))
phylo_penalty_mod2 <- mrf_penalty(data_tree_wahl2009_mod2)

data_mod3 = merge(data_meta_agg, data_ld1_agg_dv_sex)
data_mod3 = merge(data_mod3, data_mean_sim_agg_dv_sex, by=c("species", "class_dv", "sex"))
data_mod3_phylo = data_mod3[species %in% data_tree_sub1$tip.label]
data_tree_wahl2009_mod3 <- drop.tip(data_tree_sub1, setdiff(data_tree_sub1$tip.label, data_mod3_phylo$species))
phylo_penalty_mod3 <- mrf_penalty(data_tree_wahl2009_mod3)

knots = 10

# ## dv sim - non-phylo
# mod_1a = gam(dv_sim_all ~
#                s(LD1, k = knots),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod1)
# save(mod_1a, file="data/analyses_secondary/gam_regressions/mod_1a.RData")
load(file="data/analyses_secondary/gam_regressions/mod_1a.RData")
capture.output(anova(mod_1a), file = "tables/mod_1a.txt")
# mod_1b = gam(dv_sim_apo ~
#                s(LD1, k = knots),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod1)
# save(mod_1b, file="data/analyses_secondary/gam_regressions/mod_1b.RData")
load(file="data/analyses_secondary/gam_regressions/mod_1b.RData")
capture.output(anova(mod_1b), file = "tables/mod_1b.txt")

# ## dv sim - phylo
# mod_1a_phylo = gam(dv_sim_all ~
#                 s(LD1, k = knots) +
#                 s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod1),
#                 k=as.integer(length(unique(data_mod1_phylo$species))/2)),
#             method = "REML",
#             family = betar(link = "logit"),
#             data=data_mod1_phylo)
# save(mod_1a_phylo, file="data/analyses_secondary/gam_regressions/mod_1a_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_1a_phylo.RData")
capture.output(anova(mod_1a_phylo), file = "tables/mod_1a_phylo.txt")
# mod_1b_phylo = gam(dv_sim_apo ~
#                 s(LD1, k = knots) +
#                 s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod1),
#                 k=as.integer(length(unique(data_mod1_phylo$species))/2)),
#             method = "REML",
#             family = betar(link = "logit"),
#             data=data_mod1_phylo)
# save(mod_1b_phylo, file="data/analyses_secondary/gam_regressions/mod_1b_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_1b_phylo.RData")
capture.output(anova(mod_1b_phylo), file = "tables/mod_1b_phylo.txt")

## sex sim - non-phylo
# mod_2a = gam(sex_sim_all ~
#                s(LD1, k = knots) +
#                s(LD1, class_dv, k=knots, bs = "sz"),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod2)
# save(mod_2a, file="data/analyses_secondary/gam_regressions/mod_2a.RData")
load(file="data/analyses_secondary/gam_regressions/mod_2a.RData")
capture.output(anova(mod_2a), file = "tables/mod_2a.txt")
# mod_2b = gam(sex_sim_apo ~
#                s(LD1, k = knots) +
#                s(LD1, class_dv, k=knots, bs = "sz"),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod2)
# save(mod_2b, file="data/analyses_secondary/gam_regressions/mod_2b.RData")
load(file="data/analyses_secondary/gam_regressions/mod_2b.RData")
capture.output(anova(mod_2b), file = "tables/mod_2b.txt")

## sex sim - phylo
# mod_2a_phylo = gam(sex_sim_all ~
#                 s(LD1, k = knots) +
#                 s(LD1, class_dv, k=knots, bs = "sz") +
#                 s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod2),
#                 k=as.integer(length(unique(data_mod2_phylo$species))/2)),
#                method = "REML",
#                family = betar(link = "logit"),
#                data=data_mod2_phylo)
# save(mod_2a_phylo, file="data/analyses_secondary/gam_regressions/mod_2a_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_2a_phylo.RData")
capture.output(anova(mod_2a_phylo), file = "tables/mod_2a_phylo.txt")
# mod_2b_phylo = gam(sex_sim_apo ~
#                      s(LD1, k = knots) +
#                      s(LD1, class_dv, k=knots, bs = "sz") +
#                      s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod2),
#                        k=as.integer(length(unique(data_mod2_phylo$species))/2)),
#                    method = "REML",
#                    family = betar(link = "logit"),
#                    data=data_mod2_phylo)
# save(mod_2b_phylo, file="data/analyses_secondary/gam_regressions/mod_2b_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_2b_phylo.RData")
capture.output(anova(mod_2b_phylo), file = "tables/mod_2b_phylo.txt")

# ## mean sim - non-phylo
# mod_3a = gam(mean_sim_all ~
#                s(LD1, k = knots) +
#                s(LD1, sex, k=knots, bs = "sz") +
#                s(LD1, class_dv, k=knots, bs = "sz") +
#                s(LD1, sex, class_dv, k=knots, bs = "sz"),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod3)
# save(mod_3a, file="data/analyses_secondary/gam_regressions/mod_3a.RData")
load(file="data/analyses_secondary/gam_regressions/mod_3a.RData")
capture.output(anova(mod_3a), file = "tables/mod_3a.txt")

# mod_3b = gam(mean_sim_apo ~
#                s(LD1, k = knots) +
#                s(LD1, sex, k=knots, bs = "sz") +
#                s(LD1, class_dv, k=knots, bs = "sz") +
#                s(LD1, sex, class_dv, k=knots, bs = "sz"),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod3)
# save(mod_3b, file="data/analyses_secondary/gam_regressions/mod_3b.RData")
load(file="data/analyses_secondary/gam_regressions/mod_3b.RData")
capture.output(anova(mod_3b), file = "tables/mod_3b.txt")

# ## mean sim - phylo
# mod_3a_phylo = gam(mean_sim_all ~
#                 s(LD1, k = knots) +
#                 s(LD1, sex, k=knots, bs = "sz") +
#                 s(LD1, class_dv, k=knots, bs = "sz") +
#                 s(LD1, sex, class_dv, k=knots, bs = "sz") +
#                 s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod3),
#                   k=as.integer(length(unique(data_mod3_phylo$species))/2)),
#             method = "REML",
#             family = betar(link = "logit"),
#             data=data_mod3_phylo)
# save(mod_3a_phylo, file="data/analyses_secondary/gam_regressions/mod_3a_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_3a_phylo.RData")
capture.output(anova(mod_3a_phylo), file = "tables/mod_3a_phylo.txt")

# mod_3b_phylo = gam(mean_sim_apo ~
#                s(LD1, k = knots) +
#                s(LD1, sex, k=knots, bs = "sz") +
#                s(LD1, class_dv, k=knots, bs = "sz") +
#                s(LD1, sex, class_dv, k=knots, bs = "sz") +
#                s(species, bs="mrf", xt=list(penalty=phylo_penalty_mod3),
#                  k=as.integer(length(unique(data_mod3_phylo$species))/2)),
#              method = "REML",
#              family = betar(link = "logit"),
#              data=data_mod3_phylo)
# save(mod_3b_phylo, file="data/analyses_secondary/gam_regressions/mod_3b_phylo.RData")
load(file="data/analyses_secondary/gam_regressions/mod_3b_phylo.RData")
capture.output(anova(mod_3b_phylo), file = "tables/mod_3b_phylo.txt")


# figure 3/S5 - prep1 ------------------------------------------------------

## mod1 - predictions
mod_1a_pred_glob = pred_gam_newd(
  mod=mod_1a, varx = "LD1", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(LD1,sex)"))
mod_1a_pred_part = pred_gam_newd(
  mod=mod_1a, varx = "LD1", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)
mod_1b_pred_glob = pred_gam_newd(
  mod=mod_1b, varx = "LD1", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(LD1,sex)"))
mod_1b_pred_part = pred_gam_newd(
  mod=mod_1b, varx = "LD1", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)

mod_1a_phylo_pred_glob = pred_gam_newd(
  mod=mod_1a_phylo, varx = "LD1",  
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(species)", "s(LD1,sex)"))
mod_1a_phylo_pred_part = pred_gam_newd(
  mod=mod_1a_phylo, varx = "LD1",  
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")
mod_1b_phylo_pred_glob = pred_gam_newd(
  mod=mod_1b_phylo, varx = "LD1",  
  extra_x = c(0, 0), extra_y = c(0, 0), se=T, 
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(species)", "s(LD1,sex)"))
mod_1b_phylo_pred_part = pred_gam_newd(
  mod=mod_1b_phylo, varx = "LD1",  
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")

## mod2 - predictions
mod_2a_pred_glob = pred_gam_newd(
  mod=mod_2a, varx = "LD1", var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(LD1,class_dv)"))
mod_2a_pred_part = pred_gam_newd(
  mod=mod_2a, varx = "LD1", var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)
mod_2b_pred_glob = pred_gam_newd(
  mod=mod_2b, varx = "LD1", var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F,
  exclude=c("s(LD1,class_dv)"))
mod_2b_pred_part = pred_gam_newd(
  mod=mod_2b, varx = "LD1", var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)

mod_2a_phylo_pred_glob = pred_gam_newd(
  mod=mod_2a_phylo, varx = "LD1",  var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c("s(species)", "s(LD1,class_dv)"))
mod_2a_phylo_pred_part = pred_gam_newd(
  mod=mod_2a_phylo, varx = "LD1",  var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")
mod_2b_phylo_pred_glob = pred_gam_newd(
  mod=mod_2b_phylo, varx = "LD1",  var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c("s(species)", "s(LD1,class_dv)"))
mod_2b_phylo_pred_part = pred_gam_newd(
  mod=mod_2b_phylo, varx = "LD1",  var_cat="class_dv", 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")

## mod3 - predictions
mod_3a_pred_glob = pred_gam_newd(
  mod=mod_3a, varx = "LD1", var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F, exclude=c(
    "s(LD1,sex)","s(LD1,class_dv)","s(LD1,sex,class_dv)"
  ))
mod_3a_pred_part = pred_gam_newd(
  mod=mod_3a, varx = "LD1", var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)
mod_3b_pred_glob = pred_gam_newd(
  mod=mod_3b, varx = "LD1", var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F, exclude=c(
    "s(LD1,sex)","s(LD1,class_dv)","s(LD1,sex,class_dv)"
  ))
mod_3b_pred_part = pred_gam_newd(
  mod=mod_3b, varx = "LD1", var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val, 
  type="response", newdata_check = F)


mod_3a_phylo_pred_glob = pred_gam_newd(
  mod=mod_3a_phylo, varx = "LD1",  var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c(
    "s(species)", "s(LD1,sex)","s(LD1,class_dv)","s(LD1,sex,class_dv)"
  ))
mod_3a_phylo_pred_part = pred_gam_newd(
  mod=mod_3a_phylo, varx = "LD1",  var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")
mod_3b_phylo_pred_glob = pred_gam_newd(
  mod=mod_3b_phylo, varx = "LD1",  var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c(
    "s(species)", "s(LD1,sex)","s(LD1,class_dv)","s(LD1,sex,class_dv)"
  ))
mod_3b_phylo_pred_part = pred_gam_newd(
  mod=mod_3b_phylo, varx = "LD1",  var_cat=c("class_dv","sex"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude="s(species)")

# figure 3/S5 - prep2 ------------------------------------------------------

## p1
data_mod1_agg <- data_mod1[
  ,lapply(.SD, mean),
  .SDcols = c("LD1", "dv_sim_all", "dv_sim_apo"),  
  by = c("species", "unpalatability")]
data_mod1_phylo_agg <- data_mod1_phylo[
  ,lapply(.SD, mean),
  .SDcols = c("LD1", "dv_sim_all", "dv_sim_apo"),  
  by = c("species", "unpalatability")]
p1_grid_data = data_mod1_agg[order(LD1, -dv_sim_all)]
p1_grid_data[, roi_path_d := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/dorsal_male/", species, ".png")]
p1_grid_data[, roi_path_v := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/ventral_male/", species, ".png")]
p1_grid_data = p1_grid_data[species %in% c(
    ## high
    "limenitis_archippus", "oleria_onega",
    ## medium
    "castilia_castilla", 
    ## low
    "morpho_helenor", "baeotus_beotus", "archaeoprepona_demophon",
    ## unknown
    "hamadryas_arete","limenitis_lorquini"
  )]
p1_grid_data[, grid_x := rep(0:1, each = 4)]
p1_grid_data[, grid_y := rep(3:0, times = 2)]
p1_grid_data[, number:= as.numeric(factor(species, levels=unique(species)))]
p1_grid_data[, label := str_replace(as.character(species), "_", " ")]
p1_grid_data[, label := paste0("  ", str_to_sentence(label))]
p1_grid_data[species=="archaeoprepona_demophon", label:=paste0("  ", "A. demophon")]


## p2
data_mod2_agg <- data_mod2[
  ,lapply(.SD, mean),
  .SDcols = c("LD1", "sex_sim_all", "sex_sim_apo"),  
  by = c("species", "unpalatability")]
data_mod2_phylo_agg <- data_mod2_phylo[
  ,lapply(.SD, mean),
  .SDcols = c("LD1",  "sex_sim_all", "sex_sim_apo"),  
  by = c("species", "unpalatability")]
p2_grid_data = data_mod2_agg[order(LD1, -sex_sim_all)]
p2_grid_data[, roi_path_f := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/dorsal_female/", species, ".png")]
p2_grid_data[, roi_path_m := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/dorsal_male/", species, ".png")]
p2_grid_data = p2_grid_data[species %in% c(
  ## high
  "mechanitis_polymnia","heliconius_egeria",
  ## medium
  "haematera_pyrame","philaethria_wernickei",
  ## low
  "epiphile_orea", "eulaceura_osteria", #"tanaecia_julii",#
  ## unknown
  "hypolimnas_salmacis","acraea_andromacha"
)]
p2_grid_data[, grid_x := rep(0:1, each = 4)]
p2_grid_data[, grid_y := rep(3:0, times = 2)]
p2_grid_data[, number:= as.numeric(factor(species, levels=unique(species)))]
p2_grid_data[, label := str_replace(as.character(species), "_", " ")]
p2_grid_data[, label := paste0("  ", str_to_sentence(label))]
p2_grid_data[species=="gnathotriche_exclamationis", label:=paste0("  ", "G. exclamationis")]

data_mod2_phylo_agg[species %like% "heliconi"]


## p3
data_mod3_agg <- data_mod3[
  ,lapply(.SD, mean),
  .SDcols = c("LD1", "mean_sim_all", "mean_sim_apo"),  
  by = c("species", "unpalatability")]
data_mod3_phylo_agg <- data_mod3_phylo[
  ,lapply(.SD, mean),
  .SDcols = c("LD1", "mean_sim_all", "mean_sim_apo"),  
  by = c("species", "unpalatability")]
p3_grid_data = data_mod3_agg[order(LD1, -mean_sim_all)]
p3_grid_data = p3_grid_data[species %in% c(
  ## high
  "parantica_aspasia",
  ## medium
  "telenassa_teletusa","boloria_napaea",
  ## low
  "precis_octavia" ,  "cymothoe_caenis", 
  ## unknown
  "erebia_epipsodea", "coenonympha_tullia","heliconius_charithonia" 
)]
pairs_by_sex <- data_lda[
  species %in% p3_grid_data$species & class_dv=="dorsal" & !sex=="unknown", {
    ld1_diff <- as.matrix(dist(LD1))
    idx <- which(ld1_diff == max(ld1_diff), arr.ind = TRUE)[1, ]
    .(mask_name1 = mask_name[idx[1]], mask_name2 = mask_name[idx[2]], 
      LD1_1 = LD1[idx[1]], LD1_2 = LD1[idx[2]], 
      diff = abs(LD1[idx[1]] - LD1[idx[2]]))
  }, 
  by = .(species, sex)
]
pairs_dist <- pairs_by_sex[
  pairs_by_sex[
    , .I[which.max(diff)], 
    by = species
  ]$V1
]

p3_grid_data = merge(p3_grid_data, data_meta[mask_name %in% pairs_dist$mask_name1, .(species, sex, roi_path)])
setnames(p3_grid_data, "roi_path",  "roi_path_h")
p3_grid_data = merge(p3_grid_data, data_meta[mask_name %in% pairs_dist$mask_name2, .(species, roi_path)])
setnames(p3_grid_data, "roi_path",  "roi_path_l")
p3_grid_data = p3_grid_data[order(LD1, -mean_sim_all)]
p3_grid_data[, grid_x := rep(0:1, each = 4)]
p3_grid_data[, grid_y := rep(3:0, times = 2)]
p3_grid_data[, number:= as.numeric(factor(species, levels=unique(species)))]
p3_grid_data[, label := str_replace(as.character(species), "_", " ")]
p3_grid_data[, label := paste0("  ", str_to_sentence(label))]
p3_grid_data[, label := paste0(label, " (", substr(sex, 1,1), ")")]
p3_grid_data[species=="heliconius_charithonia", label := paste0("  ","H. charithonia (f)")]
p3_grid_data[species=="coenonympha_tullia", label := paste0("  ", "C. tullia (m)")]

## plotting
plot_elements = list(
  scale_color_manual(values=colors_unpala, drop=T), 
  scale_fill_manual(values=colors_unpala, drop=T),
  theme(legend.position = "none")
)


# figure 3 - plot ----------------------------------------------------------------

p1 =
  ggplot() + ggtitle("A") + 
  ylab("Dorso-ventral similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim=c(0.3,1)) +
  geom_point(data=data_mod1_phylo_agg[unpalatability=="unknown"],
             aes(x=LD1, y=dv_sim_all, color=unpalatability, label=species), size=1.2, pch=16) + 
  geom_point(data=data_mod1_phylo_agg[!unpalatability=="unknown"],
             aes(x=LD1, y=dv_sim_all, fill=unpalatability), pch=21, size=2) +
  geom_ribbon(data=mod_1a_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_1a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_ribbon(data=mod_1b_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_1b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  geom_point(data=p1_grid_data, aes(x=LD1, y=dv_sim_all, fill=unpalatability), size=6, pch=21) +
  geom_text(data=p1_grid_data, aes(x=LD1, y=dv_sim_all, label=number), size=3) +
  plot_elements + theme(
    axis.title.x = element_blank()
  )

y_lab_top = 3.65
y_lab_size = 4
text_size = 3.5
p1_g =
  ggplot(p1_grid_data) + 
  theme_void() + coord_cartesian(clip = "off") +
  geom_image(aes(image = roi_path_d, x = grid_x-0.225, y = grid_y), size = 0.22, by="width") +
  geom_image(aes(image = roi_path_v, x = grid_x+0.225, y = grid_y), size = 0.22, by="width") +
  geom_textbox(aes(x = grid_x, y = grid_y-0.38, label = label), size=text_size, width=unit(4.6, "cm"),
               height = unit(0.5, "cm"), box.r = unit(7, "pt"), lineheight = 1.4, box.padding = unit(0.2, "lines"),
               fontface = "italic", halign = 0.5, valign=1) +
  geom_point(data=p1_grid_data, aes(x=grid_x-0.42, y=grid_y-0.38, fill=unpalatability), size=6.4, pch=21, stroke = 0.4) + 
  geom_text(data=p1_grid_data, aes(x=grid_x-0.42, y=grid_y-0.38, label=number), size=3) + 
  scale_fill_manual(values=colors_unpala, guide = "none", drop=T) +  
  annotate("text", x = 0-0.225, y = y_lab_top, label = "Dorsum", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 0+0.225, y = y_lab_top, label = "Ventrum", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1-0.225, y = y_lab_top, label = "Dorsum", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1+0.225, y = y_lab_top, label = "Ventrum", size = y_lab_size, fontface = "bold") +
  theme(
    plot.margin = margin(t = 10, r = 40, b = 20, l = 10, unit = "pt"),
    plot.title = element_text(vjust=5, hjust=-0.35, size=16, face="bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )

p2 =
  ggplot() + ggtitle("B") + 
  ylab("Inter-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim=c(0.7,1)) +
  geom_point(data=data_mod2_phylo_agg[unpalatability=="unknown"],
             aes(x=LD1, y=sex_sim_all, color=unpalatability, label=species), size=1.2, pch=16) + 
  geom_point(data=data_mod2_phylo_agg[!unpalatability=="unknown"],
             aes(x=LD1, y=sex_sim_all, fill=unpalatability), pch=21, size=2) +
  geom_ribbon(data=mod_2a_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_2a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_ribbon(data=mod_2b_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_2b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  geom_point(data=p2_grid_data, aes(x=LD1, y=sex_sim_all, fill=unpalatability), size=6, pch=21) +
  geom_text(data=p2_grid_data, aes(x=LD1, y=sex_sim_all, label=number), size=3) +
  annotate("text", x = 4, y = 0.7, label = "(29 datapoints < 0.7)", size = 4) +
  plot_elements + theme(
    axis.title.x = element_blank()
  )


y_lab_top = 3.5
p2_g =
  ggplot(p2_grid_data) +
  theme_void() + coord_cartesian(clip = "off") +
  geom_image(aes(image = roi_path_f, x = grid_x-0.225, y = grid_y), size = 0.22, by="width") +
  geom_image(aes(image = roi_path_m, x = grid_x+0.225, y = grid_y), size = 0.22, by="width") +
  geom_textbox(aes(x = grid_x, y = grid_y-0.35, label = label), size=text_size, width=unit(4.6, "cm"),
               height = unit(0.5, "cm"), box.r = unit(7, "pt"), lineheight = 1.4, box.padding = unit(0.2, "lines"),
               fontface = "italic", halign = 0.5, valign=1) +
  geom_point(data=p2_grid_data, aes(x=grid_x-0.42, y=grid_y-0.35, fill=unpalatability), size=6.4, pch=21, stroke = 0.4) + 
  geom_text(data=p2_grid_data, aes(x=grid_x-0.42, y=grid_y-0.35, label=number), size=3) + 
  scale_fill_manual(values=colors_unpala, guide = "none", drop=T) +  
  annotate("text", x = 0-0.225, y = y_lab_top, label = "Female", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 0+0.225, y = y_lab_top, label = "Male", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1-0.225, y = y_lab_top, label = "Female", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1+0.225, y = y_lab_top, label = "Male", size = y_lab_size, fontface = "bold") +
  theme(
    plot.margin = margin(t = 10, r = 40, b = 20, l = 10, unit = "pt"),
    plot.title = element_text(vjust=5, hjust=-0.35, size=16, face="bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )

p3 =
  ggplot() +  ggtitle("C") +
  ylab("Intra-specific similarity") + xlab("Aposematic coloration (LD1)") +
  coord_cartesian(xlim=c(-2.5, 6), ylim=c(0.7,1)) +
  geom_point(data=data_mod3_phylo_agg[unpalatability=="unknown"],
             aes(x=LD1, y=mean_sim_all, color=unpalatability, label=species), size=1.2, pch=16) + 
  geom_point(data=data_mod3_phylo_agg[!unpalatability=="unknown"],
             aes(x=LD1, y=mean_sim_all, fill=unpalatability), pch=21, size=2) +
  geom_ribbon(data=mod_3a_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_3a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_ribbon(data=mod_3b_phylo_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_3b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  geom_point(data=p3_grid_data, aes(x=LD1, y=mean_sim_all, fill=unpalatability, label=species), size=6, pch=21) +
  geom_text(data=p3_grid_data, aes(x=LD1, y=mean_sim_all, label=number), size=3) +
  plot_elements

y_lab_top = 3.65
p3_g =
  ggplot(p3_grid_data) + 
  theme_void() + coord_cartesian(clip = "off") +
  geom_image(aes(image = roi_path_h, x = grid_x-0.225, y = grid_y), size = 0.22, by="width") +
  geom_image(aes(image = roi_path_l, x = grid_x+0.225, y = grid_y), size = 0.22, by="width") +
  geom_textbox(aes(x = grid_x, y = grid_y-0.35, label = label), size=text_size, width=unit(4.6, "cm"),
               height = unit(0.5, "cm"), box.r = unit(7, "pt"), lineheight = 1.4, box.padding = unit(0.2, "lines"),
               fontface = "italic", halign = 0.5, valign=1) +
  geom_point(data=p3_grid_data, aes(x=grid_x-0.42, y=grid_y-0.35, fill=unpalatability), size=6.4, pch=21, stroke = 0.4) + 
  geom_text(data=p3_grid_data, aes(x=grid_x-0.42, y=grid_y-0.35, label=number), size=3) + 
  scale_fill_manual(values=colors_unpala, guide = "none", drop=T) +  
  annotate("text", x = 0-0.225, y = y_lab_top, label = "High apo.", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 0+0.225, y = y_lab_top, label = "Low apo.", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1-0.225, y = y_lab_top, label = "High apo.", size = y_lab_size, fontface = "bold") +
  annotate("text", x = 1+0.225, y = y_lab_top, label = "Low apo.", size = y_lab_size, fontface = "bold") +
  theme(
    plot.margin = margin(t = 10, r = 40, b = 35, l = 10, unit = "pt"),
    plot.title = element_text(vjust=5, hjust=-0.35, size=16, face="bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )


p_leg =
  ggplot(data_mod3_phylo_agg) + ggtitle("A") + 
  geom_point(data=data_mod3_phylo_agg[unpalatability=="unknown"],
             aes(x=LD1, y=mean_sim_all, color=unpalatability, label=species), size=1.2, pch=16) + 
  geom_point(data=data_mod3_phylo_agg[!unpalatability=="unknown"],
             aes(x=LD1, y=mean_sim_all, fill=unpalatability), pch=21, size=2) +
  
  geom_ribbon(data=mod_3a_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se, alpha=sex)) +
  geom_line(data=mod_3a_pred_glob, aes(x=LD1, y=fit, linetype = sex)) +
  scale_color_manual(
    "Chemical defense", values=colors_unpala, labels=labels_unpala, drop=F,
    guide=guide_legend(order=1, title.position="top", ncol=3, byrow=T)) +
  scale_fill_manual(
    "Chemical defense", values=colors_unpala, labels=labels_unpala, drop=F,
    guide=guide_legend(order=1, title.position="top", ncol=3)) +
  scale_linetype_manual(
    "Model fit (GAM)", values=c(1,2), labels=c(
      "All pattern aspects (Dim 1-177)", "Aposematic pattern aspects (Dim 1-7)"),
    guide=guide_legend(order=2, keywidth = unit(40, "pt"), title.position="top",
                       ncol=1)) + 
  scale_alpha_manual(
    "Model fit (GAM)", values=c(0.2, 0.2), labels=c(
      "All pattern aspects (Dim 1-177)", "Aposematic pattern aspects (Dim 1-7)"),
    guide=guide_legend(order=2, keywidth = unit(40, "pt"), title.position="top",
                       ncol=1)) + 
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal", 
    legend.box.just = "left",
    legend.margin=margin(l = 1.5, unit='cm'),
    axis.title.x = element_blank()
  )

p1_grid = cowplot::plot_grid(p1, p1_g, ncol=2, rel_widths = c(0.5,0.5))
p2_grid = cowplot::plot_grid(p2, p2_g, ncol=2, rel_widths = c(0.5,0.5))
p3_grid = cowplot::plot_grid(p3, p3_g, ncol=2, rel_widths = c(0.5,0.5))
p_leg_grid <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]

p_grid = cowplot::plot_grid(p1_grid, p2_grid, p3_grid, ncol=1, rel_heights = c(0.3,0.3,0.31))
p_grid_save = cowplot::plot_grid(p_grid, p_leg_grid, ncol=1, rel_heights = c(0.9, 0.1))

width=2500
height=4000
ggsave(paste0("figures/figure3.png"), p_grid_save,
       width=width, height=height, units = c("px"), bg="white")


# figure 4 - stats1  --------------------------------------------------

# ## prep files
# species_ordered = data_tree_sub1$tip.label
# 
# data_ld1_agg = data_ld1_agg[species %in% species_ordered]
# data_ld1_agg[, species := factor(species, levels=species_ordered)]
# data_ld1_agg = data_ld1_agg[order(species)]
# 
# data_ld1_agg_dv = data_ld1_agg_dv[species %in% species_ordered]
# data_ld1_agg_dv[, species := factor(species, levels=species_ordered)]
# data_ld1_agg_dv = data_ld1_agg_dv[order(species)]
# 
# data_ld1_agg_dv_sex = data_ld1_agg_dv_sex[species %in% species_ordered]
# data_ld1_agg_dv_sex[, species := factor(species, levels=species_ordered)]
# data_ld1_agg_dv_sex = data_ld1_agg_dv_sex[order(species)]

## this takes a very long time to run (multiple hours on a HPC)
traits_list <- list()

## LD1
traits_list[["LD1"]] <- setNames(
  data_ld1_agg$LD1,
  data_ld1_agg$species
)

## LD1 dorsal/ventral
traits_list[["LD1_d"]] <- setNames(
  data_ld1_agg_dv[class_dv == "dorsal"]$LD1,
  data_ld1_agg_dv[class_dv == "dorsal"]$species
)
traits_list[["LD1_v"]] <- setNames(
  data_ld1_agg_dv[class_dv == "ventral"]$LD1,
  data_ld1_agg_dv[class_dv == "ventral"]$species
)

## LD1 dorsal/ventral / sex
traits_list[["LD1_f_d"]] <- setNames(
  data_ld1_agg_dv_sex[sex == "female" & class_dv == "dorsal"]$LD1,
  data_ld1_agg_dv_sex[sex == "female" & class_dv == "dorsal"]$species
)
traits_list[["LD1_f_v"]] <- setNames(
  data_ld1_agg_dv_sex[sex == "female" & class_dv == "ventral"]$LD1,
  data_ld1_agg_dv_sex[sex == "female" & class_dv == "ventral"]$species
)
traits_list[["LD1_m_d"]] <- setNames(
  data_ld1_agg_dv_sex[sex == "male" & class_dv == "dorsal"]$LD1,
  data_ld1_agg_dv_sex[sex == "male" & class_dv == "dorsal"]$species
)
traits_list[["LD1_m_v"]] <- setNames(
  data_ld1_agg_dv_sex[sex == "male" & class_dv == "ventral"]$LD1,
  data_ld1_agg_dv_sex[sex == "male" & class_dv == "ventral"]$species
)

for(trait_name in names(traits_list)){
  trait = traits_list[[trait_name]]
  message(trait_name, " ", length(trait))
}

subfolder = "data/analyses_secondary/phylo_models"
for (context in names(traits_list)) {
  message("Starting HRM fit for ", context)
  hrm_fit_file <- paste0(subfolder, "/hrm_fit_", context, ".RData")
  
  if (!file.exists(hrm_fit_file)) {
    hrm_fit <- NULL  
    trait <- traits_list[[context]]
    tree <- drop.tip(data_tree_sub1, setdiff(data_tree_sub1$tip.label, names(trait)))
    set.seed(42)
    hrm_fit <- try(fitmultiBM(
      tree, trait, ncat = 2, model.hrm = "ER", levs = 100, root = "nuisance",
      rand_start=TRUE, opt.method="optim",
    ), silent = TRUE)
    if (inherits(hrm_fit, "try-error") || is.null(hrm_fit)) next
    save(hrm_fit, file = hrm_fit_file)
    message("Completed HRM fit for ", context, " - saved to drive")
  } else {
    message("Skipping HRM fit for ", context, " - file already exists")
  }
}

for (context in names(traits_list)) {
  hrm_fit_file <- paste0(subfolder,"/hrm_fit_", context, ".RData")
  if (file.exists(hrm_fit_file)) {
    message("Starting ANC fit for ", context)
    anc_fit_file <- paste0(subfolder,"/anc_fit_", context, ".RData")
    if (!file.exists(anc_fit_file)) {
      load(hrm_fit_file)
      anc_hrm <- ancr(hrm_fit, tips=TRUE)
      save(anc_hrm, file=anc_fit_file)
      message("Completed ANC fit for ", context, " - saved to drive")
    } else {
      message("Skipping ANC fit for ", context, " - file already exists")
    }
  }
}


# figure 4 - stats2  --------------------------------------------------

## collect hrm fits
files_hrm_fit = list.files(subfolder)[str_detect(list.files(subfolder), "hrm")]
hrm_models = list()
for (file in files_hrm_fit){
  load(paste0(subfolder, "/", file))
  context = str_remove_all(file, ".RData")
  hrm_models[[context]] = hrm_fit
}
hrm_list_sigsq = list()
for (hrm_fit_name in names(hrm_models)){
  hrm_fit = hrm_models[[hrm_fit_name]]
  hrm_list_sigsq[[hrm_fit_name]] = data.table(
    context = hrm_fit_name,
    state=c("a", "b"),
    rate=hrm_fit$states,
    sigsq = hrm_fit$sigsq
  )
}

data_hrm_fit = rbindlist(hrm_list_sigsq)
data_hrm_fit[, context_apo := ifelse(str_detect(context, "all"), "all", "apo")]
data_hrm_fit[, context_trait := substr(context, 9,100)]
data_hrm_fit[, rate := factor(rate)]
data_hrm_fit[, rate_label := ifelse(
  sigsq == max(sigsq), "higher", "lower"
), by = context_trait]

## collect anc fits 
files_anc_fit = list.files(subfolder)[str_detect(list.files(subfolder), "anc")]
anc_hrm_models = list()
for (file in files_anc_fit){
  load(paste0(subfolder, "/", file))
  context = str_remove_all(file, ".RData")
  anc_hrm_models[[context]] = anc_hrm
}
anc_list = list()
for (anc_fit_name in names(anc_hrm_models)){
  anc_fit = anc_hrm_models[[anc_fit_name]]
  anc_fit_results = anc_fit$ace$discrete
  species_names = rownames(anc_fit_results)[is.na(as.numeric(rownames(anc_fit_results)))]
  
  anc_discrete_long=data.table::melt(
    data.table(species=species_names,
               context=anc_fit_name,
               anc_fit$ace$discrete[1:length(species_names),]),
    measure.vars = colnames(anc_fit$ace$discrete),
    variable.name="rate", value.name = "prob")
  
  anc_list[[anc_fit_name]] = anc_discrete_long
}

data_anc_fit = rbindlist(anc_list)
data_anc_fit[, context_apo := ifelse(str_detect(context, "all"), "all", "apo")]
data_anc_fit[, context_trait := substr(context, 9,100)]

## merge
data_rates = merge(
  data_anc_fit[,c("species","rate","prob","context_apo","context_trait")],
  data_hrm_fit[,c("rate", "rate_label","sigsq","context_apo","context_trait")],
  by=c("rate","context_apo","context_trait")
)
data_rates = data_rates[rate_label=="higher"]
data_rates = merge(data_meta_agg, data_rates) 
data_rates[str_detect(context_trait, "_v$"), class_dv := "ventral"]
data_rates[str_detect(context_trait, "_d$"), class_dv := "dorsal"]
data_rates[str_detect(context_trait, "_m$|_m_"), sex := "male"]
data_rates[str_detect(context_trait, "_f$|_f_"), sex := "female"]
data_rates[, rate_est := prob * sigsq]
data_rates[, variable := ifelse(str_detect(context_trait, "LD1"), "LD1",
                                  ifelse(str_detect(context_trait, "dv_sim"), "dv_sim", "sex_sim"))]
data_rates[, context_trait := gsub("_s$", "", gsub("_s_", "_", context_trait))]
# fwrite(data_rates, "data/data_secondary/ld1_rates.csv")

## LD1 
data_rates_ld1 = merge(
  data_rates[context_trait %in% c("LD1")], 
  data_ld1_agg, by=c("species"))
data_rates_ld1_dv = merge(
  data_rates[context_trait %in% c("LD1_v", "LD1_d")], 
  data_ld1_agg_dv, by=c("species", "class_dv"))
data_rates_ld1_dv_sex = merge(
  data_rates[context_trait %in% c("LD1_f_d", "LD1_f_v", "LD1_m_d", "LD1_m_v")], 
  data_ld1_agg_dv_sex, by=c("species", "class_dv", "sex"))

## ld1 
p1 = 
  ggplot(data_rates_ld1) +
  geom_smooth(aes(x=LD1, y=prob), method = "gam", 
              formula = y ~ s(x, k = 10), method.args = list(family = betar(link = "logit"))) +
  geom_point(aes(x=LD1, y=prob, fill=unpalatability), pch=21, size=2) +
  scale_fill_manual("Chemical defense", values=colors_unpala, labels=labels_unpala, drop=T) +
  theme(legend.position = "bottom")

## ld1 dv
p2 = 
  ggplot(data_rates_ld1_dv) +
  facet_wrap(~class_dv, ncol=1, scales="free") +
  geom_smooth(aes(x=LD1, y=prob), method = "gam", 
              formula = y ~ s(x, k = 10), method.args = list(family = betar(link = "logit"))) +
  geom_point(aes(x=LD1, y=prob, fill=unpalatability), pch=21, size=2) +
  scale_fill_manual("Chemical defense", values=colors_unpala, labels=labels_unpala, drop=T) +
  theme(legend.position = "bottom")

## LD1 dv / sex
p3=
  ggplot(data_rates_ld1_dv_sex) +
  facet_wrap(class_dv~sex, ncol=2, scales="free") +
  geom_smooth(aes(x=LD1, y=prob), method = "gam", 
              formula = y ~ s(x, k = 10), method.args = list(family = betar(link = "logit"))) +
  geom_point(aes(x=LD1, y=prob, fill=unpalatability), pch=21, size=2) +
  scale_fill_manual("Chemical defense", values=colors_unpala, labels=labels_unpala, drop=T) +
  theme(legend.position = "bottom")

pg1 = cowplot::plot_grid(p1,p2,p3, ncol=3)
# pg1

data_rates[,table(context_trait)]

tree = copy(data_tree_sub3) 

data_mod_4 = data_rates_ld1_dv_sex[variable == "LD1" & context_apo=="apo"]
data_mod_4[, species := factor(species, levels=tree$tip.label)]
data_mod_4[, sex := factor(sex, levels=c("female", "male"))]
data_mod_4[, class_dv := factor(class_dv, levels=c("dorsal", "ventral"))]

# rate regression (GAM 4)
phylo_penalty <- mrf_penalty(tree)
knots = 10
# mod_4 = gam(prob ~
#                s(LD1, k = knots) +
#                s(LD1, sex, k=knots, bs = "sz") +
#                s(LD1, class_dv, k=knots, bs = "sz") +
#                s(LD1, sex, class_dv, k=knots, bs = "sz") +
#                s(species, bs = "mrf", xt = list(penalty = phylo_penalty), k = 200),
#              family = betar(link = "logit"),
#              method="REML",
#              data=data_mod_4
# )
# save(mod_4, file="data/analyses_secondary/gam_regressions/mod_4.RData")
load(file="data/analyses_secondary/gam_regressions/mod_4.RData")
summary(mod_4)
# capture.output(anova(mod_4), file = "tables/mod_4.txt")
# par(mfrow=c(2,2))
# gam.check(mod_4)

# figure 4/S8 - prep ---------------------------------------------------------

## format heatmap
tree <- data_tree_sub4
data_traits_all <- dcast(
  data_ld1_agg_dv_sex[species %in% tree$tip.label],
  species ~ class_dv + sex, value.var = "LD1")
traits_all = data.frame(data_traits_all[,-"species",with=F])
rownames(traits_all) = data_traits_all$species

# keep only species present in both
tree <- drop.tip(tree, setdiff(tree$tip.label, data_traits_all$species))

## traits
trait_left <- setNames(
  data_ld1_agg[species %in% tree$tip.label]$LD1,
  data_ld1_agg[species %in% tree$tip.label]$species
  )
trait_left <- trait_left[tree$tip.label]  # reorder to tip order
trait_right <- setNames(
  droplevels(data_meta_agg[species %in% tree$tip.label]$unpalatability),
  data_meta_agg[species %in% tree$tip.label]$species
)
trait_right <- trait_right[tree$tip.label]  # reorder to tip order

## left fit
anc_fit_left <- fastAnc(tree, trait_left, vars = FALSE, CI = FALSE)
tips_left  <- data.frame(node = match(names(trait_left), tree$tip.label),
                       trait_left = as.numeric(trait_left))
nodes_left <- data.frame(node = as.numeric(names(anc_fit_left)),
                         trait_left = as.numeric(anc_fit_left))
tree_map_left <- bind_rows(tips_left, nodes_left)

## right fit 
anc_fit_right <- fastAnc(tree, trait_right, vars = FALSE, CI = FALSE)
tips_right  <- data.frame(node = match(names(trait_right), tree$tip.label),
                          trait_right = as.numeric(trait_right))
nodes_right <- data.frame(node = as.numeric(names(anc_fit_right)),
                          trait_right = as.numeric(anc_fit_right))
tree_map_right <- bind_rows(tips_right, nodes_right)

tree_phylo <- as.phylo(tree)  # strip simmap info
tree_phylo <- full_join(tree_phylo, tree_map_left, by = 'node')
tree_phylo <- full_join(tree_phylo, tree_map_right, by = 'node')
  
data_rates_ld1 = data.table(
  LD1_f_d = data_rates[context_trait=="LD1_f_d", prob],
  LD1_m_d = data_rates[context_trait=="LD1_m_d", prob],
  LD1_f_v = data_rates[context_trait=="LD1_f_v", prob],
  LD1_m_v = data_rates[context_trait=="LD1_m_v", prob],
  species = data_rates[context_trait=="LD1_m_d", species]
)

data_rates_ld1 = data_rates_ld1[species %in% tree$tip.label]
rates_left = data.frame(data_rates_ld1[,1:4])
rownames(rates_left) = data_rates_ld1$species
data_rates_ld1_long = data.table::melt(data_rates_ld1)
data_rates_ld1_long[, `:=`(x = 1, y = 1)]

## family labels
species_by_family <- data_meta_agg[species %in% tree$tip.label, .(species_list = list(species)), by = "subfamily"]
species_by_family <- setNames(species_by_family$species_list, species_by_family$subfamily)
tree_phylo = groupOTU(tree_phylo, species_by_family, group_name = "subfamily")

## tribe labels
species_by_tribe <- data_meta_agg[species %in% tree$tip.label, .(species_list = list(species)), by = "tribe"]
species_by_tribe <- setNames(species_by_tribe$species_list, species_by_tribe$tribe)
tree_phylo = groupOTU(tree_phylo, species_by_tribe, group_name = "tribe")

## mod tree data
ladderize = T
data_tree = data.table(ggtree(tree_phylo , ladderize = ladderize)$data)
data_tree = data_tree[order(y)]
data_tree[, tribe_alternate := rleid(tribe)]
data_tree[, tribe_factor := factor(
  tribe_alternate, labels = rep(c("A", "B"), length.out = max(tribe_alternate)))]
tribe_values = data.frame(data_tree$tribe_factor)
rownames(tribe_values) = data_tree$species
data_tree[, subfamily := factor(subfamily, levels=c(
  "0", "satyrinae", "charaxinae", "nymphalinae", 
  "cyrestinae", "biblidinae", "apaturinae", "pseudergolinae", "heliconiinae",
  "limenitidinae", "libytheinae", "danainae"))]
data_tree[, unpalatability := "unknown"]
data_tree[, x_tribe := ifelse(tribe_factor=="A",1,2)]

# figure 4 - plot --------------------------------------------------------------

mod_4_pred_part = pred_gam_newd(
  mod=mod_4, varx = "LD1",  var_cat=c("sex", "class_dv"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c("s(species)"))
mod_4_pred_glob = pred_gam_newd(
  mod=mod_4, varx = "LD1",  var_cat=c("sex", "class_dv"), 
  extra_x = c(0, 0), extra_y = c(0, 0), se=T,
  n_gridlines = 50, too_far_val = too_far_val,
  type="response", newdata_check = F,
  exclude=c("s(species)", "s(LD1,sex)", "s(LD1,class_dv)", "s(LD1,class_dv,sex)"))
data_mod_4[, group := paste0(class_dv, "_", sex)]
mod_4_pred_part[, group := paste0(class_dv, "_", sex)]
mod_4_pred_glob[, group := paste0(class_dv, "_", sex)]
group_labels = c( 
  "dorsal_female"="(1) Female dorsum", 
  "dorsal_male"="(2) Male dorsum", 
  "ventral_female"="(3) Female ventrum", 
  "ventral_male"="(4) Male ventrum"
  )

p1_reg =
  ggplot(data_mod_4[species_label=="label"]) + 
  coord_cartesian(ylim=c(0,1), xlim=c(-2.5, 6)) + 
  xlab("Aposematic coloration (LD1)") + ylab("Probability of higher evolutionary rate") + 
  facet_wrap(~group, ncol=1, labeller = as_labeller(group_labels)) +
  geom_point(data=data_mod_4[species_label=="no-label"],
             aes(x=LD1, y=prob, color=unpalatability), pch=16, size=1, color="gray") +
  geom_point(aes(x=LD1, y=prob, fill=unpalatability), pch=21, size=2) +
  geom_ribbon(data=mod_4_pred_part, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_4_pred_part, aes(x=LD1, y=fit), linetype=2) +
  geom_ribbon(data=mod_4_pred_glob, aes(x=LD1, ymin=fit-se, ymax=fit+se), alpha=0.3) +
  geom_line(data=mod_4_pred_glob, aes(x=LD1, y=fit)) +
  scale_fill_manual("Chemical defense", values=colors_unpala, labels=labels_unpala, drop=T) +
  scale_y_continuous(breaks=c(0,0.5,1)) + 
  theme(legend.position = "none",
        axis.title.y = element_text(margin=margin(r=10, unit="pt")),
        plot.margin = margin(r=10, t=10, b=15, unit="pt"))

p1 = ggtree(tree_phylo , ladderize = ladderize) +
  geom_tree(aes(color = trait_left),  size = 1, continuous = "color") + 
  scale_color_viridis_c(option = "plasma", guide="none") + 
  geom_tippoint(aes(x=-10, fill=subfamily), size=1.5, pch=22, stroke=NA) +
  scale_fill_manual(values=fam_cols, "Subfamily", guide=guide_legend(
    order=1, override.aes = list(linewidth=1), title.position = "top", ncol=1),
    labels=fam_labs) +
  scale_x_continuous(expand=c(0.01,0), limits=c(-10, max(data_tree$x))) +
  theme_void() + theme(
    legend.position = "none",
    plot.margin = margin(t=28,b=15,l=15,r=0, unit="pt")
    ) 

for (fam in unique(data_tree$subfamily)) {
  species <- data_tree[subfamily == fam & !is.na(label)]$label
  if (length(species) > 1) {
    mrca_node <- getMRCA(tree_phylo@phylo, species)
    p1 <- p1 + 
      geom_hilight(node = mrca_node, fill = fam_cols[[fam]], alpha = 0.2, 
                            to.bottom = T, extend = -200) +
      geom_hilight(node = mrca_node, fill = fam_cols[[fam]], alpha = 0.2, 
                   to.bottom = T, extend = 0)
  }
}

p2 <- ggtree(tree_phylo, color = NA, size = 0, ladderize = ladderize) + 
  theme_void() +  new_scale_fill() + coord_cartesian(clip="off")
p2 <- gheatmap(p2, rates_left, width = 10, offset = 0,
               colnames = T, colnames_position = "top", 
               font.size = 5, family="bold",
               colnames_angle=0, hjust=0.5, colnames_offset_y = 5,
               custom_column_labels=c("1", "2", "3", "4")) +
  scale_fill_viridis_c(guide = "none") +   
  scale_x_continuous(expand=c(0,0)) +
  theme(plot.margin = margin(t=15.5,b=15,l=-12,r=3, unit="pt"))

cmap_right <- colorRampPalette(colors_unpala[2:4])(100)
p3 <-
  ggtree(tree_phylo, aes(color = trait_right),
         size = 1, continuous = "color", ladderize = ladderize) +
  scale_x_continuous(expand=c(0.01,0), limits=c(-10, max(data_tree$x))) +
  scale_color_gradientn(colors = cmap_right, guide = "none") +
  theme_void() + 
  theme(
    legend.position = "none", 
    plot.margin = margin(t=28,b=15,l=0,r=15, unit="pt") 
    )

for (fam in unique(data_tree$subfamily)) {
  species <- data_tree[subfamily == fam & !is.na(label)]$label
  if (length(species) > 1) {
    mrca_node <- getMRCA(tree_phylo@phylo, species)
    p3 <- p3 + 
      geom_hilight(node = mrca_node, fill = fam_cols[[fam]], alpha = 0.2, 
                   to.bottom = T, extend = -200) +
      geom_hilight(node = mrca_node, fill = fam_cols[[fam]], alpha = 0.2, 
                   to.bottom = T, extend = 0)
  }
}

tribes_sub = c(
  ## danainae
  "danaini", "ithomiini", 
  ## heliconiinae
  "heliconiini", "argynnini", "acraeini",
  ## nymphalinae
  "melitaeini",
  "charaxini",
  "biblidini")
trib = "acraeini"
for (trib in tribes_sub) {
  species <- data_tree[tribe == trib & !is.na(label)]$label
  # if (length(species) > 10) {
  mrca_node <- getMRCA(tree_phylo@phylo, species)
  node_info <- data_tree[node == mrca_node]
  node_info[, tribe := str_to_title(tribe)]
  if (nrow(node_info) == 1) {
    p3 <- p3 + 
      geom_text(data = node_info,
                aes(x = x, y = y, label=tribe),
                size = 4, color = "black")
  }
  # }
}

p_leg1 =
ggplot(data_tree) +
  geom_point(aes(x=1,y=1,color=trait_left)) +
  scale_color_viridis_c(
    "Aposematic coloration (LD1)",
    option="inferno", 
    oob = scales::oob_squish,
    limits = c(-2.5,6.5),
    guide=guide_colorbar(
      order=1, keyheight=unit(0.5, "cm"), keywidth=unit(6,"cm"), 
      title.position = "top", direction = "horizontal")) +
  
  new_scale_color() + 
  geom_point(data=data_rates_ld1_long,
             aes(x=1, y=1, color=value)) +
  scale_color_viridis_c(
    "Prob. higher rate", 
    breaks=c(0,0.5,1), limits=c(0,1),
    guide=guide_colorbar(
      order=2, keyheight=unit(0.5, "cm"), keywidth=unit(4,"cm"), 
      title.position = "top", direction = "horizontal")) +
  
  new_scale_color() + 
  geom_point(aes(x=1,y=1,color=trait_right)) +
  scale_color_gradientn(
    "Chemical defense",
    colors = cmap_right, breaks=c(1,2,3), 
    labels = labels_unpala[2:4],
    guide=guide_colorbar(
      order=3, keyheight=unit(0.5, "cm"), keywidth=unit(4,"cm"), 
      title.position = "top", direction = "horizontal")) +
  
  new_scale_color() + 
  geom_point(aes(x=1,y=1,color=unpalatability)) +
  scale_color_manual(
    "",
    values=colors_unpala[1],
    labels=labels_unpala[1],
    guide=guide_legend(
      order=4, 
      theme = theme(legend.margin = margin(l=-10,unit="pt"),
                    legend.text = element_text(size=10)),
      label.position = "bottom",
      title.position = "top")) +
  
  theme(legend.position = "bottom",
        legend.margin = margin(r=20, l=20, unit="pt"))



p_leg2 =
  ggplot(data_tree[!subfamily==0]) +
  geom_ribbon(aes(x=1,y=1, ymin=1, ymax=1, fill=subfamily)) +
  
  scale_fill_manual(values=fam_cols, "Subfamily", guide=guide_legend(
      order=1, override.aes = list(linewidth=1), keywidth = unit(5, "pt"),
      title.position = "top", nrow=3, byrow=T), 
      labels=fam_labs) +
  geom_line(
    data=mod_4_pred_part,
    aes(x=LD1, y=fit, linetype=sex, group=sex)) + 
  geom_ribbon(
      data=mod_4_pred_part,
      aes(x=LD1, ymin=-1, ymax=1, linetype=sex, group=sex), alpha=0.2) + 
  scale_linetype_manual(
    "Model fit (GAM)",
    values=c(1,2), labels=c("global","partial"),
    guide=guide_legend(
      order=2, keywidth = unit(40, "pt"),
      title.position = "top", nrow=2)
  ) + 
  
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(r=20, l=20, unit="pt"))
  
pg1 = cowplot::plot_grid(p1,p2,mirror_plot(p3), ncol=3, rel_widths = c(0.4,0.2,0.4))
pg2 = cowplot::plot_grid(pg1, p1_reg, ncol=2, rel_widths = c(0.6,0.4),
                         labels=c('A', 'B'), label_size = 16)

leg1 <- get_plot_component(p_leg1, 'guide-box', return_all = TRUE)[[3]]
leg2 <- get_plot_component(p_leg2, 'guide-box', return_all = TRUE)[[3]]

pg3 = cowplot::plot_grid(leg1,leg2, ncol=1, rel_heights = c(0.4,0.6))
pg_final = cowplot::plot_grid(pg2,pg3, ncol=1, rel_heights = c(0.825,0.175))


height=4000
width=2500
ggsave(paste0("figures/figure4.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")

# figure S2 ---------------------------------------------------------------

## KDE on aggregated coords
data_plot = merge(data_meta_agg, pca_coords_butterflies_agg)

## contours
data_KDE = data_plot[species_label=="label", {
  contour_dt <- kde_2d(Dim1,Dim2,levels = c(0.25, 0.5, 0.75),
                       n_grid = 100, bw_adjust = 1)
}, by = .(unpalatability)]

## grid
data_KDE_grid <- data_plot[species_label == "label", {
  kde_grid_2d(Dim1, Dim2, n_grid = 100, bw_adjust = 1, expand=0.2)
}, by = .(unpalatability)]
probs <- c(0.05, 0.1, 0.25, 0.5, 0.75)
thr_long <- data_KDE_grid[
  , {
    o   <- order(density)
    cum <- cumsum(pmass[o])
    tot <- sum(pmass)
    cut <- sapply(probs, function(p) density[o][ which.max(cum >= p*tot) ])
    data.table(prob = paste0("p",substr(as.character(probs), 3, 20)), cutoff = as.numeric(cut))
  },
  by = .(unpalatability)
]
thr_wide <- dcast(thr_long, unpalatability ~ prob, value.var = "cutoff")
data_KDE_grid <- thr_wide[data_KDE_grid, on=.(unpalatability)]


p1 =
  ggplot(data_quanti_supp[feature_string %in% data_quanti_supp[include==T]$feature_string]) +
  coord_fixed(ylim=c(-25,17),xlim=c(-20,20)) + 
  
  xlab(paste("Dim1 - Variance explained: ", round(pca_eigenvalues[1,]$var_explained,2), "%")) +
  ylab(paste("Dim2 - Variance explained: ", round(pca_eigenvalues[2,]$var_explained,2), "%")) +
  scale_x_continuous(breaks=seq(-100,100,10), expand=c(0,0)) + 
  scale_y_continuous(breaks=seq(-100,100,10), expand=c(0,0)) + 
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  
  ## points without labels
  geom_point(
    data=data_plot[species_label=="no-label"],
    aes(x=Dim1, y=Dim2, color=unpalatability), size=1, inherit.aes = FALSE) +
  scale_color_manual("Chemical defense:", values=colors_unpala, labels=labels_unpala, drop=FALSE,
                     guide = guide_legend(keywidth = unit(1.75 , "cm"), order=1,
                                          override.aes = list(size=c(2,3,3,3)))) +
  
  ## kde grid
  geom_raster(
    data=data_KDE_grid[unpalatability == "low" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "medium" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "high" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  scale_alpha_continuous(
    range = c(0, 1),
    trans = "sqrt",
  ) +
  
  ## kde contours
  geom_path(
    data=data_KDE,
    aes(x=x, y=y, color=unpalatability, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dashed")) +
  
  ## poins with labels
  geom_point(data=data_plot[species_label=="label"],
             aes(x=Dim1, y=Dim2, fill=unpalatability), pch=21, size=2) +
  scale_fill_manual(values=colors_unpala) + 
  
  ## LDA label
  geom_segment(data=data_quanti_supp[feature_string=="LD1"],
               aes(x = 0, y = 0, xend = Dim1_s, 
                   yend = Dim2_s), size=1, 
               arrow = arrow(angle = 30, length = unit(2, "mm"),
                             ends = "last", type = "closed")) +
  geom_label(data=data_quanti_supp[feature_string=="LD1"],
             aes(x = Dim1_s, y = Dim2_s, 
                 label=feature_string), size=5, show.legend = F,
             fill=colors_unpala["high"]) +
  
  ## features
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  # geom_segment(aes(x = 0, y = 0, xend = Dim1_s, yend = Dim2_s)) +  # Arrows from (0,0)
  geom_label_repel(aes(x = Dim1_s, y = Dim2_s, label = feature_formatted), 
                   min.segment.length = 0, size = 4, 
                   force=1, direction="both",
                   # force_pull=1
  ) +  
  geom_point(aes(x = Dim1_s, y = Dim2_s), size=2, pch=16) +  
  
  ## theme settings
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14)
  )

p_leg =
  ggplot(data_plot) + 
  
  geom_point(
    data=data_plot,
    aes(x=Dim1, y=Dim2, color=unpalatability), size=1) +
  scale_color_manual(
    "Chemical defense", values=colors_unpala, labels=labels_unpala, 
    guide = guide_legend(
      order=1, keywidth = unit(40 , "pt"),
      override.aes = list(size=c(3,3,3,1)), title.position="top")) +
  ## KD lines
  geom_path(
    data=data_KDE, 
    aes(x=x, y=y, color=unpalatability, linetype=factor(level), 
        group=interaction(level_idx,unpalatability, group)), linewidth=1.25) +
  
  scale_linetype_manual(
    "Spatial density quantiles (KDE Dim1&2)", values=c("solid","longdash","dashed"),
    labels=c("25%","50%","75%"), 
    guide = guide_legend(order=2, keywidth = unit(65 , "pt"), title.position="top")) +
  
  theme(
    legend.margin = margin(l = 50, unit = "pt"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left")

leg <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]
pg_final = cowplot::plot_grid(p1, leg, nrow=2, rel_heights = c(0.85,0.15))


width = 2500
height = 3000
ggsave(paste0("figures/figureS2a.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")

# figure S3 - prep ----------------------------------------------------------------

data_plot = merge(data_meta_agg, pca_coords_butterflies_agg_dv_sex)

## contours
data_KDE = data_plot[species_label=="label", {
  contour_dt <- kde_2d(Dim1,Dim2,levels = c(0.25, 0.5, 0.75),
                       n_grid = 100, bw_adjust = 1)
}, by = .(unpalatability, class_dv, sex)]

## grid
data_KDE_grid <- data_plot[species_label == "label", {
  kde_grid_2d(Dim1, Dim2, n_grid = 100, bw_adjust = 1, expand=0.2)
}, by = .(unpalatability, class_dv, sex)]
probs <- c(0.05, 0.1, 0.25, 0.5, 0.75)
thr_long <- data_KDE_grid[
  , {
    o   <- order(density)
    cum <- cumsum(pmass[o])
    tot <- sum(pmass)
    cut <- sapply(probs, function(p) density[o][ which.max(cum >= p*tot) ])
    data.table(prob = paste0("p",substr(as.character(probs), 3, 20)), cutoff = as.numeric(cut))
  },
  by = .(unpalatability, class_dv, sex)
]
thr_wide <- dcast(thr_long, unpalatability + class_dv+sex ~ prob, value.var = "cutoff")
data_KDE_grid <- thr_wide[data_KDE_grid, on=.(unpalatability, class_dv, sex)]

data_rois_dv_sex_subset = data_rois_dv_sex[species %in% data_roi_subset$species]
data_rois_dv_sex_subset = merge(
  data_plot[, c("species","class_dv", "sex","Dim1", "Dim2")], 
  data_rois_dv_sex_subset, by=c("species", "class_dv", "sex"))

## insets
data_plot2 = merge(data_meta_agg, data_ld1_agg_dv_sex)
inset_settings = list(
  coord_cartesian(ylim=c(0,0.55), xlim=c(-2.5, 6)),
  theme_cowplot(), 
  geom_density(aes(x = LD1, fill=unpalatability), alpha=0.5), 
  scale_fill_manual(values=colors_unpala, labels=labels_unpala, guide="none"),
  scale_y_continuous(expand=c(0,0)),
  theme(
    axis.text = element_text(size=10),
    axis.title.x = element_text(size=10, margin = margin(b=-2,t=2,unit="pt")),
    axis.title.y = element_text(size=10, margin = margin(l=-2,r=2,unit="pt")),
    plot.background = element_rect(colour = "black", fill = "white", size=1),
    plot.margin = margin(t=10,b=6,l=6,r=6, unit="pt")
  ) 
)
inset_d_f <-
  ggplot(data_plot2[class_dv=="dorsal" & sex=="female" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings
inset_d_m <-
  ggplot(data_plot2[class_dv=="dorsal" & sex=="male" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings
inset_v_f <-
  ggplot(data_plot2[class_dv=="ventral" & sex=="female" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings
inset_v_m <-
  ggplot(data_plot2[class_dv=="ventral" & sex=="male" & !unpalatability=="unknown" ]) + 
  labs(x="Aposematic color score (LD1)",y="Density") + 
  inset_settings


# figure S3 - plot ----------------------------------------------------------------

## figure
p1 =
  ggplot() + 
  facet_wrap(sex~class_dv, labeller = as_labeller(c(
    "male"="Male", "female"="Female", "dorsal"="dorsum","ventral"="ventrum"
  ), multi_line=F)) + 
  xlab(paste("Dim1 - Variance explained: ", round(pca_eigenvalues[1,]$var_explained,2), "%")) +
  ylab(paste("Dim2 - Variance explained: ", round(pca_eigenvalues[2,]$var_explained,2), "%")) +
  scale_x_continuous(breaks=seq(-100,100,10)) + 
  scale_y_continuous(breaks=seq(-100,100,10)) + 
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  
  ## points without labels
  geom_point(
    data=data_plot[!species_label=="label"],
    aes(x=Dim1, y=Dim2, color=unpalatability), size=1, inherit.aes = FALSE) +
  scale_color_manual("Chemical defense:", values=colors_unpala, labels=labels_unpala, drop=FALSE,
                     guide = guide_legend(keywidth = unit(1.75 , "cm"), order=1,
                                          override.aes = list(size=c(2,3,3,3)))) +
  
  ## kde grid
  geom_raster(
    data=data_KDE_grid[unpalatability == "low" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "medium" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE_grid[unpalatability == "high" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  scale_alpha_continuous(
    range = c(0, 1),
    trans = "sqrt",
  ) +  
  
  ## kde contours
  geom_path(
    data=data_KDE,
    aes(x=x, y=y, color=unpalatability, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dashed")) +
  
  ## points with labels
  geom_point(
    data=data_plot[species_label=="label"],
    aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), size=2, pch=21) +
  scale_fill_manual("Chemical defense:", values=colors_unpala, labels=labels_unpala, guide="none") +
  
  ## pictograms
  geom_image(
    data=data_rois_dv_sex_subset,
    aes(image=roi_path, x=Dim1, y=Dim2), size=0.1, by="height",
    image_fun=function(img) {img = magick::image_background(img, "none")}) +
  
  ## points for examples
  geom_point(
    data=data_rois_dv_sex_subset,
    aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), size=3, pch=21) +
  
  ## loading arrow
  new_scale_fill() +
  geom_segment(data=pca_quanti_supp[feature_string=="LD1"],
               aes(x = 0, y = 0, xend =  Dim1*loadings_factor,
                   yend =  Dim2*loadings_factor), size=1,
               arrow = arrow(angle = 30, length = unit(2, "mm"),
                             ends = "last", type = "closed")) +
  geom_label(data=pca_quanti_supp[feature_string=="LD1"],
             aes(x =  Dim1*loadings_factor, y =  Dim2*loadings_factor,
                 label=feature_string, fill = feature_string), size=5, show.legend = FALSE) +
  scale_fill_manual(values=unname(colors_unpala["high"])) +
  
  ## model results
  geom_label(data=data.table(
    x=-26, 
    y=-29, 
    class_dv="dorsal", 
    sex="female",
    label=paste("Pillai effect size:", 0.181, "\nP-value: 0.001")),
    aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
    label.size=NA, alpha=0.85) +
  geom_label(data=data.table(
    x=-26, 
    y=-29, 
    class_dv="dorsal", 
    sex="male",
    label=paste("Pillai effect size:", 0.143, "\nP-value: 0.001")),
    aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
    label.size=NA, alpha=0.85) +
  geom_label(data=data.table(
    x=-26, 
    y=-29, 
    class_dv="ventral", 
    sex="female",
    label=paste("Pillai effect size:", 0.234, "\nP-value: 0.001")),
    aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
    label.size=NA, alpha=0.85) +
  geom_label(data=data.table(
    x=-26, 
    y=-29, 
    class_dv="ventral", 
    sex="male",
    label=paste("Pillai effect size:", 0.192, "\nP-value: 0.001")),
    aes(x=x,y=y,label=label), hjust=0, vjust=0, size = 4, fontface = "bold",
    label.size=NA, alpha=0.85) +
  
  ## theme settings
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    aspect.ratio = 1,
  )


inset_height = 0.12
inset_width = 0.2
p2 = ggdraw() +
  draw_plot(p1) +
  draw_plot(inset_d_f, x = 0.31, y = 0.53, width = inset_width, height = inset_height) +
  draw_plot(inset_v_f, x = 0.775, y = 0.53, width = inset_width, height = inset_height) +
  draw_plot(inset_d_m, x = 0.31, y = 0.05, width = inset_width, height = inset_height) +
  draw_plot(inset_v_m, x = 0.775, y = 0.05, width = inset_width, height = inset_height)



leg <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]
pg_final = cowplot::plot_grid(p2, leg, nrow=2, rel_heights = c(0.9,0.1))


width = 4000
height = 4500
ggsave(paste0("figures/figureS3.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")


# figure S4 ---------------------------------------------------------------

data_lda_results_long = data.table::melt(
  data_lda_results, id.vars = c("pcs","dataset", "class_dv"))
data_lda_results_long[, variable_f := str_to_title(variable)]
data_lda_results_long[, dataset := factor(dataset, levels=c("train", "test1_id", "test2_ood", "test3_moths"))]

p1 =
  ggplot(data_lda_results_long[pcs<=50]) +
  labs(
    x="N Dimensions",
    y="Value"
  ) + 
  facet_wrap(~variable_f, ncol=1) +
  geom_line(aes(x=pcs, y=value, color=dataset, linetype = class_dv, group=interaction(dataset,class_dv)), linewidth=1) +
  geom_vline(xintercept = 7, color="red", linetype=1, linewidth=1) + 
  scale_x_continuous(breaks=seq(0, 50, 5), limits=c(2,50), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0.2,1), expand=c(0,0)) +
  scale_linetype_manual(
    values=c("solid", "dashed"),
    labels=c("Dorsum", "Ventrum"),
    guide = guide_legend(keywidth = unit(3.5, "cm"), order=2)) + 
  scale_color_manual(
    values = c("train"="darkorchid1", "test1_id"="green", "test2_ood"="darkorange", "test3_moths"="blue"),
    labels = c("train"="Train", "test1_id"="Test (ID)", "test2_ood"="Test (OOD)", "test3_moths"="Moths"),
    guide = guide_legend(keywidth = unit(1.25, "cm"), order=1)) + 
  theme(legend.position="bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_blank(),
        plot.margin = margin(5, 15, 5, 5)
  )

width = 2500
height = 3000
ggsave(paste0("figures/figureS4.png"), p1,
       width=width, height=height, units = c("px"), bg="white")

# figure S5 - prep ---------------------------------------------------------------

data_plot1 = merge(
  data_lda_summ3[class_dv=="dorsal" & unpalatability %in% c("low", "high")], 
  pca_coords_butterflies_agg_dv[class_dv=="dorsal", c("species", "Dim1", "Dim2")], 
  by="species")

data_plot2 = merge(data_meta_moths_agg, data_lda_moths_summ1)
data_plot2 = merge(data_plot2, pca_coords_moths_agg[,c("species", "Dim1", "Dim2")], by="species")

## KDE butterflies
data_KDE1 <- data_plot1[, {
  contour_dt <- kde_2d(Dim1,Dim2,levels = c(0.25, 0.5, 0.75),  
                       n_grid = 100, bw_adjust = 1 )
}, by = .(unpalatability)]
data_KDE1_grid <- data_plot1[, {
  kde_grid_2d(Dim1, Dim2, n_grid = 100, bw_adjust = 1, expand=0.2)
}, by = .(unpalatability)]
probs <- c(0.05, 0.1, 0.25, 0.5, 0.75)
thr_long <- data_KDE1_grid[
  , {
    o   <- order(density)
    cum <- cumsum(pmass[o])
    tot <- sum(pmass)
    cut <- sapply(probs, function(p) density[o][ which.max(cum >= p*tot) ])
    data.table(prob = paste0("p",substr(as.character(probs), 3, 20)), cutoff = as.numeric(cut))
  },
  by = .(unpalatability)
]
thr_wide <- dcast(thr_long, unpalatability ~ prob, value.var = "cutoff")
data_KDE1_grid <- thr_wide[data_KDE1_grid, on=.(unpalatability)]

## KDE moths
data_KDE2 <-data_plot2[, {
  contour_dt <- kde_2d(Dim1,Dim2,levels = c(0.25, 0.5, 0.75),  
                       n_grid = 100, bw_adjust = 1)
}, by = .(unpalatability)]
data_KDE2_grid <- data_plot2[, {
  kde_grid_2d(Dim1, Dim2, n_grid = 200, bw_adjust = 1, expand=0.1)
}, by = .(unpalatability)]
probs <- c(0.05, 0.1, 0.25, 0.5, 0.75)
thr_long <- data_KDE2_grid[
  , {
    o   <- order(density)
    cum <- cumsum(pmass[o])
    tot <- sum(pmass)
    cut <- sapply(probs, function(p) density[o][ which.max(cum >= p*tot) ])
    data.table(prob = paste0("p",substr(as.character(probs), 3, 20)), cutoff = as.numeric(cut))
  },
  by = .(unpalatability)
]
thr_wide <- dcast(thr_long, unpalatability ~ prob, value.var = "cutoff")
data_KDE2_grid <- thr_wide[data_KDE2_grid, on=.(unpalatability)]

mask_coords1 = data_plot1[pred_success==F, ]
mask_coords2 = data_plot1[
  species %in% c(
    # high
    "parantica_aspasia", "argynnina_cyrila", "danaus_plexippus", 
    "poladryas_arachne", "altinote_stratonice",
    # low
    "neominois_ridingsii", "calisto_obscura","punapedaliodes_flavopunctata",
    "baeotus_beotus", "pyronia_cecilia"),]
data_roi_butterflies = rbind(mask_coords1, mask_coords2)
data_roi_butterflies[, roi_path := file.path(
  paste0("data_raw//segmentation_masks_clean//centroids_dv//dorsal//", species, ".png"))]

data_plot2[species %like% "lyge"]

mask_coords1 = data_plot2[pred_success==F, ]
mask_coords2 = data_plot2[
  species %in% c(
    # high
    "zygaena_lonicerae","rhyparia_purpurata","elophila_icciusalis", "pseudopanthera_macularia",
    "grammia_nevadensis",
    # low
    "evergestis_forficalis", "hypena_rostralis", "autographa_californica", "macrochilo_cribrumalis",
    "lygephila_pastinum"),]
data_roi_moths = rbind(mask_coords1, mask_coords2)
data_roi_moths[, roi_path := file.path(
  paste0("data_raw//segmentation_masks_moths//centroids", "//" , species, ".png"))]

## arrow length
loadings_factor = 25

# figure S5 - plot ---------------------------------------------------------------

p1 =
  ggplot(data=data_plot1) +
  coord_fixed() + 
  xlab(paste("Dim1 - Variance explained: ", round(pca_eigenvalues[1,]$var_explained,2), "%")) +
  ylab(paste("Dim2 - Variance explained: ", round(pca_eigenvalues[2,]$var_explained,2), "%")) +
  scale_x_continuous(breaks=seq(-100,100,10)) + 
  scale_y_continuous(breaks=seq(-100,100,10)) + 
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  
  ## butterflies
  geom_raster(
    data=data_KDE1_grid[unpalatability == "low" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE1_grid[unpalatability == "high" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  scale_alpha_continuous(range = c(0, 1),trans = "sqrt", guide="none") +
  geom_path(
    data=data_KDE1,
    aes(x=x, y=y, color=unpalatability, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dashed")) +
  geom_point(aes(x=Dim1, y=Dim2, color=unpalatability, label=species), size=1) +
  geom_image(data=data_roi_butterflies, aes(image=roi_path, x=Dim1, y=Dim2), size=0.1, by="height",
             image_fun=function(img) {img = magick::image_background(img, "none")}) +

  geom_point(data=data_roi_butterflies[pred_success==F],
             aes(x=Dim1, y=Dim2, label=species), size=4, color="yellow") +
  geom_point(data=data_roi_butterflies,
             aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), pch=21, size=2) +
  scale_color_manual("Butterflies - Chemical defense", values=colors_unpala, labels=labels_unpala,
                     guide = guide_legend(keywidth = unit(1.25, "cm"), order=1, override.aes = list(
                       size=3)), drop=T) +
  scale_fill_manual("Butterflies - Chemical defense", values=colors_unpala, labels=labels_unpala,
                    guide = guide_legend(keywidth = unit(1.25, "cm"), order=1, override.aes = list(
                      size=3)), drop=T) +
  scale_linetype_manual("Quantile",  values=c("solid","longdash","dashed"),  labels=c("25%","50%","75%")) + 

  
  ## moths
  new_scale_fill() + new_scale_color() + new_scale("alpha") +
  geom_raster(
    data=data_KDE2_grid[unpalatability == "low" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  geom_raster(
    data=data_KDE2_grid[unpalatability == "high" & density>p25],
    aes(x, y, alpha = density, fill=unpalatability), interpolate = TRUE) +
  scale_alpha_continuous(range = c(0, 1), trans = "sqrt", guide="none") +
  geom_path(
    data=data_KDE2,
    aes(x=x, y=y, color=unpalatability, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dashed")) +
  geom_point(data=pca_coords_moths_agg,
             aes(x=Dim1, y=Dim2, color=unpalatability, label=species), size=1) +
  geom_image(data=data_roi_moths, aes(image=roi_path, x=Dim1, y=Dim2), size=0.1, by="height",
             image_fun=function(img) {img = magick::image_background(img, "none")}) +
  geom_point(data=data_roi_moths[pred_success==F],
             aes(x=Dim1, y=Dim2, label=species), size=4, color="yellow") +
  geom_point(data=data_roi_moths,
             aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), pch=21, size=2) +
  scale_fill_manual("Moths - Strategy", values=colors_moths, labels=labels_moths,
                    guide = guide_legend(keywidth = unit(1.25, "cm"), order=2, override.aes = list(
                      size=3)), drop=T)   +
  scale_color_manual("Moths - Strategy", values=colors_moths, labels=labels_moths,
                     guide = guide_legend(keywidth = unit(1.25, "cm"), order=2, override.aes = list(
                       size=3)), drop=T) +

  ## LDA label
  new_scale_fill() +
  geom_segment(data=pca_quanti_supp[feature_string=="LD1"],
               aes(x = 0, y = 0, xend =  Dim1*loadings_factor,
                   yend =  Dim2*loadings_factor), size=1,
               arrow = arrow(angle = 30, length = unit(2, "mm"),
                             ends = "last", type = "closed")) +
  geom_label(data=pca_quanti_supp[feature_string=="LD1"],
             aes(x =  Dim1*loadings_factor, y =  Dim2*loadings_factor), size=5,
             label="LD1", fill="#d95f02") +
  
  ## theme
  theme(
    strip.text = element_text(size = 14),
    legend.position = "none"
  )


p_leg =
  ggplot(data=data_plot1) +
  coord_fixed() + 
  xlab(paste("Dim1 - Variance explained: ", round(pca_eigenvalues[1,]$var_explained,2), "%")) +
  ylab(paste("Dim2 - Variance explained: ", round(pca_eigenvalues[2,]$var_explained,2), "%")) +
  scale_x_continuous(breaks=seq(-100,100,10)) + 
  scale_y_continuous(breaks=seq(-100,100,10)) + 
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  
  ## butterflies
  geom_point(aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), pch=21, size=3) +
  scale_fill_manual("Butterflies\n(chemical defense)", values=colors_unpala, labels=labels_unpala,
                    guide = guide_legend(order=1, title.position="top"), drop=T) +
  
  
  ## moths
  new_scale_fill() + new_scale_color() + new_scale("alpha") + 
  geom_point(data=pca_coords_moths_agg,
             aes(x=Dim1, y=Dim2, fill=unpalatability, label=species), pch=21, size=3) +
  scale_fill_manual("Moths\n(strategy)", values=colors_moths, labels=labels_moths,
                    guide = guide_legend(order=2, title.position="top"), drop=T)   +
  
  ## quantiles
  geom_path(
    data=data_KDE1,
    aes(x=x, y=y, linetype=factor(level),
        group=interaction(level_idx,unpalatability, group)), linewidth=0.75, color="black") +
  scale_linetype_manual("Spatial density quantiles\n(KDE Dim1/Dim2)",  
                        values=c("solid","longdash","dashed"),  labels=c("25%","50%","75%"),
                        guide = guide_legend(
                          order=3, title.position="top",keywidth=unit(48, "pt"))) +
  
  ## theme
  theme(
    strip.text = element_text(size = 14),
    legend.position = "bottom",
    legend.margin = margin(r=10, unit="pt"),
    legend.box.margin = margin(l=100, unit="pt")
  )
leg <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]
pg_final = cowplot::plot_grid(p1, leg, nrow=2, rel_heights = c(0.9,0.1))

width = 3000
height = 2500
ggsave(paste0("figures/figureS5.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")

# figure S6 - stats -------------------------------------------------------

tree_a = data_tree_sub2

data_lda_sub_summ1_wide_sub_a = data_ld1_agg_dv_sex_wide[species %in% tree_a$tip.label]
data_lda_sub_summ1_wide_sub_a[,species := factor(species, levels = tree_a$tip.label)]
data_lda_sub_summ1_wide_sub_a = data_lda_sub_summ1_wide_sub_a[order(species)]

# 3. Create a correlation structure based on the tree
cor_struct_a <- corBrownian(form=~species, phy = tree_a)
mod1_ld1_diff_a <- gls(diff_abs ~ sign * class_dv, 
               correlation = cor_struct_a, 
               data = data_lda_sub_summ1_wide_sub_a, 
               method = "ML")
summary(mod1_ld1_diff_a)
anova(mod1_ld1_diff_a)
capture.output(anova(mod1_ld1_diff_a), file = "tables/mod_S6a.txt")

data_lda_sub_summ1_wide_sub_b = data_lda_sub_summ1_wide_sub_a[out==T]

tree_b = drop.tip(tree_a, setdiff(tree_a$tip.label, data_lda_sub_summ1_wide_sub_b$species))
data_lda_sub_summ1_wide_sub_b[,species := factor(species, levels = tree_a$tip.label)]
data_lda_sub_summ1_wide_sub_b = data_lda_sub_summ1_wide_sub_b[order(species)]

cor_struct_b <- corBrownian(form=~species, phy = tree_b)
mod1_ld1_diff_b <- gls(diff_abs ~ sign * class_dv, 
                     correlation = cor_struct_b, 
                     data = data_lda_sub_summ1_wide_sub_b, 
                     method = "ML")
summary(mod1_ld1_diff_b)
anova(mod1_ld1_diff_b)
capture.output(anova(mod1_ld1_diff_b), file = "tables/mod_S7b.txt")


# figure S6 --------------------------------------------------------------------

null_dist[, null := T]

p1 =
  ggplot(data_ld1_agg_dv_sex_wide) +
  coord_cartesian(xlim=c(-3.5,4)) + 
  labs(
    title="A",
    y="Density",
    x="Sexual differences in LD1 / male-female \n(positive=male higher, negative=female higher)") + 
  geom_density(aes(fill=class_dv, x=diff), alpha=0.5) +
  scale_fill_manual("", values=c("pink", "gray80"), labels=c("Dorsum", "Ventrum"), guide=guide_legend(order=1)) + 
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) + 
  annotate("text", x=1, y=0.8, hjust=0,size=4,
           label="Higher difference in") + 
  annotate("text", x=1, y=0.65, hjust=0,size=4,
           label="All\nMales: P<0.001\nDorsum: P<0.001") + 
  annotate("text", x=1, y=0.4, hjust=0,size=4,
           label="Outside null\nMales: P<0.001\nDorsum: P<0.001") + 
  geom_vline(xintercept = q, color = "blue", linetype = "dashed", linewidth = 0.6, show.legend=T) +
  geom_vline(xintercept=0) + 
  geom_density(
    data=null_dist, 
    aes(x=diff), color="blue", alpha=0.5, adjust = 10, linetype=2, linewidth = 0.6, show.legend = F) +
    
  geom_line(
    data=data.table(x=0.0001, y=0.0001, col="a"), 
      aes(x=x, y=y, color=col)) +  
    scale_color_manual(values="blue", label="Null\ndistribution",
                       guide=guide_legend(order=2, keywidth=unit(30,"pt"), override.aes=list(
                         linetype=2
                       ))) + 
    
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    # legend.text = element_text(size=8),
    legend.position.inside = c(0.05,0.72)
  )

p2 =
  ggplot(data_ld1_agg_dv_sex_wide[!subfamily %in% subfamily_excl]) +  ggtitle("B") + 
  scale_y_continuous(limits=c(0,2.5)) +
  facet_wrap(~class_dv, labeller = as_labeller(c("dorsal"="Dorsum", "ventral"="Ventrum"))) + ylab("Difference LD1\n") + 
  geom_boxplot(aes(x = subfamily_f, y = diff_abs, fill = sign), 
               position = position_dodge(preserve = "single"), outliers=F) +
  geom_jitter(aes(x = subfamily_f, y = diff_abs, fill = sign, color=out, group = sign), size=0.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  scale_color_manual("Null distribution", values=c("black", "gray60"),labels=c("Inside", "Outside")) + 
  guides(
    color=guide_legend(title.position="top", override.aes = list(size=3)),
    fill=guide_legend(title.position="top")
    ) + 
  scale_fill_manual("Absolute sexual diff. LD1", values=c("sienna3", "deepskyblue"),labels=c("Higher female LD1", "Higher male LD1")) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "bottom")

p3 =
  ggplot(data_ld1_diff_n[!subfamily %in% subfamily_excl]) + ggtitle("C", subtitle="All datapoints") + 
  facet_wrap(~class_dv, labeller = as_labeller(c("dorsal"="Dorsum", "ventral"="Ventrum"))) + ylab("Proportion") + 
  geom_bar(aes(x=subfamily_f, y=prop, fill=sign), position = "stack", stat = "identity") +
  geom_hline(yintercept = 0.5) +
  scale_fill_manual(values=c("sienna3", "deepskyblue")) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none")
p4 =
  ggplot(data_ld1_diff_sub_n[!subfamily %in% subfamily_excl]) + ggtitle("D", subtitle="Outside Null dist. only") + 
  facet_wrap(~class_dv, labeller = as_labeller(c("dorsal"="Dorsum", "ventral"="Ventrum"))) + ylab("Proportion") + 
  geom_bar(aes(x=subfamily_f, y=prop, fill=sign), position = "stack", stat = "identity") +
  geom_hline(yintercept = 0.5) +
  scale_fill_manual(values=c("sienna3", "deepskyblue")) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none")

pg1 = cowplot::plot_grid(p1, p2, align = "hv",rel_widths = c(0.4,0.6))
pg2 = cowplot::plot_grid(p3, p4, align = "hv", axis = "tb", rel_widths = c(0.5,0.5))
pg_f = cowplot::plot_grid(pg1, pg2, ncol=1, axis = "lr", rel_heights = c(0.55,0.45))


width=4000
height=3000
ggsave(paste0("figures/figureS6.png"), pg_f,
       width=width, height=height, units = c("px"), bg="white")


# figure S7a ---------------------------------------------------------------

data_meta_agg[, genus := tstrsplit(species, "_")[[1]]]
data_lda_summ2[, genus := tstrsplit(species, "_")[[1]]]

data_genus_prop = data_lda_summ2[, .N, by = .(genus, unpalatability_pred)]
data_genus_prop = data_genus_prop[order(genus, unpalatability_pred)]
data_genus_prop[, total := sum(N), by = .(genus)]
data_genus_prop = data_genus_prop[total>1]
data_genus_prop[, prop:=N/total]
data_genus_prop = merge(unique(data_meta_agg[,.(genus, subfamily)]), data_genus_prop)
data_genus_prop = data_genus_prop[order(-total, -prop, genus), .SD[1], by = genus]
data_genus_prop[, genus_size := ifelse(total>=15, "large", "small")]

p1 =
  ggplot(data_genus_prop) +
  labs(x="N species in genus", y="Proportion species classified as low/high") + 
  geom_jitter(aes(x=total, y=prop, color=unpalatability_pred),size=2, width=0.01, height=0.01) +
  scale_x_log10(breaks=c(0,1,2,3,4,5,10,20,50,100), limits=c(2,150)) +
  # scale_size_manual(values=c(3, 1)) + 
  geom_vline(xintercept = 18, linetype=2, size=1) + 
  scale_color_manual("Genus-level class prediction\nof LD model:", values=colors_unpala, labels=labels_unpala) +
  theme(
    legend.position = "bottom"
  )


data_genus_prop[total>=15]
data_genus_prop[,sum(total)]
data_genus_prop[, subfamily_f := str_to_title(subfamily)]
data_genus_prop[, genus_f := str_to_title(genus)]

fwrite(
  data_genus_prop[1:30, c("subfamily_f", "genus_f", "unpalatability_pred", "total", "prop")],
  "tables/figureS4b.csv"
)

height=2500
width=2000
ggsave(paste0("figures/figureS7a.png"), p1,
       width=width, height=height, units = c("px"), bg="white")

# figure S7 - plot -------------------------------------------------------

mod_1_phylo_combined <- rbindlist(list(
  mod_1a_phylo_pred_part[, `:=`(aspect = "all", model = "mod_1")],
  mod_1b_phylo_pred_part[, `:=`(aspect = "apo", model = "mod_1")]
), use.names = TRUE)

mod_2_phylo_combined <- rbindlist(list(
  mod_2a_phylo_pred_part[, `:=`(aspect = "all", model = "mod_2")],
  mod_2b_phylo_pred_part[, `:=`(aspect = "apo", model = "mod_2")]
), use.names = TRUE)

mod_3_phylo_combined <- rbindlist(list(
  mod_3a_phylo_pred_part[, `:=`(aspect = "all", model = "mod_3")],
  mod_3b_phylo_pred_part[, `:=`(aspect = "apo", model = "mod_3")]
), use.names = TRUE)

mod_1_combined <- rbindlist(list(
  mod_1a_pred_part[, `:=`(aspect = "all", model = "mod_1")],
  mod_1b_pred_part[, `:=`(aspect = "apo", model = "mod_1")]
), use.names = TRUE)

mod_2_combined <- rbindlist(list(
  mod_2a_pred_part[, `:=`(aspect = "all", model = "mod_2")],
  mod_2b_pred_part[, `:=`(aspect = "apo", model = "mod_2")]
), use.names = TRUE)

mod_3_combined <- rbindlist(list(
  mod_3a_pred_part[, `:=`(aspect = "all", model = "mod_3")],
  mod_3b_pred_part[, `:=`(aspect = "apo", model = "mod_3")]
), use.names = TRUE)

p1 =
  ggplot(mod_1_phylo_combined) + ggtitle("With phylogenetic term") + 
  ylab("Dorso-ventral differences") + 
  coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.69, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=aspect), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, linetype=aspect)) + 
  geom_line(data=mod_1a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_1b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Position:", values=c("#F8766D", "#00BFC4")) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none", 
    axis.title.x = element_blank()
  )

p4 =
  ggplot(mod_1_combined) + ggtitle("Without phylogenetic term") + 
  ylab("Dorso-ventral differences") + 
  coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.69, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=aspect), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, linetype=aspect)) + 
  geom_line(data=mod_1a_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_1b_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Position:", values=c("#F8766D", "#00BFC4")) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none", 
    axis.title = element_blank()
  )


p2 =
  ggplot(mod_2_phylo_combined) + 
  ylab("Inter-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.875, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=interaction(aspect, class_dv)), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, linetype=aspect, color=class_dv)) + 
  geom_line(data=mod_2a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_2b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Position:", values=c("#F8766D", "#7CAE00")) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none", 
    axis.title.x = element_blank()
  )

p5 =
  ggplot(mod_2_combined) + 
  ylab("Inter-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.875, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=interaction(aspect, class_dv)), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, linetype=aspect, color=class_dv)) + 
  geom_line(data=mod_2a_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_2b_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Position:", values=c("#F8766D", "#7CAE00")) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none", 
    axis.title = element_blank()
  )

p3 =
  ggplot(mod_3_phylo_combined) +
  ylab("Intra-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.825, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=interaction(aspect, sex, class_dv)), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, group=interaction(aspect, sex, class_dv), linetype=aspect, color=interaction(sex, class_dv))) + 
  geom_line(data=mod_3a_phylo_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_3b_phylo_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Position:", values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_linetype_manual("Model fit (GAM):", values=c(1,2), guide=guide_legend(order=1)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none",
    legend.key.width = unit(2, "cm"),
    legend.margin=margin(l = 2, unit='cm'),
    legend.box = "vertical",
    legend.box.just = "left"
  )
p6 =
  ggplot(mod_3_combined) +
  ylab("Intra-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.825, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=interaction(aspect, sex, class_dv)), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, group=interaction(aspect, sex, class_dv), linetype=aspect, color=interaction(sex, class_dv))) +
  geom_line(data=mod_3a_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_3b_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_linetype_manual(values=c(1,2)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
    
  )

p_leg =
  ggplot(rbind(mod_3_combined, data.table(class_dv="x",LD1=1, fit=1, aspect="apo",sex="female"), fill=T)) +
  ylab("Intra-sexual similarity") + coord_cartesian(xlim=c(-2.5, 6), ylim = c(0.8, 1)) +
  geom_ribbon(aes(x=LD1, ymin=fit-se, ymax=fit+se, group=interaction(aspect, sex, class_dv)), alpha=0.2) +
  geom_line(aes(x=LD1, y=fit, group=interaction(aspect, sex, class_dv), linetype=aspect, color=interaction(sex, class_dv))) +
  geom_line(data=mod_3a_pred_glob, aes(x=LD1, y=fit)) +
  geom_line(data=mod_3b_pred_glob, aes(x=LD1, y=fit), linetype=2) +
  scale_color_manual("Factor combinations", 
                     values=c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF", "black"),
                     labels=c("Female (Dorsum)", "Male (Dorsum)", 
                              "Female (Ventrum", "Male (Ventrum)", 
                              "Global fit (Fig. 3)"),
                     guide=guide_legend(ncol=2, byrow=T)) +
  scale_linetype_manual("Model fit (GAM)", values=c(1,2), 
                        labels=c("All pattern aspects (Dim 1-177)", "Aposematic pattern aspects (Dim 1-7)"),
                        guide=guide_legend(order=1, keywidth = unit(40, "pt"))) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(20, "pt"),
    legend.margin=margin(l = 40, unit='pt'),
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.box.just = "top",
    axis.title.y = element_blank()
    
  )

fig_leg <- get_plot_component(p_leg, 'guide-box', return_all = TRUE)[[3]]

p_grid = cowplot::plot_grid(p1,p2,p3,p4,p5,p6, ncol=2, align="v", byrow=F)
p_grid = cowplot::plot_grid(p_grid, fig_leg, ncol=1, rel_heights = c(0.9,0.1))

width=3000
height=3500
ggsave(paste0("figures/figureS8.png"), p_grid,
       width=width, height=height, units = c("px"), bg="white")


# figure S9 - prep --------------------------------------------------------

## get tree
tree <- data_tree_sub2

## format heatmaps
data_traits = data_ld1_agg[species %in% tree$tip.label]
traits = data.frame(data_traits[,-"species",with=F])
rownames(traits) = data_traits$species
data_traits_all <- dcast(
  data_ld1_agg_dv_sex[species %in% tree$tip.label],
  species ~ class_dv + sex, value.var = "LD1")
traits_all = data.frame(data_traits_all[,-"species",with=F])
rownames(traits_all) = data_traits_all$species

data_rates_ld1 = data.table(
  LD1_f_d = data_rates[context_trait=="LD1_f_d", prob],
  LD1_m_d = data_rates[context_trait=="LD1_m_d", prob],
  LD1_f_v = data_rates[context_trait=="LD1_f_v", prob],
  LD1_m_v = data_rates[context_trait=="LD1_m_v", prob],
  species = data_rates[context_trait=="LD1_m_d", species]
)

rates_all = data.frame(data_rates_ld1[,1:4])
rownames(rates_all) = data_rates_ld1$species

# keep only species present in both
tree <- drop.tip(tree, setdiff(tree$tip.label, data_traits_all$species))
tree_phylo <- as.phylo(tree)  # strip simmap info

## family labels
species_by_family <- data_meta_agg[species %in% tree$tip.label, .(species_list = list(species)), by = "subfamily"]
species_by_family <- setNames(species_by_family$species_list, species_by_family$subfamily)
tree_phylo = groupOTU(tree_phylo, species_by_family, group_name = "subfamily")

## tribe labels
species_by_tribe <- data_meta_agg[species %in% tree$tip.label, .(species_list = list(species)), by = "tribe"]
species_by_tribe <- setNames(species_by_tribe$species_list, species_by_tribe$tribe)
tree_phylo = groupOTU(tree_phylo, species_by_tribe, group_name = "tribe")

## mod tree data
ladderize = T
data_tree = data.table(ggtree(tree_phylo , ladderize = ladderize)$data)
data_tree = data_tree[order(-y)]
data_tree[, tribe_alternate := rleid(tribe)]
data_tree[, tribe_factor := factor(
  tribe_alternate, labels = rep(c("A", "B"), length.out = max(tribe_alternate)))]
tribe_values = data.frame(data_tree$tribe_factor)
rownames(tribe_values) = data_tree$species
data_tree[, subfamily := factor(subfamily, levels=c(
  "nymphalinae", "cyrestinae", "biblidinae", "apaturinae", "pseudergolinae", 
  "heliconiinae",  "limenitidinae", "satyrinae", "charaxinae", "libytheinae", "danainae"))]


# figure S9 - plot --------------------------------------------------------

## subset for images
data_tree_sub = data_tree[isTip == TRUE,]
data_tree_sub[, species := label]
data_tree_sub[, roi_path_f_d := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/", "dorsal_female", "/", species, ".png")]
data_tree_sub[, roi_path_m_d := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/", "dorsal_male", "/", species, ".png")]
data_tree_sub[, roi_path_f_v := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/", "ventral_female", "/", species, ".png")]
data_tree_sub[, roi_path_m_v := paste0("data_raw/segmentation_masks_clean/centroids_dv_sex/", "ventral_male", "/", species, ".png")]

data_tree_sub = data_tree_sub[order(y)]
y_sel = seq(6, max(data_tree_sub$y), 25)

data_tree_sub = data_tree_sub[y %in% y_sel]
data_tree_sub[, label_f := str_to_sentence(str_replace(as.character(species), "_", " "))]
data_tree_sub[, N := factor(1:.N)]
data_tree_sub[, y_img := y]

img_size = 0.04

p1 =
  ggtree(tree_phylo, ladderize = ladderize) +
  geom_tree(aes(colour = subfamily), size = 0.5) + 
  coord_cartesian(clip="off") + 
  scale_color_manual(values=fam_cols,guide="none") +
  geom_segment(data = data_tree_sub,
               aes(x=x, xend=x+105, y=y, yend=y, color=subfamily),size=0.5) +
  geom_tiplab(data=data_tree_sub, aes(x=x+60, y=y_img, image=roi_path_f_d), geom="image",
              size=img_size, inherit.aes = FALSE, angle=0) +
  geom_tiplab(data=data_tree_sub, aes(x=x+75, y=y_img, image=roi_path_m_d), geom="image",
              size=img_size, inherit.aes = FALSE, angle=0) +
  geom_tiplab(data=data_tree_sub, aes(x=x+90, y=y_img, image=roi_path_f_v), geom="image",
              size=img_size, inherit.aes = FALSE, angle=0) +
  geom_tiplab(data=data_tree_sub, aes(x=x+105, y=y_img, image=roi_path_m_v), geom="image",
              size=img_size, inherit.aes = FALSE, angle=0) +
  
  theme_void() + theme(
    legend.position = "none",
    plot.margin = margin(t=100,b=10,l=15,r=50, unit="pt")
  ) 

p1 <- p1 + new_scale_fill()
p1 = gheatmap(p1, unpalatability_labs, offset=-1.75,  width=0.06, colnames = F) + 
  scale_fill_manual(values=colors_unpala_tree, guide="none", na.value = colors_unpala["unknown"]) 

p1 <- p1 + new_scale_fill()
p1 <- gheatmap(p1, traits_all, width = 0.25, offset = 4,
               colnames = T, colnames_position = "top", 
               font.size = 4, family="bold",
               colnames_angle=90, hjust=0, colnames_offset_y = 5,
               custom_column_labels=c("Female dorsum", "Male dorsum", "Female ventrum", "Male ventrum")) +
  scale_fill_viridis_c(
    option="inferno", oob = scales::oob_squish,
    limits = c(-2.5,6.5), guide="none") 

p1 <- p1 + new_scale_fill()
p1 <- gheatmap(p1, rates_all, width = 0.25, offset = 26,
               colnames = T, colnames_position = "top", 
               font.size = 4, family="bold",
               colnames_angle=90, hjust=0, colnames_offset_y = 5,
               custom_column_labels=c("Female dorsum", "Male dorsum", "Female ventrum", "Male ventrum")) +
  scale_fill_viridis_c(guide = "none") +   
  scale_x_continuous(expand=c(0,0)) 

p_leg1 =
  ggplot(data_tree_sub) +
  
  geom_ribbon(data=data_tree[!subfamily==0],
              aes(x=1,y=1, ymin=1, ymax=1, fill=subfamily)) +
  scale_fill_manual(values=fam_cols, "Subfamily", guide=guide_legend(
    order=1, override.aes = list(linewidth=1), keywidth = unit(5, "pt"),
    title.position = "top", nrow=3, byrow=T),
    labels=fam_labs) +
  
  new_scale_color() +
  geom_line(data=data_meta_agg,
            aes(x=1,y=1,color=unpalatability), linewidth=1) +
  scale_color_manual(
    "Chemical defense labels",
    values=colors_unpala,
    labels=labels_unpala,
    na.value = colors_unpala["unknown"],
    guide=guide_legend(
      order=2,ncol=3,byrow=T,
      theme = theme(legend.margin = margin(l=-10,unit="pt"),
                    legend.text = element_text(size=10)),
      title.position = "top")) +
  
  new_scale_color() +
  geom_point(data=melt(rates_all),
             aes(x=1,y=1,color=value)) +
  scale_color_viridis_c(
    "Aposematic coloration (LD1)",
    option="inferno",
    oob = scales::oob_squish,
    limits = c(-2.5,6.5),
    guide=guide_colorbar(
      order=3, keyheight=unit(0.5, "cm"), keywidth=unit(6,"cm"),
      title.position = "top", direction = "horizontal")) +

  new_scale_color() +
  geom_point(data=melt(traits_all),
             aes(x=1, y=1, color=value)) +
  scale_color_viridis_c(
    "Prob. higher rate",
    breaks=c(0,0.5,1), limits=c(0,1),
    guide=guide_colorbar(
      order=4, keyheight=unit(0.5, "cm"), keywidth=unit(4,"cm"),
      title.position = "top", direction = "horizontal")) +

  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(r=20, l=20, unit="pt"))


leg1 <- get_plot_component(p_leg1, 'guide-box', return_all = TRUE)[[3]]
pg_final = cowplot::plot_grid(p1, leg1, ncol=1, rel_heights = c(0.925,0.075))


height=6000
width=4000
ggsave(paste0("figures/figureS9.png"), pg_final,
       width=width, height=height, units = c("px"), bg="white")


# figure S10 ---------------------------------------------------------------

data_fail = copy(data_lda_all[species %in% data_tree_sub4$tip.label & !sex=="unknown" & unpalatability=="low"])
data_fail[, pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]
data_fail_count <- data_fail[
  pred_success == FALSE, .(n_specimens_fail = uniqueN(uuid)), by = c("subfamily","species", "class_dv", "sex")]
data_fail_count = data_fail_count[order(species, class_dv, sex)]

data_fail_wide1 = dcast(data_fail_count, subfamily+species+class_dv~sex, value.var = "n_specimens_fail")
data_fail_wide1[is.na(data_fail_wide1)] = 0
data_fail_wide1[, female_prop := female / (female + male)]
data_fail_wide1[, subfamily_f := str_to_sentence(subfamily)]
data_fail_wide1[, species_f := str_to_sentence(str_replace(as.character(species), "_", " "))]

p1 = 
ggplot(data_fail_wide1) + 
  facet_wrap(~class_dv, ncol=1, labeller=as_labeller(c("dorsal"="Dorsum", "ventral"="Ventrum"))) +
  labs(y="Count", x="Female proportion of missclassified specimens") + 
  geom_histogram(aes(x=female_prop, fill=subfamily), boundary=0) +
  geom_vline(xintercept=0.5) + 
  scale_fill_manual("Subfamily", values=fam_cols, labels=fam_labs) +
  scale_y_continuous(breaks=seq(0,10,1)) +
  theme(
    legend.title.position = "top",
    legend.position = "bottom"
  )

height=2500
width=2500
ggsave(paste0("figures/figureS10.png"), p1,
       width=width, height=height, units = c("px"), bg="white")

# table S1 - stats --------------------------------------------------------

trait = setNames(
  data_ld1_agg$LD1,
  data_ld1_agg$species)

trait_d = setNames(
  data_ld1_agg_dv[class_dv=="dorsal"]$LD1,
  data_ld1_agg_dv[class_dv=="dorsal"]$species)
trait_v = setNames(
  data_ld1_agg_dv[class_dv=="ventral"]$LD1,
  data_ld1_agg_dv[class_dv=="ventral"]$species)

trait_d_f = setNames(
  data_ld1_agg_dv_sex[class_dv=="dorsal" & sex=="female"]$LD1,
  data_ld1_agg_dv_sex[class_dv=="dorsal" & sex=="female"]$species)
trait_d_m = setNames(
  data_ld1_agg_dv_sex[class_dv=="dorsal" & sex=="male"]$LD1,
  data_ld1_agg_dv_sex[class_dv=="dorsal" & sex=="male"]$species)
trait_v_f = setNames(
  data_ld1_agg_dv_sex[class_dv=="ventral" & sex=="female"]$LD1,
  data_ld1_agg_dv_sex[class_dv=="ventral" & sex=="female"]$species)
trait_v_m = setNames(
  data_ld1_agg_dv_sex[class_dv=="ventral" & sex=="male"]$LD1,
  data_ld1_agg_dv_sex[class_dv=="ventral" & sex=="male"]$species)

## labels only 
tree = data_tree_sub3

phylosig(tree,trait,test=TRUE)

phylosig(tree,trait_d,test=TRUE)
phylosig(tree,trait_v,test=TRUE)

phylosig(tree,trait_d_f,test=TRUE)
phylosig(tree,trait_d_m,test=TRUE)
phylosig(tree,trait_v_f,test=TRUE)
phylosig(tree,trait_v_m,test=TRUE)

## full trees
tree = data_tree_sub1

phylosig(tree,trait,test=TRUE)

phylosig(tree,trait_d,test=TRUE)
phylosig(tree,trait_v,test=TRUE)

phylosig(tree,trait_d_f,test=TRUE)
phylosig(tree,trait_d_m,test=TRUE)
phylosig(tree,trait_v_f,test=TRUE)
phylosig(tree,trait_v_m,test=TRUE)

# table S5 ----------------------------------------------------------------

data_fail = copy(data_lda_all)
data_fail[unpalatability %in% c("low","high"), pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]
data_fail_count <- data_fail[
  pred_success == FALSE, .(n_specimens_fail = uniqueN(uuid)), by = species]

data_fail_count[species %like% "sephisa_dichroa"]

data_fail = data_lda_summ3[pred_success==F & unpalatability %in% c("low","high")]
fail_table_d = data_lda_all[class_dv=="dorsal" & species %in% data_fail[class_dv=="dorsal"]$species, 
                         table(as.character(species), as.character(sex))]
fail_table_v = data_lda_all[class_dv=="dorsal" & species %in% data_fail[class_dv=="ventral"]$species, 
                            table(as.character(species), as.character(sex))]
data_fail_table_d <- data.table(
  class_dv = "dorsal",
  species = rownames(fail_table_d),
  female = fail_table_d[, "female"],
  male   = fail_table_d[, "male"],
  unknown   = fail_table_d[, "unknown"]
)
data_fail_table_v <- data.table(
  class_dv = "ventral",
  species = rownames(fail_table_v),
  female = fail_table_v[, "female"],
  male   = fail_table_v[, "male"],
  unknown   = fail_table_v[, "unknown"]
)
data_fail_table = rbind(data_fail_table_d, data_fail_table_v)

data_fail_table_wide = dcast(
  data_fail_table, species~class_dv, 
  value.var = c("female","male","unknown"),
  fill = 0)
data_fail_table_wide = merge(data_meta_agg[,.(subfamily, species, unpalatability, n_specimens)], data_fail_table_wide)
data_fail_table_wide = merge(data_fail_table_wide, data_fail_count)
data_fail_table_wide[, subfamily_f := str_to_sentence(subfamily)]
data_fail_table_wide[, species_f := str_to_sentence(str_replace(as.character(species), "_", " "))]
data_fail_table_wide[, ratio_fail := paste0(n_specimens, "/", n_specimens_fail)]
data_fail_table_wide[, ratio_d := paste0(female_dorsal, "/", male_dorsal, "/",  unknown_dorsal)]
data_fail_table_wide[, ratio_v := paste0(female_ventral, "/", male_ventral, "/",  unknown_ventral)]

data_fail_table_wide = data_fail_table_wide[order(unpalatability, subfamily, species)]
data_fail_table_wide_f = data_fail_table_wide[,c(
  "subfamily_f", "species_f", "unpalatability", 
  "ratio_fail", "ratio_d", "ratio_v")]

fwrite(data_fail_table_wide_f, "tables/tableS5.csv")

# table S6 ----------------------------------------------------------------

data_fail = copy(data_lda_all[species %in% data_tree_sub4$tip.label & !sex=="unknown" & unpalatability=="low"])
data_fail[, pred_success := ifelse(
  as.character(unpalatability)==as.character(unpalatability_pred), T, F)]
data_fail_count1 <- data_fail[
  pred_success == FALSE, .(n_specimens_fail = uniqueN(uuid)), by = c("subfamily","species", "class_dv", "sex")]
data_fail_count1 = data_fail_count1[order(subfamily, species, class_dv, sex)]

data_fail_count2 <- data_fail[
  pred_success == FALSE, .(n_specimens_fail = uniqueN(uuid)), by = c("subfamily","species")]
data_fail_count2 = data_fail_count2[order(subfamily, species)]

data_fail_wide2 = dcast(data_fail_count1, subfamily+species~class_dv+sex, value.var = "n_specimens_fail")
data_fail_wide2[is.na(data_fail_wide2)] = 0
data_fail_wide2[, female_prop_d := dorsal_female / (dorsal_female + dorsal_male)]
data_fail_wide2[, female_prop_v := ventral_female / (ventral_female + ventral_male)]
data_fail_wide2[, ratio_d := paste0(dorsal_female, "/", dorsal_male)]
data_fail_wide2[, ratio_v := paste0(ventral_female, "/", ventral_male)]

data_fail_wide2 = merge(data_fail_wide2, data_meta_agg[, .(species, n_specimens)])
data_fail_wide2 = merge(data_fail_wide2, data_fail_count2[, .(species, n_specimens_fail)])
data_fail_wide2[, ratio_fail := paste0(n_specimens, "/", n_specimens_fail)]
data_fail_wide2[, species_misscl := ifelse(
  species %in% data_fail_table_wide[unpalatability=="low"]$species, T, F)]

data_fail_wide2[, subfamily_f := str_to_sentence(subfamily)]
data_fail_wide2[, species_f := str_to_sentence(str_replace(as.character(species), "_", " "))]
data_fail_wide2 = data_fail_wide2[order(subfamily, species)]
data_fail_wide2[is.na(data_fail_wide2)] = "-"
fwrite(data_fail_wide2[, .(
  subfamily_f, species_f, 
  ratio_fail, species_misscl,
  ratio_d, female_prop_d , ratio_v, female_prop_v)], "tables/tableS6.csv")
