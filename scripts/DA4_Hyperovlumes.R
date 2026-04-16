#' ####################################################################### #
#' PROJECT: [DA.4 - Hypervolumes]
#' CONTENTS:
#'  - Calculating hypervolumes and comparing their sizes
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls())

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec1 <- c(
  "hypervolume",
  "car",
  "multcomp",
  "ggpubr",
  "emmeans",
  "lmerTest",
  "magrittr",
  "cowplot",
  "grid",
  "gridExtra",
  "VennDiagram",
  "pbapply",
  "parallel",
  "dplyr",
  "stringr",
  "tidyverse",
  "maditr",
  "polyclip",
  "png"
)
sapply(package_vec1, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

## Plotting Theme Set-Up --------------------------------------------------
# Basic plotting theme
my_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Color (short for 'colour') scheme
# 3x2 scale (PFGs x species -- for linear models)
# Forb species (blues): H. umbellata = "#016392"; L. orbiculata = "#A0CBE8"
# Graminoid species (orange-ish?): P. bonplandianum = "#E19825"; R. macrochaeta = "#F7C480"
# Woody species (greens): G. glomerata = "#3E8853"; V. floribundum = "#9FCD99"
pal_lm <- c("#016392", "#A0CBE8", "#E19825", "#F7C480", "#3E8853", "#9FCD99")
names(pal_lm) <- c("H. umbellata", "L. orbiculata", "P. bonplandianum", "R. macrochaeta", "G. glomerata", "V. floribundum")
# scales::show_col(pal_lm)
FG_pal <- c(
  colorRampPalette(pal_lm[1:2])(3)[2],
  colorRampPalette(pal_lm[3:4])(3)[2],
  colorRampPalette(pal_lm[5:6])(3)[2]
)
names(FG_pal) <- c("Forb", "Graminoid", "Woody")

# DATA LOADING ============================================================
source("scripts/0_data_import.R") # sourcing data import script
traits_wide$ID <- paste(traits_wide$site, traits_wide$plot_id, traits_wide$individual_nr, traits_wide$taxon, sep = "_")

# POOL ELEVATION WITHIN SITE. WAY: 3101; ACJ: 3468; TRE: 3715
traits_df <- na.omit(traits_wide[, c(-2:-5, -17:-19)])
traits_df$elevation <- ifelse(traits_df$site == "ACJ", 3468,
  ifelse(traits_df$site == "WAY", 3101,
    ifelse(traits_df$site == "TRE", 3715, NA)
  )
)
# CREATE DATA FRAME OF GROUP MEMBERSHIP OF EACH SPECIES
Grouping_df <- as.data.frame(table(traits_df[, c("taxon", "functional_group")])) # tally up each species and functional group combination
Grouping_df <- Grouping_df[Grouping_df$Freq != 0, -3] # remove every non-observed pairing
Grouping_df$taxon <- as.character(Grouping_df$taxon) # ensure that taxon column is character

# ANALYSES ================================================================

## Functionality ----------------------------------------------------------
## CALCULATE HYPERVOLUMES FOR DESIRED GROUPINGS OF DATA
FUN.Hypervolumes <- function(data = traits_df, # the data frame with traits and the grouping column
                             TraitCols = 8:12, # which columns of the data frame contain the trait values, default excludes plant height due to identical plant height for each individual and its leaves and wet mass
                             Grouping = "family" # by which column to establish groups
) {
  cl <- parallel::detectCores() / 2
  on.exit(stopCluster(cl))
  cl <- parallel::makeCluster(cl) # for parallel pbapply functions
  parallel::clusterExport(cl,
    varlist = c(
      "data",
      "Grouping",
      "TraitCols"
    ),
    envir = environment()
  )

  data_scale <- log(data[, TraitCols]) # now applying logarithmic transformation
  data_ls <- split(data_scale, data[, Grouping]) # split the data by grouping argument
  hv_ls <- pblapply(data_ls, cl = cl, FUN = function(x) {
    if (0 %in% apply(apply(x, 2, range), 2, diff)) {
      return(NA)
    } else {
      return(hypervolume::hypervolume(x))
    }
  }) # calculate hypervolume for each group separately
  return(hv_ls) # return list of volumes with one element for each group
  # volumes_ls <- lapply(hv_ls, get_volume)
  # return(volumes_ls)
}

## CALCULATE OVERLAP OF HYPERVOLUMES
FUN.Overlap <- function(data, # the list of hypervolumes from which to pull
                        names, # the names of the hypervolumes which to compare
                        what = "jaccard") {
  overlap <- hypervolume_set(data[[names[1]]], data[[names[2]]], check.memory = FALSE) # extract hypervolumes and safe as set
  if (what == "jaccard") {
    overlap <- hypervolume::hypervolume_overlap_statistics(overlap)[1] # jaccard distance calculation
  } else {
    overlap <- as.numeric(get_volume(overlap[["Intersection"]]))
  }
  return(overlap) # return numeric overlap
}

# Function to get post-hoc test stats
get_posthoc_stats <- function(model, group_var) {
  cmp <- do.call(mcp, setNames(list("Tukey"), group_var))

  stats <- glht(model, linfct = cmp) %>%
    summary(test = adjusted("bonferroni"))

  stats.df <- cbind(
    "Estimate" = stats$test$coefficients,
    "Std_Error" = stats$test$sigma,
    "z_value" = stats$test$tstat,
    "p_value" = stats$test$pvalues
  ) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Comparison")

  return(stats.df)
}

# Function to get compact letter display of post-hoc groups
get_letters <- function(model, group_var) {
  cmp <- do.call(mcp, setNames(list("Tukey"), group_var))

  letters <- glht(model, linfct = cmp) %>%
    cld() %>%
    .[[10]] %>%
    .[[1]] %>%
    toupper() %>%
    as.list()

  return(letters)
}

## Creating Hypervolumes --------------------------------------------------
message("Hypervolume calculation")
if (!file.exists(file.path("./data", "Hypervolumes.RData"))) {
  ### Individual Hypervolumes ####
  ID_hv <- FUN.Hypervolumes(Grouping = "ID")
  ID_vols <- lapply(ID_hv, get_volume)

  ### Species Hypervolumes ####
  taxon_hv <- FUN.Hypervolumes(Grouping = "taxon", TraitCols = c(6, 7:12)) # also including plant height now
  taxon_vols <- lapply(taxon_hv, get_volume)

  ### Functional Group Hypervolumes ####
  functional_group_hv <- FUN.Hypervolumes(Grouping = "functional_group", TraitCols = c(6, 7:12)) # also including plant height now
  functional_group_vols <- lapply(functional_group_hv, get_volume)

  ### Saving Hypervolume Lists ####
  save(ID_hv, ID_vols,
    taxon_hv, taxon_vols,
    functional_group_hv, functional_group_vols,
    file = file.path("./data", "Hypervolumes.RData")
  )
} else {
  load(file.path("./data", "Hypervolumes.RData"))
}

## Preparing data ----------------------------------------------------------
vols_ls <- list(
  Individuals = list(
    vols = ID_vols,
    hv = ID_hv
  ),
  Taxon = list(
    vols = taxon_vols,
    hv = taxon_hv
  ),
  Functional = list(
    vols = functional_group_vols,
    hv = functional_group_hv
  )
)
vols_ls <- lapply(vols_ls, FUN = function(x) {
  lapply(x, FUN = function(y) {
    names(y) <- str_replace_all(names(y), c(
      c(
        "Halenia umbellata" = "H. umbellata",
        "Lachemilla orbiculata" = "L. orbiculata",
        "Paspalum bonplandianum" = "P. bonplandianum",
        "Rhynchospora macrochaeta" = "R. macrochaeta",
        "Gaultheria glomerata" = "G. glomerata",
        "Vaccinium floribundum" = "V. floribundum"
      )
    ))
    y
  })
})


Vols_df <- data.frame(
  values = unlist(ID_vols),
  ID = names(unlist(ID_vols))
)
Vols_df$ID <- gsub(Vols_df$ID, pattern = ".untitled", replacement = "")
plot_df <- plyr::join(x = Vols_df, y = traits_df, by = "ID") %>%
  distinct(values, ID, .keep_all = TRUE)
# %>%
#   filter(values != max(.$values)) #  and filter out the outlier which is TRE_1_NA_Rhynchospora macrochaeta
### factorise columns that need to be factors
plot_df$taxon <- str_replace_all(plot_df$taxon, c(
  c(
    "Halenia umbellata" = "H. umbellata",
    "Lachemilla orbiculata" = "L. orbiculata",
    "Paspalum bonplandianum" = "P. bonplandianum",
    "Rhynchospora macrochaeta" = "R. macrochaeta",
    "Gaultheria glomerata" = "G. glomerata",
    "Vaccinium floribundum" = "V. floribundum"
  )
))
plot_df$taxon <- factor(plot_df$taxon, levels = c("H. umbellata", "L. orbiculata", "P. bonplandianum", "R. macrochaeta", "G. glomerata", "V. floribundum"))
plot_df$functional_group <- as.factor(plot_df$functional_group)
plot_df$site %<>% factor(levels = c("WAY", "ACJ", "TRE"))

## Sizes of Hypervolumes --------------------------------------------------
message("#### Sizes of Hypervolumes")
### Taxon -----
print("Volume by taxon")
### set comparisons for plotted boxplots
taxon.comps <- list(
  c("G. glomerata", "V. floribundum"),
  c("H. umbellata", "L. orbiculata"),
  c("P. bonplandianum", "R. macrochaeta")
)
sizes_indiv_taxon_gg <- ggplot(plot_df, aes(y = log(values), x = taxon, fill = taxon)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_lm) +
  stat_compare_means(comparisons = taxon.comps, method = "t.test", label = "p.signif") +
  my_theme +
  labs(y = "Log-Transformed Hypervolume Size", x = "", fill = "Species", title = "Hypervolumes of Individuals") +
  # facet_wrap(~site) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

sizes_total_taxon_gg <- ggplot(
  data.frame(
    Groups = names(vols_ls[["Taxon"]]$vols),
    Value = unlist(vols_ls[["Taxon"]]$vols)
  ),
  aes(x = factor(Groups, levels = names(pal_lm)), y = Value, fill = factor(Groups))
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_lm) +
  labs(x = "", y = "Hypervolume Size", title = "Total Hypervolumes") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )

sizes_taxon_gg <- plot_grid(sizes_indiv_taxon_gg + theme(legend.position = "none"),
  sizes_total_taxon_gg,
  get_legend(sizes_indiv_taxon_gg),
  rel_widths = c(1, 1, 0.5), ncol = 3
)
sizes_taxon_gg

### Functional Groups -----
print("Volume by functional group")
fg.comps <- list(c("Forb", "Graminoid"), c("Graminoid", "Woody"), c("Forb", "Woody"))
sizes_indiv_fg_gg <- ggplot(plot_df, aes(y = log(values), x = functional_group, fill = functional_group)) +
  geom_boxplot() +
  scale_fill_manual(values = FG_pal) +
  stat_compare_means(comparisons = fg.comps, method = "t.test", label = "p.signif", size = 5) +
  my_theme +
  labs(y = "Log-Transformed Hypervolume Size", x = "", fill = "Functional Group", title = "Hypervolumes of Individuals") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

sizes_total_fg_gg <- ggplot(
  data.frame(
    Groups = names(vols_ls[["Functional"]]$vols),
    Value = unlist(vols_ls[["Functional"]]$vols)
  ),
  aes(x = Groups, y = Value, fill = factor(Groups))
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = FG_pal) +
  labs(x = "", y = "Hypervolume Size", title = "Total Hypervolumes") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )

sizes_fg_gg <- plot_grid(sizes_indiv_fg_gg + theme(legend.position = "none"),
  sizes_total_fg_gg,
  get_legend(sizes_indiv_fg_gg),
  rel_widths = c(1, 1, 0.5), ncol = 3
)
sizes_fg_gg

### Manuscript Plot ----
sizes_gg <- plot_grid(sizes_fg_gg, sizes_taxon_gg, ncol = 1, labels = "AUTO")
ggsave(file.path("plots", "HV_Sizes.png"), sizes_gg, units = "in", height = 10, width = 12, dpi = 600)
ggsave(file.path("plots", "HV_Sizes.pdf"), sizes_gg, units = "in", height = 10, width = 12, dpi = 600)

### Elevation Effects ----
#### Plot ----
print("Volume by elevation")
elev_gg <- ggplot(plot_df, aes(y = values, x = elevation, color = taxon)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F) +
  scale_color_manual(values = pal_lm) +
  my_theme +
  ylab("Hypervolume Size") +
  xlab("Elevation (mASL)") +
  labs(color = "Species") +
  theme(
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
elev_gg
ggsave(file.path("plots", "Elevation.png"), elev_gg, units = "in", height = 4.7, width = 7.4, dpi = 600)
ggsave(file.path("plots", "Elevation.pdf"), elev_gg, units = "in", height = 4.7, width = 7.4, dpi = 600)

head(plot_df[order(plot_df$values, decreasing = TRUE), ])

#### Linear Mixed Effect Models ----
print("Hypervolume linear mixed effect model")
### taxonomic grouping of species
fg <- traits_wide %>%
  mutate(taxon = str_replace(taxon, " ", ".")) %>%
  distinct(taxon, functional_group, family)
### elevation at sites
el.values <- traits_wide %>%
  group_by(site) %>%
  summarise(elevation = mean(elevation))
### model data frame containing volume size, study compartments (site, plot, individual), taxonomic memembership, and elevation
temp <- Vols_df %>%
  separate(ID, c("site", "plot_id", "individual_nr", "taxon"), sep = "_") %>%
  mutate(
    taxon = str_replace(taxon, " ", "."),
    plot_id = factor(plot_id),
    individual_nr = as.numeric(individual_nr)
  ) %>%
  left_join(fg, by = "taxon") %>%
  left_join(el.values, by = "site") %>%
  na.omit() # THIS WILL REMOVE THE OUTLIERS
temp$site %<>% factor(levels = c("WAY", "ACJ", "TRE"))
temp$individual_uid <- paste(temp$site, temp$plot_id, temp$individual_nr, temp$taxon, sep = "_")

# THESE ARE THE FINAL MODELS as of April 11, 2024
print("taxon model ----")
hv.taxon.model <- lmerTest::lmer(
  values ~ 0 + scale(elevation) * taxon + (1 | site),
  data = temp
)
print(car::Anova(hv.taxon.model, type = 3))

taxon.comparisons <- emmeans(hv.taxon.model, list(pairwise ~ taxon), adjust = "tukey")
print(taxon.comparisons)

print("functional group model ----")
hv.fg.model <- lmerTest::lmer(
  values ~ 0 + scale(elevation) * functional_group + (1 | site),
  data = temp
)
print(car::Anova(hv.fg.model, type = 3))

fg.comparisons <- emmeans(hv.fg.model, list(pairwise ~ functional_group), adjust = "tukey")
print(fg.comparisons)

## Overlap of Hypervolumes ------------------------------------------------
### Functional Group -----
Overlaps_fg <- as.data.frame(do.call(rbind, combn(x = names(vols_ls[["Functional"]]$vols), m = 2, simplify = FALSE)))
colnames(Overlaps_fg) <- c("G1", "G2")
Overlap <- apply(Overlaps_fg, 1, FUN = function(y) {
  FUN.Overlap(
    data = vols_ls[["Functional"]]$hv,
    names = y
  )
})
Overlaps_fg$Value <- Overlap * 100

OV_fg <- ggplot(Overlaps_fg, aes(x = G1, y = G2, label = round(Value, 2))) +
  geom_label(fill = "white") +
  geom_tile(aes(fill = Value)) +
  scale_fill_viridis_c(direction = -1, option = "E", begin = 0.2) + # limits = c(0, 10)
  geom_label() +
  labs(x = "", y = "", fill = "Jaccard Distances", title = "Functional Groups (Jaccard Distances)") +
  my_theme +
  theme(
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(dev.size()[1] / 6, "inches")
  )
print(OV_fg)

### Taxon ----
#### Venn Diagramms -----
Overlaps_venn <- data.frame(
  G1 = c("H. umbellata", "P. bonplandianum", "G. glomerata"),
  G2 = c("L. orbiculata", "R. macrochaeta", "V. floribundum"),
  Value = c(
    FUN.Overlap(
      data = taxon_hv,
      names = Grouping_df$taxon[Grouping_df$functional_group == "Forb"],
      what = "absolute"
    ),
    FUN.Overlap(
      data = taxon_hv,
      names = Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"],
      what = "absolute"
    ),
    FUN.Overlap(
      data = taxon_hv,
      names = Grouping_df$taxon[Grouping_df$functional_group == "Woody"],
      what = "absolute"
    )
  )
)

Forbs_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`H. umbellata`), 4),
  round(as.numeric(vols_ls$Taxon$vols$`L. orbiculata`), 4),
  round(Overlaps_venn$Value[1], 4),
  c("H. umbellata", "L. orbiculata"),
  fill = pal_lm[1:2],
  alpha = 1,
  inverted = TRUE,
  ind = FALSE,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0,
  cat.dist = 0.04,
  ext.dist = -0.15,
  ext.length = 0.85
)
# Forbs_euler <- ggdraw(Forbs_euler)
Forbs_euler

Gram_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`P. bonplandianum`), 2),
  round(as.numeric(vols_ls$Taxon$vols$`R. macrochaeta`), 2),
  round(Overlaps_venn$Value[2], 4),
  c("P. bonplandianum", "R. macrochaeta"),
  fill = pal_lm[3:4],
  alpha = 1,
  ind = FALSE,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 180,
  cat.dist = 0.04,
  ext.dist = -0.15,
  ext.length = 0.85
)
# Gram_euler <- ggdraw(Gram_euler)

Wood_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`V. floribundum`), 2),
  round(as.numeric(vols_ls$Taxon$vols$`G. glomerata`), 2),
  round(Overlaps_venn$Value[3], 4),
  c("V. floribundum", "G. glomerata"),
  fill = pal_lm[6:5],
  alpha = 1,
  ind = FALSE,
  inverted = TRUE,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0,
  cat.dist = 0.04,
  ext.dist = -0.15,
  ext.length = 0.85
)
# Wood_euler <- ggdraw(Wood_euler)

# OV_venn <- plot_grid(Forbs_euler, Gram_euler, Wood_euler, ncol = 3)
# ggsave(file.path("plots", 'HV_Over_Venn.png'), OV_venn, units = 'in', height = 4.7, width = 14, dpi = 600)
png(file.path("plots", "HV_Over_Venn.png"), units = "in", height = 4.7, width = 14, res = 600)
par(mfrow = c(1, 3), mai = c(0, 0, 0, 0))
for (EulerIter in 1:3) {
  vp <- list(Forbs_euler, Gram_euler, Wood_euler)[[EulerIter]]
  Eulerpal <- pal_lm

  A <- list(list(x = as.vector(vp[[3]][[1]]), y = as.vector(vp[[3]][[2]])))
  B <- list(list(x = as.vector(vp[[4]][[1]]), y = as.vector(vp[[4]][[2]])))

  AintB <- polyclip(A, B)
  ix <- sapply(vp, function(x) grepl("text", x$name, fixed = TRUE))
  labs <- do.call(rbind.data.frame, lapply(vp[ix], `[`, c("x", "y", "label")))
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  polygon(A[[1]], col = Eulerpal[1 + (EulerIter - 1) * 2])
  polygon(B[[1]], col = Eulerpal[2 + (EulerIter - 1) * 2])
  polygon(AintB[[1]], col = "darkgrey")
  text(x = labs$x, y = labs$y, labels = labs$label)
}
dev.off()

imgurlBase <- file.path("plots", "HV_Over_Venn.png")
img <- readPNG(imgurlBase)
g <- rasterGrob(img, interpolate = TRUE)

OV_venn <- ggplot() +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()


#### Relative Overlap -----
cl <- parallel::detectCores() / 2
cl <- parallel::makeCluster(cl) # for parallel pbapply functions
parallel::clusterExport(cl,
  varlist = c(
    "vols_ls",
    "install.load.package",
    "package_vec1",
    "FUN.Overlap"
  ),
  envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec1, install.load.package))

Overlaps_taxon <- as.data.frame(expand.grid(names(pal_lm), names(pal_lm)))
colnames(Overlaps_taxon) <- c("G1", "G2")
Overlap <- pbapply(Overlaps_taxon, 1,
  cl = cl,
  FUN = function(y) {
    FUN.Overlap(
      data = vols_ls[["Taxon"]]$hv,
      names = y
    )
  }
)
stopCluster(cl)
closeAllConnections()
Overlaps_taxon$Value <- Overlap * 100

casted_OV <- dcast(Overlaps_taxon, G1 ~ G2)
G1 <- casted_OV[, 1]
casted_OV <- casted_OV[, -1]
casted_OV[lower.tri(casted_OV, diag = TRUE)] <- NA
casted_OV <- melt(casted_OV)
colnames(casted_OV) <- c("G2", "Value2")
casted_OV$G1 <- rep(G1$G1, 6)
Overlaps_taxon <- base::merge(Overlaps_taxon, casted_OV)

OV_taxon <- ggplot(Overlaps_taxon, aes(x = G1, y = G2, label = round(Value2, 2))) +
  geom_label(fill = "white") +
  geom_tile(aes(fill = Value2)) +
  scale_fill_viridis_c(direction = -1, option = "E", begin = 0.2, na.value = "white") + # , limits = c(0, 10)
  geom_label() +
  labs(x = "", y = "", fill = "Jaccard Distances", title = "Species (Jaccard Distances)") +
  my_theme +
  theme(
    text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom",
    legend.key.width = unit(dev.size()[1] / 8, "inches")
  )
print(OV_taxon)

### By Individuals in Species -----
vols_ls$Individuals$vols <- vols_ls$Individuals$vols[which(!unlist(lapply(vols_ls$Individuals$hv, is.na)))]
vols_ls$Individuals$hv <- vols_ls$Individuals$hv[which(!unlist(lapply(vols_ls$Individuals$hv, is.na)))]

spec_vec <- levels(plot_df$taxon)
cl <- parallel::detectCores() / 2
cl <- parallel::makeCluster(cl) # for parallel pbapply functions
parallel::clusterExport(cl,
  varlist = c(
    "vols_ls",
    "install.load.package",
    "package_vec1",
    "FUN.Overlap"
  ),
  envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec1, install.load.package))

Overlaps <- lapply(spec_vec, FUN = function(sp) {
  print(sp)

  IterPos <- which(lapply(
    strsplit(names(vols_ls$Individuals$hv), split = "_"),
    "[[", 4
  ) == sp)

  IterCombs <- as.data.frame(do.call(
    rbind,
    combn(
      x = names(vols_ls$Individuals$vols[IterPos]),
      m = 2, simplify = FALSE
    )
  ))
  colnames(IterCombs) <- c("ID1", "ID2")

  IterOver <- pbapply(IterCombs,
    cl = cl,
    MARGIN = 1, FUN = function(iter) {
      # print(iter)
      # sink("aux")
      over <- FUN.Overlap(
        data = vols_ls$Individuals$hv, # [which(names(vols_ls$Individuals$hv) %in% iter)]
        names = unlist(iter)
      )
      # sink(NULL)
      over
    }
  )
})
stopCluster(cl)
closeAllConnections()
names(Overlaps) <- spec_vec

OverID <- do.call(rbind, lapply(Overlaps, FUN = function(x) {
  data.frame(
    mean = mean(x),
    sd = sd(x)
  )
}))
OverID$SP <- spec_vec

OverID <- data.frame(
  Value = unlist(Overlaps),
  SP = rep(names(Overlaps), unlist(lapply(Overlaps, length)))
)
OverID$SP <- factor(OverID$SP, levels = unique(OverID$SP))

taxon.comps <- list(
  c("G. glomerata", "V. floribundum"),
  c("H. umbellata", "L. orbiculata"),
  c("P. bonplandianum", "R. macrochaeta")
)

OV_ID <- ggplot(OverID, aes(y = Value * 100, x = SP, fill = factor(SP))) +
  geom_boxplot() +
  stat_compare_means(comparisons = taxon.comps, method = "t.test", label = "p.signif") +
  scale_fill_manual(values = as.character(pal_lm)) +
  guides(fill = "none") +
  labs(x = "Species Identity", y = "Relatrive Overlap [%]", title = "Individual Hypervolumes (Jaccard Distances)") +
  theme_bw()
print(OV_ID)

### Manuscript Plot ----
overlaps_gg <- plot_grid(
  plot_grid(OV_fg + theme(legend.position = "none"), OV_taxon + theme(legend.position = "none"), labels = "AUTO"),
  get_legend(OV_taxon),
  OV_venn,
  OV_ID,
  ncol = 1, rel_heights = c(1.5, 0.1, 1, 1),
  labels = c("", "", "C", "D")
)
print(overlaps_gg)
ggsave(file.path("plots", "HV_Overlaps.png"), overlaps_gg, units = "in", height = 24 / 1.5, width = 21 / 1.5, dpi = 600)
ggsave(file.path("plots", "HV_Overlaps.pdf"), overlaps_gg, units = "in", height = 24 / 1.5, width = 21 / 1.5, dpi = 600)

## Hypervolumes by Elevation ----------------------------------------------
message("#### Hypervolumes by elevation")
elevs_vec <- unique(traits_df$elevation)

### calculate hypervolumes for each elevation separately for taxonomic and functional groupings
elev_hvs <- lapply(elevs_vec, FUN = function(elev) {
  message(elev)
  taxon_hv <- FUN.Hypervolumes(
    data = traits_df %>% filter(elevation == elev),
    Grouping = "taxon",
    TraitCols = c(6, 7:12)
  )
  fg_hv <- FUN.Hypervolumes(
    data = traits_df %>% filter(elevation == elev),
    Grouping = "functional_group",
    TraitCols = c(6, 7:12)
  )
  list(
    taxon = taxon_hv,
    fg = fg_hv
  )
})
names(elev_hvs) <- as.character(elevs_vec)

### calculate overlap between each subsequent elevation for taxonomic and functional groupings
elevs_vec <- sort(elevs_vec)
comps <- data.frame(
  first = c(elevs_vec[1], elevs_vec[2], elevs_vec[1]),
  second = c(elevs_vec[2], elevs_vec[3], elevs_vec[3])
)

OverElev_ls <- lapply(1:nrow(comps), FUN = function(i) {
  # i <- 1
  print(comps[i, ])

  ### functional group overlap
  fg_ov <- pblapply(unique(traits_df$functional_group), FUN = function(fg) {
    # fg <- unique(traits_df$functional_group)[3]

    comp_ls <- list(
      elev_hvs[[as.character(comps$first[i])]]$fg[[which(names(elev_hvs[[as.character(comps$first[i])]]$fg) == fg)]],
      elev_hvs[[as.character(comps$second[i])]]$fg[[which(names(elev_hvs[[as.character(comps$second[i])]]$fg) == fg)]]
    )
    names(comp_ls) <- c(as.character(comps$first[i]), as.character(comps$second[i]))

    ov_hv <- FUN.Overlap(
      data = comp_ls,
      names = names(comp_ls),
      what = "absolute"
    )

    data.frame(
      Elevation1 = as.character(comps$first[i]),
      Elevation2 = as.character(comps$second[i]),
      Size1 = get_volume(comp_ls[[1]]),
      Size2 = get_volume(comp_ls[[2]]),
      Functional_Group = fg,
      Overlap = ov_hv
    )
  })
  fg_ov_df <- do.call(rbind, fg_ov)

  ### taxonomic group overlap
  tx_ov <- pblapply(unique(traits_df$taxon), FUN = function(tx) {
    # tx <- unique(traits_df$taxon)[1]

    comp_ls <- list(
      elev_hvs[[as.character(comps$first[i])]]$taxon[[which(names(elev_hvs[[as.character(comps$first[i])]]$taxon) == tx)]],
      elev_hvs[[as.character(comps$second[i])]]$taxon[[which(names(elev_hvs[[as.character(comps$second[i])]]$taxon) == tx)]]
    )
    names(comp_ls) <- c(as.character(comps$first[i]), as.character(comps$second[i]))

    ov_hv <- FUN.Overlap(
      data = comp_ls,
      names = names(comp_ls),
      what = "absolute"
    )

    data.frame(
      Elevation1 = as.character(comps$first[i]),
      Elevation2 = as.character(comps$second[i]),
      Size1 = get_volume(comp_ls[[1]]),
      Size2 = get_volume(comp_ls[[2]]),
      Taxon = tx,
      Overlap = ov_hv
    )
  })
  tx_ov_df <- do.call(rbind, tx_ov)

  ### reporting back up
  list(
    taxon = tx_ov_df,
    fg = fg_ov_df
  )
})

fg_df <- do.call(rbind, lapply(OverElev_ls, function(x) x$fg))
fg_df$relativeOV <- fg_df$Overlap / (fg_df$Size1 + fg_df$Size2) * 100
fg_df$ElevDiff <- as.numeric(as.character(fg_df$Elevation2)) - as.numeric(as.character(fg_df$Elevation1))

fg_plot <- plot_grid(
  ggplot(fg_df, aes(x = Elevation1, y = Elevation2, fill = relativeOV)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, option = "E", begin = 0.2, limits = c(0, NA)) +
    facet_wrap(~Functional_Group) +
    geom_label(aes(label = paste(round(relativeOV, 2), "%")), fill = "white") +
    labs(x = "Elevation of Site 1 [m]", y = "Elevation of Site 2 [m]", fill = "Relative Overlap [%]") +
    my_theme +
    theme(
      text = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "top",
      legend.key.width = unit(dev.size()[1] / 18, "inches")
    ),
  ggplot(fg_df, aes(x = ElevDiff, y = relativeOV, fill = Functional_Group, col = Functional_Group)) +
    geom_point() +
    geom_line(linewidth = 2) +
    labs(x = "Elevation Difference [m]", y = "Relative Overlap [%]", fill = "Functional Group", color = "Functional Group") +
    scale_fill_manual(values = FG_pal) +
    scale_color_manual(values = FG_pal) +
    my_theme +
    theme(
      text = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.width = unit(dev.size()[1] / 15, "inches")
    ),
  ncol = 1
)
ggsave(file.path("plots", "HV_Over_Elevation_FG.png"), fg_plot, units = "in", height = 12, width = 12, dpi = 600)
ggsave(file.path("plots", "HV_Over_Elevation_FG.pdf"), fg_plot, units = "in", height = 12, width = 12, dpi = 600)



tx_df <- do.call(rbind, lapply(OverElev_ls, function(x) x$taxon))
tx_df$relativeOV <- tx_df$Overlap / (tx_df$Size1 + tx_df$Size2) * 100
tx_df$ElevDiff <- as.numeric(as.character(tx_df$Elevation2)) - as.numeric(as.character(tx_df$Elevation1))
tx_df$Taxon <- paste0(
  substr(unlist(lapply(strsplit(as.character(tx_df$Taxon), split = " "), "[[", 1)), 1, 1),
  ". ",
  unlist(lapply(strsplit(as.character(tx_df$Taxon), split = " "), "[[", 2))
)

tx_levels <- plot_df %>%
  dplyr::distinct(taxon, functional_group) %>%
  dplyr::mutate(
    taxon = as.character(taxon),
    functional_group = as.character(functional_group)
  ) %>%
  dplyr::arrange(
    factor(functional_group, levels = c("Forb", "Graminoid", "Woody")),
    taxon
  ) %>%
  dplyr::pull(taxon)
tx_df$Taxon <- factor(as.character(tx_df$Taxon), levels = tx_levels)

tx_plot <- plot_grid(
  ggplot(tx_df, aes(x = Elevation1, y = Elevation2, fill = relativeOV)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, option = "E", begin = 0.2, limits = c(0, NA)) +
    facet_wrap(~Taxon) +
    geom_label(aes(label = paste(round(relativeOV, 2), "%")), fill = "white") +
    labs(x = "Elevation of Site 1 [m]", y = "Elevation of Site 2 [m]", fill = "Relative Overlap [%]") +
    my_theme +
    theme(
      text = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "top",
      legend.key.width = unit(dev.size()[1] / 18, "inches")
    ),
  ggplot(tx_df, aes(x = ElevDiff, y = relativeOV, fill = Taxon, col = Taxon)) +
    geom_point() +
    geom_line(linewidth = 2) +
    labs(x = "Elevation Difference [m]", y = "Relative Overlap [%]", fill = "Species", color = "Species") +
    scale_fill_manual(values = pal_lm) +
    scale_color_manual(values = pal_lm) +
    my_theme +
    theme(
      text = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.width = unit(dev.size()[1] / 15, "inches")
    ),
  ncol = 1
)

ggsave(file.path("plots", "HV_Over_Elevation_TX.png"), tx_plot, units = "in", height = 12, width = 12, dpi = 600)
ggsave(file.path("plots", "HV_Over_Elevation_TX.pdf"), tx_plot, units = "in", height = 12, width = 12, dpi = 600)
