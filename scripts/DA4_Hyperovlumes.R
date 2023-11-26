#' ####################################################################### #
#' PROJECT: [DA.4 - Hypervolumes] 
#' CONTENTS: 
#'  - Calculating hypervolumes and comparing their sizes
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
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
  "tidyverse"
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
FG_pal <- c(colorRampPalette(pal_lm[1:2])(3)[2], 
            colorRampPalette(pal_lm[3:4])(3)[2],
            colorRampPalette(pal_lm[5:6])(3)[2])
names(FG_pal) <- c("Forb", "Graminoid", "Woody")

# DATA LOADING ============================================================
source("scripts/0_data_import.R") # sourcing data import script
traits_wide$ID <- paste(traits_wide$site, traits_wide$plot_id, traits_wide$individual_nr, traits_wide$taxon, sep="_")

# POOL ELEVATION WITHIN SITE. WAY: 3101; ACJ: 3468; TRE: 3715
traits_df <- na.omit(traits_wide[ , c(-2:-5, -17:-19)])
traits_df$elevation <- ifelse(traits_df$site == 'ACJ', 3468,
                              ifelse(traits_df$site == 'WAY', 3101,
                                     ifelse(traits_df$site == 'TRE', 3715, NA)))
# CREATE DATA FRAME OF GROUP MEMBERSHIP OF EACH SPECIES
Grouping_df <- as.data.frame(table(traits_df[,c("taxon", "functional_group")])) # tally up each species and functional group combination
Grouping_df <- Grouping_df[Grouping_df$Freq != 0, -3] # remove every non-observed pairing
Grouping_df$taxon <- as.character(Grouping_df$taxon) # ensure that taxon column is character

# ANALYSES ================================================================

## Functionality ----------------------------------------------------------
## CALCULATE HYPERVOLUMES FOR DESIRED GROUPINGS OF DATA
FUN.Hypervolumes <- function(data = traits_df, # the data frame with traits and the grouping column
                             TraitCols = 8:12, # which columns of the data frame contain the trait values, default excludes wet mass and plant height due to identical plant height for many individuals
                             Grouping = "family" # by which column to establish groups
){
  data_scale <- log(data[,TraitCols]) # now applying logarithmic transformation
  data_ls <- split(data_scale, data[,Grouping]) # split the data by grouping argument
  hv_ls <- lapply(data_ls, FUN = function(x){
    if(0 %in% apply(apply(x, 2, range), 2, diff)){
      return(NA)
    }else{
      return(hypervolume(x))  
    }
  }
  ) # calculate hypervolume for each group separately
  return(hv_ls) # return list of volumes with one element for each group
  # volumes_ls <- lapply(hv_ls, get_volume)
  # return(volumes_ls)
}

## CALCULATE OVERLAP OF HYPERVOLUMES
FUN.Overlap <- function(data, # the list of hypervolumes from which to pull
                        names, # the names of the hypervolumes which to compare
                        what = "jaccard"
){
  overlap <- hypervolume_set(data[[names[1]]], data[[names[2]]], check.memory = FALSE) # extract hypervolumes and safe as set
  if(what == "jaccard"){
    overlap <- hypervolume::hypervolume_overlap_statistics(overlap)[1] # jaccard distance calculation
  }else{
    overlap <- as.numeric(get_volume(overlap[["Intersection"]]))
  }
  return(overlap) # return numeric overlap
}

# Function to get post-hoc test stats
get_posthoc_stats <- function(model, group_var) {
  
  cmp <- do.call(mcp, setNames(list("Tukey"), group_var))
  
  stats <- glht(model, linfct = cmp) %>%
    summary(test = adjusted("bonferroni"))
  
  stats.df <- cbind('Estimate' = stats$test$coefficients,
                    'Std_Error' = stats$test$sigma,
                    'z_value' = stats$test$tstat,
                    'p_value' = stats$test$pvalues
  ) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = 'Comparison')
  
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
if(!file.exists(file.path("./data", "Hypervolumes.RData"))){
  
  ### Individual Hypervolumes ####
  ID_hv <- FUN.Hypervolumes(Grouping = "ID")
  ID_vols <- lapply(ID_hv, get_volume)
  
  ### Species Hypervolumes ####
  taxon_hv <- FUN.Hypervolumes(Grouping = "taxon", TraitCols = c(6, 8:12)) # also including plant height now
  taxon_vols <- lapply(taxon_hv, get_volume)
  
  ### Functional Group Hypervolumes ####
  functional_group_hv <- FUN.Hypervolumes(Grouping = "functional_group", TraitCols = c(6, 8:12)) # also including plant height now
  functional_group_vols <- lapply(functional_group_hv, get_volume)
  
  ### Saving Hypervolume Lists ####
  save(ID_hv, ID_vols, 
       taxon_hv, taxon_vols, 
       functional_group_hv, functional_group_vols,
       file = file.path("./data", "Hypervolumes.RData"))
}else{
  load(file.path("./data", "Hypervolumes.RData"))
}

## Comparison of Volumes --------------------------------------------------
message("Hypervolume comparison")
### Preparing data
Vols_df <- data.frame(values = unlist(ID_vols),
                      ID = names(unlist(ID_vols))
)
Vols_df$ID <- gsub(Vols_df$ID, pattern = ".untitled", replacement = "")
# create plotting data
plot_df <- plyr::join(x = Vols_df, y = traits_df, by = "ID") %>%
  distinct(values, ID, .keep_all = TRUE) 
# %>% 
#   filter(values != max(.$values)) #  and filter out the outlier which is TRE_1_NA_Rhynchospora macrochaeta
### factorise columns that need to be factors
plot_df$taxon <- factor(plot_df$taxon, levels = c('Halenia umbellata', 'Lachemilla orbiculata', 'Paspalum bonplandianum', 'Rhynchospora macrochaeta', 'Gaultheria glomerata', 'Vaccinium floribundum'))
plot_df$functional_group <- as.factor(plot_df$functional_group)
plot_df$site %<>% factor(levels = c('WAY', 'ACJ', 'TRE'))

### Volume comparison by taxon -----
print("Volume by taxon")
### set comparisons for plotted boxplots
taxon.comps <- list(c("Gaultheria glomerata", "Halenia umbellata"), c("Gaultheria glomerata", "Lachemilla orbiculata"), c("Gaultheria glomerata", "Paspalum bonplandianum"), c("Gaultheria glomerata", "Rhynchospora macrochaeta"), c("Gaultheria glomerata", "Vaccinium floribundum"),
                    c("Halenia umbellata", "Lachemilla orbiculata"), c("Halenia umbellata", "Paspalum bonplandianum"), c("Halenia umbellata", "Rhynchospora macrochaeta"), c("Halenia umbellata", "Vaccinium floribundum"),
                    c("Lachemilla orbiculata", "Paspalum bonplandianum"), c("Lachemilla orbiculata", "Rhynchospora macrochaeta"), c("Lachemilla orbiculata", "Vaccinium floribundum"),
                    c("Paspalum bonplandianum", "Rhynchospora macrochaeta"), c("Paspalum bonplandianum", "Vaccinium floribundum"),
                    c("Rhynchospora macrochaeta", "Vaccinium floribundum")
)

taxon.boxplot1 <- ggplot(plot_df, aes(y = values, x = taxon, fill = taxon)) +
  geom_boxplot() + 
  scale_fill_manual(values = pal_lm) +
  # stat_compare_means(comparisons = taxon.comps, method = 't.test', label = 'p.signif') +
  my_theme +
  ylab('Hypervolume Size') + xlab('Taxon') + labs(fill = 'Species') +
  facet_wrap(~site) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12) 
  )

print(taxon.boxplot1)
ggsave('taxon.boxplot1.png', taxon.boxplot1, units = 'in', height = 5.5, width = 11.4, dpi = 600)
ggsave('taxon.boxplot1.pdf', taxon.boxplot1, units = 'in', height = 5.5, width = 11.4, dpi = 600)

### Volume comparison by functional groups -----
print("Volume by functional group")
fg.comps <- list(c("Forb", "Graminoid"), c("Graminoid", "Woody"), c("Forb", "Woody"))
taxon.boxplot2 <- ggplot(plot_df, aes(y = values, x = functional_group, fill = taxon)) +
  geom_boxplot() + 
  scale_fill_manual(values = pal_lm) +
  # stat_compare_means(comparisons = fg.comps, method = 't.test', label = 'p.signif', size = 5) +
  # stat_compare_means(aes(group = taxon), label = 'p.signif', label.y = c(5, 9, 7.5), size = 5) +
  my_theme +
  ylab('Hypervolume Size') + xlab('Functional Group') + labs(fill = 'Species') +
  # ylim(0, 15) +
  facet_wrap(~site) +
  theme(text = element_text(size = 16),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12) 
  )

print(taxon.boxplot2)
ggsave('taxon.boxplot2.png', taxon.boxplot2, units = 'in', height = 5.5, width = 11.4, dpi = 600)
ggsave('taxon.boxplot2.pdf', taxon.boxplot2, units = 'in', height = 5.5, width = 11.4, dpi = 600)


### Volume comparison by elevation -----
print("Volume by elevation")
taxon.scatterplot <- ggplot(plot_df, aes(y = values, x = elevation, color = taxon)) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x + I(x^2), se = F) +
  scale_color_manual(values = pal_lm) +
  my_theme +
  ylab('Hypervolume Size') + xlab('Elevation (mASL)') + labs(color = 'Species') +
  theme(text = element_text(size = 16),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12) 
  )
print(taxon.scatterplot)
ggsave('taxon.scatterplot.png', taxon.scatterplot, units = 'in', height = 4.7, width = 7.4, dpi = 600)
ggsave('taxon.scatterplot.pdf', taxon.scatterplot, units = 'in', height = 4.7, width = 7.4, dpi = 600)

## Overlap of Volumes -----------------------------------------------------
print("Hypervolume overlap")
vols_ls <- list(
  Individuals = list(vols = ID_vols,
                     hv = ID_hv),
  Taxon = list(vols = taxon_vols,
               hv = taxon_hv),
  Functional = list(vols = functional_group_vols,
                    hv = functional_group_hv)
)

# plot_df <- data.frame(sp = c(Grouping_df$taxon[Grouping_df$functional_group == "Forb"], "Overlap",
#                              Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"], "Overlap",
#                              Grouping_df$taxon[Grouping_df$functional_group == "Woody"], "Overlap"),
#                       fg = rep(c("Forb", "Graminoid", "Woody"), each = 3),
#                       value = as.numeric(
#                         c(taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Forb"]],
#                           FUN.Overlap(data = taxon_hv,
#                                       names = Grouping_df$taxon[Grouping_df$functional_group == "Forb"]),
#                           taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]],
#                           FUN.Overlap(data = taxon_hv,
#                                       names = Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]),
#                           taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Woody"]],
#                           FUN.Overlap(data = taxon_hv,
#                                       names = Grouping_df$taxon[Grouping_df$functional_group == "Woody"]))
#                       ),
#                       sp2 = rep(c("sp1", "sp2", "Overlap"), 3)
# )
# plot_df$value[3] <- plot_df$value[3]*sum(plot_df$value[1:2])
# plot_df$value[6] <- plot_df$value[6]*sum(plot_df$value[4:5])
# plot_df$value[9] <- plot_df$value[9]*sum(plot_df$value[7:8])
# plot_df$sp <- c("H. umbellata", "L. orbiculata", "Overlap", 
#                 "P. bonplandianum", "R. macrochaeta", "Overlap",
#                 "G. glomerata", "V. floribundum", "Overlap"
# )
# plot_df$sp <- factor(plot_df$sp, levels=c("H. umbellata", "P. bonplandianum", "G. glomerata", 
#                                           "Overlap", 
#                                           "L. orbiculata", "R. macrochaeta", "V. floribundum"
# ))
# 
# 
# 
# c(pal_lm[c(1,3,5)], "#808080", pal_lm[c(2,4,6)])
# 
# ### plot it out
# overlap.plot <- ggplot(plot_df, aes(x = fg, y = value, fill = sp, label = sp)) + 
#   geom_bar(position="stack", stat="identity") + 
#   geom_text(size = 5, position = position_stack(vjust = 0.6), family = 'Helvetica', fontface = rep(c('italic', 'italic', 'plain'), 3)) + 
#   scale_fill_manual(values = c(pal_lm[c(1,3,5)], "#808080", pal_lm[c(2,4,6)])) +
#   my_theme +
#   labs(x = "Functional Group", y = "Hypervolume Size") +
#   theme(legend.position = "none",
#         text = element_text(size = 16)
#   ) 
# 
# ## original volume sizes
# functional_group_vols$Forb
# functional_group_vols$Graminoid
# functional_group_vols$Woody
# 
# ## overlaps
# ForbGram_ov <- FUN.Overlap(functional_group_hv, c("Forb", "Graminoid"))
# ForbWood_ov <- FUN.Overlap(functional_group_hv, c("Forb", "Woody"))
# WoodGram_ov <- FUN.Overlap(functional_group_hv, c("Woody", "Graminoid"))
# 
# ## report overlaps
# message("Forb + Graminoid")
# ForbGram_ov
# message("Forb + Woody")
# ForbWood_ov
# message("Graminoid + Woody")
# WoodGram_ov

### By Functional Group -----
Overlaps <- as.data.frame(do.call(rbind, combn(x = names(vols_ls[["Functional"]]$vols), m = 2, simplify = FALSE)))
colnames(Overlaps) <- c("G1", "G2")
Overlap <- apply(Overlaps, 1, FUN = function(y){
  FUN.Overlap(data = vols_ls[["Functional"]]$hv,
              names = y)
})
Overlaps$Value <- Overlap

Volumes <- data.frame(
  Groups = names(vols_ls[["Functional"]]$vols),
  Value = unlist(vols_ls[["Functional"]]$vols)
)

OV_Functional <- ggplot(Overlaps, aes(x = G1, y = G2, fill = Value, label = round(Value,2))) + 
  geom_tile() + 
  scale_fill_viridis_c(direction = -1, option = "E", begin = 0.2) +
  geom_label() + 
  labs(x = "", y = "", fill = "Overlap") +
  my_theme
VO_Functional <- ggplot(Volumes, aes(x = Groups, y = Value, fill = factor(Groups))) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = FG_pal) +
  labs(x = "", y = "") + 
  guides(fill = FALSE) + 
  theme_bw()
VO_FunctionalID <- ggplot(data = plot_df, aes(x = functional_group, y = values, fill = factor(functional_group))) +
  geom_violin() + 
  scale_fill_manual(values = FG_pal) +
  labs(x = "", y = "") + 
  guides(fill = FALSE) + 
  my_theme

FG_plot <- plot_grid(OV_Functional, 
          grid.arrange(arrangeGrob(plot_grid(VO_Functional, VO_FunctionalID, ncol = 2), 
                                   left = "Volume", bottom = "Functional Group")), 
          ncol = 2, rel_widths = c(1, 1.5))

### By Taxon -----
Overlaps <- data.frame(G1 = c("H. umbellata", "P. bonplandianum", "G. glomerata"),
                       G2 = c("L. orbiculata", "R. macrochaeta", "V. floribundum"),
                       Value = c(FUN.Overlap(data = taxon_hv,
                                             names = Grouping_df$taxon[Grouping_df$functional_group == "Forb"],
                                             what = "absolute"),
                                 FUN.Overlap(data = taxon_hv,
                                             names = Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"],
                                             what = "absolute"),
                                 FUN.Overlap(data = taxon_hv,
                                             names = Grouping_df$taxon[Grouping_df$functional_group == "Woody"],
                                             what = "absolute")
                                 ))

Volumes <- data.frame(
  Groups = names(vols_ls[["Taxon"]]$vols),
  Value = unlist(vols_ls[["Taxon"]]$vols),
  FG = c("Woody", "Forb", "Forb", "Graminoid", "Graminoid", "Woody")
)

Forbs_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`Halenia umbellata`), 2), 
  round(as.numeric(vols_ls$Taxon$vols$`Lachemilla orbiculata`), 2), 
  round(Overlaps$Value[1], 4),
  c("H. umbellata", "L. orbiculata"),
  fill = pal_lm[1:2],
  alpha = 1,
  ind = FALSE,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0,
  cat.dist = 0.04,
  ext.dist = -0.15,
  ext.length = 0.85
)
Forbs_euler <- ggdraw(Forbs_euler)

Gram_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`Paspalum bonplandianum`), 2), 
  round(as.numeric(vols_ls$Taxon$vols$`Rhynchospora macrochaeta`), 2), 
  round(Overlaps$Value[2], 4),
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
Gram_euler <- ggdraw(Gram_euler)

Wood_euler <- draw.pairwise.venn(
  round(as.numeric(vols_ls$Taxon$vols$`Gaultheria glomerata`), 2), 
  round(as.numeric(vols_ls$Taxon$vols$`Vaccinium floribundum`), 2), 
  round(Overlaps$Value[3], 4),
  c("G. glomerata", "V. floribundum"),
  fill = pal_lm[5:6],
  alpha = 1,
  ind = FALSE,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0,
  cat.dist = 0.04,
  ext.dist = -0.15,
  ext.length = 0.85
)
Wood_euler <- ggdraw(Wood_euler)

Taxon_plot <- plot_grid(Forbs_euler, Gram_euler, Wood_euler, ncol = 3)

### By Individuals in Species -----
vols_ls$Individuals$vols <- vols_ls$Individuals$vols[which(!unlist(lapply(vols_ls$Individuals$hv, is.na)))]
vols_ls$Individuals$hv <- vols_ls$Individuals$hv[which(!unlist(lapply(vols_ls$Individuals$hv, is.na)))]

spec_vec <- levels(plot_df$taxon)
cl <- parallel::detectCores()
cl <- parallel::makeCluster(cl) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c("vols_ls", 
                                    "install.load.package",
                                    "package_vec1",
                                    "FUN.Overlap"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec1, install.load.package))

Overlaps <- lapply(spec_vec, FUN = function(sp){
  print(sp)
  
  IterPos <- which(lapply(
      strsplit(names(vols_ls$Individuals$hv), split = "_"), 
      "[[", 4) == sp)
  
  IterCombs <-  as.data.frame(do.call(rbind, 
                            combn(x = names(vols_ls$Individuals$vols[IterPos]), 
                                  m = 2, simplify = FALSE)
  ))
  colnames(IterCombs) <- c("ID1", "ID2")
  
  IterOver <- pbapply(IterCombs, 
                      cl = cl,
                      MARGIN = 1, FUN = function(iter){
    # print(iter)
    sink("aux")
    over <- FUN.Overlap(data = vols_ls$Individuals$hv, # [which(names(vols_ls$Individuals$hv) %in% iter)]
                names = unlist(iter))
    sink(NULL)
    over
  })
})
names(Overlaps) <- spec_vec

OverID <- do.call(rbind, lapply(Overlaps, FUN = function(x){
  data.frame(mean = mean(x),
  sd = sd(x)
  )
}))
OverID$SP <- spec_vec

OverID <- data.frame(Value = unlist(Overlaps),
                     SP = rep(names(Overlaps), unlist(lapply(Overlaps, length))))
OverID$SP <- factor(OverID$SP, levels = unique(OverID$SP))

ID_plot <- ggplot(OverID, aes(y = Value, x = SP, fill = factor(SP))) +
  geom_boxplot() +
  scale_fill_manual(values = as.character(pal_lm)) + 
  guides(fill = "none") + 
  labs(x = "Species Identity", y = "Jaccard Overlap Statistic") + 
  theme_bw()

### Fusing of Plots ----
overlap.plot <- plot_grid(FG_plot, Taxon_plot, ID_plot, nrow = 3, labels = "auto")
ggsave('overlap.plot.png', overlap.plot, units = 'in', height = 10, width = 10, dpi = 600)
ggsave('overlap.plot.pdf', overlap.plot, units = 'in', height = 10, width = 10, dpi = 600)

## Linear Mixed Effect Model of Volume Size -------------------------------
print("Hypervolume linear mixed effect model")
### taxonomic grouping of species
fg <- traits_wide %>%
  mutate(taxon = str_replace(taxon, ' ', '.')) %>%
  distinct(taxon, functional_group, family)
### elevation at sites		
el.values <- traits_wide %>%
  group_by(site) %>%
  summarise(elevation = mean(elevation))
### model data frame containing volume size, study compartments (site, plot, individual), taxonomic memembership, and elevation 				
temp <- Vols_df %>%
  separate(ID, c('site', 'plot_id', 'individual_nr', 'taxon'),  sep = '_') %>%
  mutate(taxon = str_replace(taxon, ' ', '.'),
         plot_id = factor(plot_id),
         individual_nr = as.numeric(individual_nr)
  ) %>%
  left_join(fg, by = 'taxon') %>%
  left_join(el.values, by = 'site') %>%
  na.omit # THIS WILL REMOVE THE OUTLIERS
temp$site %<>% factor(levels = c('WAY', 'ACJ', 'TRE'))
temp$individual_uid <- paste(temp$site, temp$plot_id, temp$individual_nr, temp$taxon, sep = "_")

# Add elevation instead of site
# THESE ARE THE FINAL MODELS as of November 25, 2023
print("taxon model ----")
hv.taxon.model <- lmerTest::lmer(
  values ~ 0 + scale(elevation)*taxon + (1|site/plot_id), 
  data = temp)
print(car::Anova(hv.taxon.model, type = 3))

taxon.comparisons <- emmeans(hv.taxon.model, list(pairwise ~ taxon), adjust = "tukey")
print(taxon.comparisons)

print("functional group model ----")
hv.fg.model <- lmerTest::lmer(
  values ~ 0 + scale(elevation)*functional_group + (1|site/plot_id), 
  data = temp)
print(car::Anova(hv.fg.model, type = 3))

print("functional group comparison ----")
fg.comparisons <- emmeans(hv.fg.model, list(pairwise ~ functional_group), adjust = "tukey")
print(fg.comparisons)
