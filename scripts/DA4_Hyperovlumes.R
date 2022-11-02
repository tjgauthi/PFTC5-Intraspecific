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
package_vec <- c(
  "hypervolume",
  "car",
  "multcomp",
  "ggpubr",
  "emmeans",
  "lmerTest",
  "magrittr"
)
sapply(package_vec, install.load.package)

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
# scales::show_col(pal_lm) 

# DATA LOADING ============================================================
source("scripts/0_data_import.R") # sourcing data import script
traits_wide$ID <- paste(traits_wide$individual_uid, traits_wide$taxon, sep="_")

# POOL ELEVATION WITHIN SITE. WAY: 3101; ACJ: 3468; TRE: 3715
traits_df <- na.omit(traits_wide[ , c(-2:-4, -16:-18)])
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
                             TraitCols = c(6, 8:12), # which columns of the data frame contain the trait values, default excludes wet mass
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
                        names # the names of the hypervolumes which to compare
){
  overlap <- hypervolume_set(data[[names[1]]], data[[names[2]]], check.memory = FALSE) # extract hypervolumes and safe as set
  overlap <- hypervolume::hypervolume_overlap_statistics(overlap)[1] # jaccard distance calculation
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
  taxon_hv <- FUN.Hypervolumes(Grouping = "taxon")
  taxon_vols <- lapply(taxon_hv, get_volume)
  
  ### Functional Group Hypervolumes ####
  functional_group_hv <- FUN.Hypervolumes(Grouping = "functional_group")
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
# create plotting data and filter out the outlier which is TRE_1_NA_Rhynchospora macrochaeta
plot_df <- plyr::join(x = Vols_df, y = traits_df, by = "ID") %>%
  distinct(values, ID, .keep_all = TRUE) %>%
  filter(values != max(.$values)) 
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
  stat_compare_means(comparisons = taxon.comps, method = 't.test', label = 'p.signif') +
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
  stat_compare_means(comparisons = fg.comps, method = 't.test', label = 'p.signif', size = 5) +
  stat_compare_means(aes(group = taxon), label = 'p.signif', label.y = c(5, 9, 7.5), size = 5) +
  my_theme +
  ylab('Hypervolume Size') + xlab('Functional Group') + labs(fill = 'Species') +
  ylim(0, 15) +
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
### create plotting object
plot_df <- data.frame(sp = c(Grouping_df$taxon[Grouping_df$functional_group == "Forb"], "Overlap",
                             Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"], "Overlap",
                             Grouping_df$taxon[Grouping_df$functional_group == "Woody"], "Overlap"),
                      fg = rep(c("Forb", "Graminoid", "Woody"), each = 3),
                      value = as.numeric(
                        c(taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Forb"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Forb"]),
                          taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Graminoid"]),
                          taxon_vols[Grouping_df$taxon[Grouping_df$functional_group == "Woody"]],
                          FUN.Overlap(data = taxon_hv,
                                      names = Grouping_df$taxon[Grouping_df$functional_group == "Woody"]))
                      ),
                      sp2 = rep(c("sp1", "sp2", "Overlap"), 3)
)
plot_df$value[3] <- plot_df$value[3]*sum(plot_df$value[1:2])
plot_df$value[6] <- plot_df$value[6]*sum(plot_df$value[4:5])
plot_df$value[9] <- plot_df$value[9]*sum(plot_df$value[7:8])
plot_df$sp <- c("H. umbellata", "L. orbiculata", "Overlap", 
                "P. bonplandianum", "R. macrochaeta", "Overlap",
                "G. glomerata", "V. floribundum", "Overlap"
)
plot_df$sp <- factor(plot_df$sp, levels=c("H. umbellata", "P. bonplandianum", "G. glomerata", 
                                          "Overlap", 
                                          "L. orbiculata", "R. macrochaeta", "V. floribundum"
))

### plot it out
overlap.plot <- ggplot(plot_df, aes(x = fg, y = value, fill = sp, label = sp)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(size = 5, position = position_stack(vjust = 0.6), family = 'Helvetica', fontface = rep(c('italic', 'italic', 'plain'), 3)) + 
  scale_fill_manual(values = c(pal_lm[c(1,3,5)], "#808080", pal_lm[c(2,4,6)])) +
  my_theme +
  labs(x = "Functional Group", y = "Hypervolume Size") +
  theme(legend.position = "none",
        text = element_text(size = 16)
  ) 
print(overlap.plot)
ggsave('overlap.plot.png', overlap.plot, units = 'in', height = 7, width = 7, dpi = 600)
ggsave('overlap.plot.pdf', overlap.plot, units = 'in', height = 7, width = 7, dpi = 600)

## original volume sizes
functional_group_vols$Forb
functional_group_vols$Graminoid
functional_group_vols$Woody

## overlaps
ForbGram_ov <- FUN.Overlap(functional_group_hv, c("Forb", "Graminoid"))
ForbWood_ov <- FUN.Overlap(functional_group_hv, c("Forb", "Woody"))
WoodGram_ov <- FUN.Overlap(functional_group_hv, c("Woody", "Graminoid"))

## report overlaps
message("Forb + Graminoid")
ForbGram_ov
message("Forb + Woody")
ForbWood_ov
message("Graminoid + Woody")
WoodGram_ov

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
  na.omit # THIS WILL REMOVE THE OUTLIER AT INDEX = 102
temp$site %<>% factor(levels = c('WAY', 'ACJ', 'TRE'))

# Add elevation instead of site
# THESE ARE THE FINAL MODELS as of June 16, 2022
print("taxon model ----")
hv.taxon.model <- lmerTest::lmer(values ~ 0 + scale(elevation)*taxon + (1|site/plot_id), data = temp)
print(car::Anova(hv.taxon.model, type = 3))
print("functional group model ----")
hv.fg.model <- lmerTest::lmer(values ~ 0 + scale(elevation)*functional_group + (1|site/plot_id), data = temp)
print(car::Anova(hv.fg.model, type = 3))
print("taxon comparison ----")
taxon.comparisons <- emmeans(hv.taxon.model, list(pairwise ~ taxon), adjust = "tukey")
print(taxon.comparisons)
fg.comparisons <- emmeans(hv.fg.model, list(pairwise ~ functional_group), adjust = "tukey")
print(fg.comparisons)
