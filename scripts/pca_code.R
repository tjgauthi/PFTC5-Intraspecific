#### Data analysis - Ordination ####

## Authors: 

#### Load libraries ####
library(dplyr)
library(ggplot2)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(readr)
library(factoextra)

#* TO DO: 
#* PCA - use traits as  explanatory variables
#* How traits values vary between different sites?


# upload data 
data <- read_csv("data/raw/PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")

# select 2020, sites (WAY, ACJ and TRE), treatment and triats we are interest in 
data <- data %>%
  filter(site %in% c("WAY", "ACJ", "TRE") &
           year == 2020 & 
           treatment == "C" &
           trait %in% c("plant_height_cm", "leaf_area_cm2", "sla_cm2_g", "ldmc", "dry_mass_g", "leaf_thickness_mm"))

# should we log transform some of the traits?

data$site <- factor(data$site)
# data$site <- as.numeric(data$site)
data$taxon <- factor(data$taxon)
# data$taxon <- as.numeric(data$taxon)
data$functional_group <- factor(data$functional_group)
# data$functional_group <- as.numeric(data$functional_group)
# data <- data[!is.na(data$individual_nr),]

unique(data$site)
# unique(data$individual_nr)
unique(data$taxon)
unique(data$functional_group)
unique(data$trait)


# PCA - plot 1 
df1 <- data %>% 
  dplyr::select(id:value) %>% 
  filter(!id %in% c("CXX4125", "BDN3235")) %>% 
  tidyr::pivot_wider(names_from = trait, values_from = value) %>% 
  tibble::column_to_rownames("id") 


# df1.trait <- factor(data$trait)
pca_out <- prcomp(na.omit(df1[, -(1:3)]), center = TRUE, scale = TRUE)

# pca_out <- prcomp(df1, center = TRUE, scale. = TRUE)
summary(pca_out)
str(pca_out)
pca_out$rotation #look at laodings 

ggbiplot(pca_out, choices = c(1, 2)) + theme_classic()
# ggbiplot(pca_out, choices = c(3, 4)) + theme_classic()


ggbiplot(
  pca_out,
  scale = 0,
  choices = c(1, 2),
  ellipse = TRUE,
  # groups = df1.trait
) + theme_classic()


# fviz_pca_ind(pca_out, repel = TRUE,
#              col.var = "#2A2A2A", # Variables color
#              label = "none",
#              labelsize = 5, 
#              habillage = df1$functional_group,
#              addEllipses = TRUE,
#              ellipse.level = 0.95) +
#   theme_minimal(base_size = 15)
# 
# 
# 
# # PCA - plot 2
# df2 <- data %>% select(c(site, taxon, individual_nr, functional_group, value))
# df2.trait <- factor(data$trait)
# pca_out1 <- prcomp(df2, center = TRUE, scale. = TRUE)
# summary(pca_out1)
# 
# str(pca_out1)
# 
# ggbiplot(pca_out1,choices=c(1,2))+ theme_classic()
# 
# 
# ggbiplot(
#   pca_out1,
#   scale = 0,
#   choices = c(1, 2),
#   # ellipse = TRUE,
#   groups = df1.trait
# ) + theme_classic()

