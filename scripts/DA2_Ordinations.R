#### Data analysis - Ordination ####
## Authors: Fernanda, Elisa, Korina
# Others: Augustina, Fiorella

#### Load libraries ####
library(dplyr)
library(ggplot2)
library(devtools)
library(readr)
library(factoextra)
library(vegan)
# install_github("vqv/ggbiplot")
# library(ggbiplot)


# Data from data_import code
#source


# 
# # select 2020, sites (WAY, ACJ and TRE), treatment and traits we are interest in
# data <- data %>%
#   filter(
#     site %in% c("WAY", "ACJ", "TRE") &
#       year == 2020 &
#       treatment == "C" &
#       trait %in% c(
#         "plant_height_cm",
#         "leaf_area_cm2",
#         "sla_cm2_g",
#         "ldmc",
#         "dry_mass_g",
#         "leaf_thickness_mm"
#       ) &
#       taxon %in% c(
#         "Halenia umbellata",
#         "Lachemilla orbiculata",
#         "Paspalum bonplandianum",
#         "Rhynchospora macrochaeta",
#         "Vaccinium floribundum",
#         "Gaultheria glomerata"
#       )
#   )
# 
# 
# data$site <- factor(data$site)
# data$taxon <- factor(data$taxon)
# data$functional_group <- factor(data$functional_group)
# # data$functional_group <- as.numeric(data$functional_group)
# # data <- data[!is.na(data$individual_nr),]
# 
# unique(data$site)
# # unique(data$individual_nr)
# unique(data$taxon)
# unique(data$functional_group)
# unique(data$trait)

ord_traits <- traits_wide |> 
  select(site, taxon, leaf_uid, plant_height_cm, dry_mass_g, leaf_area_cm2, sla_cm2_g, ldmc, leaf_thickness_mm)


### PCA analysis -------------------------------------------------------------------

# log transforme variables
ord_traits <- ord_traits |> 
  mutate(plant_height_cm = log(plant_height_cm),
         dry_mass_g = log(dry_mass_g), 
         leaf_area_cm2 = log(leaf_area_cm2), 
         sla_cm2_g = log(sla_cm2_g), 
         ldmc = log(ldmc), 
         leaf_thickness_mm = log(leaf_thickness_mm))

#Checks for NAs
ord_traits |> 
  summarise(across(
    c(plant_height_cm, dry_mass_g, leaf_area_cm2, 
      sla_cm2_g, ldmc, leaf_thickness_mm),
    ~ sum(is.na(.))
  ))

#NO NAs, so moving forward

#df1$site <- factor(df1$site, levels = c("WAY", "ACJ", "TRE"))


# create a new DF for pca analisys
pca_out <- prcomp(ord_traits[, c(4:9)], center = TRUE, scale = TRUE)


# pca_out <- prcomp(df1, center = TRUE, scale. = TRUE)
summary(pca_out)
str(pca_out)
pca_out$rotation # look at laodings


# Plot PCA with ggplot2
scores <- as.data.frame(pca_out$x) ## getting the scores
scores.1 <- cbind(scores, ord_traits)
pca.loadings <- data.frame(Variables = rownames(pca_out$rotation), pca_out$rotation) # drawing the arrows


# Combine SITES and TAXON in the same graph
c <-
  ggplot(scores.1, aes(x = PC1, y = PC2, color=taxon, shape=site)) +  
geom_point(stat="identity", position=position_dodge(), size=3, alpha = 0.8)+ 
  scale_fill_hue(l=40) + 
  coord_fixed(ratio = 1, xlim = range(scores$PC1), ylim = range(scores$PC2))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  xlab(paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")) +
  scale_color_manual(
    breaks = c(
      "Halenia umbellata",
      "Lachemilla orbiculata",
      "Paspalum bonplandianum",
      "Rhynchospora macrochaeta",
      "Vaccinium floribundum",
      "Gaultheria glomerata"
    ),
    values = c(
      "#016392",
      "#A0CBE8",
      "#E19825",
      "#F7C480",
      "#9FCD99",
      "#3E8853"
    )
  ) +
  geom_segment(
    data = pca.loadings,
    aes(
      x = 0,
      y = 0,
      xend = (PC1 * 4),
      yend = (PC2 * 2)
    ),
    arrow = arrow(length = unit(1 / 2, "picas")),
    inherit.aes = FALSE,
    color = "black"
  ) +
  annotate(
    "text",
    x = (pca.loadings$PC1 * 4),
    y = (pca.loadings$PC2 * 3),
    label = c("H", "DM", "LA", "SLA", "LDMC", "LT"),
    size = 3.5
  ) 

c.final <-
  c + guides(
    color = guide_legend(
      title = "Plant species",
      order = 2,
      title.position = "top",
      legend.title.align = 0
    ),
    shape = guide_legend(
      "Sites",
      order = 1,
      title.position = "top",
      legend.title.align = 0
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12,  face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

# saving the plot
ggsave(
  "pca.total.png",
  width = 22,
  height = 21,
  units = "cm",
  dpi = 600,
  c_final
)


### (NOT USED) Plot taking TAXON into account
# p <- ggplot(data = scores.1, aes(x = PC1, y = PC2, color = taxon)) +
#   geom_point(size = 2) +
#   scale_fill_hue(l = 40) +
#   coord_fixed(
#     ratio = 1,
#     xlim = range(scores$PC1),
#     ylim = range(scores$PC2)
#   ) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   theme_classic() +
#   xlab("PC 1 (49.4%)") +
#   ylab("PC 2 (31.4%)") +
#   geom_segment(
#     data = pca.loadings,
#     aes(
#       x = 0,
#       y = 0,
#       xend = (PC1 * 3.5),
#       yend = (PC2 * 2)
#     ),
#     arrow = arrow(length = unit(1 / 2, "picas")),
#     color = "black"
#   ) +
#   geom_point(size = 3) +
#   annotate(
#     "text",
#     x = (pca.loadings$PC1 * 3.5),
#     y = (pca.loadings$PC2 * 2),
#     label = c("Height", "Dry mass", "leaf area", "SLA", "LDMC", "Leaf thickness")
#   ) +
#   theme(legend.position = "right") +
#   scale_color_brewer(palette = "Paired")
# p + guides(color = guide_legend(title = "Plant species"))


### (NOT USED) Plot taking SITES into account
# s <- ggplot(data = scores.1, aes(x = PC1, y = PC2, color = site)) +
#   geom_point(size = 2) +
#   scale_fill_hue(l = 40) +
#   coord_fixed(
#     ratio = 1,
#     xlim = range(scores$PC1),
#     ylim = range(scores$PC2)
#   ) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   theme_classic() +
#   xlab("PC 1 (49.4%)") +
#   ylab("PC 2 (31.4%)") +
#   geom_segment(
#     data = pca.loadings,
#     aes(
#       x = 0,
#       y = 0,
#       xend = (PC1 * 3.5),
#       yend = (PC2 * 2)
#     ),
#     arrow = arrow(length = unit(1 / 2, "picas")),
#     color = "black"
#   ) +
#   geom_point(size = 3) +
#   annotate(
#     "text",
#     x = (pca.loadings$PC1 * 3.5),
#     y = (pca.loadings$PC2 * 2),
#     label = c("Height", "Dry mass", "leaf area", "SLA", "LDMC", "Leaf thickness")
#   ) +
#   theme(legend.position = "right") +
#   scale_color_brewer(palette = "Paired")
# s + guides(color = guide_legend(title = "Sites"))


### RDA analysis ---------------------------------------------------------------------
# to see if there are significant differences in the factor

RDA_out <- rda(df1[, c(7:12)] ~ scale(elevation) * taxon, center = TRUE, scale = TRUE, data = df1)

summary(RDA_out)
print(RDA_out)

# We construct the models
RDA_1 <- rda(df1[, c(7:12)] ~ scale(elevation), center = TRUE, scale = TRUE, data = df1)
RDA_2 <- rda(df1[, c(7:12)] ~ taxon, center = TRUE, scale = TRUE, data = df1)
RDA_3 <- rda(df1[, c(7:12)] ~ scale(elevation) + taxon, center = TRUE, scale = TRUE, data = df1)
RDA_4 <- rda(df1[, c(7:12)] ~ scale(elevation) * taxon, center = TRUE, scale = TRUE, data = df1)

# then we compare the models
anova(RDA_1, RDA_3)
anova(RDA_1, RDA_4)
anova(RDA_3, RDA_4)

anova(RDA_2, RDA_3)
anova(RDA_2, RDA_4)
anova(RDA_3, RDA_4)

# when doing ANOVA the interaction between elevation and taxa
# there is an interaction between elevation and factor