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
#install_github("vqv/ggbiplot")
# library(ggbiplot)

# upload data 
data <- read_csv("data/raw/PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")

# select 2020, sites (WAY, ACJ and TRE), treatment and traits we are interest in 
data <- data %>%
  filter(site %in% c("WAY", "ACJ", "TRE") &
           year == 2020 & 
           treatment == "C" &
           trait %in% c("plant_height_cm", "leaf_area_cm2", "sla_cm2_g", "ldmc", "dry_mass_g", "leaf_thickness_mm")&
           taxon %in% c("Halenia umbellata","Lachemilla orbiculata","Paspalum bonplandianum",
                        "Rhynchospora macrochaeta","Vaccinium floribundum","Gaultheria glomerata"))


data$site <- factor(data$site)
data$taxon <- factor(data$taxon)
data$functional_group <- factor(data$functional_group)
# data$functional_group <- as.numeric(data$functional_group)
# data <- data[!is.na(data$individual_nr),]

unique(data$site)
# unique(data$individual_nr)
unique(data$taxon)
unique(data$functional_group)
unique(data$trait)


### PCA analysis -------------------------------------------------------------------
df1 <- data %>% 
  dplyr::select(site,id:value,elevation,plot_id) %>% 
  filter(!id %in% c("CXX4125", "BDN3235")) %>% # cut off two outliers 
  tidyr::pivot_wider(names_from = trait, values_from = value) %>% 
  tibble::column_to_rownames("id") 

df1 <- na.omit(df1)

# create a new DF for pca analisys
pca_out <- prcomp(df1[,c(7:12)], center = TRUE, scale = TRUE)


# pca_out <- prcomp(df1, center = TRUE, scale. = TRUE)
summary(pca_out)
str(pca_out)
pca_out$rotation #look at laodings 


# Plot PCA ggplot2
scores <- as.data.frame(pca_out$x)## getting the scores 
scores.1 <- cbind(scores, df1)
pca.loadings <- data.frame(Variables = rownames(pca_out$rotation), pca_out$rotation) # drawing the arrows


# Plot taking taxon into account
p <- ggplot(data = scores.1, aes(x = PC1, y = PC2, color=taxon)) + 
  geom_point(size=2) + 
  scale_fill_hue(l=40) + 
  coord_fixed(ratio=1, xlim=range(scores$PC1), 
              ylim=range(scores$PC2))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+theme_classic()+
  xlab("PC 1 (49.4%)")+
  ylab("PC 2 (31.4%)")+
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*3.5),
                                        yend = (PC2*2)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 3) +
  annotate("text", x = (pca.loadings$PC1*3.5), y = (pca.loadings$PC2*2),
           label = c("Height","Dry mass","leaf area","SLA","LDMC","Leaf thickness"))+
  theme(legend.position="right")+
  scale_color_brewer(palette="Paired")
p+guides(color=guide_legend(title="Plant species"))


# Plot taking sites into account
s <- ggplot(data = scores.1, aes(x = PC1, y = PC2, color=site)) + 
  geom_point(size=2) + 
  scale_fill_hue(l=40) + 
  coord_fixed(ratio=1, xlim=range(scores$PC1), 
              ylim=range(scores$PC2))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+theme_classic()+
  xlab("PC 1 (49.4%)")+
  ylab("PC 2 (31.4%)")+
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*3.5),
                                        yend = (PC2*2)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 3) +
  annotate("text", x = (pca.loadings$PC1*3.5), y = (pca.loadings$PC2*2),
           label = c("Height","Dry mass","leaf area","SLA","LDMC","Leaf thickness"))+
  theme(legend.position="right")+
  scale_color_brewer(palette="Paired")
s+guides(color=guide_legend(title="Sites"))


# Combine SITE and TAXON in the same graph 
#* need to solve color and legend 
c <- ggplot(data = scores.1, aes(x = PC1, y = PC2, color= site)) +
  stat_ellipse(data = scores.1, aes(x = PC1, y = PC2, color=taxon),  alpha = 0.5) +
  coord_fixed(ratio=1, xlim=range(scores$PC1), 
              ylim=range(scores$PC2))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+theme_classic()+
  xlab("PC 1 (49.4%)")+
  ylab("PC 2 (31.4%)")+
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*3.5),
                                        yend = (PC2*2)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 2, alpha = 0.5) +
  annotate("text", x = (pca.loadings$PC1*3.5), y = (pca.loadings$PC2*2),
           label = c("Height","Dry mass","leaf area","SLA","LDMC","Leaf thickness"))+
  theme(legend.position="right")+
  scale_color_brewer(palette="Paired")
c + guides(color=guide_legend(title="Plant species")) 
# stat_ellipse(data = scores.1, aes(x = PC1, y = PC2, color=taxon), show.legend = TRUE) +
# scale_color_brewer(palette = "Greens")

### RDA analysis ---------------------------------------------------------------------
# to see if there are significant differences in the factor 

RDA_out <- rda(df1[,c(7:12)]~scale(elevation) * taxon, center = TRUE, scale = TRUE, data=df1)

summary(RDA_out)
print(RDA_out)

#We construct the models
RDA_1<- rda(df1[,c(7:12)]~scale(elevation), center = TRUE, scale = TRUE, data=df1)
RDA_2<- rda(df1[,c(7:12)]~taxon, center = TRUE, scale = TRUE, data=df1)
RDA_3<- rda(df1[,c(7:12)]~scale(elevation) + taxon, center = TRUE, scale = TRUE, data=df1)
RDA_4<- rda(df1[,c(7:12)]~scale(elevation) * taxon, center = TRUE, scale = TRUE, data=df1)

#then we compare the models
anova(RDA_1, RDA_3)
anova(RDA_1, RDA_4)
anova(RDA_3, RDA_4)

anova(RDA_2, RDA_3)
anova(RDA_2, RDA_4)
anova(RDA_3, RDA_4)

# when doing ANOVA the interaction between elevation and taxa
# there is an interaction between elevation and factor 



### NOT USED -----------------------------------------------------
## Plot with ggbiplot
# ggbiplot(pca_out, choices = c(1, 2)) + theme_classic()
# # ggbiplot(pca_out, choices = c(3, 4)) + theme_classic()
# 
# 
# ggbiplot(
#   pca_out,
#   scale = 0,
#   choices = c(1, 2),
#   ellipse = TRUE,
#   # groups = df1.trait
# ) + theme_classic()

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

