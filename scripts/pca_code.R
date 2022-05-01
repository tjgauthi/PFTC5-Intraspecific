#### Data analysis - Ordination ####

## Authors: Fernanda, Elisa, Korina, Augustina, Fiorella 

#### Load libraries ####
library(dplyr)
library(ggplot2)
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(readr)
library(factoextra)

#* TO DO: 
#* 1- try use ggplot2 intead of bbplot
#* 2- plot sites PCA (maybe like ellipses?!) 
#* 3 plot species (colored or similar)
#* 
#* PCA - use traits as  explanatory variables
#* How traits values vary between different sites? 


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

# should we log transform some of the traits?

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


### PCA - plot 1 
df1 <- data %>% 
  dplyr::select(site,id:value) %>% 
  filter(!id %in% c("CXX4125", "BDN3235")) %>% # cut off two outliers 
  tidyr::pivot_wider(names_from = trait, values_from = value) %>% 
  tibble::column_to_rownames("id") 

df1=na.omit(df1)
#create a new DF for pca analisys

# df1.trait <- factor(data$trait)
pca_out <- prcomp(df1[,-c(1:4)], center = TRUE, scale = TRUE)

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

###PCA with ggplot2

# Plot with GGPLOT

scores <- as.data.frame(pca_out$x)## getting the scores  
#uno el df:scores a pca.com
scores.1=cbind(scores,df1)

# drawing the arrows
pca.loadings <- data.frame(Variables = rownames(pca_out$rotation), pca_out$rotation)



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


#taking sites into account


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