
# Plant functional trait course 5
# Cusco/Wayqecha, Peru - March 2020
#
# Group 1: Intraspecific Variation

#Install dev tools before running this and other codes!! 

### 0) Preamble ----
### >> a) Dependencies ----
# if(!require(skimr)){        # for quick overview of dataset
#   install.packages("skimr")
#   library(skimr)
# }

install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "tidyverse",
  "tidylog",
  "stringr",
  "here",
  "gsheet"
)
sapply(package_vec, install.load.package)

if("dataDownloader" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  devtools::install_github("Between-the-Fjords/dataDownloader")
}
library(dataDownloader)


### >> b) Data from osf ----

dir.create("data")
dir.create("data/raw")

#Download traits data from OSF
get_file(node = "gs8u6",
         file = "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv",
         path = "data/raw",
         remote_path = "traits")

### 1) Data cleaning ----

### >> Traits data ----

# traits data - complete
traits_raw <- read.csv(file.path("data", "raw", "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"),
                       header = T,
                       sep = ",") %>%
  filter(site %in% c("WAY", "ACJ", "TRE") &
           year == 2020 & treatment == "C")

#skim(traits_raw)


### 2) Data filtering ----

# Relevant species only
unique(traits_raw$taxon)
rel_sp <- c("Gaultheria glomerata","Paspalum bonplandianum",
            "Vaccinium floribundum","Rhynchospora macrochaeta",
            # These were new substitutes
            # Check other docs for consistancy
            "Halenia umbellata", "Lachemilla orbiculata") 

traits <- traits_raw %>%
  filter(taxon %in% rel_sp) 

unique(traits$taxon)

#remove raw files
rm('traits_raw')

### 3) Data Structuring ----

#removing obsolete columns after filtering
traits<-traits %>% 
  select (-c(year,season,month,treatment,burn_year,latitude,longitude,course))

#Transform from long to wide format
traits_wide<-traits %>%
  pivot_wider(names_from = trait, values_from = value)

#adding unique plot, individual, and leaf
traits_wide$plot_uid <- paste(traits_wide$site,traits_wide$plot_id, sep = "_")
traits_wide$individual_uid <- paste(traits_wide$site,traits_wide$plot_id, traits_wide$individual_nr, sep = "_")
traits_wide$leaf_uid <- paste(traits_wide$site,traits_wide$plot_id, traits_wide$individual_nr, traits_wide$id, sep = "_")




#Ordination plot: PCA
#Principal component analysis of the variance of SEVEN leaf functional traits 

#call packages
library(vegan)
library(tidyverse)
library(devtools)
library(factoextra)
library(ggbiplot)


#
traits.active <- traits_wide[, c(1,5,7,8:15,18)]
traits.active=na.omit(traits.active)
res.pca <- prcomp(traits.active[,5:11], scale = TRUE)
head(traits.active[, 1:12])
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


mc.pca <- prcomp(traits.active[,c(5:11)], center = TRUE,scale. = TRUE)


# End of script ----
#Add environmental variables 