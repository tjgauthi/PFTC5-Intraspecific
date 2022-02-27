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
library(tidyverse)
library(tidylog)
if(!require(stringr)){        # for string operations
  install.packages("stringr")
  library(stringr)
}
library(here) #uses working directory as starting point for paths
library(gsheet)
#devtools::install_github("Between-the-Fjords/dataDownloader")
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

### 4) Colour scheme ----

# 3x2 scale (PFGs x species -- for linear models)
pal_lm <- c("#016392", "#A0CBE8", "#E19825", "#F7C480", "#3E8853", "#9FCD99")
# scales::show_col(pal_lm)

# 5-level scale (levels of trait variation -- for variance partitioning)
pal_vp <- c()
# scales::show_col(pal_vp)

# End of script ----
