# Plant functional trait course 5
# Cusco/Wayqecha, Peru - March 2020
#
# Group 1: Intraspecific Variation


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

rel_sp <- c("Gaultheria glomerata", "Paspalum bonplandianum",
            "Vaccinium floribundium", "Hypericum andium", 
            "Rhynchospora macrochaeta", "Oritropium hieraciodies")

traits_ITV <- traits_raw %>%
  filter(taxon == rel_sp) 

#remove raw files
rm('traits_raw')

# End of script ----