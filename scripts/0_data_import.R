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
                       sep = ",") |> 
  filter(site %in% c("WAY", "ACJ", "TRE") &
           year == 2020 & treatment == "C")
#skim(traits_raw)


### 2) Data filtering ----

traits <- traits_raw |> 
  #Select the intraspecific species
  filter(taxon %in% c("Gaultheria glomerata", "Rhynchospora macrochaeta", "Vaccinium floribundum", "Halenia umbellata", "Lachemilla orbiculata", "Paspalum bonplandianum")) |>
  #Removing all individuals of these species that were not sampled with the ITV method (several leaves per individual)
  filter(!is.na(leaf_id)) |>  
  #Those that we can not confirm is wrong, but most likely does not belong in the ITV dataset
  filter(!id %in% c("COI1685", "BUS1756", "CMR2436", "AUB2849", "AAF7186", "BZH3536"))


unique(traits$taxon)

#remove raw files
rm('traits_raw')

### 3) Data Structuring ----

#removing obsolete columns after filtering
traits<-traits |>  
  select (-c(year,season,month,treatment,burn_year,latitude,longitude,course))

#Transform from long to wide format
traits_wide<-traits |> 
  pivot_wider(names_from = trait, values_from = value)

#adding unique plot, individual, and leaf
traits_wide <- traits_wide |> 
  mutate(plot_uid = paste(site, plot_id, sep = "_"),
         individual_uid = paste(site, taxon, plot_id, individual_nr, sep = "_"),
         leaf_uid = paste(site, taxon, plot_id, individual_nr, leaf_id, sep = "_"))

# End of script ----