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
         file = "PFTC3-Puna-PFTC5_Peru_2018-2020_FunctionalTraits_clean.csv",
         path = "data/raw",
         remote_path = "traits")

### 1) Data cleaning ----

### >> Traits data ----

# traits data - complete
traits_raw <- read.csv(file.path("data", "raw", "PFTC3-Puna-PFTC5_Peru_2018-2020_FunctionalTraits_clean.csv"),
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
  filter(!is.na(leaf_nr)) |> 
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
         leaf_uid = paste(site, taxon, plot_id, individual_nr, leaf_nr, sep = "_"))

#Code for cleaning leftover mistakes in the plant height
#This will be changed in the original cleaning code, so it will be redundant once that has been pushed and merged, and the new clean data is on OSF.
#But it doesn't cause any problems if this line of code is here then.

traits_wide <- traits_wide |> 
  mutate(
    plant_height_cm = if_else(id == "ABQ1404", 7.5, plant_height_cm,),
    plant_height_cm = if_else(id %in% c("AUP3248", "AUO6988"), 18, plant_height_cm,),
    plant_height_cm = if_else(id == "AEW0937", 16.4, plant_height_cm),
    plant_height_cm = if_else(id == "CAE6952", 13, plant_height_cm),
    plant_height_cm = if_else(id == "CER9449", 59, plant_height_cm),
    plant_height_cm = if_else(id == "AOR3155", 63.5, plant_height_cm))

rm('traits')

#we only want to use leaves with complete traits
traits_wide<-na.omit(traits_wide) #remove incomplete rows

#filter out individuals that have 1-2 leaves only
traits_wide <-subset (traits_wide, individual_uid != "TRE_Vaccinium floribundum_5_4" &
                        individual_uid !="TRE_Vaccinium floribundum_5_1" &
                        individual_uid !="TRE_Vaccinium floribundum_3_3" &
                        individual_uid !="TRE_Vaccinium floribundum_1_2" &
                        individual_uid !="TRE_Vaccinium floribundum_1_1" &
                        individual_uid !="TRE_Lachemilla orbiculata_1_12" &
                        individual_uid !="ACJ_Rhynchospora macrochaeta_5_3" &
                        individual_uid !="ACJ_Rhynchospora macrochaeta_1_3" &
                        individual_uid !="WAY_Vaccinium floribundum_5_5" &
                        individual_uid !="WAY_Lachemilla orbiculata_1_2" &
                        individual_uid !="TRE_Vaccinium floribundum_5_2" &
                        individual_uid !="TRE_Paspalum bonplandianum_4_1" &
                        individual_uid !="TRE_Lachemilla orbiculata_5_12" &
                        individual_uid !="TRE_Lachemilla orbiculata_1_1")


#removing outliers from the wetmass-drymass relationship
traits_wide <-subset (traits_wide, leaf_uid != "ACJ_Gaultheria glomerata_4_2_1" &
                        leaf_uid != "ACJ_Lachemilla orbiculata_4_3_2" &
                        leaf_uid != "TRE_Vaccinium floribundum_4_2_1" &
                        leaf_uid != "ACJ_Vaccinium floribundum_3_2_1" &
                        leaf_uid != "ACJ_Paspalum bonplandianum_5_3_3" & 
                        leaf_uid != "ACJ_Halenia umbellata_3_3_5")

#removing outlier leaf thickness
traits_wide <-subset (traits_wide, leaf_uid != "ACJ_Paspalum bonplandianum_2_3_1")

traits <- traits_wide |> 
  pivot_longer(names_to = "trait", values_to = "value", cols = 10:16)


# End of script ----

