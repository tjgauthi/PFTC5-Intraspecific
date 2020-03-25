#import excel
#loading libraries
library(readxl)
library(dplyr)
library(tidyverse)
#loading excel File with the traits
# You will need to set your own file path
raw_data <- read_excel("C:/Users/Agustina Barros/Documents/TRAITS/Traitgroup1/data/PFTC5_Peru_2020_LeafTraits_cleaned_20-03-19.xlsx",
                       col_types = c("text", "numeric", "text",
                                     "text", "text", "text", "text", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "text", "numeric", "numeric"))
raw_data
######### Data cleaning and validation #########
#need to rename vaccinium bonplandianum as vaccinium floribundum
raw_data$Species[raw_data$Genus == "Vaccinium" & raw_data$Species == "bonplandianum"] <- "floribundum"
#get rid of duplicates and triplicates based on the ID column
raw_data <-raw_data[-which(duplicated(raw_data[,1])),]
#combine genus and species into one column
raw_data$genus_species <- paste (raw_data$Genus,raw_data$Species, sep = " ")
#subset to get trait data for the trait project, control experiment, and the intraspecific species
intraspecific_data <- subset(raw_data, (Project == "T") & (Experiment == "C") &
                               (Site == "TRE" |
                                  Site == "WAY" |
                                  Site == "ACJ") &
                               (genus_species == "Paspalum bonplandianum"|
                                  genus_species == "Gaultheria glomerata"|
                                  genus_species == "Halenia umbellata"|
                                  genus_species == "Lachemilla orbiculata"|
                                  genus_species == "Rhynchospora macrochaeta"|
                                  genus_species == "Vaccinium floribundum"))
#verifying that there is only one leaf for each leaf#/individual#/species/plot/site
leaf_count <- intraspecific_data %>%
  group_by(Site, Plot_ID, genus_species, Individual_nr, Leaf_nr) %>%
  summarise(n())
#verifying the number of indiduals/species/site
individual_count <- intraspecific_data %>%
  group_by(Site, Plot_ID, genus_species, Individual_nr) %>%
  summarise(n()) %>%
  group_by(Site, genus_species) %>%
  summarise(n())
#calculate average leaf thickness
intraspecific_data$Leaf_Thickness_AVG <-rowMeans(cbind (x1 = intraspecific_data$Leaf_Thickness_1_mm, x2 = intraspecific_data$Leaf_Thickness_2_mm, x3 = intraspecific_data$Leaf_Thickness_3_mm), na.rm=TRUE)
######### Analysis #########
####### Exploratory analyses######
# histogram of species Paspalum
x<-intraspecific_data %>%
  filter(genus_species=="Paspalum bonplandianum")
x
hist(intraspecific_data$Plant_Height_cm)
hist(intraspecific_data$Wet_Mass_g)