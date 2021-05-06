#loading libraries
library(readxl)
library(dplyr)

#loading excel File with the traits
# You will need to set your own file path
raw_data <- read_excel("C:/Users/Tasha-Leigh/Dropbox/PhD/Courses/PFTC 5 - Peru/PFTC5-Intraspecific/PFTC5_Peru_2020_LeafTraits_cleaned_20-03-19.xlsx", 
                       col_types = c("text", "numeric", "text", 
                                     "text", "text", "text", "text", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "text", "numeric", "numeric"))

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



#########  Analysis #########

#calculate average leaf thickness
intraspecific_data$Leaf_Thickness_AVG <-rowMeans(cbind (x1 = intraspecific_data$Leaf_Thickness_1_mm, x2 = intraspecific_data$Leaf_Thickness_2_mm, x3 = intraspecific_data$Leaf_Thickness_3_mm), na.rm=TRUE)





#partitioning of variance
#varcomp.Height <- varcomp(lme(log(LMA),random=1|Elevation/Species/intraspecific/intraindividual, data = data, na.action = na.omit))

