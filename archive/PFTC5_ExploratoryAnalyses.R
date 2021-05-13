#loading libraries
library(tidyverse)

#verifying that there is only one leaf for each leaf#/individual#/species/plot/site
leaf_count <- intraspecific_data %>%
  group_by(Site, Plot_ID, genus_species, Individual_nr, Leaf_nr) %>%
  summarise(n())

#Calculating the number of individuals/species/site
individual_count <- intraspecific_data %>%
  group_by(Site, Plot_ID, genus_species, Individual_nr) %>%
  summarise(n()) %>%
  group_by(Site, genus_species) %>%
  summarise(n())

####### Exploratory analyses ######
# histogram of species Paspalum
x<-intraspecific_data %>%
  filter(genus_species=="Paspalum bonplandianum")
x
hist(intraspecific_data$Plant_Height_cm)
hist(intraspecific_data$Wet_Mass_g)






