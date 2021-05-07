#loading libraries
library(readxl)
library(dplyr)


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

