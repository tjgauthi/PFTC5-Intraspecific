#loading libraries
library(openxlsx)


#loading excel File
raw_data<- read.xlsx("C:/Users/Tasha-Leigh/Dropbox/PhD/Courses/PFTC 5 - Peru/PFTC5-Intraspecific/PFTC5_Peru_2020_LeafTraits_cleaned_20-03-19.xlsx", sheet = "PFTC5_Peru_2020_LeafTraits_clea")

#subset to get trait data for the control plots
trait_control <- subset(raw_data, (Project == "T") & (Experiment == "C"))

#combine genus and species into one column
trait_control$genus_species <- paste (trait_control$Genus,trait_control$Species, sep = " ")

#calculate average leaf thickness
trait_control$Leaf_Thickness_AVG <-rowMeans(cbind (x1 = trait_control$Leaf_Thickness_1_mm, x2 = trait_control$Leaf_Thickness_2_mm, x3 = trait_control$Leaf_Thickness_3_mm), na.rm=TRUE)



#partitioning of variance


#varcomp.Height <- varcomp(lme(log(LMA),random=1|Elevation/Species/intraspecific/intraindividual, data = data, na.action = na.omit))
