#loading libraries
library(openxlsx)


#loading excel File 
# You will need to set your own file path
raw_data<- read.xlsx("C:/Users/Tasha-Leigh/Dropbox/PhD/Courses/PFTC 5 - Peru/PFTC5-Intraspecific/PFTC5_Peru_2020_LeafTraits_cleaned_20-03-19.xlsx", sheet = "PFTC5_Peru_2020_LeafTraits_clea")

#combine genus and species into one column
raw_data$genus_species <- paste (raw_data$Genus,raw_data$Species, sep = " ")

#subset to get trait data for the trait project, control experiment, and the intraspecific species
intraspecific_data <- subset(raw_data, (Project == "T") & (Experiment == "C") & 
                               genus_species == "Gaultheria glomerata" |
                               genus_species == "Halenia umbellata" | 
                               genus_species == "Lachemilla orbiculata" |
                               genus_species == "Paspalum bonplandianum" |
                               genus_species == "Rhynchospora macrochaeta" |
                               genus_species == "Vaccinium floribundum" )

#calculate average leaf thickness
intraspecific_data$Leaf_Thickness_AVG <-rowMeans(cbind (x1 = intraspecific_data$Leaf_Thickness_1_mm, x2 = intraspecific_data$Leaf_Thickness_2_mm, x3 = intraspecific_data$Leaf_Thickness_3_mm), na.rm=TRUE)



#partitioning of variance


#varcomp.Height <- varcomp(lme(log(LMA),random=1|Elevation/Species/intraspecific/intraindividual, data = data, na.action = na.omit))
