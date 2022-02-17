# Task DA.1 Perform variance analysis of traits

### Call source script----

source(here::here(path = "scripts/0_data_import.R"))

#load in packages
library(tidyverse)
library(lme4)
library(reshape2)

#sets a theme
blank_theme <- theme(panel.grid.major = element_blank(), #removes major axis grid lines
                     panel.grid.minor = element_blank(), #removes minor axis grid lines
                     panel.background = element_blank(), #removes the default grey background
                     legend.key = element_blank(), #removes background behind legend keys
                     axis.line = element_line(colour = "black"), #makes axis lines black
                     text = element_text (size = 27), #sets all text size to 20 
                     axis.text = element_text(size = 24)) #sets axis text to size 15  


#check out how many leaves each individual has
leaf_count <- traits_wide %>% 
  group_by(site, plot_uid, individual_uid, taxon) %>% 
  summarize(n = n())

table(leaf_count$n)

#need to decide on how to deal with individuals with 1-2 leaves

#This code all comes from Julie Messier's web site
mod<-lmer(log(leaf_thickness_mm)~1+(1|taxon/individual_uid)+(1|site), 
          data=traits_wide, 
          na.action=na.omit)
variances<-c(unlist(lapply(VarCorr(mod),diag)), 
             attr(VarCorr(mod),"sc")^2) #get variances

#this is the same as: print(VarCorr(mod),comp="Variance")

# raw variance
variances

# % Variance

(var.comp<-variances/sum(variances))

var.comp<-as.data.frame(var.comp) #creates a dataframe from the values
var.comp<-cbind(rownames(var.comp),data.frame(var.comp,row.names=NULL)) #changes row names into a column
var.comp<-melt(var.comp,value.name="value") #makes var.comp into a variable
names(var.comp)[1]<-"Scale" #changes the first column name to "scale"

var.comp$value<-var.comp$value *100 #changes values into % 

var.comp<- var.comp %>%
  group_by(variable)%>%
  arrange(variable,desc(Scale))%>%
  mutate(labypos=cumsum(value)-0.5*value)%>%
  subset(value>1) #this line removes variance partitioning less than 1% so that there are no zero labels

VP_Plot<-ggplot(var.comp, aes(x=variable, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 0)),colour="white", size = 10)+
  #scale_fill_manual (values = scales_colours)+
  ylab("Proportion of Variance (%)")+
  xlab("")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100))+#this forces the graph to actually start at 0% and end at 100%
  #theme(legend.position = "none")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})
VP_Plot

