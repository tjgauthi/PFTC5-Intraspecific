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
                     text = element_text (size = 15), #sets all text size to 20 
                     axis.text = element_text(size = 12)) #sets axis text to size 15  


#check out how many leaves each individual has
leaf_count <- traits_wide %>% 
  group_by(site, plot_uid, individual_uid, taxon) %>% 
  summarize(n = n())

table(leaf_count$n)

#need to decide on how to deal with individuals with 1-2 leaves
traits %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

#in traits, we need to log transform the traits that should be transformed, and create individual uid column.
traits_log <- traits %>% 
  mutate(value = log(value))

output <- data.frame(NULL)
for(i in unique(traits_log$trait)){
  #This code all comes from Julie Messier's web site
  mod<-lmer(value~1+(1|functional_group/taxon/individual_nr)+(1|site), 
            data=traits_log %>% filter(trait == i), 
            na.action=na.omit)
  variances<-c(unlist(lapply(VarCorr(mod),diag)), 
               attr(VarCorr(mod),"sc")^2) #get variances
  
  var.comp<-variances/sum(variances)
  
  var.comp<-as.data.frame(var.comp) #creates a dataframe from the values
  var.comp<-cbind(rownames(var.comp),data.frame(var.comp,row.names=NULL)) #changes row names into a column
  var.comp<-melt(var.comp,value.name="value") #makes var.comp into a variable
  names(var.comp)[1]<-"Scale" #changes the first column name to "scale"
  
  var.comp$value<-var.comp$value *100 #changes values into % 
  
  var.comp<- var.comp %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "site.(Intercept)","individual_nr:(taxon:functional_group).(Intercept)", "taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
    # group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  output <- bind_rows(output, var.comp)
}

# 
# #this is the same as: print(VarCorr(mod),comp="Variance")
# 
# # raw variance
# variances
# 
# # % Variance
# 
# (var.comp<-variances/sum(variances))



VP_Plot<-ggplot(output, aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 0)),colour="white", size = 5)+
  #scale_fill_manual (values = scales_colours)+
  ylab("Proportion of Variance (%)")+
  xlab("")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  #theme(legend.position = "none")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  #scale_fill_manual(values = pal_vp) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0))+
  labs(title = "Focal Species")
VP_Plot

#### look at variance partitioning for all species
#randomly subsample from the 6 species of interest from part 1. How should it be subsampled, though? 


#in traits, we need to log transform the traits that should be transformed, and create individual uid column.
leaf_count_all <- traits_all_wide %>% 
  group_by(site, plot_uid, individual_uid, taxon) %>% 
  summarize(n = n())

table(leaf_count_all$n)


traits_all_log <- traits_all %>% 
  mutate(value = log(value)) %>% 
  separate(taxon, into = c("genus", "species"), sep = " ", remove = F) %>% 
  group_by(site, plot_id, individual_nr, taxon, trait, family, genus) %>% 
  slice_sample(n = 3) #subset down to at most 3 measurements of a species trait from each plot, not sure if this is grouped how we want for subsetting.

output_all <- data.frame(NULL)
for(i in unique(traits_all_log$trait)){
  #This code all comes from Julie Messier's web site
  mod<-lmer(value~1+(1|family/genus/taxon)+(1|site), 
            data=traits_all_log %>% filter(trait == i), 
            na.action=na.omit)
  variances<-c(unlist(lapply(VarCorr(mod),diag)), 
               attr(VarCorr(mod),"sc")^2) #get variances
  
  var.comp<-variances/sum(variances)
  
  var.comp<-as.data.frame(var.comp) #creates a dataframe from the values
  var.comp<-cbind(rownames(var.comp),data.frame(var.comp,row.names=NULL)) #changes row names into a column
  var.comp<-melt(var.comp,value.name="value") #makes var.comp into a variable
  names(var.comp)[1]<-"Scale" #changes the first column name to "scale"
  
  var.comp$value<-var.comp$value *100 #changes values into % 
  
  var.comp<- var.comp %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "site.(Intercept)","taxon:(genus:family).(Intercept)", "genus:family.(Intercept)","family.(Intercept)"))) %>% 
    # group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  output_all <- bind_rows(output_all, var.comp)
}


VP_Plot_all<-ggplot(output_all, aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 0)),colour="white", size = 5)+
  #scale_fill_manual (values = scales_colours)+
  ylab("Proportion of Variance (%)")+
  xlab("")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  #theme(legend.position = "none")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  #scale_fill_manual(values = pal_vp) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0)) +
  labs(title = "All Species")
VP_Plot_all
