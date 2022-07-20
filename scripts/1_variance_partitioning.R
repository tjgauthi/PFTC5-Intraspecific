# Task DA.1 Perform variance analysis of traits

#### Call source script----

source(here::here(path = "scripts/0_data_import.R"))


#### Setup ----

#load in packages
library(tidyverse)
library(lme4)
library(reshape2)
library(cowplot) #to arrange multiple plots in a figure
library(gcookbook)

#sets a theme
blank_theme <- theme(panel.grid.major = element_blank(), #removes major axis grid lines
                     panel.grid.minor = element_blank(), #removes minor axis grid lines
                     panel.background = element_blank(), #removes the default grey background
                     legend.key = element_blank(), #removes background behind legend keys
                     axis.line = element_line(colour = "black"), #makes axis lines black
                     text = element_text (size = 15), #sets all text size to 20 
                     axis.text = element_text(size = 12)) #sets axis text to size 15  



#### Data Organization ####
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

#### Model Structure 1 - functional group/taxon/individual + 1|site ####
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
# var.comp<-variances/sum(variances))

#Variance Partitioning Plot1
VP_Plot<-ggplot(output, aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 0)),colour="white", size = 5)+
  scale_fill_discrete(labels = c("Within Individual + Unexplained", 
                                 "Between sites",
                                 "Between individuals within taxon",
                                 "Between taxon within functional groups", 
                                 "Between functional groups"))+
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

#### Model Structure 2 - functional group/taxon/site/individual####

output2 <- data.frame(NULL)
for(i in unique(traits_log$trait)){
  #This code all comes from Julie Messier's web site
  mod2<-lmer(value~1+(1|functional_group/taxon/site/individual_nr), 
            data=traits_log %>% filter(trait == i), 
            na.action=na.omit)
  variances2<-c(unlist(lapply(VarCorr(mod2),diag)), 
               attr(VarCorr(mod2),"sc")^2) #get variances
  
  var.comp2<-variances2/sum(variances2)
  
  var.comp2<-as.data.frame(var.comp2) #creates a dataframe from the values
  var.comp2<-cbind(rownames(var.comp2),data.frame(var.comp2,row.names=NULL)) #changes row names into a column
  var.comp2<-melt(var.comp2,value.name="value") #makes var.comp into a variable
  names(var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  var.comp2$value<-var.comp2$value *100 #changes values into % 
  
  var.comp2<- var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_nr:(site:(taxon:functional_group)).(Intercept)", "site:(taxon:functional_group).(Intercept)","taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
    # group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  output2 <- bind_rows(output2, var.comp2)
}

#Variance Partitioning Plot 2 
VP_Plot2<-ggplot(output2, aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 0)),colour="white", size = 5)+
  scale_fill_manual (values = c("#358420","#53a13e","#213c67","#38537f","#899cbb"),
                     labels = c("Within Individual + Unexplained",
                                "Between individuals within sites",
                                "Between sites within taxon",
                                "Between taxon within functional groups",
                                "Between functional groups"))+
  ylab("Proportion of Variance (%)")+
  xlab("")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0))+
  labs(title = "functional group/taxon/site/individual")
VP_Plot2

#### Using model 2 to create a graph showing intra  vs interspecific variability ####

#Reclassifying the scales to be just Intraspecific and Interspecific
output2.1 <-output2 #creates output 2.1 based on model 2
output2.1[,"intrainter"]<-NA #adds a new column "intrainter
output2.1$intrainter <- ifelse (output2$Scale == "Unexplained"|
                              output2$Scale == "individual_nr:(site:(taxon:functional_group)).(Intercept)",
                              "Intraspecific","Interspecific")

#aggregates the variance scales into intraspecific or interspecific
output2.1<-aggregate(output2.1$value, by=list(trait=output2.1$trait,intrainter=output2.1$intrainter), FUN=sum)
output2.1 <- rename (output2.1,"value" = "x") #renames the output "x" to "value"

#this sets the y position of the graph labels
output2.1 <-output2.1 %>% 
  mutate(intrainter = factor(intrainter, levels = c("Intraspecific","Interspecific"))) %>% 
  arrange(trait, intrainter)%>%
  group_by(trait)%>%
  mutate(label_y=100-(cumsum(value)-0.5*value))

#intraspecific vs intraspecific variation plot
VP_Plot_IntraInter<-ggplot(output2.1, aes(x=trait, y=value))+
  geom_col(aes(fill=intrainter))+
  geom_text(aes(y=label_y,label=round(value,digits = 0)),colour="white", size = 5)+
  scale_fill_manual (values = c("#195e07","#0c244a"))+
  ylab("Proportion of Variance (%)")+
  xlab("")+
  blank_theme+
  labs(fill = "intrainter")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0))+
  labs(title = "Intraspecific vs. Interspecific Variation")
VP_Plot_IntraInter

plot_grid (VP_Plot_IntraInter + theme (legend.position = "none"),
           VP_Plot2+ theme (legend.position = "none"),
           get_legend(VP_Plot_IntraInter + 
                        theme(legend.direction = "vertical", 
                              legend.justification = "center",
                              legend.title = element_blank())),
           get_legend(VP_Plot2 + 
                        theme(legend.direction = "vertical", 
                              legend.justification = "center",
                              legend.title = element_blank())),
           ncol = 2, #assigns the # of display columns
           rel_heights =c(1,0.3)) #assigns relative row height allowing us to make the graph larger and the legend smaller

#### Comparing Model Structure 1 and 2 ####

#this pulls the 2 Vp graphs with no legend, adds appropriate titles and adds legends separately to allow resizing of legend vs. graph
plot_grid (VP_Plot + theme (legend.position = "none")+ labs(title = "functional group/taxon/individual + 1|site"),
           VP_Plot2+ theme (legend.position = "none")+labs(title = "functional group/taxon/site/individual"),
           get_legend(VP_Plot + 
                        theme(legend.direction = "vertical", 
                              legend.justification = "center",
                              legend.title = element_blank())),
           get_legend(VP_Plot2 + 
                        theme(legend.direction = "vertical", 
                              legend.justification = "center",
                              legend.title = element_blank())),
           ncol = 2, #assigns the # of display columns
           rel_heights =c(1,0.3)) #assigns relative row height allowing us to make the graph larger and the legend smaller

#### bootstrapping 95% CI for variance partitioning ####

boot_output <- data.frame(NULL)

for(i in unique(traits_log$trait)){
  
  b_output <- data.frame(NULL)
      
  dat <- traits_log %>% 
      filter(trait ==i)
  
  #create randomly-sampled dataset
  for(j in 1:500){
    samp <- dat[sample(nrow(dat), replace = T, size = nrow(dat)*0.9),]
  
  #This code all comes from Julie Messier's web site
  mod<-lmer(value~1+(1|functional_group/taxon/site/individual_nr), 
             data=samp, 
             na.action=na.omit)
  variances<-c(unlist(lapply(VarCorr(mod),diag)), 
                attr(VarCorr(mod),"sc")^2) #get variances
  
  var.comp<-variances/sum(variances)
  
  var.comp<-as.data.frame(t(var.comp)) #creates a dataframe from the values
  colnames(var.comp) <- c("Between individuals within sites",
                          "Between sites within taxon",
                          "Between taxon within functional groups", 
                          "Between functional groups",
                          "Unexplained")
  var.comp$rep <- j
  var.comp$trait <- i
  
  b_output <- bind_rows(b_output, var.comp)
  }
  boot_output <- bind_rows(boot_output, b_output)
}

boot_summary <- boot_output %>% 
  select(-rep) %>% 
  pivot_longer(cols = -trait, names_to = "part", values_to = "vals") %>% 
  group_by(trait, part) %>% 
  summarize(lower = quantile(vals, probs = 0.025), upper = quantile(vals, probs = 0.975))

varpart_real <- output2 %>% 
  mutate(part = plyr::mapvalues(Scale, from = c("Unexplained", "individual_nr:(site:(taxon:functional_group)).(Intercept)", "site:(taxon:functional_group).(Intercept)","taxon:functional_group.(Intercept)","functional_group.(Intercept)"), to = c("Unexplained","Between individuals within sites","Between sites within taxon","Between taxon within functional groups","Between functional groups"))) %>% 
  left_join(boot_summary) %>% 
  mutate(lower = lower*100,
         upper = upper*100)

#### look at variance partitioning for all species ----

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

