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
library(dplyr)#must be loaded after plyr to run summary

#define custom functions
`%notin%` <- Negate(`%in%`)

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
individual_leaf_count <- traits_wide %>% 
  group_by(site, individual_uid, taxon) %>% 
  summarize(n = n())

#check out how many leaves each site has
site_leaf_count <- traits_wide %>% 
  group_by(site, taxon) %>% 
  summarize(n = n())

#plotting a histogram for each trait to determine if the data is skewed
traits %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

#log transforming traits because of skew
traits_log <- traits %>% 
  mutate(value = log(value))

#plotting log transformed traits to confirm skew is gone
traits_log %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

#### Model Structure 2 - functional group/taxon/site/individual####

#need to remove plant height
traits_filtered <- traits_log %>% 
  filter(trait != "plant_height_cm")

output2 <- data.frame(NULL)
for(i in unique(traits_filtered$trait)){
  #This code all comes from Julie Messier's web site
  mod2<-lmer(value~1+(1|functional_group/taxon/site/individual_uid), 
            data=traits_filtered %>% filter(trait == i), 
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
    mutate(Scale = factor(Scale, levels = c("Unexplained", 
                                            "individual_uid:site:taxon:functional_group.(Intercept)", 
                                            "site:taxon:functional_group.(Intercept)",
                                            "taxon:functional_group.(Intercept)",
                                            "functional_group.(Intercept)"))) %>% 
    # group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  output2 <- bind_rows(output2, var.comp2)
} 

#### model for plant height (no intra-individual variance) ####


PH.model <-lmer(value~1+(1|functional_group/taxon/site), 
            data=traits_log %>% filter(trait == "plant_height_cm"), 
            na.action=na.omit)
PH.variances <-c(unlist(lapply(VarCorr(PH.model),diag)),
                   attr(VarCorr(PH.model),"sc")^2)
PH.var.comp<-PH.variances/sum(PH.variances)


PH.var.comp<-as.data.frame(PH.var.comp) #creates a dataframe from the values
PH.var.comp<-cbind(rownames(PH.var.comp),data.frame(PH.var.comp,row.names=NULL)) #changes row names into a column
PH.var.comp<-melt(PH.var.comp,value.name="value") #makes var.comp into a variable
names(PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

PH.var.comp$value<-PH.var.comp$value *100 #changes values into % 

PH.var.comp<- PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", 
                                          "individual_uid:(site:(taxon:functional_group)).(Intercept)", 
                                          "site:taxon:functional_group.(Intercept)",
                                          "taxon:functional_group.(Intercept)",
                                          "functional_group.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")


#### Plotting VP ####
#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
VP_legend_title <- "Ecological Scales (not plant height)"
VP_Plot2<-ggplot(output2%>% 
                   filter(trait %notin% c("plant_height_cm", "wet_mass_g")), 
                 aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
  scale_fill_manual (VP_legend_title, values = c("#5D3640","#8A5462","#C68596","#8977AA","#4F4266"),
                    #values = c("#77293e","#ae435f","#f27092","#543a83","#37245a"),
                     labels = c("Between leaves within individuals + Unexplained",
                                "Between individuals within sites",
                                "Between sites within species",
                                "Between species within functional groups",
                                "Between functional groups"))+
  ylab("Proportion of Variance (%)")+
  xlab("Log (natural) Transformed Traits")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  scale_x_discrete(limits = c(
                               #"wet_mass_g",
                              "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
                   labels = c(
                               #"wet_mass_g"="Wet Mass (g)", 
                              "dry_mass_g"="Dry Mass (g)",
                              "ldmc"="LDMC", 
                              "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
                              "sla_cm2_g"=expression("SLA "(cm^2/g)), 
                              "leaf_thickness_mm"="Leaf Thickness (mm)")) +
  theme(axis.text.x = element_text(angle = 300, hjust = 0))
  #labs(title = "functional group/taxon/site/individual")
VP_Plot2

#### plant height VP ####
PH_legend_title <- "Plant Height Ecological Scales"

PH.VP_Plot2<-ggplot(PH.var.comp,aes(x=trait, y=value))+
  geom_col(aes(fill=Scale))+
  geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
  scale_fill_manual (PH_legend_title, values = c("#8A5462","#C68596","#8977AA","#4F4266"),
                     #values = c("#77293e","#ae435f","#f27092","#543a83","#37245a"),
                     labels = c("Between individuals within sites + Unexplained",
                                "Between sites within species",
                                "Between species within functional groups",
                                "Between functional groups"))+
  ylab(NULL)+
  #xlab("Log (natural) Transformed Traits")+
  blank_theme+
  labs(fill = "Ecological Scale")+
  scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
  scale_x_discrete(labels = c("Plant Height (cm)")) +
  theme(axis.text.x = element_text(angle = 300, hjust = 0))
PH.VP_Plot2

#### combining plots ####

PH.VP_Plot2.legend<- get_legend(PH.VP_Plot2)  #gets the legend from each plot
VP_Plot2.legend<- get_legend(VP_Plot2)

temp.plots<- plot_grid(VP_Plot2 + theme(legend.position = "none"),
                       PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(VP_Plot2.legend,
                       PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)



#### bootstrapping 95% CI for variance partitioning ####

#traits_filtered is a long dataset of log transformed trait values with plant height already removed

boot_output <- data.frame(NULL)

for(i in unique(traits_filtered$trait)){
  
  b_output <- data.frame(NULL)
      
  dat <- traits_filtered %>% 
      filter(trait ==i)
  
  #create randomly-sampled dataset
  for(j in 1:500){
    samp <- dat[sample(nrow(dat), replace = T, size = nrow(dat)*0.9),]
  
  #This code all comes from Julie Messier's web site
  mod<-lmer(value~1+(1|functional_group/taxon/site/individual_uid), 
             data=samp, 
             na.action=na.omit)
  variances<-c(unlist(lapply(VarCorr(mod),diag)), 
                attr(VarCorr(mod),"sc")^2) #get variances
  
  var.comp<-variances/sum(variances)
  
  var.comp<-as.data.frame(t(var.comp)) #creates a dataframe from the values
  #colnames(var.comp) <- c("Between individuals within sites",
                         # "Between sites within taxon",
                         # "Between taxon within functional groups", 
                          #"Between functional groups",
                         # "Unexplained")
  var.comp$rep <- j
  var.comp$trait <- i
  
  b_output <- bind_rows(b_output, var.comp)
  }
  boot_output <- bind_rows(boot_output, b_output)
}

boot_summary <- boot_output %>% 
  select(-rep) %>% 
  pivot_longer(cols = -trait, names_to = "part", values_to = "vals") %>% 
  dplyr::group_by(trait, part) %>% 
  dplyr::summarise(lower = quantile(vals, probs = 0.025), upper = quantile(vals, probs = 0.975)) %>%
  mutate(lower = lower*100,
         upper = upper*100)


#### bootstrapping 95% CI for VP Plant Height #### 


dat <- traits_log %>% filter(trait == "plant_height_cm")

PH_output <- data.frame(NULL)

#create randomly-sampled dataset
  for(j in 1:500){
    samp <- dat[sample(nrow(dat), replace = T, size = nrow(dat)*0.9),] #creates a sample dataset of 90% of original data
    
    mod<- tryCatch ({lmer(value~1+(1|functional_group/taxon/site), #try catch to remove null models
              data=samp, 
              na.action=na.omit)}, error = function(e) return(NULL))
    if (is.null(mod)) next
      
    variances<-c(unlist(lapply(VarCorr(mod),diag)), 
                 attr(VarCorr(mod),"sc")^2) #get variances
    
    var.comp<-variances/sum(variances)
    
    var.comp<-as.data.frame(t(var.comp)) #creates a dataframe from the values
    var.comp$rep <- j
    var.comp$trait <- i
    
    PH_output <- bind_rows(PH_output, var.comp)
  }

PH_summary <- PH_output %>% 
  select(-rep) %>% 
  pivot_longer(cols = -trait, names_to = "part", values_to = "vals") %>% 
  na.omit() %>%
  dplyr::group_by(trait, part) %>% 
  dplyr::summarise(lower = quantile(vals, probs = 0.025), upper = quantile(vals, probs = 0.975))%>% 
  mutate(lower = lower*100,
         upper = upper*100)

#### ARCHIVE - Model Structure 2 broken down by functional group ####

#Forbs
forb_trait <- traits_log %>% 
  filter (functional_group == "Forb")

forb_output2 <- data.frame(NULL)
for(i in unique(forb_trait$trait)){
  #This code all comes from Julie Messier's web site
  forb_mod2<-lmer(value~1+(1|taxon/site/individual_uid), 
                  data=forb_trait %>% filter(trait == i), 
                  na.action=na.omit)
  forb_variances2<-c(unlist(lapply(VarCorr(forb_mod2),diag)), 
                     attr(VarCorr(forb_mod2),"sc")^2) #get variances
  
  forb_var.comp2<-forb_variances2/sum(forb_variances2)
  
  forb_var.comp2<-as.data.frame(forb_var.comp2) #creates a dataframe from the values
  forb_var.comp2<-cbind(rownames(forb_var.comp2),data.frame(forb_var.comp2,row.names=NULL)) #changes row names into a column
  forb_var.comp2<-melt(forb_var.comp2,value.name="value") #makes var.comp into a variable
  names(forb_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  forb_var.comp2$value<-forb_var.comp2$value *100 #changes values into % 
  
  forb_var.comp2<- forb_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(site:taxon).(Intercept)", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  forb_output2 <- bind_rows(forb_output2, forb_var.comp2)
} 

#model for plant height (no intra-individual variance)
forb_PH.model <-lmer(value~1+(1|taxon/site), 
                     data=forb_trait %>% filter(trait == "plant_height_cm"), 
                     na.action=na.omit)
forb_PH.variances <-c(unlist(lapply(VarCorr(forb_PH.model),diag)),
                      attr(VarCorr(forb_PH.model),"sc")^2)
forb_PH.var.comp<-forb_PH.variances/sum(forb_PH.variances)

forb_PH.var.comp<-as.data.frame(forb_PH.var.comp) #creates a dataframe from the values
forb_PH.var.comp<-cbind(rownames(forb_PH.var.comp),data.frame(forb_PH.var.comp,row.names=NULL)) #changes row names into a column
forb_PH.var.comp<-melt(forb_PH.var.comp,value.name="value") #makes var.comp into a variable
names(forb_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

forb_PH.var.comp$value<-forb_PH.var.comp$value *100 #changes values into % 

forb_PH.var.comp<- forb_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
forb_VP_legend_title <- "Forb Ecological Scales (not plant height)"
(forb_VP_Plot2<-ggplot(forb_output2%>% 
                         filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                       aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (forb_VP_legend_title, values = c("#5D3640","#8A5462","#C68596","#8977AA"),
                       labels = c("Within Individual + Unexplained",
                                  "Between individuals within sites",
                                  "Between sites within taxon",
                                  "Between taxon"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
forb_PH_legend_title <- "Forb Plant Height Ecological Scales"

(forb_PH.VP_Plot2<-ggplot(forb_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (forb_PH_legend_title, values = c("#8A5462","#C68596","#8977AA"),
                       labels = c("Between individuals within sites + Unexplained",
                                  "Between sites within taxon",
                                  "Between taxon within functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

forb_PH.VP_Plot2.legend<- get_legend(forb_PH.VP_Plot2)  #gets the legend from each plot
forb_VP_Plot2.legend<- get_legend(forb_VP_Plot2)

temp.plots<- plot_grid(forb_VP_Plot2 + theme(legend.position = "none"),
                       forb_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(forb_VP_Plot2.legend,
                       forb_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)




##Graminoid

graminoid_trait <- traits_log %>% 
  filter (functional_group == "Graminoid")

graminoid_output2 <- data.frame(NULL)
for(i in unique(graminoid_trait$trait)){
  #This code all comes from Julie Messier's web site
  graminoid_mod2<-lmer(value~1+(1|taxon/site/individual_uid), 
                       data=graminoid_trait %>% filter(trait == i), 
                       na.action=na.omit)
  graminoid_variances2<-c(unlist(lapply(VarCorr(graminoid_mod2),diag)), 
                          attr(VarCorr(graminoid_mod2),"sc")^2) #get variances
  
  graminoid_var.comp2<-graminoid_variances2/sum(graminoid_variances2)
  
  graminoid_var.comp2<-as.data.frame(graminoid_var.comp2) #creates a dataframe from the values
  graminoid_var.comp2<-cbind(rownames(graminoid_var.comp2),data.frame(graminoid_var.comp2,row.names=NULL)) #changes row names into a column
  graminoid_var.comp2<-melt(graminoid_var.comp2,value.name="value") #makes var.comp into a variable
  names(graminoid_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  graminoid_var.comp2$value<-graminoid_var.comp2$value *100 #changes values into % 
  
  graminoid_var.comp2<- graminoid_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(site:taxon).(Intercept)", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  graminoid_output2 <- bind_rows(graminoid_output2, graminoid_var.comp2)
} 

#model for plant height (no intra-individual variance)
graminoid_PH.model <-lmer(value~1+(1|taxon/site), 
                          data=graminoid_trait %>% filter(trait == "plant_height_cm"), 
                          na.action=na.omit)
graminoid_PH.variances <-c(unlist(lapply(VarCorr(graminoid_PH.model),diag)),
                           attr(VarCorr(graminoid_PH.model),"sc")^2)
graminoid_PH.var.comp<-graminoid_PH.variances/sum(graminoid_PH.variances)

graminoid_PH.var.comp<-as.data.frame(graminoid_PH.var.comp) #creates a dataframe from the values
graminoid_PH.var.comp<-cbind(rownames(graminoid_PH.var.comp),data.frame(graminoid_PH.var.comp,row.names=NULL)) #changes row names into a column
graminoid_PH.var.comp<-melt(graminoid_PH.var.comp,value.name="value") #makes var.comp into a variable
names(graminoid_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

graminoid_PH.var.comp$value<-graminoid_PH.var.comp$value *100 #changes values into % 

graminoid_PH.var.comp<- graminoid_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
graminoid_VP_legend_title <- "graminoid Ecological Scales (not plant height)"
(graminoid_VP_Plot2<-ggplot(graminoid_output2%>% 
                              filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                            aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (graminoid_VP_legend_title, values = c("#5D3640","#8A5462","#C68596","#8977AA"),
                       labels = c("Within Individual + Unexplained",
                                  "Between individuals within sites",
                                  "Between sites within taxon",
                                  "Between taxon"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
graminoid_PH_legend_title <- "graminoid Plant Height Ecological Scales"

(graminoid_PH.VP_Plot2<-ggplot(graminoid_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (graminoid_PH_legend_title, values = c("#8A5462","#C68596","#8977AA"),
                       labels = c("Between individuals within sites + Unexplained",
                                  "Between sites within taxon",
                                  "Between taxon within functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

graminoid_PH.VP_Plot2.legend<- get_legend(graminoid_PH.VP_Plot2)  #gets the legend from each plot
graminoid_VP_Plot2.legend<- get_legend(graminoid_VP_Plot2)

temp.plots<- plot_grid(graminoid_VP_Plot2 + theme(legend.position = "none"),
                       graminoid_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(graminoid_VP_Plot2.legend,
                       graminoid_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)


#woody

woody_trait <- traits_log %>% 
  filter (functional_group == "Woody")

woody_output2 <- data.frame(NULL)
for(i in unique(woody_trait$trait)){
  #This code all comes from Julie Messier's web site
  woody_mod2<-lmer(value~1+(1|taxon/site/individual_uid), 
                   data=woody_trait %>% filter(trait == i), 
                   na.action=na.omit)
  woody_variances2<-c(unlist(lapply(VarCorr(woody_mod2),diag)), 
                      attr(VarCorr(woody_mod2),"sc")^2) #get variances
  
  woody_var.comp2<-woody_variances2/sum(woody_variances2)
  
  woody_var.comp2<-as.data.frame(woody_var.comp2) #creates a dataframe from the values
  woody_var.comp2<-cbind(rownames(woody_var.comp2),data.frame(woody_var.comp2,row.names=NULL)) #changes row names into a column
  woody_var.comp2<-melt(woody_var.comp2,value.name="value") #makes var.comp into a variable
  names(woody_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  woody_var.comp2$value<-woody_var.comp2$value *100 #changes values into % 
  
  woody_var.comp2<- woody_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(site:taxon).(Intercept)", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  woody_output2 <- bind_rows(woody_output2, woody_var.comp2)
} 

#model for plant height (no intra-individual variance)
woody_PH.model <-lmer(value~1+(1|taxon/site), 
                      data=woody_trait %>% filter(trait == "plant_height_cm"), 
                      na.action=na.omit)
woody_PH.variances <-c(unlist(lapply(VarCorr(woody_PH.model),diag)),
                       attr(VarCorr(woody_PH.model),"sc")^2)
woody_PH.var.comp<-woody_PH.variances/sum(woody_PH.variances)

woody_PH.var.comp<-as.data.frame(woody_PH.var.comp) #creates a dataframe from the values
woody_PH.var.comp<-cbind(rownames(woody_PH.var.comp),data.frame(woody_PH.var.comp,row.names=NULL)) #changes row names into a column
woody_PH.var.comp<-melt(woody_PH.var.comp,value.name="value") #makes var.comp into a variable
names(woody_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

woody_PH.var.comp$value<-woody_PH.var.comp$value *100 #changes values into % 

woody_PH.var.comp<- woody_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "site:taxon.(Intercept)","taxon.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
woody_VP_legend_title <- "woody Ecological Scales (not plant height)"
(woody_VP_Plot2<-ggplot(woody_output2%>% 
                          filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                        aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (woody_VP_legend_title, values = c("#5D3640","#8A5462","#C68596","#8977AA"),
                       labels = c("Within Individual + Unexplained",
                                  "Between individuals within sites",
                                  "Between sites within taxon",
                                  "Between taxon"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
woody_PH_legend_title <- "woody Plant Height Ecological Scales"

(woody_PH.VP_Plot2<-ggplot(woody_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (woody_PH_legend_title, values = c("#8A5462","#C68596","#8977AA"),
                       labels = c("Between individuals within sites + Unexplained",
                                  "Between sites within taxon",
                                  "Between taxon within functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

woody_PH.VP_Plot2.legend<- get_legend(woody_PH.VP_Plot2)  #gets the legend from each plot
woody_VP_Plot2.legend<- get_legend(woody_VP_Plot2)

temp.plots<- plot_grid(woody_VP_Plot2 + theme(legend.position = "none"),
                       woody_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(woody_VP_Plot2.legend,
                       woody_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

#### ARCHIVE - Model Structure 2 broken down by site ####

ACJ_trait <- traits_log %>% 
  filter (site == "ACJ")

ACJ_output2 <- data.frame(NULL)
for(i in unique(ACJ_trait$trait)){
  #This code all comes from Julie Messier's web site
  ACJ_mod2<-lmer(value~1+(1|functional_group/taxon/individual_uid),
                 data=ACJ_trait %>% filter(trait == i), 
                 na.action=na.omit)
  ACJ_variances2<-c(unlist(lapply(VarCorr(ACJ_mod2),diag)), 
                    attr(VarCorr(ACJ_mod2),"sc")^2) #get variances
  
  ACJ_var.comp2<-ACJ_variances2/sum(ACJ_variances2)
  
  ACJ_var.comp2<-as.data.frame(ACJ_var.comp2) #creates a dataframe from the values
  ACJ_var.comp2<-cbind(rownames(ACJ_var.comp2),data.frame(ACJ_var.comp2,row.names=NULL)) #changes row names into a column
  ACJ_var.comp2<-melt(ACJ_var.comp2,value.name="value") #makes var.comp into a variable
  names(ACJ_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  ACJ_var.comp2$value<-ACJ_var.comp2$value *100 #changes values into % 
  
  ACJ_var.comp2<- ACJ_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(taxon:functional_group).(Intercept)","taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  ACJ_output2 <- bind_rows(ACJ_output2, ACJ_var.comp2)
} 

#model for plant height (no intra-individual variance)
ACJ_PH.model <-lmer(value~1+(1|functional_group/taxon), 
                    data=ACJ_trait %>% filter(trait == "plant_height_cm"), 
                    na.action=na.omit)
ACJ_PH.variances <-c(unlist(lapply(VarCorr(ACJ_PH.model),diag)),
                     attr(VarCorr(ACJ_PH.model),"sc")^2)
ACJ_PH.var.comp<-ACJ_PH.variances/sum(ACJ_PH.variances)

ACJ_PH.var.comp<-as.data.frame(ACJ_PH.var.comp) #creates a dataframe from the values
ACJ_PH.var.comp<-cbind(rownames(ACJ_PH.var.comp),data.frame(ACJ_PH.var.comp,row.names=NULL)) #changes row names into a column
ACJ_PH.var.comp<-melt(ACJ_PH.var.comp,value.name="value") #makes var.comp into a variable
names(ACJ_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

ACJ_PH.var.comp$value<-ACJ_PH.var.comp$value *100 #changes values into % 

ACJ_PH.var.comp<- ACJ_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
ACJ_VP_legend_title <- "ACJ Ecological Scales (not plant height)"
(ACJ_VP_Plot2<-ggplot(ACJ_output2%>% 
                        filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                      aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (ACJ_VP_legend_title, values = c("#5D3640","#8A5462","#8977AA","#4F4266"),
                       labels = c("Between leaves within individuals + Unexplained",
                                  "Between individuals within species",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
ACJ_PH_legend_title <- "ACJ Plant Height Ecological Scales"

(ACJ_PH.VP_Plot2<-ggplot(ACJ_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (ACJ_PH_legend_title, values = c("#8A5462","#8977AA","#4F4266"),
                       labels = c("Between individuals within species + Unexplained",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

ACJ_PH.VP_Plot2.legend<- get_legend(ACJ_PH.VP_Plot2)  #gets the legend from each plot
ACJ_VP_Plot2.legend<- get_legend(ACJ_VP_Plot2)

temp.plots<- plot_grid(ACJ_VP_Plot2 + theme(legend.position = "none"),
                       ACJ_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(ACJ_VP_Plot2.legend,
                       ACJ_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)

### WAY ###

WAY_trait <- traits_log %>% 
  filter (site == "WAY")

WAY_output2 <- data.frame(NULL)
for(i in unique(WAY_trait$trait)){
  #This code all comes from Julie Messier's web site
  WAY_mod2<-lmer(value~1+(1|functional_group/taxon/individual_uid),
                 data=WAY_trait %>% filter(trait == i), 
                 na.action=na.omit)
  WAY_variances2<-c(unlist(lapply(VarCorr(WAY_mod2),diag)), 
                    attr(VarCorr(WAY_mod2),"sc")^2) #get variances
  
  WAY_var.comp2<-WAY_variances2/sum(WAY_variances2)
  
  WAY_var.comp2<-as.data.frame(WAY_var.comp2) #creates a dataframe from the values
  WAY_var.comp2<-cbind(rownames(WAY_var.comp2),data.frame(WAY_var.comp2,row.names=NULL)) #changes row names into a column
  WAY_var.comp2<-melt(WAY_var.comp2,value.name="value") #makes var.comp into a variable
  names(WAY_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  WAY_var.comp2$value<-WAY_var.comp2$value *100 #changes values into % 
  
  WAY_var.comp2<- WAY_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(taxon:functional_group).(Intercept)","taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  WAY_output2 <- bind_rows(WAY_output2, WAY_var.comp2)
} 

#model for plant height (no intra-individual variance)
WAY_PH.model <-lmer(value~1+(1|functional_group/taxon), 
                    data=WAY_trait %>% filter(trait == "plant_height_cm"), 
                    na.action=na.omit)
WAY_PH.variances <-c(unlist(lapply(VarCorr(WAY_PH.model),diag)),
                     attr(VarCorr(WAY_PH.model),"sc")^2)
WAY_PH.var.comp<-WAY_PH.variances/sum(WAY_PH.variances)

WAY_PH.var.comp<-as.data.frame(WAY_PH.var.comp) #creates a dataframe from the values
WAY_PH.var.comp<-cbind(rownames(WAY_PH.var.comp),data.frame(WAY_PH.var.comp,row.names=NULL)) #changes row names into a column
WAY_PH.var.comp<-melt(WAY_PH.var.comp,value.name="value") #makes var.comp into a variable
names(WAY_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

WAY_PH.var.comp$value<-WAY_PH.var.comp$value *100 #changes values into % 

WAY_PH.var.comp<- WAY_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
WAY_VP_legend_title <- "WAY Ecological Scales (not plant height)"
(WAY_VP_Plot2<-ggplot(WAY_output2%>% 
                        filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                      aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (WAY_VP_legend_title, values = c("#5D3640","#8A5462","#8977AA","#4F4266"),
                       labels = c("Between leaves within individuals + Unexplained",
                                  "Between individuals within species",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
WAY_PH_legend_title <- "WAY Plant Height Ecological Scales"

(WAY_PH.VP_Plot2<-ggplot(WAY_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (WAY_PH_legend_title, values = c("#8A5462","#8977AA","#4F4266"),
                       labels = c("Between individuals within species + Unexplained",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

WAY_PH.VP_Plot2.legend<- get_legend(WAY_PH.VP_Plot2)  #gets the legend from each plot
WAY_VP_Plot2.legend<- get_legend(WAY_VP_Plot2)

temp.plots<- plot_grid(WAY_VP_Plot2 + theme(legend.position = "none"),
                       WAY_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(WAY_VP_Plot2.legend,
                       WAY_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)

### TRE ###

TRE_trait <- traits_log %>% 
  filter (site == "TRE")

TRE_output2 <- data.frame(NULL)
for(i in unique(TRE_trait$trait)){
  #This code all comes from Julie Messier's web site
  TRE_mod2<-lmer(value~1+(1|functional_group/taxon/individual_uid),
                 data=TRE_trait %>% filter(trait == i), 
                 na.action=na.omit)
  TRE_variances2<-c(unlist(lapply(VarCorr(TRE_mod2),diag)), 
                    attr(VarCorr(TRE_mod2),"sc")^2) #get variances
  
  TRE_var.comp2<-TRE_variances2/sum(TRE_variances2)
  
  TRE_var.comp2<-as.data.frame(TRE_var.comp2) #creates a dataframe from the values
  TRE_var.comp2<-cbind(rownames(TRE_var.comp2),data.frame(TRE_var.comp2,row.names=NULL)) #changes row names into a column
  TRE_var.comp2<-melt(TRE_var.comp2,value.name="value") #makes var.comp into a variable
  names(TRE_var.comp2)[1]<-"Scale" #changes the first column name to "scale"
  
  TRE_var.comp2$value<-TRE_var.comp2$value *100 #changes values into % 
  
  TRE_var.comp2<- TRE_var.comp2 %>%
    mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
    mutate(Scale = factor(Scale, levels = c("Unexplained", "individual_uid:(taxon:functional_group).(Intercept)","taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
    #group_by(variable)%>%
    arrange(variable, Scale)%>%
    mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
    #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
    mutate(trait = i)
  
  TRE_output2 <- bind_rows(TRE_output2, TRE_var.comp2)
} 

#model for plant height (no intra-individual variance)
TRE_PH.model <-lmer(value~1+(1|functional_group/taxon), 
                    data=TRE_trait %>% filter(trait == "plant_height_cm"), 
                    na.action=na.omit)
TRE_PH.variances <-c(unlist(lapply(VarCorr(TRE_PH.model),diag)),
                     attr(VarCorr(TRE_PH.model),"sc")^2)
TRE_PH.var.comp<-TRE_PH.variances/sum(TRE_PH.variances)

TRE_PH.var.comp<-as.data.frame(TRE_PH.var.comp) #creates a dataframe from the values
TRE_PH.var.comp<-cbind(rownames(TRE_PH.var.comp),data.frame(TRE_PH.var.comp,row.names=NULL)) #changes row names into a column
TRE_PH.var.comp<-melt(TRE_PH.var.comp,value.name="value") #makes var.comp into a variable
names(TRE_PH.var.comp)[1]<-"Scale" #changes the first column name to "scale"

TRE_PH.var.comp$value<-TRE_PH.var.comp$value *100 #changes values into % 

TRE_PH.var.comp<- TRE_PH.var.comp %>%
  mutate(Scale = plyr::mapvalues(Scale, from = c(""), to = c("Unexplained"))) %>% 
  mutate(Scale = factor(Scale, levels = c("Unexplained", "taxon:functional_group.(Intercept)","functional_group.(Intercept)"))) %>% 
  # group_by(variable)%>%
  arrange(variable, Scale)%>%
  mutate(labypos=100-(cumsum(value)-0.5*value)) %>%
  #subset(value>1) %>% #this line removes variance partitioning less than 1% so that there are no zero labels
  mutate(trait = "plant_height_cm")



#Colour palette for graph: https://paletton.com/#uid=54n140kjJo3hfJliyuOl7gxlT9k
#Colour palette is in low saturation 
#Variance Partitioning Plot 2 
TRE_VP_legend_title <- "TRE Ecological Scales (not plant height)"
(TRE_VP_Plot2<-ggplot(TRE_output2%>% 
                        filter(trait %notin% c("wet_mass_g","plant_height_cm")), 
                      aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (TRE_VP_legend_title, values = c("#5D3640","#8A5462","#8977AA","#4F4266"),
                       labels = c("Between leaves within individuals + Unexplained",
                                  "Between individuals within species",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab("Proportion of Variance (%)")+
    xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(limits = c(
      # "wet_mass_g",
      "dry_mass_g", "ldmc", "leaf_area_cm2","sla_cm2_g","leaf_thickness_mm"),
      labels = c(
        # "wet_mass_g"="Wet Mass (g)", 
        "dry_mass_g"="Dry Mass (g)",
        "ldmc"="LDMC", 
        "leaf_area_cm2"=expression(paste("Leaf Area "(cm^2))),
        "sla_cm2_g"=expression("SLA "(cm^2/g)), 
        "leaf_thickness_mm"="Leaf Thickness (mm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))
#labs(title = "functional group/taxon/site/individual")


#plant height VP
TRE_PH_legend_title <- "TRE Plant Height Ecological Scales"

(TRE_PH.VP_Plot2<-ggplot(TRE_PH.var.comp,aes(x=trait, y=value))+
    geom_col(aes(fill=Scale))+
    geom_text(aes(y=labypos, label=round(value,digits = 1)),colour="white", size = 5)+
    scale_fill_manual (TRE_PH_legend_title, values = c("#8A5462","#8977AA","#4F4266"),
                       labels = c("Between individuals within species + Unexplained",
                                  "Between species within functional groups",
                                  "Between functional groups"))+
    ylab(NULL)+
    #xlab("Log (natural) Transformed Traits")+
    blank_theme+
    labs(fill = "Ecological Scale")+
    scale_y_continuous(expand=c(0,0),limits=c(0,100.1))+#this forces the graph to actually start at 0% and end at 100%
    scale_x_discrete(labels = c("Plant Height (cm)")) +
    theme(axis.text.x = element_text(angle = 300, hjust = 0)))


#combining plots

TRE_PH.VP_Plot2.legend<- get_legend(TRE_PH.VP_Plot2)  #gets the legend from each plot
TRE_VP_Plot2.legend<- get_legend(TRE_VP_Plot2)

temp.plots<- plot_grid(TRE_VP_Plot2 + theme(legend.position = "none"),
                       TRE_PH.VP_Plot2 + theme(legend.position = "none") + xlab(NULL),
                       ncol=2,
                       align = 'h',
                       rel_widths = c(1,0.3))

temp.legend<-plot_grid(TRE_VP_Plot2.legend,
                       TRE_PH.VP_Plot2.legend,
                       align = 'v',
                       rel_heights = c(0.5,1),
                       nrow = 2)
plot_grid(temp.plots,
          temp.legend,
          rel_widths =  c(1,0.5),
          ncol = 2)

ggsave ("VP.jpeg",
        width = 6.5, height=5, units = "in",
        dpi = 800)

#### ARCHIVE - Model Structure 1 - functional group/taxon/individual + 1|site ####
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
  #scale_fill_discrete(labels = c("Within Individual + Unexplained", 
                                 #"Between sites",
                                # "Between individuals within taxon",
                                # "Between taxon within functional groups", 
                                # "Between functional groups"))+
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


#### ARCHIVE - Comparing Model Structure 1 and 2 ####

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

#### ARCHIVE - Using model 2 to create a graph showing intra  vs interspecific variability ####

#Reclassifying the scales to be just Intraspecific and Interspecific
output2.1 <-output2 #creates output 2.1 based on model 2
output2.1[,"intrainter"]<-NA #adds a new column "intrainter
output2.1$intrainter <- ifelse (output2$Scale == "Unexplained"|
                                  output2$Scale == "individual_uid:(site:(taxon:functional_group)).(Intercept)"|
                                  output2$Scale == "site:(taxon:functional_group).(Intercept)",
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


#### ARCHIVE - look at variance partitioning for all species ----

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




#### ARCHIVE - Exporting Data ####
library(writexl)
write_xlsx(traits_wide,"C:/Users/tjgau/Desktop/Work/PFTC5/test.xlsx")

#### ARCHIVE - checking for outliers ####
traits_log <- traits %>% 
  mutate(value = log(value))

#removing outliers from the wetmass-drymass relationship
traits_log <-subset (traits_log, leaf_uid != "ACJ_Gaultheria glomerata_4_2_1" &
                        leaf_uid != "ACJ_Lachemilla orbiculata_4_3_2" &
                        leaf_uid != "TRE_Vaccinium floribundum_4_2_1" &
                        leaf_uid != "ACJ_Vaccinium floribundum_3_2_1" &
                        leaf_uid != "ACJ_Paspalum bonplandianum_5_3_3" & 
                        leaf_uid != "ACJ_Halenia umbellata_3_3_5")

traits_log <-subset (traits_log, leaf_uid != "ACJ_Paspalum bonplandianum_2_3_1") #removed due to outlier leaf thickness

#checking wet vs dry mass
ggplot(data = traits_wide, aes(x = dry_mass_g, y = leaf_area_cm2)) +
  geom_point()+
  facet_wrap(~ taxon, scales = "free") +
  blank_theme



ggplot(data = traits_log, aes(x = functional_group, y = value)) +
  geom_jitter(alpha=0.2, aes(colour=taxon))+
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ trait, scales = "free_y") +
  blank_theme


#### ARCHIVE - selecting individuals with only 4 leaves ####


#filter out only vaccinium that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "TRE_Vaccinium floribundum_3_1")

#removing outliers from the wetmass-drymass relationship
traits_log <-subset (traits_log, leaf_uid != "ACJ_Gaultheria glomerata_4_2_1" &
                       leaf_uid != "ACJ_Lachemilla orbiculata_4_3_2" &
                       leaf_uid != "TRE_Vaccinium floribundum_4_2_1" &
                       leaf_uid != "ACJ_Vaccinium floribundum_3_2_1" &
                       leaf_uid != "ACJ_Paspalum bonplandianum_5_3_3" & 
                       leaf_uid != "ACJ_Halenia umbellata_3_3_5")

#removing outlier leaf thickness
traits_log <-subset (traits_log, leaf_uid != "ACJ_Paspalum bonplandianum_2_3_1")

#check out how many leaves each individual has
individual_leaf_count <- traits_log %>% 
  group_by(site, individual_uid, taxon,trait) %>% 
  summarize(
    leaf_count = n(),
    .groups = "drop"
  )


#subset_data_4 <- traits_wide %>%
#group_by(individual_uid) %>%
# filter individuals with exactly 4 leaves
# filter(n() == 4) 

#selecting individuals with 5 leaves and randomly removing 1 leaf
#subset_data_5 <- traits_wide %>%
#group_by(individual_uid) %>%
# filter individuals with exactly 5 leaves
# filter(n() == 5) 


#%>%
# randomly select 4 leaves per individual
#slice_sample(n = 4)

#Traits_Wide_4 <- bind_rows(subset_data_4 , subset_data_5)

#Traits_long_4 <- Traits_Wide_4 %>%
#pivot_longer(names_to = "trait", values_to = "value", cols = 10:16)

#traits_log <- Traits_long_4 %>% 
#mutate(value = log(value))%>%
# filter(trait != "wet_mass_g" & trait != "plant_height_cm" )

#### ARCHIVE - filtering out specific samples that only have 3 leaves ####
# to see which may be contributing to zero variability in plant height for fuctional groups

#create log traits (delete after testing)
traits_log <- traits %>% 
  mutate(value = log(value))

#filter out only gautheria that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "ACJ_Gaultheria glomerata_4_1" &
                       individual_uid !="TRE_Gaultheria glomerata_4_2" &
                       individual_uid !="WAY_Gaultheria glomerata_4_1")

#filter out only halenia that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "ACJ_Halenia umbellata_1_3" &
                       individual_uid !="WAY_Halenia umbellata_1_1")

#filter out only lachemilla that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "TRE_Lachemilla orbiculata_3_2" &
                       individual_uid !="TRE_Lachemilla orbiculata_4_2" &
                       individual_uid !="TRE_Lachemilla orbiculata_3_1"&
                       individual_uid !="TRE_Lachemilla orbiculata_2_2"&
                       individual_uid !="WAY_Lachemilla orbiculata_1_1"&
                       individual_uid !="WAY_Lachemilla orbiculata_1_12"&
                       individual_uid !="WAY_Gaultheria glomerata_4_1")

#filter out only paspalum that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "WAY_Paspalum bonplandianum_4_1" &
                       individual_uid !="ACJ_Paspalum bonplandianum_3_2" &
                       individual_uid !="ACJ_Paspalum bonplandianum_5_2"&
                       individual_uid !="ACJ_Paspalum bonplandianum_2_1")

#filter out only rhynchospora that has 3 leaves
traits_log <-subset (traits_log, individual_uid != "TRE_Rhynchospora macrochaeta_1_2" &
                       individual_uid !="ACJ_Rhynchospora macrochaeta_4_2")