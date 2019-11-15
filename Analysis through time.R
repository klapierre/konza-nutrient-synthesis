library(vegan)
library(gridExtra)
library(doBy)
library(grid)
library(devtools)
install_github("NCEAS/codyn", ref = "anderson")
library(codyn)
library(tidyverse)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's laptop:
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')

#kim's desktop:
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')



###data
#community data
community<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  #create a replicate variable unique to each experiment
  mutate(replicate=paste(project_name, treatment, plot_id, sep='::'))%>%
  mutate(project_trt=paste(project_name, treatment, sep='::'))%>%
  #filter out 0s
  filter(abundance>0)%>%
  #filter out pre-treatment data
  mutate(pretrt=ifelse(project_name=='pplots'&calendar_year==2002, 1, ifelse(project_name=='nutnet'&calendar_year==2007, 1, 0)))%>%
  filter(pretrt==0)%>%
  select(-pretrt)
  #need to add in ghost fire?

#list of experiments, treatments, and plots
plots <- community%>%
  select(project_name, treatment, plot_id, replicate, calendar_year)%>%
  unique()

#list of experiments, treatments
trt <- community%>%
  select(project_name, treatment, project_trt)%>%
  unique()

#list of experiments
proj <- community%>%
  select(project_name)%>%
  unique()
  



###calculate change metrics
###RAC change
#year to year variation
yearlyRACchange <- RAC_change(community, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', replicate.var='replicate')%>%
  left_join(plots)%>%
  mutate(comparison='yearly')


###compositional change

#makes an empty dataframe
compChange=data.frame(row.names=1) 

###first: composition_change is the change from year to year (yearly) or to the first year of each experiment (cumulative)
####second: dispersion is the average dispersion of plots within a treatment to treatment centriod
for(i in 1:length(trt$project_trt)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- community[community$project_trt==as.character(trt$project_trt[i]),]
  
  #calculating composition change from year to year
  yearly <- multivariate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'abundance', replicate.var = 'plot_id')%>%
    mutate(project_trt=trt$project_trt[i], comparison='yearly')

  #calculating composition change from year to year
  cumulative <- multivariate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'abundance', replicate.var = 'plot_id', reference.time=min(subset$calendar_year))%>%
    mutate(project_trt=trt$project_trt[i], comparison='cumulative')
  
  all <- rbind(yearly, cumulative)
  
  #pasting into the dataframe made for this analysis
  compChange=rbind(all, compChange)  
}

compChange <- compChange%>%
  separate(project_trt, c('project_name', 'treatment'), sep='::')


##figures

# #yearly RAC changes
# ggplot(data=yearlyRACchange, aes(x=calendar_year2, y=richness_change, group=treatment))+
#   geom_point(aes(color=treatment))+
#   geom_line()+
#   geom_hline(yintercept = 0)+
#   #theme(legend.position = "none")+
#   facet_wrap(~project_name, scales="free")

##yearly change
# red= n only
# blue=p only
# purple=n+p
#control=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compChange, project_name=="pplots"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control", "N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("black", "blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)

BGP_ub<-
  ggplot(data=subset(compChange, project_name=="BGP unburned"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P0", "N10 P0", "N0 P1"), values=c("black","red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
  ggplot(data=subset(compChange, project_name=="BGP burned"&treatment!='u_u_n+p'&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control", "N10 P0", "N10 P1", "N0 P1"), values=c("black","red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=0, max=0.8)

nutnet<-
  ggplot(data=subset(compChange, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="K"&treatment!="PK"&treatment!="NK"&treatment!="NPK"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control","N10 P0 K0","N10 P10 K0","N0 P10 K0"), values=c("black", "red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change")+
  ggtitle("Nutrient Network")+
  ylim(min=0, max=0.8)

invert<-
  ggplot(data=subset(compChange, project_name=="invert"&treatment!="x_caged_x"&treatment!="x_caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged_x"&treatment!="x_x_insecticide"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide", 'Control'), values=c("purple","purple","black"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)


grid.arrange(BGP_ub, BGP_b, pplots, nutnet, invert)
#export at 1600x1000



##pretty graphs - cumulative change
# red= n only
# blue=p only
# purple=n+p
#control=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compChange, project_name=="pplots"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control", "N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("black", "blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)

BGP_ub<-
ggplot(data=subset(compChange, project_name=="BGP unburned"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P0", "N10 P0", "N0 P1"), values=c("black","red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
ggplot(data=subset(compChange, project_name=="BGP burned"&treatment!='u_u_n+p'&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control", "N10 P0", "N10 P1", "N0 P1"), values=c("black","red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=0, max=0.8)

nutnet<-
ggplot(data=subset(compChange, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="K"&treatment!="PK"&treatment!="NK"&treatment!="NPK"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("control","N10 P0 K0","N10 P10 K0","N0 P10 K0"), values=c("black", "red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change")+
  ggtitle("Nutrient Network")+
  ylim(min=0, max=0.8)

invert<-
ggplot(data=subset(compChange, project_name=="invert"&treatment!="x_caged_x"&treatment!="x_caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged_x"&treatment!="x_x_insecticide"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide", 'Control'), values=c("purple","purple","black"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)


grid.arrange(BGP_ub, BGP_b, pplots, nutnet, invert)
#export at 1600x1000



#calculate difference in compositional changes
compChangeCtl <- compChange%>%
  filter(treatment=='control'|treatment=='b_u_c'|treatment=='N1P0'|treatment=='u_u_c'|treatment=='x_x_x'|treatment==0)%>%
  rename(composition_change_ctl=composition_change)%>%
  select(-treatment, -dispersion_change)

compChangeDiff <- compChange%>%
  filter(treatment!='control'&treatment!='b_u_c'&treatment!='N1P0'&treatment!='u_u_c'&treatment!='x_x_x'&treatment!=0)%>%
  left_join(compChangeCtl)%>%
  mutate(composition_change_diff=composition_change-composition_change_ctl)


##pretty graphs - cumulative change difference
# red= n only
# blue=p only
# purple=n+p
#control=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compChangeDiff, project_name=="pplots"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change Difference")+
  ggtitle("Phosphorus Plots")+
  ylim(min=-0.1, max=0.45)

BGP_ub<-
  ggplot(data=subset(compChangeDiff, project_name=="BGP unburned"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N0 P1"), values=c("red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change Difference")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=-0.1, max=0.45)

BGP_b<-
  ggplot(data=subset(compChangeDiff, project_name=="BGP burned"&treatment!='u_u_n+p'&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N10 P1", "N0 P1"), values=c("red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change Difference")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=-0.1, max=0.45)

nutnet<-
  ggplot(data=subset(compChangeDiff, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="K"&treatment!="PK"&treatment!="NK"&treatment!="NPK"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0 K0","N10 P10 K0","N0 P10 K0"), values=c("red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change Difference")+
  ggtitle("Nutrient Network")+
  ylim(min=-0.1, max=0.45)

invert<-
  ggplot(data=subset(compChangeDiff, project_name=="invert"&treatment!="x_caged_x"&treatment!="x_caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged_x"&treatment!="x_x_insecticide"&comparison=='cumulative'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide"), values=c("purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Cumulative Community Change Difference")+
  ggtitle("Invertebrate Removal")+
  ylim(min=-0.1, max=0.45)


grid.arrange(BGP_ub, BGP_b, pplots, nutnet, invert)
#export at 1600x1000




##pretty graphs - yearly change difference
# red= n only
# blue=p only
# purple=n+p
#control=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compChangeDiff, project_name=="pplots"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change Difference")+
  ggtitle("Phosphorus Plots")+
  ylim(min=-0.1, max=0.45)

BGP_ub<-
  ggplot(data=subset(compChangeDiff, project_name=="BGP unburned"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N0 P1"), values=c("red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Community Change Difference")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=-0.1, max=0.45)

BGP_b<-
  ggplot(data=subset(compChangeDiff, project_name=="BGP burned"&treatment!='u_u_n+p'&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N10 P1", "N0 P1"), values=c("red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Cumulative Community Change Difference")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=-0.1, max=0.45)

nutnet<-
  ggplot(data=subset(compChangeDiff, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="K"&treatment!="PK"&treatment!="NK"&treatment!="NPK"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0 K0","N10 P10 K0","N0 P10 K0"), values=c("red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Cumulative Community Change Difference")+
  ggtitle("Nutrient Network")+
  ylim(min=-0.1, max=0.45)

invert<-
  ggplot(data=subset(compChangeDiff, project_name=="invert"&treatment!="x_caged_x"&treatment!="x_caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged_x"&treatment!="x_x_insecticide"&comparison=='yearly'), aes(x=as.numeric(calendar_year2), y=composition_change_diff, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide"), values=c("purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Yearly Cumulative Community Change Difference")+
  ggtitle("Invertebrate Removal")+
  ylim(min=-0.1, max=0.45)


grid.arrange(BGP_ub, BGP_b, pplots, nutnet, invert)
#export at 1600x1000




###compositional difference from control plots

#makes an empty dataframe
compDiff=data.frame(row.names=1) 

###first: composition_diff is the distance between trt and controls
####second: dispersion is the average dispersion of plots within a treatment to treatment centriod
for(i in 1:length(proj$project_name)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- community[community$project_name==as.character(proj$project_name[i]),]%>%
    mutate(treatment2=as.character(ifelse(treatment=='control', 'control', ifelse(treatment=='b_u_c', 'control', ifelse(treatment=='N1P0', 'control', ifelse(treatment=='u_u_c', 'control', ifelse(treatment=='x_x_x', 'control', ifelse(treatment==0, 'control', as.character(treatment)))))))))
  
  #calculating composition diffrence from controls
  all <- multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'abundance', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='control')%>%
    mutate(project_name=proj$project_name[i])
  
  #pasting into the dataframe made for this analysis
  compDiff=rbind(all, compDiff)  
}



##figures

##compositional difference
# red= n only
# blue=p only
# purple=n+p
#control=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compDiff, project_name=="pplots"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment22))+
  geom_point(aes(color=treatment22), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)

BGP_ub<-
  ggplot(data=subset(compDiff, project_name=="BGP unburned"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment22))+
  geom_point(aes(color=treatment22), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N0 P1"), values=c("red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
  ggplot(data=subset(compDiff, project_name=="BGP burned"&treatment22!='u_u_n+p'), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment22))+
  geom_point(aes(color=treatment22), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N10 P1", "N0 P1"), values=c("red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=0, max=0.8)

nutnet<-
  ggplot(data=subset(compDiff, project_name=="nutnet"&treatment22!="fence"&treatment22!="NPKfence"&treatment22!="K"&treatment22!="PK"&treatment22!="NK"&treatment22!="NPK"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment22))+
  geom_point(aes(color=treatment22), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0 K0","N10 P10 K0","N0 P10 K0"), values=c("red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Nutrient Network")+
  ylim(min=0, max=0.8)

invert<-
  ggplot(data=subset(compDiff, project_name=="invert"&treatment22!="x_caged_x"&treatment22!="x_caged_insecticide"&treatment22!="NPK_caged_insecticide"&treatment22!="NPK_caged_x"&treatment22!="x_x_insecticide"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment22))+
  geom_point(aes(color=treatment22), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide"), values=c("purple","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)


grid.arrange(BGP_ub, BGP_b, pplots, nutnet, invert)
#export at 1600x1000



#generating figure of compositional change by experiment year with all experiments included in one panel

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

#experiment year information
experimentYear <- read.csv('experiment years.csv')%>%
  mutate(calendar_year2=calendar_year)

compChangeTime <- compChange%>%
  filter(project_name!='ChANGE'&project_name!='ghost fire')%>%
  merge(experimentYear, by=c('project_name', 'calendar_year2'))

#cumulative comp change
compChangeTimeSubset <- compChangeTime%>%
  filter(treatment=='b_u_n'|treatment=='u_u_n'|treatment=='N2P0'|treatment=='N'|treatment=='NPK_x_x')

ggplot(data=subset(compChangeTimeSubset,experiment_year>0&comparison=='cumulative'), aes(x=experiment_year, y=composition_change, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment))) +
  xlab('Experiment Year') +
  ylab('Cumulative Community Change') +
  scale_color_discrete(name='Experiment') +
  # ylim(0,0.8) +
  xlim(0,30)
#export at 1000x800

#comp diff
compDiffTime <- compDiff%>%
  filter(project_name!='ChANGE'&project_name!='ghost fire')%>%
  merge(experimentYear, by=c('project_name', 'calendar_year'))

compDiffTimeSubset <- compDiffTime%>%
  filter(treatment22=='b_u_n'|treatment22=='u_u_n'|treatment22=='N2P0'|treatment22=='N'|treatment22=='NPK_x_x')

ggplot(data=subset(compDiffTimeSubset,experiment_year>0), aes(x=experiment_year, y=composition_diff, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment22))) +
  xlab('Experiment Year') +
  ylab('Community Difference') +
  scale_color_discrete(name='Experiment') +
  # ylim(0,0.8) +
  xlim(0,30)
#export at 1000x800


###comparing yearly change to covariates (herbivores and weather)

#comparing to herbivore outbreaks
herbivores <- read.csv('konza_herbivore_site mean.csv')%>%
  rename(calendar_year=year)%>%
  left_join(compChange)%>%
  filter(comparison=='yearly')%>%
  filter(treatment=='b_u_n'|treatment=='u_u_n'|treatment=='N2P0'|treatment=='N'|treatment=='NPK_x_x')

ggplot(herbivores, aes(x=grasshopper_abund, y=composition_change)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Change') + xlab('Grasshopper Abundance')
#export 1000x800

ggplot(herbivores, aes(x=small_mammal_abund, y=composition_change)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Change') + xlab('Small Mammal Abundance')
#export 1000x800


#comparing to weather data
precip<-read.csv("WETDRY.csv")%>%
  rename(calendar_year=YEAR)%>%
  select(calendar_year, AprSeptPrecip, AnnualPrecip)%>%
  left_join(compChange)%>%
  filter(comparison=='yearly')%>%
  filter(treatment=='b_u_n'|treatment=='u_u_n'|treatment=='N2P0'|treatment=='N'|treatment=='NPK_x_x')

ggplot(precip, aes(x=AnnualPrecip, y=composition_change)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Change') + xlab('Annual Precipitation (mm)')
#export 1000x800

ggplot(precip, aes(x=AprSeptPrecip, y=composition_change)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Change') + xlab('Apr-Sep Precipitation (mm)')
#export 1000x800




###comparing difference to covariates (herbivores, weather)

herbivoresDiff <- read.csv('konza_herbivore_site mean.csv')%>%
  rename(calendar_year=year)%>%
  left_join(compDiff)%>%
  filter(treatment22=='b_u_n'|treatment22=='u_u_n'|treatment22=='N2P0'|treatment22=='N'|treatment22=='NPK_x_x')

ggplot(herbivoresDiff, aes(x=grasshopper_abund, y=composition_diff)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Difference') + xlab('Grasshopper Abundance')
#export 1000x800

ggplot(herbivoresDiff, aes(x=small_mammal_abund, y=composition_diff)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Difference') + xlab('Small Mammal Abundance')
#export 1000x800


#comparing to weather data
precipDiff <- read.csv("WETDRY.csv")%>%
  rename(calendar_year=YEAR)%>%
  select(calendar_year, AprSeptPrecip, AnnualPrecip)%>%
  left_join(compDiff)%>%
  filter(treatment22=='b_u_n'|treatment22=='u_u_n'|treatment22=='N2P0'|treatment22=='N'|treatment22=='NPK_x_x')

ggplot(precipDiff, aes(x=AnnualPrecip, y=composition_diff)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Difference') + xlab('Annual Precipitation (mm)')
#export 1000x800

ggplot(precipDiff, aes(x=AprSeptPrecip, y=composition_diff)) +
  geom_point(size=5) +
  facet_wrap(~project_name) +
  ylab('Community Difference') + xlab('Apr-Sep Precipitation (mm)')
#export 1000x800

