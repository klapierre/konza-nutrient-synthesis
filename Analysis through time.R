library(vegan)
library(gridExtra)
library(doBy)
library(grid)
library(codyn)
library(tidyverse)


#kim's laptop:
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')

#kim's desktop:
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}



theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

`%notin%` <- Negate(`%in%`)

###data
#community data
community<-read.csv("Konza_nutrient synthesis_spp comp_04142020.csv")%>%
  #consistently name control plots
  mutate(treatment2=ifelse(treatment %in% c('x_x_x', 'u_u_c', 'control', 'N1P0', 'b_u_c', '0'), 'control', ifelse(treatment=='A C'&project_name=='GF Burned', 'control', ifelse(treatment=='P C'&project_name=='GF Unburned', 'control', as.character(treatment)))))%>%
  select(-treatment)%>%
  rename(treatment=treatment2)%>%
  #create a replicate variable unique to each experiment
  mutate(replicate=paste(project_name, treatment, plot_id, sep='::'))%>%
  mutate(project_trt=paste(project_name, treatment, sep='::'))%>%
  #filter out 0s
  filter(abundance>0)%>%
<<<<<<< HEAD
  mutate(treatment=as.character(as.factor(treatment)))%>%
  #rename control treatments for consistency across experiments
  mutate(treatment2=ifelse(treatment=="N1P0"|treatment=="b_u_c"|treatment==0|treatment=="control"|treatment=="control_control"|treatment=="u_u_c", 'control', treatment))%>%
  select(-treatment)%>%
  mutate(treatment=treatment2)%>%
  #create project_trt column
  mutate(project_trt=paste(project_name, treatment, sep='::'))
  # #drop ChANGE, ghost fire, and Ukulinga data because time series too short
  # #drop restoration plots because planted communities
  # filter(!(project_name %in% c('ukulinga annual', 'ukulinga four', 'ukulinga unburned', 'restoration', 'ChANGE', 'ghost fire')))
=======
  #filter out sugar addition in Ghost Fire
  filter(treatment %notin% c('A S', 'P S'))%>%
  #filter out pre-treatment data
  mutate(pretrt=ifelse(project_name=='pplots'&calendar_year==2002, 1, ifelse(project_name=='nutnet'&calendar_year==2007, 1, ifelse(project_name=='GF Burned'&calendar_year==2014, 1, ifelse(project_name=='GF Unburned'&calendar_year==2014, 1, ifelse(project_name=='ChANGE'&calendar_year==2013, 1, 0))))))%>%
  filter(pretrt==0)%>%
  select(-pretrt)
>>>>>>> f11d0d389b5b4b3d5d35c73af6f7c958a106761a

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



<<<<<<< HEAD
#multivariate change
multivariateChange <- multivariate_change(community, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', replicate.var='replicate', treatment.var=)

=======
###calculate change metrics
###RAC change
#year to year variation
yearlyRACchange <- RAC_change(community, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', replicate.var='replicate')%>%
  left_join(plots)%>%
  mutate(comparison='yearly')
>>>>>>> f11d0d389b5b4b3d5d35c73af6f7c958a106761a


###compositional change

#makes an empty dataframe
compChange2=data.frame(row.names=1) 

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
  compChange2=rbind(all, compChange2)  
}

compChange <- compChange2%>%
  separate(project_trt, c('project_name', 'treatment'), sep='::')




###calculate difference metrics
###compositional difference

#makes an empty dataframe
compDiff=data.frame(row.names=1) 

for(i in 1:length(proj$project_name)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- community[community$project_name==as.character(proj$project_name[i]),]
  
  #calculating composition difference to control plots
  diff <- multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'abundance', replicate.var = 'plot_id', treatment.var='treatment', reference.treatment='control')%>%
    mutate(project_name=proj$project_name[i])
  
  #pasting into the dataframe made for this analysis
  compDiff=rbind(diff, compDiff)  
}

# write.csv(compDiff, 'Konza_nutrient synthesis_comp difference_05222020.csv')


##figures

##compositional difference
# red= n only
# blue=p only
# purple=n+p
# black=other combo
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(compDiff, project_name=="pplots"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
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
  ggplot(data=subset(compDiff, project_name=="BGP unburned"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0", "N0 P1"), values=c("red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
  ggplot(data=subset(compDiff, project_name=="BGP burned"&treatment2!='u_u_n+p'), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
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
  ggplot(data=subset(compDiff, project_name=="nutnet"&treatment2!="fence"&treatment2!="NPKfence"&treatment2!="K"&treatment2!="PK"&treatment2!="NK"&treatment2!="NPK"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
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
  ggplot(data=subset(compDiff, project_name=="invert"&treatment2!="x_caged_x"&treatment2!="x_caged_insecticide"&treatment2!="NPK_caged_insecticide"&treatment2!="NPK_caged_x"&treatment2!="x_x_insecticide"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "N10 P10 K10\n& Insecticide"), values=c("purple","black"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)

change<-
  ggplot(data=subset(compDiff, project_name=="ChANGE"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",values=c("red","red","red","red","red","red","red"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("ChANGE")+
  ylim(min=0, max=0.8)

GFburned<-
  ggplot(data=subset(compDiff, project_name=="GF Burned"&treatment2 %in% c('A U','P U')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c('N','N + litter'), values=c("red","black"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("GF Burned")+
  ylim(min=0, max=0.8)

GFunburned<-
  ggplot(data=subset(compDiff, project_name=="GF Unburned"&treatment2 %in% c('A U','P U')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c('N - litter','N'), values=c("black","red"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("GF unburned")+
  ylim(min=0, max=0.8)

pushViewport(viewport(layout=grid.layout(3,3)))
print(BGP_b, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(GFburned, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(change, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(pplots, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(nutnet, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(invert, vp=viewport(layout.pos.row=2, layout.pos.col=3))
print(BGP_ub, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(GFunburned, vp=viewport(layout.pos.row=3, layout.pos.col=2))
#export at 1400x800



#generating figure of compositional difference by experiment year with all experiments included in one panel

#experiment year information
experimentYear <- read.csv('experiment years.csv')

compDiffTime <- compDiff%>%
  left_join(experimentYear)

#comp difference across all experiments
compDiffTimeSubset <- compDiffTime%>%
  mutate(keep=ifelse(treatment2 %in% c('b_u_n','u_u_n','N2P0','N','NPK_x_x','10'), 1, ifelse(treatment2=='A U'&project_name=='GF Burned', 1, ifelse(treatment2=='P U'&project_name=='GF Unburned', 1, 0))))%>%
  filter(keep==1)

ggplot(data=subset(compDiffTimeSubset,experiment_year>0), aes(x=experiment_year, y=composition_diff, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment2))) +
  xlab('Experiment Year') +
  ylab('Composition Difference') +
  scale_color_discrete(name='Experiment') +
  # ylim(0,0.8) +
  xlim(0,30)
#export at 1200x800


#calculating change in difference from year x to x-1 for each experiment
compDiffChange <- compDiffTime%>%
  group_by(project_name, treatment2)%>%
  mutate(comp_diff_change=(composition_diff-lag(composition_diff, order_by=experiment_year)))%>%
  ungroup()%>%
  mutate(burn_regime=ifelse(project_name %in% c('BGP burned', 'ChANGE', 'GF Burned'), 'annual', ifelse(project_name %in% c('BGP unburned', 'GF Unburned'), 'unburned', 'two-year')))

# write.csv(compDiffChange, 'Konza_nutrient synthesis_change in yearly diff_05222020.csv')

#burn effect
ggplot(data=subset(compDiffChange, experiment_year<7), aes(x=burn_regime, y=comp_diff_change, fill=as.factor(burned))) +
  geom_boxplot() +
  ylab('Yearly Change in Compositional Difference') + xlab('Burn Regime')
#export at 800x800


###comparing difference to covariates (herbivores, weather)

herbivores <- read.csv('konza_herbivore_site mean.csv')

quantile(herbivores$grasshopper_abund, probs = c(0.05,0.95), na.rm=T) #grasshopper quantiles
quantile(herbivores$small_mammal_abund, probs = c(0.05,0.95), na.rm=T) #small mammal quantiles

compDiffEarly <- compDiffChange%>%
  filter(experiment_year<7)

quantile(compDiffEarly$comp_diff_change, probs = c(0.05,0.95), na.rm=T) #composition difference quantiles

herbivoresDiff <- herbivores%>%
  rename(calendar_year=year)%>%
  left_join(compDiffChange)%>%
  mutate(pre_change=ifelse(experiment_year<7, 1, 0))%>%
  filter(!is.na(comp_diff_change))


ggplot(subset(herbivoresDiff, pre_change==1), aes(x=grasshopper_abund, y=comp_diff_change, color=project_name)) +
  geom_point(size=5) +
  ylab('Yearly Change in Compositional Difference') + xlab('Grasshopper Abundance') +
  xlim(c(35,85)) +
  scale_color_discrete(name="Experiment") +
  theme(legend.position=c(0.1,0.85)) +
  geom_hline(yintercept=-0.08781136, linetype='dashed') + #5th quantile composition difference
  geom_hline(yintercept=0.19865817, linetype='dashed') + #95th quantile composition difference
  geom_vline(xintercept=51.29563, linetype='dashed') + #5th quantile grasshopper abundance
  geom_vline(xintercept=76.40751, linetype='dashed') #95th quantile grasshopper abundance
#export 1200x800
ggplot(subset(herbivoresDiff, pre_change==1), aes(x=grasshopper_abund)) +
  geom_density() +
  ylab('Density') + xlab('Grasshopper Abundance') +
  xlim(c(35,85))
#export 1200x200

ggplot(subset(herbivoresDiff, pre_change==1), aes(x=small_mammal_abund, y=comp_diff_change, color=project_name)) +
  geom_point(size=5) +
  ylab('Yearly Change in Compositional Difference') + xlab('Small Mammal Abundance') +
  xlim(c(5,70)) +
  # scale_color_discrete(name="Experiment") +
  theme(legend.position='none') +
  geom_hline(yintercept=-0.08781136, linetype='dashed') + #5th quantile composition difference
  geom_hline(yintercept=0.19865817, linetype='dashed') + #95th quantile composition difference
  geom_vline(xintercept=9.0500, linetype='dashed') + #5th quantile small mammal abundance
  geom_vline(xintercept=50.8375, linetype='dashed') #95th quantile small mammal abundance
#export 1200x800
ggplot(subset(herbivoresDiff, pre_change==1), aes(x=small_mammal_abund)) +
  geom_density() +
  ylab('Density') + xlab('Small Mammal Abundance') +
  xlim(c(5,70))
#export 1200x200


#comparing to weather data
precip <- read.csv("WETDRY.csv")

quantile(precip$AnnualPrecip, probs = c(0.05,0.95), na.rm=T) #annual precp quantiles
quantile(precip$AprSeptPrecip, probs = c(0.05,0.95), na.rm=T) #growing season precip quantiles

precipDiff <- precip%>%
  rename(calendar_year=YEAR)%>%
  select(calendar_year, AprSeptPrecip, AnnualPrecip)%>%
  left_join(compDiffChange)%>%
  filter(!is.na(comp_diff_change))%>%
  mutate(pre_change=ifelse(experiment_year<7, 1, 0))


ggplot(subset(precipDiff, pre_change==1), aes(x=AnnualPrecip, y=comp_diff_change, color=project_name)) +
  geom_point(size=5) +
  ylab('Yearly Change in Compositional Difference') + xlab('Annual Precipitation (mm)') +
  theme(legend.position='none') +
  xlim(c(500,1200)) +
  geom_hline(yintercept=-0.08781136, linetype='dashed') + #5th quantile composition difference
  geom_hline(yintercept=0.19865817, linetype='dashed') + #95th quantile composition difference
  geom_vline(xintercept=596.50, linetype='dashed') + #5th quantile annual precip
  geom_vline(xintercept=1094.15, linetype='dashed') #95th quantile annual precip
#export 1200x800
ggplot(subset(precipDiff, pre_change==1), aes(x=AnnualPrecip)) +
  geom_density() +
  ylab('Density') + xlab('Annual Precipitation (mm)') +
  xlim(c(500,1200))
#export 1200x200

ggplot(subset(precipDiff, pre_change==1), aes(x=AprSeptPrecip, y=comp_diff_change, color=project_name)) +
  geom_point(size=5) +
  ylab('Yearly Change in Compositional Difference') + xlab('Growing Season Precipitation (mm)') +
  theme(legend.position='none') +
  xlim(c(300,1000)) +
  geom_hline(yintercept=-0.08462066, linetype='dashed') + #5th quantile composition difference
  geom_hline(yintercept=0.19948423, linetype='dashed') + #95th quantile composition difference
  geom_vline(xintercept=400.650, linetype='dashed') + #5th quantile growing season precip
  geom_vline(xintercept=861.175, linetype='dashed') #95th quantile growing season precip
#export 1200x800
ggplot(subset(precipDiff, pre_change==1), aes(x=AprSeptPrecip)) +
  geom_density() +
  ylab('Density') + xlab('Growing Season Precipitation (mm)') +
  xlim(c(300,1000))
#export 1200x200