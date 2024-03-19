################################################################################
##  2_temporal trends and drivers.R: Identifying patterns in community difference
##  through time and possible drivers.
##
##  Authors: Kimberly Komatsu, Meghan Avolio
################################################################################

##### Workspace Set-Up #####

library(vegan)
library(gridExtra)
library(doBy)
library(grid)
library(codyn)
library(ggthemes)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\data') #kim's
setwd("~/Dropbox/Konza Nutrient Synthesis") #meghan's

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.35), axis.text.x=element_text(size=25),
             axis.title.y=element_text(size=30, angle=90, vjust=0.5), axis.text.y=element_text(size=25),
             plot.title = element_text(size=30, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=25))

# bar graph summary statistics function
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

#negation function
`%notin%` <- Negate(`%in%`)



##### Read in Data #####

#community data

community<-read.csv("Konza_nutrient synthesis_spp comp_20240304.csv")%>%
  #consistently name control plots
  mutate(treatment2=ifelse(treatment %in% c('x_x_x', 'u_u_c', 'control', 'N1P0', 'b_u_c', '0'), 'control', ifelse(treatment=='A C'&project_name=='GF Burned', 'control', ifelse(treatment=='P C'&project_name=='GF Unburned', 'control', as.character(treatment)))))%>%
  select(-treatment)%>%
  rename(treatment=treatment2)%>%
  #create a replicate variable unique to each experiment
  mutate(replicate=paste(project_name, treatment, plot_id, sep='::'))%>%
  #filter out 0s
  filter(abundance>0)%>%
  mutate(treatment=as.character(as.factor(treatment)))%>%
  #create project_trt column
  mutate(project_trt=paste(project_name, treatment, sep='::')) %>% 
  # #drop Ukulinga data because time series too short and restoration plots because planted communities
  filter(!(project_name %in% c('ukulinga annual', 'ukulinga four', 'ukulinga unburned', 'restoration'))) %>% 
  #filter out sugar addition in Ghost Fire
  filter(treatment %notin% c('A S', 'P S'))%>%
  #filter out pre-treatment data
  mutate(pretrt=ifelse(project_name=='pplots'&calendar_year==2002, 1, ifelse(project_name=='nutnet'&calendar_year==2007, 1, ifelse(project_name=='GF Burned'&calendar_year==2014, 1, ifelse(project_name=='GF Unburned'&calendar_year==2014, 1, ifelse(project_name=='ChANGE'&calendar_year==2013, 1, 0))))))%>%
  filter(pretrt==0)%>%
  select(-pretrt)

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



##### Calculate Community Change Metrics #####

###RAC change
#year to year variation
yearlyRACchange <- RAC_change(community, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', replicate.var='replicate')%>%
  left_join(plots)%>%
  mutate(comparison='yearly')

### Compositional Change

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



##### Calculate Community Difference Metrics #####

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

# write.csv(compDiff, 'Konza_nutrient synthesis_comp difference_20240305.csv')



##### Compositional Difference Figure -- All Experiments Together #####

# Experiment year information
experimentYear <- read.csv('experiment years.csv')

compDiffTime <- compDiff%>%
  left_join(experimentYear)

# Compositional difference across all experiments
compDiffTimeSubset <- compDiffTime%>%
  mutate(keep=ifelse(treatment2 %in% c('b_u_n','u_u_n','N2P0','N','NPK_x_x','10', '5'), 1, ifelse(treatment2=='A_U'&project_name=='GF Burned', 1, ifelse(treatment2=='P_U'&project_name=='GF Unburned', 1, 0))))%>%
  filter(keep==1)

# Figure by expt year
exptYearFig <- ggplot(data=subset(compDiffTimeSubset,experiment_year>0 & !(project_name %in% c('GF Burned', 'GF Unburned'))), aes(x=experiment_year, y=composition_diff, color=project_name)) +
  geom_rect(aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill="grey", color=NA, alpha=0.01) + #all
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment2)), size=2) +
  xlab('Experiment Year') +
  ylab('Community Difference') +
  scale_color_manual(values=c('#f5892a', '#f2cc3a', 'grey', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  ylim(0,0.8) +
  scale_x_continuous(limits = c(1, 10), breaks = seq(from=2, to=10, by=2)) +
  theme(legend.position='none')
# ggsave(exptYearFig, file='C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 1b_temporalTrajectories_20240306.png', width=7.5, height=7.5, units='in', dpi=300, bg='white')

# Figure by calendar year
calYearFig <- ggplot(data=subset(compDiffTimeSubset,experiment_year>0 & !(project_name %in% c('GF Burned', 'GF Unburned'))), aes(x=calendar_year, y=composition_diff, color=project_name)) +
  geom_rect(aes(xmin = 1988, xmax = 1995, ymin = -Inf, ymax = Inf), fill="#f2cc3a", color=NA, alpha=0.006) + #bgp
  geom_rect(aes(xmin = 2004.5, xmax = 2007.5, ymin = -Inf, ymax = Inf), fill="#db4c23", color=NA, alpha=0.006) + #pplots
  geom_rect(aes(xmin = 2010.5, xmax = 2012.5, ymin = -Inf, ymax = Inf), fill="#54c4b7", color=NA, alpha=0.006) + #nutnet
  geom_rect(aes(xmin = 2011.5, xmax = 2013.5, ymin = -Inf, ymax = Inf), fill="#39869e", color=NA, alpha=0.006) + #invert
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment2)), size=2) +
  xlab('Calendar Year') +
  ylab('Community Difference') +
  scale_color_manual(name='Experiment',
                     values=c('#f5892a', '#f2cc3a', 'grey', '#39869e', '#54c4b7', '#db4c23'), 
                     labels=c("BGP Burned", 'BGP Unburned', 'ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  ylim(0,0.8) 
# ggsave(calYearFig, file='C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 1a_temporalTrajectories_20240306.png', width=17, height=7.5, units='in', dpi=300, bg='white')



##### Compositional Difference Figure -- Interactions with Other Treatments #####

# red= n only
# blue=p only
# purple=n+p
# black=other combo

pplots<-
  ggplot(data=subset(compDiff, project_name=="pplots" & treatment2 %in% c('N1P3', 'N2P0', 'N2P3')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",
                     labels=c("P", "N", "NP"),
                     values=c("blue","red","purple"))+
  geom_line()+
  # geom_vline(xintercept = 2004, linetype="longdash")+
  # geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.88))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)

BGP_ub<-
  ggplot(data=subset(compDiff, project_name=="BGP unburned"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("N", "P"), values=c("red","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.9))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
  ggplot(data=subset(compDiff, project_name=="BGP burned"&treatment2!='u_u_n+p'), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("N", "NP", "P"), values=c("red","purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 1989,linetype="longdash")+
  # geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.88))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=0, max=0.8)

nutnet<-
  ggplot(data=subset(compDiff, project_name=="nutnet"&treatment2!="fence"&treatment2!="NPK_fence"&treatment2!="K"&treatment2!="PK"&treatment2!="NK"&treatment2!="NPK"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("N","NP","P"), values=c("red", "purple","blue"))+
  geom_line()+
  # geom_vline(xintercept = 2010,linetype="longdash")+
  # geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.88))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Nutrient Network")+
  ylim(min=0, max=0.8)

invert<-
  ggplot(data=subset(compDiff, project_name=="invert"&treatment2!="x_x_caged"&treatment2!="x_insect_caged"&treatment2!="NPK_insect_caged"&treatment2!="NPK_x_caged"&treatment2!="x_insect_x"), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",labels=c("NPK", "NPK+I"), values=c("purple","black"))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.9))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)

change<-
  ggplot(data=subset(compDiff, project_name=="ChANGE" & treatment2 %in% c('3', '5')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
  geom_point(aes(color=treatment2), size=5)+
  scale_color_manual(name="Treatment",values=c("pink","red"),
                     labels=c('N5', 'N10'))+
  geom_line()+
  # geom_vline(xintercept = 2011,linetype="longdash")+
  # geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.1,0.9))+
  xlab("Year")+
  ylab("Community Difference")+
  ggtitle("ChANGE")+
  ylim(min=0, max=0.8)

# GFburned<-
#   ggplot(data=subset(compDiff, project_name=="GF Burned"&treatment2 %in% c('A_U','P_U')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
#   geom_point(aes(color=treatment2), size=5)+
#   scale_color_manual(name="Treatment",labels=c('N','N + litter'), values=c("pink","black"))+
#   geom_line()+
#   # geom_vline(xintercept = 2011,linetype="longdash")+
#   # geom_vline(xintercept = 2013)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Difference")+
#   ggtitle("GF Burned")+
#   ylim(min=0, max=0.8)
# 
# GFunburned<-
#   ggplot(data=subset(compDiff, project_name=="GF Unburned"&treatment2 %in% c('A_U','P_U')), aes(x=as.numeric(calendar_year), y=composition_diff, group=treatment2))+
#   geom_point(aes(color=treatment2), size=5)+
#   scale_color_manual(name="Treatment",labels=c('N - litter','N'), values=c("black","pink"))+
#   geom_line()+
#   # geom_vline(xintercept = 2011,linetype="longdash")+
#   # geom_vline(xintercept = 2013)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Difference")+
#   ggtitle("GF Unburned")+
#   ylim(min=0, max=0.8)

pushViewport(viewport(layout=grid.layout(2,3)))
print(BGP_b, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(BGP_ub, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(change, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(pplots, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(nutnet, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(invert, vp=viewport(layout.pos.row=2, layout.pos.col=3))
# print(GFburned, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(GFunburned, vp=viewport(layout.pos.row=3, layout.pos.col=2))
#export at 2800x1400



##### Calculating Change in Difference (from year x to x-1 for each experiment) #####

compDiffChange <- compDiffTime%>%
  group_by(project_name, treatment2)%>%
  mutate(comp_diff_change=(composition_diff-lag(composition_diff, order_by=experiment_year)))%>%
  ungroup()%>%
  mutate(burn_regime=ifelse(project_name %in% c('BGP burned', 'ChANGE', 'GF Burned'), 'annual', ifelse(project_name %in% c('BGP unburned', 'GF Unburned'), 'unburned', 'two-year')))

# write.csv(compDiffChange, 'Konza_nutrient synthesis_change in yearly diff_20240306.csv')

#burn effect
ggplot(data=subset(compDiffChange, experiment_year<7), aes(x=burn_regime, y=comp_diff_change, fill=as.factor(burned))) +
  geom_boxplot() +
  ylab('Yearly Change in Compositional Difference') + xlab('Burn Regime')
#export at 800x800



###### Comparing Difference to Covariates (herbivores, weather) #####

# import drivers data
drivers <- read.csv('Konza_nutrient synthesis_drivers_20240304.csv')

# calculate quantiles of change for early years of focal experiments
compDiffEarly <- compDiffChange%>%
  filter(!(project_name %in% c('GF Burned', 'GF Unburned'))) %>% 
  filter(treatment2 %in% c('5', 'b_u_n', 'u_u_n', 'N', 'N2P0', 'NPK_x_x')) %>% 
  filter(experiment_year<7)

quantile(compDiffEarly$comp_diff_change, probs = c(0.05,0.95), na.rm=T) #composition difference quantiles



# comparing to herbivores

quantile(drivers$grasshopper, probs = c(0.05,0.95), na.rm=T) #grasshopper quantiles
quantile(drivers$mammal, probs = c(0.05,0.95), na.rm=T) #small mammal quantiles

herbivoresDiff <- drivers %>%
  left_join(compDiffEarly, multiple='all') %>%
  filter(!is.na(comp_diff_change))

#models
with(herbivoresDiff, cor.test(grasshopper, comp_diff_change, method = "pearson", use = "complete.obs"))
with(herbivoresDiff, cor.test(mammal, comp_diff_change, method = "pearson", use = "complete.obs"))


#figures

ggplot(herbivoresDiff, aes(x=grasshopper, y=comp_diff_change, color=project_name)) +
  ylab('Yearly Change in\nCompositional Difference') + xlab('Grasshopper Abundance') +
  theme(legend.position='none') +
  scale_color_manual(name='Experiment',
                     values=c('grey', '#39869e', '#54c4b7', '#db4c23'), 
                     labels=c('ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  xlim(c(0,450)) +
  geom_hline(yintercept=-0.1139534  , linetype='dashed', size=2) + #5th quantile composition difference
  geom_hline(yintercept=0.2442750 , linetype='dashed', size=2) + #95th quantile composition difference
  geom_vline(xintercept=38.1500, linetype='dashed', size=2) + #5th quantile grasshopper abundance
  geom_vline(xintercept=411.0875, linetype='dashed', size=2) + #95th quantile grasshopper abundance 
  geom_point(size=10)
#export 1200x800

ggplot(drivers, aes(x=grasshopper)) +
  geom_density(size=2) +
  ylab('Density') + xlab('Grasshopper Abundance') +
  xlim(c(0,450))
#export 1200x250

ggplot(herbivoresDiff, aes(x=mammal, y=comp_diff_change, color=project_name)) +
  ylab('Yearly Change in\nCompositional Difference') + xlab('Small Mammal Abundance') +
  theme(legend.position='none') +
  scale_color_manual(name='Experiment',
                     values=c('grey', '#39869e', '#54c4b7', '#db4c23'), 
                     labels=c('ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  xlim(c(0,35)) +
  geom_hline(yintercept=-0.1139534  , linetype='dashed',size=2) + #5th quantile composition difference
  geom_hline(yintercept=0.2442750 , linetype='dashed',size=2) + #95th quantile composition difference
  geom_vline(xintercept=5.5750, linetype='dashed',size=2) + #5th quantile small mammal abundance
  geom_vline(xintercept=29.0375, linetype='dashed',size=2) + #95th quantile small mammal abundance
  geom_point(size=10)
#export 1200x800

ggplot(drivers, aes(x=mammal)) +
  geom_density(size=2) +
  ylab('Density') + xlab('Small Mammal Abundance') +
  xlim(c(0,35))
#export 1200x250


#comparing to weather data

quantile(drivers$year_precip, probs = c(0.05,0.95), na.rm=T) #annual precp quantiles
quantile(drivers$grow_precip, probs = c(0.05,0.95), na.rm=T) #growing season precip quantiles

precipDiff <- drivers %>%
  left_join(compDiffEarly, multiple='all') %>%
  filter(!is.na(comp_diff_change))

#models
with(precipDiff, cor.test(year_precip, comp_diff_change, method = "pearson", use = "complete.obs"))
with(precipDiff, cor.test(grow_precip, comp_diff_change, method = "pearson", use = "complete.obs"))

#figures

ggplot(precipDiff, aes(x=year_precip, y=comp_diff_change, color=project_name)) +
  ylab('Yearly Change in\nCompositional Difference') + xlab('Annual Precipitation (mm)') +
  theme(legend.position='none') +
  scale_color_manual(name='Experiment',
                     values=c('grey', '#39869e', '#54c4b7', '#db4c23'), 
                     labels=c('ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  xlim(c(400,1350)) +
  geom_hline(yintercept=-0.1139534  , linetype='dashed',size=2) + #5th quantile composition difference
  geom_hline(yintercept=0.2442750 , linetype='dashed',size=2) + #95th quantile composition difference
  geom_vline(xintercept=516, linetype='dashed',size=2) + #5th quantile annual precip
  geom_vline(xintercept=1119, linetype='dashed',size=2) + #95th quantile annual precip 
  geom_point(size=10)
#export 1200x800

ggplot(precipDiff, aes(x=year_precip)) +
  geom_density(size=2) +
  ylab('Density') + xlab('Annual Precipitation (mm)') +
  xlim(c(400,1350))
#export 1200x250

ggplot(precipDiff, aes(x=grow_precip, y=comp_diff_change, color=project_name)) +
  ylab('Yearly Change in\nCompositional Difference') + xlab('Growing Season Precipitation (mm)') +
  scale_color_manual(name='Experiment',
                     values=c('grey', '#39869e', '#54c4b7', '#db4c23'), 
                     labels=c('ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  theme(legend.position='none') +
  xlim(c(300,800)) +
  geom_hline(yintercept=-0.1139534  , linetype='dashed',size=2) + #5th quantile composition difference
  geom_hline(yintercept=0.2442750 , linetype='dashed',size=2) + #95th quantile composition difference
  geom_vline(xintercept=305, linetype='dashed',size=2) + #5th quantile growing season precip
  geom_vline(xintercept=690, linetype='dashed',size=2) + #95th quantile growing season precip
  geom_point(size=10) 
#export 1200x800

ggplot(precipDiff, aes(x=grow_precip)) +
  geom_density(size=2) +
  ylab('Density') + xlab('Growing Season Precipitation (mm)') +
  xlim(c(300,800))
#export 1200x250