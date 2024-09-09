################################################################################
##  4_ordinations.R: Identifying patterns community responses
##  through time and possible drivers.
##
##  Authors: Kimberly Komatsu
################################################################################

##### Workspace Set-Up #####

library(vegan)
library(grid)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

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

setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\data') #kim's

##### Import data #####

data <- read.csv("Konza_nutrient synthesis_spp comp_20240304.csv") %>% 
  filter(!(genus_species %in% c('litter', 'litter ', 'bare_ground ', 'NA_NA'))) %>% 
  mutate(treatment=ifelse(treatment=='N1P0', 'control', treatment))

exyrs <- read.csv("experiment years.csv")


##### Calculate Relative Abundances #####

totAbundundance <- data %>%
  group_by(project_name, plot_id, calendar_year) %>%
  summarise(tot=sum(abundance)) %>% 
  ungroup()

relabund <- data %>%
  left_join(totAbundundance) %>%
  mutate(relabund=abundance/tot) %>%
  filter(relabund!="0")


##### Determine experiment start date and number of plots #####

yrs <- data %>%
  select(project_name, calendar_year) %>%
  unique %>%
  group_by(project_name) %>%
  summarise(start_year=min(calendar_year)) %>%
  ungroup() %>%
  mutate(pretreatment=ifelse(project_name %in% c('ChANGE', 'pplots', 'nutnet'), 1, 0)) %>% 
  mutate(start_year=start_year+pretreatment) %>% 
  select(project_name, start_year)


##### Merge cover and year data #####
dataYear <- relabund %>% 
  left_join(yrs) %>% 
  mutate(experiment_year=calendar_year-start_year+1) %>% 
  filter(experiment_year>0,
         project_name %in% c('pplots', 'nutnet', 'BGP burned', 'BGP unburned', 'ChANGE', 'invert'),
         treatment %in% c('control', 'NPK_x_x', 'b_u_n', 'u_u_n', 'N', 'N2P0', '5')) %>% 
  mutate(treatment2=ifelse(treatment=='control', 'control', 'N')) %>% 
  select(project_name, plot_id, experiment_year, treatment2, genus_species, relabund) %>% 
  mutate(genus_species=ifelse(genus_species=='andropogon_scoparius', 'schizachyrium_scoparium',
                       ifelse(genus_species=='kuhnia_eupatorioides', 'brickellia_eupatorioides',
                       ifelse(genus_species=='carex_heliophila', 'carex_inops',
                       ifelse(genus_species=='sporobolus_asper', 'sporobolus_compositus',
                       ifelse(genus_species=='aster_ericoides', 'symphyotrichum_ericoides',
                       ifelse(genus_species %in% c('aster_oblongifolia', 'aster_oblongifolius'), 'symphyotrichum_oblongifolium',
                       ifelse(genus_species=='solidago_canadensis', 'solidago_altissima',
                       genus_species)))))))) %>% 
  mutate(genus_species=str_to_lower(genus_species)) %>% 
  mutate(genus_species=str_replace(genus_species, " ", "_")) %>% 
  pivot_wider(names_from=genus_species, values_from=relabund, values_fill=0)



##### NMDS #####
dataYearNMDS <- metaMDS(dataYear[,5:183])

#gathers NMDS scores
coreScores <- data.frame(scores(dataYearNMDS, display='sites'))

#merge NMDS scores with treatments
coreScoresTrt <- dataYear[,1:4] %>%
  cbind(coreScores)

coreScoresTrtMean <- coreScoresTrt %>% 
  group_by(project_name, treatment2, experiment_year) %>% 
  summarize(NMDS1_mean=mean(NMDS1), x_err=(sd(NMDS1)/sqrt(length(NMDS1))), NMDS2_mean=mean(NMDS2), y_err=(sd(NMDS2)/sqrt(length(NMDS2)))) %>% 
  ungroup()

#plot NMDS

pplots<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='pplots'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                     label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('NMDS1') + ylab('NMDS2') +
  theme(legend.position='none') +
  ggtitle("(d) Phosphorus Plots") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)

BGP_ub<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='BGP unburned'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                     label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('') + ylab('') +
  theme(legend.position='none') +
  ggtitle("(b) Belowground Plots Unburned") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)

BGP_b<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='BGP burned'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                           label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('') + ylab('NMDS2') +
  theme(legend.position=c(0.2,0.8), legend.title=element_blank()) +
  ggtitle("(a) Belowground Plots Burned") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)

nutnet<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='nutnet'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                         label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('NMDS1') + ylab('') +
  theme(legend.position='none') +
  ggtitle("(e) Nutrient Network") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)

invert<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='invert'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                         label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('NMDS1') + ylab('') +
  theme(legend.position='none') +
  ggtitle("(f) Invertebrate Removal") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)

change<-
  ggplot(data=subset(coreScoresTrtMean, project_name=='ChANGE'), aes(x=NMDS1_mean, y=NMDS2_mean, color=treatment2, group=treatment2, 
                                                                         label=experiment_year)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err)) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err)) +
  scale_color_manual(values=c('black', 'red3'), labels=c('control', 'N')) +
  geom_text(color='white', size=3) +
  xlab('') + ylab('') +
  theme(legend.position='none') +
  ggtitle("(c) ChANGE") +
  xlim(min=-1, max=2.5) +
  ylim(min=-1.5, max=1.5)


pushViewport(viewport(layout=grid.layout(2,3)))
print(BGP_b, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(BGP_ub, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(change, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(pplots, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(nutnet, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(invert, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 1600x800


