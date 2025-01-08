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
  select(project_name, plot_id, experiment_year, calendar_year, treatment2, genus_species, relabund) %>% 
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
dataYearNMDS <- metaMDS(dataYear[,6:183])

#gathers NMDS scores
coreScores <- data.frame(scores(dataYearNMDS, display='sites'))

#merge NMDS scores with treatments
coreScoresTrt <- dataYear[,1:5] %>%
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


##### single plot for first and final year, N addition only #####
coreScoresTrtMeanNonly <- coreScoresTrtMean %>% 
  group_by(project_name) %>% 
  mutate(keep=ifelse(experiment_year==min(experiment_year), 'start',
              ifelse(experiment_year==max(experiment_year), 'end', 'middle'))) %>% 
  ungroup() %>% 
  filter(keep=='end')

NMDScombinedFig <- ggplot(data=coreScoresTrtMeanNonly, 
       aes(x=NMDS1_mean, y=NMDS2_mean, color=project_name, group=treatment2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=NMDS2_mean-y_err, ymax=NMDS2_mean+y_err), width=0.02) +
  geom_errorbar(aes(xmin=NMDS1_mean-x_err, xmax=NMDS1_mean+x_err), width=0.02) +
  scale_color_manual(values=c('#f5892a', '#f2cc3a', '#686868', '#39869e', '#54c4b7', '#db4c23'),
                     labels=c('BGP Burned', 'BGP Unburned', 'ChANGE', 'Invert', 'NutNet', 'PPlots')) +
  # geom_text(color='black', size=3) +
  xlab('NMDS1') + ylab('NMDS2') +
  theme(legend.position=c(0.17,0.88), legend.title=element_blank()) +
  # xlim(min=-1, max=2.5) +
  # ylim(min=-1.5, max=1.5) 
  geom_segment(x=0.3107672, xend=2.13234452, y=-0.29681054, yend=0.43714202, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black') + #BGP burned #f5892a
  geom_segment(x=1.5638310, xend=1.80415726, y=0.55707143, yend=0.33578529, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black') + #BGP unburned
  geom_segment(x=-0.1589214, xend=-0.54187134, y=-0.39746504, yend=-0.51749001, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black') + #ChANGE
  geom_segment(x=-0.3142156, xend=0.45248747, y=0.60356896, yend=-0.07327285, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black') + #Invert
  geom_segment(x=-0.4321114, xend=-0.46987265, y=0.62175717, yend=-0.05017150, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black') + #NutNet
  geom_segment(x=-0.2496861, xend=-0.21767465, y=0.10334450, yend=-0.39144485, 
               arrow=arrow(type='closed', length=unit(0.02, "npc")), color='black')  #PPlots
#export at 800x800
# ggsave(NMDScombinedFig, file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig3_NMDS_all_20250108.png', width=8, height=8, units='in', dpi=300, bg='white')
