################################################################################
##  3_species responses.R: Identifying patterns in individual species responses
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
library(RColorBrewer)
library(ggthemes)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\data') #kim's
# setwd("~/Dropbox/Konza Nutrient Synthesis") #meghan's

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.35), axis.text.x=element_text(size=25),
             axis.title.y=element_text(size=30, angle=90, vjust=0.5), axis.text.y=element_text(size=25),
             plot.title = element_text(size=30, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=25),
             strip.text.x=element_text(size=25), strip.text.y=element_text(size=25))

# bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count]<- byFactorNames
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

data <- read.csv("Konza_nutrient synthesis_spp comp_20240304.csv") %>% 
  filter(genus_species %notin% c('litter', 'litter '))

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

nplots <- data %>%
  mutate(treatment2=ifelse(treatment %in% c('x_x_x', 'b_u_c', 'u_u_c', 'control', 'N1P0', '0'), 'control',
                    ifelse(treatment %in% c('NPK_x_x', 'b_u_n', 'u_u_n', 'N', 'N2P0', '5'), 'nitrogen', 'ignore'))) %>%
  filter(treatment2=="control") %>%
  select(project_name, plot_id) %>%
  unique() %>%
  group_by(project_name) %>%
  summarize(nplot=length(plot_id)) %>% 
  ungroup()


##### Calculate DCi and rank #####

domsp <- relabund %>%
  mutate(treatment2=ifelse(treatment %in% c('x_x_x', 'b_u_c', 'u_u_c', 'control', 'N1P0', '0'), 'control',
                    ifelse(treatment %in% c('NPK_x_x', 'b_u_n', 'u_u_n', 'N', 'N2P0', '5'), 'nitrogen', 'ignore'))) %>%
  filter(treatment2!="ignore") %>% 
  select(-treatment) %>% 
  rename(treatment=treatment2) %>% 
  left_join(yrs) %>% 
  mutate(experiment_year=calendar_year-start_year+1) %>%  
  group_by(project_name, experiment_year, treatment, genus_species) %>%
  summarize(mabund=mean(relabund), n=length(relabund)) %>%
  ungroup() %>% 
  left_join(nplots) %>%
  group_by(project_name, experiment_year, treatment) %>% 
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi)) %>%
  ungroup() %>% 
  mutate(genus_species=str_to_lower(genus_species)) %>% 
  separate(genus_species, into=c('genus', 'species1'), sep=' ') %>% 
  separate (genus, into=c('genus2', 'species2'), sep='_') %>% 
  mutate(species=ifelse(is.na(species2), species1, species2),
         genus_species=paste(genus2, species, sep='_')) %>% 
  select(-genus2, -species1, -species2, -species) %>% 
  mutate(genus_species=ifelse(genus_species=='andropogon_scoparius', 'schizachyrium_scoparium',
                       ifelse(genus_species=='kuhnia_eupatorioides', 'brickellia_eupatorioides',
                       ifelse(genus_species=='carex_heliophila', 'carex_inops',
                       ifelse(genus_species=='sporobolus_asper', 'sporobolus_compositus',
                       ifelse(genus_species=='aster_ericoides', 'symphyotrichum_ericoides',
                       ifelse(genus_species %in% c('aster_oblongifolia', 'aster_oblongifolius'), 'symphyotrichum_oblongifolium',
                       ifelse(genus_species=='solidago_canadensis', 'solidago_altissima',
                              genus_species)))))))) %>% 
  mutate(genus_species2=ifelse(genus_species=='ambrosia_psilostachya', 'Ambrosia psilostachya',
                        ifelse(genus_species=='amorpha_canescens', 'Amorpha canescens',
                        ifelse(genus_species=='andropogon_gerardii', 'Andropogon gerardii',
                        ifelse(genus_species=='bouteloua_curtipendula', 'Bouteloua curtipendula',
                        ifelse(genus_species=='dichanthelium_oligosanthes', 'Dichanthelium oligosanthes',
                        ifelse(genus_species=='panicum_virgatum', 'Panicum virgatum',
                        ifelse(genus_species=='physalis_pumila', 'Physalis pumila',
                        ifelse(genus_species=='schizachyrium_scoparium', 'Schizachyrium scoparium',
                        ifelse(genus_species=='sorghastrum_nutans', 'Sorghastrum nutans', genus_species)))))))))) %>% 
  mutate(project2=ifelse(project_name=='invert', 'Invert Removal',
                  ifelse(project_name=='nutnet', 'NutNet',
                  ifelse(project_name=='pplots', 'PPlots', project_name))))

domsp$genus_species2 = factor(domsp$genus_species2, levels=c('Schizachyrium scoparium', 'Sorghastrum nutans', 'Andropogon gerardii',
                                                             'Bouteloua curtipendula', 'Panicum virgatum', 'Dichanthelium oligosanthes', 
                                                             'Ambrosia psilostachya', 'Amorpha canescens', 'Physalis pumila', 'NA'))


#DCi figure
ggplot(data=subset(domsp, experiment_year<7 & experiment_year!=0 & treatment=='nitrogen' &
                          !(project_name %in% c('GF Burned', 'GF Unburned', 'ukulinga_annual', 'ukulinga_four', 'ukulinga_unburned', 'ChANGE')) &
                          genus_species %in% c('ambrosia_psilostachya', 'amorpha_canescens', 'andropogon_gerardii',
                                               'bouteloua_curtipendula', 'dichanthelium_oligosanthes', 'panicum_virgatum', 'physalis_pumila',
                                               'schizachyrium_scoparium', 'sorghastrum_nutans')), 
       aes(x=experiment_year, y=DCi, color=project_name)) +
  geom_point(size=5) +
  geom_line(linewidth=2) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(from=1, to=6, by=1)) +
  xlab('Experiment Year') +
  ylab('DCi') +
  scale_color_manual(values=c('#f5892a', '#f2cc3a', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  facet_wrap(~genus_species2, scales='free_y') +
  theme(legend.position='none')


#rel abundance figure
ggplot(data=subset(domsp, experiment_year<7 & experiment_year!=0 & treatment=='nitrogen' &
                     !(project_name %in% c('GF Burned', 'GF Unburned', 'ukulinga_annual', 'ukulinga_four', 'ukulinga_unburned', 'ChANGE')) &
                     genus_species %in% c('ambrosia_psilostachya', 'amorpha_canescens', 'andropogon_gerardii',
                                          'bouteloua_curtipendula', 'dichanthelium_oligosanthes', 'panicum_virgatum', 'physalis_pumila',
                                          'schizachyrium_scoparium', 'sorghastrum_nutans')), 
       aes(x=experiment_year, y=100*mabund, color=project_name)) +
  geom_point(size=5) +
  geom_line(linewidth=2) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(from=1, to=6, by=1)) +
  xlab('Experiment Year') +
  ylab('Relative Abundance') +
  scale_color_manual(values=c('#f5892a', '#f2cc3a', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  facet_wrap(~genus_species2, scales='free_y') +
  theme(legend.position='none')
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 5_sppTrajectories_20240703.png', width=17, height=17, units='in', dpi=300, bg='white')






#absolute abundance figures

abund <- relabund %>%
  mutate(treatment2=ifelse(treatment %in% c('x_x_x', 'b_u_c', 'u_u_c', 'control', 'N1P0', '0'), 'control',
                           ifelse(treatment %in% c('NPK_x_x', 'b_u_n', 'u_u_n', 'N', 'N2P0', '5'), 'nitrogen', 'ignore'))) %>%
  filter(treatment2!="ignore") %>% 
  select(-treatment) %>% 
  rename(treatment=treatment2) %>% 
  left_join(yrs) %>% 
  mutate(experiment_year=calendar_year-start_year+1) %>%  
  group_by(project_name, experiment_year, treatment, genus_species) %>%
  summarize(mabund=mean(cov), n=length(cov)) %>%
  ungroup() %>% 
  mutate(genus_species=str_to_lower(genus_species)) %>% 
  separate(genus_species, into=c('genus', 'species1'), sep=' ') %>% 
  separate (genus, into=c('genus2', 'species2'), sep='_') %>% 
  mutate(species=ifelse(is.na(species2), species1, species2),
         genus_species=paste(genus2, species, sep='_')) %>% 
  select(-genus2, -species1, -species2, -species) %>% 
  mutate(genus_species=ifelse(genus_species=='andropogon_scoparius', 'schizachyrium_scoparium',
                       ifelse(genus_species=='kuhnia_eupatorioides', 'brickellia_eupatorioides',
                       ifelse(genus_species=='carex_heliophila', 'carex_inops',
                       ifelse(genus_species=='sporobolus_asper', 'sporobolus_compositus',
                       ifelse(genus_species=='aster_ericoides', 'symphyotrichum_ericoides',
                       ifelse(genus_species %in% c('aster_oblongifolia', 'aster_oblongifolius'), 'symphyotrichum_oblongifolium',
                       ifelse(genus_species=='solidago_canadensis', 'solidago_altissima',
                       genus_species)))))))) %>% 
  mutate(genus_species2=ifelse(genus_species=='ambrosia_psilostachya', 'Ambrosia psilostachya',
                        ifelse(genus_species=='amorpha_canescens', 'Amorpha canescens',
                        ifelse(genus_species=='andropogon_gerardii', 'Andropogon gerardii',
                        ifelse(genus_species=='bouteloua_curtipendula', 'Bouteloua curtipendula',
                        ifelse(genus_species=='dichanthelium_oligosanthes', 'Dichanthelium oligosanthes',
                        ifelse(genus_species=='panicum_virgatum', 'Panicum virgatum',
                        ifelse(genus_species=='physalis_pumila', 'Physalis pumila',
                        ifelse(genus_species=='schizachyrium_scoparium', 'Schizachyrium scoparium',
                        ifelse(genus_species=='sorghastrum_nutans', 'Sorghastrum nutans', genus_species)))))))))) %>% 
  mutate(project2=ifelse(project_name=='invert', 'Invert Removal',
                  ifelse(project_name=='nutnet', 'NutNet',
                  ifelse(project_name=='pplots', 'PPlots', project_name))))

abund$genus_species2 = factor(abund$genus_species2, levels=c('Schizachyrium scoparium', 'Sorghastrum nutans', 'Andropogon gerardii',
                                                             'Bouteloua curtipendula', 'Panicum virgatum', 'Dichanthelium oligosanthes', 
                                                             'Ambrosia psilostachya', 'Amorpha canescens', 'Physalis pumila', 'NA'))

ggplot(data=subset(abund, experiment_year<7 & experiment_year!=0 & treatment=='nitrogen' &
                     !(project_name %in% c('GF Burned', 'GF Unburned', 'ukulinga_annual', 'ukulinga_four', 'ukulinga_unburned', 'ChANGE')) &
                     genus_species %in% c('ambrosia_psilostachya', 'amorpha_canescens', 'andropogon_gerardii',
                                          'bouteloua_curtipendula', 'dichanthelium_oligosanthes', 'panicum_virgatum',
                                          'schizachyrium_scoparium', 'sorghastrum_nutans')), 
       aes(x=experiment_year, y=mabund, color=project_name)) +
  geom_point(size=5) +
  geom_line(linewidth=2) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(from=1, to=6, by=1)) +
  xlab('Experiment Year') +
  ylab('Relative Abundance') +
  scale_color_manual(values=c('#f5892a', '#f2cc3a', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  facet_wrap(~genus_species2, scales='free_y') +
  theme(legend.position='none')
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 5_sppTrajectoriesAbsAbund_20240703.png', width=17, height=17, units='in', dpi=300, bg='white')



# abs abundance by expt
abund$project_name = factor(abund$project_name, levels=c('BGP burned', 'BGP unburned', 'ChANGE', 'pplots', 'nutnet', 'invert'))


ggplot(data=subset(abund, experiment_year<7 & experiment_year!=0 & treatment=='nitrogen' &
                     !(project2 %in% c('GF Burned', 'GF Unburned', 'ukulinga_annual', 'ukulinga_four', 'ukulinga_unburned')) &
                     genus_species %in% c('ambrosia_psilostachya', 'amorpha_canescens', 'andropogon_gerardii',
                                          'bouteloua_curtipendula', 'dichanthelium_oligosanthes', 'panicum_virgatum',
                                          'schizachyrium_scoparium', 'sorghastrum_nutans')), 
       aes(x=experiment_year, y=mabund, color=genus_species2)) +
  geom_point(size=5) +
  geom_line(linewidth=2) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(from=1, to=6, by=1)) +
  xlab('Experiment Year') +
  ylab('Percent Abundance') +
  scale_color_brewer(palette = "Paired") +
  # scale_color_manual(values=c('#f5892a', '#f2cc3a', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  facet_wrap(~project2, scales='free_y') +
  theme(legend.position='bottom')
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 5_sppTrajectoriesByExptAbsAbund_20240703.png', width=20, height=15, units='in', dpi=300, bg='white')


# relative abundance by expt
domsp$project_name = factor(domsp$project_name, levels=c('BGP burned', 'BGP unburned', 'ChANGE', 'pplots', 'nutnet', 'invert'))


ggplot(data=subset(domsp, experiment_year<7 & experiment_year!=0 & treatment=='nitrogen' &
                     !(project2 %in% c('GF Burned', 'GF Unburned', 'ukulinga_annual', 'ukulinga_four', 'ukulinga_unburned')) &
                     genus_species %in% c('ambrosia_psilostachya', 'amorpha_canescens', 'andropogon_gerardii',
                                          'bouteloua_curtipendula', 'dichanthelium_oligosanthes', 'panicum_virgatum',
                                          'schizachyrium_scoparium', 'sorghastrum_nutans')), 
       aes(x=experiment_year, y=100*mabund, color=genus_species2)) +
  geom_point(size=5) +
  geom_line(linewidth=2) +
  scale_x_continuous(limits = c(1, 6), breaks = seq(from=1, to=6, by=1)) +
  xlab('Experiment Year') +
  ylab('Relative Abundance (%)') +
  scale_color_brewer(palette = "Paired") +
  # scale_color_manual(values=c('#f5892a', '#f2cc3a', '#39869e', '#54c4b7', '#db4c23'), name=element_blank()) +
  facet_wrap(~project2, scales='free_y') +
  theme(legend.position='bottom')
# ggsave(file='C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\figures\\Fig 5_sppTrajectoriesByExptRelAbund_20240703.png', width=20, height=15, units='in', dpi=300, bg='white')












# domsp_change <- relabund %>%
#   right_join(trt_CN) %>%
#   right_join(n_domsp) %>%
#   filter(trt!="C") %>%
#   group_by(project_name, genus_species, calendar_year) %>%
#   summarize(mabund=mean(relabund)) %>%
#   ungroup()
# 
# ggplot(data=domsp_change, aes(x=calendar_year, y=mabund, group=sp, color=sp))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~project_name, scales="free")
# 
# 
# N_yr7 <- relabund %>%
#   right_join(trt_CN) %>%
#   right_join(exyrs) %>%
#   filter(trt=="N") %>%
#   filter(experiment_year==7) %>%
#   group_by(project_name, genus_species) %>%
#   summarize(mabund=mean(relabund), n=length(relabund)) %>%
#   left_join(nplots) %>%
#   mutate(freq=n/nplot,
#          DCi=(mabund*freq)/2,
#          rank=rank(-DCi)) %>%
#   filter(rank<4)
#   
# N_yr14 <- relabund %>%
#   right_join(trt_CN) %>%
#   right_join(exyrs) %>%
#   filter(trt=="N") %>%
#   filter(experiment_year==14) %>%
#   group_by(project_name, genus_species) %>%
#   summarize(mabund=mean(relabund), n=length(relabund)) %>%
#   left_join(nplots) %>%
#   mutate(freq=n/nplot,
#          DCi=(mabund*freq)/2,
#          rank=rank(-DCi)) %>%
#   filter(rank<4)
# 
# N_yr9 <- relabund %>%
#   right_join(trt_CN) %>%
#   right_join(exyrs) %>%
#   filter(trt=="N") %>%
#   filter(experiment_year==9) %>%
#   group_by(project_name, genus_species) %>%
#   summarize(mabund=mean(relabund), n=length(relabund)) %>%
#   left_join(nplots) %>%
#   mutate(freq=n/nplot,
#          DCi=(mabund*freq)/2,
#          rank=rank(-DCi)) %>%
#   filter(rank<4)
# 
# lyrs <- data %>%
#   select(project_name, calendar_year) %>%
#   unique %>%
#   group_by(project_name) %>%
#   summarise(end=max(calendar_year)) %>%
#   ungroup() %>%
#   select(project_name, end) %>%
#   rename(calendar_year=end)
# 
# N_yrlast <- relabund %>%
#   right_join(trt_CN) %>%
#   right_join(lyrs) %>%
#   filter(trt=="N") %>%
#   group_by(project_name, genus_species) %>%
#   summarize(mabund=mean(relabund), n=length(relabund)) %>%
#   left_join(nplots) %>%
#   mutate(freq=n/nplot,
#          DCi=(mabund*freq)/2,
#          rank=rank(-DCi)) %>%
#   filter(rank<4)