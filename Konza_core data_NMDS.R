library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis')

core <- read.csv('Konza_core data_spp comp_PVC021.csv')

coreClean <- core%>%
  unite(genus_species, AB_GENUS, AB_SPECIES)%>%
  mutate(year=RECYEAR, month=RECMONTH, watershed=WATERSHED, soil=SOILTYPE, cover=Cover)%>%
  select(year, month, watershed, soil, genus_species, Transect, Plot, cover)%>%
  group_by(year, watershed, soil, genus_species, Transect)%>%
  summarise(cover=mean(cover))%>%
  ungroup()%>%
  mutate(plot_id=Transect)%>%
  select(-Transect)%>%
  group_by(year, watershed, soil, genus_species, plot_id)%>%
  summarise(cover=max(cover))%>%
  ungroup()

coreFirstLast <- coreClean%>%
  group_by(watershed)%>%
  summarise(max=max(year), min=min(year))%>%
  merge(coreClean, by=c('watershed'))%>%
  mutate(max_year=ifelse(max==year, 1, 0), min_year=ifelse(min==year, 1, 0))%>%
  filter(max_year==1|min_year==1)%>%
  spread(key=genus_species, value=cover, fill=0)


coreNMDS <- metaMDS(coreFirstLast[,9:288])

coreScores <- data.frame(scores(coreNMDS, display='sites'))

coreScoresTrt <- coreFirstLast[,1:9]%>%
  cbind(coreScores)%>%
  filter(soil!='s'&watershed!='00fa'&watershed!='00fb'&watershed!='00wa'&watershed!='00wb'&watershed!='0sub'&watershed!='0sua'&watershed!='0spa'&watershed!='0spb'&watershed!='r01a'&watershed!='r01b'&watershed!='r20a'&watershed!='r20b')%>%
  mutate(burn=ifelse(watershed=='001d'|watershed=='n01a', 1, ifelse(watershed=='002c'|watershed=='002d', 2, ifelse(watershed=='004a'|watershed=='004b'|watershed=='n04a'|watershed=='n04d', 4, 20))))%>%
  mutate(graze=ifelse(watershed=='n01a'|watershed=='n01b'|watershed=='n04a'|watershed=='n04d'|watershed=='n20a'|watershed=='n20b', 'grazed', 'ungrazed'))%>%
  unite(trt, graze, burn, remove=F)%>%
  mutate(time=ifelse(max_year==1, 'last', 'first'))


ggplot(data=coreScoresTrt, aes(x=NMDS1, y=NMDS2, color=trt, shape=soil)) +
  geom_point(size=5) +
  scale_color_manual(breaks=c('grazed_1', 'grazed_4', 'grazed_20', 'ungrazed_1', 'ungrazed_2', 'ungrazed_4', 'ungrazed_20'),
                    values=c('#FBDF41', '#F37F28', '#EB200F', '#66FC5F', '#49B644', '#2C7029', '#102A0E'),
                    name='Grazing, Burn',
                    labels=c('Grazed, 1yr', 'Grazed, 20yr', 'Grazed, 4yr', 'Ungrazed, 1yr', 'Ungrazed, 2yr', 'Ungrazed, 20yr', 'Ungrazed, 4yr')) +
  scale_shape_discrete(name='Soil Type') +
  facet_wrap(~time)

  
# coreScoresMean <- coreScoresTrt%>%
#   group_by(time, trt, soil)%>%
#   summarise(NMDS1_mean=mean(NMDS1), NMDS2_mean=mean(NMDS2))
# 
# ggplot(data=coreScoresMean, aes(x=NMDS1_mean, y=NMDS2_mean, color=trt, shape=interaction(soil,time))) +
#   geom_point(size=5) +
#   scale_color_manual(breaks=c('grazed_1', 'grazed_4', 'grazed_20', 'ungrazed_1', 'ungrazed_2', 'ungrazed_4', 'ungrazed_20'),
#                      values=c('#FBDF41', '#F37F28', '#EB200F', '#66FC5F', '#49B644', '#2C7029', '#102A0E'),
#                      name='Grazing, Burn',
#                      labels=c('Grazed, 1yr', 'Grazed, 20yr', 'Grazed, 4yr', 'Ungrazed, 1yr', 'Ungrazed, 2yr', 'Ungrazed, 20yr', 'Ungrazed, 4yr')) +
#   scale_shape_manual(name='Soil Type, Year',
#                      values=c(2,1,17,16))




  
  