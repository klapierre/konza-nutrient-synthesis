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

setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis')

#get data
core <- read.csv('Konza_core data_spp comp_PVC021.csv')

#note -- set first year as 1993 for all watersheds!
coreClean <- core%>%
  #get a species column
  unite(genus_species, AB_GENUS, AB_SPECIES)%>%
  #clean up column names
  mutate(year=RECYEAR, month=RECMONTH, watershed=WATERSHED, soil=SOILTYPE, cover=Cover)%>%
  select(year, month, watershed, soil, genus_species, Transect, Plot, cover)%>%
  #max cover within a year in a plot
  group_by(year, watershed, soil, genus_species, Transect, Plot)%>%
  summarise(cover=max(cover))%>%
  ungroup()%>%
  #mean cover across plots within a transect
  group_by(year, watershed, soil, genus_species, Transect)%>%
  summarise(cover=mean(cover))%>%
  ungroup()%>%
  mutate(plot_id=Transect)%>%
  select(-Transect)%>%
  filter(year>1992)

#get first and last year of data for each watershed
coreFirstLast <- coreClean%>%
  group_by(watershed)%>%
  summarise(max=max(year), min=min(year))%>%
  merge(coreClean, by=c('watershed'))%>%
  mutate(max_year=ifelse(max==year, 1, 0), min_year=ifelse(min==year, 1, 0))%>%
  filter(max_year==1|min_year==1)%>%
  spread(key=genus_species, value=cover, fill=0)


###NMDS for first and last years
coreNMDS <- metaMDS(coreFirstLast[,9:286])

#gathers NMDS scores
coreScores <- data.frame(scores(coreNMDS, display='sites'))

#merge NMDS scores with treatments and filter down to just the watersheds of interest
coreScoresTrt <- coreFirstLast[,1:9]%>%
  cbind(coreScores)%>%
  filter(soil!='s'&watershed!='00fa'&watershed!='00fb'&watershed!='00wa'&watershed!='00wb'&watershed!='0sub'&watershed!='0sua'&watershed!='0spa'&watershed!='0spb'&watershed!='r01a'&watershed!='r01b'&watershed!='r20a'&watershed!='r20b'&watershed!='001a'&watershed!='004f'&watershed!='020a'&watershed!='020d')%>%
  mutate(burn=ifelse(watershed=='001d'|watershed=='n01a', 1, ifelse(watershed=='002c'|watershed=='002d', 2, ifelse(watershed=='004a'|watershed=='004b'|watershed=='n04a'|watershed=='n04d', 4, 20))))%>%
  mutate(graze=ifelse(watershed=='n01a'|watershed=='n01b'|watershed=='n04a'|watershed=='n04d'|watershed=='n20a'|watershed=='n20b', 'grazed', 'ungrazed'))%>%
  unite(trt, graze, burn, remove=F)%>%
  mutate(time=ifelse(max_year==1, 'last', 'first'))

#plot NMDS
ggplot(data=coreScoresTrt, aes(x=NMDS1, y=NMDS2, color=trt, shape=soil)) +
  geom_point(size=5) +
  scale_color_manual(breaks=c('grazed_1', 'grazed_4', 'grazed_20', 'ungrazed_1', 'ungrazed_2', 'ungrazed_4', 'ungrazed_20'),
                    values=c('#FBDF41', '#F37F28', '#EB200F', '#66FC5F', '#49B644', '#2C7029', '#102A0E'),
                    name='Grazing, Burn',
                    labels=c('Grazed, 1yr', 'Grazed, 20yr', 'Grazed, 4yr', 'Ungrazed, 1yr', 'Ungrazed, 2yr', 'Ungrazed, 20yr', 'Ungrazed, 4yr')) +
  scale_shape_discrete(name='Soil Type') +
  facet_wrap(~time)
#export at 1600x1000

  
# #plot NMDS with all years in one panel
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



###Bray-Curtis Dissimilarity from year 1 through time
bc_mean_change=data.frame(watershed=c(), soil=c(), year=c(), mean_change=c()) 

coreTrt <- coreClean%>%
  mutate(treatment=paste(watershed,soil, sep='::'))

coreClean <- coreClean%>%
    mutate(treatment=paste(watershed,soil, sep='::'))

coreTrtList<-unique(coreTrt$treatment)

###first, get bray curtis dissimilarity values for each year within each experiment between all combinations of plots
###second, get distance of each plot within a trt to the trt centroid 
###third: mean_change is the distance between trt and control centriods
for(i in 1:length(coreTrtList)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=coreClean%>%
    filter(treatment==coreTrtList[i], cover!=0)%>%
    select(watershed, soil, treatment, year, genus_species, cover, plot_id)
  
  #get the name of the first year of watershed
  minYear<-subset%>%
    mutate(min=min(year))%>%
    select(min)%>%
    unique
  
  #transpose data
  species=subset%>%
    spread(genus_species, cover, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,6:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$year, type="centroid")
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(treatment=coreTrtList[i],
                      year=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==minYear$min,]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  mc<-cent_C_T%>%
    separate(treatment, c("watershed","soil"), sep="::")
  
  #pasting dispersions into the dataframe made for this analysis
  bc_mean_change=rbind(mc, bc_mean_change)  
}


#filter down to watersheds we are interested in
mean_change <- bc_mean_change%>%
  filter(soil!='s'&watershed!='00fa'&watershed!='00fb'&watershed!='00wa'&watershed!='00wb'&watershed!='0sub'&watershed!='0sua'&watershed!='0spa'&watershed!='0spb'&watershed!='r01a'&watershed!='r01b'&watershed!='r20a'&watershed!='r20b'&watershed!='001a'&watershed!='004f'&watershed!='020a'&watershed!='020d')%>%
  mutate(burn=ifelse(watershed=='001d'|watershed=='n01a', 1, ifelse(watershed=='002c'|watershed=='002d', 2, ifelse(watershed=='004a'|watershed=='004b'|watershed=='n04a'|watershed=='n04d', 4, 20))))%>%
  mutate(graze=ifelse(watershed=='n01a'|watershed=='n01b'|watershed=='n04a'|watershed=='n04d'|watershed=='n20a'|watershed=='n20b', 'grazed', 'ungrazed'))%>%
  unite(trt, graze, burn, remove=F)%>%
  mutate(year=as.numeric(as.character(year)))

#dissimilarity through time (for all watersheds)
ggplot(data=mean_change, aes(x=year, y=mean_change, color=soil)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=16)) +
  scale_x_continuous(breaks=seq(1983, 2016, 3)) +
  ylab('Community Change') +
  xlab('Year') +
  facet_wrap(~watershed)

#dissimilarity through time (avg response for trts)
ggplot(data=barGraphStats(data=mean_change, variable="mean_change", byFactorNames=c("trt", "graze", "burn", "year", "soil")), aes(x=year, y=mean, color=soil)) +
  geom_point() +
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  geom_line() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=16)) +
  scale_x_continuous(breaks=seq(1983, 2016, 3)) +
  ylab('Community Change') +
  xlab('Year') +
  scale_color_manual(values=c('#FF3300', '#330033')) +
  facet_grid(graze~burn)
#export at 1500x1000
  