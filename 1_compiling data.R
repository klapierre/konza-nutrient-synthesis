################################################################################
##  1_compiling data.R: Gathering community and abundance data across experiments
##  and Konza long-term observational records.
##
##  Authors: Kimberly Komatsu, Meghan Avolio
################################################################################

##### Workspace Set-Up #####

library(tidyverse)

### Setting working directories
setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Konza Nutrient Synthesis\\data') #kim's
setwd("~/Dropbox/Konza Nutrient Synthesis") #meghan's

### Functions
#not in function
`%notin%` <- Negate(`%in%`)



##### Read in Data #####

#konza spp lists
spp <- read.csv('species_list\\PPS011_new KNZ spp list.csv')

#treatments - for experiments where treatments designations are not in the species composition file
bgp_trt <- read.csv("treatments\\belowground_plots_anpp_1989-2015.csv")%>%
  filter(MOW!="m")%>%
  mutate(treat_other_name=paste(BURN, MOW, NUTRIENT, sep="_"),
         plot_id=PLOT)%>%
  filter(treat_other_name!="__")%>%
  mutate(project_name=ifelse(treat_other_name %in% c('u_u_n','u_u_c','u_u_p','u_u_b'), 'BGP unburned', 'BGP burned'),
         treatment=ifelse(treat_other_name %in% c('u_u_c','b_u_c'), 'control', treat_other_name))%>%
  select(project_name,plot_id,treatment)%>%
  unique()

change_trt <- read.csv("treatments\\ChANGE_treatments.csv")%>%
  mutate(plot_id=plot, treatment=ifelse(N==0, 'control', as.factor(N)), 
         project_name='ChANGE')%>%
  select(project_name,plot_id, treatment)

nutnet_trt <- read.csv('treatments\\KNZ_NutNet_trt.csv')%>%
  mutate(project_name='nutnet')%>%
  rename(treatment=treat_other_name, plot_id=plot)%>%
  select(project_name, plot_id, treatment)

invert_trt <- read.csv('treatments\\invert_trt.csv')%>%
  mutate(project_name='invert',
         treat_other_name=paste(NPK,insecticide,exclose,sep='_'))%>%
  rename(plot_id=plot)%>%
  mutate(treatment=ifelse(treat_other_name=='x_x_x', 'control', treat_other_name))%>%
  select(project_name, plot_id, treatment)

ghostfire_trt <- read.csv('treatments\\ghost fire_trt.csv')%>%
  mutate(project_name=ifelse(Burn.Trt=='Annual', 'GF Burned', 'GF Unburned'), 
         plot_id=Plot,
         treat_other_name=paste(Litter,Nutrient, sep='_'),
         treatment=ifelse(treat_other_name=='A_C'&project_name=='GF Burned', 'control', ifelse(treat_other_name=='P_C'&project_name=='GF Unburned', 'control', treat_other_name)))%>%
  select(project_name,plot_id, treatment)

ukulinga_trt <- read.csv('treatments\\ukulinga_trt.csv')%>%
  mutate(project_name=ifelse(site=='1D', 'ukulinga_annual', ifelse(site=='4F', 'ukulinga_four', 'ukulinga_unburned')),
         treatment=ifelse(fert==1, 'control', 'N'))

#species comp data with species numbers
sp_invert <- read.csv('species_comp\\VIR011_invertRemoval.csv')%>%
  mutate(project_name='invert')%>%
  rename(calendar_year=RecYear, plot_id=Plot)%>%
  mutate(code=ifelse(SpeCode %in% c(900, 957), 517, ifelse(SpeCode==795, 246, ifelse(SpeCode==825, 518, ifelse(SpeCode==960, 158, ifelse(SpeCode==956, 152, ifelse(SpeCode==251, 16, SpeCode)))))))%>%
  left_join(spp)%>%
  left_join(invert_trt)%>%
  mutate(genus_species=paste(genus, species))%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, Cover)%>%
  filter(genus_species!='bare ground')%>%
  filter(calendar_year<2019)%>% #cessation begins
  group_by(project_name, calendar_year, plot_id, treatment, genus_species)%>%
  summarise(abundance=max(Cover))%>%
  ungroup()

sp_nutnet <- read.csv('species_comp\\NUT011_nutnet.csv')%>%
  mutate(project_name="nutnet")%>%
  rename(calendar_year=RecYear, plot_id=Plot)%>%
  mutate(code=ifelse(Sppnum %in% c(900, 957), 517, ifelse(Sppnum==795, 246, ifelse(Sppnum %in% c(951,825), 518, ifelse(Sppnum==960, 158, ifelse(Sppnum==956, 152, ifelse(Sppnum==251, 16, ifelse(Sppnum==267, 52, ifelse(Sppnum==726, 246, Sppnum)))))))))%>%
  left_join(spp)%>%
  left_join(nutnet_trt)%>%
  filter(Taxa!='bare ground'&Taxa!='litter')%>%
  mutate(genus_species=paste(genus, species))%>%
  group_by(project_name, calendar_year, plot_id, treatment, genus_species)%>%
  summarise(abundance=max(Cover))%>%
  ungroup()

sp_bgp <- read.csv("species_comp\\BGPVC_belowgroundPlots.csv")%>%
  rename(calendar_year=RecYear, code=SpeciesCode, plot_id=Plot)%>%
  left_join(spp)%>%
  mutate(genus_species=paste(genus,species,sep='_'))%>%
  mutate(abundance=ifelse(CoverClass==1, 1, ifelse(CoverClass==2, 2.5, ifelse(CoverClass==3, 10, ifelse(CoverClass==4, 37.5, ifelse(CoverClass==5, 62.5, ifelse(CoverClass==6, 85, 97.5)))))))%>%
  group_by(genus_species, plot_id, calendar_year)%>%
  summarize(abundance=mean(abundance))%>% #average over subplots
  ungroup()%>%
  right_join(bgp_trt)%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)%>%
  filter(calendar_year<2017) #remove data after cessation of treatmetns

sp_ghostfire <- read.csv("species_comp\\GFE011_ghostFire.csv")%>%
  mutate(plot_id=paste(Block,Plot,sep=''),
         project_name=ifelse(BurnTrt=='Annual', 'GF Burned', 'GF Unburned'))%>%
  rename(code=Spnum, calendar_year=RecYear)%>%
  left_join(spp)%>%
  left_join(ghostfire_trt)%>%
  mutate(genus_species=paste(genus,species,sep='_'))%>%
  group_by(project_name, calendar_year, treatment, genus_species, plot_id)%>%
  summarise(abundance=max(Cover))%>%
  ungroup()

#without species numbers
sp_pplots <- read.csv("species_comp\\PPL011_pplots.csv")%>%
  mutate(project_name="pplots", 
         genus_species=paste(Genus, Species, sep='_'),
         treatment=as.factor(ifelse(Treatment=='n1p0', 'control', as.character(Treatment))))%>%
  rename(calendar_year=RecYear, plot_id=PlotID)%>%
  group_by(project_name, calendar_year, plot_id, treatment, genus_species)%>%
  summarise(abundance=max(Abundance))%>%
  ungroup()

sp_change_1 <- read.csv("species_comp\\KNZ_change_sppComp.csv")%>%
  rename(calendar_year=year, plot_id=plot, genus_species=species)%>%
  left_join(change_trt)%>%
  group_by(project_name, calendar_year, plot_id, treatment, genus_species)%>%
  summarise(abundance=max(cover))%>%
  ungroup()

sp_change_2 <- read.csv("species_comp\\ChANGE_sppcomp_2020-2023.csv")%>%
  mutate(project_name='ChANGE') %>% 
  rename(calendar_year=year, plot_id=plot, genus_species=species)%>%
  select(-treatment) %>% 
  left_join(change_trt)%>%
  group_by(project_name, calendar_year, plot_id, treatment, genus_species)%>%
  summarise(abundance=max(cover))%>%
  ungroup()

sp_change <- rbind(sp_change_1, sp_change_2)

sp_ukulinga <- read.csv('species_comp\\KNZ_ukulinga_sppComp.csv')%>%
  left_join(ukulinga_trt)%>%
  gather(key='year', value='cover', cov2005:cov2010)%>%
  separate(col=year, into=c('trash','calendar_year'), sep='v')%>%
  rename(code=spnum)%>%
  left_join(spp)%>%
  mutate(calendar_year=as.numeric(calendar_year),
         genus_species=paste(genus,species,sep='_'))%>%
  group_by(project_name, calendar_year, treatment, genus_species, plot, subplot, ssubplot)%>%
  summarise(max_cov=max(cover))%>%
  ungroup()%>%
  group_by(project_name, calendar_year, treatment, genus_species, plot, subplot)%>%
  summarise(abundance=mean(max_cov))%>% #averages across sub-subplots
  ungroup()%>%
  mutate(plot_id=paste(plot,subplot,sep='_'))%>%
  select(project_name, calendar_year, treatment, genus_species, abundance, plot_id)


##### Merging Species Composition Data #####

sp_all <- sp_pplots%>%
  rbind(sp_bgp)%>%
  rbind(sp_change)%>%
  rbind(sp_invert)%>%
  rbind(sp_nutnet)%>%
  rbind(sp_ghostfire)%>%
  rbind(sp_ukulinga)%>%
  filter(abundance>0)



##### Calculating Relative Cover #####

totCov <- sp_all%>%
  group_by(project_name, calendar_year, treatment, plot_id)%>%
  summarise(tot_cover=sum(abundance))%>%
  ungroup()

relCov <- sp_all%>%
  rename(cov=abundance)%>%
  left_join(totCov)%>%
  mutate(abundance=cov/tot_cover)

# write.csv(relCov, 'Konza_nutrient synthesis_spp comp_20240703.csv', row.names=F)



##### Read in Covaraite Data #####

#grasshoppers
grasshopper <- read.csv('drivers\\CGR022_grasshoppers.csv')%>%
  filter(!is.na(TOTAL), SPECIES %notin% c('Gryllidae spp.','Oecanthinae spp.','Tettigoniidae spp.','Orchelimum spp.','Neoconocephalus ensiger','Neoconocephalus robustus','Scudderia texensis'))%>% #remove non-Acrididae species, which were not counted in all years
  group_by(RECYEAR,RECMONTH,RECDAY,WATERSHED,SOILTYPE,REPSITE)%>%
  summarise(total_count=sum(TOTAL))%>% #sum of grasshoppers across all species for each rep
  ungroup()%>%
  filter(WATERSHED %notin% c('010d','0sub','n00b','n01a','n01b','n04a','n04d','n20a','n20b'))%>% #drop watersheds not in our experimental framework (e.g., 10 yr burn, grazed)
  filter(WATERSHED %notin% c('000b','004d','004g'))%>% #drop watersheds with too little data
  group_by(RECYEAR, RECMONTH, RECDAY, WATERSHED, SOILTYPE)%>%
  summarise(trans_count=sum(total_count))%>% #sum across reps within a transect
  ungroup()%>%
  group_by(RECYEAR, WATERSHED, SOILTYPE)%>%
  summarise(max_count=max(trans_count))%>% #max across two sampling time points (early vs late season)
  ungroup()%>%
  group_by(RECYEAR, WATERSHED)%>%
  summarise(ws_count=mean(max_count))%>% #average over transects within a watershed
  ungroup()%>%
  mutate(burn=ifelse(WATERSHED %in% c('001d','0spb'), 'annual', ifelse(WATERSHED %in% c('002c','002d'), 'two', ifelse(WATERSHED %in% c('004b','004f'), 'four', 'twenty'))))%>%
  group_by(RECYEAR, burn)%>%
  summarise(burn_count=mean(ws_count))%>% #average by burn treatment
  ungroup()%>%
  group_by(RECYEAR)%>%
  summarise(grasshopper=mean(burn_count))%>% #average across all burn types
  ungroup()%>%
  rename(calendar_year=RECYEAR)

# ggplot(data=grasshopper, aes(x=calendar_year, y=grasshopper)) +
#   geom_point() + geom_smooth(method='lm', se=F)

#small mammals
mammal_early <- read.csv('drivers\\CSM01_small mammals_1981-2013.csv')%>%
  separate(col=WATERSHED.LINE, into=c('WATERSHED','LINE'), sep='-')%>%
  filter(WATERSHED %in% c('001D','004B','004F','020B'))%>%  #drop watersheds not in our experimental framework (e.g., grazed)
  filter(SEASON=='AU')%>% #only keep fall sampling dates for comparability to newest records
  mutate(count=Pm+Rmeg+Sh+Bh+Rmon+St+Mo+Pl+Ch+Mm+Nf+Sc+Zh+Cp)%>% #calculate sum across species
  group_by(RECYEAR, WATERSHED, LINE)%>%
  summarise(trans_count=sum(count))%>% #sum across four sampling dates
  ungroup()%>%
  group_by(RECYEAR, WATERSHED)%>%
  summarise(ws_count=mean(trans_count))%>% #mean across transects within watershed
  ungroup()%>%
  mutate(burn=ifelse(WATERSHED %in% c('001D'), 'annual', ifelse(WATERSHED %in% c('004B','004F'), 'four', 'twenty')))%>%
  group_by(RECYEAR, burn)%>%
  summarise(burn_count=mean(ws_count))%>% #mean by burn treatment
  ungroup()%>%
  group_by(RECYEAR)%>%
  summarise(mammal=mean(burn_count))%>% #mean across all burn types
  ungroup()%>%
  rename(calendar_year=RECYEAR)
  
mammal_late <- read.csv('drivers\\CSM08_small mammals_2016-2021.csv')%>%
  filter(Watershed %in% c('1D','4B','4F','20B'))%>%  #drop watersheds not in our experimental framework (e.g., grazed)
  filter(RecaptureStatus=='N')%>% #filter out recaptures
  group_by(RecYear, Watershed, Transect)%>%
  summarise(trans_count=length(Species))%>% #sum across four sampling dates
  ungroup()%>%
  group_by(RecYear, Watershed)%>%
  summarise(ws_count=mean(trans_count))%>% #mean across transects within watershed
  ungroup()%>%
  mutate(burn=ifelse(Watershed %in% c('1D'), 'annual', ifelse(Watershed %in% c('4B','4F'), 'four', 'twenty')))%>%
  group_by(RecYear, burn)%>%
  summarise(burn_count=mean(ws_count))%>% #mean by burn treatment
  ungroup()%>%
  group_by(RecYear)%>%
  summarise(mammal=mean(burn_count))%>% #mean across all burn types
  ungroup()%>%
  rename(calendar_year=RecYear)

mammal <- mammal_early%>%
  rbind(mammal_late)

# ggplot(data=mammal, aes(x=calendar_year, y=mammal)) +
#   geom_point() + geom_smooth(method='lm', se=F)

#weather
precip <- read.csv('drivers\\AWE012_weather.csv')%>%
  filter(DPPT!='.') %>%
  mutate(month=paste('a',RECMONTH, sep=''),
         DPPT=as.integer(DPPT))%>%
  mutate(year=ifelse(RECMONTH %in% c(9, 10, 11, 12), RECYEAR+1, RECYEAR)) %>% 
  group_by(year, month)%>%
  summarise(precip=sum(DPPT))%>%
  ungroup()%>%
  spread(key=month, value=precip, fill=0)%>%
  mutate(grow_precip=(a4+a5+a6+a7+a8),
         year_precip=(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12))%>%
  filter(year>1982, year<2024) %>% 
  rename(calendar_year=year)%>%
  select(calendar_year, grow_precip, year_precip)

growTemp <- read.csv('drivers\\AWE012_weather.csv')%>%
  filter(TAVE!='.', TAVE>-30) %>% 
  mutate(temp=as.numeric(TAVE))%>%
  filter(RECMONTH %in% c(4,5,6,7,8))%>%
  group_by(RECYEAR)%>%
  summarise(grow_temp=mean(temp))%>%
  ungroup()%>%
  rename(calendar_year=RECYEAR)

yearTemp <- read.csv('drivers\\AWE012_weather.csv')%>%
  filter(TAVE!='.', TAVE>-30) %>% 
  mutate(temp=as.numeric(TAVE))%>%
  mutate(year=ifelse(RECMONTH %in% c(9, 10, 11, 12), RECYEAR+1, RECYEAR)) %>% 
  group_by(year)%>%
  summarise(year_temp=mean(temp))%>%
  ungroup()%>%
  filter(year>1982, year<2024)%>%
  rename(calendar_year=year) 

#combine drivers
drivers <- precip%>%
  full_join(growTemp)%>%
  full_join(yearTemp)%>%
  full_join(grasshopper)%>%
  full_join(mammal)

# write.csv(drivers, 'Konza_nutrient synthesis_drivers_20240304.csv', row.names=F)