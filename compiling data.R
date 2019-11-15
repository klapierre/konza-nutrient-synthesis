library(tidyverse)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's:
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')

#konza spp list
spp <- read.csv('konza_spplist.csv')
nutnetSp <- read.csv('nutnet_spplist.csv')%>%
  mutate(sppnum=SppNum)
viSp <- read.csv('invert_spplist.csv')%>%
  mutate(sppnum=SppNum)


#treatments
bgp_trt <- read.csv("belowground_plots_anpp_1989-2015.csv")%>%
  filter(MOW!="m")%>%
  mutate(treatment=paste(BURN, MOW, NUTRIENT, sep="_"),
         plot_id=PLOT)%>%
  select(plot_id,treatment)%>%
  unique()%>%
  filter(treatment!="__")

change_trt <- read.csv("ChANGE_treatments.csv")%>%
  mutate(plot_id=plot, treatment=as.factor(N))%>%
  select(plot_id, treatment)

nutnet_trt <- read.csv('KNZ_NutNet_trt.csv')%>%
  mutate(project_name='nutnet')%>%
  select(project_name, plot, treat_other_name)

ghostfire_trt <- read.csv('ghost fire_trt.csv')%>%
  mutate(plot_id=paste(Burn.Trt,Plot, sep='_'))%>%
  mutate(treatment=paste(Litter,Nutrient))%>%
  select(plot_id, treatment)


###species data
sp_pplots <- read.csv("pplots_sppcomp_2002-2018.csv")%>%
  mutate(project_name="pplots")%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)

sp_bgp <- read.csv("belowground_plots_sppcomp_1989-2017.csv")%>%
  mutate(project_name="BGP", 
         genus_species=paste(Ab_genus, Ab_species, sep=" "),
         plot_id=Plot,
         calendar_year=RecYear,
         abundance=ifelse(CoverClass==1, 1, ifelse(CoverClass==2, 2.5, ifelse(CoverClass==3, 10, ifelse(CoverClass==4, 37.5, ifelse(CoverClass==5, 62.5, ifelse(CoverClass==6, 85, 97.5)))))))%>%
  group_by(project_name, genus_species, plot_id, calendar_year)%>%
  summarize(abundance=mean(abundance))%>%
  ungroup()%>%
  #check with john about why some data has no spp id, genus, or species
  filter(genus_species!=' ')%>%
  right_join(bgp_trt)%>%
  mutate(project_name=ifelse(treatment=='u_u_n'|treatment=='u_u_c'|treatment=='u_u_p'|treatment=='u_u_b', 'BGP unburned', 'BGP burned'))%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)

sp_change <- read.csv("ChANGE_sppcomp_2013-2018.csv")%>%
  mutate(project_name='ChANGE',
         calendar_year=Year,
         plot_id=as.numeric(Plot), 
         genus_species=Species)%>%
  mutate(abundance=pmax(June,August))%>%
  select(project_name, calendar_year, plot_id, genus_species, abundance)%>%
  right_join(change_trt)

sp_invert <- read.csv('VI_sppcomp_2008-2019.csv')%>%
  mutate(project_name='invert', calendar_year=year, plot_id=plot, treatment=paste(NPK, exclose, insecticide, sep='_'), abundance=rel_cover, spnum=sppnum)%>%
  left_join(viSp)%>%
  mutate(genus_species=paste(genus, species))%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)%>%
  filter(genus_species!='bare ground')%>%
  filter(calendar_year<2019) #cessation begins

sp_nutnet <- read.csv('nutnet_sppcomp_2007-2019.csv')%>%
  select(-date, -note_cover)%>%
  spread(key=season, value=cover, fill=0)%>%
  group_by(year, plot, taxa)%>%
  mutate(cover=pmax(fall,spring))%>%
  ungroup()%>%
  left_join(nutnet_trt)%>%
  select(-taxa)%>%
  left_join(nutnetSp)%>%
  mutate(genus_species=ifelse(taxa=='Muhlenbergia cuspidata', 'Muhlenbergia racemosa', ifelse(taxa=='Euphorbia serpens', 'Euphorbia nutans', as.character(taxa))))%>%
  mutate(project_name='nutnet', calendar_year=year, plot_id=plot, treatment=treat_other_name, abundance=cover)%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)%>%
  filter(genus_species!='bare ground', treatment!='fence'&treatment!='NPK_fence')%>%
  filter(calendar_year<2019) #for consistency with other datasets

sp_ghostfire<-read.csv("ghost fire_spp comp_2014-2018.csv")%>%
  mutate(genus_species=paste(genus,species), project_name=ifelse(Burn.Trt=='Annual', 'GF Burned', 'GF Unburned'), calendar_year=Year, blockplot=paste(Block,Plot, sep=''), plot_id=paste(Burn.Trt, blockplot, sep='_'))%>%
  left_join(ghostfire_trt)%>%
  group_by(project_name, calendar_year, treatment, genus_species, plot_id)%>%
  summarise(abundance=max(cover))%>%
  ungroup()%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)

# sp_restoration <- read.csv('restoration plots_spp comp_1999-2012.csv')%>%
#   filter(DEPTH==1)%>%
#   select(-OBS, -BLOCK, -DEPTH, -NUTRIENT, -RESTORE_YR, -RESIN_NO3, -H, -R, -ANPP, -WPTRT)%>%
#   gather(key=genus_species, value=cover, ANGE:TRRE3)%>%
#   mutate(cover=as.numeric(cover))%>%
#   filter(cover>0)%>%
#   group_by(PLOT, TRTCOMB, YEAR, genus_species)%>%
#   summarise(abundance=mean(cover))%>%
#   ungroup()%>%
#   mutate(project_name='restoration', calendar_year=YEAR, plot_id=PLOT, treatment=TRTCOMB)%>%
#   select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)
# 
# sp_ukulinga_all <- read.csv('KNZ UK_spp comp_2005-2010.csv')%>%
#   unite(plot_label, plot, fert, sep='_', remove=F)
# 
# ukulinga_plots <- sp_ukulinga_all%>%
#   select(plot_label)%>%
#   unique()%>%
#   arrange(plot_label)%>%
#   mutate(plot_id=seq(1,6, by=1))
# 
# sp_ukulinga <- sp_ukulinga_all%>%
#   merge(ukulinga_plots, by=c('plot_label'))%>%
#   gather(key=covyear, value=cover, cov2005:cov2010)%>%
#   filter(cover>0)%>%
#   separate(covyear, c('cov', 'year'), sep='v')%>%
#   group_by(site, fert, plot_id, subplot, spnum, year)%>%
#   summarise(cover=mean(cover))%>%
#   ungroup()%>%
#   group_by(site, fert, plot_id, spnum, year)%>%
#   summarise(cover=mean(cover))%>%
#   ungroup()%>%
#   merge(knz_spplist, by=c('spnum'))%>%
#   mutate(treatment=ifelse(fert==1, 'control', 'N'), project_name=ifelse(site=='1D', 'ukulinga annual', ifelse(site=='4F', 'ukulinga four', 'ukulinga unburned')), abundance=cover, calendar_year=year)%>%
#   select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)


###merging
#species data
sp_all <- sp_pplots%>%
  rbind(sp_bgp)%>%
  rbind(sp_change)%>%
  rbind(sp_invert)%>%
  rbind(sp_nutnet)%>%
  rbind(sp_ghostfire)
  # rbind(sp_restoration)%>%
  # rbind(sp_ukulinga)

# write.csv(sp_all, 'Konza_nutrient synthesis_spp comp_11152019.csv', row.names=F)

# #anpp data
# anpp_all <- anpp_pplots%>%
#   rbind(anpp_bgp)%>%
# #  rbind(anpp_change)%>%
#   rbind(anpp_invert)%>%
#   rbind(anpp_nutnet)%>%
#   rbind(anpp_ghostfire)%>%
#   rbind(anpp_restoration)





