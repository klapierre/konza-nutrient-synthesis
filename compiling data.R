library(tidyverse)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's:
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data\\2018 analysis')

#konza spp list
spp <- read.csv('konza_spplist.csv')


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
  mutate(site='nutnet')


###species data
sp_pplots <- read.csv("pplots_sppcomp_2002-2016.csv")%>%
  mutate(project_name="pplots")

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
  mutate(project_name=ifelse(treatment=='u_u_n'|treatment=='u_u_c'|treatment=='u_u_p'|treatment=='u_u_b', 'BGP unburned', 'BGP burned'))

sp_change <- read.csv("ChANGE_sppcomp_2013-2017.csv")%>%
  mutate(project_name='ChANGE',
         calendar_year=Year,
         plot_id=as.numeric(Plot), 
         genus_species=Species)%>%
  group_by(project_name, calendar_year, plot_id, genus_species)%>%
  mutate(abundance=as.numeric(max(June,August)))%>%
  ungroup()%>%
  select(project_name, calendar_year, plot_id, genus_species, abundance)%>%
  right_join(change_trt)

sp_invert <- read.csv('VI_sppcomp_2008-2017.csv')%>%
  mutate(project_name='invert', calendar_year=year, plot_id=plot, treatment=paste(NPK, exclose, insecticide, sep='_'), abundance=rel_cover, spnum=sppnum)%>%
  left_join(spp)%>%
  mutate(genus_species=paste(genus, spp))%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)%>%
  filter(genus_species!='bare ground')

sp_nutnet <- read.csv('nutnet_sppcomp_2007-2017.csv')%>%
  select(-date)%>%
  spread(key=season, value=cover, fill=0)%>%
  group_by(year, plot, taxa)%>%
  mutate(cover=max(fall,spring))%>%
  ungroup()%>%
  left_join(nutnet_trt)%>%
  mutate(project_name='nutnet', calendar_year=year, plot_id=plot, treatment=treat_other_name, genus_species=taxa, abundance=cover)%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)%>%
  filter(genus_species!='bare ground')

# sp_ghostfire_clean<-read.csv("ghost fire_spp comp_2014-2015.csv")%>%
#   mutate(genus_species=Species,
#   project_name='ghost fire',
#   calendar_year=Year)%>%
#   group_by(project_name, calendar_year, genus_species)%>%
#   mutate(abundance=max(June,August))%>%
#   select(project_name, calendar_year, genus_species, abundance, Burn.Trt, Block, Plot)
# 
# sp_ghostfire<-merge(sp_ghostfire_clean, ghostfire_trt, by=c('Burn.Trt', 'Block', 'Plot'))%>%
#   select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)

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

knz_spplist <- read.csv('konza_spplist.csv')

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

  


# ###anpp data
# anpp_bgp_raw<-read.csv("BGPE_ANPP_1986-2015.csv")%>%
#   group_by(RecYear)%>%
#   mutate(maxmonth=max(RecMonth))%>%
#   filter(maxmonth==RecMonth)
# 
# anpp_bgp_cleaned<-anpp_bgp_raw%>%
#   mutate(calendar_year=RecYear,plot_id=Plot, 
#          anpp=Lvgrass+Forbs+Cuyrdd+Woody,
#          project_name="BGP")%>%
#   group_by(calendar_year, plot_id, project_name)%>%
#   summarize(anpp=mean(anpp))
# 
# anpp_bgp<-merge(anpp_bgp_cleaned, bgp_trt, by="plot_id")%>%
#   mutate(project_name=ifelse(treatment=='u_u_n'|treatment=='u_u_c'|treatment=='u_u_p'|treatment=='u_u_b', 'BGP unburned', 'BGP burned'), anpp=10*anpp)
# 
# anpp_pplots<-read.csv("pplots_anpp_2002-2015.csv")%>%
#   mutate(project_name="pplots")
# 
# anpp_invert <- read.csv('Vert Invert_anpp_2009-2015.csv')%>%
#   mutate(project_name='invert', calendar_year=date, plot_id=plot, treatment=trt_other_name, anpp=10*anpp)%>%
#   select(project_name, calendar_year, plot_id, treatment, anpp)
# 
# anpp_nutnet <- read.csv('NutNet_anpp_2007-2015.csv')%>%
#   mutate(project_name='nutnet', calendar_year=year, plot_id=plot, treatment=treat_other_name, anpp=total, anpp=10*anpp)%>%
#   select(project_name, calendar_year, plot_id, treatment, anpp)
# 
# anpp_ghostfire <- read.csv('ghost fire_anpp_2014-2015.csv')%>%
#   mutate(project_name='ghost fire', calendar_year=Year, Burn.Trt=ifelse(BurnFreq==1, 'Annual', 'Unburned'), anpp=(Grass+Forb+Woody), anpp=10*anpp)%>%
#   merge(ghostfire_trt, by=c('Burn.Trt', 'Block', 'Plot'))%>%
#   select(project_name, calendar_year, plot_id, treatment, anpp)
# 
# anpp_restoration <- read.csv('restoration plots_spp comp_1999-2012.csv')%>%
#   filter(DEPTH==1)%>%
#   select(PLOT, SUBPLOT, TRTCOMB, YEAR, ANPP)%>%
#   filter(ANPP!='.')%>%
#   mutate(ANPP=as.numeric(ANPP))%>%
#   group_by(PLOT, TRTCOMB, YEAR)%>%
#   summarise(anpp=mean(ANPP))%>%
#   ungroup()%>%
#   mutate(project_name='restoration', calendar_year=YEAR, plot_id=PLOT, treatment=TRTCOMB)%>%
#   select(project_name, calendar_year, plot_id, treatment, anpp)


###merging
#species data
sp_all <- sp_pplots%>%
  rbind(sp_bgp)%>%
  rbind(sp_change)%>%
  rbind(sp_invert)%>%
  rbind(sp_nutnet)
  # rbind(sp_ghostfire)%>%
  # rbind(sp_restoration)%>%
  # rbind(sp_ukulinga)

# #anpp data
# anpp_all <- anpp_pplots%>%
#   rbind(anpp_bgp)%>%
# #  rbind(anpp_change)%>%
#   rbind(anpp_invert)%>%
#   rbind(anpp_nutnet)%>%
#   rbind(anpp_ghostfire)%>%
#   rbind(anpp_restoration)





