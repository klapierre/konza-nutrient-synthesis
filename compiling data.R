setwd("~/Dropbox/Konza Nutrient Synthesis")

#treatments
bgp_trt<-read.csv("BGPE_ANPP_1986-2015.csv")%>%
  filter(Mow!="m")%>%
  tbl_df()%>%
  mutate(treatment=paste(Burn, Mow, Nutrient, sep="_"),
         plot_id=Plot)%>%
  select(plot_id,treatment)%>%
  unique()%>%
  filter(treatment!="__")

###species data
sp_pplots<-read.csv("pplots_spp comp_2002-2015.csv")%>%
  mutate(project_name="pplots")

sp_bgp_clean<-read.csv("BGPE_spp comp_1986-2012.csv")%>%
  mutate(project_name="BGP", 
         genus_species=paste(Ab_genus, Ab_species, sep=" "),
         plot_id=Plot,
         calendar_year=RecYear,
         abundance=ifelse(CoverClass==1, 1, ifelse(CoverClass==2, 2.5, ifelse(CoverClass==3, 10, ifelse(CoverClass==4, 37.5, ifelse(CoverClass==5, 62.5, ifelse(CoverClass==6, 85, 97.5)))))))%>%
  group_by(project_name, genus_species, plot_id, calendar_year)%>%
  summarize(abundance=mean(abundance))
  
sp_bgp<-merge(bgp_trt, sp_bgp_clean, by="plot_id")           
  

sp_change<-read.csv("ChANGE_spp comp_2013-2016.csv")%>%
  tbl_df%>%
  mutate(project_name=Experiment,
         calendar_year=Year,
         plot_id=Plot, 
         genus_species=Species)%>%
  group_by(project_name, calendar_year, plot_id, genus_species)%>%
  mutate(abundance=max(June,August))%>%
  select(project_name, calendar_year, plot_id, genus_species, abundance)
  

###anpp data
anpp_bgp_raw<-read.csv("BGPE_ANPP_1986-2015.csv")%>%
  group_by(RecYear)%>%
  mutate(maxmonth=max(RecMonth))%>%
  filter(maxmonth==RecMonth)
  
anpp_bgp_cleaned<-anpp_bgp_raw%>%
  mutate(calendar_year=RecYear,plot_id=Plot, 
         anpp=Lvgrass+Forbs+Cuyrdd+Woody,
         project_name="BGP")%>%
  group_by(calendar_year, plot_id, project_name)%>%
  summarize(anpp=mean(anpp))
         
anpp_bgp<-merge(anpp_bgp_cleaned, bgp_trt, by="plot_id")

anpp_pplots<-read.csv("pplots_anpp_2002-2015.csv")%>%
  mutate(project_name="pplots")
