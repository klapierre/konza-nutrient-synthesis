library(tidyr)
library(dplyr)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's:
setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis')

#treatments
bgp_trt<-read.csv("BGPE_ANPP_1986-2015.csv")%>%
  tbl_df()%>%
  mutate(treatment=paste(Burn, Mow, Nutrient, sep="_"),
         plot_id=Plot)%>%
  select(plot_id,treatment)%>%
  unique()

###species data
sp_pplots<-read.csv("pplots_spp comp_2002-2015.csv")

sp_bgp<-read.csv("BGPE_spp comp_1986-2012.csv")%>%
  group_by()

sp_invert <- read.csv('Vert Invert_spp comp_2009-2016.csv')%>%
  mutate(project_name='invert', calendar_year=year, plot_id=plot, treatment=trt_other_name, genus_species=taxa, abundance=cover)%>%
  select(project_name, calendar_year, plot_id, treatment, genus_species, abundance)


###anpp data
anpp_bgp<-read.csv("BGPE_ANPP_1986-2015.csv")%>%
  group_by(RecYear)%>%
  mutate(maxmonth=max(RecMonth))%>%
  filter(maxmonth==RecMonth)
  
anpp_bgp_cleaned<-anpp_bgp%>%
  mutate(experiment_year=RecYear,plot_id=Plot, 
         anpp=Lvgrass+Forbs+Cuyrdd+Woody,
         project_name="BGP",
         community_type=0)%>%
  group_by(experiment_year, plot_id, project_name, community_type)%>%
  summarize(anpp=mean(anpp))
         
anpp_pplots<-read.csv("pplots_anpp_2002-2015.csv")

anpp_invert <- read.csv('Vert Invert_anpp_2009-2015.csv')%>%
  mutate(project_name='invert', calendar_year=date, plot_id=plot, treatment=trt_other_name)%>%
  select(project_name, calendar_year, plot_id, treatment, anpp)







