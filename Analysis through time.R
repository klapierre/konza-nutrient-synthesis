setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(ggplot2)

####subset konza for sally
reldata<-read.csv("SpeciesRawAbundance_Dec2016.csv")
konza<-reldata%>%
  filter(site_code=="KNZ")
write.csv(konza, "~/Dropbox/Konza Nutrient Synthesis/Konza_subset.csv")

######take corre data metrics and look at bgp and pplots through time

sp_dat<-read.csv("ForBayesianAnalysis_Dec2016.csv")%>%
  filter(project_name=="pplots"|project_name=="BGP")%>%
  filter(treatment=="b_u_b"|treatment=="b_u_c"|treatment=="b_u_n"|treatment=="b_u_p"|treatment=="N1P1"|treatment=="N1P2"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P1"|treatment=="N2P2"|treatment=="N2P3")

anpp_dat<-read.csv("ForBayesianAnalysisANPP_Dec2016.csv")%>%
  filter(project_name=="pplots"|project_name=="BGP")%>%
  filter(treatment=="b_u_b"|treatment=="b_u_c"|treatment=="b_u_n"|treatment=="b_u_p"|treatment=="N1P1"|treatment=="N1P2"|treatment=="N1P3"|treatment=="N2P0"|treatment=="N2P1"|treatment=="N2P2"|treatment=="N2P3")

#richness
ggplot(data=dat, aes(x=calendar_year, y=S_PC, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_line()+
  facet_wrap(~project_name, scales="free")
#mean change
ggplot(data=sp_dat, aes(x=calendar_year, y=mean_change, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_line()+
  facet_wrap(~project_name, scales="free")
#ANPP
ggplot(data=anpp_dat, aes(x=calendar_year, y=anpp_PC, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_line()+
  facet_wrap(~project_name, scales="free")
