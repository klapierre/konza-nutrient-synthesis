library(tidyverse)
library(codyn)

setwd("C://Users/mavolio2/Dropbox/Konza Research/Nutrient synthesis")

data<-read.csv("Konza_nutrient synthesis_spp comp_04142020.csv")
exyrs<-read.csv("experiment years.csv")

theme_set(theme_bw(12))

totabund<-data%>%
  group_by(project_name, plot_id, calendar_year)%>%
  summarise(tot=sum(abundance))

relabund<-data%>%
  left_join(totabund)%>%
  mutate(relabund=abundance/tot)%>%
  filter(relabund!="0")

yrs<-data%>%
  select(project_name, calendar_year)%>%
  unique%>%
  group_by(project_name)%>%
  summarise(start=min(calendar_year))%>%
  ungroup()%>%
  select(project_name, start)%>%
  rename(calendar_year=start)

trt_CN<-data%>%
  select(project_name, treatment)%>%
  unique()%>%
  mutate(trt=ifelse(treatment=="N1P0", "C", ifelse(treatment=="N2P0", "N", ifelse(treatment=="u_u_c", "C", ifelse(treatment=='u_u_n', "N", ifelse(treatment=="b_u_c", "C", ifelse(treatment=="b_u_n", "N", ifelse(treatment==0, "C", ifelse(treatment==10, "N", ifelse(treatment=="x_x_x", "C", ifelse(treatment=="NPK_x_x", "N", ifelse(treatment=="control", "C", ifelse(treatment=="N", "N", ifelse(treatment=="A C"&project_name=="GF Burned", "C", ifelse(treatment=="A U"&project_name=="GF Burned", "N", ifelse(treatment=="P C"&project_name=="GF Unburned", "C", ifelse(treatment=="P U"&project_name=="GF Unburned", "N", "ign")))))))))))))))))%>%
  filter(trt!="ign")

nplots<-data%>%
  right_join(trt_CN)%>%
  right_join(yrs)%>%
  filter(trt=="C")%>%
  select(project_name, plot_id)%>%
  unique()%>%
  group_by(project_name)%>%
  summarize(nplot=length(plot_id))

c_domsp<-relabund%>%
  right_join(trt_CN)%>%
  right_join(yrs)%>%
  filter(trt=="C")%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)%>%
  select(project_name, genus_species)
#using the control plots doesn't work well b/c sometimes they are not in the N treated plots and trackign thier change through time is odd. Also most of these experiments the first year is pre-treatment

n_domsp<-relabund%>%
  right_join(trt_CN)%>%
  right_join(yrs)%>%
  filter(trt=="N")%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)%>%
  select(project_name, genus_species)

domsp_change<-relabund%>%
  right_join(trt_CN)%>%
  right_join(n_domsp)%>%
  filter(trt!="C")%>%
  group_by(project_name, genus_species, calendar_year)%>%
  summarize(mabund=mean(relabund))%>%
  mutate(sp=ifelse(genus_species=="androp gerar"|genus_species=="andropogon gerardii"|genus_species=="Andropogon gerardii"|genus_species=="andropogon_gerardii", "A. gerardii", ifelse(genus_species=="andropogon scoparius"|genus_species=="andropogon_scoparius"|genus_species=="schiza scopa"|genus_species=="Schizachyrium scoparius","S. scoparius", ifelse(genus_species=="artemisia ludoviciana", "A. ludoviciana", ifelse(genus_species=="bothri bladh", "B. bladhii", ifelse(genus_species=="bouteloua curtipendula", "B. curipendula", ifelse(genus_species=="Panicum virgatum", "P. virgatum", ifelse(genus_species=="poa pratensis", "P. pratensis", ifelse(genus_species=="Solidago missouriensis", "S. missouriensis", ifelse(genus_species=="sorgha nutan"|genus_species=="Sorghastrum nutans"|genus_species=="sorghastrum_nutans", "S. nutans", "na"))))))))))

ggplot(data=domsp_change, aes(x=calendar_year, y=mabund, group=sp, color=sp))+
  geom_point()+
  geom_line()+
  facet_wrap(~project_name, scales="free")


N_yr7<-relabund%>%
  right_join(trt_CN)%>%
  right_join(exyrs)%>%
  filter(trt=="N")%>%
  filter(experiment_year==7)%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)
  
N_yr14<-relabund%>%
  right_join(trt_CN)%>%
  right_join(exyrs)%>%
  filter(trt=="N")%>%
  filter(experiment_year==14)%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)

N_yr9<-relabund%>%
  right_join(trt_CN)%>%
  right_join(exyrs)%>%
  filter(trt=="N")%>%
  filter(experiment_year==9)%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)

lyrs<-data%>%
  select(project_name, calendar_year)%>%
  unique%>%
  group_by(project_name)%>%
  summarise(end=max(calendar_year))%>%
  ungroup()%>%
  select(project_name, end)%>%
  rename(calendar_year=end)

N_yrlast<-relabund%>%
  right_join(trt_CN)%>%
  right_join(lyrs)%>%
  filter(trt=="N")%>%
  group_by(project_name, genus_species)%>%
  summarize(mabund=mean(relabund), n=length(relabund))%>%
  left_join(nplots)%>%
  mutate(freq=n/nplot,
         DCi=(mabund*freq)/2,
         rank=rank(-DCi))%>%
  filter(rank<4)


         