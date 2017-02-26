library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)
library(doBy)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's:
setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis')

#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

##first get a list of controls
controls<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(control=ifelse(treatment=="N1P0", 1, ifelse(treatment=="b_u_c", 1, ifelse(treatment==0, 1, ifelse(treatment=="control", 1, ifelse(treatment=="control_control", 1, ifelse(treatment=="u_u_c", 1, 0)))))))%>%
  select(project_name, treatment, control)%>%
  unique()

##select N10 treatments
nitrogen<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(nitrogen=ifelse(treatment=="N2P0", 1, ifelse(treatment=="N2P3",1, ifelse(treatment=="b_u_n", 1, ifelse(treatment=="b_u_b",1, ifelse(treatment==10, 1, ifelse(treatment==5, 1, ifelse(treatment=="NPK", 1, ifelse(treatment=="N",1, ifelse(treatment=="control_nitrogen", 1, ifelse(treatment=="NP", 1,0)))))))))))%>%
  select(project_name, treatment, nitrogen)%>%
  unique()

#read in climate data
precip<-read.csv("WETDRY.csv")%>%
  mutate(calendar_year=YEAR, gprecip=AprSeptPrecip)%>%
  select(calendar_year, gprecip)

##do species data
spdata<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(id=paste(project_name, calendar_year, sep="::"))

##richness
sp_rich <- group_by(spdata, project_name, treatment, calendar_year, plot_id) %>% 
  summarize(S=S(abundance))%>%
  tbl_df()%>%
  group_by(project_name, treatment, calendar_year)%>%
  summarize(S=mean(S))

sp_rich_control<-merge(sp_rich, controls, by=c("treatment","project_name"))%>%
  filter(control==1)%>%
  mutate(c_S=S)%>%
  select(-S, -control,-treatment)

sp_rich_rel1<-merge(sp_rich, sp_rich_control, by=c("project_name","calendar_year"))%>%
  mutate(PC_S=(S-c_S)/c_S)

sp_rich_rel<-merge(sp_rich_rel1, controls, by=c("treatment","project_name"))%>%
  filter(control!=1)

nit_sp<-merge(sp_rich_rel, nitrogen, by=c("project_name","treatment"))%>%
  filter(nitrogen==1)

#richness graph
ggplot(data=nit_sp, aes(x=calendar_year, y=PC_S, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_line()+
  geom_hline(yintercept = 0)+
  #theme(legend.position = "none")+
  facet_wrap(~project_name, scales="free")


###ANPP data
anppdata<-read.csv("Konza_nutrient synthesis_anpp.csv")%>%
  tbl_df()%>%
  group_by(project_name, treatment, calendar_year)%>%
  summarize(anpp=mean(anpp))

anpp_control<-merge(anppdata, controls, by=c("treatment","project_name"))%>%
  filter(control==1)%>%
  mutate(c_anpp=anpp)%>%
  select(-anpp, -control,-treatment)

anpp_rel1<-merge(anppdata, anpp_control, by=c("project_name","calendar_year"))%>%
  mutate(PC_anpp=(anpp-c_anpp)/c_anpp)

anpp_rel<-merge(anpp_rel1, controls, by=c("treatment","project_name"))%>%
  filter(control!=1)

write.csv(anpp_rel,"anpp_percent_change.csv", row.names = F)

nit_anpp<-merge(anpp_rel, nitrogen, by=c("project_name","treatment"))%>%
  filter(nitrogen==1)
#anpp graph
# ggplot(data=anpp_rel, aes(x=calendar_year, y=PC_anpp, group=treatment))+
#   geom_point(aes(color=treatment))+
#   geom_line()+
#   geom_hline(yintercept = 0)+
#   #theme(legend.position = "none")+
#   facet_wrap(~project_name, scales="free")
# 
# nit_anpp_precip<-merge(nit_anpp, precip, by="calendar_year")

##graphing this
ggplot(data=anppdata, aes(x=calendar_year, y=anpp, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = "lm", se=F, color="black")+
  facet_wrap(~project_name, scales="free")

###looking at bgp only
bgp_precip<-nit_anpp_precip%>%
  filter(project_name=="BGP")%>%
  mutate(shift=as.factor(ifelse(calendar_year>1991, 1, 0)))

ggplot(data=bgp_precip, aes(x=gprecip, y=PC_anpp, group=treatment))+
    geom_point(aes(color=shift, shape=treatment), size=5)

###looking with precip and not percent change
nitrogen_only<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(nitrogen=ifelse(treatment=="N2P0", 1, ifelse(treatment=="b_u_n", 1,ifelse(treatment==10, 1, ifelse(treatment=="N",1, ifelse(treatment=="control_nitrogen", 1, 0))))))%>%
  select(project_name, treatment, nitrogen)%>%
  unique()

anpp_precip<-merge(anppdata, precip, by="calendar_year")

anpp_precip1<-merge(anpp_precip, nitrogen_only, by=c('project_name','treatment'))
anpp_precip2<-merge(anpp_precip1, controls, by=c("project_name","treatment"))%>%
  filter(nitrogen==1|control==1)

ggplot(anpp_precip2, aes(x=gprecip, y=anpp, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=F, color="black")+
  facet_wrap(~project_name, scales="free")



##Do BC community shifts


spdata2<-merge(spdata, controls, by=c("project_name","treatment"))

bc_mean_change=data.frame(project_name=c(), calendar_year=c(), treatment=c(),mean_change=c()) 

pjt_nm<-unique(spdata2$id)

###first, get bray curtis dissimilarity values for each year within each experiment between all combinations of plots
###second, get distance of each plot within a trt to the trt centroid 
###third: mean_change is the distance between trt and control centriods
for(i in 1:length(pjt_nm)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=spdata2%>%
    filter(id==pjt_nm[i], abundance!=0)%>%
    select(project_name, calendar_year, treatment, genus_species, abundance, plot_id, control)
  
  test<-subset%>%
    group_by(plot_id, calendar_year, genus_species)%>%
    summarize(dup=length(abundance))
  
  #get the name of the control treatment
  cont<-subset%>%
    select(treatment, control)%>%
    unique%>%
    filter(control==1)%>%
    select(treatment)
  
  #transpose data
  species=subset%>%
    spread(genus_species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,6:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(id=pjt_nm[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==cont$treatment,]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  mc<-cent_C_T%>%
    separate(id, c("project_name","calendar_year"), sep="::")
  
  #pasting dispersions into the dataframe made for this analysis
  bc_mean_change=rbind(mc, bc_mean_change)  
}

mean_change<-merge(bc_mean_change, controls, by=c("treatment","project_name"))%>%
  filter(control!=1)

write.csv(mean_change, "community_change_nut_experiments.csv", row.names = F)

###graphing BC
##overall
ggplot(data=mean_change, aes(x=calendar_year, y=mean_change, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_line()+
  theme(legend.position="none")+
  facet_wrap(~project_name, scales = "free_x")

##pretty graphs
# red= n only
# blue=p only
# purple=n+p
#other trt=black
theme_set(theme_bw(12))

pplots<-
  ggplot(data=subset(mean_change, project_name=="pplots"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N0 P2.5", "N0 P5", "N0 P10","N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("blue","blue","blue","red","purple","purple","purple"))+
  geom_line()+
  geom_vline(xintercept = 2004, linetype="longdash")+
  geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)

BGP_ub<-
ggplot(data=subset(mean_change, project_name=="BGP unburned"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P1", "N10 P0", "N0 P1"), values=c("purple","red","blue"))+
  geom_line()+
  geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Belowground Plots Unburned")+
  ylim(min=0, max=0.8)

BGP_b<-
ggplot(data=subset(mean_change, project_name=="BGP burned"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P1", "N10 P0", "N0 P1"), values=c("purple","red","blue"))+
  geom_line()+
  geom_vline(xintercept = 1989,linetype="longdash")+
  geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Belowground Plots Burned")+
  ylim(min=0, max=0.8)

nutnet<-
ggplot(data=subset(mean_change, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P0","N0 P0 K10" ,"N10 P0 K10","N10 P10 K0" ,"N10 P10 K10","N0 P10 K0", "N0 P10 K10"), values=c("red","black", "red", "purple","purple","blue","blue"))+
  geom_line()+
  geom_vline(xintercept = 2010,linetype="longdash")+
  geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Nutrient Network")+
  ylim(min=0, max=0.8)

invert<-
ggplot(data=subset(mean_change, project_name=="invert"&treatment!="caged"&treatment!="caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10 P10 K10", "Insecticide", "N10 P10 K10\n& Insecticide"), values=c("purple","black","purple"))+
  geom_line()+
  geom_vline(xintercept = 2011,linetype="longdash")+
  geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)

restoration<-
ggplot(data=subset(mean_change, project_name=="restoration"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N5", "Carbon"), values=c("red","black"))+
  geom_line()+
  geom_vline(xintercept = 2005, linetype="longdash")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Restoration Plots")+
  ylim(min=0, max=0.8)

uk_annual<-
  ggplot(data=subset(mean_change, project_name=="ukulinga annual"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Ukulinga Annual Burn")+
  ylim(min=0, max=0.8)

uk_four<-
  ggplot(data=subset(mean_change, project_name=="ukulinga four"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Ukulinga 4-Yr Burn")+
  ylim(min=0, max=0.8)

uk_ub<-
  ggplot(data=subset(mean_change, project_name=="ukulinga unburned"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Ukulinga Unburned")+
  ylim(min=0, max=0.8)

grid.arrange(pplots, BGP_ub, BGP_b, nutnet, invert, restoration, uk_annual, uk_four, uk_ub)


#####
##### doing NMDS with all experiments
#####

sp_names<-spdata2%>%
  select(genus_species)%>%
  unique()
##export this and reimport corrected names.
#write.csv(sp_names, "species_list_allexp.csv", row.names=F)
sp_names_clean<-read.csv("species_list_allexp_cleaned.csv")

spdata3<-merge(sp_names_clean, spdata, by="genus_species")

####this is not working for BGP! why???
###export and reimport with BGP fixed
# bgp<-spdata%>%
#   filter(project_name=="BGP unburned"|project_name=="BGP burned")
# 
# write.csv(bgp,"species_data_BGP_problems.csv", row.names = F)

bgpfix<-read.csv("species_data_BGP_problems_fix.csv")

bgpmerge<-merge(bgpfix, sp_names_clean, by="genus_species")

spdata4<-rbind(bgpmerge, spdata3)

###finally do NMDS
lst_yr<-spdata4%>%
  group_by(project_name)%>%
  mutate(lastyear=max(calendar_year))%>%
  filter(calendar_year==lastyear)%>%
  select(project_name, treatment, plot_id, abundance, species)%>%
  filter(project_name!="ChANGE"&project_name!="ghost fire")%>%
  filter(treatment!="caged"&treatment!="caged_insecticide"&treatment!="fence"&treatment!="NPK_caged"&treatment!="NPK_caged_insecticide"&treatment!="NPKfence")

sp_wide<-lst_yr%>%
  spread(species, abundance, fill=0)

mds<-metaMDS(sp_wide[,4:131],autotransform=FALSE, shrink=FALSE, trymax = 1000)
mds

###plotting this
theme_set(theme_bw(12))

plots<-as.data.frame(sp_wide[,1:3])
scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores)

funcs<-function(x)c(mn=mean(x),se=sd(x)/sqrt(length(x)))
means<-summaryBy(NMDS1+NMDS2~treatment, data=scores2, FUN=funcs)

ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=treatment))+
  geom_point(size=5,aes(shape=project_name))+
  scale_shape_manual(name="Experiment", values=c(0,1,2,5,6,7,9,10,12))+
  scale_color_manual(name="Treatment", values=c("purple","green","red","blue","purple","green","red","blue", "black","green","black","black","red","green","blue","blue","blue","red","purple","purple","purple","red","purple","purple","purple","purple","blue","blue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

