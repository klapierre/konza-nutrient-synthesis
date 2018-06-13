library(tidyverse)
library(vegan)
library(gridExtra)
library(doBy)
library(grid)
library(codyn)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's laptop:
setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis')

#kim's desktop:
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\konza projects\\Konza Nutrient Synthesis\\Threshold project\\data')



###data
#community data
community<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  #create a replicate variable unique to each experiment
  mutate(replicate=paste(project_name, treatment, plot_id, sep='::'))%>%
  #filter out 0s
  filter(abundance>0)
  # #drop ChANGE, ghost fire, and Ukulinga data because time series too short
  # #drop restoration plots because planted communities
  # filter(!(project_name %in% c('ukulinga annual', 'ukulinga four', 'ukulinga unburned', 'restoration', 'ChANGE', 'ghost fire')))

#list of experiments, treatments, and plots
plots <- community%>%
  select(project_name, treatment, plot_id, replicate, calendar_year)%>%
  unique()
  
#climate data
precip<-read.csv("WETDRY.csv")%>%
  mutate(calendar_year=YEAR, gprecip=AprSeptPrecip)%>%
  select(calendar_year, gprecip)


###calculate difference metrics


###calculate change metrics
#RAC change
RACchange <- RAC_change(community, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', replicate.var='replicate')%>%
  left_join(plots)



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

# write.csv(anpp_rel,"anpp_percent_change.csv", row.names = F)

nit_anpp<-merge(anpp_rel, nitrogen, by=c("project_name","treatment"))%>%
  filter(nitrogen==1)
#anpp graph
# ggplot(data=anpp_rel, aes(x=calendar_year, y=PC_anpp, group=treatment))+
#   geom_point(aes(color=treatment))+
#   geom_line()+
#   geom_hline(yintercept = 0)+
#   #theme(legend.position = "none")+
#   facet_wrap(~project_name, scales="free")

nit_anpp_precip<-merge(nit_anpp, precip, by="calendar_year")

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

# write.csv(mean_change, "community_change_nut_experiments.csv", row.names = F)

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
##### doing NMDS with all experiments (only keep N only treatments)
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
  filter(project_name!="ChANGE"&project_name!="ghost fire"&project_name!='restoration')%>%
  filter(treatment!="caged"&treatment!="caged_insecticide"&treatment!="fence"&treatment!="NPK_caged"&treatment!="NPK_caged_insecticide"&treatment!="NPKfence"&treatment!='u_u_b'&treatment!='u_u_p'&treatment!='b_u_b'&treatment!='b_u_p'&treatment!='insecticide'&treatment!='K'&treatment!='N1P1'&treatment!='N1P2'&treatment!='N1P3'&treatment!='N2P1'&treatment!='N2P2'&treatment!='N2P3'&treatment!='NK'&treatment!='NP'&treatment!='NPK_insecticide'&treatment!='P'&treatment!='PK')%>%
  mutate(treatment2=ifelse(treatment=='b_u_n', 'N', ifelse(treatment=='u_u_n', 'N', ifelse(treatment=='NPK', 'N', ifelse(treatment=='N', 'N', ifelse(treatment=='N2P0', 'N', 'control'))))))%>%
  select(-treatment)

sp_wide<-lst_yr%>%
  spread(species, abundance, fill=0)

mds<-metaMDS(sp_wide[,4:110],autotransform=FALSE, shrink=FALSE, trymax = 1000)
mds

###plotting this
theme_set(theme_bw(12))

plots<-as.data.frame(sp_wide[,1:3])
scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores)

funcs<-function(x)c(mn=mean(x),se=sd(x)/sqrt(length(x)))
means<-summaryBy(NMDS1+NMDS2~treatment2+project_name, data=scores2, FUN=funcs)

ggplot(means, aes(x=NMDS1.mn, y=NMDS2.mn, color=project_name))+
  geom_point(size=8,aes(shape=treatment2))+
  geom_errorbar(aes(ymin=NMDS2.mn-NMDS2.se, ymax=NMDS2.mn+NMDS2.se),color="black")+
  geom_errorbarh(aes(xmin=NMDS1.mn-NMDS1.se, xmax=NMDS1.mn+NMDS1.se),color="black")+
  scale_shape_manual(name="Treatment", values=c(16,17))+
  scale_color_manual(name="Experiment", values=c("purple","green","red","blue","orange","black","dark green","pink"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("NMDS1")+
  ylab("NMDS2")





###getting species appearances and disappearances to compare to the mean change through time
# codyn function modification ---------------------------------------------
#modifying codyn functions to output integer numbers of species appearing and disappearing, plus total spp number over two year periods, rather than ratios
#modifying codyn functions to output turnover from first year to all other years, rather than one year to the next
turnover_allyears <- function(df, 
                              time.var, 
                              species.var, 
                              abundance.var, 
                              metric=c("total", "disappearance","appearance")) {
  
  # allows partial argument matching
  metric = match.arg(metric) 
  
  # sort and remove 0s
  df <- df[order(df[[time.var]]),]
  df <- df[which(df[[abundance.var]]>0),]
  
  ## split data by year
  templist <- split(df, df[[time.var]])
  
  ## create two time points (first year and each other year)
  t1 <- templist[1]
  t2 <- templist[-1]
  
  ## calculate turnover for across all time points
  out <- Map(turnover_twoyears, t1, t2, species.var, metric)
  output <- as.data.frame(unlist(out))
  names(output)[1] = metric
  
  ## add time variable column
  alltemp <- unique(df[[time.var]])
  output[time.var] =  alltemp[2:length(alltemp)]
  
  # results
  return(output)
}

turnover_twoyears <- function(d1, d2, 
                              species.var, 
                              metric=c("total", "disappearance","appearance")){
  
  # allows partial argument matching
  metric = match.arg(metric)
  
  # create character vectors of unique species from each df
  d1spp <- as.character(unique(d1[[species.var]]))
  d2spp <- as.character(unique(d2[[species.var]]))
  
  # ID shared species
  commspp <- intersect(d1spp, d2spp)
  
  # count number not present in d2
  disappear <- length(d1spp)-length(commspp)
  
  # count number that appear in d2
  appear <- length(d2spp)-length(commspp)
  
  # calculate total richness
  totrich <- sum(disappear, appear, length(commspp))
  
  # output based on metric 
  if(metric == "total"){
    output <- totrich
  } else {
    if(metric == "appearance"){
      output <- appear
    } else {
      if(metric == "disappearance"){
        output <- disappear
      }
    }
  }
  
  # results
  return(output)
}


# generating appearances and disappearances for each experiment ---------------------------------------------
#make a new dataframe with just the label
spdata3 <- spdata2%>%
  mutate(id2=paste(project_name, treatment, sep='::'))

exp_code=spdata3%>%
  select(id2)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1)

for(i in 1:length(exp_code$id2)) {
  
  #creates a dataset for each unique trt, exp combo
  subset=spdata3[spdata3$id2==as.character(exp_code$id2[i]),]%>%
    select(id2, calendar_year, genus_species, abundance, plot_id)%>%
    group_by(id2, calendar_year, plot_id, genus_species)%>%
    summarise(abundance=max(abundance))%>%
    ungroup()
  
  #need this to keep track of experiment labels
  labels=subset%>%
    select(id2, calendar_year)%>%
    unique()
  
  #calculating appearances and disappearances
  appear<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', metric='appearance')
  disappear<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', metric='disappearance')
  total<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='abundance', metric='total')
  
  #merging back with labels to get back experiment labels
  turnover<-merge(appear, disappear, by=c('calendar_year'))
  turnoverAll<-merge(turnover, total, by=c('calendar_year'))
  turnoverLabel<-merge(turnoverAll, labels, by=c('calendar_year'), all=T)
  turnoverLabel[is.na(turnoverLabel)] <- 0
  
  #pasting into the dataframe made for this analysis
  for.analysis=rbind(turnoverLabel, for.analysis)  
}

for.analysis2 <- for.analysis%>%
  gather(key=variable, value=num_spp, appearance:total)%>%
  separate(col=id2, into=c('project_name', 'treatment'), sep='::')%>%
  mutate(num_spp2=ifelse(variable=='disappearance', -1*num_spp, num_spp))
  

#generating figure of appearances/disappearances for each experiment
pplots<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="pplots"&treatment!="N1P1"&treatment!="N1P2"&treatment!="N1P3"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("N1P0", "N2P0", "N2P1", "N2P2", "N2P3"), labels=c("N0 P0", "N10 P0","N10 P2.5","N10 P5","N10 P10"), values=c("blue", "red","purple","purple","purple"))+
  geom_line()+
  geom_vline(xintercept = 2004, linetype="longdash")+
  geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Phosphorus Plots")

BGP_ub<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="BGP unburned"&treatment!="u_u_p"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("u_u_c", "u_u_n", "u_u_b"), labels=c("N0 P0", "N10 P0", "N10 P1"), values=c("blue","red","purple"))+
  geom_line()+
  geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Belowground Plots Unburned")

BGP_b<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="BGP burned"&treatment!="b_u_p"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("b_u_c", "b_u_n", "b_u_b"), labels=c("N0 P0", "N10 P0", "N10 P1"), values=c("blue","red","purple"))+
  geom_line()+
  geom_vline(xintercept = 1989,linetype="longdash")+
  geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Belowground Plots Burned")

nutnet<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="P"&treatment!="K"&treatment!="PK"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("control", "N", "NP", "NK", "NPK"), labels=c("N0 P0 K0","N10 P0 K0" ,"N10 P10 K0","N10 P0 K10" ,"N10 P10 K10"), values=c("blue","red", "purple", "purple","purple"))+
  geom_line()+
  geom_vline(xintercept = 2010,linetype="longdash")+
  geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Nutrient Network")

invert<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="invert"&treatment!="caged"&treatment!="caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged"&treatment!="insecticide"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("control", "NPK", "NPK_insecticide"), labels=c("N0 P0 K0", "N10 P10 K10", "N10 P10 K10\n& Insecticide"), values=c("blue","purple","purple"))+
  geom_line()+
  geom_vline(xintercept = 2011,linetype="longdash")+
  geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Invertebrate Removal")

restoration<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="restoration"&treatment!='C'), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c('control', 'N'), labels=c("N0", "N5"), values=c("blue","red"))+
  geom_line()+
  geom_vline(xintercept = 2005, linetype="longdash")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Restoration Plots")

uk_annual<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="ukulinga annual"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("control", "N"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Change")+
  ggtitle("Ukulinga Annual Burn")

uk_four<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="ukulinga four"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("control", "N"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Change")+
  ggtitle("Ukulinga 4-Yr Burn")

uk_ub<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="ukulinga unburned"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("control", "N"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Change")+
  ggtitle("Ukulinga Unburned")

grid.arrange(pplots, BGP_ub, BGP_b, nutnet, invert, restoration, uk_annual, uk_four, uk_ub)
#export at 2400x1600






#generating figure of mean change by experiment year with all experiments included in one panel
experimentYear <- read.csv('experiment years.csv')

meanChangeTime <- mean_change%>%
  filter(project_name!='ChANGE'&project_name!='ghost fire')%>%
  merge(experimentYear, by=c('project_name', 'calendar_year'))

meanChangeTimeSubset <- meanChangeTime%>%
  filter(treatment=='b_u_n'|treatment=='b_u_b'|treatment=='u_u_n'|treatment=='u_u_b'|treatment=='N2P0'|treatment=='N2P1'|treatment=='N2P2'|treatment=='N2P3'|treatment=='N'|treatment=='NP'|treatment=='NK'|treatment=='NPK')%>%
  filter(project_name!='ukulinga annual'&project_name!='ukulinga four'&project_name!='ukulinga unburned'&project_name!='restoration')

meanChangeTimeSubset_nochange <- meanChangeTime%>%
  filter(project_name=='ukulinga annual'|project_name=='ukulinga four'|project_name=='ukulinga unburned'|project_name=='restoration')%>%
  filter(treatment!='C')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))


changePlot <- ggplot(data=subset(meanChangeTimeSubset,experiment_year>0), aes(x=experiment_year, y=mean_change, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment))) +
  xlab('Experiment Year') +
  ylab('Community Change') +
  scale_color_discrete(name='Experiment') +
  ylim(0,0.8) +
  xlim(0,30)

nochangePlot <- ggplot(data=meanChangeTimeSubset_nochange, aes(x=experiment_year, y=mean_change, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment))) +
  xlab('Experiment Year') +
  ylab('Community Change') +
  scale_color_discrete(name='Experiment') +
  ylim(0,0.8) +
  xlim(0,30)

pushViewport(viewport(layout=grid.layout(2,1)))
print(changePlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(nochangePlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))





###BGP unburned only NMDS spread
bgpUnburned<-spdata4%>%
  filter(project_name=='BGP unburned')%>%
  select(project_name, treatment, calendar_year, plot_id, abundance, species)%>%
  filter(treatment!="caged"&treatment!="caged_insecticide"&treatment!="fence"&treatment!="NPK_caged"&treatment!="NPK_caged_insecticide"&treatment!="NPKfence"&treatment!='u_u_b'&treatment!='u_u_p'&treatment!='b_u_b'&treatment!='b_u_p'&treatment!='insecticide'&treatment!='K'&treatment!='N1P1'&treatment!='N1P2'&treatment!='N1P3'&treatment!='N2P1'&treatment!='N2P2'&treatment!='N2P3'&treatment!='NK'&treatment!='NP'&treatment!='NPK_insecticide'&treatment!='P'&treatment!='PK')%>%
  mutate(treatment2=ifelse(treatment=='b_u_n', 'N', ifelse(treatment=='u_u_n', 'N', ifelse(treatment=='NPK', 'N', ifelse(treatment=='N', 'N', ifelse(treatment=='N2P0', 'N', 'control'))))))%>%
  select(-treatment)

bgpwide1<-bgpUnburned%>%
  spread(species, abundance, fill=0)

mds1<-metaMDS(bgpwide1[,5:115],autotransform=FALSE, shrink=FALSE, trymax = 1000)
mds1

###plotting this
theme_set(theme_bw(12))

plots<-as.data.frame(bgpwide1[,1:4])
scores1 <- data.frame(scores(mds1, display="sites"))  # Extracts NMDS scores for year "i" #
scores2a<- cbind(plots, scores1)


ggplot(scores2a, aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=8,aes(shape=treatment2))+
  scale_shape_manual(name="Treatment", values=c(16,17))+
  scale_color_manual(name="Year", values=c("purple","green","red","blue","orange","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("NMDS1")+
  ylab("NMDS2")









