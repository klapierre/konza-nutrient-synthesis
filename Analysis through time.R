library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(gridExtra)
library(doBy)
library(grid)
library(codyn)

#meghan's:
setwd("~/Dropbox/Konza Nutrient Synthesis")

#kim's:
setwd('C:\\Users\\Kim\\Dropbox\\konza projects\\Konza Nutrient Synthesis\\Threshold Project\\data')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

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

##first get a list of controls
controls<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(control=ifelse(treatment=="N1P0"|treatment=="b_u_c"|treatment==0|treatment=="control"|treatment=="control_control"|treatment=="u_u_c", 1, 0))%>%
  select(project_name, treatment, control)%>%
  unique()

##select N only treatments
nitrogen<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(nitrogen=ifelse(treatment=="N2P0"|treatment=="b_u_n"|treatment==10|treatment=="NPK"|treatment=="N"|treatment=="control_nitrogen", 1, 0))%>%
  select(project_name, treatment, nitrogen)%>%
  unique()

#read in climate data
precip<-read.csv("WETDRY.csv")%>%
  mutate(calendar_year=YEAR, gprecip=AprSeptPrecip)%>%
  select(calendar_year, gprecip)

##do species data
spdata<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(id=paste(project_name, calendar_year, sep="::"))%>%
  #filter pre-treatment data
  filter(id!='ChANGE::2013'&id!='ghost fire::2014'&id!='nutnet::2007'&id!='pplots::2002'&id!='ukulinga annual::2005'&id!='ukulinga four::2005'&id!='ukulinga unburned::2005')

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

# ###looking at bgp only
# bgp_precip<-nit_anpp_precip%>%
#   filter(project_name=="BGP")%>%
#   mutate(shift=as.factor(ifelse(calendar_year>1991, 1, 0)))
# 
# ggplot(data=bgp_precip, aes(x=gprecip, y=PC_anpp, group=treatment))+
#     geom_point(aes(color=shift, shape=treatment), size=5)

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

###graphs for N treatments only
#pplots only as an example
  ggplot(data=subset(mean_change, project_name=="pplots"&treatment!='N1P0'&treatment!='N1P1'&treatment!='N1P2'&treatment!='N1P3'&treatment!='N2P1'&treatment!='N2P2'&treatment!='N2P3'), aes(x=as.numeric(calendar_year), y=mean_change))+
  geom_point(size=5)+
  geom_line()+
  geom_vline(xintercept = 2004, linetype="longdash")+
  geom_vline(xintercept = 2007)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Phosphorus Plots")+
  ylim(min=0, max=0.8)+
  xlim(min=2002.5, max=2015.5)
  #export at 1000x600



###evidence that multiple nutrient limitation does not limit change (just N alone is as much change as N plus other nutrients)
##pretty graphs
# red= n only
# blue=p only
# purple=n+p
#other trt=black

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
  ylim(min=0, max=0.8)+
  xlim(min=2002.5, max=2015.5)
#export at 1000x600

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
ggplot(data=subset(mean_change, project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!='K'&treatment!='PK'&treatment!='NK'), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",labels=c("N10P0","N10P10K0" ,"N10P10K10","N0P10K0"), values=c("red","purple", "purple","blue"))+
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
  scale_color_manual(name="Treatment",labels=c("N10P10K10", "Insecticide", "N10P10K10\n& Insecticide"), values=c("purple","black","purple"))+
  geom_line()+
  geom_vline(xintercept = 2011,linetype="longdash")+
  geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Community Change")+
  ggtitle("Invertebrate Removal")+
  ylim(min=0, max=0.8)

# restoration<-
# ggplot(data=subset(mean_change, project_name=="restoration"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
#   geom_point(aes(color=treatment), size=5)+
#   scale_color_manual(name="Treatment",labels=c("N5", "Carbon"), values=c("red","black"))+
#   geom_line()+
#   geom_vline(xintercept = 2005, linetype="longdash")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Change")+
#   ggtitle("Restoration Plots")+
#   ylim(min=0, max=0.8)
# 
# uk_annual<-
#   ggplot(data=subset(mean_change, project_name=="ukulinga annual"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
#   geom_point(aes(color=treatment), size=5)+
#   scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
#   geom_line()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Change")+
#   ggtitle("Ukulinga Annual Burn")+
#   ylim(min=0, max=0.8)
# 
# uk_four<-
#   ggplot(data=subset(mean_change, project_name=="ukulinga four"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
#   geom_point(aes(color=treatment), size=5)+
#   scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
#   geom_line()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Change")+
#   ggtitle("Ukulinga 4-Yr Burn")+
#   ylim(min=0, max=0.8)
# 
# uk_ub<-
#   ggplot(data=subset(mean_change, project_name=="ukulinga unburned"), aes(x=as.numeric(calendar_year), y=mean_change, group=treatment))+
#   geom_point(aes(color=treatment), size=5)+
#   scale_color_manual(name="Treatment",labels=c("N10"), values=c("red"))+
#   geom_line()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   xlab("Year")+
#   ylab("Community Change")+
#   ggtitle("Ukulinga Unburned")+
#   ylim(min=0, max=0.8)

grid.arrange(pplots, BGP_b, BGP_ub, nutnet, invert, ncol=3)
#export at 1600x600

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
means<-summaryBy(NMDS1+NMDS2~treatment+project_name, data=scores2, FUN=funcs)

ggplot(means, aes(x=NMDS1.mn, y=NMDS2.mn, color=treatment))+
  geom_point(size=8,aes(shape=project_name))+
  geom_errorbar(aes(ymin=NMDS2.mn-NMDS2.se, ymax=NMDS2.mn+NMDS2.se),color="black")+
  geom_errorbarh(aes(xmin=NMDS1.mn-NMDS1.se, xmax=NMDS1.mn+NMDS1.se),color="black")+
  scale_shape_manual(name="Experiment", values=c(15,16,4,17,18,11,9,10,13))+
  scale_color_manual(name="Treatment", values=c("purple","green","red","blue","purple","green","red","blue", "black","green","black","black","red","green","blue","blue","blue","red","purple","purple","purple","red","purple","purple","purple","purple","blue","blue"))+
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
  mutate(id2=paste(project_name, treatment, sep='::'))%>%
  filter(project_name!='ghost fire'&project_name!='restoration'&project_name!='ChANGE'&project_name!='ukulinga four')

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
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="pplots"&treatment!="N1P1"&treatment!="N1P2"&treatment!="N1P3"&treatment!="N2P1"&treatment!='N2P2'&treatment!='N2P3'), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("N1P0", "N2P0"), labels=c("control", "N10"), values=c("blue", "red"))+
  geom_line()+
  geom_vline(xintercept = 2004, linetype="longdash")+
  geom_vline(xintercept = 2006)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Phosphorus Plots")

BGP_ub<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="BGP unburned"&treatment!="u_u_p"&treatment!='u_u_b'), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("u_u_c", "u_u_n"), labels=c("N0", "N10"), values=c("blue","red"))+
  geom_line()+
  geom_vline(xintercept = 1989)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Belowground Plots Unburned")

BGP_b<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="BGP burned"&treatment!="b_u_p"&treatment!='b_u_b'), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("b_u_c", "b_u_n"), labels=c("N0", "N10"), values=c("blue","red"))+
  geom_line()+
  geom_vline(xintercept = 1989,linetype="longdash")+
  geom_vline(xintercept = 1994)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Belowground Plots Burned")

nutnet<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="nutnet"&treatment!="fence"&treatment!="NPKfence"&treatment!="P"&treatment!="K"&treatment!="PK"&treatment!="NP"&treatment!="PK"&treatment!="NK"&treatment!="NPK"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("control", "N"), labels=c("N0 P0 K0","N10"), values=c("blue","red"))+
  geom_line()+
  geom_vline(xintercept = 2010,linetype="longdash")+
  geom_vline(xintercept = 2012)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Nutrient Network")

invert<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="invert"&treatment!="caged"&treatment!="caged_insecticide"&treatment!="NPK_caged_insecticide"&treatment!="NPK_caged"&treatment!="insecticide"&treatment!="NPK_insecticide"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment",breaks=c("control", "NPK"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  geom_vline(xintercept = 2011,linetype="longdash")+
  geom_vline(xintercept = 2013)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Species Change")+
  ggtitle("Invertebrate Removal")

uk_annual<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="ukulinga annual"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("control", "N"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Change")+
  ggtitle("Ukulinga Annual Burn")

uk_ub<-
  ggplot(data=subset(for.analysis2, variable!='total'&project_name=="ukulinga unburned"), aes(x=as.numeric(calendar_year), y=num_spp2, group=interaction(treatment,variable)))+
  geom_point(aes(color=treatment), size=5)+
  scale_color_manual(name="Treatment", breaks=c("control", "N"), labels=c("N0", "N10"), values=c("blue", "red"))+
  geom_line()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Change")+
  ggtitle("Ukulinga Unburned")

grid.arrange(pplots, BGP_b, BGP_ub, nutnet, uk_annual, uk_ub, invert, ncol=3)
#export at 2400x1600






#generating figure of mean change by experiment year with all experiments included in one panel
experimentYear <- read.csv('experiment years.csv')

meanChangeTime <- mean_change%>%
  filter(project_name!='ChANGE'&project_name!='ghost fire')%>%
  merge(experimentYear, by=c('project_name', 'calendar_year'))

meanChangeTimeSubset <- meanChangeTime%>%
  filter(treatment=='b_u_n'|treatment=='u_u_n'|treatment=='N2P0'|treatment=='N'|treatment=='NPK')%>%
  filter(project_name!='ukulinga four'&project_name!='restoration')%>%
  mutate(remove=ifelse(treatment=='NPK'&project_name=='nutnet', 1, 0))%>%
  filter(remove!=1)


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))


ggplot(data=subset(meanChangeTimeSubset,experiment_year>0), aes(x=experiment_year, y=mean_change, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment))) +
  scale_color_manual(name="Experiment",labels=c("BGP burned", "BGP unburned", "Invert Removals", "NutNet", "PPlots", "Ukulinga burned", "Ukulinga unburned"), values=c("black", "green", "red", "orange", "blue", "dark grey", "purple")) +
  xlab('Experiment Year') +
  ylab('Community Change') +
  ylim(0,0.8) +
  xlim(0,30) +
  geom_vline(xintercept = 6)
#export at 1400x700

#burn regime
ggplot(data=subset(meanChangeTimeSubset,experiment_year>0), aes(x=experiment_year, y=mean_change, color=project_name)) +
  geom_point(size=5) +
  geom_line(aes(group=interaction(project_name,treatment))) +
  scale_color_manual(name="Experiment",labels=c("BGP burned", "BGP unburned", "Invert Removals", "NutNet", "PPlots", "Ukulinga burned", "Ukulinga unburned"), values=c("red", "brown", "orange", "orange", "orange", "red", "brown")) +
  xlab('Experiment Year') +
  ylab('Community Change') +
  ylim(0,0.8) +
  xlim(0,30) +
  geom_vline(xintercept = 6)
#export at 1400x700




###change in abundance
abund <- spdata4%>%
  left_join(experimentYear)%>%
  select(-X, -calendar_year, -id)%>%
  mutate(experiment_year_num=paste("yr", experiment_year, sep=''))%>%
  select(-experiment_year)%>%
  group_by(genus_species, plot_id, treatment, project_name, species, experiment_year_num)%>%
  summarise(abundance=mean(abundance))%>%
  ungroup()%>%
  group_by(genus_species, plot_id, treatment, project_name, species)%>%
  spread(key=experiment_year_num, value=abundance, fill=0)%>%
  ungroup()%>%
  mutate(change_yr2=yr2-yr1, change_yr3=yr3-yr1, change_yr4=yr4-yr1, change_yr5=yr5-yr1, change_yr6=yr6-yr1, change_yr7=yr7-yr1, change_yr8=yr8-yr1, change_yr9=yr9-yr1, change_yr10=yr10-yr1, change_yr11=yr11-yr1, change_yr12=yr12-yr1, change_yr13=yr13-yr1, change_yr14=yr14-yr1, change_yr15=yr15-yr1, change_yr20=yr20-yr1, change_yr25=yr25-yr1, change_yr30=yr30-yr1)

# yr3 <- ggplot(data=abund, aes(x=yr1,y=change_yr3)) +
#   geom_point() +
#   xlab('Initial Abundance') +
#   ylab('Abundance Change')
# 
# yr5 <- ggplot(data=abund, aes(x=yr1, y=change_yr5)) +
#   geom_point() +
#   xlab('Initial Abundance') +
#   ylab('Abundance Change')
# 
# yr7 <- ggplot(data=abund, aes(x=yr1, y=change_yr7)) +
#   geom_point() +
#   xlab('Initial Abundance') +
#   ylab('Abundance Change')
# 
# yr10 <- ggplot(data=abund, aes(x=yr1, y=change_yr10)) +
#   geom_point() +
#   xlab('Initial Abundance') +
#   ylab('Abundance Change')
# 
# grid.arrange(yr3, yr5, yr10, ncol=3)
# #print at 1000x500



#loss of dominant spp
domCtl <- spdata4%>%
  left_join(experimentYear)%>%
  select(-X, -calendar_year, -id)%>%
  filter(species=='andropogon gerardii'|species=='sorghastrum nutans'|species=='schizachyrium scoparium')%>%
  group_by(genus_species, treatment, project_name, species, experiment_year)%>%
  summarise(abundance=mean(abundance))%>%
  ungroup()%>%
  filter(treatment=='control'|treatment=='N1P0'|treatment=='b_u_c'|treatment=='u_u_c')%>%
  select(-treatment)
names(domCtl)[names(domCtl)=="abundance"] <- "abundance_ctl"

dom <- spdata4%>%
  left_join(experimentYear)%>%
  select(-X, -calendar_year, -id)%>%
  filter(species=='andropogon gerardii'|species=='sorghastrum nutans'|species=='schizachyrium scoparium')%>%
  group_by(genus_species, plot_id, treatment, project_name, species, experiment_year)%>%
  summarise(abundance=mean(abundance))%>%
  ungroup()%>%
  filter(treatment=='N'|treatment=='N2P0'|treatment=='NPK'|treatment=='b_u_n'|treatment=='u_u_n')%>%
  mutate(remove=ifelse(project_name=='nutnet'&treatment=='NPK', 1, 0))%>%
  filter(remove==0&experiment_year!=1)%>%
  merge(domCtl, by=c('genus_species', 'project_name', 'species', 'experiment_year'))%>%
  mutate(abund_change=log(abundance/abundance_ctl))

andro <- ggplot(barGraphStats(data=subset(dom, species=='andropogon gerardii'&project_name!='restoration'&project_name!='ukulinga four'), variable="abund_change", byFactorNames=c("experiment_year", "project_name")), aes(x=experiment_year, y=mean, group=project_name, color=project_name)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  theme(legend.position = 'none') +
  scale_color_manual(name="Experiment",labels=c("BGP burned", "BGP unburned", "Invert Removals", "NutNet", "PPlots", "Ukulinga burned", "Ukulinga unburned"), values=c("black", "green", "red", "orange", "blue", "dark grey", "purple")) +
  xlab('Experiment Year') +
  ylab('Change in Relative Abundance')

sorg <- ggplot(barGraphStats(data=subset(dom, species=='sorghastrum nutans'&project_name!='restoration'&project_name!='ukulinga four'), variable="abundance", byFactorNames=c("experiment_year", "project_name")), aes(x=experiment_year, y=mean, group=project_name, color=project_name)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  theme(legend.position = 'none') +
  scale_color_manual(name="Experiment",labels=c("BGP burned", "BGP unburned", "Invert Removals", "NutNet", "PPlots", "Ukulinga burned", "Ukulinga unburned"), values=c("black", "green", "red", "orange", "blue", "dark grey", "purple")) +
  xlab('Experiment Year') +
  ylab('Relative Abundance')

schiz <- ggplot(barGraphStats(data=subset(dom, species=='schizachyrium scoparium'&project_name!='restoration'&project_name!='ukulinga four'), variable="abundance", byFactorNames=c("experiment_year", "project_name")), aes(x=experiment_year, y=mean, group=project_name, color=project_name)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  theme(legend.position = 'none') +
  scale_color_manual(name="Experiment",labels=c("BGP burned", "BGP unburned", "Invert Removals", "NutNet", "PPlots", "Ukulinga burned", "Ukulinga unburned"), values=c("black", "green", "red", "orange", "blue", "dark grey", "purple")) +
  xlab('Experiment Year') +
  ylab('Relative Abundance')

grid.arrange(andro, sorg, schiz, ncol=3)
#export at 1000x500

#mean dom spp loss across all experiments
andro <- ggplot(barGraphStats(data=subset(dom, species=='andropogon gerardii'&project_name!='restoration'&project_name!='ukulinga four'), variable="abund_change", byFactorNames=c("experiment_year")), aes(x=experiment_year, y=mean)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  xlab('Experiment Year') +
  ylab('ln Change in Relative Abundance') +
  xlim(1,10) +
  ylim(-3,1) +
  geom_hline(yintercept=0) +
  geom_text()

sorg <- ggplot(barGraphStats(data=subset(dom, species=='sorghastrum nutans'&project_name!='restoration'&project_name!='ukulinga four'), variable="abund_change", byFactorNames=c("experiment_year")), aes(x=experiment_year, y=mean)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  xlab('Experiment Year') +
  ylab('ln Change in Relative Abundance') +
  xlim(1,10) +
  ylim(-3,1) +
  geom_hline(yintercept=0)

schiz <- ggplot(barGraphStats(data=subset(dom, species=='schizachyrium scoparium'&project_name!='restoration'&project_name!='ukulinga four'), variable="abund_change", byFactorNames=c("experiment_year")), aes(x=experiment_year, y=mean)) +
  geom_point(size=5) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  xlab('Experiment Year') +
  ylab('ln Change in Relative Abundance') +
  xlim(1,10) +
  ylim(-3,1) +
  geom_hline(yintercept=0)

grid.arrange(andro, sorg, schiz, ncol=3)
#export at 1000x500



###precip data
climate <- read.csv('Konza_Daily_Met_data_summary.csv')

hist(climate$Total_precip_mean, breaks=10)
hist(climate$Temp_mean, breaks=10)
hist(climate$Total_dry_days, breaks=10)


