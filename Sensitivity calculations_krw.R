### Calculating sensitivity of ANPP ###
### By: Kevin Wilcox (wilcoxkr@gmail.com) ###
### Last modified: Feb 27th, 2017 ###


# Prep workspace ----------------------------------------------------------

rm(list)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
#setwd("C:\\Users\\WilcoxKR.WILCOXKR-LAPTOP\\Dropbox\\Konza Nutrient Synthesis\\")


# Read in data ------------------------------------------------------------

# ID control plots #
controls<-read.csv("Konza_nutrient synthesis_spp comp.csv")%>%
  mutate(control=ifelse(treatment=="N1P0", 1, ifelse(treatment=="b_u_c", 1, ifelse(treatment==0, 1, ifelse(treatment=="control", 1, ifelse(treatment=="control_control", 1, ifelse(treatment=="u_u_c", 1, 0)))))))%>%
  select(project_name, treatment, control)%>%
  unique()

# Mean anpp by year #
anpp.means <- read.csv("annual means_anpp_with comm change info.csv")

# weather data #
climate <- read.csv("KNZ_weather data.csv")
climate <- climate[,1:3]

# Calculating sensitivity -------------------------------------------------

# Calculate average anpp in each treatment, experiment, and community type (pre, post) #
baselineANPP <- ddply(anpp.means, .(project_name, treatment, comm_type), summarize,
                      base.anpp = mean(anpp, na.rm=T))
anpp.means <- merge(anpp.means, baselineANPP, by=c("project_name", "treatment", "comm_type"), all=T)

# Create column for difference between current year anpp and average anpp #
anpp.means$anpp.diff <- with(anpp.means, (anpp-base.anpp))

# Combine anpp data with annual precipitation data #
anpp.ppt <- merge(anpp.means, climate, by.x="calendar_year",
                  by.y="YEAR", all.x=T)

# Calculate mean annual and growing season precipitation from 1983-2016 #
  # KW note Feb27,2017: I think a potential issue with taking the long term mean is that,
  # depending on which years a certain experiment/community type occurs in, you might
  # have a large difference in precipitation amount, but not much difference in anpp difference
  # simply because you are comparing precip to the long-term average and anpp to the average
  # within a smaller time period. I'm guessing this won't be an issue for experiments/community
  # types that exist over longer time periods, but it will likely be an issue for most pre-community
  # shifts... One solution to this may be to use the average precip for each exp/community type
  # for precipitation as well.
anpp.ppt$MAP <- mean(climate$AnnualPrecip,na.rm=T)
anpp.ppt$MGSP <- mean(climate$AprSeptPrecip,na.rm=T)
anpp.ppt$AP.diff <- with(anpp.ppt, (AnnualPrecip-MAP))
anpp.ppt$GSP.diff <- with(anpp.ppt, (AprSeptPrecip-MGSP))

# Calculate sensitivity
anpp.ppt$sensitivity <- with(anpp.ppt, (anpp.diff/AP.diff))
anpp.ppt$GSPsensitivity <- with(anpp.ppt, (anpp.diff/AP.diff))

# Identify wet and dry years (<25th percentile and >75th percentile) and exclude "average" years #
dry.cutoff <- quantile(climate$AnnualPrecip, probs=0.25, na.rm=T)
wet.cutoff <- quantile(climate$AnnualPrecip, probs=0.75, na.rm=T)
anpp.ppt$ppt_year <- with(anpp.ppt, 
                          (ifelse(AnnualPrecip>wet.cutoff,
                                            "wet", ifelse(AnnualPrecip<dry.cutoff,
                                                          "dry","average"))))
anpp.ppt.noavg <- subset(anpp.ppt, ppt_year != "average")

# Calculate means and st.err for each project, treatment, and community type #
sens.means <- ddply(anpp.ppt.noavg, .(project_name, comm_type),summarize,
               sensitivity = mean(sensitivity,na.rm=T),
               sens.se = sd(sensitivity,na.rm=T)/sqrt(length(sensitivity)))
sens.means <- subset(sens.means, comm_type=="a"|comm_type=="b"|comm_type=="n")
sens.means$sens.se[sens.means$sens.se=="Inf"] <- NA

# order factors for plotting #
sens.means$comm_type <- gsub("a","after",sens.means$comm_type)
sens.means$comm_type <- gsub("b","before",sens.means$comm_type)
sens.means$comm_type <- gsub("n","control",sens.means$comm_type)
sens.means$comm_type <- factor(sens.means$comm_type,levels=c("control","before","after"))

# plot overall senstiivity in control plots, before and after community change
sens.plot.all <- ggplot(sens.means, aes(x=comm_type, y=sensitivity))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=sensitivity-sens.se, ymax=sensitivity+sens.se))+
  facet_wrap(~project_name)+
  theme_few()

# split sensitivity into sensitivity to wet and dry years, plot error bars aren't working correctly # 
wet.dry.means <- ddply(anpp.ppt.noavg, .(project_name, comm_type, ppt_year),summarize,
               sensitivity = mean(sensitivity,na.rm=T),
               sens.se = sd(sensitivity)/sqrt(length(sensitivity)))
wet.dry.means <- subset(wet.dry.means, comm_type=="a"|comm_type=="b"|comm_type=="n")

wet.dry.means$comm_type <- gsub("a","after",wet.dry.means$comm_type)
wet.dry.means$comm_type <- gsub("b","before",wet.dry.means$comm_type)
wet.dry.means$comm_type <- gsub("n","control",wet.dry.means$comm_type)
wet.dry.means$comm_type <- factor(wet.dry.means$comm_type,levels=c("control","before","after"))

sens.plot.wetdry <- ggplot(wet.dry.means, aes(x=comm_type, y=sensitivity, fill=ppt_year, colour=ppt_year))+
  geom_bar(stat="identity", position="dodge")+
  geom_errorbar(position="dodge", aes(ymin=sensitivity-sens.se, ymax=sensitivity+sens.se, width=0))+
  facet_wrap(~project_name)+
  theme_few()

#ggsave("Figures//sensitivity by community type.png", sens.plot.all, width=5, height=4, units="in")
#ggsave("Figures//sensitivity by community type_wet and dry.png", sens.plot.wetdry, width=5, height=4, units="in")



