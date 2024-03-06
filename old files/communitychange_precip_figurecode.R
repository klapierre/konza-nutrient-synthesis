setwd("C:/Users/Ellen/Desktop/Konza Synthesis Workshop")

comnut=read.csv(file="community_change_nut_experiments_withGSprecip.csv",header=TRUE,check.names=FALSE)
attach(comnut)
head(comnut)

unique(project_name)

par(mfrow=c(3,3))
#############################################################
#############################################################
#list trt levels
trt <- treatment[project_name=="pplots"]
unique(trt)
range(mean_change[project_name=="pplots"])

#####barplot
x_N1P1 <- calendar_year[project_name=="pplots" & treatment=="N1P1"]
y1_N1P1 <- mean_change[project_name=="pplots" & treatment=="N1P1"]
y2_N1P1 <- AprSeptPrecip[project_name=="pplots" & treatment=="N1P1"]

par(mar=c(5,4,4,5)+.1)

pr = barplot(y2_N1P1, border = NA, axes = FALSE, plot=FALSE, names.arg = x_N1P1)
##precip barplot
barplot(y2_N1P1,axes = FALSE,xlim=range(pr), main="Phosphorus plots")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900,1000))

par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2002,2015), ylim=c(0,0.6), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6))
box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation",side=4,line=3, cex=1)

####N2P0
x_N2P0 <- calendar_year[project_name=="pplots" & treatment=="N2P0"]
y1_N2P0 <- mean_change[project_name=="pplots" & treatment=="N2P0"]
y2_N2P0 <- AprSeptPrecip[project_name=="pplots" & treatment=="N2P0"]
lines(x_N2P0, y1_N2P0, xlim=range(x_N2P0), ylim=range(y1_N2P0), pch=16,col=1)
points (y1_N2P0~x_N2P0, cex = 2.5, col = 1, bg=1, pch=21)

####N2P1
x_N2P1 <- calendar_year[project_name=="pplots" & treatment=="N2P1"]
y1_N2P1 <- mean_change[project_name=="pplots" & treatment=="N2P1"]
y2_N2P1 <- AprSeptPrecip[project_name=="pplots" & treatment=="N2P1"]
lines(x_N2P1, y1_N2P1, xlim=range(x_N1P2), ylim=range(y1_N1P2), pch=16,col=1)
points (y1_N2P1~x_N2P1, pch = 21, cex = 2.5, col = 2, bg=2)

####N2P2
x_N2P2 <- calendar_year[project_name=="pplots" & treatment=="N2P2"]
y1_N2P2 <- mean_change[project_name=="pplots" & treatment=="N2P2"]
y2_N2P2 <- AprSeptPrecip[project_name=="pplots" & treatment=="N2P2"]
lines(x_N2P2, y1_N2P2, xlim=range(x_N2P2), ylim=range(y1_N2P2), pch=16,col=1)
points (y1_N2P2~x_N2P2, pch = 21, cex = 2.5, col =3, bg=3)

####N2P3
x_N2P3 <- calendar_year[project_name=="pplots" & treatment=="N2P3"]
y1_N2P3 <- mean_change[project_name=="pplots" & treatment=="N2P3"]
y2_N2P3 <- AprSeptPrecip[project_name=="pplots" & treatment=="N2P3"]
lines(x_N2P3, y1_N2P3, xlim=range(x_N2P3), ylim=range(y1_N2P3), pch=16,col=1)
points (y1_N2P3~x_N2P3, pch = 21, cex = 2.5, col =4, bg=4)

#legend("topleft",legend=c("N2P0","N2P1","N2P2","N2P3"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2,3,4), pt.bg=c(1,2,3,4,5,6,7))

##########################################################################################
############these trt dont have N###################################################################
#N1P1 points
lines(x_N1P1, y1_N1P1, xlim=range(x_N1P1), ylim=range(y1_N1P1), pch=16,col=1)
points (y1_N1P1~x_N1P1, pch = 21, cex = 2.5, col = 5, bg=5)

####N1P2
x_N1P2 <- calendar_year[project_name=="pplots" & treatment=="N1P2"]
y1_N1P2 <- mean_change[project_name=="pplots" & treatment=="N1P2"]
y2_N1P2 <- AprSeptPrecip[project_name=="pplots" & treatment=="N1P2"]
lines(x_N1P2, y1_N1P2, xlim=range(x_N1P2), ylim=range(y1_N1P2), pch=16,col=1)
points (y1_N1P2~x_N1P2, pch = 21, cex = 2.5, col =6, bg=6)

####N1P3
x_N1P3 <- calendar_year[project_name=="pplots" & treatment=="N1P3"]
y1_N1P3 <- mean_change[project_name=="pplots" & treatment=="N1P3"]
y2_N1P3 <- AprSeptPrecip[project_name=="pplots" & treatment=="N1P3"]
lines(x_N1P3, y1_N1P3, xlim=range(x_N1P3), ylim=range(y1_N1P3), pch=16,col=1)
points (y1_N1P3~x_N1P3, pch = 21, cex = 2.5, col = 7, bg=7)

legend("topleft",legend=c("N2P0","N2P1","N2P2","N2P3","N1P1","N1P2","N1P3"), bty="n", text.font=6, cex=1, pch=c(21), col=c(1,2,3,4,5,6,7), pt.bg=c(1,2,3,4,5,6,7))


################################
###############################################################
#########BGP burned#####################################################################
#list trt levels
trt <- treatment[project_name=="BGP burned"]
unique(trt)
range(mean_change[project_name=="BGP burned"])
range(calendar_year[project_name=="BGP burned"])

#####b_u_b
x_b_u_b <- calendar_year[project_name=="BGP burned" & treatment=="b_u_b"]
y1_b_u_b <- mean_change[project_name=="BGP burned" & treatment=="b_u_b"]
y2 <- c(478,433,402.7,569.1,1051,477.8,824.1,564.7,490.1,584.1,751.1,401.7,711.8,508.3,571.7,738.6,685.2,560.6,634.5,888.7,720.6,466,466,405.6,604.4,480.3,747.1)
yr<- c(1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

#par(mfrow=c(1,2))
par(mar=c(5,4,4,5)+.1)

pr = barplot(y2, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2,axes = FALSE,xlim=range(pr), main="Belowground plots burned")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900,1000))

#b_u_b points
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(1989,2015), ylim=c(0,0.8), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_b_u_b, y1_b_u_b, xlim=range(x_b_u_b), ylim=range(y1_b_u_b), pch=16,col=1)
points (y1_b_u_b~x_b_u_b, pch = 21, cex = 2.5, col = 1, bg=1)

####b_u_n
x_b_u_n <- calendar_year[project_name=="BGP burned" & treatment=="b_u_n"]
y1_b_u_n <- mean_change[project_name=="BGP burned" & treatment=="b_u_n"]
y2_b_u_n <- AprSeptPrecip[project_name=="BGP burned" & treatment=="b_u_n"]
lines(x_b_u_n, y1_b_u_n, xlim=range(x_b_u_n), ylim=range(y1_b_u_n), pch=16,col=1)
points (y1_b_u_n~x_b_u_n, pch = 21, cex = 2.5, col = 2, bg=2)

#legend("topleft",legend=c("bub","bun"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2), pt.bg=c(1,2))

#######################this trt doesn't have N
####b_u_p
x_b_u_p <- calendar_year[project_name=="BGP burned" & treatment=="b_u_p"]
y1_b_u_p <- mean_change[project_name=="BGP burned" & treatment=="b_u_p"]
y2_b_u_p <- AprSeptPrecip[project_name=="BGP burned" & treatment=="b_u_p"]
lines(x_b_u_p, y1_b_u_p, xlim=range(x_b_u_p), ylim=range(y1_b_u_p), pch=16,col=1)
points (y1_b_u_p~x_b_u_p, pch = 21, cex = 2.5, col = 3, bg=3)

legend("topleft",legend=c("bub","bun","bup"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2,3), pt.bg=c(1,2,3))

################################
###############################################################
#########BGP unburned#####################################################################
#list trt levels
head(comnut)
unique(project_name)
trt <- treatment[project_name=="BGP unburned"]
unique(trt)
range(mean_change[project_name=="BGP unburned"])
range(calendar_year[project_name=="BGP unburned"])

#####u_u_b
x_u_u_b <- calendar_year[project_name=="BGP unburned" & treatment=="u_u_b"]
y1_u_u_b <- mean_change[project_name=="BGP unburned" & treatment=="u_u_b"]
y2 <- c(478,433,402.7,569.1,1051,477.8,824.1,564.7,490.1,584.1,751.1,401.7,711.8,508.3,571.7,738.6,685.2,560.6,634.5,888.7,720.6,692.1,466,405.6,604.4,480.3,747.1)
yr<- c(1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)

#par(mfrow=c(1,2))
par(mar=c(5,4,4,5)+.1)

pr = barplot(y2, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2,axes = FALSE,xlim=range(pr),main="Belowground plots unburned")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900,1000))

#u_u_b points
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(1989,2015), ylim=c(0,0.8), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_u_u_b, y1_u_u_b, xlim=range(x_u_u_b), ylim=range(y1_u_u_b), pch=16,col=1)
points (y1_u_u_b~x_u_u_b, pch = 21, cex = 2.5, col = 1, bg=1)

####u_u_n
x_u_u_n <- calendar_year[project_name=="BGP burned" & treatment=="b_u_n"]
y1_u_u_n <- mean_change[project_name=="BGP burned" & treatment=="b_u_n"]
y2_u_u_n <- AprSeptPrecip[project_name=="BGP burned" & treatment=="b_u_n"]
lines(x_u_u_n, y1_u_u_n, xlim=range(x_u_u_n), ylim=range(y1_u_u_n), pch=16,col=1)
points (y1_u_u_n~x_u_u_n, pch = 21, cex = 2.5, col = 2, bg=2)

#legend("topleft",legend=c("uub","uun"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2), pt.bg=c(1,2))

#######################this trt doesn't have N
####u_u_p
x_u_u_p <- calendar_year[project_name=="BGP unburned" & treatment=="u_u_p"]
y1_u_u_p <- mean_change[project_name=="BGP unburned" & treatment=="u_u_p"]
y2_u_u_p <- AprSeptPrecip[project_name=="BGP unburned" & treatment=="u_u_p"]
lines(x_u_u_p, y1_u_u_p, xlim=range(x_u_u_p), ylim=range(y1_u_u_p), pch=16,col=1)
points (y1_u_u_p~x_u_u_p, pch = 21, cex = 2.5, col = 3, bg=3)

legend("topleft",legend=c("uub","uun","uup"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2,3), pt.bg=c(1,2,3))

#######################not needed now############################
##########
####u_u_p
x_u_u_p <- calendar_year[project_name=="BGP burned" & treatment=="u_u_p"]
y1_u_u_p <- mean_change[project_name=="BGP burned" & treatment=="u_u_p"]
y2_u_u_p <- AprSeptPrecip[project_name=="BGP burned" & treatment=="u_u_p"]
lines(x_u_u_p, y1_u_u_p, xlim=range(x_u_u_p), ylim=range(y1_u_u_p), pch=16,col=1)
points (y1_u_u_p~x_u_u_p, pch = 21, cex = 2.5, col = 3, bg=3)

legend("topleft",legend=c("uub","uun","uup"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2,3), pt.bg=c(1,2,3))
##################################################################
##############################################
#########################################################
#############NUTNET
#list trt levels
head(comnut)
unique(project_name)
trt <- treatment[project_name=="nutnet"]
unique(trt)
range(mean_change[project_name=="nutnet"])
range(AprSeptPrecip[project_name=="nutnet"])
range(calendar_year[project_name=="nutnet"])

#####N
x_N <- calendar_year[project_name=="nutnet" & treatment=="N"]
y1_N <- mean_change[project_name=="nutnet" & treatment=="N"]
y2_N <- AprSeptPrecip[project_name=="nutnet" & treatment=="N"]
yr<- c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016)

#par(mfrow=c(1,2))
par(mar=c(5,4,4,5)+.1)

pr = barplot(y2_N, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2_N,axes = FALSE,xlim=range(pr),main="Nutrient network")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))

#N points
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2007,2016), ylim=c(0,0.7), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_N, y1_N, xlim=range(x_N), ylim=range(y1_N), pch=16,col=1)
points (y1_N~x_N, pch = 21, cex = 2.5, col = 1, bg=1)

####NP
x_NP <- calendar_year[project_name=="nutnet" & treatment=="NP"]
y1_NP <- mean_change[project_name=="nutnet" & treatment=="NP"]
lines(x_NP, y1_NP, xlim=range(x_NP), ylim=range(y1_NP), pch=16,col=1)
points (y1_NP~x_NP, pch = 21, cex = 2.5, col = 2, bg=2)

####NK
x_NK <- calendar_year[project_name=="nutnet" & treatment=="NK"]
y1_NK <- mean_change[project_name=="nutnet" & treatment=="NK"]
lines(x_NK, y1_NK, xlim=range(x_NK), ylim=range(y1_NK), pch=16,col=1)
points (y1_NK~x_NK, pch = 21, cex = 2.5, col = 3, bg=3)

####NPK
x_NPK <- calendar_year[project_name=="nutnet" & treatment=="NPK"]
y1_NPK <- mean_change[project_name=="nutnet" & treatment=="NPK"]
lines(x_NPK, y1_NPK, xlim=range(x_NPK), ylim=range(y1_NPK), pch=16,col=1)
points (y1_NPK~x_NPK, pch = 21, cex = 2.5, col = 4, bg=4)

#legend("topleft",legend=c("N","NP","NK","NPK"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2,3,4), pt.bg=c(1,2,3,4))
#######################not needed now############################
##########
####K
x_K <- calendar_year[project_name=="nutnet" & treatment=="K"]
y1_K <- mean_change[project_name=="nutnet" & treatment=="K"]
lines(x_K, y1_K, xlim=range(x_K), ylim=range(y1_K), pch=16,col=1)
points (y1_K~x_K, pch = 21, cex = 2.5, col = 5, bg=5)

####P
x_P <- calendar_year[project_name=="nutnet" & treatment=="P"]
y1_P <- mean_change[project_name=="nutnet" & treatment=="P"]
lines(x_P, y1_P, xlim=range(x_P), ylim=range(y1_P), pch=16,col=1)
points (y1_P~x_P, pch = 21, cex = 2.5, col = 6, bg=6)

####PK
x_PK <- calendar_year[project_name=="nutnet" & treatment=="PK"]
y1_PK <- mean_change[project_name=="nutnet" & treatment=="PK"]
lines(x_PK, y1_PK, xlim=range(x_PK), ylim=range(y1_PK), pch=16,col=1)
points (y1_PK~x_PK, pch = 21, cex = 2.5, col = 7, bg=7)

####fence
x_fence <- calendar_year[project_name=="nutnet" & treatment=="fence"]
y1_fence <- mean_change[project_name=="nutnet" & treatment=="fence"]
lines(x_fence, y1_fence, xlim=range(x_fence), ylim=range(y1_fence), pch=16,col=1)
points (y1_fence~x_fence, pch = 21, cex = 2.5, col = 10, bg=10)

####NPKfence
x_NPKfence <- calendar_year[project_name=="nutnet" & treatment=="NPKfence"]
y1_NPKfence <- mean_change[project_name=="nutnet" & treatment=="NPKfence"]
lines(x_NPKfence, y1_NPKfence, xlim=range(x_NPKfence), ylim=range(y1_NPKfence), pch=16,col=1)
points (y1_NPKfence~x_NPKfence, pch = 21, cex = 2.5, col = 9, bg=9)

legend("topleft",legend=c("N","NP","NK","NPK","K","P"), bty="n", text.font=6, cex=1, pch=c(21), 
col=c(1,2,3,4,5,6), pt.bg=c(1,2,3,4,5,6))
legend(x=2008,y=0.74,legend=c("PK","fence","NPKfence"), bty="n", text.font=6, cex=1, pch=c(21), 
col=c(7,10,9), pt.bg=c(7,10,9))
##################################################################
##############################################
#########################################################
#############invert
trt <- treatment[project_name=="invert"]
unique(trt)
range(mean_change[project_name=="invert"])
range(AprSeptPrecip[project_name=="invert"])
range(calendar_year[project_name=="invert"])

#####NPK
x_NPK <- calendar_year[project_name=="invert" & treatment=="NPK"]
y1_NPK <- mean_change[project_name=="invert" & treatment=="NPK"]
y2 <- c(720.6,692.1,466,405.6,604.4,480.3,747.1,852)
yr<- c(2009,2010,2011,2012,2013,2014,2015,2016)

#par(mfrow=c(1,2))
par(mar=c(5,4,4,5)+.1)

pr = barplot(y2, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2,axes = FALSE,xlim=range(pr),main="Invertebrate removal")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))

#NPK points
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2009,2016), ylim=c(0,0.8), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2009,2010,2011,2012,2013,2014,2015,2016))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_NPK, y1_NPK, xlim=range(x_NPK), ylim=range(y1_NPK), pch=16,col=1)
points (y1_NPK~x_NPK, pch = 21, cex = 2.5, col = 1, bg=1)
#legend("topleft",legend=c("NPK"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1), pt.bg=c(1))

####NPK_caged
x_NPK_caged <- calendar_year[project_name=="invert" & treatment=="NPK_caged"]
y1_NPK_caged <- mean_change[project_name=="invert" & treatment=="NPK_caged"]
lines(x_NPK_caged, y1_NPK_caged, xlim=range(x_NPK_caged), ylim=range(y1_NPK_caged), pch=16,col=1)
points (y1_NPK_caged~x_NPK_caged, pch = 21, cex = 2.5, col = 2, bg=2)

####NPK_insecticide
x_NPK_insecticide <- calendar_year[project_name=="invert" & treatment=="NPK_insecticide"]
y1_NPK_insecticide <- mean_change[project_name=="invert" & treatment=="NPK_insecticide"]
lines(x_NPK_insecticide, y1_NPK_insecticide, xlim=range(x_NPK_insecticide), ylim=range(y1_NPK_insecticide), pch=16,col=1)
points (y1_NPK_insecticide~x_NPK_insecticide, pch = 21, cex = 2.5, col = 3, bg=3)

####NPK_caged_insecticide
x_NPK_caged_insecticide <- calendar_year[project_name=="invert" & treatment=="NPK_caged_insecticide"]
y1_NPK_caged_insecticide <- mean_change[project_name=="invert" & treatment=="NPK_caged_insecticide"]
lines(x_NPK_caged_insecticide, y1_NPK_caged_insecticide, xlim=range(x_NPK_caged_insecticide), ylim=range(y1_NPK_caged_insecticide), pch=16,col=1)
points (y1_NPK_caged_insecticide~x_NPK_caged_insecticide, pch = 21, cex = 2.5, col =4, bg=4)

####caged
x_caged <- calendar_year[project_name=="invert" & treatment=="caged"]
y1_caged <- mean_change[project_name=="invert" & treatment=="caged"]
lines(x_caged, y1_caged, xlim=range(x_caged), ylim=range(y1_caged), pch=16,col=1)
points (y1_caged~x_caged, pch = 21, cex = 2.5, col =5, bg=5)

####insecticide
x_insecticide <- calendar_year[project_name=="invert" & treatment=="insecticide"]
y1_insecticide <- mean_change[project_name=="invert" & treatment=="insecticide"]
lines(x_insecticide, y1_insecticide, xlim=range(x_insecticide), ylim=range(y1_insecticide), pch=16,col=1)
points (y1_insecticide~x_insecticide, pch = 21, cex = 2.5, col = 6, bg=6)

####caged_insecticide
x_caged_insecticide <- calendar_year[project_name=="invert" & treatment=="caged_insecticide"]
y1_caged_insecticide <- mean_change[project_name=="invert" & treatment=="caged_insecticide"]
lines(x_caged_insecticide, y1_caged_insecticide, xlim=range(x_caged_insecticide), ylim=range(y1_caged_insecticide), pch=16,col=1)
points (y1_caged_insecticide~x_caged_insecticide, pch = 21, cex = 2.5, col = 7, bg=7)

legend("topleft",legend=c("NPK","NPK caged","NPK insecticide","NPK caged insecticide","caged","insecticide","caged insecticide"), bty="n", text.font=6, cex=1, pch=c(21), 
col=c(1,2,3,4,5,6,7), pt.bg=c(1,2,3,4,5,6,7))
##################################################################
##############################################
#########################################################
#############restoration
trt <- treatment[project_name=="restoration"]
unique(trt)

range(mean_change[project_name=="restoration"])
range(AprSeptPrecip[project_name=="restoration"])
range(calendar_year[project_name=="restoration"])

#####N
x_Nr <- calendar_year[project_name=="restoration" & treatment=="N"]
y1_Nr <- mean_change[project_name=="restoration" & treatment=="N"]
y2r <- c(584.1,751.1,401.7,711.8,508.3,571.7,738.6,685.2,560.6,634.5,888.7,720.6,692.1,466,405.6)
yr<- c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012)

#par(mfrow=c(1,2))
par(mar=c(5,4,4,5)+.1)

pr = barplot(y2r, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2r,axes = FALSE,xlim=range(pr),main="Restoration plots")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))

#N
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(1998,2012), ylim=c(0,0.8), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012))
axis(side = 2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_Nr, y1_Nr, xlim=range(x_N), ylim=range(y1_N), pch=16,col=1)
points (y1_Nr~x_Nr, pch = 21, cex = 2.5, col = 1, bg=1)

####C
x_C <- calendar_year[project_name=="restoration" & treatment=="C"]
y1_C <- mean_change[project_name=="restoration" & treatment=="C"]
lines(x_C, y1_C, xlim=range(x_C), ylim=range(y1_C), pch=16,col=1)
points (y1_C~x_C, pch = 21, cex = 2.5, col = 2, bg=2)

legend("topleft",legend=c("N","C"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1,2), pt.bg=c(1,2))

##################################################################
##############################################
#########################################################
#############ukulinga annual
trt <- treatment[project_name=="ukulinga annual"]
unique(trt)
range(mean_change[project_name=="ukulinga annual"])
range(AprSeptPrecip[project_name=="ukulinga annual"])
range(calendar_year[project_name=="ukulinga annual"])

#####N
x_N <- calendar_year[project_name=="ukulinga annual" & treatment=="N"]
y1_N <- mean_change[project_name=="ukulinga annual" & treatment=="N"]
y2_N <- c(560.6,634.5,888.7,720.6,692.1)
yr<- c(2006,2007,2008,2009,2010)

pr = barplot(y2_N, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2_N,axes = FALSE,xlim=range(pr),main="Ukulinga annual burn")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))
#N
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2006,2010), ylim=c(0,0.35), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2006,2007,2008,2009,2010))
axis(side = 2, at = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_N, y1_N, xlim=range(x_N), ylim=range(y1_N), pch=16,col=1)
points (y1_N~x_N, pch = 21, cex = 2.5, col = 1, bg=1)

legend("topleft",legend=c("N"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1), pt.bg=c(1))

##############################################
#########################################################
#############ukulinga four
trt <- treatment[project_name=="ukulinga four"]
unique(trt)
range(mean_change[project_name=="ukulinga four"])
range(AprSeptPrecip[project_name=="ukulinga four"])
range(calendar_year[project_name=="ukulinga four"])

#####N
x_N <- calendar_year[project_name=="ukulinga four" & treatment=="N"]
y1_N <- mean_change[project_name=="ukulinga four" & treatment=="N"]
y2_N <- c(560.6,634.5,888.7,720.6,692.1)
yr<- c(2006,2007,2008,2009,2010)

pr = barplot(y2_N, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2_N,axes = FALSE,xlim=range(pr),main="Ukulinga 4yr burn")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))
#N
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2006,2010), ylim=c(0,0.2), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2006,2007,2008,2009,2010))
axis(side = 2, at = c(0,0.05,0.1,0.15,0.2))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_N, y1_N, xlim=range(x_N), ylim=range(y1_N), pch=16,col=1)
points (y1_N~x_N, pch = 21, cex = 2.5, col = 1, bg=1)

legend("topleft",legend=c("N"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1), pt.bg=c(1))

##############################################
#########################################################
#############ukulinga unburned
trt <- treatment[project_name=="ukulinga unburned"]
unique(trt)

#####N
x_N <- calendar_year[project_name=="ukulinga unburned" & treatment=="N"]
y1_N <- mean_change[project_name=="ukulinga unburned" & treatment=="N"]
y2_N <- c(560.6,634.5,888.7,720.6,692.1)
yr<- c(2006,2007,2008,2009,2010)

par(mar=c(5,4,4,5)+0.1)

pr = barplot(y2_N, border = NA, axes = FALSE, plot=FALSE, names.arg = yr)
##precip barplot
barplot(y2_N,axes = FALSE,xlim=range(pr),main="Ukulinga unburned")
axis(side = 4, at = c(0,100,200,300,400,500,600,700,800,900))
#N
par(new=TRUE)
plot(1, type="n", xlab="", ylab="", xlim=c(2006,2010), ylim=c(0,0.25), axes = FALSE, main="", cex.lab=1.6,cex=2.5, pch=21)
axis(side = 1, at = c(2006,2007,2008,2009,2010))
axis(side = 2, at = c(0,0.05,0.1,0.15,0.2,0.25))

box()
title(ylab="Community change", line=2, cex.lab=1.6)
mtext("Growing season precipitation (mm)",side=4,line=3, cex=1)
lines(x_N, y1_N, xlim=range(x_N), ylim=range(y1_N), pch=16,col=1)
points (y1_N~x_N, pch = 21, cex = 2.5, col = 1, bg=1)

legend("topleft",legend=c("N"), bty="n", text.font=6, cex=1.5, pch=c(21), col=c(1), pt.bg=c(1))

