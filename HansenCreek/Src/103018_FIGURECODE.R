rm(list = ls())

######################################################################################################################
######################################        Importing Data        #####################################
#####################################################################################################################


data.calc<- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Min_Nit_calc.csv")
colnames(data.calc) <- c("x", "ID", "Bank", "Distance", "Site", "Reach", "GW", "NetMin", "NetNit", "NO3Conc", "NH4Conc")
data.calc <- cbind(data.calc, log_dist=log(data.calc$Distance))[-69,]

data.exp<- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Figures/Calculatd_Data.csv")
data.exp <- cbind(data.exp, log_dist=log(data.exp$Distance))

data4 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Results/OrganicNitrogen.csv")
data4 <- cbind(data4, log_dist=log(data4$Distance))

data3 <- merge(data4, data.exp, by = c("Reach", "Site", "Distance", "Bank"))

data5 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Results/NO3.StableIsotope.2.csv")
data5 <- cbind(data5, log_dist=log(data5$Distance))

data6 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Results/NH4.StableIsotope.2.csv")
data6 <- cbind(data6, log_dist=log(data6$Distance))


setwd("~/Documents/Graduate Work 2/Hansen Creek")
rawData <- Hansen_bulkSoilIsotopes <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Hansen_bulkSoilIsotopes.csv")
rawData <- cbind(rawData, log_dist=log(rawData$Distance))
colNames <- names(rawData)

######################################################################################################################
######################################        Plot Functions          #####################################
#####################################################################################################################



FunctionRemoveReplactes <- function(dataFrame, groupingColumns, dataColumns, FUN) {
  
  colNames <- names(dataFrame)
  index <-  colNames %in% dataColumns
  
  listo <- vector("list", length(groupingColumns))
  for (i in 1:sum(index)) listo[[i]] <- dataFrame[,groupingColumns[i]]
  
  out <- dataFrame[,index]
  
  #names(data) <- c(groupingColumns, dataColumns)
  
  return (out)
}


######################################################################################################################
######################################        Average Replicates       #####################################
#####################################################################################################################

groupingVars1 <- c("Reach", "Site", "Bank", "Distance")
groupingVars2 <- c("Bank", "Distance")
dataVars <- c("pctN", "d15N", "pctC","d13C")
dataVars2 <- c("GW", "NetMin", "NetNit")
dataVars3 <- c("NH4_Conc", "N03_Conc")
dataVars4 <- c("OrgN", "massN", "GW")
dataVars5 <- c("d15N", "fit", "lwr", "upr")

names <- c(groupingVars1, dataVars)
names2 <- c("Reach", "Distance", "Site", "Bank", "GW", "NetMin", "NetNit")


out <- FunctionRemoveReplactes(rawData, groupingVars1, dataVars, mean)
data <- aggregate(out, by = list(rawData$Reach, rawData$Site, rawData$Bank, rawData$Distance, rawData$log_dist), mean)
colnames(data) <- c(groupingVars1, "log_dist", dataVars)


######################################################################################################################
###################################### Creating a dataset for each variable of interest #####################################
#####################################################################################################################
data.subsets <- c("Reach", "Site", "Bank", "Distance", "log_dist")

massN <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Results/MassNitrogen.csv") 
#massN <- cbind(massN, log_dist=log10(massN$Distance))

C.N <- (data$pctC/data$pctN)*(12.0107/14.0057)
data <- cbind(data, C.N)
data.d15N <- data[c(data.subsets,"d15N")]
data.d15N <- merge(data.d15N, massN, by = c("Reach", "Site", "Distance", "Bank"))
data.C.N <- data[c(data.subsets,"C.N")]
data.pctN2 <- data[c(data.subsets,"pctN")]
data.pctN <- data[c(data.subsets,"pctN")]
data.pctN <- merge(data.d15N, massN, by = c("Reach", "Site", "Distance", "Bank"))
data.d13C <- data[c(data.subsets,"d13C")]
data.GW <- data.exp[c(data.subsets,"GW")]
data.Min <- merge(data.exp[c(data.subsets,"NetMin", "GW")], data4[c(data.subsets,"OrgN")], by = data.subsets)
data.Nit <- na.omit(merge(data.exp[c(data.subsets,"NetNit", "GW")], data.exp[c(data.subsets,"NH4_Conc","GW")], by = c(data.subsets, "GW")))
data.NH4 <- na.omit(data.exp[c(data.subsets,"NH4_Conc","GW")])
data.NO3 <- na.omit(data.exp[c(data.subsets,"NO3_Conc", "GW")])
data.OrgN <- merge(data.exp[c(data.subsets, "GW")], data4[c(data.subsets,"OrgN")], by = data.subsets)
#data.massN <-data4[c(data.subsets,"massN", "GW")]
data.d15NNO3 <- na.omit(data5[c(data.subsets,"d15N", "N.ug")])
data.d15NNH4 <- na.omit(data6[c(data.subsets,"d15N")])
data.mass <- massN[c(data.subsets,"massN", "GW")]

data.NO3 <- merge(data.NO3, data.Nit, by = c("Reach", "Site", "Distance", "Bank"))
data.NH4 <- merge(data.NH4, data.Min, by = c("Reach", "Site", "Distance", "Bank"))



######################################################################################################################
######################################        Supported models       #####################################
#####################################################################################################################

library(nlme)
library(lme4)
library(lmerTest)
library(AICcmodavg)


model.d15N <- lm(d15N.x~Bank*log_dist+I(log_dist^2)*Bank, data.d15N)
model.d13C <- lm(d13C~log_dist, data.d13C)
model.pctN <- lm(pctN.x~Bank*log_dist, data.pctN)
model.Min <- lm(NetMin~OrgN, data.Min)
model.Nit <- lm(NetNit~Bank+NH4_Conc+GW, data.Nit)
model.NH4 <- lm(NH4_Conc~Bank, data.NH4)
model.NO3 <- lm(NO3_Conc~Bank+GW.x, data.NO3)
model.OrgN <- lm(OrgN~Bank*log_dist+I(log_dist^2)*Bank+GW, data.OrgN)
model.d15NNO3 <- lm(d15N~Bank*log_dist+N.ug, data.d15NNO3)
model.d15NNH4 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data.d15NNH4)
model.massN <- lm(massN~Bank*log_dist+GW, data.mass)
model.C.N <- lm(C.N~Bank*log_dist+I(log_dist^2)*Bank, data.C.N)
model.GW <- lm(GW~log_dist, data.GW)

par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(predict(model.Min), data.Min$NetMin, ylim = c(-2.5, 2.5), pch=19, cex=0.5)
plot(predict(model.Min), residuals(model.Min))
coefplot(model.Min, main = "Net Mineralization")

#qqnorm(resid(model.Min))
#qqline(resid(model.Min))
#shapiro.test(data.Min$NetMin)



data.C.N <- cbind(data.C.N, predict(model.C.N, data.C.N, interval="confidence"))
data.d15N <- cbind(data.d15N, predict(model.d15N, data.d15N, interval="confidence"))
data.d13C <- cbind(data.d13C, predict(model.d13C, data.d13C, interval="confidence"))
data.pctN <- cbind(data.pctN, predict(model.pctN, data.pctN, interval="confidence"))
data.GW <- cbind(data.GW, predict(model.GW, data.GW, interval="confidence"))
data.Min <- cbind(data.Min, predict(model.Min, data.Min, interval="confidence"))
data.Nit <- cbind(data.Nit, predict(model.Nit, data.Nit, interval="confidence"))
data.NH4 <- na.omit(cbind(data.NH4, predict(model.NH4, data.NH4, interval="confidence")))
data.NO3 <- cbind(data.NO3, predict(model.NO3, data.NO3, interval="confidence"))
data.OrgN <-  cbind(data.OrgN, predict(model.OrgN, data.OrgN, interval="confidence"))
data.d15NNO3 <- cbind(data.d15NNO3, predict(model.d15NNO3, data.d15NNO3, interval="confidence"))
data.d15NNH4 <- na.omit(cbind(data.d15NNH4, predict(model.d15NNH4, data.d15NNH4, interval="confidence")))
data.massN <- cbind(data.mass, predict(model.massN, data.mass, interval="confidence"))



######################################################################################################################
######################################        Calculate average and SD for       #####################################
######################################    each bank and distance across sites    #####################################
#####################################################################################################################
dataVars5 <- c("fit", "lwr", "upr")

mean.d15N <- aggregate(data.d15N[c(dataVars5, "d15N.x")], by= list(Bank = data.d15N$Bank, Distance= data.d15N$Distance), FUN=mean)
mean.d13C <- aggregate(data.d13C[c(dataVars5, "d13C")], by= list(Bank = data.d13C$Bank, Distance= data.d13C$Distance), FUN=mean)
mean.pctN <- aggregate(data.pctN[c(dataVars5, "pctN.x")], by= list(Bank = data.pctN$Bank, Distance= data.pctN$Distance), FUN=mean)
mean.C.N <- aggregate(data.C.N[c(dataVars5, "C.N")], by= list(Bank = data.C.N$Bank, Distance= data.C.N$Distance), FUN=mean)

mean.GW <- aggregate(data.GW[c(dataVars5,"GW")], by= list(Bank = data.GW$Bank, Distance= data.GW$Distance), FUN=mean)

#pred.GW <- 
mean.Min <- aggregate(data.Min[c(dataVars5, "NetMin")], by= list(Bank = data.Min$Bank, Distance= data.Min$Distance), FUN=mean)
mean.Nit <- aggregate(data.Nit[c(dataVars5, "NetNit")], by= list(Bank = data.Nit$Bank, Distance= data.Nit$Distance), FUN=mean)
mean.NH4 <- aggregate(data.NH4[c(dataVars5, "NH4_Conc")], by= list(Bank = data.NH4$Bank, Distance= data.NH4$Distance), FUN=mean)
mean.NO3 <- aggregate(data.NO3[c(dataVars5, "NO3_Conc")], by= list(Bank = data.NO3$Bank, Distance= data.NO3$Distance), FUN=mean)
mean.OrgN <- aggregate(data.OrgN[c(dataVars5, "OrgN")], by= list(Bank = data.OrgN$Bank, Distance= data.OrgN$Distance), FUN=mean)
mean.d15NNO3 <- aggregate(data.d15NNO3[c(dataVars5, "d15N")], by= list(Bank = data.d15NNO3$Bank, Distance= data.d15NNO3$Distance), FUN=mean)
mean.d15NNH4 <-aggregate(data.d15NNH4[c(dataVars5, "d15N")], by= list(Bank = data.d15NNH4$Bank, Distance= data.d15NNH4$Distance), FUN=mean)
mean.massN <-aggregate(data.massN[c(dataVars5, "massN")], by= list(Bank = data.massN$Bank, Distance= data.massN$Distance), FUN=mean)


sd.d15N <- aggregate(data.d15N[c(dataVars5, "d15N.x")], by= list(Bank = data.d15N$Bank, Distance= data.d15N$Distance), FUN=sd)
sd.d13C <- aggregate(data.d13C[c(dataVars5, "d13C")], by= list(Bank = data.d13C$Bank, Distance= data.d13C$Distance), FUN=sd)
sd.pctN <- aggregate(data.pctN[c(dataVars5, "pctN.x")], by= list(Bank = data.pctN$Bank, Distance= data.pctN$Distance), FUN=sd)

sd.Min <- aggregate(data.Min[c(dataVars5, "NetMin")], by= list(Bank = data.Min$Bank, Distance= data.Min$Distance), FUN=sd)
sd.Nit <- aggregate(data.Nit[c(dataVars5, "NetNit")], by= list(Bank = data.Nit$Bank, Distance= data.Nit$Distance), FUN=sd)
sd.NH4 <- aggregate(data.NH4[c(dataVars5, "NH4_Conc")], by= list(Bank = data.NH4$Bank, Distance= data.NH4$Distance), FUN=sd)
sd.NO3 <- aggregate(data.NO3[c(dataVars5, "NO3_Conc")], by= list(Bank = data.NO3$Bank, Distance= data.NO3$Distance), FUN=sd)
sd.OrgN <- aggregate(data.OrgN[c(dataVars5, "OrgN")], by= list(Bank = data.OrgN$Bank, Distance= data.OrgN$Distance), FUN=sd)
sd.d15NNO3 <- aggregate(data.d15NNO3[c(dataVars5, "d15N")], by= list(Bank = data.d15NNO3$Bank, Distance= data.d15NNO3$Distance), FUN=sd)
sd.d15NNH4 <- aggregate(data.d15NNH4[c(dataVars5, "d15N")], by= list(Bank = data.d15NNH4$Bank, Distance= data.d15NNH4$Distance), FUN=sd)
sd.massN <-aggregate(data.massN[c(dataVars5, "massN")], by= list(Bank = data.massN$Bank, Distance= data.massN$Distance), FUN=sd)
sd.C.N <- aggregate(data.C.N[c(dataVars5, "C.N")], by= list(Bank = data.C.N$Bank, Distance= data.C.N$Distance), FUN=sd)
sd.GW <- aggregate(data.GW[c(dataVars5, "GW")], by= list(Bank = data.GW$Bank, Distance= data.GW$Distance), FUN=sd)

######################################################################################################################
######################################        Defining axes and river scales       ##################################
#####################################################################################################################

#define y axes 
axd15N <- seq(2,12,2)
axd13C <- seq(-29, -25, 1)
axpctN <- seq(0.5, 3.5, 0.5)
axGW <- seq(0, 6, 1)
axMin <- seq(0, 10, 2)
axNit <- seq(0, 10, 2)
axNO3 <- seq(0, 25, 5)
axNH4 <- seq(0, 300, 50)
axOrgN <- seq(1000, 7000, 1000)
axd15NNO3 <- seq(-14, 10, 4)
axd15NNH4 <- seq(0, 40, 10)
axdmassN <- seq(0.0, 6, 1)
axdC.N <- seq(10, 18, 2)



stream.0 <- 0
stream.1 <- 1
labels.1 <- c("Salmon Enhanced", "Salmon Depleted")
labels.0 <- c("", "")
xaxis.0 <- c("", 50)
xaxis.1 <- c("Meters", 1,3,6,10,20)

#Make plots
mat <- matrix(data = c(1,0, 2, 3, 0, 4, 
                       5, 0, 6, 7, 0, 8,
                       9, 0, 10, 11, 0, 12, 
                       13, 0, 14, 15, 0, 16,
                       17, 0, 18, 19, 0, 20,
                       21, 22, 23, 24, 25, 26), ncol=6, nrow=6, byrow=TRUE)
pdf(file="HansenSoil_Plots2.pdf", width=11, height=12)
lay <- layout(mat=mat, widths=c(6.25,1,5, 6.25,1,5,
                                6.25,1,5, 6.25,1,5, 
                                6.25,1,5, 6.25,1,5,
                                6.25,1,5, 6.25,1,5,
                                6.25,1,5, 6.25,1,5,
                                6.25,1,5, 6.25,1,5),
              heights=c(3,3, 3, 3,3, 1.5))


#layout.show(26)

#par(mfrow=c(1,1))
xlim.R<-c(0,25)
xlim.L<-c(25,0)

######################################################################################################################
######################################        d15N      ##################################
#####################################################################################################################
par(mar=c(0,5,1,0), oma=rep(0,4))

tempData <- subset(data.d15N, Bank =="Left")
tempData.1 <-subset(mean.d15N, Bank =="Left") 
tempDataSD <- subset(sd.d15N, Bank =="Left")

plot(x=tempData$Distance, y=tempData$d15N.x, type="p", pch = 20, cex=1.25, col='gray58', 
      bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N", sep="")),
     xlab="", xlim=xlim.L, ylim<- c(1, 14), lty=stream.0, cex.lab=1.25)



points(x=tempData$Distance, y=tempData$d15N.x, pch = 20, cex=1.25, col='gray58')
lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
mtext("a.", 3, line=-1, adj=0, font=2)
axis(2, at=axd15N, cex=2, line=-1)
axis(1, at=xaxis.0[2:6])
text(x=10, y=13, labels=labels.1[1], cex=1.5, font=2)
mtext("per mil vs. air", 2, cex=0.7, line = 1.5)

par(mar=c(0,0,1,0))
tempData <- subset(data.d15N, Bank =="Right")
tempData.1 <-subset(mean.d15N, Bank =="Right")
tempDataSD <- subset(sd.d15N, Bank =="Right")
 

plot(x=tempData$Distance, y=tempData$d15N.x, pch = 20, col='gray58',type="p",  cex=1.25, 
       bty="n", xaxs="i", axes=F, ylab="", xlab="", xlim=xlim.R,ylim<- c(1, 14),lty=stream.0)
lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])
text(x=10, y=13, labels=labels.1[2], cex=1.5, font=2)

######################################################################################################################
######################################        d13C   ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.d13C, Bank =="Left")
tempData.1 <-subset(mean.d13C, Bank =="Left") 
tempDataSD <- subset(sd.d13C, Bank =="Left")

plot(x=tempData$Distance, y=tempData$d13C, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= expression(paste(delta^13, "C", sep="")),  xlab="", 
     xlim=xlim.L, ylim = c(-29.25, -24.5), lty=stream.0, cex.lab=1.25)
mtext("b.", 3, line=-1, adj=0, font=2)
mtext("per mil vs. VPDB", 2, cex=0.7, line = 1.5)


lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axd13C, line=-1)
axis(1, at=xaxis.0[2:6])
text(x=10, y=-24.85, labels=labels.1[1], cex=1.5, font=2)

par(mar=c(0,0,1,0))

tempData <- subset(data.d13C, Bank =="Right")
tempData.1 <-subset(mean.d13C, Bank =="Right")
tempDataSD <- subset(sd.d13C, Bank =="Right")

plot(x=tempData$Distance, y=tempData$d13C,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab="",
     xlab="", xlim=xlim.R, ylim = c(-29.25, -24.5), lty=stream.0)



lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])
axis(1, at=xaxis.0[2:6])
text(x=10, y=-24.85, labels=labels.1[2], cex=1.5, font=2)

######################################################################################################################
######################################        d15N NH4   ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.d15NNH4, Bank =="Left")
tempData.1 <-subset(mean.d15NNH4, Bank =="Left") 
tempDataSD <- subset(sd.d15NNH4, Bank =="Left")

plot(x=tempData$Distance, y=tempData$d15N, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= expression(paste(delta^15, "N ", NH[4]^" +")),  xlab="", 
     xlim=xlim.L, ylim = c(-1, 47), lty=stream.0, cex.lab=1.25)
mtext("c.", 3, line=-1, adj=0, font=2)
mtext(expression(paste("(per mil vs. air)")), 2, cex=0.7, line = 1.5)


lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axd15NNH4, line=-1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.d15NNH4, Bank =="Right")
tempData.1 <-subset(mean.d15NNH4, Bank =="Right")
tempDataSD <- subset(sd.d15NNH4, Bank =="Right")

plot(x=tempData$Distance, y=tempData$d15N,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab="",
     xlab="", xlim=xlim.R, ylim = c(-1, 47), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])
axis(1, at=xaxis.0[2:6])





######################################################################################################################
######################################        NH4      ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.NH4, Bank =="Left")
tempData.1 <-subset(mean.NH4, Bank =="Left") 
tempDataSD <- subset(sd.NH4, Bank =="Left")

plot(x=tempData$Distance, y=tempData$NH4_Conc, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= expression(paste(NH[4]^" +")),  xlab="", 
     xlim=xlim.L, ylim = c(-20,350), lty=stream.0, cex.lab=1.25)
mtext("d.", 3, line=-1, adj=0, font=2)
mtext(expression(paste( mu,"g  N ", g^-1, " dry soil")), 2, cex=0.7, line = 1.5)


lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axNH4, line=-1)
axis(1, at=xaxis.0[2:6])
text(x=10, y=75, labels=labels.0[1], cex=1.5, font=2)

par(mar=c(0,0,1,0))

tempData <- subset(data.NH4, Bank =="Right")
tempData.1 <-subset(mean.NH4, Bank =="Right")
tempDataSD <- subset(sd.NH4, Bank =="Right")

plot(x=tempData$Distance, y=tempData$NH4_Conc,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N (per mil vs. air)", sep="")),
     xlab="", xlim=xlim.R, ylim = c(-20,350), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])
text(x=10, y=12700, labels=labels.0[2], cex=1.5, font=2)

legend(x=1, y=300, legend=c("Data", "Predicted Value", "95% Confidence Interval"), 
       col=c("gray58", "gray8", "gray8"), lty= c(0, 1, 5), pch=c(16,1,1), bty="o")

######################################################################################################################
######################################        NO3      ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.NO3, Bank =="Left")
tempData.1 <-subset(mean.NO3, Bank =="Left") 
tempDataSD <- subset(sd.NO3, Bank =="Left")

plot(x=tempData$Distance, y=tempData$NO3_Conc, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= expression(paste(NO[3]^"- ")),  xlab="", 
     xlim=xlim.L, ylim = c(-5,30), lty=stream.0, cex.lab=1.25)

mtext("e.", 3, line=-1, adj=0, font=2)
lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axNO3, line=-1)
axis(1, at=xaxis.0[2:6])
mtext(expression(paste( mu,"g  N ", g^-1, " dry soil")), 2, cex=0.7, line = 1.5)

par(mar=c(0,0,1,0))

tempData <- subset(data.NO3, Bank =="Right")
tempData.1 <-subset(mean.NO3, Bank =="Right")
tempDataSD <- subset(sd.NO3, Bank =="Right")

plot(x=tempData$Distance, y=tempData$NO3_Conc,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N (per mil vs. air)", sep="")),
     xlab="", xlim=xlim.R, ylim = c(-5,30), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])


######################################################################################################################
######################################        NetMin     ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.Min, Bank =="Left")
tempData.1 <-subset(mean.Min, Bank =="Left") 
tempDataSD <- subset(sd.Min, Bank =="Left")

plot(x=tempData$Distance, y=tempData$NetMin, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= "Net Mineralization",  xlab="", 
     xlim=xlim.L, ylim = c(-1,11), lty=stream.0, cex.lab=1.25)
mtext("f.", 3, line=-1, adj=0, font=2)
mtext(expression(paste(mu,"g  N ", g^-1, " dry soil ", d^-1)), 2, cex=0.7, line = 1.5)

lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axMin, line=-1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.Min, Bank =="Right")
tempData.1 <-subset(mean.Min, Bank =="Right")
tempDataSD <- subset(sd.Min, Bank =="Right")

plot(x=tempData$Distance, y=tempData$NetMin,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab="",
     xlab="", xlim=xlim.R, ylim = c(-1,11), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])



######################################################################################################################
######################################        NetNit   ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.Nit, Bank =="Left")
tempData.1 <-subset(mean.Nit, Bank =="Left") 
tempDataSD <- subset(sd.Nit, Bank =="Left")

plot(x=tempData$Distance, y=tempData$NetNit, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= "Net Nitrification ",  xlab="", 
     xlim=xlim.L, ylim = c(-1,11), lty=stream.0, cex.lab=1.25)
mtext("g.", 3, line=-1, adj=0, font=2)
mtext(expression(paste(mu,"g  N ", g^-1, " dry soil ", d^-1)), 2, cex=0.7, line = 1.5)

lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axMin, line=-1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.Nit, Bank =="Right")
tempData.1 <-subset(mean.Nit, Bank =="Right")
tempDataSD <- subset(sd.Nit, Bank =="Right")

plot(x=tempData$Distance, y=tempData$NetNit,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab="",
     xlab="", xlim=xlim.R, ylim = c(-1,11), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])




######################################################################################################################
######################################        OrgN     ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.OrgN, Bank =="Left")
tempData.1 <-subset(mean.OrgN, Bank =="Left") 
tempDataSD <- subset(sd.OrgN, Bank =="Left")
ylab <- expression(paste("Organic Nitrogen ", mu,"g  N ", mg^-1, " dry soil"))


plot(x=tempData$Distance, y=tempData$OrgN, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= "Organic Nitrogen",  xlab="", 
     xlim=xlim.L, ylim = c(0,8000), lty=stream.0, cex.lab=1.25)
mtext("h.", 3, line=-1, adj=0, font=2)
mtext(expression(paste(mu,"g  N ", mg^-1, " dry soil")), 2, cex=0.7, line = 1.5)

lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axOrgN, line = -1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.OrgN, Bank =="Right")
tempData.1 <-subset(mean.OrgN, Bank =="Right")
tempDataSD <- subset(sd.OrgN, Bank =="Right")

plot(x=tempData$Distance, y=tempData$OrgN,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N (per mil vs. air)", sep="")),
     xlab="", xlim=xlim.R, ylim = c(0,8000), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])

######################################################################################################################
######################################     GW    ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.GW, Bank =="Left")
tempData.1 <-subset(mean.GW, Bank =="Left") 
tempDataSD <- subset(sd.GW, Bank =="Left")


axGW
plot(x=tempData$Distance, y=tempData$GW, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= "Gravimetric Water",  xlab="", 
     xlim=xlim.L, ylim = c(0,7), lty=stream.0, cex.lab=1.25)
mtext("i.", 3, line=-1, adj=0, font=2)
mtext("Content", 2, cex=0.85, line = 1.5)

lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axGW, line = -1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.GW, Bank =="Right")
tempData.1 <-subset(mean.GW, Bank =="Right")
tempDataSD <- subset(sd.GW, Bank =="Right")

plot(x=tempData$Distance, y=tempData$GW,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N (per mil vs. air)", sep="")),
     xlab="", xlim=xlim.R, ylim = c(0,7), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])

######################################################################################################################
######################################       C:N     ##################################
#####################################################################################################################
par(mar=c(0,5,1,0))

tempData <- subset(data.C.N, Bank =="Left")
tempData.1 <-subset(mean.C.N, Bank =="Left") 
tempDataSD <- subset(sd.C.N, Bank =="Left")



plot(x=tempData$Distance, y=tempData$C.N, type='p',pch = 20, cex=1.25, col='gray58',bty="n", xaxs="i", axes=F,
     ylab= "C:N",  xlab="", 
     xlim=xlim.L, ylim = c(9,20), lty=stream.0, cex.lab=1.25)
mtext("j.", 3, line=-1, adj=0, font=2)


lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col= 'gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type= "n" , lty=5, errbar.col = 'gray8')
axis(2, at=axdC.N, line = -1)
axis(1, at=xaxis.0[2:6])


par(mar=c(0,0,1,0))

tempData <- subset(data.C.N, Bank =="Right")
tempData.1 <-subset(mean.C.N, Bank =="Right")
tempDataSD <- subset(sd.C.N, Bank =="Right")

plot(x=tempData$Distance, y=tempData$C.N,type="p", pch = 20, cex=1.25, col='gray58', 
     bty="n", xaxs="i", axes=F,  ylab=expression(paste(delta^15, "N (per mil vs. air)", sep="")),
     xlab="", xlim=xlim.R, ylim = c(9,20), lty=stream.0)




lines(sort(tempData.1$Distance), tempData.1$fit[order(tempData.1$Distance)], col='gray8', type='o')
errbar(x=tempData.1$Distance, y=tempData.1$fit, yplus=tempData.1$upr, 
       yminus=tempData.1$lwr, add=T, type="n", lty=5, errbar.col = 'gray8')
axis(1, at=xaxis.0[2:6])


######################################################################################################################
######################################        River left   ##################################
#####################################################################################################################
par(mar=c(4,5,0,0))

plot(x=seq(0,20,1), y=seq(0,20,1), xlim = xlim.L, ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab=xaxis.1[1], col="brown", lwd=3, cex.lab=1.25, font.lab=2)
axis(1, at=xaxis.1[2:6])

par(mar=c(4,0,0,0))

plot(x=c(0,1,5,6), y=c(0,-20,-20,0), xlim = c(0,6), ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab="", col="brown", lwd=3)
axis(1, at=xaxis.0[2:6])
lines(x=c(1,5), y=c(-20,-20), col="blue", lwd=4)


par(mar=c(4,0,0,0))

plot(x=seq(0,20,1), y=seq(0,20,1), xlim = xlim.R, ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab=xaxis.1[1], col="brown", lwd=3, cex.lab=1.25, font.lab=2)
axis(1, at=xaxis.1[2:6])

par(mar=c(4,5,0,0))

plot(x=seq(0,20,1), y=seq(0,20,1), xlim = xlim.L, ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab=xaxis.1[1], col="brown", lwd=3, cex.lab=1.25, font.lab=2)
axis(1, at=xaxis.1[2:6])

par(mar=c(4,0,0,0))

plot(x=c(0,1,5,6), y=c(0,-20,-20,0), xlim = c(0,6), ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab="", col="brown", lwd=3)
axis(1, at=xaxis.0[2:6])
lines(x=c(1,5), y=c(-20,-20), col="blue", lwd=4)


par(mar=c(4,0,0,0))

plot(x=seq(0,20,1), y=seq(0,20,1), xlim = xlim.R, ylim=c(-20,30), type='l', bty="n", xaxs="i",
     axes=F, ylab="", xlab=xaxis.1[1], col="brown", lwd=3, cex.lab=1.25, font.lab=2)
axis(1, at=xaxis.1[2:6])


dev.off()


######################################################################################################################
######################################        Supplementary Material  ##################################
#####################################################################################################################
library(dotwhisker)
library(broom)
library(dplyr)

pdf(file="SuppMat.pdf", width=11, height=12)

######################################   d15N
par(mfrow=c(5,4), mar=c(5,5,2,2))
plot(predict(model.d15N), data.d15N$d15N.x, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="d15N")

plot(predict(model.d15N), residuals(model.d15N), ylab = "Residuals", xlab = "Predicted",
     main="d15N",  pch=19, cex=0.5)
######################################   d13C
plot(predict(model.d13C), data.d13C$d13C, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="d13C")

plot(predict(model.d13C), residuals(model.d13C), ylab = "Residuals", xlab = "Predicted",
     main="d13C",  pch=19, cex=0.5, ylim=c(-2.5, 2.5))

######################################   d15NNH4
plot(predict(model.d15NNH4), data.d15NNH4$d15N, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="d15N NH4")

plot(predict(model.d15NNH4), residuals(model.d15NNH4), ylab = "Residuals", xlab = "Predicted",
     main="d15N NH4",  pch=19, cex=0.5)


######################################   NH4
plot(predict(model.NH4), data.NH4$NH4Conc, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="NH4")

plot(predict(model.NH4), residuals(model.NH4), ylab = "Residuals", xlab = "Predicted",
     main="NH4",  pch=19, cex=0.5)

######################################   NO3
plot(predict(model.NO3), data.NO3$NO3Conc, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="NO3")

plot(predict(model.NH4), residuals(model.NH4), ylab = "Residuals", xlab = "Predicted",
     main="NO3",  pch=19, cex=0.5)

######################################   NetMin
plot(predict(model.Min), data.Min$NetMin,, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="Net Mineralization")

plot(predict(model.Min), residuals(model.Min), ylab = "Residuals", xlab = "Predicted",
     main="Net Mineralization",  pch=19, cex=0.5, ylim=c(-2.5,2.5))

######################################   NetNit
plot(predict(model.Nit), data.Nit$NetNit,, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="Net Nitrification")

plot(predict(model.Nit), residuals(model.Nit),ylab = "Residuals", xlab = "Predicted",
     main="Net Nitrification",  pch=19, cex=0.5, ylim=c(-2.5,2.5))

######################################  OrgN

plot(predict(model.OrgN), data.OrgN$OrgN, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="Organic Nitrogen")

plot(predict(model.OrgN), residuals(model.OrgN), ylab = "Residuals", xlab = "Predicted",
     main="Organic Nitrogen",  pch=19, cex=0.5)

######################################  GW

plot(predict(model.GW), data.GW$GW[1:87], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="GW")

plot(predict(model.GW), residuals(model.GW), ylab = "Residuals", xlab = "Predicted",
     main="GW",  pch=19, cex=0.5)

######################################  GW

plot(predict(model.C.N), data.C.N$C.N, pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="C:N")

plot(predict(model.C.N), residuals(model.C.N), ylab = "Residuals", xlab = "Predicted",
     main="C:N",  pch=19, cex=0.5)


dev.off()

######################################################################################################################
######################################        Supplementary Material 2  ##################################
#####################################################################################################################
library(dotwhisker)
library(broom)
library(dplyr)

########################################d15N

model2.d15N <- lm(d15N.x~Bank*log_dist+I(log_dist^2)*Bank+massN, data.d15N)

dwplot(list(model.d15N, model2.d15N),
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 1))+ 
         theme_test() + xlab("Coefficient Estimate") + ylab("") +
         geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
         ggtitle("d15N") +
         theme(plot.title = element_text(face="bold"),
               legend.position = c(0.007, 0.01),
               legend.justification = c(0, 0), 
               legend.background = element_rect(colour="grey80"),
               legend.title = element_blank())

########################################d13C
model2.d13C <- lm(d13C~log_dist+massN, data.d13C)

########################################d15NNH4
dwplot(list(model.d15NNH4),
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 1))+ 
  theme_test() + xlab("Coefficient Estimate") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("d15N NH4") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.007, 0.01),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank())

########################################NH4

model2.NH4 <- lm(NH4Conc~Bank+GW.x, data.NH4)
model3.NH4 <- lm(NH4Conc~Bank*log_dist.x, data.NH4)
model4.NH4 <- lm(NH4Conc~Bank*log_dist.x+GW.x, data.NH4)
model5.NH4 <- lm(NH4Conc~Bank*log_dist.x+I(log_dist.x^2)*Bank, data.NH4)
model6.NH4 <- lm(NH4Conc~log_dist.x, data.NH4)

dwplot(list(model.d15NNH4, model2.NH4, model3.NH4, model4.NH4),
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 1))+ 
  theme_test() + xlab("Coefficient Estimate") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("NH4") +
  theme(plot.title = element_text(face="bold"),
       
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank())



dwplot(list(model.C.N), show_intercept = TRUE)
model.d15N <- lm(d15N.x~Bank*log_dist+I(log_dist^2)*Bank, data.d15N)
model.d13C <- lm(d13C~log_dist, data.d13C)
model.pctN <- lm(pctN.x~Bank*log_dist, data.pctN)
model.Min <- lm(NetMin~OrgN, data.Min)
model.Nit <- lm(NetNit~Bank+NH4Conc+GW, data.Nit)
model.NH4 <- lm(NH4Conc~Bank, data.NH4)
model.NO3 <- lm(NO3Conc~Bank+GW.x, data.NO3)
model.OrgN <- lm(OrgN~Bank*log_dist+I(log_dist^2)*Bank+GW, data.OrgN)
model.d15NNO3 <- lm(d15N~Bank*log_dist+N.ug, data.d15NNO3)
model.d15NNH4 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data.d15NNH4)
model.massN <- lm(massN~Bank*log_dist+GW, data.mass)
model.C.N <- lm(C.N~Bank*log_dist+I(log_dist^2)*Bank, data.C.N)
model.GW <- lm(GW~log_dist, data.GW)


coefplot(model.Min, main = "GW")