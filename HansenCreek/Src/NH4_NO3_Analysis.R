rm(list = ls())
setwd("~/Documents/Graduate Work 2/Hansen Creek/Results")

library(readr)
Modified_HansenSoil_KClData_24Aug2017 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Modified_HansenSoil_KClData_24Aug2017.csv") 

data1 <- Modified_HansenSoil_KClData_24Aug2017
colnames(data1) <- c("Sample", "ID", "Lab", "NH4", "NO3", "Side", "Distance", "Transect", "Incubation", "Location", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin")

tin_weights <- c(3.081, 3.098, 3.118, 3.123, 3.078, 3.107, 3.091, 3.129, 3.108)
mean_tin <- mean(tin_weights)
days <- 15 #number of days incubated

BlankI <- subset(data1, Location=="Blank" & Incubation=="I")
BlankI_NH4 <- mean(BlankI$NH4)
BlankI_NO3 <- mean(BlankI$NO3)

BlankF <- subset(data1, Location=="Blank" & Incubation=="F")
BlankF_NH4 <- mean(BlankF$NH4)
BlankF_NO3 <- mean(BlankF$NO3)

NH4_Corrected <- ifelse(data1$Incubation == "I", data1$NH4-BlankI_NH4, ifelse(data1$Incubation=="F", data1$NH4-BlankF_NH4, NA))
NO3_Corrected <- ifelse(data1$Incubation == "I", data1$NO3-BlankI_NO3, ifelse(data1$Incubation=="F", data1$NO3-BlankF_NO3, NA))

Wet_Weight <- data1$Soil.Wet.Tin - mean_tin
Dry_Weight <- data1$Soil.Dry.Tin - mean_tin
Water_Weight <- Wet_Weight - Dry_Weight
GW <- Water_Weight/Dry_Weight
NH4_Conc <- NH4_Corrected*1000*(100/1000)*(data1$Initial.Wet*(1+GW))  ###concentrations corrected for gravimetric water content
NO3_Conc <- NO3_Corrected*1000*(100/1000)*(data1$Initial.Wet*(1+GW)) 
 

data2<- cbind(data1, NH4_Corrected, NO3_Corrected, Wet_Weight, Dry_Weight, Water_Weight, GW, NH4_Conc, NO3_Conc)

data2.1 <- subset(data2, !(ID == "Blank1 " | ID == "Blank2" | ID == "Blank3" | ID == "Blank1F" | ID == "Blank2F" | ID == "Blank3F"))
length(data2.1$Sample)
Initial <- subset(data2.1, Incubation ==  "I")
Final <- subset(data2.1, Incubation ==  "F")
length(Initial$Sample)
length(Final$Sample)
Net_Mineralization <- ((Final$NH4_Conc+Final$NO3_Conc)-(Initial$NH4_Conc+Initial$NO3_Conc))/days
Net_Nitrification <- (Final$NO3_Conc-Initial$NO3_Conc)/days

data3 <- cbind(Initial, Net_Mineralization, Net_Nitrification)
colnames(data3)<- c("Sample", "ID", "Lab", "NH4", "NO3", "Bank", "Distance", "Site", "Incubation","Reach", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin", "NH4_Corr", "N03_Corr", "Wet_Weight", "Dry_Weight", "Water_Weight", "GW", "NH4_Conc", "N03_Conc", "NetMin", "NetNit")
write.csv(data3, file = "Calculatd_Data.csv")

data3 <-data3[-79,]
data3 <-data3[-11,]

##################################################################################################################################
#####################################       Linear Mixed Effects Model      #######################################################
#################################################################################################################################
library(nlme)
library(lme4)
library(lmerTest)
library(AICcmodavg)

#########################################################################################################################
######################################            Net Mineralization                #####################################
#########################################################################################################################

summary(nest.Full.6m.d15N)
data1m <- data3[which(data3$Distance==1),]
nest.Full.1m.min <- lmer(NetMin~Bank+GW+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.min  <- lmer(NetMin~GW+(1|Reach/Site), data1m, REML=0)

data3m <- data3[which(data3$Distance==3),]
nest.Full.3m.min <- lmer(NetMin~Bank+GW+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.min <- lmer(NetMin~GW+(1|Reach/Site), data3m, REML=0)

data6m <- data3[which(data3$Distance==6),]
nest.Full.6m.min <- lmer(NetMin~Bank+GW+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.min <- lmer(NetMin~GW+(1|Reach/Site), data6m, REML=0)

data10m <- data3[which(data3$Distance==10),]
nest.Full.10m.min <- lmer(NetMin~Bank+GW+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.min  <- lmer(NetMin~GW+(1|Reach/Site), data10m, REML=0)

data20m <- data3[which(data3$Distance==20),]
nest.Full.20m.min <- lmer(NetMin~Bank+GW+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.min  <- lmer(NetMin~GW+(1|Reach/Site), data20m, REML=0)

AIC.min.full <- c(AIC(logLik(nest.Full.1m.min)), AIC(logLik(nest.Full.3m.min)), AIC(logLik(nest.Full.6m.min)), AIC(logLik(nest.Full.10m.min)), AIC(logLik(nest.Full.20m.min)))
AIC.min.null <- c(AIC(logLik(nest.Null.1m.min)), AIC(logLik(nest.Null.3m.min)), AIC(logLik(nest.Null.6m.min)), AIC(logLik(nest.Null.10m.min)), AIC(logLik(nest.Null.20m.min)))

AICc.min.full <- c(AICc(nest.Full.1m.min), AICc(nest.Full.3m.min), AICc(nest.Full.6m.min), AICc(nest.Full.10m.min), AICc(nest.Full.20m.min))
AICc.min.null <- c(AICc(nest.Null.1m.min), AICc(nest.Null.3m.min), AICc(nest.Null.6m.min), AICc(nest.Null.10m.min), AICc(nest.Null.20m.min))


Mineralization <- data.frame(AIC.Full=AIC.min.full, AIC.Null=AIC.min.null, AICc.Full=AICc.min.full, AICc.Null=AICc.min.null)
row.names(Mineralization)<- c("1m", "3m", "6m", "10m", "20m")

DelAIC.Full<-rep(1, 5)
for(i in 1:5){
  if(Mineralization$AIC.Full[i]>Mineralization$AIC.Null[i]){
    DelAIC.Full[i]<- Mineralization$AIC.Full[i]- Mineralization$AIC.Null[i]
  }else {DelAIC.Full[i] <-0}
  print(DelAIC.Full)
}


DelAIC.Null<-rep(1, 5)
for(i in 1:5){
  if(Mineralization$AIC.Null[i]>Mineralization$AIC.Full[i]){
    DelAIC.Null[i]<- Mineralization$AIC.Null[i]- Mineralization$AIC.Full[i]
  }else {DelAIC.Null[i] <-0}
  print(DelAIC.Null)
}

DelAICc.Full<-rep(1, 5)
for(i in 1:5){
  if(Mineralization$AICc.Full[i]>Mineralization$AICc.Null[i]){
    DelAICc.Full[i]<- Mineralization$AICc.Full[i]- Mineralization$AICc.Null[i]
  }else {DelAICc.Full[i] <-0}
  print(DelAICc.Full)
}


DelAICc.Null<-rep(1, 5)
for(i in 1:5){
  if(Mineralization$AICc.Null[i]>Mineralization$AICc.Full[i]){
    DelAICc.Null[i]<- Mineralization$AICc.Null[i]- Mineralization$AICc.Full[i]
  }else {DelAICc.Null[i] <-0}
  print(DelAICc.Null)
}


library(lmerTest)

EstimateBankL <- rbind(coef(summary(nest.Full.1m.min))[1,1], coef(summary(nest.Full.3m.min))[1,1], coef(summary(nest.Full.6m.min))[1,1], coef(summary(nest.Full.10m.min))[1,1], coef(summary(nest.Full.20m.min))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.min))[1,2], coef(summary(nest.Full.3m.min))[1,2], coef(summary(nest.Full.6m.min))[1,2], coef(summary(nest.Full.10m.min))[1,2], coef(summary(nest.Full.20m.min))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.min))[1,5], coef(summary(nest.Full.3m.min))[1,5], coef(summary(nest.Full.6m.min))[1,5], coef(summary(nest.Full.10m.min))[1,5], coef(summary(nest.Full.20m.min))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.min))[2,1], coef(summary(nest.Full.3m.min))[2,1], coef(summary(nest.Full.6m.min))[2,1], coef(summary(nest.Full.10m.min))[2,1], coef(summary(nest.Full.20m.min))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.min))[2,2], coef(summary(nest.Full.3m.min))[2,2], coef(summary(nest.Full.6m.min))[2,2], coef(summary(nest.Full.10m.min))[2,2], coef(summary(nest.Full.20m.min))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.min))[2,5], coef(summary(nest.Full.3m.min))[2,5], coef(summary(nest.Full.6m.min))[2,5], coef(summary(nest.Full.10m.min))[2,5], coef(summary(nest.Full.20m.min))[2,5])

Mineralization <- cbind(Mineralization, DelAIC.Full, DelAIC.Null, DelAICc.Full, DelAICc.Null)
Mineralization <- cbind(Mineralization, cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
names <- c("AIC.Full", "AIC.Null", "AICc.Full", "AICc.Null", "DelAIC.Full", "DelAIC.Null", "DelAICc.Full", "DelAICc.Null", "IntBankL", "StdBankL", "pBankL", "EstimateBankR", "StdBankR", "pBankR")
colnames(Mineralization) <- names

####testing the residuals
qqnorm(resid(nest.Full.1m.min), pch=16, main="Min")
qqline(resid(nest.Full.1m.min))

qqnorm(resid(nest.Full.3m.min), pch=16, main="Min")
qqline(resid(nest.Full.3m.min))

qqnorm(resid(nest.Full.6m.min), pch=16, main="Min")
qqline(resid(nest.Full.6m.min))

qqnorm(resid(nest.Full.10m.min), pch=16, main="Min")
qqline(resid(nest.Full.10m.min))

qqnorm(resid(nest.Full.20m.min), pch=16, main="Min")
qqline(resid(nest.Full.20m.min))

Min.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.min))$p.value, shapiro.test(resid(nest.Full.3m.min))$p.value,
                   shapiro.test(resid(nest.Full.6m.min))$p.value, shapiro.test(resid(nest.Full.10m.min))$p.value,
                   shapiro.test(resid(nest.Full.20m.min))$p.value)
                   
qqnorm(resid(nest.Null.1m.min), pch=16, main="Min")
qqline(resid(nest.Null.1m.min))

qqnorm(resid(nest.Null.3m.min), pch=16, main="Min")
qqline(resid(nest.Null.3m.min))

qqnorm(resid(nest.Null.6m.min), pch=16, main="Min")
qqline(resid(nest.Null.6m.min))

qqnorm(resid(nest.Null.10m.min), pch=16, main="Min")
qqline(resid(nest.Null.10m.min))

qqnorm(resid(nest.Null.20m.min), pch=16, main="Min")
qqline(resid(nest.Null.20m.min))

Min.Shap.Null <- c(shapiro.test(resid(nest.Null.1m.min))$p.value, shapiro.test(resid(nest.Null.3m.min))$p.value,
                   shapiro.test(resid(nest.Null.6m.min))$p.value, shapiro.test(resid(nest.Null.10m.min))$p.value,
                   shapiro.test(resid(nest.Null.20m.min))$p.value)




#########################################################################################################################
######################################            Net Nitrification                  #####################################
#########################################################################################################################

nest.Full.1m.nit <- lmer(NetNit~Bank+GW+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.nit  <- lmer(NetNit~GW+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.nit <- lmer(NetNit~Bank+GW+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.nit <- lmer(NetNit~GW+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.nit <- lmer(NetNit~Bank+GW+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.nit <- lmer(NetNit~GW+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.nit <- lmer(NetNit~Bank+GW+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.nit  <- lmer(NetNit~GW+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.nit <- lmer(NetNit~Bank+GW+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.nit  <- lmer(NetNit~GW+(1|Reach/Site), data20m, REML=0)

AIC.nit.full <- c(AIC(logLik(nest.Full.1m.nit)), AIC(logLik(nest.Full.3m.nit)), AIC(logLik(nest.Full.6m.nit)), AIC(logLik(nest.Full.10m.nit)), AIC(logLik(nest.Full.20m.nit)))
AIC.nit.null <- c(AIC(logLik(nest.Null.1m.nit)), AIC(logLik(nest.Null.3m.nit)), AIC(logLik(nest.Null.6m.nit)), AIC(logLik(nest.Null.10m.nit)), AIC(logLik(nest.Null.20m.nit)))

AICc.nit.full <- c(AICc(nest.Full.1m.nit), AICc(nest.Full.3m.nit), AICc(nest.Full.6m.nit), AICc(nest.Full.10m.nit), AICc(nest.Full.20m.nit))
AICc.nit.null <- c(AICc(nest.Null.1m.nit), AICc(nest.Null.3m.nit), AICc(nest.Null.6m.nit), AICc(nest.Null.10m.nit), AICc(nest.Null.20m.nit))


Nitrification <- data.frame(AIC.Full=AIC.nit.full, AIC.Null=AIC.nit.null, AICc.Full=AICc.nit.full, AICc.Null=AICc.nit.null)
row.names(Nitrification)<- c("1m", "3m", "6m", "10m", "20m")

DelAIC.Full<-rep(1, 5)
for(i in 1:5){
  if(Nitrification$AIC.Full[i]>Nitrification$AIC.Null[i]){
    DelAIC.Full[i]<- Nitrification$AIC.Full[i]- Nitrification$AIC.Null[i]
  }else {DelAIC.Full[i] <-0}
  print(DelAIC.Full)
}


DelAIC.Null<-rep(1, 5)
for(i in 1:5){
  if(Nitrification$AIC.Null[i]>Nitrification$AIC.Full[i]){
    DelAIC.Null[i]<- Nitrification$AIC.Null[i]- Nitrification$AIC.Full[i]
  }else {DelAIC.Null[i] <-0}
  print(DelAIC.Null)
}

DelAICc.Full<-rep(1, 5)
for(i in 1:5){
  if(Nitrification$AICc.Full[i]>Nitrification$AICc.Null[i]){
    DelAICc.Full[i]<- Nitrification$AICc.Full[i]- Nitrification$AICc.Null[i]
  }else {DelAICc.Full[i] <-0}
  print(DelAICc.Full)
}


DelAICc.Null<-rep(1, 5)
for(i in 1:5){
  if(Nitrification$AICc.Null[i]>Nitrification$AICc.Full[i]){
    DelAICc.Null[i]<- Nitrification$AICc.Null[i]- Nitrification$AICc.Full[i]
  }else {DelAICc.Null[i] <-0}
  print(DelAICc.Null)
}

EstimateBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,1], coef(summary(nest.Full.3m.nit))[1,1], coef(summary(nest.Full.6m.nit))[1,1], coef(summary(nest.Full.10m.nit))[1,1], coef(summary(nest.Full.20m.nit))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,2], coef(summary(nest.Full.3m.nit))[1,2], coef(summary(nest.Full.6m.nit))[1,2], coef(summary(nest.Full.10m.nit))[1,2], coef(summary(nest.Full.20m.nit))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,5], coef(summary(nest.Full.3m.nit))[1,5], coef(summary(nest.Full.6m.nit))[1,5], coef(summary(nest.Full.10m.nit))[1,5], coef(summary(nest.Full.20m.nit))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,1], coef(summary(nest.Full.3m.nit))[2,1], coef(summary(nest.Full.6m.nit))[2,1], coef(summary(nest.Full.10m.nit))[2,1], coef(summary(nest.Full.20m.nit))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,2], coef(summary(nest.Full.3m.nit))[2,2], coef(summary(nest.Full.6m.nit))[2,2], coef(summary(nest.Full.10m.nit))[2,2], coef(summary(nest.Full.20m.nit))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,5], coef(summary(nest.Full.3m.nit))[2,5], coef(summary(nest.Full.6m.nit))[2,5], coef(summary(nest.Full.10m.nit))[2,5], coef(summary(nest.Full.20m.nit))[2,5])

Nitrification <- cbind(Nitrification, DelAIC.Full, DelAIC.Null, DelAICc.Full, DelAICc.Null)
Nitrification <- cbind(Nitrification, cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
colnames(Nitrification) <- names

####testing the residuals
qqnorm(resid(nest.Full.1m.nit), pch=16, main="Nit")
qqline(resid(nest.Full.1m.nit))

qqnorm(resid(nest.Full.3m.nit), pch=16, main="Nit")
qqline(resid(nest.Full.3m.nit))

qqnorm(resid(nest.Full.6m.nit), pch=16, main="Nit")
qqline(resid(nest.Full.6m.nit))

qqnorm(resid(nest.Full.10m.nit), pch=16, main="Nit")
qqline(resid(nest.Full.10m.nit))

qqnorm(resid(nest.Full.20m.nit), pch=16, main="Nit")
qqline(resid(nest.Full.20m.nit))

nit.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.nit))$p.value, shapiro.test(resid(nest.Full.3m.nit))$p.value,
                   shapiro.test(resid(nest.Full.6m.nit))$p.value, shapiro.test(resid(nest.Full.10m.nit))$p.value,
                   shapiro.test(resid(nest.Full.20m.nit))$p.value)

qqnorm(resid(nest.Null.1m.nit), pch=16, main="Nit")
qqline(resid(nest.Null.1m.nit))

qqnorm(resid(nest.Null.3m.nit), pch=16, main="Nit")
qqline(resid(nest.Null.3m.nit))

qqnorm(resid(nest.Null.6m.nit), pch=16, main="Nit")
qqline(resid(nest.Null.6m.nit))

qqnorm(resid(nest.Null.10m.nit), pch=16, main="Nit")
qqline(resid(nest.Null.10m.nit))

qqnorm(resid(nest.Null.20m.nit), pch=16, main="Nit")
qqline(resid(nest.Null.20m.nit))

nit.Shap.Null <- c(shapiro.test(resid(nest.Null.1m.nit))$p.value, shapiro.test(resid(nest.Null.3m.nit))$p.value,
                   shapiro.test(resid(nest.Null.6m.nit))$p.value, shapiro.test(resid(nest.Null.10m.nit))$p.value,
                   shapiro.test(resid(nest.Null.20m.nit))$p.value)



#########################################################################################################################
######################################            NH4 Concentration                #####################################
#########################################################################################################################

nest.Full.1m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.NH4  <- lmer(NH4_Conc~GW+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.NH4 <- lmer(NH4_Conc~GW+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.NH4 <- lmer(NH4_Conc~GW+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.NH4  <- lmer(NH4_Conc~GW+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.NH4  <- lmer(NH4_Conc~GW+(1|Reach/Site), data20m, REML=0)

AIC.NH4.full <- c(AIC(logLik(nest.Full.1m.NH4)), AIC(logLik(nest.Full.3m.NH4)), AIC(logLik(nest.Full.6m.NH4)), AIC(logLik(nest.Full.10m.NH4)), AIC(logLik(nest.Full.20m.NH4)))
AIC.NH4.null <- c(AIC(logLik(nest.Null.1m.NH4)), AIC(logLik(nest.Null.3m.NH4)), AIC(logLik(nest.Null.6m.NH4)), AIC(logLik(nest.Null.10m.NH4)), AIC(logLik(nest.Null.20m.NH4)))

NH4Concentration <- data.frame(AIC.Full=AIC.NH4.full, AIC.Null=AIC.NH4.null)
row.names(NH4Concentration)<- c("1m", "3m", "6m", "10m", "20m")

DelAIC.Full<-rep(1, 5)
for(i in 1:5){
  if(NH4Concentration$AIC.Full[i]>NH4Concentration$AIC.Null[i]){
    DelAIC.Full[i]<- NH4Concentration$AIC.Full[i]- NH4Concentration$AIC.Null[i]
  }else {DelAIC.Full[i] <-0}
  print(DelAIC.Full)
}


DelAIC.Null<-rep(1, 5)
for(i in 1:5){
  if(NH4Concentration$AIC.Null[i]>NH4Concentration$AIC.Full[i]){
    DelAIC.Null[i]<- NH4Concentration$AIC.Null[i]- NH4Concentration$AIC.Full[i]
  }else {DelAIC.Null[i] <-0}
  print(DelAIC.Null)
}


EstimateBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,1], coef(summary(nest.Full.3m.NH4))[1,1], coef(summary(nest.Full.6m.NH4))[1,1], coef(summary(nest.Full.10m.NH4))[1,1], coef(summary(nest.Full.20m.NH4))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,2], coef(summary(nest.Full.3m.NH4))[1,2], coef(summary(nest.Full.6m.NH4))[1,2], coef(summary(nest.Full.10m.NH4))[1,2], coef(summary(nest.Full.20m.NH4))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,5], coef(summary(nest.Full.3m.NH4))[1,5], coef(summary(nest.Full.6m.NH4))[1,5], coef(summary(nest.Full.10m.NH4))[1,5], coef(summary(nest.Full.20m.NH4))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,1], coef(summary(nest.Full.3m.NH4))[2,1], coef(summary(nest.Full.6m.NH4))[2,1], coef(summary(nest.Full.10m.NH4))[2,1], coef(summary(nest.Full.20m.NH4))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,2], coef(summary(nest.Full.3m.NH4))[2,2], coef(summary(nest.Full.6m.NH4))[2,2], coef(summary(nest.Full.10m.NH4))[2,2], coef(summary(nest.Full.20m.NH4))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,5], coef(summary(nest.Full.3m.NH4))[2,5], coef(summary(nest.Full.6m.NH4))[2,5], coef(summary(nest.Full.10m.NH4))[2,5], coef(summary(nest.Full.20m.NH4))[2,5])

NH4Concentration <- cbind(NH4Concentration, DelAIC.Full, DelAIC.Null)
NH4Concentration  <- cbind(NH4Concentration , cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
colnames(NH4Concentration ) <- names


####testing the residuals
qqnorm(resid(nest.Full.1m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Full.1m.NH4))

qqnorm(resid(nest.Full.3m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Full.3m.NH4))

qqnorm(resid(nest.Full.6m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Full.6m.NH4))

qqnorm(resid(nest.Full.10m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Full.10m.NH4))

qqnorm(resid(nest.Full.20m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Full.20m.NH4))

NH4.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.NH4))$p.value, shapiro.test(resid(nest.Full.3m.NH4))$p.value,
                   shapiro.test(resid(nest.Full.6m.NH4))$p.value, shapiro.test(resid(nest.Full.10m.NH4))$p.value,
                   shapiro.test(resid(nest.Full.20m.NH4))$p.value)

qqnorm(resid(nest.Null.1m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Null.1m.NH4))

qqnorm(resid(nest.Null.3m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Null.3m.NH4))

qqnorm(resid(nest.Null.6m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Null.6m.NH4))

qqnorm(resid(nest.Null.10m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Null.10m.NH4))

qqnorm(resid(nest.Null.20m.NH4), pch=16, main="Nh4")
qqline(resid(nest.Null.20m.NH4))

NH4.Shap.Null <- c(shapiro.test(resid(nest.Null.1m.NH4))$p.value, shapiro.test(resid(nest.Null.3m.NH4))$p.value,
                   shapiro.test(resid(nest.Null.6m.NH4))$p.value, shapiro.test(resid(nest.Null.10m.NH4))$p.value,
                   shapiro.test(resid(nest.Null.20m.NH4))$p.value)




#########################################################################################################################
######################################            N03 Concentration                 #####################################
#########################################################################################################################

nest.Full.1m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.N03  <- lmer(N03_Conc~GW+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.N03 <- lmer(N03_Conc~GW+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.N03 <- lmer(N03_Conc~GW+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.N03  <- lmer(N03_Conc~GW+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.N03  <- lmer(N03_Conc~GW+(1|Reach/Site), data20m, REML=0)

AIC.N03.full <- c(AIC(logLik(nest.Full.1m.N03)), AIC(logLik(nest.Full.3m.N03)), AIC(logLik(nest.Full.6m.N03)), AIC(logLik(nest.Full.10m.N03)), AIC(logLik(nest.Full.20m.N03)))
AIC.N03.null <- c(AIC(logLik(nest.Null.1m.N03)), AIC(logLik(nest.Null.3m.N03)), AIC(logLik(nest.Null.6m.N03)), AIC(logLik(nest.Null.10m.N03)), AIC(logLik(nest.Null.20m.N03)))

N03Concentration <- data.frame(AIC.Full=AIC.N03.full, AIC.Null=AIC.N03.null)
row.names(N03Concentration)<- c("1m", "3m", "6m", "10m", "20m")

DelAIC.Full<-rep(1, 5)
for(i in 1:5){
  if(N03Concentration$AIC.Full[i]>N03Concentration$AIC.Null[i]){
    DelAIC.Full[i]<- N03Concentration$AIC.Full[i]- N03Concentration$AIC.Null[i]
  }else {DelAIC.Full[i] <-0}
  print(DelAIC.Full)
}


DelAIC.Null<-rep(1, 5)
for(i in 1:5){
  if(N03Concentration$AIC.Null[i]>N03Concentration$AIC.Full[i]){
    DelAIC.Null[i]<- N03Concentration$AIC.Null[i]- N03Concentration$AIC.Full[i]
  }else {DelAIC.Null[i] <-0}
  print(DelAIC.Null)
}



EstimateBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,1], coef(summary(nest.Full.3m.N03))[1,1], coef(summary(nest.Full.6m.N03))[1,1], coef(summary(nest.Full.10m.N03))[1,1], coef(summary(nest.Full.20m.N03))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,2], coef(summary(nest.Full.3m.N03))[1,2], coef(summary(nest.Full.6m.N03))[1,2], coef(summary(nest.Full.10m.N03))[1,2], coef(summary(nest.Full.20m.N03))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,5], coef(summary(nest.Full.3m.N03))[1,5], coef(summary(nest.Full.6m.N03))[1,5], coef(summary(nest.Full.10m.N03))[1,5], coef(summary(nest.Full.20m.N03))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,1], coef(summary(nest.Full.3m.N03))[2,1], coef(summary(nest.Full.6m.N03))[2,1], coef(summary(nest.Full.10m.N03))[2,1], coef(summary(nest.Full.20m.N03))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,2], coef(summary(nest.Full.3m.N03))[2,2], coef(summary(nest.Full.6m.N03))[2,2], coef(summary(nest.Full.10m.N03))[2,2], coef(summary(nest.Full.20m.N03))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,5], coef(summary(nest.Full.3m.N03))[2,5], coef(summary(nest.Full.6m.N03))[2,5], coef(summary(nest.Full.10m.N03))[2,5], coef(summary(nest.Full.20m.N03))[2,5])

N03Concentration <- cbind(N03Concentration, DelAIC.Full, DelAIC.Null)
N03Concentration  <- cbind(N03Concentration , cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
colnames(N03Concentration ) <- names




####testing the residuals
qqnorm(resid(nest.Full.1m.N03), pch=16, main="N03")
qqline(resid(nest.Full.1m.N03))

qqnorm(resid(nest.Full.3m.N03), pch=16, main="N03")
qqline(resid(nest.Full.3m.N03))

qqnorm(resid(nest.Full.6m.N03), pch=16, main="N03")
qqline(resid(nest.Full.6m.N03))

qqnorm(resid(nest.Full.10m.N03), pch=16, main="N03")
qqline(resid(nest.Full.10m.N03))

qqnorm(resid(nest.Full.20m.N03), pch=16, main="N03")
qqline(resid(nest.Full.20m.N03))

N03.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.N03))$p.value, shapiro.test(resid(nest.Full.3m.N03))$p.value,
                   shapiro.test(resid(nest.Full.6m.N03))$p.value, shapiro.test(resid(nest.Full.10m.N03))$p.value,
                   shapiro.test(resid(nest.Full.20m.N03))$p.value)

qqnorm(resid(nest.Null.1m.N03), pch=16, main="N03")
qqline(resid(nest.Null.1m.N03))

qqnorm(resid(nest.Null.3m.N03), pch=16, main="N03")
qqline(resid(nest.Null.3m.N03))

qqnorm(resid(nest.Null.6m.N03), pch=16, main="N03")
qqline(resid(nest.Null.6m.N03))

qqnorm(resid(nest.Null.10m.N03), pch=16, main="N03")
qqline(resid(nest.Null.10m.N03))

qqnorm(resid(nest.Null.20m.N03), pch=16, main="N03")
qqline(resid(nest.Null.20m.N03))

N03.Shap.Null <- c(shapiro.test(resid(nest.Null.1m.N03))$p.value, shapiro.test(resid(nest.Null.3m.N03))$p.value,
                   shapiro.test(resid(nest.Null.6m.N03))$p.value, shapiro.test(resid(nest.Null.10m.N03))$p.value,
                   shapiro.test(resid(nest.Null.20m.N03))$p.value)






write.csv(Mineralization, file = "Mineralization.csv")
write.csv(Nitrification, file = "Nitrification.csv")
write.csv(N03Concentration, file = "N03Concentration.csv")
write.csv(NH4Concentration, file = "NH4Concentration.csv")

############################################################################################################
############################          Bulk Isotope Models             ######################################
############################################################################################################

data4 <- Hansen_bulkSoilIsotopes <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Hansen_bulkSoilIsotopes.csv")
colnames(data4) <- c("ID", "Reach", "Site", "Bank", "Distance", "Replicate", "Weight", "pctN", "d15N", "pctC", "d13C")
data4.1 <- subset(data4, !(Replicate == 2))
length(data4.1)

###############################################################################
################### Calculating Org N #########################################
###############################################################################

massN <- data4$Weight * data4$pctN /100
data4.2 <- cbind(data4, massN)

data4.3 <- merge(data4.2, data3, by = c("Reach", "Site", "Distance", "Bank"))
massNH4 <- data4.3$NH4_Conc/1000000*data4.3$Weight
massNO3 <- data4.3$N03_Conc/1000000*data4.3$Weight
OrgN <- (data4.3$massN-massNH4-massNO3)*1000

data4.4 <- cbind(data4.3, massNH4, massNO3, OrgN)

data5 <- data.frame(data4.4$ID.x, data4.4$Reach, data4.4$Site, data4.4$Distance, data4.4$Bank, 
          data4.4$massN, data4.4$massNH4, data4.4$massNO3, data4.4$OrgN)
colnames(data5) <- c("ID", "Reach", "Site", "Distance", "Bank", "massN", "massNH4", "massNO3", "OrgN")

write.csv(data5, file = "OrganicNitrogen.csv")

data5.1 <- merge(data4.1, data5, by = c("Reach", "Site", "Distance", "Bank", "ID"))

###########################################         d15N          ##########################################
data1m <- data4.1[which(data5.1$Distance==1),]
nest.Full.1m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.d15N  <- lmer(d15N~(1|Reach/Site), data1m, REML=0)

data3m <- data4.1[which(data5.1$Distance==3),]
nest.Full.3m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.d15N <- lmer(d15N~(1|Reach/Site), data3m, REML=0)
anova(nest.Full.3m.d15N, nest.Null.3m.d15N)
data3mRR <- data3m[which(data3m$Bank=='Right'),]
data3mRR <- data3m[which(data3m$Bank=='Right'),]

data6m <- data4.1[which(data5.1$Distance==6),]
nest.Full.6m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.d15N <- lmer(d15N~(1|Reach/Site), data6m, REML=0)

data10m <- data4.1[which(data5.1$Distance==10),]
nest.Full.10m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.d15N  <- lmer(d15N~(1|Reach/Site), data10m, REML=0)

data20m <- data4.1[which(data5.1$Distance==20),]
nest.Full.20m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.d15N  <- lmer(d15N~(1|Reach/Site), data20m, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,1], coef(summary(nest.Full.3m.d15N))[1,1], coef(summary(nest.Full.6m.d15N))[1,1], coef(summary(nest.Full.10m.d15N))[1,1], coef(summary(nest.Full.20m.d15N))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,2], coef(summary(nest.Full.3m.d15N))[1,2], coef(summary(nest.Full.6m.d15N))[1,2], coef(summary(nest.Full.10m.d15N))[1,2], coef(summary(nest.Full.20m.d15N))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,5], coef(summary(nest.Full.3m.d15N))[1,5], coef(summary(nest.Full.6m.d15N))[1,5], coef(summary(nest.Full.10m.d15N))[1,5], coef(summary(nest.Full.20m.d15N))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,1], coef(summary(nest.Full.3m.d15N))[2,1], coef(summary(nest.Full.6m.d15N))[2,1], coef(summary(nest.Full.10m.d15N))[2,1], coef(summary(nest.Full.20m.d15N))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,2], coef(summary(nest.Full.3m.d15N))[2,2], coef(summary(nest.Full.6m.d15N))[2,2], coef(summary(nest.Full.10m.d15N))[2,2], coef(summary(nest.Full.20m.d15N))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,5], coef(summary(nest.Full.3m.d15N))[2,5], coef(summary(nest.Full.6m.d15N))[2,5], coef(summary(nest.Full.10m.d15N))[2,5], coef(summary(nest.Full.20m.d15N))[2,5])

d15N <- cbind(d15N, DelAIC.Full, DelAIC.Null, DelAICc.Full, DelAICc.Null)
d15N  <- cbind(d15N , cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
colnames(d15N ) <- names
write.csv(d15N, file = "d15N.csv")



###########################################         Organic N          ##########################################
data1m <- data4.1[which(data5.1$Distance==1),]
nest.Full.1m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data1m, REML=0)
nest.Null.1m.d15N  <- lmer(d15N~(1|Reach/Site), data1m, REML=0)

data3m <- data4.1[which(data5.1$Distance==3),]
nest.Full.3m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data3m, REML=0)
nest.Null.3m.d15N <- lmer(d15N~(1|Reach/Site), data3m, REML=0)
anova(nest.Full.3m.d15N, nest.Null.3m.d15N)
data3mRR <- data3m[which(data3m$Bank=='Right'),]
data3mRR <- data3m[which(data3m$Bank=='Right'),]

data6m <- data4.1[which(data5.1$Distance==6),]
nest.Full.6m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data6m, REML=0)
nest.Null.6m.d15N <- lmer(d15N~(1|Reach/Site), data6m, REML=0)

data10m <- data4.1[which(data5.1$Distance==10),]
nest.Full.10m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data10m, REML=0)
nest.Null.10m.d15N  <- lmer(d15N~(1|Reach/Site), data10m, REML=0)

data20m <- data4.1[which(data5.1$Distance==20),]
nest.Full.20m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data20m, REML=0)
nest.Null.20m.d15N  <- lmer(d15N~(1|Reach/Site), data20m, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,1], coef(summary(nest.Full.3m.d15N))[1,1], coef(summary(nest.Full.6m.d15N))[1,1], coef(summary(nest.Full.10m.d15N))[1,1], coef(summary(nest.Full.20m.d15N))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,2], coef(summary(nest.Full.3m.d15N))[1,2], coef(summary(nest.Full.6m.d15N))[1,2], coef(summary(nest.Full.10m.d15N))[1,2], coef(summary(nest.Full.20m.d15N))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,5], coef(summary(nest.Full.3m.d15N))[1,5], coef(summary(nest.Full.6m.d15N))[1,5], coef(summary(nest.Full.10m.d15N))[1,5], coef(summary(nest.Full.20m.d15N))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,1], coef(summary(nest.Full.3m.d15N))[2,1], coef(summary(nest.Full.6m.d15N))[2,1], coef(summary(nest.Full.10m.d15N))[2,1], coef(summary(nest.Full.20m.d15N))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,2], coef(summary(nest.Full.3m.d15N))[2,2], coef(summary(nest.Full.6m.d15N))[2,2], coef(summary(nest.Full.10m.d15N))[2,2], coef(summary(nest.Full.20m.d15N))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,5], coef(summary(nest.Full.3m.d15N))[2,5], coef(summary(nest.Full.6m.d15N))[2,5], coef(summary(nest.Full.10m.d15N))[2,5], coef(summary(nest.Full.20m.d15N))[2,5])

d15N <- cbind(d15N, DelAIC.Full, DelAIC.Null, DelAICc.Full, DelAICc.Null)
d15N  <- cbind(d15N , cbind(EstimateBankL, StdBankL, pBankL), cbind(EstimateBankR, StdBankR, pBankR))
colnames(d15N ) <- names
write.csv(d15N, file = "d15N.csv")












####testing the residuals
qqnorm(resid(nest.Full.1m.d15N), pch=16, main="d15N")
qqline(resid(nest.Full.1m.d15N))

qqnorm(resid(nest.Full.3m.d15N), pch=16, main="d15N")
qqline(resid(nest.Full.3m.d15N))

qqnorm(resid(nest.Full.6m.d15N), pch=16, main="d15N")
qqline(resid(nest.Full.6m.d15N))

qqnorm(resid(nest.Full.10m.d15N), pch=16, main="d15N")
qqline(resid(nest.Full.10m.d15N))

qqnorm(resid(nest.Full.20m.d15N), pch=16, main="d15N")
qqline(resid(nest.Full.20m.d15N))

d15N.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.d15N))$p.value, shapiro.test(resid(nest.Full.3m.d15N))$p.value,
                   shapiro.test(resid(nest.Full.6m.d15N))$p.value, shapiro.test(resid(nest.Full.10m.d15N))$p.value,
                   shapiro.test(resid(nest.Full.20m.d15N))$p.value)

qqnorm(resid(nest.Null.1m.d15N), pch=16, main="d15N")
qqline(resid(nest.Null.1m.d15N))

qqnorm(resid(nest.Null.3m.d15N), pch=16, main="d15N")
qqline(resid(nest.Null.3m.d15N))

qqnorm(resid(nest.Null.6m.d15N), pch=16, main="d15N")
qqline(resid(nest.Null.6m.d15N))

qqnorm(resid(nest.Null.10m.d15N), pch=16, main="d15N")
qqline(resid(nest.Null.10m.d15N))

qqnorm(resid(nest.Null.20m.d15N), pch=16, main="d15N")
qqline(resid(nest.Null.20m.d15N))

d15N.Shap.Null <- c(shapiro.test(resid(nest.Null.1m.d15N))$p.value, shapiro.test(resid(nest.Null.3m.d15N))$p.value,
                   shapiro.test(resid(nest.Null.6m.d15N))$p.value, shapiro.test(resid(nest.Null.10m.d15N))$p.value,
                   shapiro.test(resid(nest.Null.20m.d15N))$p.value)




Mineralization
Nitrification
d15N
N03Concentration
NH4Concentration




######################################################################################################################
######################################        Shapiro testing Normality          #####################################
#####################################################################################################################
Shapiro.Residuals <- cbind(d15N.Shap.Full, d15N.Shap.Null, nit.Shap.Full, nit.Shap.Null, Min.Shap.Full, Min.Shap.Null,
                           NH4.Shap.Full, NH4.Shap.Null, N03.Shap.Full, N03.Shap.Null)
rownames(Shapiro.Residuals) <- c("1m", "3m", "6m", "10m", "20m")

write.csv(Shapiro.Residuals, file = "ShapiroResiduals.csv")

qqnorm(data4.1$d15N, pch=16, main="d15N")
qqline(data4.1$d15N)
shapiro.test(data4.1$d15N)$p.value


qqnorm(data3$NetMin, pch=16, main="Net Mineralization")
qqline(data3$NetMin)
shapiro.test(data3$NetMin)$p.value

qqnorm(data3$NetNit, pch=16, main="Net Nitrification")
qqline(data3$NetNit)
shapiro.test(data3$NetNit)$p.value

qqnorm(data3$NH4_Conc, pch=16, main="NH4")
qqline(data3$NH4_Conc)
shapiro.test(data3$NH4_Conc)$p.value

qqnorm(data3$N03_Conc, pch=16, main="N03")
qqline(data3$N03_Conc)
shapiro.test(data3$N03_Conc)$p.value


######################################################################################################################
######################################        Salmon contribution to soil        #####################################
#####################################################################################################################

######################################################################################################################
######################################        Bulk Data          #####################################
#####################################################################################################################

data3mRR <- data3m[which(data3m$Bank=='Right'),]
data3mRL <- data3m[which(data3m$Bank=='Left'),]
mean(data3mRR$d15N)  #7.800356
mean(data3mRL$d15N) #9.158011
