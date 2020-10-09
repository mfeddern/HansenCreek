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
colnames(data3)<- c("Sample", "ID", "Lab", "NH4", "NO3", "Bank", "Distance", "Site", "Incubation","Reach", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin", "NH4_Corr", "NO3_Corr", "Wet_Weight", "Dry_Weight", "Water_Weight", "GW", "NH4_Conc", "NO3_Conc", "NetMin", "NetNit")
write.csv(data3, file = "Calculatd_Data.csv")

data3 <-data3[-79,]
data3 <-data3[-11,]

###############################################################################
################### Calculating Org N #########################################
###############################################################################

############################################################################################################
############################          Bulk Isotope Models             ######################################
############################################################################################################

names <- c("Analysis","IntBankL", "StdBankL", "pBankL", "EstimateBankR", "StdBankR", "pBankR")
names2 <- c("1m", "3m", "6m", "10m", "20m")

data4 <- Hansen_bulkSoilIsotopes <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Hansen_bulkSoilIsotopes.csv")
colnames(data4) <- c("ID", "Reach", "Site", "Bank", "Distance", "Replicate", "Weight", "pctN", "d15N", "pctC", "d13C")
data4.1 <- subset(data4, !(Replicate == 2))
length(data4.1)
data4.11 <- subset(data4.1, !(pctN > 100))

massN <- (data4.11$Weight * data4.11$pctN) #/data4.11$Weight
massN2 <- data4.11$pctN*10

data4.2 <- cbind(data4.11, massN, massN2)
data4.3 <- Initial
colnames(data4.3) <- c("Sample", "ID", "Lab", "NH4", "NO3", "Bank", "Distance", "Site", "Incubation", "Reach", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin", "NH4Corrected", "NO3Corrected", "Wet_Weight", "Dry_Weight", "Water_Weight", "GW", "NH4_Conc", "NO3_Conc")
data4.4 <- merge(data4.3, data4.2, by = c("Reach", "Site", "Distance", "Bank"))



mass.norm <- (data4.4$pctN)*(data4.4$Dry_Weight)
data4.5 <- cbind(mass.norm, data4.4)





write.csv(data4.5, file="MassNitrogen.csv")

data4.3 <- merge(data4.2, data3, by = c("Reach", "Site", "Distance", "Bank"))
massNH4 <- data4.3$NH4_Conc/1000000*data4.3$Weight
massNO3 <- data4.3$NO3_Conc/1000000*data4.3$Weight
OrgN <- (data4.3$massN-massNH4-massNO3)*1000

#Gravimetri
data4.4 <- cbind(data4.3, massNH4, massNO3, OrgN)
data3 <-data3[-74,]
data5 <- data.frame(data4.4$ID.x, data4.4$Reach, data4.4$Site, data4.4$Distance, data4.4$Bank, 
                    data4.4$massN, data4.4$massNH4, data4.4$massNO3, data4.4$OrgN)
colnames(data5) <- c("ID", "Reach", "Site", "Distance", "Bank", "massN", "massNH4", "massNO3", "OrgN")

write.csv(data5, file = "OrganicNitrogen.csv")

data5.1 <- merge(data4.1, data5, by = c("Reach", "Site", "Distance", "Bank", "ID"))
merge(data5.1, data3, by = c("Reach", "Site", "Distance", "Bank", "ID"))



##################################################################################################################################
#####################################       Linear Mixed Effects Model      #######################################################
#################################################################################################################################
library(nlme)
library(lme4)
library(lmerTest)
library(AICcmodavg)


#OrgN <- cbind(OrgN, DelAIC.Full, DelAICc.Full)
Analysis <- rep("Organic Nitrogen", 5)
OrgN  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(OrgN) <- names
row.names(OrgN) <- names2
write.csv(OrgN, file = "OrgN.csv")

log_dist <- log10(data5.1$Distance)
data5.1 <- cbind(data5.1, log_dist)

data3 <-merge(data3, data5.1, by = c("Reach", "Site", "Distance", "Bank"))
log_dist <- log10(data3$Distance)
data3 <- cbind(data3, log_dist)

############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################ ############################################ ############################################ 

data6 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Inorganic N Soil Extraction Data.csv") 
data6 <- data.frame(data6$Inorganic, data6$Test, data6$SampleID, data6$Reach, data6$ReachID, data6$Bank, data6$Distance, data6$PeakArea.N..Vs., data6$N.quantity..ug., data6$d15N.vs.Air.N2..permil.)
colnames(data6) <- c("Batch", "Test", "ID", "Reach", "Site", "Bank", "Distance", "Peak.Area", "N ug", "d15N")
#data6 <- cbind(data6, log_dist)
data6.1 <- subset(data6, Bank=="Right" | Bank=="Left")
data6.NH4 <- subset(data6.1, Test=="NH4")
data6.NO3 <- subset(data6.1, Test=="NO3")
write.csv(data6.NO3, file = "NO3.StableIsotope.csv")
write.csv(data6.NH4, file = "NH4.StableIsotope.csv")

data6.NH4 <- merge(data6.NH4, data3, by = c("Reach", "Site", "Distance", "Bank"))
data6.NO3 <- merge(data6.NO3, data3, by = c("Reach", "Site", "Distance", "Bank"))

#########################################################################################################################
######################################           Model Selection               #####################################
#########################################################################################################################
Function.MixedEffects <- Function.ModelSelection <- function(dataframe, variable, extra, extra2) {
  
  aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+extra2+OrgN+(1|Reach/Site), dataframe, REML=TRUE)), #1
                      AIC(lmer(variable~log_dist*Bank+extra2+extra+(1|Reach), dataframe, REML=TRUE)), #2
                      AIC(lmer(variable~Bank*log_dist+extra2+extra+(1|Site), dataframe, REML=TRUE))) #3
  return(aic.output)
                      
}

Function.MixedEffects(data3, data3$NetMin, data3$OrgN, data3$GW)
Function.MixedEffects(data3, data3$NetNit, data3$NH4_Conc,data3$GW)
Function.MixedEffects(data3, data3$NO3_Conc, data3$OrgN, data3$GW)
Function.MixedEffects(data3, data3$NH4_Conc, data3$OrgN, data3$GW)
Function.MixedEffects(data5.1, data5.1$d15N, data5.1$OrgN, data5.1$massN)
Function.MixedEffects(data5.1, data5.1$d13C, data5.1$OrgN, data5.1$massN)
Function.MixedEffects(data6.NO3, data6.NO3$NO3_Conc, data6.NO3$OrgN, data6.NO3$GW)
Function.MixedEffects(data6.NH4, data6.NH4$NH4_Conc, data6.NH4$OrgN, data6.NH4$GW)

Function.ModelSelection <- function(dataframe, variable, extra, extra2, mixed) {
  
  aic.output <- rbind(AIC(lm(variable~Bank*log_dist+extra2+extra, dataframe)), #1
  AIC(lm(variable~Bank*log_dist+extra, dataframe)), #2
  AIC(lm(variable~Bank*log_dist+extra2, dataframe)), #3
  AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank+extra2+extra, dataframe)), #4
  
  AIC(lmer(variable~Bank*log_dist+extra2+extra+(1|mixed), dataframe, REML=FALSE)), #5
  AIC(lmer(variable~Bank*log_dist+I(log_dist^2)*Bank+extra2+extra+(1|mixed), dataframe, REML=FALSE)), #6
  AIC(lmer(variable~Bank*log_dist+extra+(1|mixed), dataframe, REML=FALSE)), #7
  AIC(lmer(variable~Bank*log_dist+extra2+(1|mixed), dataframe, REML=FALSE)), #8
  AIC(lmer(variable~Bank+log_dist+extra2+(1|mixed), dataframe, REML=FALSE)), #9
  AIC(lmer(variable~I(log_dist^2)*Bank+extra+extra2+(1|mixed), dataframe, REML=FALSE)), #10
  AIC(lmer(variable~log_dist+extra+extra2+(1|mixed), dataframe, REML=FALSE)),#11
  AIC(lmer(variable~Bank+extra+extra2+(1|mixed), dataframe, REML=FALSE)),#12
  
  AIC(lm(variable~Bank*log_dist, dataframe)), #13
  AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank, dataframe)), #14
  AIC(lm(variable~log_dist+extra, dataframe)), #15
  AIC(lm(variable~log_dist+extra2, dataframe)), #16
  AIC(lm(variable~Bank+extra, dataframe)), #17
  AIC(lm(variable~Bank+extra2, dataframe)),#18
  AIC(lm(variable~extra, dataframe)), #19
  AIC(lm(variable~extra2, dataframe)) #20
  )
 
   names <- seq(1,20,1)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[4:20,])
  aic.output <- cbind(aic.output, delaic)
 colnames(aic.output)<- c("AIC", "delAIC")
   return(aic.output)
}



NetMin <- Function.ModelSelection(data3, data3$NetMin, data3$OrgN, data3$GW, data3$Reach)
NetNit <- Function.ModelSelection(data3, data3$NetNit, data3$NH4_Conc, data3$GW, data3$Reach)
NO3Conc <- Function.ModelSelection(data3, data3$NO3_Conc, data3$OrgN, data3$GW, data3$Reach)
NH4Conc <- Function.ModelSelection(data3, data3$NH4_Conc, data3$OrgN, data3$GW, data3$Reach)
d15N <- Function.ModelSelection(data5.1, data5.1$d15N, data5.1$OrgN, data5.1$massN, data5.1$Reach)
d13C <- Function.ModelSelection(data5.1, data5.1$d13C, data5.1$OrgN, data5.1$massN, data5.1$Site)
d15N.NO3 <- Function.ModelSelection(data6.NO3, data6.NO3$NO3_Conc, data6.NO3$OrgN, data6.NO3$GW, data6.NO3$Site)
d15N.NH4 <- Function.ModelSelection(data6.NH4, data6.NH4$NH4_Conc, data6.NH4$OrgN, data6.NH4$GW, data6.NH4$Site)

#########################################################################################################################
######################################            Net Mineralization                #####################################
#########################################################################################################################

data3 <-merge(data3, data5.1, by = c("Reach", "Site", "Distance", "Bank"))
log_dist <- log10(data3$Distance)
data3 <- cbind(data3, log_dist)

AIC(lmer(NetMin~Bank*log_dist+GW+OrgN+(1|Reach/Site), data3, REML=TRUE)) #1160.92
AIC(lmer(NetMin~log_dist*Bank+GW+OrgN+(1|Reach), data3, REML=TRUE)) #1159.014
AIC(lmer(NetMin~Bank*log_dist+GW+OrgN+(1|Site), data3, REML=TRUE)) #1159.859

AIC(lm(NetMin~Bank*log_dist+GW+OrgN, data3)) #1209.917
AIC(lm(NetMin~Bank*log_dist+OrgN, data3)) #1211.413***
AIC(lm(NetMin~Bank*log_dist+GW, data3)) #1213.996
AIC(lm(NetMin~Bank*log_dist+I(log_dist^2)*Bank+GW+OrgN, data3)) #1213.02

AIC(lmer(NetMin~Bank*log_dist+GW+OrgN+(1|Reach), data3, REML=FALSE)) #1211.553
AIC(lmer(NetMin~Bank*log_dist+I(log_dist^2)*Bank+GW+OrgN+(1|Reach), data3, REML=FALSE)) #1214.664
AIC(lmer(NetMin~Bank*log_dist+OrgN+(1|Reach), data3, REML=FALSE)) #1213.254
AIC(lmer(NetMin~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1215.914
AIC(lmer(NetMin~Bank+log_dist+GW+(1|Reach), data3, REML=FALSE)) #1217.068

Min.Model <- lm(NetMin~Bank*log_dist+OrgN, data3)
summary(Min.Model)

df <- data.frame("log_dist" = log(seq(0.1, 10, 0.1)), 
                "Bank" = c("Left", "Right"))
df.pred <- predict(Min.Model, type="response")
df.pred.se <- predict(Min.Model, type="response", se.fit=TRUE)$se.fit
plot(df.pred, df$log_dist)

qqnorm(resid(Min.Model), pch=16, main="Min")
qqline(resid(Min.Model))
shapiro.test(resid(Min.Model))$p.value


#########################################################################################################################
######################################            Net Nitrification                  #####################################
#########################################################################################################################

AIC(lmer(NetNit~Bank*log_dist+GW+NH4_Conc+(1|Reach/Site), data3, REML=TRUE)) #1112.985
AIC(lmer(NetNit~log_dist*Bank+GW+NH4_Conc+(1|Reach), data3, REML=TRUE)) #1114.196
AIC(lmer(NetNit~Bank*log_dist+GW+NH4_Conc+(1|Site), data3, REML=TRUE)) #11114.323

AIC(lm(NetNit~Bank*log_dist+GW+NH4_Conc, data3)) #1150.352***
AIC(lm(NetNit~Bank*log_dist+NH4_Conc, data3)) #1154.353
AIC(lm(NetNit~Bank*log_dist+GW, data3)) #1168.244
AIC(lm(NetNit~Bank*log_dist+I(log_dist^2)*Bank+GW+NH4_Conc, data3)) #1149.903

AIC(lmer(NetNit~Bank*log_dist+GW+NH4_Conc+(1|Reach), data3, REML=FALSE)) #1151.013
AIC(lmer(NetNit~Bank*log_dist+I(log_dist^2)*Bank+GW+NH4_Conc+(1|Reach), data3, REML=FALSE)) #1150.337
AIC(lmer(NetNit~Bank*log_dist+NH4_Conc+(1|Reach), data3, REML=FALSE)) #1155.067
AIC(lmer(NetNit~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1169.571



Nit.Model <- lm(NetNit~Bank*log_dist+GW+NH4_Conc, data3)
summary(Nit.Model)

qqnorm(resid(Nit.Model), pch=16, main="Nit")
qqline(resid(Nit.Model))
shapiro.test(resid(Nit.Model))$p.value



#########################################################################################################################
######################################            NH4 Concentration                #####################################
#########################################################################################################################


AIC(lmer(NH4_Conc~Bank*log_dist+GW+(1|Reach/Site), data3, REML=TRUE)) #1633.679
AIC(lmer(NH4_Conc~log_dist*Bank+GW+(1|Reach), data3, REML=TRUE)) #1631.683
AIC(lmer(NH4_Conc~Bank*log_dist+GW+(1|Site), data3, REML=TRUE)) #1631.683

AIC(lm(NH4_Conc~Bank*log_dist+GW, data3)) #1707.82
AIC(lm(NH4_Conc~Bank*log_dist, data3)) #1706.709***
AIC(lm(NH4_Conc~Bank*log_dist+GW, data3)) #1707.82
AIC(lm(NH4_Conc~Bank*log_dist+I(log_dist^2)*Bank+GW, data3)) #1709.758
AIC(lm(NH4_Conc~Bank*log_dist+I(log_dist^2)*Bank, data3)) ##1708.561
AIC(lm(NH4_Conc~Bank*log_dist+I(log_dist^2)*Bank+I(log_dist^3)*Bank, data3)) ##1712.163

AIC(lmer(NH4_Conc~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1709.82
AIC(lmer(NH4_Conc~Bank*log_dist+I(log_dist^2)*Bank+GW+(1|Reach), data3, REML=FALSE)) #1711.758
AIC(lmer(NH4_Conc~Bank*log_dist+(1|Reach), data3, REML=FALSE)) #1708.709
AIC(lmer(NH4_Conc~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1709.82



NH4.Model <- lm(NH4_Conc~Bank*log_dist+GW, data3)
summary(NH4.Model)

qqnorm(resid(NH4.Model), pch=16, main="Nit")
qqline(resid(NH4.Model))
shapiro.test(resid(NH4.Model))$p.value


#########################################################################################################################
######################################            N03 Concentration                 #####################################
#########################################################################################################################

AIC(lmer(NO3_Conc~Bank*log_dist+GW+(1|Reach/Site), data3, REML=TRUE)) #1313.044
AIC(lmer(NO3_Conc~log_dist*Bank+GW+(1|Reach), data3, REML=TRUE)) #1311.081
AIC(lmer(NO3_Conc~Bank*log_dist+GW+(1|Site), data3, REML=TRUE)) #1311.081

AIC(lm(NO3_Conc~Bank*log_dist+GW, data3)) #1365.088
AIC(lm(NO3_Conc~Bank*log_dist, data3)) #1384.036
AIC(lm(NO3_Conc~Bank*log_dist+GW, data3)) #1365.088
AIC(lm(NO3_Conc~Bank*log_dist+I(log_dist^2)*Bank+GW, data3)) #1364.88
AIC(lm(NO3_Conc~Bank*log_dist+I(log_dist^2)*Bank, data3)) #1353.367****

AIC(lmer(NO3_Conc~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1367.088
AIC(lmer(NO3_Conc~Bank*log_dist+I(log_dist^2)*Bank+GW+(1|Reach), data3, REML=FALSE)) #1366.88
AIC(lmer(NO3_Conc~Bank*log_dist+(1|Reach), data3, REML=FALSE)) #1386.036
AIC(lmer(NO3_Conc~Bank*log_dist+GW+(1|Reach), data3, REML=FALSE)) #1367.088



NO3.Model <- lm(NO3_Conc~Bank*log_dist+GW, data3)
summary(NO3.Model)

qqnorm(resid(NO3.Model), pch=16, main="Nit")
qqline(resid(NO3.Model))
shapiro.test(resid(NO3.Model))$p.value


############################################ ############################################ ############################################ 
############################################        d15N   ############################################ 
############################################ ############################################ ############################################ 
log_dist <- log10(data5.1$Distance)
data5.1 <- cbind(data5.1, log_dist)

AIC(lmer(d15N~Bank*log_dist+(1|Reach/Site), data5.1, REML=TRUE)) #335.0971
AIC(lmer(d15N~log_dist*Bank+(1|Reach), data5.1, REML=TRUE)) #333.0971
AIC(lmer(d15N~Bank*log_dist+(1|Site), data5.1, REML=TRUE)) #336.378

AIC(lm(d15N~Bank*log_dist, data5.1)) #334.3029
AIC(lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data5.1)) #327.384***

AIC(lmer(d15N~Bank*log_dist+(1|Reach), data5.1, REML=FALSE)) #333.9585
AIC(lmer(d15N~Bank*log_dist+I(log_dist^2)*Bank+(1|Reach), data5.1, REML=FALSE)) #325.946


d15N.Model <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data5.1)
summary(d15N.Model)

qqnorm(resid(NO3.Model), pch=16, main="Nit")
qqline(resid(NO3.Model))
shapiro.test(resid(NO3.Model))$p.value

############################################ ############################################ ############################################ 
############################################        d13C   ############################################ 
############################################ ############################################ ############################################ 
log_dist <- log10(data5.1$Distance)
data5.1 <- cbind(data5.1, log_dist)

AIC(lmer(d13C~Bank*log_dist+(1|Reach/Site), data5.1, REML=TRUE)) #132.6292
AIC(lmer(d13C~log_dist*Bank+(1|Reach), data5.1, REML=TRUE)) #133.2152
AIC(lmer(d13C~Bank*log_dist+(1|Site), data5.1, REML=TRUE)) #127.5222

AIC(lm(d15N~Bank*log_dist, data5.1)) #344.1004
AIC(lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data5.1)) #334.9355***

AIC(lmer(d15N~Bank*log_dist+(1|Reach), data5.1, REML=FALSE)) #343.5799
AIC(lmer(d15N~Bank*log_dist+I(log_dist^2)*Bank+(1|Reach), data5.1, REML=FALSE)) #333.544


d15N.Model <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data5.1)
summary(d15N.Model)

qqnorm(resid(NO3.Model), pch=16, main="Nit")
qqline(resid(NO3.Model))
shapiro.test(resid(NO3.Model))$p.value



############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################ ############################################ ############################################ 

data6 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Inorganic N Soil Extraction Data.csv") 
data6 <- data.frame(data6$Inorganic, data6$Test, data6$SampleID, data6$Reach, data6$ReachID, data6$Bank, data6$Distance, data6$PeakArea.N..Vs., data6$N.quantity..ug., data6$d15N.vs.Air.N2..permil.)
colnames(data6) <- c("Batch", "Test", "ID", "Reach", "Site", "Bank", "Distance", "Peak.Area", "N ug", "d15N")
log_dist <- log10(data6$Distance)
data6 <- cbind(data6, log_dist)
data6.1 <- subset(data6, Bank=="Right" | Bank=="Left")
data6.NH4 <- subset(data6.1, Test=="NH4")
data6.NO3 <- subset(data6.1, Test=="NO3")
write.csv(data6.NO3, file = "NO3.StableIsotope.csv")
write.csv(data6.NH4, file = "NH4.StableIsotope.csv")

############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################               NH4                            ############################################ 
############################################ ############################################ ############################################ 


AIC(lmer(d15N~Bank*log_dist+(1|Reach/Site), data6.NH4, REML=TRUE)) #512.8744
AIC(lmer(d15N~log_dist*Bank+(1|Reach), data6.NH4, REML=TRUE)) #510.8744
AIC(lmer(d15N~Bank*log_dist+(1|Site), data6.NH4, REML=TRUE)) #510.8533

AIC(lm(d15N~Bank*log_dist, data6.NH4)) #518.2149
AIC(lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NH4)) #512.1655**

AIC(lmer(d15N~Bank*log_dist+(1|Reach), data6.NH4, REML=FALSE)) #520.2149
AIC(lmer(d15N~Bank*log_dist+I(log_dist^2)*Bank+(1|Reach), data6.NH4, REML=FALSE)) #514.1655


d15N.NH4.Model <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NH4)
summary(d15N.NH4.Model)

qqnorm(resid(NO3.Model), pch=16, main="Nit")
qqline(resid(NO3.Model))
shapiro.test(resid(NO3.Model))$p.value

############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################               NO3                            ############################################ 
############################################ ############################################ ############################################ 



AIC(lmer(d15N~Bank*log_dist+(1|Reach/Site), data6.NO3, REML=TRUE)) #535.2171
AIC(lmer(d15N~log_dist*Bank+(1|Reach), data6.NO3, REML=TRUE)) #533.2171
AIC(lmer(d15N~Bank*log_dist+(1|Site), data6.NO3, REML=TRUE)) #533.6621

AIC(lm(d15N~Bank*log_dist, data6.NO3)) #541.3865
AIC(lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NO3)) #539.1965**
#AIC(lm(d15N~Bank*Distance+I(Distance^2)*Bank, data6.NO3)) #542.2418

AIC(lmer(d15N~Bank*log_dist+(1|Reach), data6.NO3, REML=FALSE)) #543.2973
AIC(lmer(d15N~Bank*log_dist+I(log_dist^2)*Bank+(1|Reach), data6.NO3, REML=FALSE)) #541.0478


d15N.NO3.Model <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NO3)
summary(d15N.NO3.Model)

qqnorm(resid(NO3.Model), pch=16, main="Nit")
qqline(resid(NO3.Model))
shapiro.test(resid(NO3.Model))$p.value

d15N.NO3 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NO3)
plot(d15N~log_dist, data=data6.NO3)
lines(sort(log_dist), fitted(d15N.NO3)[order(log_dist)], col='red', type='b')
lines(log_dist*Bank, fitted(d15N.NO3), col='red', type='b', data6.NO3)


#########################################################################################################################
######################################           Mass Nitrogen                #####################################
#########################################################################################################################
data.subsets <- c("Reach", "Site", "Bank", "Distance", "GW")
data.mass <- data4.5[c(data.subsets,"massN")]
massN <- cbind(data.mass, log_dist=log10(data.mass$Distance))

AIC(lmer(massN~Bank*log_dist+GW+(1|Reach/Site), massN, REML=TRUE)) #261.7137
AIC(lmer(massN~log_dist*Bank+GW+(1|Reach), massN, REML=TRUE)) #259.7137
AIC(lmer(massN~Bank*log_dist+GW+(1|Site), massN, REML=TRUE)) #259.7137

AIC(lm(massN~Bank*log_dist+GW, massN)) #249.5836****
AIC(lm(massN~Bank*log_dist, massN)) #307.495
AIC(lm(massN~Bank*log_dist+I(log_dist^2)*Bank+GW, massN)) #251.0279
AIC(lm(massN~Bank*log_dist+I(log_dist^2)*Bank, massN)) #309.7083

AIC(lmer(massN~Bank*log_dist+GW+(1|Reach), massN, REML=FALSE)) #1367.088
AIC(lmer(massN~Bank*log_dist+I(log_dist^2)*Bank+GW+(1|Reach), massN, REML=FALSE)) #1366.88
AIC(lmer(massN~Bank*log_dist+(1|Reach), massN, REML=FALSE)) #1386.036
AIC(lmer(massN~Bank*log_dist+GW+(1|Reach), massN, REML=FALSE)) #1367.088



NO3.Model <- lm(NO3_Conc~Bank*log_dist+GW, data3)
summary(NO3.Model)



######################################################################################################################
######################################        Supported models       #####################################
#####################################################################################################################

model.d15N <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data5.1)
model.d13C <- lm(d13C~Bank*log_dist+I(log_dist^2)*Bank, data5.1)
model.Min <- lm(NetMin~Bank*log_dist+OrgN, data3)
model.Nit <- lm(NetNit~Bank*log_dist+GW+NH4_Conc, data3)
model.NH4 <- lm(NH4_Conc~Bank*log_dist, data3)
model.NO3 <- lm(NO3_Conc~Bank*log_dist+I(log_dist^2)*Bank, data3)
model.d15NNO3 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NO3)
model.d15NNH4 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NH4)
model.massN <- lm(massN~Bank*log_dist+GW, massN)

summary(model.d15N)
anova(model.d15N)
res.d15N <-cbind(Estimate = coef(summary(model.d15N))[2:6,1], anova(model.d15N)[1:5,], test=rep("Bulk d15N",5)) ##distance effect, and bank effect
write.csv(res.d15N, file = "Bulk.d15N.csv")

summary(model.d13C)
anova(model.d13C)
res.d13C <-cbind(Estimate = coef(summary(model.d13C))[2:6,1], anova(model.d13C)[1:5,], test=rep("Bulk d13C",5))
write.csv(res.d13C, file = "Bulk.d13C.csv")

summary(model.Min)
anova(model.Min) ##no distance or bank effect
res.Min <-cbind(Estimate = coef(summary(model.Min))[2:5,1], anova(model.Min)[1:4,], test=rep("Net Mineralization",4))
write.csv(res.Min, file = "NetMineralization.csv")

summary(model.Nit)
anova(model.Nit)  ###no bank or distance effect (GW an dNH4 conc effect)
res.Nit <-cbind(Estimate = coef(summary(model.Nit))[2:6,1], anova(model.Nit)[1:5,], test=rep("Net Nitrification",5))
write.csv(res.Nit, file = "NetNitrification.csv")

summary(model.NH4)
anova(model.NH4) ##bank effect, no distance effect
res.NH4 <-cbind(Estimate = coef(summary(model.NH4))[2:4,1], anova(model.NH4)[1:3,], test=rep("NH4 Concentration",3))
write.csv(res.NH4, file = "NH4.Concentration.csv")

summary(model.NO3)
anova(model.NO3) ## bank effect, no distance effect
res.NO3 <-cbind(Estimate = coef(summary(model.NO3))[2:6,1], anova(model.NO3)[1:5,], test=rep("NO3 Concentration",5))
write.csv(res.NO3, file = "NO3.Concentration.csv")

summary(model.d15NNH4)
anova(model.d15NNH4) ##bank and distance effect
res.d15NNH4 <-cbind(Estimate = coef(summary(model.d15NNH4))[2:6,1], anova(model.d15NNH4)[1:5,], test=rep("d15N NH4",5))
write.csv(res.d15NNH4, file = "d15N.NH4.csv")

summary(model.d15NNO3)
anova(model.d15NNO3) ##bank and distance effect
res.d15NNO3 <-cbind(Estimate = coef(summary(model.d15NNO3))[2:6,1], anova(model.d15NNH4)[1:5,], test=rep("d15N NO3",5))
write.csv(res.d15NNO3, file = "d15N.NO3.csv")

summary(model.massN)
anova(model.massN) ##bank and distance effect
res.massN <-cbind(Estimate = coef(summary(model.massN))[2:6,1], anova(model.massN)[1:5,], test=rep("d15N NO3",5))
write.csv(res.massN, file = "massNres.csv")

model.output <- rbind(res.NH4, res.NO3, res.Nit, res.Min, res.d13C, res.d15N, res.d15NNH4, res.d15NNO3)
write.csv(model.output, file = "Model.Output.csv")



######################################################################################################################
######################################        Salmon contribution to soil        #####################################
#####################################################################################################################
subset(data6.NO3, d15N>12)
subset(data6.NH4, d15N>12)
