rm(list = ls())

library(readr)
Modified_HansenSoil_KClData_24Aug2017 <- read.csv("Data/Raw/Modified_HansenSoil_KClData_24Aug2017.csv") 

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
write.csv(data3, file = "Data/Clean/Calculatd_Data.csv")

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




##################################################################################################################################
#####################################       Linear Mixed Effects Model      #######################################################
#################################################################################################################################
library(nlme)
library(lme4)
library(lmerTest)
library(AICcmodavg)

############################################ ############################################ ############################################ 
############################################        Organic Nitrogen   ############################################ 
############################################ ############################################ ############################################ 


data1m.O <- data5.1[which(data5.1$Distance==1),]
nest.Full.1m.OrgN <- lmer(OrgN~Bank+(1|Reach/Site), data1m.O, REML=0)

data3m.O <- data5.1[which(data5.1$Distance==3),]
nest.Full.3m.OrgN <- lmer(OrgN~Bank+(1|Reach/Site), data3m.O, REML=0)

data6m.O <- data5.1[which(data5.1$Distance==6),]
nest.Full.6m.OrgN <- lmer(OrgN~Bank+(1|Reach/Site), data6m.O, REML=0)

data10m.O <- data5.1[which(data5.1$Distance==10),]
nest.Full.10m.OrgN <- lmer(OrgN~Bank+(1|Reach/Site), data10m.O, REML=0)

data20m.O <- data5.1[which(data5.1$Distance==20),]
nest.Full.20m.OrgN <- lmer(OrgN~Bank+(1|Reach/Site), data20m.O, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.OrgN))[1,1], coef(summary(nest.Full.3m.OrgN))[1,1], coef(summary(nest.Full.6m.OrgN))[1,1], coef(summary(nest.Full.10m.OrgN))[1,1], coef(summary(nest.Full.20m.OrgN))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.OrgN))[1,2], coef(summary(nest.Full.3m.OrgN))[1,2], coef(summary(nest.Full.6m.OrgN))[1,2], coef(summary(nest.Full.10m.OrgN))[1,2], coef(summary(nest.Full.20m.OrgN))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.OrgN))[1,5], coef(summary(nest.Full.3m.OrgN))[1,5], coef(summary(nest.Full.6m.OrgN))[1,5], coef(summary(nest.Full.10m.OrgN))[1,5], coef(summary(nest.Full.20m.OrgN))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.OrgN))[2,1], coef(summary(nest.Full.3m.OrgN))[2,1], coef(summary(nest.Full.6m.OrgN))[2,1], coef(summary(nest.Full.10m.OrgN))[2,1], coef(summary(nest.Full.20m.OrgN))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.OrgN))[2,2], coef(summary(nest.Full.3m.OrgN))[2,2], coef(summary(nest.Full.6m.OrgN))[2,2], coef(summary(nest.Full.10m.OrgN))[2,2], coef(summary(nest.Full.20m.OrgN))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.OrgN))[2,5], coef(summary(nest.Full.3m.OrgN))[2,5], coef(summary(nest.Full.6m.OrgN))[2,5], coef(summary(nest.Full.10m.OrgN))[2,5], coef(summary(nest.Full.20m.OrgN))[2,5])

#OrgN <- cbind(OrgN, DelAIC.Full, DelAICc.Full)
Analysis <- rep("Organic Nitrogen", 5)
OrgN  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(OrgN) <- names
row.names(OrgN) <- names2
write.csv(OrgN, file = "OrgN.csv")



#########################################################################################################################
######################################            Net Mineralization                #####################################
#########################################################################################################################

data3 <-merge(data3, data5.1, by = c("Reach", "Site", "Distance", "Bank"))

data1m <- data3[which(data3$Distance==1),]
nest.Full.1m.min <- lmer(NetMin~Bank+GW+OrgN+(1|Reach/Site), data1m, REML=0)


data3m <- data3[which(data3$Distance==3),]
nest.Full.3m.min <- lmer(NetMin~Bank+GW+OrgN(1|Reach/Site), data3m, REML=0)


data6m <- data3[which(data3$Distance==6),]
nest.Full.6m.min <- lmer(NetMin~Bank+GW+OrgN(1|Reach/Site), data6m, REML=0)


data10m <- data3[which(data3$Distance==10),]
nest.Full.10m.min <- lmer(NetMin~Bank+GW+OrgN(1|Reach/Site), data10m, REML=0)


data20m <- data3[which(data3$Distance==20),]
nest.Full.20m.min <- lmer(NetMin~Bank+GW+OrgN(1|Reach/Site), data20m, REML=0)


library(lmerTest)

EstimateBankL <- rbind(coef(summary(nest.Full.1m.min))[1,1], coef(summary(nest.Full.3m.min))[1,1], coef(summary(nest.Full.6m.min))[1,1], coef(summary(nest.Full.10m.min))[1,1], coef(summary(nest.Full.20m.min))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.min))[1,2], coef(summary(nest.Full.3m.min))[1,2], coef(summary(nest.Full.6m.min))[1,2], coef(summary(nest.Full.10m.min))[1,2], coef(summary(nest.Full.20m.min))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.min))[1,5], coef(summary(nest.Full.3m.min))[1,5], coef(summary(nest.Full.6m.min))[1,5], coef(summary(nest.Full.10m.min))[1,5], coef(summary(nest.Full.20m.min))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.min))[2,1], coef(summary(nest.Full.3m.min))[2,1], coef(summary(nest.Full.6m.min))[2,1], coef(summary(nest.Full.10m.min))[2,1], coef(summary(nest.Full.20m.min))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.min))[2,2], coef(summary(nest.Full.3m.min))[2,2], coef(summary(nest.Full.6m.min))[2,2], coef(summary(nest.Full.10m.min))[2,2], coef(summary(nest.Full.20m.min))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.min))[2,5], coef(summary(nest.Full.3m.min))[2,5], coef(summary(nest.Full.6m.min))[2,5], coef(summary(nest.Full.10m.min))[2,5], coef(summary(nest.Full.20m.min))[2,5])

#Mineralization <- cbind(Mineralization, DelAIC.Full, DelAICc.Full)
Analysis <- rep("Net Mineralization", 5)
Mineralization <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(Mineralization) <- names
row.names(Mineralization) <- names2

############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 


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


#########################################################################################################################
######################################            Net Nitrification                  #####################################
#########################################################################################################################

nest.Full.1m.nit <- lmer(NetNit~Bank+GW+NH4_Conc+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.nit <- lmer(NetNit~Bank+GW+NH4_Conc+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.nit <- lmer(NetNit~Bank+GW+NH4_Conc+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.nit <- lmer(NetNit~Bank+GW+NH4_Conc+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.nit <- lmer(NetNit~Bank+GW+NH4_Conc+(1|Reach/Site), data20m, REML=0)


EstimateBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,1], coef(summary(nest.Full.3m.nit))[1,1], coef(summary(nest.Full.6m.nit))[1,1], coef(summary(nest.Full.10m.nit))[1,1], coef(summary(nest.Full.20m.nit))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,2], coef(summary(nest.Full.3m.nit))[1,2], coef(summary(nest.Full.6m.nit))[1,2], coef(summary(nest.Full.10m.nit))[1,2], coef(summary(nest.Full.20m.nit))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.nit))[1,5], coef(summary(nest.Full.3m.nit))[1,5], coef(summary(nest.Full.6m.nit))[1,5], coef(summary(nest.Full.10m.nit))[1,5], coef(summary(nest.Full.20m.nit))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,1], coef(summary(nest.Full.3m.nit))[2,1], coef(summary(nest.Full.6m.nit))[2,1], coef(summary(nest.Full.10m.nit))[2,1], coef(summary(nest.Full.20m.nit))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,2], coef(summary(nest.Full.3m.nit))[2,2], coef(summary(nest.Full.6m.nit))[2,2], coef(summary(nest.Full.10m.nit))[2,2], coef(summary(nest.Full.20m.nit))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.nit))[2,5], coef(summary(nest.Full.3m.nit))[2,5], coef(summary(nest.Full.6m.nit))[2,5], coef(summary(nest.Full.10m.nit))[2,5], coef(summary(nest.Full.20m.nit))[2,5])

#Nitrification <- cbind(AIC.Full, DelAIC.Full, DelAICc.Full)
Analysis <- rep("Net Nitrification", 5)
Nitrification <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(Nitrification) <- names
row.names(Nitrification) <- names2

############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 

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



#########################################################################################################################
######################################            NH4 Concentration                #####################################
#########################################################################################################################

nest.Full.1m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.NH4 <- lmer(NH4_Conc~Bank+GW+(1|Reach/Site), data20m, REML=0)


EstimateBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,1], coef(summary(nest.Full.3m.NH4))[1,1], coef(summary(nest.Full.6m.NH4))[1,1], coef(summary(nest.Full.10m.NH4))[1,1], coef(summary(nest.Full.20m.NH4))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,2], coef(summary(nest.Full.3m.NH4))[1,2], coef(summary(nest.Full.6m.NH4))[1,2], coef(summary(nest.Full.10m.NH4))[1,2], coef(summary(nest.Full.20m.NH4))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.NH4))[1,5], coef(summary(nest.Full.3m.NH4))[1,5], coef(summary(nest.Full.6m.NH4))[1,5], coef(summary(nest.Full.10m.NH4))[1,5], coef(summary(nest.Full.20m.NH4))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,1], coef(summary(nest.Full.3m.NH4))[2,1], coef(summary(nest.Full.6m.NH4))[2,1], coef(summary(nest.Full.10m.NH4))[2,1], coef(summary(nest.Full.20m.NH4))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,2], coef(summary(nest.Full.3m.NH4))[2,2], coef(summary(nest.Full.6m.NH4))[2,2], coef(summary(nest.Full.10m.NH4))[2,2], coef(summary(nest.Full.20m.NH4))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.NH4))[2,5], coef(summary(nest.Full.3m.NH4))[2,5], coef(summary(nest.Full.6m.NH4))[2,5], coef(summary(nest.Full.10m.NH4))[2,5], coef(summary(nest.Full.20m.NH4))[2,5])

#NH4Concentration <- cbind(NH4Concentration, DelAIC.Full)
Analysis <- rep("NH4 Concentration")
NH4Concentration  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(NH4Concentration ) <- names
row.names(NH4Concentration) <- names2

############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 

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



#########################################################################################################################
######################################            N03 Concentration                 #####################################
#########################################################################################################################

nest.Full.1m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data1m, REML=0)

nest.Full.3m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data3m, REML=0)

nest.Full.6m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data6m, REML=0)

nest.Full.10m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data10m, REML=0)

nest.Full.20m.N03 <- lmer(N03_Conc~Bank+GW+(1|Reach/Site), data20m, REML=0)


EstimateBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,1], coef(summary(nest.Full.3m.N03))[1,1], coef(summary(nest.Full.6m.N03))[1,1], coef(summary(nest.Full.10m.N03))[1,1], coef(summary(nest.Full.20m.N03))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,2], coef(summary(nest.Full.3m.N03))[1,2], coef(summary(nest.Full.6m.N03))[1,2], coef(summary(nest.Full.10m.N03))[1,2], coef(summary(nest.Full.20m.N03))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.N03))[1,5], coef(summary(nest.Full.3m.N03))[1,5], coef(summary(nest.Full.6m.N03))[1,5], coef(summary(nest.Full.10m.N03))[1,5], coef(summary(nest.Full.20m.N03))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,1], coef(summary(nest.Full.3m.N03))[2,1], coef(summary(nest.Full.6m.N03))[2,1], coef(summary(nest.Full.10m.N03))[2,1], coef(summary(nest.Full.20m.N03))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,2], coef(summary(nest.Full.3m.N03))[2,2], coef(summary(nest.Full.6m.N03))[2,2], coef(summary(nest.Full.10m.N03))[2,2], coef(summary(nest.Full.20m.N03))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.N03))[2,5], coef(summary(nest.Full.3m.N03))[2,5], coef(summary(nest.Full.6m.N03))[2,5], coef(summary(nest.Full.10m.N03))[2,5], coef(summary(nest.Full.20m.N03))[2,5])

#N03Concentration <- cbind(N03Concentration, DelAIC.Full)
Analysis <- rep("NO3 Concentration", 5)
NO3Concentration  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(NO3Concentration ) <- names
row.names(NO3Concentration) <- names2


############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 

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


write.csv(Mineralization, file = "Mineralization.csv")
write.csv(Nitrification, file = "Nitrification.csv")
write.csv(N03Concentration, file = "N03Concentration.csv")
write.csv(NH4Concentration, file = "NH4Concentration.csv")



############################################ ############################################ ############################################ 
############################################        d15N   ############################################ 
############################################ ############################################ ############################################ 


data1m <- data5.1[which(data5.1$Distance==1),]
nest.Full.1m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data1m, REML=0)

data3m <- data5.1[which(data5.1$Distance==3),]
nest.Full.3m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data3m, REML=0)

data6m <- data5.1[which(data5.1$Distance==6),]
nest.Full.6m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data6m, REML=0)

data10m <- data5.1[which(data5.1$Distance==10),]
nest.Full.10m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data10m, REML=0)

data20m <- data4.1[which(data5.1$Distance==20),]
nest.Full.20m.d15N <- lmer(d15N~Bank+(1|Reach/Site), data20m, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,1], coef(summary(nest.Full.3m.d15N))[1,1], coef(summary(nest.Full.6m.d15N))[1,1], coef(summary(nest.Full.10m.d15N))[1,1], coef(summary(nest.Full.20m.d15N))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,2], coef(summary(nest.Full.3m.d15N))[1,2], coef(summary(nest.Full.6m.d15N))[1,2], coef(summary(nest.Full.10m.d15N))[1,2], coef(summary(nest.Full.20m.d15N))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.d15N))[1,5], coef(summary(nest.Full.3m.d15N))[1,5], coef(summary(nest.Full.6m.d15N))[1,5], coef(summary(nest.Full.10m.d15N))[1,5], coef(summary(nest.Full.20m.d15N))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,1], coef(summary(nest.Full.3m.d15N))[2,1], coef(summary(nest.Full.6m.d15N))[2,1], coef(summary(nest.Full.10m.d15N))[2,1], coef(summary(nest.Full.20m.d15N))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,2], coef(summary(nest.Full.3m.d15N))[2,2], coef(summary(nest.Full.6m.d15N))[2,2], coef(summary(nest.Full.10m.d15N))[2,2], coef(summary(nest.Full.20m.d15N))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.d15N))[2,5], coef(summary(nest.Full.3m.d15N))[2,5], coef(summary(nest.Full.6m.d15N))[2,5], coef(summary(nest.Full.10m.d15N))[2,5], coef(summary(nest.Full.20m.d15N))[2,5])

#d15N <- cbind(d15N, DelAIC.Full, DelAICc.Full)
Analysis <- rep("Bulk d15N", 5)
d15N  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(d15N) <- names
row.names(d15N) <- names2
write.csv(d15N, file = "d15N.csv")





############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 

qqnorm(resid(nest.Full.1m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.1m.OrgN))

qqnorm(resid(nest.Full.3m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.3m.OrgN))

qqnorm(resid(nest.Full.6m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.6m.OrgN))

qqnorm(resid(nest.Full.10m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.10m.OrgN))

qqnorm(resid(nest.Full.20m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.20m.OrgN))

OrgN.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.OrgN))$p.value, shapiro.test(resid(nest.Full.3m.OrgN))$p.value,
                   shapiro.test(resid(nest.Full.6m.OrgN))$p.value, shapiro.test(resid(nest.Full.10m.OrgN))$p.value,
                   shapiro.test(resid(nest.Full.20m.OrgN))$p.value)




Mineralization
Nitrification
d15N
N03Concentration
NH4Concentration
OrgN



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




############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################ ############################################ ############################################ 

data6 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Inorganic N Soil Extraction Data.csv") 
data6 <- data.frame(data6$Inorganic, data6$Test, data6$SampleID, data6$Reach, data6$ReachID, data6$Bank, data6$Distance, data6$PeakArea.N..Vs., data6$N.quantity..ug., data6$d15N.vs.Air.N2..permil.)
colnames(data6) <- c("Batch", "Test", "ID", "Reach", "Site", "Bank", "Distance", "Peak.Area", "N ug", "d15N")
data6.1 <- subset(data6, Bank=="Right" | Bank=="Left")
data6.NH4 <- subset(data6.1, Test=="NH4")
data6.NO3 <- subset(data6.1, Test=="NO3")
write.csv(data6.NO3, file = "NO3.StableIsotope.csv")
write.csv(data6.NH4, file = "NH4.StableIsotope.csv")

############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################               NH4                            ############################################ 
############################################ ############################################ ############################################ 

data1m <- data6.NH4[which(data6.NH4$Distance==1),]
nest.Full.1m.NH4d15N <- lmer(d15N~Bank+(1|Reach/Site), data1m, REML=0)

data3m <- data6.NH4[which(data6.NH4$Distance==3),]
nest.Full.3m.NH4d15N <- lmer(d15N~Bank+(1|Reach/Site), data3m, REML=0)

data6m <- data6.NH4[which(data6.NH4$Distance==6),]
nest.Full.6m.NH4d15N <- lmer(d15N~Bank+(1|Reach/Site), data6m, REML=0)

data10m <- data6.NH4[which(data6.NH4$Distance==10),]
nest.Full.10m.NH4d15N <- lmer(d15N~Bank+(1|Reach/Site), data10m, REML=0)

data20m <- data6.NH4[which(data6.NH4$Distance==20),]
nest.Full.20m.NH4d15N <- lmer(d15N~Bank+(1|Reach/Site), data20m, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.NH4d15N))[1,1], coef(summary(nest.Full.3m.NH4d15N))[1,1], coef(summary(nest.Full.6m.NH4d15N))[1,1], coef(summary(nest.Full.10m.NH4d15N))[1,1], coef(summary(nest.Full.20m.NH4d15N))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.NH4d15N))[1,2], coef(summary(nest.Full.3m.NH4d15N))[1,2], coef(summary(nest.Full.6m.NH4d15N))[1,2], coef(summary(nest.Full.10m.NH4d15N))[1,2], coef(summary(nest.Full.20m.NH4d15N))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.NH4d15N))[1,5], coef(summary(nest.Full.3m.NH4d15N))[1,5], coef(summary(nest.Full.6m.NH4d15N))[1,5], coef(summary(nest.Full.10m.NH4d15N))[1,5], coef(summary(nest.Full.20m.NH4d15N))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.NH4d15N))[2,1], coef(summary(nest.Full.3m.NH4d15N))[2,1], coef(summary(nest.Full.6m.NH4d15N))[2,1], coef(summary(nest.Full.10m.NH4d15N))[2,1], coef(summary(nest.Full.20m.NH4d15N))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.NH4d15N))[2,2], coef(summary(nest.Full.3m.NH4d15N))[2,2], coef(summary(nest.Full.6m.NH4d15N))[2,2], coef(summary(nest.Full.10m.NH4d15N))[2,2], coef(summary(nest.Full.20m.NH4d15N))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.NH4d15N))[2,5], coef(summary(nest.Full.3m.NH4d15N))[2,5], coef(summary(nest.Full.6m.NH4d15N))[2,5], coef(summary(nest.Full.10m.NH4d15N))[2,5], coef(summary(nest.Full.20m.NH4d15N))[2,5])

#NH4d15N <- cbind(NH4d15N, DelAIC.Full, DelAICc.Full)
Analysis <- rep("NH4 d15N", 5)
NH4d15N  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(NH4d15N) <- names
row.names(NH4d15N) <- names2
write.csv(NH4d15N, file = "NH4d15N.csv")


############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################               NO3                            ############################################ 
############################################ ############################################ ############################################ 

data1m <- data6.NO3[which(data6.NO3$Distance==1),]
nest.Full.1m.NO3d15N <- lmer(d15N~Bank+(1|Reach/Site), data1m, REML=0)

data3m <- data6.NO3[which(data6.NO3$Distance==3),]
nest.Full.3m.NO3d15N <- lmer(d15N~Bank+(1|Reach/Site), data3m, REML=0)

data6m <- data6.NO3[which(data6.NO3$Distance==6),]
nest.Full.6m.NO3d15N <- lmer(d15N~Bank+(1|Reach/Site), data6m, REML=0)

data10m <- data6.NO3[which(data6.NO3$Distance==10),]
nest.Full.10m.NO3d15N <- lmer(d15N~Bank+(1|Reach/Site), data10m, REML=0)

data20m <- data6.NO3[which(data6.NO3$Distance==20),]
nest.Full.20m.NO3d15N <- lmer(d15N~Bank+(1|Reach/Site), data20m, REML=0)



EstimateBankL <- rbind(coef(summary(nest.Full.1m.NO3d15N))[1,1], coef(summary(nest.Full.3m.NO3d15N))[1,1], coef(summary(nest.Full.6m.NO3d15N))[1,1], coef(summary(nest.Full.10m.NO3d15N))[1,1], coef(summary(nest.Full.20m.NO3d15N))[1,1])
StdBankL <- rbind(coef(summary(nest.Full.1m.NO3d15N))[1,2], coef(summary(nest.Full.3m.NO3d15N))[1,2], coef(summary(nest.Full.6m.NO3d15N))[1,2], coef(summary(nest.Full.10m.NO3d15N))[1,2], coef(summary(nest.Full.20m.NO3d15N))[1,2])
pBankL <- rbind(coef(summary(nest.Full.1m.NO3d15N))[1,5], coef(summary(nest.Full.3m.NO3d15N))[1,5], coef(summary(nest.Full.6m.NO3d15N))[1,5], coef(summary(nest.Full.10m.NO3d15N))[1,5], coef(summary(nest.Full.20m.NO3d15N))[1,5])

EstimateBankR <- rbind(coef(summary(nest.Full.1m.NO3d15N))[2,1], coef(summary(nest.Full.3m.NO3d15N))[2,1], coef(summary(nest.Full.6m.NO3d15N))[2,1], coef(summary(nest.Full.10m.NO3d15N))[2,1], coef(summary(nest.Full.20m.NO3d15N))[2,1])
StdBankR <- rbind(coef(summary(nest.Full.1m.NO3d15N))[2,2], coef(summary(nest.Full.3m.NO3d15N))[2,2], coef(summary(nest.Full.6m.NO3d15N))[2,2], coef(summary(nest.Full.10m.NO3d15N))[2,2], coef(summary(nest.Full.20m.NO3d15N))[2,2])
pBankR <- rbind(coef(summary(nest.Full.1m.NO3d15N))[2,5], coef(summary(nest.Full.3m.NO3d15N))[2,5], coef(summary(nest.Full.6m.NO3d15N))[2,5], coef(summary(nest.Full.10m.NO3d15N))[2,5], coef(summary(nest.Full.20m.NO3d15N))[2,5])

#NO3d15N <- cbind(NO3d15N, DelAIC.Full, DelAICc.Full)
Analysis <- rep("NO3 d15N", 5)
NO3d15N  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
colnames(NO3d15N) <- names
row.names(NO3d15N) <- names2
write.csv(NO3d15N, file = "NO3d15N.csv")




############################################ ############################################ ############################################ 
############################################        testing the residuals   ############################################ 
############################################ ############################################ ############################################ 

qqnorm(resid(nest.Full.1m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.1m.OrgN))

qqnorm(resid(nest.Full.3m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.3m.OrgN))

qqnorm(resid(nest.Full.6m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.6m.OrgN))

qqnorm(resid(nest.Full.10m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.10m.OrgN))

qqnorm(resid(nest.Full.20m.OrgN), pch=16, main="OrgN")
qqline(resid(nest.Full.20m.OrgN))

OrgN.Shap.Full <- c(shapiro.test(resid(nest.Full.1m.OrgN))$p.value, shapiro.test(resid(nest.Full.3m.OrgN))$p.value,
                    shapiro.test(resid(nest.Full.6m.OrgN))$p.value, shapiro.test(resid(nest.Full.10m.OrgN))$p.value,
                    shapiro.test(resid(nest.Full.20m.OrgN))$p.value)




Mineralization
Nitrification
d15N
N03Concentration
NH4Concentration
OrgN
NO3d15N
NH4d15N

model.output <- rbind(Mineralization, Nitrification, d15N, NO3Concentration, NH4Concentration, OrgN, NO3d15N, NH4d15N)
colnames(model.output) <- names
write.csv(model.output, file = "Model Output.csv")

######################################################################################################################
######################################        Salmon contribution to soil        #####################################
#####################################################################################################################


