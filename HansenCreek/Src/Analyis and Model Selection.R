rm(list = ls())
setwd("~/Documents/Graduate Work 2/Hansen Creek/Results")

library(readr)
Modified_HansenSoil_KClData_24Aug2017 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Modified_HansenSoil_KClData_24Aug2017.csv") 

data1 <- Modified_HansenSoil_KClData_24Aug2017
colnames(data1) <- c("Sample", "ID", "Lab", "NH4", "NO3", "Side", "Distance", "Transect", "Incubation", "Location", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin")

tin_weights <- c(3.081, 3.098, 3.118, 3.123, 3.078, 3.107, 3.091, 3.129, 3.108)
mean_tin <- mean(tin_weights)
days <- 15 #number of days incubated

KCL_mL <- 100  #100 mL of 2N KCl added to all soil samples. 

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

# GWC = gravimetric water content = (wet - dry)/dry
GW <- (Wet_Weight - Dry_Weight)/Dry_Weight

# WetToDry = Mass_wet / Mass_dry
WetToDry <- Wet_Weight / Dry_Weight

# Calculation taken from method 4D9 of Kellogg Soil Survey Laboratory Methods Manual (pg. 364)
NH4_Conc <- NH4_Corrected * 1000 * (KCL_mL/1000) * WetToDry / data1$Initial.Wet 
NO3_Conc <- NO3_Corrected*1000*(KCL_mL/1000) * WetToDry / data1$Initial.Wet

data2<- cbind(data1, NH4_Corrected, NO3_Corrected, Wet_Weight, Dry_Weight,  GW, NH4_Conc, NO3_Conc)

data2.1 <- subset(data2, !(ID == "Blank1 " | ID == "Blank2" | ID == "Blank3" | ID == "Blank1F" | ID == "Blank2F" | ID == "Blank3F"))
length(data2.1$Sample)
Initial <- subset(data2.1, Incubation ==  "I")
Final <- subset(data2.1, Incubation ==  "F")
length(Initial$Sample)
length(Final$Sample)
Net_Mineralization <- ((Final$NH4_Conc+Final$NO3_Conc)-(Initial$NH4_Conc+Initial$NO3_Conc))/days
Net_Nitrification <- (Final$NO3_Conc-Initial$NO3_Conc)/days

data3 <- cbind(Initial, Net_Mineralization, Net_Nitrification)
colnames(data3)<- c("Sample", "ID", "Lab", "NH4", "NO3", "Bank", "Distance", "Site", "Incubation","Reach", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin", "NH4_Corr", "NO3_Corr", "Wet_Weight", "Dry_Weight", "GW", "NH4_Conc", "NO3_Conc", "NetMin", "NetNit")
write.csv(data3, file = "Calculatd_Data.csv")

data3 <-data3[-79,]
data3 <-data3[-11,]
write.csv(data3, file = "Calculatd_Data.csv")

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
C.N <- (data4.11$pctC/data4.11$pctN)*(12.0107/14.0057)
data4.11 <- cbind(data4.11, C.N)

massN <- (data4.11$Weight * data4.11$pctN) #/data4.11$Weight
massN2 <- data4.11$pctN*10
pctN <- data4.11

data4.2 <- cbind(data4.11, massN, massN2)

data4.3 <- Initial
colnames(data4.3) <- c("Sample", "ID", "Lab", "NH4", "NO3", "Bank", "Distance", "Site", "Incubation", "Reach", "Initial.Wet", "Soil.Wet.Tin", "Soil.Dry.Tin", "NH4Corrected", "NO3Corrected", "Dry_Weight", "Water_Weight", "GW", "NH4_Conc", "NO3_Conc")
data4.4 <- merge(data4.3, data4.2, by = c("Reach", "Site", "Distance", "Bank"))



mass.norm <- (data4.4$pctN)*(data4.4$Dry_Weight)
data4.5 <- cbind(mass.norm, data4.4)





write.csv(data4.5, file="MassNitrogen.csv")

data4.3 <- merge(data4.2, data3, by = c("Reach", "Site", "Distance", "Bank"))
massNH4 <- data4.3$NH4_Conc/100000*data4.3$Weight
massNO3 <- data4.3$NO3_Conc/100000*data4.3$Weight
OrgN <- (data4.3$massN-massNH4-massNO3)*1000

#Gravimetri
data4.4 <- cbind(data4.3, massNH4, massNO3, OrgN)
data3 <-data3[-74,]
data5 <- data.frame(data4.4$ID.x, data4.4$Reach, data4.4$Site, data4.4$Distance, data4.4$Bank, 
                    data4.4$massN, data4.4$massNH4, data4.4$massNO3, data4.4$OrgN, data4.4$C.N)
colnames(data5) <- c("ID", "Reach", "Site", "Distance", "Bank", "massN", "massNH4", "massNO3", "OrgN", "C.N")

write.csv(data5, file = "OrganicNitrogen.csv")

data5.1 <- merge(data4.1, data5, by = c("Reach", "Site", "Distance", "Bank", "ID"))
#merge(data3[, c("pctN", "d15N","OrgN")], data3, by = c("Reach", "Site", "Distance", "Bank", "ID"))



##################################################################################################################################
#####################################       Linear Mixed Effects Model      #######################################################
#################################################################################################################################
library(nlme)
library(lme4)
library(lmerTest)
library(AICcmodavg)


#OrgN <- cbind(OrgN, DelAIC.Full, DelAICc.Full)
#Analysis <- rep("Organic Nitrogen", 5)
#OrgN  <- cbind(Analysis, EstimateBankL, StdBankL, pBankL, EstimateBankR, StdBankR, pBankR)
#colnames(OrgN) <- names
#row.names(OrgN) <- names2
#write.csv(OrgN, file = "OrgN.csv")

log_dist <- log(data3$Distance)
data3 <- cbind(data3, log_dist)


data3 <-merge(data5.1, data3, by = c("Reach", "Site", "Distance", "Bank"))



############################################ ############################################ ############################################ 
############################################        Inorganic Nitrogen Stable Isotope   ############################################ 
############################################ ############################################ ############################################ 

data6 <- read.csv("~/Documents/Graduate Work 2/Hansen Creek/Data/Raw Data/Inorganic N Soil Extraction Data 2.csv") 
NH4.std <- mean(-0.13851, -0.88937, -1.0352, -0.87689, -0.99886) #-0.1385
NO3.std <- mean(1.0494, 1.3199, 1.0955, 1.3511, 0.88867) #1.0494


############################################ ############################################ ############################################ 
############################################        Mixing Model                         ############################################ 
############################################ ############################################ ############################################

blankNH4 <- subset(data6, Reach=="BLANK" & Test=="NH4")

Batch1.NH4.Blank <- mean(subset(blankNH4, Batch==1)$d15N.vs.Air.N2..permil.) # 3.3648
Batch2.NH4.Blank <- mean(subset(blankNH4, Batch==2)$d15N.vs.Air.N2..permil.) # 9.2983
Batch3.NH4.Blank <- mean(subset(blankNH4, Batch==3)$d15N.vs.Air.N2..permil.) # -0.4890
Batch4.NH4.Blank <- mean(subset(blankNH4, Batch==4)$d15N.vs.Air.N2..permil.) # 5.3567

Batch1.NH4.BlankMass <- mean(subset(blankNH4, Batch==1)$N.quantity..ug.) # 3.3648
Batch2.NH4.BlankMass <- mean(subset(blankNH4, Batch==2)$N.quantity..ug.) # 9.2983
Batch3.NH4.BlankMass <- mean(subset(blankNH4, Batch==3)$N.quantity..ug.) # -0.4890
Batch4.NH4.BlankMass <- mean(subset(blankNH4, Batch==4)$N.quantity..ug.) # 5.3567

data6.0 <- data.frame(data6$Batch, data6$Test, data6$SampleID, data6$Reach, data6$ReachID, data6$Bank, data6$Distance, data6$PeakArea.N..Vs., data6$N.quantity..ug., data6$d15N.vs.Air.N2..permil.)
colnames(data6.0) <- c("Batch", "Test", "ID", "Reach", "Site", "Bank", "Distance", "Peak.Area", "N ug", "d15N")
#data6 <- cbind(data6, log_dist)

data6.NH4 <- subset(data6.0, Test=="NH4" & Reach != "BLANK")

d15NH4.Blank <- ifelse(data6.NH4$Batch==1, (data6.NH4$d15N*(Batch1.NH4.BlankMass+data6.NH4$'N ug')-(Batch1.NH4.Blank*Batch1.NH4.BlankMass))/(data6.NH4$'N ug'), 
               ifelse(data6.NH4$Batch==2, (data6.NH4$d15N*(Batch2.NH4.BlankMass+data6.NH4$'N ug')-(Batch2.NH4.Blank*Batch2.NH4.BlankMass))/(data6.NH4$'N ug'), 
                             ifelse(data6.NH4$Batch==3, (data6.NH4$d15N*(Batch3.NH4.BlankMass+data6.NH4$'N ug')-(Batch3.NH4.Blank*Batch3.NH4.BlankMass))/(data6.NH4$'N ug'), 
                                    ifelse(data6.NH4$Batch==4, (data6.NH4$d15N*(Batch4.NH4.BlankMass+data6.NH4$'N ug')-(Batch4.NH4.Blank*Batch4.NH4.BlankMass))/(data6.NH4$'N ug'), 0))))

data6.1NH4 <- cbind(data6.NH4, d15NH4.Blank) 
##################NO3 Data ########

blankNO3 <- subset(data6, Reach=="BLANK" & Test=="NO3")

Batch1.NO3.Blank <- mean(subset(blankNO3, Batch==1)$d15N.vs.Air.N2..permil.) # 3.3648
Batch2.NO3.Blank <- mean(subset(blankNO3, Batch==2)$d15N.vs.Air.N2..permil.) # 9.2983
Batch3.NO3.Blank <- mean(subset(blankNO3, Batch==3)$d15N.vs.Air.N2..permil.) # -0.4890
Batch4.NO3.Blank <- mean(subset(blankNO3, Batch==4)$d15N.vs.Air.N2..permil.) # 5.3567

Batch1.NO3.BlankMass <- mean(subset(blankNO3, Batch==1)$N.quantity..ug.) # 3.3648
Batch2.NO3.BlankMass <- mean(subset(blankNO3, Batch==2)$N.quantity..ug.) # 9.2983
Batch3.NO3.BlankMass <- mean(subset(blankNO3, Batch==3)$N.quantity..ug.) # -0.4890
Batch4.NO3.BlankMass <- mean(subset(blankNO3, Batch==4)$N.quantity..ug.) # 5.3567

data6.0 <- data.frame(data6$Batch, data6$Test, data6$SampleID, data6$Reach, data6$ReachID, data6$Bank, data6$Distance, data6$PeakArea.N..Vs., data6$N.quantity..ug., data6$d15N.vs.Air.N2..permil.)
colnames(data6.0) <- c("Batch", "Test", "ID", "Reach", "Site", "Bank", "Distance", "Peak.Area", "N ug", "d15N")
#data6 <- cbind(data6, log_dist)

data6.NO3 <- subset(data6.0, Test=="NO3" & Reach != "BLANK")

d15NO3.Blank <- ifelse(data6.NO3$Batch==1, (data6.NO3$d15N*(Batch1.NO3.BlankMass+data6.NO3$'N ug')-(Batch1.NO3.Blank*Batch1.NO3.BlankMass))/(data6.NO3$'N ug'), 
                       ifelse(data6.NO3$Batch==2, (data6.NO3$d15N*(Batch2.NO3.BlankMass+data6.NO3$'N ug')-(Batch2.NO3.Blank*Batch2.NO3.BlankMass))/(data6.NO3$'N ug'), 
                              ifelse(data6.NO3$Batch==3, (data6.NO3$d15N*(Batch3.NO3.BlankMass+data6.NO3$'N ug')-(Batch3.NO3.Blank*Batch3.NO3.BlankMass))/(data6.NO3$'N ug'), 
                                     ifelse(data6.NO3$Batch==4, (data6.NO3$d15N*(Batch4.NO3.BlankMass+data6.NO3$'N ug')-(Batch4.NO3.Blank*Batch4.NO3.BlankMass))/(data6.NO3$'N ug'), 0))))


data6.1NO3 <- cbind(data6.NO3, d15NO3.Blank) 

############################################ ############################################ ############################################ 
############################################        Step Correction                      ############################################ 
############################################ ############################################ ############################################
standardsNH4 <- subset(data6.1NH4, Reach=="STD" & Test=="NH4")
plot(standardsNH4$Batch, standardsNH4$d15N.vs.Air.N2..permil.)
abline(a=-0.1385, b=0)

standardsNO3 <- subset(data6.1NO3, Reach=="STD" & Test=="NO3")
plot(standardsNO3$Batch, standardsNO3$d15N.vs.Air.N2..permil.)
abline(a=1.0494, b=0)

Batch1.NH4 <- mean(subset(standardsNH4, Batch==1)$d15NH4.Blank) #-1.042766
Batch2.NH4 <- mean(subset(standardsNH4, Batch==2)$d15NH4.Blank) #-0.7439264
Batch3.NH4 <- mean(subset(standardsNH4, Batch==3)$d15NH4.Blank) #0.4839145
Batch4.NH4 <- mean(subset(standardsNH4, Batch==4)$d15NH4.Blank) #-0.5850174


Batch1.NO3 <- mean(subset(standardsNO3, Batch==1)$d15NO3.Blank) #-14.6874
Batch2.NO3 <- mean(subset(standardsNO3, Batch==2)$d15NO3.Blank) #-9.593663
Batch3.NO3 <- mean(subset(standardsNO3, Batch==3)$d15NO3.Blank) #1.295889
Batch4.NO3 <- mean(subset(standardsNO3, Batch==4)$d15NO3.Blank) #-6.44805



#data6.1 <- subset(data6.0, Bank=="Right" | Bank=="Left")
#data6.NH4 <- subset(data6.1, Test=="NH4")
#data6.NO3 <- subset(data6.1, Test=="NO3")


d15NH4.Corr <- ifelse(data6.1NH4$Batch==1, data6.1NH4$d15NH4.Blank-(Batch1.NH4-NH4.std), 
                      ifelse(data6.1NH4$Batch==2, data6.1NH4$d15NH4.Blank-(Batch2.NH4-NH4.std),
                             ifelse(data6.1NH4$Batch==3, data6.1NH4$d15NH4.Blank - (Batch3.NH4-NH4.std),
                                    ifelse(data6.1NH4$Batch==4, data6.1NH4$d15NH4.Blank - (Batch4.NH4-NH4.std), 0))))

data6.1NH4 <- cbind(data6.1NH4, d15NH4.Corr)  
diff.NH4.d15N <- data6.1NH4$d15N-data6.1NH4$d15NH4.Corr
NH4.CorrectedData <- cbind(data6.1NH4, diff.NH4.d15N)

d15NO3.Corr <- ifelse(data6.1NO3$Batch==1, data6.1NO3$d15NO3.Blank-(Batch1.NO3-NO3.std), 
                      ifelse(data6.1NO3$Batch==2, data6.1NO3$d15NO3.Blank-(Batch2.NO3-NO3.std),
                             ifelse(data6.1NO3$Batch==3, data6.1NO3$d15NO3.Blank - (Batch3.NO3-NO3.std),
                                    ifelse(data6.1NO3$Batch==4, data6.1NO3$d15NO3.Blank - (Batch4.NO3-NO3.std), 0))))

data6.1NO3 <- cbind(data6.1NO3, d15NO3.Corr)  
diff.NO3.d15N <- data6.1NO3$d15N-data6.1NO3$d15NO3.Corr
NO3.CorrectedData <- cbind(data6.1NO3, diff.NO3.d15N)

diff.NO3 <- data.frame(diff.NO3.d15N, na.omit = TRUE)

min(na.omit(diff.NO3.d15N))
max(na.omit(diff.NO3.d15N))

write.csv(NO3.CorrectedData, file = "NO3.StableIsotope.csv")
write.csv(NH4.CorrectedData, file = "NH4.StableIsotope.csv")

data6.NH4 <- merge(data6.1NH4, data3, by = c("Reach", "Site", "Distance", "Bank"))
data6.NO3 <- merge(data6.1NO3, data3, by = c("Reach", "Site", "Distance", "Bank"))
data6.NH4 <- data6.NH4[-c(81),]
#########################################################################################################################
######################################           Model Selection               #####################################
#########################################################################################################################
Function.MixedEffects2 <- Function.ModelSelection <- function(dataframe, variable, extra, extra2) {
  
  aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+extra2+OrgN+(1|Reach/Site), dataframe, REML=TRUE)), #1
                      AIC(lmer(variable~log_dist*Bank+extra2+extra+(1|Reach), dataframe, REML=TRUE)), #2
                      AIC(lmer(variable~Bank*log_dist+extra2+extra+(1|Site), dataframe, REML=TRUE))) #3
  return(aic.output)
                      
}

Function.MixedEffects1 <- Function.ModelSelection <- function(dataframe, variable, extra) {
  
  aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+extra+(1|Reach/Site), dataframe, REML=TRUE)), #1
                      AIC(lmer(variable~log_dist*Bank+extra+(1|Reach), dataframe, REML=TRUE)), #2
                      AIC(lmer(variable~Bank*log_dist+extra+(1|Site), dataframe, REML=TRUE))) #3
  return(aic.output)
  
}

Function.MixedEffects2(data3, data3$NetMin, data3$OrgN, data3$GW) #Reach
Function.MixedEffects2(data3, data3$NetNit, data3$NH4_Conc,data3$GW) #Reach
Function.MixedEffects1(data3, data3$C.N,data3$GW) #Reach
Function.MixedEffects1(data3, data3$OrgN, data3$GW)
Function.MixedEffects1(data3, data3$NO3_Conc, data3$GW) #Reach
Function.MixedEffects1(data3, data3$NH4_Conc, data3$GW) #Reach
Function.MixedEffects1(data3, data3$d15N, data3$massN) #Site
Function.MixedEffects1(data3, data3$d13C, data3$massN) #Reach
Function.MixedEffects1(data6.NO3, data6.NO3$d15NO3.Corr, data6.NO3$GW) #Reach
Function.MixedEffects1(data6.NH4, data6.NH4$d15NH4.Corr, data6.NH4$GW) #Reach


################### two extra covariates
#dataframe <- data3
#variable <- data3$NetMin
#extra <- data3$OrgN
#extra2 <- data3$GW
#mixed <- data3$Reach


Function.ModelSelection.Fixed <- function(dataframe, variable) {
  
  aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+(1|Reach/Site), dataframe, REML=TRUE)), #1
                      AIC(lmer(variable~Bank*log_dist+(1|Reach), dataframe, REML=TRUE)), #2
                      AIC(lmer(variable~Bank*log_dist+(1|Site), dataframe, REML=TRUE)) #3
                        )
  
  names <- seq(1,3,1)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:3,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight)
  colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight")
  return(aic.output)
 
}

NetMin.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$NetMin)) #2
NetNit.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$NetNit)) #2
NO3Conc.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$NO3_Conc)) #2
NH4Conc.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$NH4_Conc)) #3
d15N.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$d15N)) #2
OrgN.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$OrgN)) #2
d13C.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$d13C)) #3
d15N.NO3.Fix <- data.frame(Function.ModelSelection.Fixed(data6.NO3, data6.NO3$d15NO3.Corr)) #2
d15N.NH4.Fix <- data.frame(Function.ModelSelection.Fixed(data6.NH4, data6.NH4$d15NH4.Corr)) #3
PercentN.Fix <- data.frame(Function.ModelSelection.Fixed(data3, data3$pctN)) #2
C.NFix <- data.frame(Function.ModelSelection.Fixed(data3, data3$C.N)) #2

###none of the models have a delAIC that is more than 2 between model 2 and 3 

Function.ModelSelection.Mixed <- function(dataframe, variable) {
  
  aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+(1|Reach), dataframe, REML=FALSE)), #1
                      AIC(lmer(variable~Bank*log_dist+(1|Site), dataframe, REML=FALSE)), #2
                      AIC(lm(variable~Bank*log_dist, dataframe)), #3
                      AIC(lm(variable~Bank*log_dist, dataframe)) #4
  )
  
  names <- seq(1,4,1)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:4,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight)
  colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight")
  return(aic.output)
  
}

NetMin.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$NetMin)) # 3,4
NetNit.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$NetNit)) #3,4
NO3Conc.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$NO3_Conc)) #3,4
NH4Conc.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$NH4_Conc)) #3,4
OrgN.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$OrgN)) #2 BUT NOT SIGNIFICANT
d15N.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$d15N)) #2 BUT NOT SIGNIFICANT
d13C.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$d13C)) #2 BUT NOT SIGNIFICANT
d15N.NO3.Mix <- data.frame(Function.ModelSelection.Mixed(data6.NO3, data6.NO3$d15NO3.Corr)) #1*****
d15N.NH4.Mix <- data.frame(Function.ModelSelection.Mixed(data6.NH4, data6.NH4$d15NH4.Corr)) #3, 4
PercentN.Mix <- data.frame(Function.ModelSelection.Mixed(data3, data3$pctN)) #3, 4
C.NMix <- data.frame(Function.ModelSelection.Mixed(data3, data3$C.N)) #1


########################## 1 extra covariate
  Function.ModelSelection1 <- function(dataframe, variable, extra, mixed) {
    
    aic.output <- rbind(AIC(lmer(variable~Bank*log_dist+(1|Reach), dataframe, REML=FALSE)), #1
                        AIC(lmer(variable~Bank*log_dist+extra+(1|Reach), dataframe, REML=FALSE)), #2
                        AIC(lmer(variable~Bank*log_dist+I(log_dist^2)*Bank+extra+(1|Reach), dataframe, REML=FALSE)), #3
                        AIC(lmer(variable~extra+(1|Reach), dataframe, REML=FALSE)), #4
                        AIC(lmer(variable~log_dist+extra+(1|Reach), dataframe, REML=FALSE)), #5
                        
                        AIC(lm(variable~Bank*log_dist+extra, dataframe)), #6
                        AIC(lm(variable~Bank*log_dist, dataframe)), #7
                        AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank, dataframe)), #8
                        AIC(lm(variable~log_dist+extra, dataframe)), #9
                        AIC(lm(variable~log_dist, dataframe)), #10
                        AIC(lm(variable~Bank+extra, dataframe)), #11
                        AIC(lm(variable~extra, dataframe)), #12
                        AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank+extra, dataframe)) #13
    )
  
   names <- seq(1,13,1)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:13,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight)
 colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight")
   return(aic.output)
  }
  

dataframe<- data3
variable<-data3$NetMin
extra<- data3$OrgN
extra2 <- data3$GW
mixed<- data3$Reach

Function.ModelSelection2 <- function(dataframe, variable, extra, extra2, mixed) {
    
    aic.output <- rbind(AIC(lm(variable~Bank*log_dist+extra2+extra, dataframe)), #1
                        AIC(lm(variable~Bank*log_dist+extra, dataframe)), #2
                        AIC(lm(variable~Bank*log_dist+extra2, dataframe)), #3
                        AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank+extra2+extra, dataframe)), #4
                  
                        AIC(lm(variable~Bank*log_dist, dataframe)), #5
                        AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank, dataframe)), #6
                        AIC(lm(variable~log_dist+extra, dataframe)), #7
                        AIC(lm(variable~log_dist+extra2, dataframe)), #8
                        AIC(lm(variable~Bank+extra, dataframe)), #9
                        AIC(lm(variable~Bank+extra2, dataframe)),#10
                        AIC(lm(variable~extra, dataframe)), #11
                        AIC(lm(variable~extra2, dataframe)), #12
                        AIC(lm(variable~Bank+extra+extra2, dataframe)), #13
                        AIC(lm(variable~extra+extra2, dataframe)), #14
                      
                        AIC(lm(variable~log_dist^2+extra2, dataframe)), #15
                        AIC(lm(variable~Bank, dataframe)), #16
                        #AIC(lm(variable~log_dist^2, dataframe)) #26
                        
                        AIC(lm(variable~Bank+log_dist+extra2+extra, dataframe)), #17
                        AIC(lm(variable~Bank+log_dist+extra, dataframe)), #18
                        AIC(lm(variable~Bank+log_dist+extra2, dataframe)), #19
                        AIC(lm(variable~Bank+log_dist, dataframe)) #20

    )
    n<-20
    names <- seq(1,n,1)
    effect <- c(3, 3, 3, 3, 3, 3, 2, 2, 1, 1, 4, 4, 1, 4, 2, 2, 1, 3, 3, 3)
    delaic <- aic.output-min(aic.output[1:n,])
    aic.weight1 <- exp(-0.5*delaic)
    aic.weight <- aic.weight1/sum(aic.weight1)
    aic.output <- cbind(aic.output, delaic, aic.weight, effect, names)
    colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight", "Effect", "names")
    return(aic.output)
    
  }


#NH4Conc <- data.frame(Function.ModelSelection3(data3, data3$NH4_Conc, data3$GW, data3$Reach))
#dataframe <- data3
#variable <- data3$OrgN
#extra <- data3$GW
#mixed <- data3$Reach
  
  
Function.ModelSelection3 <- function(dataframe, variable, extra, extra2, mixed) {
  
  aic.output <- rbind(
                      AIC(lm(variable~Bank*log_dist+extra, dataframe)), #1
                      AIC(lm(variable~Bank*log_dist, dataframe)),  #2
                      AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank, dataframe)),  #3
                      AIC(lm(variable~log_dist+extra, dataframe)),  #4
                      AIC(lm(variable~Bank+extra, dataframe)),  #5
                      AIC(lm(variable~extra, dataframe)),  #6
                      #AIC(lm(variable~log_dist^2+extra, dataframe)), 
                      AIC(lm(variable~Bank, dataframe)),#7
                      
                      AIC(lm(variable~log_dist, dataframe)), #8
                      #AIC(lm(variable~log_dist^2, dataframe)),
                      AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank+extra, dataframe)), #9
                      
                      
                      AIC(lm(variable~Bank+log_dist, dataframe)), #10
                      AIC(lm(variable~Bank+log_dist+extra, dataframe)) #11

  )
  
  n <- 11
  names <- seq(1,n,1)
  effect <- c(3, 3, 3, 2, 1, 4, 1, 2, 3, 3, 3)
  #row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:n,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight, effect, names)
  colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight", "Effect", "names")
  
  return(aic.output)
  
}

Function.ModelSelection4 <- function(dataframe, variable) {
  
  aic.output <- rbind(
    AIC(lm(variable~Bank*log_dist, dataframe)), #13 #1
    AIC(lm(variable~Bank*log_dist+I(log_dist^2)*Bank, dataframe)), #14 #2
    AIC(lm(variable~Bank, dataframe)),#3
    
    AIC(lm(variable~log_dist, dataframe)), #4
    AIC(lm(variable~log_dist+Bank, dataframe))#5

  )
  
  names <- seq(1,5,1)
  effect <- c(3, 3, 1, 2,3)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:5,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight, effect)
  colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight", "Effect")
  return(aic.output)
  
}

Function.ModelSelection5 <- function(dataframe, variable, extra) {
  
  aic.output <- rbind(
    AIC(lmer(variable~Bank*log_dist+extra+(1|Reach), dataframe)), #1
    AIC(lmer(variable~Bank*log_dist+(1|Reach), dataframe)), #13 #2
    AIC(lmer(variable~Bank*log_dist+I(log_dist^2)*Bank+(1|Reach), dataframe)), #14 #3
    AIC(lmer(variable~log_dist+extra+(1|Reach), dataframe)), #15 #4
    AIC(lmer(variable~Bank+extra+(1|Reach), dataframe)), #17 #5
    AIC(lmer(variable~extra+(1|Reach), dataframe)), #19 #6
 
    AIC(lmer(variable~Bank+(1|Reach), dataframe)),#7
    
    AIC(lmer(variable~log_dist+(1|Reach), dataframe)), #8
    
    AIC(lmer(variable~Bank*log_dist+I(log_dist^2)*Bank+extra+(1|Reach), dataframe)) #9
    #AIC(lmer(variable~log_dist^2+(1|Reach), dataframe)),#10
   # AIC(lmer(variable~log_dist^2+extra+(1|Reach), dataframe)), #23 #7
  )
  
  names <- seq(9,1)
  effect <- c(3, 3, 3, 2, 1, 4, 1, 2, 3)
  row.names(aic.output) <- names
  delaic <- aic.output-min(aic.output[1:9,])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output <- cbind(aic.output, delaic, aic.weight, effect)
  colnames(aic.output)<- c("AIC", "delAIC", "AIC Weight", "Effect")
  
  return(aic.output)
  
}


ModelID1 <- seq(1,9,1)

#dataframe <- NetMin

Function.Models.Supported <-function(dataframe){
  
  dataframe.2 <- subset (dataframe[,c("AIC", "delAIC", "Effect", "names")], dataframe[,'delAIC'] <= 2)
  delaic <- dataframe.2[,'AIC']-min(dataframe.2[,'AIC'])
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  aic.output.SUPP <- cbind(dataframe.2[,'AIC'], delaic, aic.weight, dataframe.2[,'Effect'], dataframe.2[,'names'])
  colnames(aic.output.SUPP)<- c("AIC", "delAIC", "AIC Weight", "Effect", "names")
  return(aic.output.SUPP)
} 

NetMin <- data.frame(Function.ModelSelection2(data3, data3$NetMin, data3$OrgN, data3$GW, data3$Reach))
NetMin.Hyp <- cbind(sum(subset(NetMin, Effect==1)$AIC.Weight), sum(subset(NetMin, Effect==2)$AIC.Weight), sum(subset(NetMin, Effect==3)$AIC.Weight), sum(subset(NetMin, Effect==4)$AIC.Weight))
NetMin.Supp <- data.frame(Function.Models.Supported(NetMin))
NetMin.Hyp.Supp <- cbind(sum(subset(NetMin.Supp, Effect==1)$AIC.Weight), sum(subset(NetMin.Supp, Effect==2)$AIC.Weight), sum(subset(NetMin.Supp, Effect==3)$AIC.Weight), sum(subset(NetMin.Supp, Effect==4)$AIC.Weight))

NetNit <- data.frame(Function.ModelSelection2(data3, data3$NetNit, data3$NH4_Conc, data3$GW, data3$Reach))
NetNit.Hyp <- cbind(sum(subset(NetNit, Effect==1)$AIC.Weight), sum(subset(NetNit, Effect==2)$AIC.Weight), sum(subset(NetNit, Effect==3)$AIC.Weight), sum(subset(NetNit, Effect==4)$AIC.Weight))
NetNit.Supp <- data.frame(Function.Models.Supported(NetNit))
NetNit.Hyp.Supp <- cbind(sum(subset(NetNit.Supp, Effect==1)$AIC.Weight), sum(subset(NetNit.Supp, Effect==2)$AIC.Weight), sum(subset(NetNit.Supp, Effect==3)$AIC.Weight), sum(subset(NetNit.Supp, Effect==4)$AIC.Weight))

summary(lm(NetMin~I(log_dist^2)*Bank, data3))
summary(lm(NetMin~log_dist*Bank+I(log_dist^2)*Bank, data3))
summary(lm(NetMin~Bank+log_dist+I(log_dist^2)*Bank, data3))

NO3Conc <- data.frame(Function.ModelSelection3(data3, data3$NO3_Conc, data3$GW, data3$Reach))
NO3Conc.Hyp <- cbind(sum(subset(NO3Conc, Effect==1)$AIC.Weight), sum(subset(NO3Conc, Effect==2)$AIC.Weight), sum(subset(NO3Conc, Effect==3)$AIC.Weight), sum(subset(NO3Conc, Effect==4)$AIC.Weight))
NO3Conc.Supp <- data.frame(Function.Models.Supported(NO3Conc))
NO3Conc.Hyp.Supp <- cbind(sum(subset(NO3Conc.Supp, Effect==1)$AIC.Weight), sum(subset(NO3Conc.Supp, Effect==2)$AIC.Weight), sum(subset(NO3Conc.Supp, Effect==3)$AIC.Weight), sum(subset(NO3Conc.Supp, Effect==4)$AIC.Weight))


NH4Conc <- data.frame(Function.ModelSelection3(data3, data3$NH4_Conc, data3$GW, data3$Reach))
NH4Conc.Hyp <- cbind(sum(subset(NH4Conc, Effect==1)$AIC.Weight), sum(subset(NH4Conc, Effect==2)$AIC.Weight), sum(subset(NH4Conc, Effect==3)$AIC.Weight), sum(subset(NH4Conc, Effect==4)$AIC.Weight))
NH4Conc.Supp <- data.frame(Function.Models.Supported(NH4Conc))
NH4Conc.Hyp.Supp <- cbind(sum(subset(NH4Conc.Supp, Effect==1)$AIC.Weight), sum(subset(NH4Conc.Supp, Effect==2)$AIC.Weight), sum(subset(NH4Conc.Supp, Effect==3)$AIC.Weight), sum(subset(NH4Conc.Supp, Effect==4)$AIC.Weight))


OrgN. <- data.frame(Function.ModelSelection3(data3, data3$OrgN, data3$GW, data3$Site))
OrgN.Hyp <- cbind(sum(subset(OrgN., Effect==1)$AIC.Weight), sum(subset(OrgN., Effect==2)$AIC.Weight), sum(subset(OrgN., Effect==3)$AIC.Weight), sum(subset(OrgN., Effect==4)$AIC.Weight))
OrgN.Supp <- data.frame(Function.Models.Supported(OrgN.))
OrgN.Hyp.Supp <- cbind(sum(subset(OrgN.Supp, Effect==1)$AIC.Weight), sum(subset(OrgN.Supp, Effect==2)$AIC.Weight), sum(subset(OrgN.Supp, Effect==3)$AIC.Weight), sum(subset(OrgN.Supp, Effect==4)$AIC.Weight))
subset(OrgN., OrgN.[,'delAIC'] <= 2)

d15N <- data.frame(Function.ModelSelection3(data3, data3$d15N, data3$massN, data3$Site))
d15N.Hyp <- cbind(sum(subset(d15N, Effect==1)$AIC.Weight), sum(subset(d15N, Effect==2)$AIC.Weight), sum(subset(d15N, Effect==3)$AIC.Weight), sum(subset(d15N, Effect==4)$AIC.Weight))
d15N.Supp <- data.frame(Function.Models.Supported(d15N))
d15N.Hyp.Supp <- cbind(sum(subset(d15N.Supp, Effect==1)$AIC.Weight), sum(subset(d15N.Supp, Effect==2)$AIC.Weight), sum(subset(d15N.Supp, Effect==3)$AIC.Weight), sum(subset(d15N.Supp, Effect==4)$AIC.Weight))

d13C <- data.frame(Function.ModelSelection3(data3, data3$d13C, data3$massN, data3$Reach))
d13C.Hyp <- cbind(sum(subset(d13C, Effect==1)$AIC.Weight), sum(subset(d13C, Effect==2)$AIC.Weight), sum(subset(d13C, Effect==3)$AIC.Weight), sum(subset(d13C, Effect==4)$AIC.Weight))
d13C.Supp <- data.frame(Function.Models.Supported(d13C))
d13C.Hyp.Supp <- cbind(sum(subset(d13C.Supp, Effect==1)$AIC.Weight), sum(subset(d13C.Supp, Effect==2)$AIC.Weight), sum(subset(d13C.Supp, Effect==3)$AIC.Weight), sum(subset(d13C.Supp, Effect==4)$AIC.Weight))


d15N.NO3 <- data.frame(Function.ModelSelection3(data6.NO3, data6.NO3$d15NO3.Corr,  data6.NO3$GW, data6.NO3$Reach))
d15N.NO3.Hyp <- cbind(sum(subset(d15N.NO3, Effect==1)$AIC.Weight), sum(subset(d15N.NO3, Effect==2)$AIC.Weight), sum(subset(d15N.NO3, Effect==3)$AIC.Weight), sum(subset(d15N.NO3, Effect==4)$AIC.Weight))
d15N.NO3.Supp <- data.frame(Function.Models.Supported(d15N.NO3))
d15N.NO3.Hyp.Supp <- cbind(sum(subset(d13C.Supp, Effect==1)$AIC.Weight), sum(subset(d13C.Supp, Effect==2)$AIC.Weight), sum(subset(d13C.Supp, Effect==3)$AIC.Weight), sum(subset(d13C.Supp, Effect==4)$AIC.Weight))


d15N.NH4 <- data.frame(Function.ModelSelection3(data6.NH4, data6.NH4$d15NH4.Corr,  data6.NH4$GW, data6.NH4$Reach))
d15N.NH4.Hyp <- cbind(sum(subset(d15N.NH4, Effect==1)$AIC.Weight), sum(subset(d15N.NH4, Effect==2)$AIC.Weight), sum(subset(d15N.NH4, Effect==3)$AIC.Weight), sum(subset(d15N.NH4, Effect==4)$AIC.Weight))
d15N.NH4.Supp <- data.frame(Function.Models.Supported(d15N.NH4))
d15N.NH4.Hyp.Supp <- cbind(sum(subset(d15N.NH4.Supp, Effect==1)$AIC.Weight), sum(subset(d15N.NH4.Supp, Effect==2)$AIC.Weight), sum(subset(d15N.NH4.Supp, Effect==3)$AIC.Weight), sum(subset(d15N.NH4.Supp, Effect==4)$AIC.Weight))


PercentN <- data.frame(Function.ModelSelection4(data3, data3$pctN))
PercentN.Hyp <- cbind(sum(subset(PercentN, Effect==1)$AIC.Weight), sum(subset(PercentN, Effect==2)$AIC.Weight), sum(subset(PercentN, Effect==3)$AIC.Weight), sum(subset(PercentN, Effect==4)$AIC.Weight))
PercentN.Supp <- data.frame(Function.Models.Supported(PercentN))
PercentN.Hyp.Supp <- cbind(sum(subset(PercentN.Supp, Effect==1)$AIC.Weight), sum(subset(PercentN.Supp, Effect==2)$AIC.Weight), sum(subset(PercentN.Supp, Effect==3)$AIC.Weight), sum(subset(PercentN.Supp, Effect==4)$AIC.Weight))

C.N <- data.frame(ModelID1, Function.ModelSelection5(data3, data3$C.N, data3$GW))
C.NHyp <- cbind(sum(subset(C.N, Effect==1)$AIC.Weight), sum(subset(C.N, Effect==2)$AIC.Weight), sum(subset(C.N, Effect==3)$AIC.Weight), sum(subset(C.N, Effect==4)$AIC.Weight))
C.NSupp <- data.frame(Function.Models.Supported(C.N))
C.NHyp.Supp <- cbind(sum(subset(C.NSupp, Effect==1)$AIC.Weight), sum(subset(C.NSupp, Effect==2)$AIC.Weight), sum(subset(C.NSupp, Effect==3)$AIC.Weight), sum(subset(C.NSupp, Effect==4)$AIC.Weight))
subset(C.N, C.N[,'delAIC'] <= 2)

GW <- data.frame(Function.ModelSelection4(data3, data3$GW))
GW.Hyp <- cbind(sum(subset(GW, Effect==1)$AIC.Weight), sum(subset(GW, Effect==2)$AIC.Weight), sum(subset(GW, Effect==3)$AIC.Weight), sum(subset(GW, Effect==4)$AIC.Weight))
GW.Supp <- data.frame(Function.Models.Supported(GW))
GW.Hyp.Supp <- cbind(sum(subset(GW.Supp, Effect==1)$AIC.Weight), sum(subset(GW.Supp, Effect==2)$AIC.Weight), sum(subset(GW.Supp, Effect==3)$AIC.Weight), sum(subset(GW.Supp, Effect==4)$AIC.Weight))


d15N.NO3 <- data.frame(Function.ModelSelection2(data6.NO3, data6.NO3$d15NO3.Corr,  data6.NO3$GW, data6.NO3$massN, data6.NO3$Reach))
d15N.NH4 <- data.frame(Function.ModelSelection2(data6.NH4, data6.NH4$d15NH4.Corr,  data6.NH4$GW,data6.NO3$massN, data6.NH4$Reach))

Rownames <- c("d15N", "d13C", "d15NH4", "d15NO3", "Percent Nitrogen", "NO3 Concentration", "NH4 Concentration", "Net Nitrification", "Net Mineralization", "Organic Nitrgen")
Hyp <- data.frame(rbind(d15N.Hyp, d13C.Hyp, d15N.NH4.Hyp, d15N.NO3.Hyp, PercentN.Hyp, NO3Conc.Hyp, NH4Conc.Hyp, NetNit.Hyp, NetMin.Hyp, OrgN.Hyp))
Hyp <- cbind(Rownames, Hyp)
colnames(Hyp) <- c("Test","Hyp 1 Bank", "Hyp 2 Distance", "Hyp 3 Bank & Distance", "Hyp 4 Other Variables")
write.csv(Hyp, file="HypothesisTable.csv")

Rownames <- c("d15N", "d13C", "d15NH4", "d15NO3", "Percent Nitrogen", "NO3 Concentration", "NH4 Concentration", "Net Nitrification", "Net Mineralization", "Organic Nitrgen", "GW", "C:N")
Hyp.Supp <- data.frame(rbind(d15N.Hyp.Supp, d13C.Hyp.Supp, d15N.NH4.Hyp.Supp, d15N.NO3.Hyp.Supp, PercentN.Hyp.Supp, NO3Conc.Hyp.Supp, NH4Conc.Hyp.Supp, NetNit.Hyp.Supp, NetMin.Hyp.Supp, OrgN.Hyp.Supp, GW.Hyp.Supp, C.NHyp.Supp))
Hyp.Supp <- cbind(Rownames, Hyp.Supp)
colnames(Hyp.Supp) <- c("Test","Hyp 1 Bank", "Hyp 2 Distance", "Hyp 3 Bank & Distance", "Hyp 4 Other Variables")
write.csv(Hyp.Supp, file="HypothesisTableSupported.csv")




######################################################################################################################
######################################        Supported models       #####################################
#####################################################################################################################

model.d15N <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data3)
model.d13C <- lm(d13C~Bank*log_dist+I(log_dist^2)*Bank, data3)
model.Min <- lm(NetMin~Bank*log_dist+OrgN, data3)
model.Nit <- lm(NetNit~Bank*log_dist+GW+NH4_Conc, data3)
model.NH4 <- lm(NH4_Conc~Bank*log_dist, data3)
model.NO3 <- lm(NO3_Conc~Bank*log_dist+I(log_dist^2)*Bank, data3)
model.d15N.NO3 <- lm(d15N~Bank*log_dist+I(log_dist^2)*Bank, data6.NO3)
model.d15N.NH4 <- lm(d15N.x~Bank*log_dist+I(log_dist^2)*Bank, data6.NH4)
#model.massN <- lm(massN~Bank*log_dist+GW, massN)
model.OrgN <- lm(OrgN~Bank*log_dist+I(log_dist^2)*Bank+GW, data3)
model.GW <- lm(GW~log_dist, data3)
model.C.N <- lm(C.N~Bank*log_dist+I(log_dist^2)*Bank, data3)
lm(C.N~Bank, data3)

summary(model.NH4)
summary(model.NO3)
summary(model.d15N)
summary(model.d13C)



Left <- subset(data3,  Bank == 'Left')
Right <- subset(data3, Bank == 'Right')

mean(Left$NO3_Conc)
mean(Right$NO3_Conc)

subset(NetMin, delAIC<2)
subset(NetNit, delAIC<2)
subset(d15N, delAIC<2)
subset(d13C, delAIC<2)
subset(NO3Conc, delAIC<2)
subset(NH4Conc, delAIC<2)
subset(d15N.NH4, delAIC<2)
subset(OrgN., delAIC<2)
subset(C.N, delAIC<2)
subset(GW, delAIC<2)

######################################################################################################################
######################################        Salmon contribution to soil        #####################################
#####################################################################################################################
subset(data6.NO3, d15NO3.Corr>12.62)
subset(data6.NH4, d15NH4.Corr>12.62)
length(subset(data6.NH4, d15NH4.Corr>12.62)$d15NH4.Corr) #15
length(subset(data6.NO3, d15NO3.Corr>12.62)$d15NO3.Corr) #16
length(subset(data3, d15N>12)$d15N)#0

plot(data3$d15N)


Function.MixingModel <- Function.MixingModel <- function(dataframe) {
  
  mixing.output <- ((dataframe["SAM",]-dataframe["TEM",])/(dataframe["MEM",]-dataframe["TEM",]))*100
  
  return(mixing.output)
  
}

max.NO3 <- max(data6.NO3$d15NO3.Corr) 
max.NH4 <- max(data6.NH4$d15NH4.Corr, na.rm = TRUE) 
mean.NO3 <- mean(data6.NO3$d15NO3.Corr) 
mean.NH4 <- 15
#<- mean(subset(data6.NH4, Distance==3 & Bank=='Left')$d15NH4.Corr,na.rm = TRUE) 
min.NH4 <- min(subset(data6.NH4, Bank == 'Left')$d15NH4.Corr, na.rm = TRUE) 
max(subset(data6.NH4, Bank == 'Left')$d15NH4.Corr, na.rm = TRUE) 

names <- c("SAM", "TEM", "MEM")
names2 <- c("Quinn et al. Enhanced", "Quinn et al. Depleted", "Observed Max NH4", "Enhanced Max NH4 End Member", "Depleted Max NH4 End Member",
            "Observed Mean NH4","Enhanced Mean NH4 End Member", "Depleted Mean NH4 End Member")
  
  
quinn.ENH1 <- data.frame(c(10.7, -1.74, 12.62), row.names = names)
quinn.DEP1 <- data.frame(c(6.4, -1.74, 12.62), row.names = names)

DEP.MAX <- data.frame(c(6.4, -1.74, max.NH4), row.names = names)
ENH.MAX <- data.frame(c(10.7, -1.74, max.NH4), row.names = names)
NH4.max <- data.frame(c(max.NH4, -1.74, 12.62), row.names = names)

DEP.MEAN <- data.frame(c(6.4, -1.74, mean.NH4), row.names = names)
ENH.MEAN <- data.frame(c(10.7, -1.74, mean.NH4), row.names = names)
NH4.mean <- data.frame(c(mean.NH4, -1.74, 12.62), row.names = names)

Percent.MDN <- data.frame(c(Function.MixingModel(quinn.ENH1),
Function.MixingModel(quinn.DEP1),

Function.MixingModel(NH4.max),
Function.MixingModel(ENH.MAX),
Function.MixingModel(DEP.MAX),

Function.MixingModel(NH4.mean),
Function.MixingModel(ENH.MEAN),
Function.MixingModel(DEP.MEAN)), row.names = names2)

colnames(Percent.MDN) <- c("% MDN")

Percent.MDN

Function.MixingModel2 <- Function.MixingModel <- function(dataframe) {
  
  mixing.output <- ((dataframe["SAM",]-dataframe["MEM",])/(dataframe["TEM",]-dataframe["MEM",]))*100
  
  return(mixing.output)
  
}

Percent.TDN <- data.frame(c(Function.MixingModel2(quinn.ENH1),
                            Function.MixingModel2(quinn.DEP1),
                            
                            Function.MixingModel2(NH4.max),
                            Function.MixingModel2(ENH.MAX),
                            Function.MixingModel2(DEP.MAX),
                            
                            Function.MixingModel2(NH4.mean),
                            Function.MixingModel2(ENH.MEAN),
                            Function.MixingModel2(DEP.MEAN)), row.names = names2)

colnames(Percent.TDN) <- c("% MDN")

write.csv(Percent.TDN, file = "Percent.TDN.csv")
write.csv(Percent.MDN, file="Percent.MDN.csv")          


max(Left$d15N)
max(Right$d15N)
Left

par(mar=c(5,5,0.25,0.25), oma=rep(0,4))
plot(data6.NH4$d15N.y, data6.NH4$d15NH4.Corr, pch=16,  ylab= expression(paste(delta^15, "N ", NH[4]^" +")), xlab=expression(paste(delta^15, "N", sep="")))
abline(a=1.6822, b=1.1010, lty=2, col='red')
summary(lm(data6.NH4$d15NH4.Corr~data6.NH4$d15N.y))
text( 10,40, expression(paste(R^2," = 0.08")))
text( 10,37, "p < 0.05")


plot(data6.NH4$pctN, data6.NH4$d15NH4.Corr, pch=16,  ylab= expression(paste(delta^15, "N ", NH[4]^" +")), xlab=expression(paste(delta^15, "N", sep="")))
summary(lm(data6.NH4$d15NH4.Corr~data6.NH4$massN))

nrow(subset(data6.NH4, Bank == 'Right'))#44
nrow(subset(data6.NH4, Bank == 'Left'))#39

ex.12 <- subset(data6.NH4, d15NH4.Corr>12.62)

ex.12.R <- subset(ex.12, Bank == 'Right')
nrow(ex.12.R) #4
4/44 #9%
ex.12.L <- subset(ex.12, Bank == 'Left')
nrow(ex.12.L) #10
10/39 #26%

sum(na.omit(subset(data6.NH4, Bank =='Left')$d15NH4.Corr - subset(data6.NH4, Bank =='Left')$d15N.x) > 0) #35
sum(na.omit(subset(data6.NH4, Bank =='Left')$d15NH4.Corr - subset(data6.NH4, Bank =='Left')$d15N.x) < 0) # 2
35/37 #95%

sum(na.omit(subset(data6.NH4, Bank =='Right')$d15NH4.Corr - subset(data6.NH4, Bank =='Right')$d15N.x) > 0) #37
sum(na.omit(subset(data6.NH4, Bank =='Right')$d15NH4.Corr - subset(data6.NH4, Bank =='Right')$d15N.x) < 0) # 7
37/44 #84%

