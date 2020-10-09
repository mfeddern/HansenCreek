data1 <- Holtgrieve_HansenSoil_Samplelist_base

RL_I<-subset(data1, Side ==  "RL" & Incubation == "I")
RL_F <- subset(data1, Side ==  "RL" & Incubation == "F")
RR_I<-subset(data1, Side ==  "RR" & Incubation == "I")
RR_F <- subset(data1, Side ==  "RR" & Incubation == "F")
RR <- subset(data1, Side ==  "RR")
RL <- subset(data1, Side ==  "RL")
Initial <- rbind(RL_I, RR_I)
Final <- rbind(RL_F, RR_F)
NH4_RR_Change <- RR_F$`NH4-N mg/L`- RR_I$`NH4-N mg/L`
NH4_RL_Change <- RL_F$`NH4-N mg/L`- RL_I$`NH4-N mg/L`
NO3_RR_Change <- RR_F$`NO3-N mg/L`- RR_I$`NO3-N mg/L`
NO3_RL_Change <- RL_F$`NO3-N mg/L`- RL_I$`NO3-N mg/L`


WaterContentRR <- RR$`wet weight without tin`-RR$`dry weight without tin`
WaterContentRL <- RL$`wet weight without tin`-RL$`dry weight without tin`

#############################################################################################################
############################## paired t-tests comparison ####################################################
#############################################################################################################


####Initial NH4 side comparison
t.test(RL_I$`NH4-N mg/L`, RR_I$`NH4-N mg/L`, paired=TRUE, conf.level=0.95)
mean(RL_I$`NH4-N mg/L`) #1.288222
sd(RL_I$`NH4-N mg/L`) #2.51702
mean(RR_I$`NH4-N mg/L`) #0.4157778
sd(RR_I$`NH4-N mg/L`) #0.2884164
#t = 2.2866, df = 44, p-value = 0.02709

###Initial NO3 side comparison
t.test(RL_I$`NO3-N mg/L`, RR_I$`NO3-N mg/L`, paired=TRUE, conf.level=0.95)
#t = 0.2627, df = 44, p-value = 0.794

###Final NH4 side comparison
t.test(RL_F$`NO3-N mg/L`, RR_F$`NO3-N mg/L`, paired=TRUE, conf.level=0.95)
#t = 0.55047, df = 44, p-value = 0.5848

####Final NO3 side comparison
t.test(RL_F$`NH4-N mg/L`, RR_F$`NH4-N mg/L`, paired=TRUE, conf.level=0.95)
#t = 1.4568, df = 44, p-value = 0.1523





###RL NH4 change
t.test(RL_I$`NH4-N mg/L`, RL_F$`NH4-N mg/L`, paired=TRUE, conf.level=0.95)
#t = -1.4068, df = 44, p-value = 0.1665

###RR NH4 change
t.test(RR_I$`NH4-N mg/L`, RR_F$`NH4-N mg/L`, paired=TRUE, conf.level=0.95)
#t = -2.7824, df = 44, p-value = 0.007919
mean(RR_I$`NH4-N mg/L`) #0.4158
mean(RR_F$`NH4-N mg/L`) #0.794
sd(RR_I$`NH4-N mg/L`) #0.288
sd(RR_F$`NH4-N mg/L`) #0.987

###RL NO3 change
t.test(RL_I$`NO3-N mg/L`, RL_F$`NO3-N mg/L`, paired=TRUE, conf.level=0.95)
#t = -5.9504, df = 44, p-value = 3.991e-07
mean(RL_I$`NO3-N mg/L`) #0.2224
mean(RL_F$`NO3-N mg/L`) #1.148
sd(RL_I$`NO3-N mg/L`) #0.1743
sd(RL_F$`NO3-N mg/L`) #1.1339



###RR N03 change
t.test(RR_I$`NO3-N mg/L`, RR_F$`NO3-N mg/L`, paired=TRUE, conf.level=0.95)
#t = -4.3688, df = 44, p-value = 7.501e-05
mean(RR_I$`NO3-N mg/L`) #0.2046
mean(RR_F$`NO3-N mg/L`) #0.9815
sd(RR_I$`NO3-N mg/L`) #0.4157
sd(RR_F$`NO3-N mg/L`) #1.5357





##Comparison of water content on either side
t.test(WaterContentRL, WaterContentRR, paired=TRUE, conf.level=0.95)
#t = -1.7335, df = 44, p-value = 0.09001

##Comparing NH4 Change
t.test(NH4_RR_Change, NH4_RL_Change, paired=TRUE, conf.level=0.95)
#t = 0.63368, df = 44, p-value = 0.5296

##Comparing NO3 Change
t.test(NO3_RR_Change, NO3_RL_Change, paired=TRUE, conf.level=0.95)
#t = -0.59333, df = 44, p-value = 0.556




#############################################################################################################
####################### Plotting and including distances ####################################################
#############################################################################################################

plot(RR_I$`River Distance`, RR_I$`NH4-N mg/L`, pch=16, ylab = "NH4 mg/L", xlab="Distance", main = "Initial")
points(RL_I$`River Distance`, RL_I$`NH4-N mg/L`, pch=16, col="red")
legend("topleft", legend=c("RR (no salmon)", "RL (salmon)"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)


plot(RR_F$`River Distance`, RR_F$`NH4-N mg/L`, pch=16, ylab = "NH4 mg/L", xlab="Distance", main = "Final")
points(RL_F$`River Distance`, RL_F$`NH4-N mg/L`, pch=16, col="red")
legend("topleft", legend=c("RR (no salmon)", "RL (salmon)"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)

plot(RR_I$`River Distance`, RR_I$`NO3-N mg/L`, pch=16, ylab = "NO3 mg/L", xlab="Distance", main = "Initial")
points(RL_I$`River Distance`, RL_I$`NO3-N mg/L`, pch=16, col="red")
legend("topleft", legend=c("RR (no salmon)", "RL (salmon)"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)


plot(RR_F$`River Distance`, RR_F$`NO3-N mg/L`, pch=16, ylab = "NO3 mg/L", xlab="Distance", main = "Final")
points(RL_F$`River Distance`, RL_F$`NO3-N mg/L`, pch=16, col="red")
legend("topleft", legend=c("RR (no salmon)", "RL (salmon)"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)

plot(RR_I$`River Distance`, RR_I$`NH4-N mg/L`, pch=16, ylab = "NH4 mg/L", xlab="Distance", main = "RR (no salmon)")
points(RR_F$`River Distance`, RR_F$`NH4-N mg/L`, pch=16, col="red")
legend("topleft", legend=c("RR Initial", "RR Final"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)

plot(RL_I$`River Distance`, RL_I$`NH4-N mg/L`, pch=16, ylab = "NH4 mg/L", xlab="Distance", main = "RL (salmon)")
points(RL_F$`River Distance`, RL_F$`NH4-N mg/L`, pch=16, col="red")
legend("topright", legend=c("RL Initial", "RL Final"), col=c("black", "red"),lty=1:1, cex=0.8, box.lty=0)

#############################################################################################################
#################################### box plot comparison ####################################################
#############################################################################################################

boxplot(`NH4-N mg/L`~Incubation, data=RR)
boxplot(`NH4-N mg/L`~Incubation, data=RL)
boxplot(`NO3-N mg/L`~Incubation, data=RR)
boxplot(`NO3-N mg/L`~Incubation, data=RL)

boxplot(`NO3-N mg/L`~`River Distance`, data=RR_F)
boxplot(`NO3-N mg/L`~`River Distance`, data=RR_I)

boxplot(`NH4-N mg/L`~`River Distance`, data=RR_F)
boxplot(`NH4-N mg/L`~`River Distance`, data=RR_I)

boxplot(`NO3-N mg/L`~`River Distance`, data=RL_F)
boxplot(`NO3-N mg/L`~`River Distance`, data=RL_I)

boxplot(`NH4-N mg/L`~`River Distance`, data=RL_F)
boxplot(`NH4-N mg/L`~`River Distance`, data=RL_I)


#############################################################################################################
################################### anova comparison ####################################################
#############################################################################################################
model1.1 <- lm(`NH4-N mg/L`~Incubation+Side+`River Distance`+`Transect Number`, data=data)
summary(model1.1)
anova(model1.1)

model1.2 <- lm(`NH4-N mg/L`~`River Distance`+Side, data=data)
summary(model1.2)
anova(model1.2)


model2.1 <- lm(`NO3-N mg/L`~Incubation+Side+`River Distance`, data=data)
summary(model2.1)
anova(model2.1)

model2.2 <- lm(`NO3-N mg/L`~Incubation+`River Distance`, data=data)
summary(model2.2)
anova(model2.2)



