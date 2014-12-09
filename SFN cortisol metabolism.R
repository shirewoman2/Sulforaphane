# SFN cortisol metabolism

setwd("D:/My Documents/Work/Lin Lab/Brocco rif cortisol project")
SFN <- read.csv("All SFN cortisol data.csv", header=T, skip=1)

SFN$FRatio <- SFN$SixBOHF/SFN$Cortisol
SFN$ERatio <- SFN$SixBOHE/SFN$Cortisone
SFN$TotalRatio <- (SFN$SixBOHF+SFN$SixBOHE)/(SFN$Cortisol+SFN$Cortisone)

Colors <- rep(NA, nrow(SFN))
Colors[which(SFN$Effector == "rif")] <- "red"
Colors[which(SFN$Effector == "brocco")] <- "limegreen"
Colors[which(SFN$Effector == "brocco/rif")] <- "blue"

Symbols <- rep(NA, nrow(SFN))
Symbols[which(SFN$ClinDay==1)] <- 1
Symbols[which(SFN$ClinDay==8)] <- 16

SFN$Induced <- rep("baseline", nrow(SFN))

SFN$Induced[which(SFN$ClinDay==8 & 
                    (SFN$Effector=="rif" | SFN$Effector=="brocco/rif"))] <- "induced"

Induced <- SFN[which(SFN$Induced=="induced"),]
Baseline <- SFN[which(SFN$Induced=="baseline"),]



# F metabolic ratio
plot(SFN$FRatio, SFN$MDZCL, col=Colors, pch=Symbols,
     xlab="Cortisol Metabolic Ratio", ylab="MDZ CL/F",
     main="Cortisol")

ttest.F <- t.test(Baseline$FRatio, Induced$FRatio, paired=F, var.equal=F) # I know that these are actually paired samples, but I haven't set up the pairs yet. Just using this to get an idea of the differences. 

boxplot(Baseline$FRatio, Induced$FRatio, names=c("Baseline", "Induced"),
        ylab="Cortisol Metabolic Ratio", main="Cortisol", 
        sub=paste("p =", signif(ttest.F$p.value, digits=2)))



# E metabolic ratio
plot(SFN$ERatio, SFN$MDZCL, col=Colors, pch=Symbols,
     xlab="Cortisone Metabolic Ratio", ylab="MDZ CL/F",
     main="Cortisone")

ttest.E <- t.test(Baseline$ERatio, Induced$ERatio, paired=F, var.equal=F) # I know that these are actually paired samples, but I haven't set up the pairs yet. Just using this to get an idea of the differences. 


boxplot(Baseline$ERatio, Induced$ERatio, names=c("Baseline", "Induced"),
        ylab="Cortisone Metabolic Ratio", main="Cortisone",
        sub=paste("p =", signif(ttest.E$p.value, digits=2)))



# Total metabolic ratio
plot(SFN$TotalRatio, SFN$MDZCL, col=Colors, pch=Symbols,
     xlab="Cortisol + Cortisone Metabolic Ratio", ylab="MDZ CL/F",
     main="Cortisone + Cortisol")

ttest.tot <- t.test(Baseline$TotalRatio, Induced$TotalRatio, paired=F, var.equal=F) # I know that these are actually paired samples, but I haven't set up the pairs yet. Just using this to get an idea of the differences. 


boxplot(Baseline$ERatio, Induced$ERatio, names=c("Baseline", "Induced"),
        ylab="Total Metabolic Ratio", main="Cortisone + Cortisol",
        sub=paste("p =", signif(ttest.tot$p.value, digits=2)))


