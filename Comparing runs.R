# COmparing sulforaphane ESI- plasma runs to see which files I should use.
# 1/30/13

setwd("D:/My Documents/Work/Lin Lab/Brocco rif cortisol project/Sulforaphane metabolomics/ESI- samples/Plasma")
Data <- read.csv("Comparing runs.csv", header=T)

Run1114 <- Data[Data$Run=="20111114",]
Run1115 <- Data[Data$Run=="20111115",]
Run1130 <- Data[Data$Run=="20111130",]
Run1201 <- Data[Data$Run=="20111201",]
Run1205 <- Data[Data$Run=="20111205",]

Colors <- c("black", "gray", "red", "darkgreen", "blue")
Runs <- c("1114", "1115", "1130", "1201", "1205")

###########################
# RTs
par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,7], Run1115[,7], Run1130[,7], Run1201[,7], Run1205[,7],
        col=Colors, names=Runs, main="RT of d4-salicylic acid",
        xlab="Run", ylab="RT (min)")


par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,9], Run1115[,9], Run1130[,9], Run1201[,9], Run1205[,9],
        col=Colors, names=Runs, main="RT of prednisolone-21-sulfate",
        xlab="Run", ylab="RT (min)")


par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,11], Run1115[,11], Run1130[,11], Run1201[,11], Run1205[,11],
        col=Colors, names=Runs, main="RT of d7-stearic acid",
        xlab="Run", ylab="RT (min)")



###############################
# Abundances
par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,8], Run1115[,8], Run1130[,8], Run1201[,8], Run1205[,8],
        col=Colors, names=Runs, main="Abundance of d4-salicylic acid",
        xlab="Run", ylab="Peak Area")


par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,10], Run1115[,10], Run1130[,10], Run1201[,10], Run1205[,10],
        col=Colors, names=Runs, main="Abundance of prednisolone-21-sulfate",
        xlab="Run", ylab="Peak Area")


par(mfrow=c(1,1), mar=c(3,3,2.5,0.5), mgp=c(1.5,0.5,0))
boxplot(Run1114[,12], Run1115[,12], Run1130[,12], Run1201[,12], Run1205[,12],
        col=Colors, names=Runs, main="Abundance of d7-stearic acid",
        xlab="Run", ylab="Peak Area")
