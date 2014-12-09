# SFN checking cortisol values

library(reshape2)
library(plyr)

setwd("D:/Users/Laura/Documents/Work/Lin Lab/Sulforaphane project")
SFN <- read.csv("SFN checking cortisol values.csv", skip=1, 
                colClasses=c(rep("factor", 5), "numeric", "factor", 
                             "factor", rep("numeric", 13)))
# CV <- read.csv("SFN CV.csv", colClasses=c("factor", rep("numeric", 4)))
# Something's screwy with the values for CV. Use the recalculated values below: CV.recalc

SampList <- melt(SFN[,c(3,7:13)], measure=c("C6bOHF", "C6bOHE", "CF", "CE"))
names(SampList)[names(SampList)=="variable"] <- "Analyte1"
names(SampList)[names(SampList)=="value"] <- "Concentration"

Analyte2 <- colsplit(SampList$Analyte1, "C", c("X", "Analyte"))
SampList$Analyte2 <- paste("X",Analyte2$Analyte, sep="")

SampList$Analyte1 <- NULL
names(SampList)[names(SampList)=="Analyte2"] <- "Analyte"
SampList$Analyte <- as.factor(SampList$Analyte)

Max <- melt(SFN[,c(3,7:9,14:17)], measure=c("Max6bOHF", "Max6bOHE", "MaxF", "MaxE"))
names(Max)[names(Max)=="variable"] <- "Analyte1"
names(Max)[names(Max)=="value"] <- "Max"

Analyte2 <- colsplit(Max$Analyte1, "Max", c("X", "Analyte"))
Max$Analyte2 <- paste("X",Analyte2$Analyte, sep="")

Max$Analyte1 <- NULL
names(Max)[names(Max)=="Analyte2"] <- "Analyte"
Max$Analyte <- as.factor(Max$Analyte)

SampList <- merge(SampList, Max, by=c("Analyte", "SubjSamp", "RunDate", "Scientist", "Vol"))

Min <- melt(SFN[,c(3,7:9,18:21)], measure=c("Min6bOHF", "Min6bOHE", "MinF", "MinE"))
names(Min)[names(Min)=="variable"] <- "Analyte1"
names(Min)[names(Min)=="value"] <- "Min"

Analyte2 <- colsplit(Min$Analyte1, "Min", c("X", "Analyte"))
Min$Analyte2 <- paste("X",Analyte2$Analyte, sep="")

Min$Analyte1 <- NULL
names(Min)[names(Min)=="Analyte2"] <- "Analyte"
Min$Analyte <- as.factor(Min$Analyte)

SampList <- merge(SampList, Min, by=c("Analyte", "SubjSamp", "RunDate", "Scientist", "Vol"))

SampList <- SampList[complete.cases(SampList$Concentration),]
SampList$Use <- rep("unknown", nrow(SampList))

# CV.melt <- melt(CV, id="SubjSamp")
# names(CV.melt)[names(CV.melt)=="variable"] <- "Analyte1"
# names(CV.melt)[names(CV.melt)=="value"] <- "CV"
# Analyte2 <- colsplit(CV.melt$Analyte1, "CV", c("X", "Analyte"))
# CV.melt$Analyte2 <- paste("X",Analyte2$Analyte, sep="")
# CV.melt <- CV.melt[,c("SubjSamp", "Analyte2", "CV")]
# names(CV.melt)[names(CV.melt)=="Analyte2"] <- "Analyte"
# SampList <- merge(SampList, CV.melt, by=c("SubjSamp", "Analyte"))

CV <- ddply(SampList, 
            c("SubjSamp", "Analyte"), summarise,
            Mean=mean(Concentration, na.rm=T),
            SD=sd(Concentration, na.rm=T),
            CV.recalc=SD/Mean)

CV.recalc <- melt(CV[,c("SubjSamp", "Analyte", "CV.recalc")], id=c("SubjSamp", "Analyte"))
names(CV.recalc)[names(CV.recalc)=="value"] <- "CV.recalc"
CV.recalc$variable <- NULL
  
Compare <-  merge(CV.melt, CV.recalc, by=c("SubjSamp", "Analyte"))

SampList <- merge(SampList, CV.recalc, by=c("SubjSamp", "Analyte"))
SampList$Use[SampList$CV.recalc < 0.15] <- "yes"
SampList$Use[is.na(SampList$CV.recalc)] <- "yes"

SampList <- SampList[order(SampList$Analyte, SampList$SubjSamp, SampList$RunDate),]
row.names(SampList) <- 1:nrow(SampList)

SampList.orig <- SampList

CheckedRows <- list()
CheckedRows[["X6bOHF"]] <- c()
CheckedRows[["X6bOHF"]] <- c()
CheckedRows[["XF"]] <- c()
CheckedRows[["XE"]] <- c()


# Already checked these samples and know I want to keep all values.
SubjSamp.keep <- c("2371 180", "2378 210", "2434 110", "2434 310",  # 6bOHE
                   "2371 310", "2438 310") # F
Analyte.keep <- c(rep("X6bOHE", 4), rep("XF", 2))
Keep <- data.frame("SubjSamp"=SubjSamp.keep, "Analyte"=Analyte.keep)

for (i in 1:nrow(Keep)){
  SampList$Use[SampList$Analyte==as.character(Keep$Analyte[i]) & 
                 SampList$SubjSamp==as.character(Keep$SubjSamp[i])] <- "yes"  
}


SampList$Use[SampList$Max > 1.1] <- "no"
SampList$Use[SampList$Min < 0.9] <- "no"


# 6b-OH-F check
Check <- SampList[(SampList$Analyte=="X6bOHF" & SampList$Use=="unknown"),]
CheckedRows[["X6bOHF"]] <- data.frame("Rows"=row.names(Check))

S <- as.character(unique(Check$SubjSamp))
length(S)
# 6

YesRows <- c(264,265,     254,     439,
             271,272,     298,     
             303,304,     254,
             306,307,     299,
             343,344,     412,
             440,442)
NoRows <- c(263,     255,
            270,
            302,
            305,
            342,
            441)
Check[Check$SubjSamp==S[4],]

SampList$Use[YesRows] <- "yes"
SampList$Use[NoRows] <- "no"

CheckedRows[["X6bOHF"]]$Use <- SampList$Use[as.numeric(row.names(Check))]


SampList.final <- SampList


# 6b-OH-E check
Check <- SampList[(SampList$Analyte=="X6bOHE" & SampList$Use=="unknown"),]
CheckedRows[["X6bOHE"]] <- data.frame("Rows"=row.names(Check))

S <- as.character(unique(Check$SubjSamp))
length(S)
# 21
# 4

YesRows <- c(12,       68,69,     167,169,
             15,16,    71,72,     185,
             21,22,    74,75,     207,208,
             24,25,    87,88,     213,214,
             32,33,    97,98,     221,222,
             50,51,    109,110,     
             59,       150,151,
             63,       153,154
             )
NoRows <- c(11,13,     67,     168,
            14,        70,     184,
            20,        73,     209,
            23,        86,     215,
            31,        96,     220,
            49,        108,
            58,60,     152,
            62,64,     155
            )

Check[Check$SubjSamp==S[4],]

YesRows <- c(29,30, 104,187, 246)
NoRows <- c(247)

SampList$Use[YesRows] <- "yes"
SampList$Use[NoRows] <- "no"

CheckedRows[["X6bOHE"]]$Use <- SampList$Use[as.numeric(row.names(Check))]


SampList.final <- SampList




# F check
Check <- SampList[(SampList$Analyte=="XF" & SampList$Use=="unknown"),]
CheckedRows[["XF"]] <- data.frame("Rows"=row.names(Check))

S <- as.character(unique(Check$SubjSamp))
length(S)
# 11
# 7

YesRows <- c(734,          772,773,
             737,738,          818,820,
             746,747,          831,832,
             751,752,          875,876,
             754,755,          889,
             949,950)
NoRows <- c(733,735,          771,
            736,          819,
            745,          830,
            750,          874,
            753,          888,890,
            948)

Check[Check$SubjSamp==S[7],]

YesRows <- c(742, 760,762,768,769,790,816, 919)

SampList$Use[YesRows] <- "yes"
SampList$Use[NoRows] <- "no"

CheckedRows[["XF"]]$Use <- SampList$Use[as.numeric(row.names(Check))]


SampList.final <- SampList



# E check
Check <- SampList[(SampList$Analyte=="XE" & SampList$Use=="unknown"),]
CheckedRows[["XE"]] <- data.frame("Rows"=row.names(Check))

S <- as.character(unique(Check$SubjSamp))
length(S)
# 11
# 3

YesRows <- c(487,         525,526,
             496,497,         533,534,
             499,500,         543,544,
             504,505,         584,585,
             507,508,         642,
             521,522
             )
NoRows <- c(486,488,         524,
            495,         535,
            498,         542,
            503,         583,
            506,         641,643,
            520
            )

Check[Check$SubjSamp==S[3],]

YesRows <- c(609,629,660)
NoRows <- c(661)

SampList$Use[YesRows] <- "yes"
SampList$Use[NoRows] <- "no"

CheckedRows[["XE"]]$Use <- SampList$Use[as.numeric(row.names(Check))]


SampList.final <- SampList


#############
# Recalculating CV
SFN.recalc <- ddply(SampList[SampList$Use=="yes",], 
                    c("SubjSamp", "Analyte"), summarise,
                    Mean=mean(Concentration, na.rm=T),
                    SD=sd(Concentration, na.rm=T),
                    CV=SD/Mean)
SFN.recalc <- merge(SFN.recalc, CV.melt, by=c("SubjSamp", "Analyte"))
names(SFN.recalc)[names(SFN.recalc)=="CV.x"] <- "CV.recalc"
names(SFN.recalc)[names(SFN.recalc)=="CV.y"] <- "CV.orig"
SFN.recalc$HigherCV <- rep("no", nrow(SFN.recalc))
SFN.recalc$HigherCV[SFN.recalc$CV.recalc > 1.1*(SFN.recalc$CV.orig)] <- "yes"

Higher <- SFN.recalc[SFN.recalc$HigherCV=="yes",]
Higher <- Higher[order(Higher$Analyte),]

# SampList[c(SampList$Analyte==Higher$Analyte[3] &
#              SampList$SubjSamp==Higher$SubjSamp[3]),]
# 
# SampList$Use[905] <- "no"


###################
Meta <- unique(SFN[,c("Subject", "Sample", "SubjSamp", "Tx", "ClinDay", "MDZCL")])
Meta <- Meta[order(Meta$Subject, Meta$Tx, Meta$ClinDay),]
row.names(Meta) <- 1:nrow(Meta)

SFN.final <- dcast(SFN.recalc[,c("SubjSamp", "Analyte", "Mean")], 
                   SubjSamp ~ Analyte, value.var="Mean")
SFN.final <- merge(Meta, SFN.final, by="SubjSamp")
write.csv(SFN.final, "SFN final cortisol data.csv", row.names=F)

SFN.long <- merge(Meta, SFN.recalc[,c("SubjSamp", "Analyte", "Mean", "SD", "CV.recalc")], by="SubjSamp")
names(SFN.long)[names(SFN.long)=="CV.recalc"] <- "CV"
write.csv(SFN.long, "SFN final cortisol data - long format with mean sd CV.csv", row.names=F)

save.image("20140423 SFN checking cortisol values.RData")

