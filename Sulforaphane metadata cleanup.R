# Sulforaphane metadata cleanup
# 3/27/14 LS

# Housekeeping =========================================================
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(lubridate)
library(xlsx)

MainDir <- "C:/Users/Laura/Documents/Sulforaphane project"


# SFN urine: Files it appears Tauri used ===================================
# Working through the myriad SFN urine injections and making sure
# that I'm using the correct ones. 
setwd("H:/Data copied 09_2011_B/Sulphuraphane")

TauriCefs <- list.files(pattern = ".cef$")
TauriCefs2 <- gsub("_1.cef", ".cef", TauriCefs)
TauriCefs <- data.frame(File2 = TauriCefs2, 
                        TauriUsed = "yes")
TauriCefs$File <- sub(".cef", "", TauriCefs$File2)
TauriCefs <- unique(TauriCefs)

# write.csv(TauriCefs, "cefs Tauri made.csv", row.names = F)

# Samples and Files ====================================================
# All .d files we have
setwd("G:/Data/Metabolomics/Laura/Sulforaphane project/All SFN d files")
Files <- read.xlsx("All SFN d files known 20150327.xlsx", sheetName = "Sheet1", 
               startRow = 2)
names(Files) <- c("Sample", "File", "DateTime", "Method")
Files <- unique(Files)

# Removing .d to match file names elsewhere
Files$File <- sub(".d", "", Files$File)

# Sample type
Files$SampType <- "clinical"
Files$SampType[str_detect(Files$Sample, "QC")] <- "QC"
Files$SampType[str_detect(Files$Sample, "[Mm]aster QC")] <- "Master QC"
Files$SampType[str_detect(Files$Sample, "blank")] <- "blank"
Files$SampType[str_detect(Files$Sample, "IS")] <- "blank"
Files$SampType[str_detect(Files$Sample, "[Ww]ater")] <- "water"

# Subject
Files$Subject <- str_extract(Files$Sample, " [0-9]{4} |^[0-9]{4} |^[0-9]{4}$")
Files$Subject[Files$SampType != "clinical"] <- NA
Files$Subject <- str_trim(Files$Subject)

# Sample
Files$SampCode <- str_extract(Files$Sample, " [0-9]{3} | [0-9]{3}$|[0-9]{3}-")
Files$SampCode <- sub("-", "", Files$SampCode)
Files$SampCode <- str_trim(Files$SampCode)

setwd(MainDir)
Wklst20111201 <- read.xlsx("20111201 sulforaphane plasma ESI- worklist.xlsx",
                           sheetName = "Sheet1")

# Matching up the missing sample codes from this worklist.
Wklst20111201 <- plyr::rename(Wklst20111201, c("File.name" = "File"))
Files <- join(Files, Wklst20111201[, c("File", "Sample.Code")], by = "File")
rm(Wklst20111201)

for (i in 1:nrow(Files)){
      if(is.na(Files$SampCode[i])){
            Files$SampCode[i] <- Files$Sample.Code[i]
      }
}

Files$Sample.Code <- NULL

# Any missing sample codes for clinical samples?
NoSampCode <- subset(Files, is.na(SampCode) & SampType == "clinical")
# I have no information anywhere that I can find for the file "20110526U-_001".
# This was one that Tauri ran, so I'll need to ask him about it later. The 
# sample column just says "2". 

Files$Mode[Files$Method %in% c("metabolomicsbotA.m", "metabolomicsbotAdivB.m",
                               "metabolomics+ESI.m")] <- "Epos"
Files$Mode[Files$Method == "metabolomicsnegESI.m"] <- "Eneg"
Files$Mode[Files$Method == "metabolomicsAPCI.m"] <- "Apos"
# 1 file was missing method, but it's in the middle of Epos injections.
Files$Mode[Files$File == "20111110W+_001"] <- "Epos"

# Runs
Files$Run <- paste(Files$Mode, str_sub(Files$File, 1, 8), sep = ".")
unique(Files$Run)

# Matrix
Files$Matrix[Files$Run %in% c("Epos.20100929", "Epos.20101109", 
                              "Epos.20110419", "Epos.20110420",
                              "Epos.20110421", "Eneg.20110525",
                              "Eneg.20110526", "Eneg.20110524",
                              "Epos.20110311", "Epos.20110422",
                              "Eneg.20110531", "Apos.20110531",
                              "Apos.20110601", "Apos.20110602",
                              "Apos.20110060")] <- "urine"
Files$Matrix[Files$Run %in% c("Epos.20111108", "Epos.20111109",
                              "Epos.20111110", "Eneg.20111114",
                              "Eneg.20111130", "Eneg.20111201",
                              "Eneg.20111205")] <- "plasma"
Files$Matrix[Files$SampType == "Master QC"] <- "urine"
Files$Matrix[Files$SampType %in% c("blank", "water")] <- "water"

# Adding sample ID
Files$SampleID <- paste0("S", Files$Subject, ".", Files$SampCode, ".", 
                         toupper(str_sub(Files$Matrix, 1,1)))
Files$SampleID[Files$SampType != "clinical"] <- 
      paste(toupper(str_sub(Files$SampType[Files$SampType != "clinical"], 
                                1, 1)), 
            rank(Files$DateTime[Files$SampType != "clinical"]), sep = ".")

# Adding file and mode unique ID
Files$SampModeMat <- paste(Files$SampleID, Files$Mode, Files$Matrix, sep = ".")

# Noting which files Tauri used
Files <- join(Files, TauriCefs[, c("TauriUsed", "File")], by = "File")
setdiff(TauriCefs$File, Files$File)
# Adding note about which files we should use
Files$Use <- "use"

Files <- arrange(Files, DateTime)

# Runs df
Runs <- unique(Files[Files$SampType == "clinical", c("Run", "Mode", "Matrix")])

# Subjects & samples ---------------------------------------------------
setwd(MainDir)
Subjects <- read.xlsx("All sulforaphane sample metadata v3.xlsx", 
                  sheetName = "lookup", 
                  rowIndex = 8:32, colIndex = 8:17)
Subjects <- plyr::rename(Subjects, c("Height.cm" = "Height",
                                     "Weight..Kg" = "Weight"))


EffectCode <- read.xlsx("All sulforaphane sample metadata v3.xlsx", 
                        sheetName = "lookup", 
                        rowIndex = 2:6, colIndex = 8:9)
EffectCode <- plyr::rename(EffectCode, c("Effector.code" = "TxCode"))
EffectCode$EffectAbb <- c("BR", "B", "R")


MDZ <- read.xlsx("All sulforaphane sample metadata v3.xlsx", 
                 sheetName = "lookup", 
                 rowIndex = 36:180, colIndex = 8:12)
MDZ <- plyr::rename(MDZ, c("MDZ.Cl.calculated.by.Yvonne" = "MDZCL"))
MDZ$MDZCL <- as.numeric(MDZ$MDZCL)
MDZ$Sample <- NULL

# Making samples data.frame
Samples <- data.frame(Subject = rep(Subjects$Subject),
                      SampCode = rep(paste0(rep(1:3, each = 2), rep(c(1,8), 3), 
                                            rep(c(0,1), each = 6)),
                                     each = length(Subjects$Subject)))

# Extracting out which treatment each sample was based on the codes
Samples$TxNum <- as.numeric(str_sub(Samples$SampCode, 1, 1))
Samples <- join(Samples, Subjects[, c("Subject", "Order")], by = "Subject")
Samples$TxCode <- str_sub(Samples$Order, Samples$TxNum, Samples$TxNum)
Samples <- join(Samples, EffectCode, by = "TxCode")

# Matrix
Samples$Matrix[str_detect(Samples$SampCode, "0$")] <- "plasma"
Samples$Matrix[str_detect(Samples$SampCode, "1$")] <- "urine"

# Treatment day
Samples$Day <- str_sub(Samples$SampCode, 2, 2)

# Sometimes didn't have the right sample and had to use a different one.
# Noting that.
Samples$SampCode[Samples$Subject == 2233 & Samples$TxNum == 1 &
                       Samples$Day == 1 & Samples$Matrix == "plasma"] <- "111"

# Time point
Samples$TimePt.code <- str_sub(Samples$SampCode, 3,3)
Samples$TimePt[Samples$TimePt.code == "0"] <- "baseline"
Samples$TimePt[Samples$TimePt.code == "1" & 
                     Samples$Matrix == "urine"] <- "0-6 hr"
Samples$TimePt[Samples$TimePt.code == "1" & 
                     Samples$Matrix == "plasma"] <- "15 min"

# Adding MDZ CL info
Samples <- join(Samples, MDZ, by = c("Subject", "Effector", "Day"))

# Unique identifier: SampleID
Samples$SampleID <- paste0("S", Samples$Subject, ".", Samples$SampCode,
                           ".", toupper(str_sub(Samples$Matrix, 1,1)))


# Checking data integrity ------------------------------------------
# There should be 1 file in ESI+ and 1 in ESI- for
# each of these samples. Checking. 
ClinFiles <- arrange(Files[Files$SampType == "clinical" & Files$Use == "use", ], 
                     SampleID, Matrix, Mode)
FileCheck <- count(ClinFiles, SampModeMat)
PossProbs <- subset(Files, SampModeMat %in% 
                           FileCheck$SampModeMat[FileCheck$n != 1])
PossProbs <- arrange(PossProbs, SampModeMat)
PossProbs <- subset(PossProbs, Use == "use")

# Runs that appear to have been redone because they were problematic:
BadRuns <- c("Epos.20100929", "Epos.20110311", "Eneg.20111114")
Files$Use[Files$Run %in% BadRuns] <- "do not use"
Files$Note[Files$Run %in% BadRuns] <- "appears to be a bad run"

# Files acquired with the method "metabolomicsbotA.m" seem to be off. Removing them.
Files$Use[Files$Method %in% 
                c("metabolomicsbotA.m", "metabolomicsbotAdivB.m")] <- "do not use"
Files$Note[Files$Method %in% 
                 c("metabolomicsbotA.m", "metabolomicsbotAdivB.m")] <- "bad method"

# Files where there are replicates
ProbFiles <- c("20111201P-_005", "20111201P-_007", "20100929 ESI+ u_26",
               "20111110P+_005", "20111110P+_007", "20111205P-_036", 
               "20111205P-_038", "20111130P-_005", "20111130P-_007",
               "20111130P-_012", "20111110P+_075", "20111110P+_077", 
               "20111108P+_07", "20111108P+_09", "20111201P-_037",
               "20111201P-_039", "20111110P+_037", "20111110P+_039",
               "20111109P+_06", "20111109P+_08", "20111205P-_004", 
               "20111205P-_006", "20101109 ESI+ b-04", "20110419U+_046b",
               "20110601U+_040b", "20110525U-_006b", "20110601U+_001a",
               "20110419U+_031b", "20110419U+_020b", "20110419U+_021b",
               "20110524U-_039b", "20110531U-_001b", "20110524Q-_039",
               "20110526U-_023b", "20110421U+_011b")
Files$Use[Files$File %in% ProbFiles] <- "do not use"
Files$Note[Files$File %in% ProbFiles] <- "replicate injection"

# Files that were failed or weird injections 
BadInj <- c("20110419U+_046", "20110420U+_042", "20110602U+_018")

Files$Use[Files$File %in% BadInj] <- "do not use"
Files$Note[Files$File %in% BadInj] <- "failed injection"


# Checking again, now that I've removed problem injections.
ClinFiles <- arrange(Files[Files$SampType == "clinical" & Files$Use == "use", ], 
                     SampleID, Matrix, Mode)
FileCheck <- count(ClinFiles, SampModeMat)
PossProbs <- subset(Files, SampModeMat %in% 
                           FileCheck$SampModeMat[FileCheck$n != 1])
PossProbs <- arrange(PossProbs, SampModeMat)
PossProbs <- subset(PossProbs, Use == "use")
# No problems remain. Could still be missing files, though.

# Adding sample info to Files df to check for missing files
AllPossFiles <- join(Files, Samples[, c("Subject", "SampCode", "EffectAbb", "Day", 
                                        "TimePt.code", "Matrix", "SampleID")], 
                     by = c("Subject", "SampCode", "Matrix", "SampleID"), 
                     type = "full")
AllPossFiles <- arrange(AllPossFiles, SampleID)

# Samples that have no files
MissingFiles <- subset(AllPossFiles, is.na(File))

# Checking how many samples have 2 files that we SHOULD use in Epos and Eneg
SampCheck <- count(subset(Files, Use == "use" & SampType == "clinical" &
                                Mode %in% c("Epos", "Eneg")), 
                   SampleID)
# Missing some samples still. Will have to look into this later.

# Saving main metadata --------------------------------------------
setwd(MainDir)
save(Files, Samples, Subjects, Runs, file = "SFN metadata.RData")


