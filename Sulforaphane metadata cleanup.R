# Sulforaphane metadata cleanup
# 10/9/14 LS

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
# OrigDir <- "H:/Data copied 09_2011_B/Sulphuraphane"
# PosDir <- "F:/Sulforaphane/20110419 sulforaphane urine ESI+ and APCI+"
# NegDir <- "F:/Sulforaphane/20110524 sulforaphane urine ESI-"



# SFN urine: Files it appears Tauri used ===================================
# Working through the myriad SFN urine injections and making sure
# that I'm using the correct ones. 
setwd("H:/Data copied 09_2011_B/Sulphuraphane")

TauriCefs <- list.files(pattern = ".cef$")
TauriCefs2 <- gsub("_1.cef", ".cef", TauriCefs)
TauriCefs <- data.frame(File2 = TauriCefs2, 
                        Used = "yes")

# write.csv(TauriCefs, "cefs Tauri made.csv", row.names = F)


# Samples and Files ====================================================
setwd(MainDir)

# All possible files that I know of (checking on whether these match
# what Tauri has)
Samples <- read.csv("SFN metadata.csv", na.strings = c("Has not run", "ND"),
                    skip = 1)
Samples <- Samples[, c("Subject", "SampleCode", "Effector", "Day", 
                       "TimePoint", "Matrix", "MDZCL", "Gender",
                       "DatePrepped", "File.Epos", "Date.Epos", "File.Eneg",
                       "Date.Eneg", "File.Apos", "Date.Apos")]
Samples$TimePoint[Samples$TimePoint == "" & Samples$Matrix == "plasma"] <- 
      "5 min"

Samples$SampleID <- paste(paste0("S", Samples$Subject), 
                          Samples$Effector,
                          Samples$Day, 
                          Samples$TimePoint, 
                          Samples$Matrix, sep = ".")
Samples$SampleID <- sub("brocco/rif", "BR", Samples$SampleID)
Samples$SampleID <- sub("brocco", "B", Samples$SampleID)
Samples$SampleID <- sub("rif", "R", Samples$SampleID)

Samples$SampleID <- sub("baseline", "BL", Samples$SampleID)
Samples$SampleID <- sub("0-6 hrs", "06", Samples$SampleID)
Samples$SampleID <- sub("6-24 hrs", "624", Samples$SampleID)
Samples$SampleID <- sub("5 min", "5", Samples$SampleID)

Samples$SampleID <- sub("plasma", "P", Samples$SampleID)
Samples$SampleID <- sub("urine", "U", Samples$SampleID)

Samples$SampleID.code <- paste("S", Samples$Subject, Samples$SampleCode, 
                               Samples$Matrix, sep = ".")

# Removing samples we never ran at all
NeverRan <- Samples[is.na(Samples$File.Epos) & 
                          is.na(Samples$File.Eneg) &
                          is.na(Samples$File.Apos), ]

Samples <- Samples[!(Samples$SampleID.code %in% NeverRan$SampleID.code), ]

Files <- Samples[, c("Subject", "Matrix", "DatePrepped", "File.Epos",
                     "File.Eneg", "File.Apos", "Date.Epos", "Date.Eneg", 
                     "Date.Apos", "SampleID", "SampleID.code")] %>% 
      gather(key, value, -Subject, -Matrix, -DatePrepped, -SampleID,
             -SampleID.code) %>% 
      # Separate the key into components type and mode
      separate(key, c("type", "Mode"), "\\.") %>%
      # Spread the type back into the columns
      spread(type, value)
Files <- Files[complete.cases(Files$File), ]

ddply(Files, c("Mode", "Matrix"), function(x) nrow(x))
# Good. Have all the modes and matrices expected.

# Checking that all the files exist.
setwd("F:/Sulforaphane/All SFN mzdata files")
Files$Exists <- file.exists(paste0(Files$File, ".mzdata.xml"))

Missing <- subset(Files, Exists == FALSE & Mode != "Apos")


# Keeping in Samples df only columns that make sense to be there.
Samples <- Samples[, c("Subject", "SampleCode", "Effector", "Day", "TimePoint",
                       "Matrix", "SampleID", "SampleID.code", "MDZCL", "Gender")]

# Checking for replicates
which(duplicated(Samples$SampleID))
# Good. No replicates.


# Saving main metadata --------------------------------------------
setwd(MainDir)
save(Files, Samples, file = "SFN metadata.RData")

# Making DAReprocessor worklist ===========================================
Worklist <- read.xlsx2("Sulforaphane urine worklists.xlsx", 
                       sheetName = "worklists", 
                       colClasses = rep("character", 5),
                       stringsAsFactors = FALSE)

Worklist$File2 <- gsub("d", "cef", basename(Worklist$File))

Worklist <- join(Worklist, TauriCefs, by = "File2", type = "left")

Worklist$Subject <- rep(NA, nrow(Worklist))

Clin <- which(str_detect(Worklist$Sample, "Subject"))

for (i in Clin){
      Worklist$Subject[i] <- str_sub(Worklist$Sample[i],
                                     str_locate(Worklist$Sample[i], 
                                                "Subject")[2]+2,
                                     str_locate(Worklist$Sample[i],
                                                "Subject")[2]+5)
}

Worklist$SampleType <- rep(NA, nrow(Worklist))
Worklist$SampleType[Clin] <- "clinical"

QC <- which(str_detect(Worklist$Sample, "QC"))
Worklist$SampleType[QC] <- "QC"

MasterQC <- which(str_detect(Worklist$Sample, "Master QC"))
Worklist$SampleType[MasterQC] <- "Master QC"

Worklist$Sample2 <- rep(NA, nrow(Worklist))

for (i in Clin){
      Worklist$Sample2[i] <- str_sub(Worklist$Sample[i],
                                     str_locate(Worklist$Sample[i], 
                                                "Subject")[2]+9,
                                     str_locate(Worklist$Sample[i],
                                                "Subject")[2]+11)
}

Worklist$Dataset <- rep(NA, nrow(Worklist))
Worklist$Dataset[Worklist$Method == "metabolomics+ESI.m"] <- "EposU"
Worklist$Dataset[Worklist$Method == "metabolomicsnegESI.m"] <- "EnegU"
Worklist$Dataset[Worklist$Method == "metabolomicsAPCI.m"] <- "AposU"


# Checking to see if there are any clinical samples that were NOT used before 
# that are in the worklists. 
which(is.na(Worklist$Used) & Worklist$SampleType == "Clin")
# good!

Worklist$SampleID <- paste(Worklist$Subject, Worklist$Sample2)

Worklist$File3 <- rep(NA, nrow(Worklist))
Worklist$File3[Worklist$Dataset == "EposU" | Worklist$Dataset == "AposU"] <- 
      paste0("F:\\Sulforaphane\\20110419 sulforaphane urine ESI+ and APCI+\\", 
                         basename(Worklist$File[Worklist$Dataset == "EposU" | 
                                                      Worklist$Dataset == 
                                                      "AposU"]))

Worklist$File3[Worklist$Dataset == "EnegU"] <- 
                     paste0("F:\\Sulforaphane\\20110524 sulforaphane urine ESI-\\", 
             basename(Worklist$File[Worklist$Dataset == "EnegU"]))

DAReproc <- Worklist[Worklist$SampleType %in% c("clinical", "Master QC", "QC"), 
                     c("Sample", "Vial.position", "Method", "File3")]
DAReproc$Method <- rep("C:\\Users\\Laura\\Documents\\MassHunter methods\\General mzData export - 500 count peak height.m",
                       nrow(DAReproc))

write.csv(DAReproc, "20141010 DAReprocessor worklist for SFN urine.csv", 
          row.names = FALSE)



