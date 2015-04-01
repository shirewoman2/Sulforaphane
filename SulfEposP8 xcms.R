# SulfEposP8 xcms
# 10/1/14 LS

# Notes to the users of this script ---------------------------------


###       ITEMS YOU WILL WANT TO ADJUST EACH TIME YOU RUN THIS SCRIPT:       ###
# 1. Set the working directory (RawDataDir).

# 2. If you have more than 25 or so samples, split up your samples when you run
# xcmsSet() so that you save your progress occasionally. For example, set up 
# the code like this:
#     xs1 <- xcmsSet(Samples[1:25], method = "centWave",  ppm=15, peakwidth=c(4,12), 
#     snthresh = 20, mzCenterFun="apex", prefilter=c(10, 10000),
#     integrate = 1, fitgauss= TRUE)
#     save(xs1, file="xcmsSet objects.RData")
#
#     xs2 <- xcmsSet(Samples[26:50], method = "centWave",  ppm=15, peakwidth=c(4,12), 
#     snthresh = 20, mzCenterFun="apex", prefilter=c(10, 10000),
#     integrate = 1, fitgauss= TRUE)
#     save(c(xs1, xs2), file="xcmsSet objects.RData")
# etc. until you've done peak picking on all of your samples. That way, if the
# computer shuts down or crashes, you've saved your progress and don't have to 
# start over. To load those files into memory, here is the code:
#     load("xcmsSet objects.RData")

# 3. Parameters I commonly tweak to get things to look right:
#     xcmsSet() -- snthresh and prefilter
#     group() -- bw

# 4. Find and replace "SulfEposP8" with some single-word name that has meaning to
# you and doesn't start with a number or symbol (those are the rules for naming
# objects in R). Now, your objects will be clearly named and the output files
# will be similarly named.

# 5. I added a requirement for the quality checking step that any mass feature
# we look at must be present in at least 25% of all samples. You can change that
# number to whatever you want, though. See the "Quality control" section. 

# 6. Set the ionization mode (see the section "Housekeeping") as appropriate 
# for your data for CAMERA to work properly and for the appropriate set of 
# reference ions to be removed from your data.



###                          GENERAL TIPS                                   ###
# Do yourself a favor and make sure that ALL the files you want to extract 
# data from and ONLY the files you want to extract data from are all in the 
# same folder. It will make your life vastly less complicated. 



# Housekeeping -------------------------------------------------
# This script uses the package xcms to generate a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled. It then performs
# a check on the quality of the data extraction, collapses ions that are likely
# isotopic peaks, and checks on the quality of the selection of isotopes.

# Dataset: SulfEposP8

IonizationMode <- "positive" # Change to "positive" or "negative" as appropriate.

library(RANN)
library(Biobase)
library(seqinr)
library(plyr)
library(lubridate)
library(timeDate)
library(ggplot2)
library(stringr)
library(reshape2)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(xcms)
library(CAMERA)


# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
RawDataDir <- "G:/Data/Metabolomics/Laura/Sulforaphane project/20130604 Sulf EposP mzData files"

MainDir <- "C:/Users/Laura/Documents/Sulforaphane project/8 SulfEposP"

# Set the working directory
setwd(RawDataDir) 

# Getting samples and peak picking ------------------------------------------
# Looking for the appropriate files
Samples <- list.files( getwd(), pattern="mzdata.xml", full.names=F, 
                       recursive = TRUE) 

# Getting the names of the output tables' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(Samples))

# Making a note of start time for this step
tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
xs1 <- xcmsSet(Samples[1:20], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs1, file="xs1.RData")


xs2 <- xcmsSet(Samples[21:40], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs2, file="xs2.RData")

xs3 <- xcmsSet(Samples[41:60], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs3, file="xs3.RData")

xs4 <- xcmsSet(Samples[61:80], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs4, file="xs4.RData")

xs5 <- xcmsSet(Samples[81:100], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs5, file="xs5.RData")

xs6 <- xcmsSet(Samples[101:120], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs6, file="xs6.RData")

xs7 <- xcmsSet(Samples[121:140], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs7, file="xs7.RData")

xs8 <- xcmsSet(Samples[141:160], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs8, file="xs8.RData")

xs9 <- xcmsSet(Samples[161:169], method = "centWave",  ppm=PPM, 
               peakwidth=c(4,12), 
               snthresh = 15, mzCenterFun="apex", prefilter=c(10, 1000),
               integrate = 1, fitgauss= TRUE)
save(xs9, file="xs9.RData")


# Tracking how long this step takes
tpick.final <- Sys.time()
tpick <- as.numeric(difftime(tpick.final, tpick.init, units="mins")) 

write.csv(tpick, "tpick.csv")


# Initial peak alignment -------------------------------------------------------
tgroup.init <- Sys.time()

# Peak alignment
SulfEposP8.grouped <- group(c(xs1, xs2, xs3, xs4, xs5, xs6, xs7, xs8, xs9), 
                            method = "density", bw = 4, minfrac = 0,
                        minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(xs1, xs2, .. xsn), method = ...)

tgroup.final <- Sys.time()
tgroup <- as.numeric(difftime(tgroup.final, tgroup.init, units="mins"))
write.csv(tgroup, "tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(xs1, xs2, xs3, xs4, xs5, xs6, xs7, xs8, xs9)


# RT correction -------------------------------------------------------
tretcor.init <- Sys.time() 

# Performing retention time correction
SulfEposP8.RTcor <- retcor(SulfEposP8.grouped, method = "peakgroups", 
                       missing = 10, extra = 50, smooth = "loess", 
                       family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("SulfEposP8 RTcor plot.png", width = 6, height = 8, units = "in", res = 300)
plotrt(SulfEposP8.RTcor, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

tretcor.final <- Sys.time()
tretcor <- as.numeric(difftime(tretcor.final, tretcor.init, units="mins"))
write.csv(tretcor, "tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
SulfEposP8.grouped2 <- group(SulfEposP8.RTcor, method = "density", 
                         minsamp = 1, minfrac = 0, mzwid = 0.007, 
                         bw = 2, max=100)

# Making a table of the data before the recursive peak filling step
SulfEposP8.unfilled <- peakTable(SulfEposP8.grouped2)

setwd(MainDir)
write.csv(SulfEposP8.unfilled, 
          "SulfEposP8 peak table - unfilled.csv")
setwd(RawDataDir)

tgroup2.final <- Sys.time()
tgroup2 <- as.numeric(difftime(tgroup2.final, tgroup2.init, units="mins"))
write.csv(tgroup2, "tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
SulfEposP8.filledpeaks <- fillPeaks(SulfEposP8.grouped2)

tfillPeaks.final <- Sys.time()
tfillPeaks <- as.numeric(difftime(tfillPeaks.final, tfillPeaks.init, 
                                  units="mins"))
write.csv(tfillPeaks, "tfillPeaks.csv")

# Generate a table with all the peaks ----------------------------------------
tPeakTable.init <- Sys.time()

# Making a table of the recursively filled data
SulfEposP8.allpeaks <- peakTable(SulfEposP8.filledpeaks, 
                             filebase="SulfEposP8 peak table")

# Making a column with the mass feature name
SulfEposP8.allpeaks$MassFeature <- paste("I", round((
      SulfEposP8.allpeaks$mz),digits=4), "R", round(
            SulfEposP8.allpeaks$rt/60, digits=2), sep="")

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
SulfEposP8.allpeaks$RT <- SulfEposP8.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "positive",
       RefIons <- c(121.0509, 922.0098), # ESI+       
       RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(SulfEposP8.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 SulfEposP8.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
SulfEposP8.noref <- SulfEposP8.allpeaks[-unlist(RefMFs), ]

setwd(MainDir)
write.csv(SulfEposP8.noref, "SulfEposP8 peak table without ref ions.csv")

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
SulfEposP8.after2 <- subset(SulfEposP8.noref, SulfEposP8.noref$rt 
                        > 120)
write.csv(SulfEposP8.after2, 
          "SulfEposP8 peak table - RT above 2 min.csv")

# Counting how many times a mass feature was detected BEFORE recursive peak
# filling step.
Count.df <- SulfEposP8.unfilled[as.numeric(row.names(
      SulfEposP8.after2)), SampCol]

SulfEposP8.after2$Count <- rep(NA, nrow(SulfEposP8.after2))

for (m in 1:nrow(Count.df)){
      SulfEposP8.after2$Count[m] <- length(which(!is.na(as.numeric(
            Count.df[m, ]))))
}

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
SulfEposP8.filter <- SulfEposP8.after2[SulfEposP8.after2$Count >= 
                                     0.25*length(SampCol), ]
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
rm(Count.df)
write.csv(SulfEposP8.filter, 
          "SulfEposP8 peak table - mass features must be present in at least 25% of all samples.csv")

setwd(MainDir)
save(SulfEposP8.after2, SulfEposP8.allpeaks, SulfEposP8.filter, SulfEposP8.noref,
     SulfEposP8.unfilled, file = "SulfEposP8 - all main dataframes.RData")
save(SulfEposP8.filter, file = "SulfEposP8 - filtered dataframe only.RData")


tPeakTable.final <- Sys.time()
tPeakTable <- as.numeric(difftime(tPeakTable.final, tPeakTable.init, 
                                  units="mins"))
write.csv(tPeakTable, "tPeakTable.csv")

# Quality control ----------------------------------------
tQC.init <- Sys.time()

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
MFs <- as.numeric(sample(row.names(SulfEposP8.filter), 30))
RandSamp <- as.numeric(sample(length(Samples), 10))
write.csv(MFs, paste(Sys.Date(), 
                     "SulfEposP8 randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "SulfEposP8 randomly selected samples.csv"))

EIC.uncorrected <- list()
EIC.corrected <- list()

# Looking at the extracted ion chromatograms of each sample and mass feature.
# This will take some time to process, so don't expect instant results. 
for (m in MFs){
      EIC.uncorrected[[m]] <- getEIC(SulfEposP8.filledpeaks, rt = "raw", 
                                     groupidx = m, sampleidx = RandSamp)
      EIC.corrected[[m]] <- getEIC(SulfEposP8.filledpeaks, rt = "corrected", 
                                   groupidx = m, sampleidx = RandSamp)  
}

# Setting the palette for drawing the EICs. The last sample in "RandSamp"
# will be colored red and will be the one used to generate the raw mass 
# spectral data.
ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp) - 1), "red")

xset.raw <- xcmsRaw(Samples[RandSamp[10]], profstep = 0.01, profmethod = "bin")

setwd(MainDir)
pdf("SulfEposP8 peak picking and alignment quality check.pdf", 8.5, 11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the last sample in RandSamp. A dashed 
# horizontal line shows where the calculated m/z is; a dashed vertical line 
# shows the median RT for that compound. A semi-transparent gray rectangle shows 
# the compound's m/z within a window equal to the value of PPM.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:30){
      m <- MFs[i]
      plot(EIC.uncorrected[[m]], SulfEposP8.filledpeaks, groupidx = 1, 
           rtrange = 60, col = MyColors, main = MFs[m])
      mtext(paste(i, SulfEposP8.allpeaks$MassFeature[m]), side = 3, line = -1, 
            adj = 0, padj = 0, cex = 0.8)
      plot(EIC.corrected[[m]], SulfEposP8.filledpeaks, groupidx = 1, 
           rtrange = 60, col = MyColors)
      
      RT <- SulfEposP8.allpeaks$rt[m]
      RTRange <- c(RT - 30, RT + 30)
      
      mz <- SulfEposP8.allpeaks$mz[m]
      mzRange <- c(mz - 0.02, mz + 0.02)
      mzRange.poly.low <- mz - mz * 7.5/1e6
      mzRange.poly.up <- mz * 7.5/1e6 + mz
      
      plotRaw(xset.raw, mzrange = mzRange, rtrange = RTRange, log = FALSE)
      abline(h = mz, lty = 2, col = "gray35")
      mtext(paste("abund =", round(SulfEposP8.allpeaks[m, (length(RandSamp))], 
                                   digits = 0)), 
            side = 3, line = -1, adj = 0, padj = 0, cex = 0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, mzRange.poly.low, 
                mzRange.poly.low), 
              col = col2alpha("blue", alpha = 0.1), border=NA)
      abline(v = RT, lty = 2, col = "gray35")
      
}

dev.off()

setwd(RawDataDir)

tQC.final <- Sys.time()
tQC <- as.numeric(difftime(tQC.final, tQC.init, units = "mins"))
write.csv(tQC, "tQC.csv")

# Processing the data using CAMERA -----------------------------
tCAMERA.init <- Sys.time()

# Create a CAMERA object
SulfEposP8.xsa <- xsAnnotate(SulfEposP8.filledpeaks)

# Find the compound groups based on RT
SulfEposP8.xsaF <- groupFWHM(SulfEposP8.xsa, perfwhm = 1.5)

# Check that peaks in the same group correlate well enough to continue to be 
# included. Generate pseudospectra.
SulfEposP8.xsaC <- groupCorr(SulfEposP8.xsaF, cor_eic_th = 0.85, 
                         pval=0.0001, calcIso=TRUE, calcCiS=TRUE)

# Find the isotopes within each pseudospectrum.
SulfEposP8.xsaFI <- findIsotopes(SulfEposP8.xsaC, maxcharge = 3, maxiso = 8, 
                             ppm = 10, mzabs = 0.015, intval = "maxo", 
                             minfrac = 0.25)

# Find the adducts within each pseudospectrum.
SulfEposP8.xsaFA <- findAdducts(SulfEposP8.xsaFI, ppm = 10, mzabs = 0.015, 
                            multiplier = 3, polarity = IonizationMode, 
                            rules = NULL, max_peaks = 100)

# Make a table of the annotated peaks.
SulfEposP8.annot <- getPeaklist(SulfEposP8.xsaFA)

# Once again, remove reference ions from the peak table.
MF.refs <- list()

for (i in 1:length(RefIons)){
      MF.refs[[i]] <- which(SulfEposP8.annot$mz < 
                                  (RefIons[i] + PPM/1e6*RefIons[i]) &
                                  SulfEposP8.annot$mz > 
                                  RefIons[i] - PPM/1e6*RefIons[i])
}

SulfEposP8.annot <- SulfEposP8.annot[-unlist(MF.refs), ]

# Save the original annotated peak list.
setwd(MainDir)
write.csv(SulfEposP8.annot, file = "SulfEposP8 annotated peaklist.csv")
setwd(RawDataDir)

# Make the column "isotopes" be character data instead of the default, factor.
# Replace empty cells with "NA".
SulfEposP8.annot$isotopes <- as.character(SulfEposP8.annot$isotopes)
SulfEposP8.annot$isotopes[SulfEposP8.annot$isotopes == ""] <- NA

# Processing CAMERA data to remove redundant, isotopic mass features.
# Making a column with the mass feature name.
SulfEposP8.annot$MassFeature <- paste0("I", round((SulfEposP8.annot$mz),
                                              digits = 4), "R", 
                                   round(SulfEposP8.annot$rt/60, 
                                         digits = 2))

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
SulfEposP8.annot$RT <- round(SulfEposP8.annot$rt/60, digits = 2)

# Take a subset of the data that is only the mass features that appear
# to have isotopes. 
Iso <- SulfEposP8.annot[, c("MassFeature", "mz", "RT", "isotopes", 
                        "pcgroup")]
Iso <- Iso[complete.cases(Iso$isotopes), ]

# Split the isotopes column into three pieces:
# 1. Iso$IsoGroup (numeric isotope group)
# 2. Iso$IonType (the isotope of that particular ion s/a "M" or "M+1")
# 3. Iso$Charge (the charge of the isotope)
IsoString <- str_split(Iso$isotopes, "\\]")
String1 <- sapply(IsoString, function(x) x[[1]])
String2 <- sapply(IsoString, function(x) x[[2]])
String3 <- sapply(IsoString, function(x) x[[3]])
IsoString.df <- data.frame("IsoGroup" = String1, 
                           "IonType" = String2, 
                           "Charge" = String3)
IsoString.df$IsoGroup <- gsub("\\[", "", IsoString.df$IsoGroup)
IsoString.df$IonType <- gsub("\\[", "", IsoString.df$IonType)

# Put the three new columns back into the subsetted data
Iso <- data.frame(Iso, IsoString.df)

# Add back the other fields to the data.frame
Iso <- merge(Iso, SulfEposP8.annot, by=c("MassFeature", "mz", "RT", 
                                     "isotopes", "pcgroup"), 
             type = "left")
rm(IsoString, String1, String2, String3, IsoString.df)

# For each detected isotope group, calculate the sum of all MS peaks. Sort the
# data.frame by the isotope group number. This data.frame, IsoCollapse, 
# includes only one summed value for each isotope group and it's THAT value
# that I want to use in the final data -- not the value in each of the mass
# features that were isotopes of each other. 
# IsoCollapse has one row per isotope group.
IsoCollapse <- ddply(Iso[, c("IsoGroup", SampCol)], "IsoGroup", 
                     function(x) colSums(x[, SampCol]))
IsoCollapse$IsoGroup <- as.numeric(IsoCollapse$IsoGroup)
IsoCollapse <- arrange(IsoCollapse, IsoGroup)

# Make a data.frame that lists each mass feature with its isotope group.
IsoGroup.df <- data.frame("MassFeature" = rep(NA, nrow(IsoCollapse)),
                          "IsoGroup" = IsoCollapse$IsoGroup)

# For each isotope group, sort the data by m/z and make the name of the 
# new mass feature (really, the group of mass features that are isotopes of one
# another) be the name of the major ion (the lowest m/z ion).
for (i in 1:nrow(IsoCollapse)){
      
      df <- Iso[Iso$IsoGroup == IsoCollapse$IsoGroup[i], ]
      df <- arrange(df, mz)
      IsoGroup.df$MassFeature[i] <- df$MassFeature[1]
      rm(df)
}

# Merge the two data.frames so that we've got all the data we need.
IsoCollapse <- join(IsoGroup.df, IsoCollapse, by="IsoGroup")
IsoCollapse <- join(SulfEposP8.annot[, c("MassFeature", "mz", "RT", 
                                     "pcgroup",
                                     "isotopes", "adduct")], 
                    IsoCollapse, 
                    by="MassFeature", type="right")

# Get all the mass features that weren't isotopes. Assign their isotope group 
# as "NA".
SulfEposP8.collapse <- SulfEposP8.annot[is.na(SulfEposP8.annot$isotopes), ]
SulfEposP8.collapse$IsoGroup <- rep(NA, nrow(SulfEposP8.collapse))

# Add rows to the SulfEposP8.collapse data.frame that are the collapsed 
# isotopic mass features.
SulfEposP8.collapse <- rbind(SulfEposP8.collapse[, c(
      "MassFeature", "mz", "RT", "pcgroup", "isotopes", "adduct", "IsoGroup",
      SampCol)], IsoCollapse)

# Remove referece ions using the "Count" data calculated earlier in this script.
# (Only non-reference mass features will have a value for "Count".)
SulfEposP8.filter2 <- join(SulfEposP8.after2[, c("MassFeature", 
                                         "Count")],
                       SulfEposP8.collapse, by = "MassFeature", 
                       type = "right")
SulfEposP8.filter2 <- SulfEposP8.filter2[
      complete.cases(SulfEposP8.filter2$Count), ]

# Retain only mass features with RT >= 2 minutes.
SulfEposP8.filter2 <- SulfEposP8.filter2[
      SulfEposP8.filter2$Count >= 0.25*length(SampCol) & 
            SulfEposP8.filter2$RT >= 2, ]

setwd(MainDir)
write.csv(SulfEposP8.filter2, "SulfEposP8 - isotopes collapsed.csv")
setwd(RawDataDir)

tCAMERA.final <- Sys.time()
tCAMERA <- as.numeric(difftime(tCAMERA.final, tCAMERA.init, units="mins"))
write.csv(tCAMERA, "tCAMERA.csv")


# CAMERA quality check ------------------------------------------
tQCCAMERA.init <- Sys.time()

# List all the pseudospectra. (See user's manual for the package CAMERA for a 
# definition of "pseudospectra".)
PS <- as.numeric(unique(SulfEposP8.filter2$pcgroup[
      complete.cases(SulfEposP8.filter2$isotopes)]))

# Randomly select 32 pseudospectra to scrutinize.
PS.rand <- sample(PS, 32)

setwd(MainDir)
pdf("SulfEposP8 CAMERA quality check.pdf", height=11, width=8.5)
par(mfrow=c(4,2))

# For each pseudospectrum, plot the EICs and the mass spectrum.
for (p in PS.rand){
      
      plotEICs(SulfEposP8.xsaFA, pspec=p, maxlabel=5)      
      plotPsSpectrum(SulfEposP8.xsaFA, pspec=p, maxlabel=5) 
      
}

dev.off()
setwd(RawDataDir)

tQCCAMERA.final <- Sys.time()
tQCCAMERA <- as.numeric(difftime(tQCCAMERA.final, tQCCAMERA.init, units="mins"))
write.csv(tQCCAMERA, "tQCCAMERA.csv")

# Calculating processing times for each step --------------------------------------------
tpick <- read.csv("tpick.csv")
tgroup <- read.csv("tgroup.csv")
tretcor <- read.csv("tretcor.csv")
tgroup2 <- read.csv("tgroup2.csv")
tfillPeaks <- read.csv("tfillPeaks.csv")
tPeakTable <- read.csv("tPeakTable.csv")
tQC <- read.csv("tQC.csv")
tCAMERA <- read.csv("tCAMERA.csv")
tQCCAMERA <- read.csv("tQCCAMERA.csv")

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- rbind(tpick, tgroup, tretcor, tgroup2, tfillPeaks, tPeakTable, tQC, 
               tCAMERA, tQCCAMERA)
Times$Step <- c("pick peaks", # tpick
                "group 1", # tgroup
                "retcor", # tretcor
                "group 2", # tgroup2
                "fill peaks", # tfillPeaks
                "peak table", # tPeakTable
                "QC", # tQC
                "CAMERA", # tCAMERA
                "QC CAMERA") # tQCCAMERA

Times <- rename(Times, c("x" = "Duration"))

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "SulfEposP8 processing times.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SulfEposP8") +
      guides(fill=FALSE)
ggsave("SulfEposP8 processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SulfEposP8") +
      guides(fill=FALSE)
ggsave("SulfEposP8 processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                                 # "Number of MFs after collapsing isotopic peaks"),
                      "Count" = c(nrow(SulfEposP8.unfilled),
                                  nrow(SulfEposP8.noref),
                                  nrow(SulfEposP8.after2),
                                  nrow(SulfEposP8.filter)))
                                  # nrow(SulfEposP8.filter2)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "SulfEposP8 counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("SulfEposP8 bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(RawDataDir)
save(SulfEposP8.filledpeaks, file = "SulfEposP8 xcmsSet object.RData")
save.image("SulfEposP8 workspace.RData") # This saves EVERYTHING that is currently
# in your workspace, which is a pretty huge file. Skip this step if you don't 
# think you'll need that. 
