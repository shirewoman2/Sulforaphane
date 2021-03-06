# SFNEnegU11 xcms
# 3/30/15 LS

# Notes to the users of this script ---------------------------------


###       ITEMS YOU WILL WANT TO ADJUST EACH TIME YOU RUN THIS SCRIPT:       ###
# 1. Set the working directory (RawDataDir).

# 2. If you have more than 25 or so samples, split up your samples when you run
# xcmsSet() so that you save your progress occasionally. For example, set up 
# the code like this:
#     xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=15, peakwidth=c(4,12), 
#     snthresh = 20, mzCenterFun="apex", prefilter=c(10, 10000),
#     integrate = 1, fitgauss= TRUE)
#     save(xs1, file="xcmsSet objects.RData")
#
#     xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=15, peakwidth=c(4,12), 
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

# 4. Find and replace "SFNEnegU11" with some single-word name that has meaning to
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
# Do yourself a favor and make a metadata file that includes information on
# each sample. Also include any info you might later want to sort or filter by
# such as matrix, date prepped, date run, what kind of sample is it (clinical,
# QC, Master QC, water), treatment group, whether that injection was a good one
# if some of them were not, etc. Then, use that file throughout your analyses, 
# including here, where you can use it to pick which samples you want to 
# work with. 


# Housekeeping -------------------------------------------------
# This script uses the package xcms to generate a list of peaks that have been
# aligned, RT corrected, realigned, and recursively filled. It then performs
# a check on the quality of the data extraction, collapses ions that are likely
# isotopic peaks, and checks on the quality of the selection of isotopes.

# Dataset: SFNEnegU11

IonizationMode <- "negative" # Change to "positive" or "negative" as appropriate.

library(plyr)
library(ggplot2)
library(gridExtra)
library(xcms)
library(XLConnect)
library(dplyr)
library(seqinr)
library(lubridate)

# Set path of the working directory, i.e. the folder where your files are. 
# Note that you need forward slashes!
MainDir <- "C:/Users/Laura/Documents/Sulforaphane project/SFNEnegU11"
RawDataDir <- "G:/Data/Metabolomics/Laura/Sulforaphane project/All SFN mzdata files"

# Loading metadata
setwd("C:/Users/Laura/Documents/Sulforaphane project")
load("SFN metadata.RData")

# Samples to use
GoodSamples <- paste0(Files$File[Files$Mode == "Eneg" & 
                                       Files$Matrix == "urine" &
                                       Files$Use == "use" &
                                       Files$SampType == "clinical"],
                  ".mzdata.xml")

# Getting the names of the output data.frames' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(GoodSamples))


# General functions ----------------------------------------------
ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=90),
                  panel.background = element_rect(fill="white", color=NA),
                  panel.grid.minor.y = element_line(color = NA),
                  panel.grid.minor.x = element_line(color = NA),
                  panel.grid.major = element_line(colour = NA),
                  plot.background = element_rect(fill="white", colour=NA),
                  panel.border = element_rect(color="black", fill=NA),
                  strip.background = element_rect(color=NA, fill="white"),
                  legend.background = element_rect(color=NA, fill=NA),
                  legend.key = element_rect(color=NA, fill=NA)
            )   
}

# Call up that theme before plotting graphs.
theme_set(ThemeLaura())


# Peak picking ------------------------------------------
setwd(RawDataDir)

# Making a note of start time for this step
SFNEnegU11.tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 15
Prefilter <- c(10, 5000)

SFNEnegU11.xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs1, file="SFNEnegU11 xs1.RData")

SFNEnegU11.xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs2, file="SFNEnegU11 xs2.RData")

SFNEnegU11.xs3 <- xcmsSet(GoodSamples[51:75], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs3, file="SFNEnegU11 xs3.RData")

SFNEnegU11.xs4 <- xcmsSet(GoodSamples[76:100], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs4, file="SFNEnegU11 xs4.RData")

SFNEnegU11.xs5 <- xcmsSet(GoodSamples[101:125], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs5, file="SFNEnegU11 xs5.RData")

SFNEnegU11.xs6 <- xcmsSet(GoodSamples[126:136], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SFNEnegU11.xs6, file="SFNEnegU11 xs6.RData")


# Tracking how long this step takes
SFNEnegU11.tpick.final <- Sys.time()
SFNEnegU11.tpick <- as.numeric(difftime(SFNEnegU11.tpick.final, 
                                        SFNEnegU11.tpick.init, units = "mins"))
write.csv(SFNEnegU11.tpick, "SFNEnegU11 tpick.csv")



# Initial peak alignment -------------------------------------------------------
SFNEnegU11.tgroup.init <- Sys.time()

# Peak alignment
SFNEnegU11 <- group(c(SFNEnegU11.xs1, SFNEnegU11.xs2, SFNEnegU11.xs3,
                    SFNEnegU11.xs4, SFNEnegU11.xs5, SFNEnegU11.xs6), 
                    method = "density", bw = 4, minfrac = 0,
                  minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(SFNEnegU11.xs1, SFNEnegU11.xs2, .. SFNEnegU11.xsn), method = ...)

SFNEnegU11.tgroup.final <- Sys.time()
SFNEnegU11.tgroup <- as.numeric(difftime(SFNEnegU11.tgroup.final, 
                                       SFNEnegU11.tgroup.init, units="mins"))
write.csv(SFNEnegU11.tgroup, "SFNEnegU11 tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(SFNEnegU11.xs1, SFNEnegU11.xs2, SFNEnegU11.xs3,
   SFNEnegU11.xs4, SFNEnegU11.xs5, SFNEnegU11.xs6)


# RT correction -------------------------------------------------------
SFNEnegU11.tretcor.init <- Sys.time() 

# Performing retention time correction
SFNEnegU11 <- retcor(SFNEnegU11, method = "peakgroups", 
                   missing = 0.1*length(GoodSamples), 
                   extra = 2*length(GoodSamples), smooth = "loess", 
                   family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("SFNEnegU11 RTcor plot.png", width = 4, height = 5, units = "in", res = 300)
plotrt(SFNEnegU11, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

SFNEnegU11.tretcor.final <- Sys.time()
SFNEnegU11.tretcor <- as.numeric(
      difftime(SFNEnegU11.tretcor.final, SFNEnegU11.tretcor.init, units="mins"))
write.csv(SFNEnegU11.tretcor, "SFNEnegU11 tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
SFNEnegU11.tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
SFNEnegU11 <- group(SFNEnegU11, method = "density", 
                  minsamp = 1, minfrac = 0, mzwid = 0.007, 
                  bw = 2, max=100)

# Making a data.frame of the data before the recursive peak filling step
SFNEnegU11.unfilled <- peakTable(SFNEnegU11)
SFNEnegU11.unfilled$MassFeature <- paste("I", round((
      SFNEnegU11.unfilled$mz),digits=4), "R", round(
            SFNEnegU11.unfilled$rt/60, digits=2), sep="")

SFNEnegU11.tgroup2.final <- Sys.time()
SFNEnegU11.tgroup2 <- as.numeric(difftime(SFNEnegU11.tgroup2.final, 
                                        SFNEnegU11.tgroup2.init, units="mins"))
write.csv(SFNEnegU11.tgroup2, "SFNEnegU11 tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
SFNEnegU11.tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
SFNEnegU11 <- fillPeaks(SFNEnegU11)

# Saving main object
save(SFNEnegU11, file = "SFNEnegU11 xcmsSet object.RData")

SFNEnegU11.tfillPeaks.final <- Sys.time()
SFNEnegU11.tfillPeaks <- as.numeric(difftime(SFNEnegU11.tfillPeaks.final, 
                                           SFNEnegU11.tfillPeaks.init, 
                                           units="mins"))
write.csv(SFNEnegU11.tfillPeaks, "SFNEnegU11 tfillPeaks.csv")

# Generate a data.frame with all the peaks ----------------------------------
SFNEnegU11.tPeakTable.init <- Sys.time()

# Making a data.frame of the recursively filled data
SFNEnegU11.allpeaks <- peakTable(SFNEnegU11)

# Setting up a function to count length of columns with 0 or NA in the 
# data before recursive peak filling step
Detected <- function(x) {length(x) - length(which(is.na(x) | x == 0))}
# Counting
SFNEnegU11.allpeaks$Count <- apply(SFNEnegU11.unfilled[, SampCol], 
                                   MARGIN = 1, FUN = Detected)

# Making a column with the mass feature name
SFNEnegU11.allpeaks$MassFeature <- paste("I", round((
      SFNEnegU11.allpeaks$mz),digits=4), "R", round(
            SFNEnegU11.allpeaks$rt/60, digits=2), sep="")

# Making a column with the mass feature name as xcms sets it
SFNEnegU11.allpeaks$groupname <- groupnames(SFNEnegU11)

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
SFNEnegU11.allpeaks$RT <- SFNEnegU11.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "positive",
       RefIons <- c(121.0509, 922.0098), # ESI+       
       RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(SFNEnegU11.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 SFNEnegU11.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
SFNEnegU11.noref <- SFNEnegU11.allpeaks[-unlist(RefMFs), ]

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
SFNEnegU11.after2 <- subset(SFNEnegU11.noref, SFNEnegU11.noref$rt 
                          > 120)

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
SFNEnegU11.filter <- subset(SFNEnegU11.after2, Count >= 0.25*length(SampCol))
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
setwd(MainDir)
write.csv(SFNEnegU11.filter, 
          "SFNEnegU11 peak table - mass features must be present in at least 25% of all samples.csv")

save(SFNEnegU11.after2, SFNEnegU11.allpeaks, SFNEnegU11.filter, SFNEnegU11.noref,
     SFNEnegU11.unfilled, file = "SFNEnegU11 - all main dataframes.RData")
save(SFNEnegU11.filter, file = "SFNEnegU11 - filtered dataframe only.RData")

setwd(RawDataDir)
SFNEnegU11.tPeakTable.final <- Sys.time()
SFNEnegU11.tPeakTable <- as.numeric(difftime(SFNEnegU11.tPeakTable.final, 
                                           SFNEnegU11.tPeakTable.init, 
                                           units="mins"))
write.csv(SFNEnegU11.tPeakTable, "SFNEnegU11 tPeakTable.csv")

# Quality control ----------------------------------------
SFNEnegU11.tQC.init <- Sys.time()

setwd(MainDir)

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
set.seed(253)

MFs <- sample(SFNEnegU11.filter$groupname, 30)
RandSamp <- sample(1:length(GoodSamples), 10)
write.csv(MFs, paste(Sys.Date(), 
                     "SFNEnegU11 randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "SFNEnegU11 randomly selected samples.csv"))

LastSamp <- RandSamp[length(RandSamp)]

EIC.uncorrected <- list()
EIC.corrected <- list()

# This next step will take some time to process, so don't expect instant results. 
for (i in 1:length(MFs)){
      EIC.uncorrected[[i]] <- getEIC(SFNEnegU11, rt="raw", 
                                     groupidx=MFs[i], sampleidx=RandSamp)
      EIC.corrected[[i]] <- getEIC(SFNEnegU11, rt="corrected", 
                                   groupidx=MFs[i], sampleidx=RandSamp)  
}

ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp)-1), "red")

setwd(RawDataDir)
xset.raw <- xcmsRaw(GoodSamples[LastSamp], profstep=0.01, profmethod="bin")

setwd(MainDir)
pdf("SFNEnegU11 quality check.pdf", 8.5,11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the 1st sample for that compound with a 
# dashed horizontal line where the calculated m/z is.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:length(MFs)){
      
      plot(EIC.uncorrected[[i]], SFNEnegU11, groupidx=1, rtrange=60, 
           col=MyColors, main=MFs[i])
      mtext(paste(i, SFNEnegU11.filter$MassFeature[
            SFNEnegU11.filter$groupname == MFs[i]]), side=3, line=-1, 
            adj=0, padj=0, cex=0.8)
      plot(EIC.corrected[[i]], SFNEnegU11, groupidx=1, rtrange=60, 
           col=MyColors)
      
      RT <- SFNEnegU11.filter$rt[SFNEnegU11.filter$groupname == MFs[i]]
      RTRange <- c(RT-30, RT+30)
      
      mz <- SFNEnegU11.filter$mz[SFNEnegU11.filter$groupname == MFs[i]]
      mzRange <- c(mz-0.02, mz+0.02)
      mzRange.poly.low <- mz- mz*(0.5*PPM)/1e6
      mzRange.poly.up <- mz*(0.5*PPM)/1e6 + mz
      
      plotRaw(xset.raw, mzrange=mzRange, rtrange=RTRange, log=FALSE)
      abline(h=mz, lty=2, col="gray35")
      mtext(paste("abund =", round(SFNEnegU11.filter[
            SFNEnegU11.filter$groupname == MFs[i], SampCol[LastSamp]], digits=0)), 
            side=3, line=-1, adj=0, padj=0, cex=0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, 
                mzRange.poly.low, mzRange.poly.low), 
              col=col2alpha("blue", alpha=0.1), border=NA)
      abline(v=RT, lty=2, col="gray35")
      
}

dev.off()

save(EIC.corrected, EIC.uncorrected, xset.raw, 
     file = "SFNEnegU11 QC data.RData")

setwd(RawDataDir)
SFNEnegU11.tQC.final <- Sys.time()
SFNEnegU11.tQC <- as.numeric(difftime(SFNEnegU11.tQC.final, 
                                    SFNEnegU11.tQC.init, units = "mins"))
write.csv(SFNEnegU11.tQC, "SFNEnegU11 tQC.csv")



# Calculating processing times for each step --------------------------------
setwd(RawDataDir)

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- data.frame(rbind(SFNEnegU11.tpick, SFNEnegU11.tgroup, 
                          SFNEnegU11.tretcor, SFNEnegU11.tgroup2, 
                          SFNEnegU11.tfillPeaks, 
                          SFNEnegU11.tPeakTable, SFNEnegU11.tQC))
Times$Step <- c("pick peaks", # SFNEnegU11.tpick
                "group 1", # SFNEnegU11.tgroup
                "retcor", # SFNEnegU11.tretcor
                "group 2", # SFNEnegU11.tgroup2
                "fill peaks", # SFNEnegU11.tfillPeaks
                "peak table", # SFNEnegU11.tPeakTable
                "QC") # SFNEnegU11.tQC
names(Times) <- c("Duration", "Step")
row.names(Times) <- 1:nrow(Times)

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "SFNEnegU11 processing times.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SFNEnegU11") +
      guides(fill=FALSE)
ggsave("SFNEnegU11 processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SFNEnegU11") +
      guides(fill=FALSE)
ggsave("SFNEnegU11 processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(SFNEnegU11.unfilled),
                                  nrow(SFNEnegU11.noref),
                                  nrow(SFNEnegU11.after2),
                                  nrow(SFNEnegU11.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "SFNEnegU11 counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("SFNEnegU11 bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(MainDir)
# save.image("SFNEnegU11 workspace.RData") # This saves EVERYTHING that is 
# currently in your workspace, which is a pretty huge file. Skip this step if 
# you don't think you'll need that. 
