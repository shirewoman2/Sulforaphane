# SulfEposU10 xcms
# 3/28/15 LS

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

# 4. Find and replace "SulfEposU10" with some single-word name that has meaning to
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

# Dataset: SulfEposU10

IonizationMode <- "positive" # Change to "positive" or "negative" as appropriate.

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
MainDir <- "C:/Users/Laura/Documents/Sulforaphane project/SFNEposU10"
RawDataDir <- "G:/Data/Metabolomics/Laura/Sulforaphane project/All SFN mzdata files"

# Loading metadata
setwd("C:/Users/Laura/Documents/Sulforaphane project")
load("SFN metadata.RData")

# Samples to use
GoodSamples <- paste0(Files$File[Files$Mode == "Epos" & 
                                       Files$Matrix == "urine" &
                                       Files$Use == "use" &
                                       Files$SampType == "clinical"],
                  ".mzdata.xml")

# Getting the names of the output data.frames' sample columns
SampCol <- sub("mzdata.xml", "mzdata", make.names(GoodSamples))

# Getting the name of the column that will contain the number of times
# an ion was detected.
CountCol <- make.names(basename(RawDataDir))



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
SulfEposU10.tpick.init <- Sys.time() 

# Setting mass accuracy in ppm. Set this to whatever you think is appropriate 
# for your instrument. Err on the higher side.
PPM <- 15 

# Peak picking
SNthresh <- 15
Prefilter <- c(10, 5000)

SulfEposU10.xs1 <- xcmsSet(GoodSamples[1:25], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs1, file="SulfEposU10 xs1.RData")

SulfEposU10.xs2 <- xcmsSet(GoodSamples[26:50], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs2, file="SulfEposU10 xs2.RData")

SulfEposU10.xs3 <- xcmsSet(GoodSamples[51:75], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs3, file="SulfEposU10 xs3.RData")

SulfEposU10.xs4 <- xcmsSet(GoodSamples[76:100], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs4, file="SulfEposU10 xs4.RData")

SulfEposU10.xs5 <- xcmsSet(GoodSamples[101:125], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs5, file="SulfEposU10 xs5.RData")

SulfEposU10.tpick.5 <- mdy_hms("3/29/15 2:23:00")

SulfEposU10.tpick.6 <- Sys.time()

SulfEposU10.xs6 <- xcmsSet(GoodSamples[126:137], method = "centWave",  ppm=PPM, 
                        peakwidth=c(4,12), 
                        snthresh = SNthresh, mzCenterFun="apex", prefilter = 
                              Prefilter,
                        integrate = 1, fitgauss= TRUE)
save(SulfEposU10.xs6, file="SulfEposU10 xs6.RData")


# Tracking how long this step takes
SulfEposU10.tpick.final <- Sys.time()
SulfEposU10.tpick <- as.numeric(difftime(SulfEposU10.tpick.final, 
                                      SulfEposU10.tpick.6, units="mins")) +
      as.numeric(difftime(SulfEposU10.tpick.5,
                          SulfEposU10.tpick.init, units = "mins"))
write.csv(SulfEposU10.tpick, "SulfEposU10 tpick.csv")



# Initial peak alignment -------------------------------------------------------
tgroup.init <- Sys.time()

# Peak alignment
SulfEposU10 <- group(c(SulfEposU10.xs1, SulfEposU10.xs2, SulfEposU10.xs3,
                    SulfEposU10.xs4, SulfEposU10.xs5, SulfEposU10.xs6), 
                    method = "density", bw = 4, minfrac = 0,
                  minsamp = 1, mzwid = 0.007, max = 100)
# If you split your samples into multiple xcmsSet steps, change the first bit
# of the above command to:
# group (c(SulfEposU10.xs1, SulfEposU10.xs2, .. SulfEposU10.xsn), method = ...)

SulfEposU10.tgroup.final <- Sys.time()
SulfEposU10.tgroup <- as.numeric(difftime(SulfEposU10.tgroup.final, 
                                       SulfEposU10.tgroup.init, units="mins"))
write.csv(SulfEposU10.tgroup, "SulfEposU10 tgroup.csv")

# Removing xcmsSet objects b/c they're such RAM eaters. 
rm(SulfEposU10.xs1, SulfEposU10.xs2, SulfEposU10.xs3,
   SulfEposU10.xs4, SulfEposU10.xs5, SulfEposU10.xs6,
   SulfEposU10.xs7)


# RT correction -------------------------------------------------------
SulfEposU10.tretcor.init <- Sys.time() 

# Performing retention time correction
SulfEposU10 <- retcor(SulfEposU10, method = "peakgroups", 
                   missing = 0.1*length(GoodSamples), 
                   extra = 2*length(GoodSamples), smooth = "loess", 
                   family = "symmetric", plottype = NULL)

# Making a plot of the retention time deviations for each sample
setwd(MainDir)
png("SulfEposU10 RTcor plot.png", width = 4, height = 5, units = "in", res = 300)
plotrt(SulfEposU10, leg = F, densplit = T)
dev.off()
setwd(RawDataDir)

SulfEposU10.tretcor.final <- Sys.time()
SulfEposU10.tretcor <- as.numeric(difftime(SulfEposU10.tretcor.final, SulfEposU10.tretcor.init, units="mins"))
write.csv(SulfEposU10.tretcor, "SulfEposU10 tretcor.csv")

# Peak align the RT-corrected data -----------------------------------------
SulfEposU10.tgroup2.init <- Sys.time()

# Refining peak alignment after RT correction step
SulfEposU10 <- group(SulfEposU10, method = "density", 
                  minsamp = 1, minfrac = 0, mzwid = 0.007, 
                  bw = 2, max=100)

# Making a data.frame of the data before the recursive peak filling step
SulfEposU10.unfilled <- peakTable(SulfEposU10)
SulfEposU10.unfilled$MassFeature <- paste("I", round((
      SulfEposU10.unfilled$mz),digits=4), "R", round(
            SulfEposU10.unfilled$rt/60, digits=2), sep="")

SulfEposU10.tgroup2.final <- Sys.time()
SulfEposU10.tgroup2 <- as.numeric(difftime(SulfEposU10.tgroup2.final, 
                                        SulfEposU10.tgroup2.init, units="mins"))
write.csv(SulfEposU10.tgroup2, "SulfEposU10 tgroup2.csv") 

# Recursive peak filling -------------------------------------------------------
SulfEposU10.tfillPeaks.init <- Sys.time()

# Recursively filling all detected peaks
SulfEposU10 <- fillPeaks(SulfEposU10)

SulfEposU10.tfillPeaks.final <- Sys.time()
SulfEposU10.tfillPeaks <- as.numeric(difftime(SulfEposU10.tfillPeaks.final, 
                                           SulfEposU10.tfillPeaks.init, 
                                           units="mins"))
write.csv(SulfEposU10.tfillPeaks, "SulfEposU10 tfillPeaks.csv")

# Generate a data.frame with all the peaks ----------------------------------
SulfEposU10.tPeakTable.init <- Sys.time()

# Making a data.frame of the recursively filled data
SulfEposU10.allpeaks <- peakTable(SulfEposU10)

# Changing the name of the count column
names(SulfEposU10.allpeaks)[names(SulfEposU10.allpeaks) == CountCol] <- "Count"

# Making a column with the mass feature name
SulfEposU10.allpeaks$MassFeature <- paste("I", round((
      SulfEposU10.allpeaks$mz),digits=4), "R", round(
            SulfEposU10.allpeaks$rt/60, digits=2), sep="")

# Making a column with the mass feature name as xcms sets it
SulfEposU10.allpeaks$groupname <- groupnames(SulfEposU10)

# Making a column with the RT in minutes. Note that this is different
# from the column "rt", which is the RT in seconds. 
SulfEposU10.allpeaks$RT <- SulfEposU10.allpeaks$rt/60

# Removing reference ions
# Checking on which reference ions are in these data
ifelse(IonizationMode == "positive",
       RefIons <- c(121.0509, 922.0098), # ESI+       
       RefIons <- c(112.9856, 119.0363, 980.015)) # ESI-

RefMFs <- list()

# Finding mass features that are really just reference ions
for (m in 1:length(RefIons)){
      
      RefMFs[[m]] <- which(SulfEposU10.allpeaks$mz < 
                                 (RefIons[m] + PPM/1e6*RefIons[m]) &
                                 SulfEposU10.allpeaks$mz > 
                                 RefIons[m] - PPM/1e6*RefIons[m])
      
}

# Removing reference ions from the data
SulfEposU10.noref <- SulfEposU10.allpeaks[-unlist(RefMFs), ]

# Retaining only mass features that elute after 2 min since our RT isn't 
# very reproducible before 2 minutes.
SulfEposU10.after2 <- subset(SulfEposU10.noref, SulfEposU10.noref$rt 
                          > 120)

# Only retaining mass features that were detected in the initial peak-picking 
# step in at least 25% of all samples.
SulfEposU10.filter <- subset(SulfEposU10.after2, Count >= 0.25*length(SampCol))
# If you want to be more or less stringent in what fraction of samples you're
# requiring something to be detected in initially, change the "0.25" to
# something higher or lower. The "0.25" means that a mass feature had to be
# found in at least 25% of all samples to be retained for consideration.
setwd(MainDir)
write.csv(SulfEposU10.filter, 
          "SulfEposU10 peak table - mass features must be present in at least 25% of all samples.csv")

save(SulfEposU10.after2, SulfEposU10.allpeaks, SulfEposU10.filter, SulfEposU10.noref,
     SulfEposU10.unfilled, file = "SulfEposU10 - all main dataframes.RData")
save(SulfEposU10.filter, file = "SulfEposU10 - filtered dataframe only.RData")

setwd(RawDataDir)
SulfEposU10.tPeakTable.final <- Sys.time()
SulfEposU10.tPeakTable <- as.numeric(difftime(SulfEposU10.tPeakTable.final, 
                                           SulfEposU10.tPeakTable.init, 
                                           units="mins"))
write.csv(SulfEposU10.tPeakTable, "SulfEposU10 tPeakTable.csv")

# Quality control ----------------------------------------
SulfEposU10.tQC.init <- Sys.time()

setwd(MainDir)

# Selecting some random mass features and samples to scrutinize and then
# saving the names of those mass features and samples.
set.seed(253)

MFs <- sample(SulfEposU10.filter$groupname, 30)
RandSamp <- sample(1:length(GoodSamples), 10)
write.csv(MFs, paste(Sys.Date(), 
                     "SulfEposU10 randomly selected mass features.csv"))
write.csv(RandSamp, paste(Sys.Date(), 
                          "SulfEposU10 randomly selected samples.csv"))

LastSamp <- RandSamp[length(RandSamp)]

EIC.uncorrected <- list()
EIC.corrected <- list()

# This next step will take some time to process, so don't expect instant results. 
for (i in 1:length(MFs)){
      EIC.uncorrected[[i]] <- getEIC(SulfEposU10, rt="raw", 
                                     groupidx=MFs[i], sampleidx=RandSamp)
      EIC.corrected[[i]] <- getEIC(SulfEposU10, rt="corrected", 
                                   groupidx=MFs[i], sampleidx=RandSamp)  
}

ColRainbow <- colorRampPalette(c("green", "blue", "purple"))
MyColors <- c(ColRainbow(length(RandSamp)-1), "red")

setwd(RawDataDir)
xset.raw <- xcmsRaw(GoodSamples[LastSamp], profstep=0.01, profmethod="bin")

setwd(MainDir)
pdf("SulfEposU10 quality check.pdf", 8.5,11)

# 1st column shows the uncorrected EICs.
# 2nd column shows the RT-corrected EICs.
# 3rd column shows the m/z vs. RT for the 1st sample for that compound with a 
# dashed horizontal line where the calculated m/z is.

par(mfrow=c(4,3), mar=c(3,3,3,0.5))
for(i in 1:length(MFs)){
      
      plot(EIC.uncorrected[[i]], SulfEposU10, groupidx=1, rtrange=60, 
           col=MyColors, main=MFs[i])
      mtext(paste(i, SulfEposU10.filter$MassFeature[
            SulfEposU10.filter$groupname == MFs[i]]), side=3, line=-1, 
            adj=0, padj=0, cex=0.8)
      plot(EIC.corrected[[i]], SulfEposU10, groupidx=1, rtrange=60, 
           col=MyColors)
      
      RT <- SulfEposU10.filter$rt[SulfEposU10.filter$groupname == MFs[i]]
      RTRange <- c(RT-30, RT+30)
      
      mz <- SulfEposU10.filter$mz[SulfEposU10.filter$groupname == MFs[i]]
      mzRange <- c(mz-0.02, mz+0.02)
      mzRange.poly.low <- mz- mz*(0.5*PPM)/1e6
      mzRange.poly.up <- mz*(0.5*PPM)/1e6 + mz
      
      plotRaw(xset.raw, mzrange=mzRange, rtrange=RTRange, log=FALSE)
      abline(h=mz, lty=2, col="gray35")
      mtext(paste("abund =", round(SulfEposU10.filter[
            SulfEposU10.filter$groupname == MFs[i], SampCol[LastSamp]], digits=0)), 
            side=3, line=-1, adj=0, padj=0, cex=0.8)
      polygon(c(RTRange[2], RTRange[1], RTRange[1], RTRange[2]), 
              c(mzRange.poly.up, mzRange.poly.up, 
                mzRange.poly.low, mzRange.poly.low), 
              col=col2alpha("blue", alpha=0.1), border=NA)
      abline(v=RT, lty=2, col="gray35")
      
}

dev.off()

save(EIC.corrected, EIC.uncorrected, xset.raw, 
     file = "SulfEposU10 QC data.RData")

setwd(RawDataDir)
SulfEposU10.tQC.final <- Sys.time()
SulfEposU10.tQC <- as.numeric(difftime(SulfEposU10.tQC.final, 
                                    SulfEposU10.tQC.init, units = "mins"))
write.csv(SulfEposU10.tQC, "SulfEposU10 tQC.csv")



# Calculating processing times for each step --------------------------------
setwd(RawDataDir)

# Make a data.frame with all the times that each step required. 
# Units for Times are in minutes.
Times <- data.frame(rbind(SulfEposU10.tpick, SulfEposU10.tgroup, 
                          SulfEposU10.tretcor, SulfEposU10.tgroup2, 
                          SulfEposU10.tfillPeaks, 
                          SulfEposU10.tPeakTable, SulfEposU10.tQC))
Times$Step <- c("pick peaks", # SulfEposU10.tpick
                "group 1", # SulfEposU10.tgroup
                "retcor", # SulfEposU10.tretcor
                "group 2", # SulfEposU10.tgroup2
                "fill peaks", # SulfEposU10.tfillPeaks
                "peak table", # SulfEposU10.tPeakTable
                "QC") # SulfEposU10.tQC
names(Times) <- c("Duration", "Step")
row.names(Times) <- 1:nrow(Times)

# Calculate the total time in hours. 
TotalTime <- c("Total time in hours", format(sum(Times$Duration)/60, digits=2))
Times$Duration <- round(Times$Duration, digits=2)
Times <- Times[, c("Step", "Duration")]

setwd(MainDir)
write.csv(rbind(Times, TotalTime), "SulfEposU10 processing times.csv", row.names=F)

Times$Step <- factor(Times$Step, levels=Times$Step)

# Make a bar graph showing how long each step took.
ggplot(Times, aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SulfEposU10") +
      guides(fill=FALSE)
ggsave("SulfEposU10 processing times bar plot.png", height=6, width=6.5, dpi=300)

# Since peak picking usually takes SOOOOOO long compared to everything else,
# make another graph that doesn't show the peak picking step.
ggplot(Times[-1,], aes(x=Step, y=Duration, fill=Step)) + 
      geom_bar(stat="identity") +
      ylab("Duration (min)") + ggtitle("Processing times for SulfEposU10") +
      guides(fill=FALSE)
ggsave("SulfEposU10 processing times without peak picking step - bar plot.png", 
       height = 6, width = 6.5, dpi = 300)

# Number of mass features at each step ------------------------------
# Make a data.frame listing the number of mass features at each step.
MFCount <- data.frame("Step" = c("Number of MFs initially", 
                                 "Number of MFs after removing reference masses",
                                 "Number of MFs eluting after 2 min", 
                                 "Number of MFs present in >= 25% of samples at peak-picking step"),
                      "Count" = c(nrow(SulfEposU10.unfilled),
                                  nrow(SulfEposU10.noref),
                                  nrow(SulfEposU10.after2),
                                  nrow(SulfEposU10.filter)))


MFCount$Step <- factor(MFCount$Step, levels = MFCount$Step)

write.csv(MFCount, "SulfEposU10 counts of mass features at each step.csv", 
          row.names=F)


# Make a bar graph showing how many mass features there were at each step.
ggplot(MFCount, aes(x = Step, y = Count, fill = Step)) + 
      geom_bar(stat = "identity") + guides(fill = FALSE) +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
ggsave("SulfEposU10 bar chart of numbers of mass features at each step.png")


# Saving final workspace ------------------------------
setwd(MainDir)
save(SulfEposU10, file = "SulfEposU10 xcmsSet object.RData")
# save.image("SulfEposU10 workspace.RData") # This saves EVERYTHING that is 
# currently in your workspace, which is a pretty huge file. Skip this step if 
# you don't think you'll need that. 
