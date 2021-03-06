# SFN plasma processing

# This script processes the data from xcms.


# Housekeeping, including data tidying ---------------------------------------
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(lubridate)
library(stringr)

Dataset <- c("SFNEposP8", "SFNEnegP9")

Directory <- list("C:/Users/Laura/Documents/Sulforaphane project/8 SulfEposP",
                  "C:/Users/Laura/Documents/Sulforaphane project/9 SulfEnegP")
names(Directory) <- Dataset

MainDir <- c("C:/Users/Laura/Documents/Sulforaphane project")

File <- c("SulfEposP8 peak table - mass features must be present in at least 25% of all samples.csv",
          "SulfEnegP9 peak table - mass features must be present in at least 25% of all samples.csv")
names(File) <- Dataset

setwd(MainDir)
load("SFN metadata.RData")
Files$SampCol <- make.names(paste0(Files$File, ".mzdata"))
Files$Dataset[Files$Mode == "Epos" & Files$Matrix == "plasma"] <- Dataset[1]
Files$Dataset[Files$Mode == "Eneg" & Files$Matrix == "plasma"] <- Dataset[2]

# Lists:
Data.filtered <- list()
Data.log <- list()

for (j in Dataset){
      setwd(as.character(Directory[j]))
      Data.filtered[[j]] <- read.csv(File[j])      
}

names(Data.filtered) <- Dataset

Samples <- list()
ClinSamples <- list()

# Selecting only the columns of interest in the data and renaming them sensibly
for (j in Dataset){
      
      FilesToUse <- subset(Files, Dataset == j & Use == "use" &
                                 SampType %in% c("clinical", "QC"))
      FilesToUse$Check <- FilesToUse$SampCol %in% names(Data.filtered[[j]])
      FilesToUse <- subset(FilesToUse, Check == TRUE)
      
      Samples[[j]] <- FilesToUse$SampleID
      ClinSamples[[j]] <- FilesToUse$SampleID[FilesToUse$SampType == "clinical"]
      
      Data.filtered[[j]] <- Data.filtered[[j]][, c("MassFeature", "mz", 
                                                   "RT", 
                                                   FilesToUse$SampCol)]
      names(Data.filtered[[j]]) <- c("MassFeature", "mz", "RT", 
                                     FilesToUse$SampleID)
      
      
}



# Functions and themes ----------------------------------
# Theme for graphs made using ggplot2
ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=0),
                  panel.background = element_rect(fill="white", color=NA),
                  panel.grid.minor.y = element_line(color="white"),
                  panel.grid.minor.x = element_line(color="white"),
                  panel.grid.major = element_line(colour = "white"),
                  plot.background = element_rect(fill="white", color=NA),
                  panel.border = element_rect(color="black", fill=NA),
                  strip.background = element_rect(color=NA, fill="white"),
                  legend.background = element_rect(color=NA, fill="white"),
                  legend.key = element_rect(color=NA, fill="white")
            )   
}

theme_set(ThemeLaura())


# Function for writing files with a note on the first line. 
my.write <- function(x, file, header, f = write.csv, ...){
      # create and open the file connection
      datafile <- file(file, open = 'wt')
      # close on exit
      on.exit(close(datafile))
      # if a header is defined, write it to the file
      if(!missing(header)) writeLines(header,con=datafile)
      # write the file using the defined function and required addition arguments  
      f(x, datafile,...)
}



# Processing the data ----------------------------------------------

# Starting from data that were already filtered by frequency to retain MFs 
# detected initially in >= 25% of all samples. 

Hist <- list()
Data.TCCnorm <- list()

for (j in Dataset){
      
      # Adding 1 to everything to avoid taking log(0)
      Data.filtered[[j]][, Samples[[j]]] <- Data.filtered[[j]][, Samples[[j]]] + 1
      
      # Calculate the TCC sum.
      TCCsum <- apply(Data.filtered[[j]][ , Samples[[j]]], 
                      2, sum, na.rm=TRUE)
      
      
      # Normalize each sample by the TCC sum. 
      Data.TCCnorm[[j]] <- Data.filtered[[j]][ , Samples[[j]]]
      for (s in 1:length(Samples[[j]])){
            Data.TCCnorm[[j]][, Samples[[j]][s]] <- 
                  (Data.filtered[[j]][, Samples[[j]][s]]/ 
                         TCCsum[s])*1e6
            rm(s)
      }  
      
      
      # log10 transforming
      Data.log[[j]] <- data.frame(MassFeature = Data.filtered[[j]]$MassFeature, 
                                  mz = Data.filtered[[j]]$mz, 
                                  RT = Data.filtered[[j]]$RT,
                                  log10(Data.TCCnorm[[j]]))
      
      
      # Saving the files
      setwd(as.character(Directory[j]))
      FileNote <- paste("Dataset:", j,
                        "This file was generated on",Sys.Date(), 
                        "using the script SCOR MDZ preprocessing.R.")
      my.write(Data.log[[j]], paste(j,"- preprocessed.csv"), 
               header=FileNote, row.names=F)
      
      
      # Checking distributions. 
      Data.melt <- gather(Data.log[[j]], Sample, Abundance, -MassFeature, 
                          -mz, -RT)
      
      
      Hist[[j]] <- ggplot(Data.melt, aes(x=Abundance)) + 
            geom_density() +
            xlab("log10(Abundance)") + ggtitle(j)
      
      rm(Data.melt, TCCsum)
      
}


setwd(MainDir)
png("Kernel density plots of MF abundances.png", height = 6, width = 8,
    units = "in", res=600)
grid.arrange(Hist[[1]], Hist[[2]], ncol=2)
dev.off()

save(Data.log, Directory, Dataset, 
     Samples, ClinSamples, file="SCOR MDZ main data.RData")
save(Data.TCCnorm, Data.filtered, 
     file = "SCOR MDZ intermediate steps to final data.RData")
