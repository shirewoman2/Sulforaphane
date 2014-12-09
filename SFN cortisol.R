# SFN cortisol
# 1/8/14 LS

library(ggplot2)
library(reshape2)


setwd("D:/My Documents/Work/Lin Lab/Brocco rif cortisol project")
SFN.all <- read.csv("SFN cortisol.csv", header=T, skip=1)
# SFN.all includes all the samples, including subj 2246, for whom we don't have MDZ CL

# Removes subject 2246
SFN <- subset(SFN.all, MDZCL != "NA")


# Calculating the molar ratios for metabolite to parent
SFN$Ratio.6bOHE.E <- SFN$Conc.6bOHE/SFN$Conc.E
SFN$Ratio.6bOHF.F <- SFN$Conc.6bOHF/SFN$Conc.F
SFN$Ratio.all <- (SFN$Conc.6bOHE+SFN$Conc.6bOHF)/(SFN$Conc.E+SFN$Conc.F)
SFN$ClinDay <- as.factor(SFN$ClinDay)

# Linear regression of all samples against MDZ CL
fit.F <- lm(SFN$MDZCL ~ SFN$Ratio.6bOHF.F)
fit.E <- lm(SFN$MDZCL ~ SFN$Ratio.6bOHE.E)
fit.all <- lm(SFN$MDZCL ~ SFN$Ratio.all)

cor.F <- cor(SFN$MDZCL, SFN$Ratio.6bOHF.F)
cor.E <- cor(SFN$MDZCL, SFN$Ratio.6bOHE.E)
cor.all <- cor(SFN$MDZCL, SFN$Ratio.all)



# Simple plot of E ratios
plot(x=SFN$Ratio.6bOHE.E, y=SFN$MDZCL)
abline(fit.E)


# Setting up theme for ggplot2
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
      plot.background = element_rect(fill="white"),
      panel.border = element_rect(color="black", fill=NA)
    )   
}

theme_set(ThemeLaura())


# Cortisol
Xtext <- max(SFN$Ratio.6bOHF.F)*0.9
Ytext <- max(SFN$MDZCL)*0.95

ggplot(SFN, aes(x=SFN$Ratio.6bOHF.F, y=SFN$MDZCL, 
                color=SFN$Tx, shape=SFN$ClinDay, 
                label=rep(1,nrow(SFN)), group=rep(1,nrow(SFN)))) +
  geom_point() + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("6b-OH-F/F molar ratio") + ylab("MDZ CL/F") +
  stat_smooth(method=lm, col="black") + 
  ggtitle("Cortisol") +
  annotate("text", x=Xtext, y=Ytext, label= paste("p =",
             signif(summary(fit.F)$coef[2,"Pr(>|t|)"], digits=2),
             "\nr =",signif(cor.F, digits=2)), 
           fontface="italic", size=3.5)     
ggsave("SNF cortisol correlations.png", width=7, height=4.5)

# Cortisone
Xtext <- max(SFN$Ratio.6bOHE.E)*0.9
Ytext <- max(SFN$MDZCL)*0.95

ggplot(SFN, aes(x=SFN$Ratio.6bOHE.E, y=SFN$MDZCL, 
                color=SFN$Tx, shape=SFN$ClinDay, 
                label=rep(1,nrow(SFN)), group=rep(1,nrow(SFN)))) +
  geom_point() + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("6b-OH-E/E molar ratio") + ylab("MDZ CL/F") +
  stat_smooth(method=lm, col="black") + 
  ggtitle("Cortisone") +
  annotate("text", x=Xtext, y=Ytext, label= 
             paste("p =",
                   signif(summary(fit.E)$coef[2,"Pr(>|t|)"], digits=2),
                   "\nr =",signif(cor.E, digits=2)), 
           fontface="italic", size=3.5)     
ggsave("SNF cortisone correlations.png", width=7, height=4.5)



# Cortisol + cortisone
Xtext <- max(SFN$Ratio.all)*0.9
Ytext <- max(SFN$MDZCL)*0.95

ggplot(SFN, aes(x=SFN$Ratio.all, y=SFN$MDZCL, 
                color=SFN$Tx, shape=SFN$ClinDay, 
                label=rep(1,nrow(SFN)), group=rep(1,nrow(SFN)))) +
  geom_point() + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("Metabolites-to-parent molar ratio") + ylab("MDZ CL/F") +
  stat_smooth(method=lm, col="black") + 
  ggtitle("Cortisol + Cortisone") +
  annotate("text", x=Xtext, y=Ytext, label= 
             paste("p =",
                   signif(summary(fit.E)$coef[2,"Pr(>|t|)"], digits=2),
                   "\nr =",signif(cor.E, digits=2)), 
           fontface="italic", size=3.5)     
ggsave("SNF cortisol + cortisone correlations.png", width=7, height=4.5)




# All 3 plots faceted.
SFN.melt <- melt(SFN[,c(1,4,5,6,11,12,13)], 
                 id=c("SubjSamp", "ClinDay", "Tx", "MDZCL"))
names(SFN.melt)[5:6] <- c("Analyte", "Ratio")


ggplot(SFN.melt, aes(x=SFN.melt$Ratio, y=SFN.melt$MDZCL, 
                color=SFN.melt$Tx, shape=SFN.melt$ClinDay)) + 
  geom_point() + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("Metabolite-to-parent molar ratio") + ylab("MDZ CL/F") +
  facet_wrap( ~ Analyte, nrow=1)
# Nah. It doesn't work b/c the x axis scales differ so much.

################################################################
# GEE, before-and-after plots
################################################################
library(geepack)

SFN.gee.F <- geeglm(SFN$MDZCL ~ SFN$Ratio.6bOHF.F, 
                              id=interaction(SFN$Subject),
                              family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.F)
p.gee.F <- signif(coef(summary(SFN.gee.F))[2,4], digits=2)


SFN.gee.E <- geeglm(SFN$MDZCL ~ SFN$Ratio.6bOHE.E, 
                    id=interaction(SFN$Subject),
                    family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.E)
p.gee.E <- signif(coef(summary(SFN.gee.E))[2,4], digits=2)


SFN.gee.all <- geeglm(SFN$MDZCL ~ SFN$Ratio.all, 
                    id=interaction(SFN$Subject),
                    family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.all)
p.gee.all <- signif(coef(summary(SFN.gee.all))[2,4], digits=2)



################################################################
save.image("SFN cortisol.RData")


