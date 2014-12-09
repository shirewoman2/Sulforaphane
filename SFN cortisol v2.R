# SFN cortisol v2
# 4/23/14 LS


library(ggplot2)
library(reshape2)

###############################################
# Housekeeping
########################
setwd("D:/Users/Laura/Documents/Work/Lin Lab/Sulforaphane project")
SFN.all <- read.csv("SFN cortisol.csv", header=T, skip=1,
                    colClasses=c(rep("factor", 5),rep("numeric",5)))
# SFN.all includes all the samples, including subj 2246, for whom we don't have MDZ CL

# Removes subject 2246
SFN <- subset(SFN.all, MDZCL != "NA")

# Calculating the molar ratios for metabolite to parent
SFN$Ratio.6bOHE.E <- SFN$C6bOHE/SFN$CE
SFN$Ratio.6bOHF.F <- SFN$C6bOHF/SFN$CF
SFN$Ratio.all <- (SFN$C6bOHE+SFN$C6bOHF)/(SFN$CE+SFN$CF)

# Adding a categorical variable to determine whether induced or baseline
SFN$Induced <- rep(NA, nrow(SFN))
SFN$Induced[(SFN$Tx=="rif" | SFN$Tx=="brocco/rif") & SFN$ClinDay==8] <- "induced"
SFN$Induced[SFN$ClinDay==1] <- "baseline"
SFN$Induced[SFN$Tx=="brocco"] <- "baseline"
SFN$Induced <- factor(SFN$Induced, levels=c("baseline", "induced"))

summary(SFN)

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
      plot.background = element_rect(fill="white", colour=NA),
      panel.border = element_rect(color="black", fill=NA),
      strip.background = element_rect(color=NA, fill="white"),
      legend.background = element_rect(color=NA, fill="white"),
      legend.key = element_rect(color=NA, fill="white")
    )   
}

theme_set(ThemeLaura())




##################################################
# Scatter and density plots
#########################
# Getting data into shape
SFN.melt <- melt(SFN[,-c(7:10)], 
                 id=c("SubjSamp", "Subject", "Sample", "ClinDay", "Tx", "MDZCL", "Induced"))
names(SFN.melt)[names(SFN.melt)=="variable"] <- "Analyte"
names(SFN.melt)[names(SFN.melt)=="value"] <- "Ratio"

levels(SFN.melt$Analyte)[levels(SFN.melt$Analyte)=="Ratio.6bOHF.F"] <- "Cortisol molar ratio"
levels(SFN.melt$Analyte)[levels(SFN.melt$Analyte)=="Ratio.6bOHE.E"] <- "Cortisone molar ratio"
levels(SFN.melt$Analyte)[levels(SFN.melt$Analyte)=="Ratio.all"] <- "Combined molar ratio"

levels(SFN.melt$Analyte) <- c("Cortisol molar ratio", 
                              "Cortisone molar ratio",
                              "Combined molar ratio")
sapply(SFN.melt, class)



ggplot(SFN.melt, aes(x=Ratio, y=MDZCL, group=Analyte,
                color=Tx, shape=ClinDay)) + 
  geom_point() + stat_smooth(method=lm, col="black", group="Analyte") + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("Metabolite-to-parent molar ratio") + ylab("MDZ CL/F") +
  facet_grid(. ~ Analyte, scales="free")

ggsave("SNF all scatter plots.png", width=12, height=4.5, dpi=300)


# Density plots
ggplot(SFN.melt, aes(x=Ratio, group=Analyte)) +
  geom_density(fill="darkgray") + 
  facet_wrap(~ Analyte, scales="free", ncol=3) +
  ggtitle("Density plots of molar ratios")

ggsave("Density plots of molar ratios.png", width=8, height=4.5, dpi=300)


ggplot(SFN.melt, aes(x=log10(Ratio), group=Analyte)) +
  geom_density(fill="darkgray") + 
  facet_wrap(~ Analyte, scales="free", ncol=3) +
  ggtitle("Density plots of log-transformed molar ratios")

ggsave("Density plots of log-transformed molar ratios.png", width=8, height=4.5, dpi=300)



# log transformed plots
ggplot(SFN.melt, aes(x=log10(Ratio), y=log10(MDZCL), group=Analyte,
                     color=Tx, shape=ClinDay)) + 
  geom_point() + stat_smooth(method=lm, col="black", group="Analyte") + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("log10(Metabolite-to-parent molar ratio)") + ylab("log10(MDZ CL/F)") +
  facet_grid(. ~ Analyte, scales="free")

ggsave("SNF all scatter plots - log transformed.png", width=12, height=4.5, dpi=300)



ggplot(SFN.melt, aes(x=log10(Ratio), y=log10(MDZCL), group=Induced,
                     color=Tx, shape=ClinDay)) + 
  geom_point() + stat_smooth(method=lm, col="black", group="Induced") + 
  scale_shape_manual(values=c(1,16)) +
  scale_color_manual(values=c("limegreen", "blue", "red")) +
  labs(color="Treatment", shape="Clinical Day") +
  xlab("log10(Metabolite-to-parent molar ratio)") + ylab("log10(MDZ CL/F)") +
  facet_grid(. ~ Analyte, scales="free")

ggsave("SNF all scatter plots - log transformed - regression line by induction status.png", width=12, height=4.5, dpi=300)




################################################################
# GEE
################################################################
library(geepack)

SFN.gee.F <- geeglm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHF.F), 
                              id=interaction(SFN$Subject),
                              family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.F)
p.gee.F <- signif(coef(summary(SFN.gee.F))[2,4], digits=2)


SFN.gee.E <- geeglm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHE.E), 
                    id=interaction(SFN$Subject),
                    family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.E)
p.gee.E <- signif(coef(summary(SFN.gee.E))[2,4], digits=2)


SFN.gee.all <- geeglm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.all), 
                    id=interaction(SFN$Subject),
                    family=gaussian, corstr="ar1", std.err="jack")
summary(SFN.gee.all)
p.gee.all <- signif(coef(summary(SFN.gee.all))[2,4], digits=2)

gee.pvals <- c(p.gee.F, p.gee.E, p.gee.all)

################################################################
# Multiple linear regression with subject as covariate
################################################################

SFN.lm.F <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHF.F) + SFN$Subject)
summary(SFN.lm.F)
p.lm.F <- signif(coef(summary(SFN.lm.F))[2,4], digits=2)


SFN.lm.E <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHE.E) + SFN$Subject)
summary(SFN.lm.E)
p.lm.E <- signif(coef(summary(SFN.lm.E))[2,4], digits=2)


SFN.lm.all <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.all) + SFN$Subject)
summary(SFN.lm.all)
p.lm.all <- signif(coef(summary(SFN.lm.all))[2,4], digits=2)

lm.pvals <- c(p.lm.F, p.lm.E, p.lm.all)



################################################################
# Simple linear regression
################################################################

SFN.slm.F <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHF.F))
summary(SFN.slm.F)
p.slm.F <- signif(coef(summary(SFN.slm.F))[2,4], digits=2)


SFN.slm.E <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.6bOHE.E))
summary(SFN.slm.E)
p.slm.E <- signif(coef(summary(SFN.slm.E))[2,4], digits=2)


SFN.slm.all <- lm(log10(SFN$MDZCL) ~ log10(SFN$Ratio.all))
summary(SFN.slm.all)
p.slm.all <- signif(coef(summary(SFN.slm.all))[2,4], digits=2)

slm.pvals <- c(p.slm.F, p.slm.E, p.slm.all)


################################################################
# Multiple linear regression with subject + induction as covariates
################################################################

SFN.lm2.F <- lm(log10(MDZCL) ~ log10(Ratio.6bOHF.F) + Subject + Induced, data=SFN)
summary(SFN.lm2.F)
p.lm2.F <- signif(coef(summary(SFN.lm2.F))[2,4], digits=2)


SFN.lm2.E <- lm(log10(MDZCL) ~ log10(Ratio.6bOHE.E) + Subject + Induced, data=SFN)
summary(SFN.lm2.E)
p.lm2.E <- signif(coef(summary(SFN.lm2.E))[2,4], digits=2)


SFN.lm2.all <- lm(log10(MDZCL) ~ log10(Ratio.all) + Subject + Induced, data=SFN)
summary(SFN.lm2.all)
p.lm2.all <- signif(coef(summary(SFN.lm2.all))[2,4], digits=2)

lm2.pvals <- c(p.lm2.F, p.lm2.E, p.lm2.all)
lm2.pvals


################################################################
# Multiple linear regression with subject as covariate -- baseline only
################################################################
SFN.baseline <- SFN[SFN$Induced=="baseline",]


SFN.lm3.F <- lm(log10(MDZCL) ~ log10(Ratio.6bOHF.F) + Subject, 
                data=SFN.baseline)
summary(SFN.lm3.F)
p.lm3.F <- signif(coef(summary(SFN.lm3.F))[2,4], digits=2)


SFN.lm3.E <- lm(log10(MDZCL) ~ log10(Ratio.6bOHE.E) + Subject, 
                data=SFN.baseline)
summary(SFN.lm3.E)
p.lm3.E <- signif(coef(summary(SFN.lm3.E))[2,4], digits=2)


SFN.lm3.all <- lm(log10(MDZCL) ~ log10(Ratio.all) + Subject, 
                  data=SFN.baseline)
summary(SFN.lm3.all)
p.lm3.all <- signif(coef(summary(SFN.lm3.all))[2,4], digits=2)

lm3.pvals <- c(p.lm3.F, p.lm3.E, p.lm3.all)
lm3.pvals


################################################################
# Full model multiple linear regression with subject as covariate -- baseline only
################################################################
SFN.base.melt <- melt(SFN.baseline[,-c(7:10,14)], 
                      id=c("SubjSamp", "Subject", "Sample",
                           "Tx", "ClinDay", "MDZCL"))
names(SFN.base.melt)[names(SFN.base.melt)=="variable"] <- "Analyte"
names(SFN.base.melt)[names(SFN.base.melt)=="value"] <- "Ratio"

SFN.lm3 <- lm(log10(MDZCL) ~ log10(Ratio) + Subject + Analyte, 
              data=SFN.base.melt)
summary(SFN.lm3)
p.lm3 <- signif(coef(summary(SFN.lm3))[2,4], digits=2)

# test whether to include subject
ReducedModel <- update(SFN.lm3, .~. -Subject)
anova(SFN.lm3, ReducedModel)

# Analysis of Variance Table
# 
# Model 1: log10(MDZCL) ~ log10(Ratio) + Subject + Analyte
# Model 2: log10(MDZCL) ~ log10(Ratio) + Analyte
# Res.Df  RSS  Df Sum of Sq    F Pr(>F)    
# 1    250 1.40                              
# 2    272 5.18 -22     -3.78 30.8 <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Have to include Subject b/c ANOVA p < 2e-16.


################################################################
save.image("SFN cortisol v2.RData")


