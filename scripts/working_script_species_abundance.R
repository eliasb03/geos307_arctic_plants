# USE THIS - 2024 Version 
# Point Frame Data for GEOS 307 Biogeography and Global Change 2024
# Created by Zoe Feb 2020; Adapted by N Hewitt 2021-2; Pedro Gonzalez 2020; Jenna Loesberg 2022

# Step 1: Set your working directory:
### follow along with TA's instructions in lab

getwd() # check and make sure it's where you want it (i.e. where the csv is)

# Step 2: Install and load packages
## Install any packages if you have not previously (remove hashtag and run)

# install.packages("plyr")
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("ggplot2")

## Load Packages:

library(tidyverse)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)

#Step 3: Read in data and make some graphs
#===============================================#
# EXAMPLE: Arctic poppy (PAPRAD) at WILL site  #
#===============================================#

# Load the dataset:
data_path <- "data/SpeciesAbundance_1995-2019.csv"
AlexSpeciesAbund_1995_2019 <- read.csv(data_path, stringsAsFactors = F, strip.white = T, na.strings = c("","NA"))
# read csv into R; species are columns, sites are rows, presences = 0, absences = 1


# This is an example species. You will repeat this with your group's assigned species
# Create a new data frame with only PAPRAD species abundance at CASS site (subset)

DRYAS_SpeciesAbund_1995_2019 <- subset(AlexSpeciesAbund_1995_2019, AlexSpeciesAbund_1995_2019$SUBSITE=="DRYAS", select = c(SUBSITE,TRTMT,PLOT,SPP,YEAR,Abundance))
CASTET_DRYAS_SpeciesAbund_1995_2019 <- subset(DRYAS_SpeciesAbund_1995_2019, DRYAS_SpeciesAbund_1995_2019$SPP=="CASTET", select = c(SUBSITE,TRTMT,PLOT,SPP,YEAR,Abundance))

# Plot PAPRAD abundance difference between control and OTC plots for all years
boxplot(Abundance~TRTMT, data = PAPRAD_WILL_SpeciesAbund_1995_2019, yaxt="n", col=(c("cyan3","darkorange")), boxwex = 0.5, main = "PAPRAD abundance WILL site ")
axis(2,las=2)
# Test for significant difference between control and OTC plots using TukeyHSD test
TukeyHSD(aov(PAPRAD_WILL_SpeciesAbund_1995_2019$Abundance~PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT), 'PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT')

# Plot PAPRAD abundance difference between control and OTC plots for each year
boxplot(Abundance~TRTMT*YEAR, data = PAPRAD_WILL_SpeciesAbund_1995_2019, yaxt="n", xaxt="n", col=(c("cyan3","darkorange")), boxwex = 0.5, main = "PAPRAD abundance WILL site ")
axis(2,las=2)
text(x=1:16, y=-1, labels = c("CTL", "OTC"), xpd = NA )
text(x=1:16, y=-2.5, labels = c("1995","1995","1996","1996","2000","2000","2007","2007","2010","2010","2015","2015","2019","2019"), xpd = NA )

# Plot CTL and OTC abundance trend
plot(PAPRAD_WILL_SpeciesAbund_1995_2019$YEAR, PAPRAD_WILL_SpeciesAbund_1995_2019$Abundance, yaxt="n",type = "p", pch=16, col= ifelse(PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT == "CTL", "cyan3","darkorange"), main = "WILL site PAPRAD abundance trend over time", xlab = "Year", ylab = "Abundance")
axis(2,las=2)
CTL_Abund<-subset(PAPRAD_WILL_SpeciesAbund_1995_2019, PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT=="CTL", select = c(SUBSITE,TRTMT,PLOT,SPP,YEAR,Abundance))
OTC_Abund<-subset(PAPRAD_WILL_SpeciesAbund_1995_2019, PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT=="OTC", select = c(SUBSITE,TRTMT,PLOT,SPP,YEAR,Abundance))
CTLfit<-lm(Abundance~YEAR, data=CTL_Abund)
OTCfit<-lm(Abundance~YEAR, data=OTC_Abund)
abline(CTLfit, col="cyan3")
abline(OTCfit, col="darkorange")
stats.text<-paste0("R^2=", round(summary(CTLfit)$r.squared, digits = 2),", p=",round(cor.test(CTL_Abund$Abundance, CTL_Abund$YEAR)$p.value, digits = 4), ", b=",round(coef(CTLfit)[2], digits = 4))
text(1994.5, 17, stats.text, col="cyan3", cex = 1, pos=4, offset = 0)
stats.text<-paste0("R^2=", round(summary(OTCfit)$r.squared, digits = 2),", p=",round(cor.test(OTC_Abund$Abundance, OTC_Abund$YEAR)$p.value, digits = 4), ", b=",round(coef(OTCfit)[2], digits = 4))
text(1994.5, 19, stats.text, col="darkorange", cex = 1, pos=4, offset = 0)
# Statistis of trend lines
summary(OTCfit)
summary(CTLfit)


# Test for significant difference between control and OTC plots

# A) Using TukeyHSD test
TukeyHSD(aov(PAPRAD_WILL_SpeciesAbund_1995_2019$Abundance~PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT), 'PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT')

# B) Using student's t test
CTL_Abund<-subset(PAPRAD_WILL_SpeciesAbund_1995_2019$Abundance, PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT=="CTL")
OTC_Abund<-subset(PAPRAD_WILL_SpeciesAbund_1995_2019$Abundance, PAPRAD_WILL_SpeciesAbund_1995_2019$TRTMT=="OTC")
t.test(CTL_Abund,OTC_Abund,var.equal=TRUE, paired=FALSE)

### Your turn! Just copy and paste the above analysis below, and change the site and species names to your group's assigned site and species
### This will mean you have to change "WILLOW" and "PAPRAD" in lines 33 and 34 to your assigned species, as well as renaming the new data frames




