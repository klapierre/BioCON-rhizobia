################################################################################
##  BioCON_plant soil feedback_Toomer_2023.R: Greenhouse study with three BioCON legume species, inoculated with a mix of rhizobial strains.
##
##  Authors: Kimberly Komatsu
################################################################################

library(nlme)
library(emmeans)
library(performance)
library(readxl)
library(tidyverse)

#kim's working directory
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\UNCG undergraduate students\\Toomer_2023')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=20),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=20),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=16), legend.text=element_text(size=16))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


options(contrasts=c('contr.sum','contr.poly'))

#not in function
`%!in%` = Negate(`%in%`)

###########################################################################
###########################################################################

#### Read in data ####
data <- read_excel('Toomer_legume_coexistence_2023.xlsx', sheet = "compiled_data") %>% 
  na.omit() %>% 
  filter(rhizobia==1) %>% 
  separate(treatment, into=c('sp1', 'sp2'), sep='/', remove=F)


#### Checking model assumptions and cleaning data ####

with(subset(data, sp2=='LP'), hist(aboveground_introduced))
with(subset(data, sp2=='LC'), hist(aboveground_introduced))
with(subset(data, sp2=='AC'), hist(aboveground_introduced))

with(data, hist(aboveground_introduced))
with(data, hist(log(aboveground_introduced)))

with(subset(data, sp2=='LP'), hist(belowground_introduced))
with(subset(data, sp2=='LC'), hist(belowground_introduced))
with(subset(data, sp2=='AC'), hist(belowground_introduced))

with(data, hist(belowground_introduced))
with(data, hist(log(belowground_introduced)))


#### Aboveground Biomass Models and Figures ####

# LP - no effect of sp 1 on aboveground biomass
summary(LPabovegroundIntroduced <- lme(aboveground_introduced ~ as.factor(sp1),
                                   data=subset(data, sp2=='LP' & aboveground_introduced<100), 
                                   random=~1|rep))
check_model(LPabovegroundIntroduced)
anova.lme(LPabovegroundIntroduced, type='sequential')
emmeans(LPabovegroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")

summary(aov(aboveground_introduced ~ as.factor(sp1),
            data=subset(data, sp2=='LP' & aboveground_introduced<100)))

# LC - no effect of sp 1 on aboveground biomass
summary(LCabovegroundIntroduced <- lme(aboveground_introduced ~ as.factor(sp1),
                                       data=subset(data, sp2=='LC'), 
                                       random=~1|rep))
check_model(LCabovegroundIntroduced)
anova.lme(LCabovegroundIntroduced, type='sequential')
emmeans(LCabovegroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")

summary(aov(aboveground_introduced ~ as.factor(sp1),
            data=subset(data, sp2=='LC')))

# AC - no effect of sp 1 on aboveground biomass
summary(ACabovegroundIntroduced <- lme(log(aboveground_introduced) ~ as.factor(sp1),
                                       data=subset(data, sp2=='AC' & aboveground_introduced<30), 
                                       random=~1|rep))
check_model(ACabovegroundIntroduced)
anova.lme(ACabovegroundIntroduced, type='sequential')
emmeans(ACabovegroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")

summary(aov(aboveground_introduced ~ as.factor(sp1),
            data=subset(data, sp2=='AC' & aboveground_introduced<30)))

# all together
summary(abovegroundIntroduced <- lme(aboveground_introduced ~ as.factor(sp1)*as.factor(sp2),
                                       data=data, 
                                       random=~1|rep))
check_model(abovegroundIntroduced)
anova.lme(abovegroundIntroduced, type='sequential')
emmeans(abovegroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")


#### Belowground Biomass Models and Figures ####

# LP - effect of sp 1 on belowground biomass
# highest with AC, middle with self, lowest with LC
summary(LPbelowgroundIntroduced <- lme(belowground_introduced ~ as.factor(sp1),
                                       data=subset(data, sp2=='LP'), 
                                       random=~1|rep))
check_model(LPbelowgroundIntroduced)
anova.lme(LPbelowgroundIntroduced, type='sequential')
emmeans(LPbelowgroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")

# LC - no effect of sp 1 on belowground biomass
summary(LCbelowgroundIntroduced <- lme(log(belowground_introduced) ~ as.factor(sp1),
                                       data=subset(data, sp2=='LC'), 
                                       random=~1|rep))
check_model(LCbelowgroundIntroduced)
anova.lme(LCbelowgroundIntroduced, type='sequential')
emmeans(LCbelowgroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")

# AC - no effect of sp 1 on belowground biomass
summary(ACbelowgroundIntroduced <- lme(belowground_introduced ~ as.factor(sp1),
                                       data=subset(data, sp2=='AC' & belowground_introduced<30), 
                                       random=~1|rep))
check_model(ACbelowgroundIntroduced)
anova.lme(ACbelowgroundIntroduced, type='sequential')
emmeans(ACbelowgroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")


# all together
summary(belowgroundIntroduced <- lme(belowground_introduced ~ as.factor(sp1)*as.factor(sp2),
                                     data=data, 
                                     random=~1|rep))
check_model(belowgroundIntroduced)
anova.lme(belowgroundIntroduced, type='sequential')
emmeans(belowgroundIntroduced, pairwise~as.factor(sp1), adjust="tukey")


#### Nodules Models and Figures ####

# LP - no effect of sp 1 on belowground biomass
summary(LPnodulesIntroduced <- lme(nodules_introduced ~ as.factor(sp1),
                                   data=subset(data, sp2=='LP'), 
                                    random=~1|rep))
check_model(LPnodulesIntroduced)
anova.lme(LPnodulesIntroduced, type='sequential')
emmeans(LPnodulesIntroduced, pairwise~as.factor(sp1), adjust="tukey")

# LC - no effect of sp 1 on belowground biomass
summary(LCnodulesIntroduced <- lme(nodules_introduced ~ as.factor(sp1),
                                   data=subset(data, sp2=='LC'), 
                                   random=~1|rep))
check_model(LCnodulesIntroduced)
anova.lme(LCnodulesIntroduced, type='sequential')
emmeans(LCnodulesIntroduced, pairwise~as.factor(sp1), adjust="tukey")

# AC - no effect of sp 1 on belowground biomass
summary(ACnodulesIntroduced <- lme(nodules_introduced ~ as.factor(sp1),
                                   data=subset(data, sp2=='AC' & nodules_introduced<30), 
                                   random=~1|rep))
check_model(ACnodulesIntroduced)
anova.lme(ACnodulesIntroduced, type='sequential')
emmeans(ACnodulesIntroduced, pairwise~as.factor(sp1), adjust="tukey")

# # all together
# summary(nodulesIntroduced <- lme(nodules_introduced ~ as.factor(sp1)*as.factor(sp2),
#                                  data=data, 
#                                  random=~1|rep))
# check_model(nodulesIntroduced)
# anova.lme(nodulesIntroduced, type='sequential')
# emmeans(nodulesIntroduced, pairwise~as.factor(sp1), adjust="tukey")






