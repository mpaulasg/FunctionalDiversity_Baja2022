################################################################################
##
## Script for FD analysis 
##
## Paper: Fish functional diversity is modulated by 
## small-scale habitat complexity in a temperate ecosystem 
## (to submit in Hydrobiologia)
##
## Authors: Sgarlatta, M. Paula, Ramirez-Valdez, Arturo,
## Ladah, Lydia B., Calderon-Aguilera, Luis E.
## 
## Code by Paula Sgarlatta 
##
##
##  Scripts to run:
##
## 1/ Statistical analysis
## 
###############################################################################

rm(list=ls()) # cleaning memory

#load packages

library(here)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(car)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(dplyr)
library(DHARMa)

# Load data

fd_values <- read.csv(here::here("data", "fd_values_replicates.csv"))

metadata <- read.csv(here::here("data", "kelp&rocky_metadata_replicate.csv"))

habitat_kelp <- read.csv(here::here("data", "habitat_kelp.csv"))

habitat_RR <- read.csv(here::here("data", "habitat_RR.csv"))


#Prepare data for stats

fd_stats <- fd_values %>% 
  left_join(metadata, by=c("X"="Code")) %>% 
  select(-X, -HabitatID, -Replicate) %>% 
  relocate(Site, Habitat, .before = "sp_richn")


##### t test between habitats (RR and KF) #####

##Species richness

#Normality - QQ plots

sp_kelp <- fd_stats %>% 
  filter(Habitat =="Kelp")

sp_rocky<- fd_stats %>% 
  filter(Habitat =="Rocky")

ggqqplot(sp_kelp$sp_richn)   
ggqqplot(sp_rocky$sp_richn)    # Both good

t.test(sp_kelp$sp_richn,sp_rocky$sp_richn, var.equal=TRUE)

#No differences between habitats - p=0.46

##Functional richness

#Normality
ggqqplot(sp_kelp$fric)  
ggqqplot(sp_rocky$fric)   # Both good 

t.test(sp_kelp$fric,sp_rocky$fric, var.equal=TRUE)

#No differences between habitats - p=0.74


##Functional evenness

#Normality
ggqqplot(sp_kelp$feve)   
ggqqplot(sp_rocky$feve)  # This one is not good, so we will do
#a Wilcoxon Test 

wilcox.test(sp_kelp$feve,sp_rocky$feve, paired = TRUE)

#No Differences between habitats - p=0.23

##Functional divergence

#Normality
ggqqplot(sp_kelp$fdiv)  
ggqqplot(sp_rocky$fdiv)   #Both good  

t.test(sp_kelp$fdiv,sp_rocky$fdiv, var.equal=TRUE)

### Differences between habitats - p=0.02

ggboxplot(fd_stats, x = "Habitat", y = "fdiv", 
          color = "Habitat",palette = c("#00AFBB", "#E7B800"),
          ylab = "Functional Divergence", xlab = "")


##### ANOVA between sites (within each habitat) #####


############################KELP FORESTS #############################


## S

s_KF <- aov(sp_richn ~ Site, data=sp_kelp)

summary(s_KF)

#Not significant - p=0.937

#Normality

ggqqplot(s_KF$residuals) #Good

#Homogeneity of variances
bartlett.test(sp_richn ~ Site, data=sp_kelp)  #p=0.6148


## Fric

fric_KF <- aov(fric ~ Site, data=sp_kelp)

summary(fric_KF)

#No significant - p=0.427

#Normality
ggqqplot(fric_KF$residuals) #Good

#Homogeneity of variances
bartlett.test(fric ~ Site, data=sp_kelp)  #p=0.1778


## FEve

feve_KF <- aov(feve ~ Site, data=sp_kelp)

summary(feve_KF)

# Significant - p=0.0337

#Normality
ggqqplot(feve_KF$residuals) #Good

#Homogeneity of variances
bartlett.test(feve ~ Site, data=sp_kelp)  #p=0.42


## FDiv


fdiv_KF <- aov(fdiv ~ Site, data=sp_kelp)

summary(fdiv_KF)

###Not significant - p=0.5

#Normality
ggqqplot(fdiv_KF$residuals)

#Homogeneity of variances
bartlett.test(fdiv ~ Site, data=sp_kelp)  #p=0.78

####################ROCKY REEFS ######################

## S

s_RR <- aov(sp_richn ~ Site, data=sp_rocky)
summary(s_RR)

#Not significant - p=0.91

#Normality

ggqqplot(s_RR$residuals) #Good

#Homogeneity of variances
bartlett.test(sp_richn ~ Site, data=sp_rocky)  #p=0.88

## FRic

fric_RR <- aov(fric ~ Site, data=sp_rocky)
summary(fric_RR)

#Not significant - p=0.96

#Normality

ggqqplot(fric_RR$residuals) #Good


#Homogeneity of variances
bartlett.test(fric ~ Site, data=sp_rocky)  #p=0.09


##########NOW GLMM and HABITAT VARIABLES ######


##########KELP FORESTS ##########


### make tables to see if the variables are correlated

#Needs to eliminate site

correlation_KF <- habitat_kelp[-c(1)]


cor(correlation_KF) #stipes and diameter highly correlated - run 2 models with variables separated

## Join data

glmm_KF <- left_join(habitat_kelp, fd_values, by = c("transect"="X"))


## GLMM models 


S_stipe <- glmmTMB (sp_richn ~ RI + Substrate + Stipes + 
          Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_stipe) 
Anova(S_stipe)#No significant variables


S_stipe1 <- glmmTMB (sp_richn ~ RI + Stipes 
  + Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_stipe1) 
Anova(S_stipe1)#No significant variables


S_stipe2 <- glmmTMB (sp_richn ~  Substrate + Stipes + 
      Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_stipe2) #No significant variables
Anova(S_stipe2)


S_stipe3 <- glmmTMB (sp_richn ~  Stipes + Ind_kelp+ 
                       (1|transect),data=glmm_KF, family = poisson())

summary (S_stipe3) #No significant variables
Anova(S_stipe3)


#check assumptions

qqnorm(residuals(S_stipe))
simulationOutput <- simulateResiduals(fittedModel =S_stipe, n = 250)
fittedModel=S_stipe
plot(simulationOutput)
testDispersion(simulationOutput)


##check-deviation significant

S_diameter <- glmmTMB (sp_richn ~ RI + Substrate + 
      Diameter + Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_diameter) #No significant variables
Anova(S_diameter)


S_diameter1 <- glmmTMB (sp_richn ~ Substrate + Diameter + 
                Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_diameter1) #No significant variables

Anova(S_diameter1)

S_diameter2 <- glmmTMB (sp_richn ~ RI + Diameter + 
                          Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_diameter2) #No significant variables

Anova(S_diameter2)

S_diameter3 <- glmmTMB (sp_richn ~ Diameter + 
          Ind_kelp+ (1|transect),data=glmm_KF, family = poisson())

summary (S_diameter3) #No significant variables

Anova(S_diameter3)


S_diameter4 <- glmmTMB (sp_richn ~ Ind_kelp+ 
                (1|transect),data=glmm_KF, family = poisson())

summary (S_diameter4) #No significant variables

Anova(S_diameter4)


#check assumptions

qqnorm(residuals(S_diameter))

simulationOutput <- simulateResiduals(fittedModel =S_diameter, n = 250)
fittedModel=S_diameter
plot(simulationOutput)
testDispersion(simulationOutput)

##check-deviation significant


###########FRic###############


mod1_stipe <- glmmTMB (fric ~ RI + Stipes + Substrate + 
  Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod1_stipe)
Anova(mod1_stipe)

########Ind kelp p=0.02, Stipes=0.001

#check assumptions

qqnorm(residuals(mod1_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod1_stipe, n = 250)
fittedModel=mod1_stipe
plot(simulationOutput)
testDispersion(simulationOutput) #Good


mod2_diameter <- glmmTMB (fric ~ RI + Substrate + Diameter + 
            Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod2_diameter)

Anova(mod2_diameter)

########Ind kelp p=0.03, Diameter= 0.03

#check assumptions

qqnorm(residuals(mod2_diameter))


simulationOutput <- simulateResiduals(fittedModel =mod2_diameter, n = 250)
fittedModel=mod2_diameter
plot(simulationOutput)
testDispersion(simulationOutput) #Good


##################FEve#########################



mod_feve_stipe <- glmmTMB (feve ~ RI + Stipes + Substrate + 
            Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod_feve_stipe)
Anova(mod_feve_stipe)

# No significance

#check assumptions

qqnorm(residuals(mod_feve_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod_feve_stipe, n = 250)
fittedModel=mod_feve_stipe
plot(simulationOutput) ######## Not good - check
testDispersion(simulationOutput)# This one - good

mod_feve_diameter <- glmmTMB (feve ~ RI + Substrate + Diameter + 
        Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod_feve_diameter) #No significant

#check assumptions

qqnorm(residuals(mod_feve_diameter))

simulationOutput <- simulateResiduals(fittedModel =mod_feve_diameter, n = 250)
fittedModel=mod_feve_diameter
plot(simulationOutput) # This one, not good
testDispersion(simulationOutput) #This one, good



##################Fdiv#########################



mod_fdiv_stipe <- glmmTMB (fdiv ~ RI + Stipes + Substrate + 
              Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod_fdiv_stipe) #No significant 
Anova(mod_fdiv_stipe)

#check assumptions

qqnorm(residuals(mod_fdiv_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod_fdiv_stipe, n = 250)
fittedModel=mod_fdiv_stipe
plot(simulationOutput)
testDispersion(simulationOutput)#Good


mod_fdiv_diameter <- glmmTMB (fdiv ~ RI + Substrate + Diameter + 
                    Ind_kelp+ (1|transect),data=glmm_KF, family = beta_family())

summary (mod_fdiv_diameter) #No significant variables
Anova(mod_fdiv_diameter)

#check assumptions

qqnorm(residuals(mod_fdiv_diameter))

simulationOutput <- simulateResiduals(fittedModel =mod_fdiv_diameter, n = 250)
fittedModel=mod_fdiv_diameter
plot(simulationOutput)
testDispersion(simulationOutput)
#Good


#######################ROCKY REEFS ####################

correlation_RR <- habitat_RR[-c(1)]



cor(correlation_RR) #no high correlations

## Join data

glmm_RR <- left_join(habitat_RR, fd_values, by = c("transect"="X"))

## Because it's all rock (substrate) we are not going to add this one 
# to the models

## GLMM models 


mod1 <- glmmTMB (sp_richn ~ Depth_caves + RI + (1|transect),
                 data=glmm_RR, family = poisson())

summary (mod1) #No significant variables
Anova(mod1)

mod3 <- glmmTMB (sp_richn ~ Depth_caves + (1|transect),data=glmm_RR, 
                 family = poisson())

summary (mod3)#No significant variables


mod4 <- glmmTMB (sp_richn ~ RI + (1|transect),data=glmm_RR, 
                 family = poisson)

summary (mod4)#No significant variables


#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1, n = 250)
fittedModel=mod1
plot(simulationOutput)
testDispersion(simulationOutput)

##check-deviation significant

##########FRic###################

mod1_fric <- glmmTMB (fric ~ Depth_caves + RI + 
                        (1|transect),data=glmm_RR, family = beta_family())

summary(mod1_fric)#No significant variables

mod3_fric <- glmmTMB (fric ~ Depth_caves +(1|transect),
                      data=glmm_RR, family = beta_family())

summary (mod3_fric)#No significant variables


mod4_fric <- glmmTMB (fric ~ RI + (1|transect),data=glmm_RR, 
                      family = beta_family())

summary (mod4_fric)#No significant variables


#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1_fric, n = 250)
fittedModel=mod1_fric
plot(simulationOutput)#Not good
testDispersion(simulationOutput)#Good


##########Feve###################

mod1_feve <- glmmTMB (feve ~ Depth_caves + RI + (1|transect),
                      data=glmm_RR, family = beta_family())

summary(mod1_feve)#No significant variables

mod3_feve <- glmmTMB (feve ~ Depth_caves +(1|transect),
                      data=glmm_RR, family = beta_family())

summary (mod3_feve)#No significant variables


mod4_feve <- glmmTMB (feve ~ RI + (1|transect),
                      data=glmm_RR, family = beta_family())

summary (mod4_feve)#No significant variables

#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1_feve, n = 250)
fittedModel=mod1_feve
plot(simulationOutput)
testDispersion(simulationOutput)# Good


##########Fdiv###################

mod1_fdiv <- glmmTMB (fdiv ~ Depth_caves + RI + (1|transect),
                      data=glmm_RR, family = beta_family())

summary(mod1_fdiv)#No significant variables

mod3_fdiv <- glmmTMB (fdiv ~ Depth_caves +(1|transect),
                      data=glmm_RR, family = beta_family())

summary (mod3_fdiv)#No significant variables


mod4_fdiv <- glmmTMB (fdiv ~ RI + (1|transect),
                      data=glmm_RR, family = beta_family())

summary (mod4_fdiv)#No significant variables

#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1_fdiv, n = 250)
fittedModel=mod1_fdiv
plot(simulationOutput)#Not good
testDispersion(simulationOutput)# Good


################ end of code ######################