################################################################################
##
## Script for FD analysis 
##
## Paper: Fish functional diversity is modulated by 
## small-scale habitat complexity in a temperate ecosystem 
## (to submit in Hydrobiologia)
##
## Authors: Sgarlatta, M. Paula, Ramírez-Valdez, Arturo,
## Gómez-Gómez, Antonio, Ladah, Lydia B., 
## Calderon-Aguilera, Luis E.
## 
## Code by Paula Sgarlatta 
##
##
##  Scripts to run:
##
## 1/ Statistical analysis
## 
###############################################################################

#rm(list=ls()) # cleaning memory

#load packages

library(ggpubr)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(dplyr)
library(DHARMa)

##### t test between habitats (RR and KF) #####

##Species richness

#Normality - QQ plots
ggqqplot(fd_kelp$sp_richn)   
ggqqplot(fd_RR$sp_richn)    

t.test(fd_kelp$sp_richn,fd_RR$sp_richn, var.equal=TRUE)

#No differences between habitats - p=0.45

kelp <- data.frame(KF=fd_kelp$sp_richn)

row.names(kelp)<-row.names(rep(c("K"), each=35))

two_habitats <- as.matrix(c(fd_kelp$sp_richn,fd_RR$sp_richn))

## This is not working. I'll do it manually-try again to do it in R



##Functional richness

#Normality
ggqqplot(fd_kelp$fric)  
ggqqplot(fd_RR$fric)    

t.test(fd_kelp$fric,fd_RR$fric, var.equal=TRUE)

#No differences between habitats - p=0.387

kelp <- data.frame(KF=fd_kelp$fric)

row.names(kelp)<-row.names(rep(c("K"), each=35))

two_habitats <- as.matrix(c(fd_kelp$fric,fd_RR$fric))

## This is not working. I'll do it manually-try again to do it in R


##Functional evenness

#Normality
ggqqplot(fd_kelp$feve)   
ggqqplot(fd_RR$feve)     

t.test(fd_kelp$feve,fd_RR$feve, var.equal=TRUE)

#Differences between habitats - p=0.002

##Functional divergence

#Normality
ggqqplot(fd_kelp$fdiv)  ## only few points outside 
ggqqplot(fd_RR$fdiv)     ## only few points outside 

t.test(fd_kelp$fdiv,fd_RR$fdiv, var.equal=TRUE)

### Differences between habitats - p=0.002


##### ANOVA between sites (within each habitat) #####

# Adding transect to first column

kelp <- tibble::rownames_to_column(fd_kelp, "transect")
RR <- tibble::rownames_to_column(fd_RR, "transect")

# Adding site to data

kelp$site <- substr(kelp$transect,1,2)
RR$site <- substr(RR$transect,1,2)

############################KELP FORESTS #############################


## S

s_KF <- aov(sp_richn ~ site, data=kelp)

summary(s_KF)

#Not significant - p=0.614

#Normality

ggqqplot(s_KF$residuals)

#Homogeneity of variances
bartlett.test(sp_richn ~ site, data=kelp)  #p=0.717


## Fric

fric_KF <- aov(fric ~ site, data=kelp)

summary(fric_KF)

###Significant - p=0.014

TukeyHSD(fric_KF) #RB-LR p=0.008

#Normality
ggqqplot(fric_KF$residuals)

#Homogeneity of variances
bartlett.test(fric ~ site, data=kelp)  #p=0.49


## FEve

feve_KF <- aov(feve ~ site, data=kelp)

summary(feve_KF)

###Not significant - p=0.1

#Normality
ggqqplot(feve_KF$residuals)

#Homogeneity of variances
bartlett.test(feve ~ site, data=kelp)  #p=0.08


## FDiv


fdiv_KF <- aov(fdiv ~ site, data=kelp)

summary(fdiv_KF)

###Not significant - p=0.5

#Normality
ggqqplot(fdiv_KF$residuals)

#Homogeneity of variances
bartlett.test(fdiv ~ site, data=kelp)  #p=0.67

####################ROCKY REEFS ######################

## S

s_RR <- aov(sp_richn ~ site, data=RR)
summary(s_RR)

#Not significant - p=0.3

#Normality

ggqqplot(s_RR$residuals)

#Homogeneity of variances
bartlett.test(sp_richn ~ site, data=RR)  #p=0.574

## FRic

fric_RR <- aov(fric ~ site, data=RR)
summary(fric_RR)

#Not significant - p=0.2

#Normality

ggqqplot(fric_RR$residuals)


#Homogeneity of variances
bartlett.test(fric ~ site, data=RR)  #p=0.15


##########NOW GLMM and HABITAT VARIABLES ######


##########KELP FORESTS ##########


#loading data

habitat_KF <- read.csv("data/habitat_kelp.csv")



### make tables to see if the variables are correlated

#Needs to eliminate site

correlation_KF <- habitat_KF[-c(1)]


cor(correlation_KF) #stipes and diameter highly correlated - run 2 models with variables separated

## Join data

glmm_KF <- left_join(habitat_KF, kelp, by = "transect" )


## GLMM models 


S_stipe <- glmmTMB (sp_richn ~ RI + Substrate + Stipes + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_stipe) #No significant variables


S_stipe1 <- glmmTMB (sp_richn ~ RI + Stipes + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_stipe1) #No significant variables


S_stipe2 <- glmmTMB (sp_richn ~  Substrate + Stipes + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_stipe2) #No significant variables


S_stipe3 <- glmmTMB (sp_richn ~  Stipes + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_stipe3) #No significant variables



#check assumptions

qqnorm(residuals(S_stipe))
simulationOutput <- simulateResiduals(fittedModel =S_stipe, n = 250)
fittedModel=S_stipe
plot(simulationOutput)
testDispersion(simulationOutput)


##check-deviation significant

S_diameter <- glmmTMB (sp_richn ~ RI + Substrate + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_diameter) #No significant variables


S_diameter1 <- glmmTMB (sp_richn ~ Substrate + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_diameter1) #No significant variables

S_diameter2 <- glmmTMB (sp_richn ~ RI + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_diameter2) #No significant variables

S_diameter3 <- glmmTMB (sp_richn ~ Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_diameter3) #No significant variables


S_diameter4 <- glmmTMB (sp_richn ~ Ind_kelp+ (1|site),data=glmm_KF, family = poisson())

summary (S_diameter4) #No significant variables


#check assumptions

qqnorm(residuals(S_diameter))

simulationOutput <- simulateResiduals(fittedModel =S_diameter, n = 250)
fittedModel=S_diameter
plot(simulationOutput)
testDispersion(simulationOutput)

##check-deviation significant


###########FRic###############


mod1_stipe <- glmmTMB (fric ~ RI + Stipes + Substrate + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod1_stipe)

########Ind kelp p=0.0016, RI=0.005

#check assumptions

qqnorm(residuals(mod1_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod1_stipe, n = 250)
fittedModel=mod1_stipe
plot(simulationOutput)
testDispersion(simulationOutput)


mod2_diameter <- glmmTMB (fric ~ RI + Substrate + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod2_diameter)

########Ind kelp p=0.001, RI= 0.005

#check assumptions

qqnorm(residuals(mod2_diameter))


simulationOutput <- simulateResiduals(fittedModel =mod2_diameter, n = 250)
fittedModel=mod2_diameter
plot(simulationOutput)
testDispersion(simulationOutput)


##################FEve#########################



mod_feve_stipe <- glmmTMB (feve ~ RI + Stipes + Substrate + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod_feve_stipe)

# RI= 0.02

#check assumptions

qqnorm(residuals(mod_feve_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod_feve_stipe, n = 250)
fittedModel=mod_feve_stipe
plot(simulationOutput)


mod_feve_diameter <- glmmTMB (feve ~ RI + Substrate + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod_feve_diameter) #RI=0.04

#check assumptions

qqnorm(residuals(mod_feve_diameter))

simulationOutput <- simulateResiduals(fittedModel =mod_feve_diameter, n = 250)
fittedModel=mod_feve_diameter
plot(simulationOutput)
testDispersion(simulationOutput)



##################Fdiv#########################



mod_fdiv_stipe <- glmmTMB (fdiv ~ RI + Stipes + Substrate + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod_fdiv_stipe) #No significant 

#check assumptions

qqnorm(residuals(mod_fdiv_stipe))

simulationOutput <- simulateResiduals(fittedModel =mod_fdiv_stipe, n = 250)
fittedModel=mod_fdiv_stipe
plot(simulationOutput)
testDispersion(simulationOutput)


mod_fdiv_diameter <- glmmTMB (fdiv ~ RI + Substrate + Diameter + Ind_kelp+ (1|site),data=glmm_KF, family = beta_family())

summary (mod_fdiv_diameter) #No significant variables

#check assumptions

qqnorm(residuals(mod_fdiv_diameter))

simulationOutput <- simulateResiduals(fittedModel =mod_fdiv_diameter, n = 250)
fittedModel=mod_fdiv_diameter
plot(simulationOutput)
testDispersion(simulationOutput)



#######################ROCKY REEFS ####################

habitat_RR <- read.csv("data/habitat_RR.csv")



correlation_RR <- habitat_RR[-c(1)]



cor(correlation_RR) #no high correlations



## Join data

glmm_RR <- left_join(habitat_RR, RR, by = "transect" )


## GLMM models 


mod1 <- glmmTMB (sp_richn ~ Depth_caves + RI + Substrate + (1|site),data=glmm_RR, family = poisson())

summary (mod1) #No significant variables

mod2 <- glmmTMB (sp_richn ~ Depth_caves + RI + (1|site),data=glmm_RR, family = poisson())

summary (mod2)#No significant variables

mod3 <- glmmTMB (sp_richn ~ Depth_caves +Substrate+ (1|site),data=glmm_RR, family = poisson())

summary (mod3)#No significant variables

mod4 <- glmmTMB (sp_richn ~ RI + Substrate + (1|site),data=glmm_RR, family = poisson)

summary (mod4)


#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1, n = 250)
fittedModel=mod1
plot(simulationOutput)
testDispersion(simulationOutput)

##check-deviation significant

##########FRic###################

mod1_fric <- glmmTMB (fric ~ Depth_caves + RI + Substrate + (1|site),data=glmm_RR, family = beta_family())

summary(mod1_fric)#No significant variables

##check warnings?

mod2_fric <- glmmTMB (fric ~ Depth_caves + RI + (1|site),data=glmm_RR, family = beta_family())

summary (mod2_fric)#No significant variables

mod3_fric <- glmmTMB (fric ~ Depth_caves +Substrate+ (1|site),data=glmm_RR, family = beta_family())

summary (mod3_fric)#No significant variables


mod4_fric <- glmmTMB (fric ~ RI + Substrate + (1|site),data=glmm_RR, family = beta_family())

summary (mod4_fric)#No significant variables


#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1_fric, n = 250)
fittedModel=mod1_fric
plot(simulationOutput)
testDispersion(simulationOutput)


##########Feve###################

mod1_feve <- glmmTMB (feve ~ Depth_caves + RI + Substrate + (1|site),data=glmm_RR, family = beta_family())

summary(mod1_feve)#No significant variables

mod2_feve <- glmmTMB (feve ~ Depth_caves + RI + (1|site),data=glmm_RR, family = beta_family())

summary (mod2_feve)#No significant variables

mod3_feve <- glmmTMB (feve ~ Depth_caves +Substrate+ (1|site),data=glmm_RR, family = beta_family())

summary (mod3_feve)#No significant variables


mod4_feve <- glmmTMB (feve ~ RI + Substrate + (1|site),data=glmm_RR, family = beta_family())

summary (mod4_feve)#No significant variables

#check assumptions

simulationOutput <- simulateResiduals(fittedModel =mod1_feve, n = 250)
fittedModel=mod1_feve
plot(simulationOutput)
testDispersion(simulationOutput)


################ end of code ######################