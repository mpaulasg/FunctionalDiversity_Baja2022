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
## Code by Paula Sgarlatta - following code from Camille Magneville 
## (https://github.com/CmlMagneville/mFD) 
##
##
##  Scripts to run:
##
## 1/ Load and transform dataset - kelp and rocky reefs
## 2/ Basic FD analysis
## 3/ Alpha indices (FRic, FEve, FDiv)
## 4/ Graph of functional space (comparing habitats)
## 5/ FD indices per site and transect 
## (from each habitat separately) for statistical analysis

###############################################################################

rm(list=ls()) # cleaning memory

#load packages

library(here)
library(mFD)
library(ggplot2)
library(tidyverse)
library(dplyr)

###############################################################

#### 1 - Load (and transform) datasets - kelp and rocky reefs

traits <- read.csv(here::here("data", "Baja_traits_fish.csv"))

# transform columns type:

#size as numeric
traits$maximum_length <- as.numeric(traits$maximum_length)

# diet as ordinal

traits$diet <- as.factor(traits$diet) 

#habitat as ordinal

traits$habitat <-factor(traits[,"habitat"], levels=c("B", "S-B","FRO", "PEL", 
                                                     "ROV" ), ordered = TRUE )

# #gregariousness as ordinal
# Change the order of the groups with Solitary < In Pairs or 
#sometimes aggregating < Forming Schools

traits$gregariousness <- factor(traits[, "gregariousness"],
                    levels = c("SOL", "PAIR","SCHO"),
                    ordered = TRUE)

#morphology as ordinal

traits$morphology <- factor(traits[,"morphology"], 
                            levels=c("COMH", "COML", "EEL", "RAY"), ordered = TRUE ) 

# add species as row names:

rownames(traits) <- traits$species_id
traits <- tibble::column_to_rownames(traits, var = "species")

# load trait types data:
traits_type <- read.csv(here::here("data", "traits_type.csv"))

################################################################################



#### 2 - Basic FD analysis - for plotting functional space ####


## Species traits summary:

fish_traits_summ <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                       sp_tr  = traits)
# repartition of traits:

fish_traits_summ$tr_summary_list


### Compute functional distance:

# compute functional distance based on traits:
sp_dist <- mFD::funct.dist(
  sp_tr  = traits,
  tr_cat = traits_type,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)

# is there a lot of distances = 0? (group into FEs or not?)
sp_dist2 <- as.matrix(sp_dist)
inds <- which(sp_dist2 == min(sp_dist2), arr.ind=TRUE)
inds
# It is ok: only between Sebastes.chrysomelas and Sebastes.carnatus 
# Don't need to gather into FEs

sp_dist

### Compute functional spaces and their quality:

# compute the quality of functional spaces:
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")
# Same warning message saying that some distance equal 0 
# but have checked and it is ok

# retrieve table of quality metric:
fspaces_quality$quality_fspaces

# We can see that the functional space with the lowest quality metric 
#is the one with 4 dimensions (pcoa_4d=0.037),

# illustrate quality of functional space (for Supp Material)

fs_qual <- mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average","pcoa_1d", "pcoa_2d", "pcoa_3d", "pcoa_4d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

fs_qual

faxes_coord_habitat <-fspaces_quality$details_fspaces$sp_pc_coord

# export figure

ggsave(plot = fs_qual, file  = here::here("figures", "figure1_SM.png"), 
       height = 12, width = 25, unit = "cm")

### Test correlation between traits and functional axes 
## (also for Supp Material):

# retrieve coordinates of species:
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

# test correlation between traits and axes:
cor_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = traits, 
  sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# get the table of correlation:
cor_tr_faxes$tr_faxes_stat

# get the plot:
traits_axes <- cor_tr_faxes$tr_faxes_plot
traits_axes   

# export figure

ggsave(plot = traits_axes, filename  = here::here("figures", "figure2_SM.png"), 
       height = 12, width = 25, unit = "cm")

################################################################################

#### 3 - Compute functional alpha diversity indices

## Load fish biomass data and metadata

fish <- read.csv(here::here("data", "fish_kelp&rocky_reef.csv"))  
  
metadata   <- read.csv(here::here("data", "kelp&rocky_metadata.csv"))

#Prepare data to sum/average

fish_tosum <- fish %>% 
  left_join(metadata, by="Code") 

fish_replicates <-fish_tosum %>% 
  pivot_longer(cols= contains("."), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, HabitatID, Replicate, species) %>%
  summarize(total=sum(value)) %>% 
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Replicate, HabitatID, sep="_"), .before="Site") %>% 
  select(-Site, -HabitatID,-Replicate) %>% 
    column_to_rownames("Code") %>% 
  as.matrix()
   

### Compute indices:
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = fish_replicates,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_values_replicates <- alpha_fd_indices$functional_diversity_indices

## saving as csv

write.csv(fd_values_replicates, file=here::here("data", 
  "fd_values_replicates.csv") )


###############################################################

#### 4 - Graph functional space 
#### Only for kelp forests and rocky reef  in general

#Metadata for replicates

metadata_replicates <- read.csv(here::here("data", "kelp&rocky_metadata_replicate.csv"))

# computing biomass of species in each habitat
habitat_fric <- rbind( 
  Kelp = apply(fish_replicates[metadata_replicates[which(metadata_replicates$HabitatID=="K"),"Code"],],2,max ),
  Rocky = apply(fish_replicates[metadata_replicates[which(metadata_replicates$HabitatID=="R"),"Code"],],2,max )
) 


### Calculate fric for kelp and rocky reef

fd_habitat <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = habitat_fric,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_habitat_replicates <- fd_habitat$functional_diversity_indices

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = fd_habitat,
  plot_asb_nm              = c("Kelp", "Rocky"),
  ind_nm                   = c("fric"),
  faxes                    = paste0("PC", 1:4),
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "mediumseagreen", asb2 = "sienna"),
  color_vert               = c(pool = "grey50", asb1 = "mediumseagreen", asb2 = "sienna"),
  fill_sp                  = c(pool = NA, asb1 = "mediumseagreen", asb2 = "sienna"),
  fill_vert                = c(pool = NA, asb1 = "mediumseagreen", asb2 = "sienna"),
  color_ch                 = c(pool = NA, asb1 = "mediumseagreen", asb2 = "sienna"),
  fill_ch                  = c(pool = "white", asb1 = "mediumseagreen", asb2 = "sienna"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  #color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  #fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 


plot_fric <- plots_alpha$"fric"$"patchwork"
plot_fric

# export figure

ggsave(plot = plot_fric, file  = here::here("figures", "figure3.png"),
       height = 25, width = 25, unit = "cm")

############## end of code 