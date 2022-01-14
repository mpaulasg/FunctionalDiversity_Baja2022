################################################################################
##
## Script for FD analysis 
##
## Paper: Fish functional diversity is modulated by 
## small-scale habitat complexity in a temperate ecosystem 
## (to submit in Hydrobiologia)
##
## Authors: Sgarlatta, M. Paula, Ramirez-Valdez, Arturo,
## Gomez-Gomez, Antonio, Ladah, Lydia B., 
## Calderon-Aguilera, Luis E.
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
#library(remotes)
#remotes::install_github("CmlMagneville/mFD")
library(mFD)
library(ggplot2)


#### 1 - Load (and transform) datasets ####

# load traits data - Kelp and rocky reefs:
traits <- read.csv("data/Baja_traits_fish.csv")

# transform columns type:

#size as numeric
traits$maximum_length <- as.numeric(traits$maximum_length)

# diet as factor

traits$Diet <- as.factor(traits$Diet) #I think this is the correct way

#traits$Diet <- factor(traits[,"Diet"], levels=c("Herbivorous", "Invertivores/zooplanktivores", 
      #"Macrocarnivores", "Mobile benthic invertivores", "Zooplanktivores" ), ordered = TRUE)

#habitat as ordinal

traits$Habitat <-factor(traits[,"Habitat"], levels=c("Bottom", "Slightly-bottom","In fronds", "Pelagic", 
                                                     "Roving" ), ordered = TRUE )
#gregariousness as ordinal

traits$Gregariousness <- as.ordered(traits$Gregariousness)
ordered (traits$Gregariousness)

# Change the order of the groups with Solitary < In Pairs or 
#sometimes aggregating < Forming Schools

traits[, "Gregariousness"] <- factor(traits[, "Gregariousness"],
                                          levels = c("Solitary", "In pairs or sometimes aggregating","Forming schools"), ordered = TRUE)

#morphology as ordinal

traits$Morphology <- factor(traits[,"Morphology"], levels=c("COMA", "COMB", "ANG", "DEP"), ordered = TRUE ) ###Need to change this to English

# add species as row names:

rownames(traits) <- traits$species_id
traits <- tibble::column_to_rownames(traits, var = "Species")

# load fish composition data - both kelp forests 
#and rocky reefs

fish <- read.csv("data/fish_kelp&rocky_reef.csv")

# group between habitats:

fish$hab <- rep(c("Kelp", "Reef"), each=36) 

# create a new dataframe for habitat (and not transects):

fish_habitat <- as.data.frame(matrix(ncol = ncol(fish), nrow = 2))
colnames(fish_habitat) <- colnames(fish)
fish_habitat$site_id <- c("Kelp", "Reef")

rownames(fish_habitat) <- fish_habitat$site_id

# fill this new df:
for (c in colnames(fish[, - c(1, ncol(fish))])) {
  hab_kelp <- mean(fish[which(fish$hab == "Kelp"), c])
  fish_habitat["Kelp", c] <- hab_kelp
  hab_rocky <- mean(fish[which(fish$hab == "Reef"), c])
  fish_habitat["Reef", c] <- hab_rocky
}
# remove Site, site id and habitat:

fish_habitat <- fish_habitat[, - c(1,42,43, ncol(fish_habitat))]

# load trait types data:
traits_type <- read.csv("data/traits_type.csv")

################################################################################



#### 2 - Basic FD analysis ####


## Species traits summary:

fish_traits_summ <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                 sp_tr  = traits)
# repartition of traits:

fish_traits_summ$tr_summary_list


## summary of the assemblages * species data frame:

fish_habitat <- as.matrix(fish_habitat)

fish_traits_summ <- mFD::asb.sp.summary(asb_sp_w = fish_habitat) 

# retrieve occurrence matrix:
occ_matrix <- fish_traits_summ$asb_sp_occ

# total biomass of each species:
fish_traits_summ$sp_tot_w

# total biomass for each asb (ie habitat):
fish_traits_summ$asb_tot_w

# species richness per habitat:
fish_traits_summ$asb_sp_richn

# list of species present in each assemblage (ie habitat):
fish_traits_summ$asb_sp_nm


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

# We can see that the functional space with the lowest quality metric is the one with 4 dimensions (pcoa_4d=0.42),
# Will use 3D to build the functional space

# illustrate quality of functional space (for Supp Material)

fs_qual <- mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", "pcoa_4d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

fs_qual

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
traits_axes   # Will need to change the labels 

# export figure

ggsave(plot = traits_axes, filename  = here::here("figures", "figure2_SM.png"), 
       height = 12, width = 25, unit = "cm")

################################################################################

#### 3 - Compute functional alpha diversity indices ####

### Compute indices:
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = fish_habitat,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)


fd_values <- alpha_fd_indices$functional_diversity_indices


###############################################################

#### 4 - Graph functional space ####

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("Kelp", "Reef"),
  ind_nm                   = c("fric"),
  faxes                    = paste0("PC", 1:3),
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


###############################################################

#### 5 - FD indices per site and transect for 
####              statistical analysis   

#### Load (and transform) datasets for kelp ####

rm(list=ls()) # cleaning memory

#load packages

#library(remotes)
#remotes::install_github("CmlMagneville/mFD")
library(mFD)
library(ggplot2)

traits_kelp <- read.csv("data/kelp_traits.csv")

# transform columns type:

#size as numeric:

traits_kelp$maximum_length <- as.numeric(traits_kelp$maximum_length)


#diet as factor

traits_kelp$Diet <- as.factor(traits_kelp$Diet)

#traits_kelp$Trophic.groups <- factor(traits_kelp[,"Trophic.groups"], levels=c("Herbivorous", "Invertivores/zooplanktivores", 
                            #"Macrocarnivores", "Mobile benthic invertivores", "Zooplanktivores" ), ordered = TRUE)

#habitat as ordinal

traits_kelp$Habitat <-factor(traits_kelp[,"Habitat"], levels=c("Bottom", "Slightly-bottom", "In fronds", 
                      "Pelagic", "Roving" ), ordered = TRUE )

#gregariousness as ordinal

traits_kelp$Gregariousness <- as.ordered(traits_kelp$Gregariousness)
ordered (traits_kelp$Gregariousness)

# Change the order of the groups with Solitary < In Pairs or 
#sometimes aggregating < Forming Schools

traits_kelp[, "Gregariousness"] <- factor(traits_kelp[, "Gregariousness"],
                                     levels = c("Solitary", "In pairs or sometimes aggregating","Forming schools"), ordered = TRUE)

#morphology as ordinal:

traits_kelp$Morphology <- factor(traits_kelp[,"Morphology"], levels=c("COMA", "COMB", "DEP"), ordered = TRUE )

# add species as row names:

rownames(traits_kelp) <- traits_kelp$species_id
traits_kelp <- tibble::column_to_rownames(traits_kelp, var = "Species")

# load fish composition data for kelp forests 

fish_kelp <- read.csv("data/fish_kelp.csv")

# add transects as row names:
rownames(fish_kelp) <- fish_kelp$species_id
fish_kelp <- tibble::column_to_rownames(fish_kelp, var = "Site")

# load trait types data:
traits_type <- read.csv("data/traits_type.csv")

###### Load (and transform) datasets for rocky reef ######

traits_RR <- read.csv("data/RR_traits.csv")

# transform columns type:

#size as numerical
traits_RR$maximum_length <- as.numeric(traits_RR$maximum_length)

#Diet as factor
traits_RR$Diet <- as.factor(traits_RR$Diet)

#habitat as ordinal

traits_RR$Habitat <-factor(traits_RR[,"Habitat"], levels=c("Bottom", "Slightly-bottom", 
                                                               "Pelagic", "Roving" ), ordered = TRUE )
#gregariousness as ordinal

traits_RR$Gregariousness <- as.ordered(traits_RR$Gregariousness)

# Change the order of the groups with Solitary < In Pairs or 
#sometimes aggregating < Forming Schools

traits_RR[, "Gregariousness"] <- factor(traits_RR[, "Gregariousness"],
        levels = c("Solitary", "In pairs or sometimes aggregating","Forming schools"), ordered = TRUE)

#morphology as ordinal 

traits_RR$Morphology <- factor(traits_RR[,"Morphology"], levels=c("COMA", "COMB", "ANG"), ordered = TRUE )

# add species as row names:

rownames(traits_RR) <- traits_RR$species_id
traits_RR <- tibble::column_to_rownames(traits_RR, var = "Species")

# load fish composition data for rocky reefs 

fish_RR <- read.csv("data/fish_RR.csv")

# add transects as row names:
rownames(fish_RR) <- fish_RR$species_id
fish_RR <- tibble::column_to_rownames(fish_RR, var = "Site")

################################################################################



#### 6 - Basic FD analysis for kelp forest and RR separated ####


## Species traits summary:


kelp_traits_summ <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                       sp_tr  = traits_kelp)



RR_traits_summ <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                       sp_tr  = traits_RR)

# Distribution of traits:

kelp_traits_summ$tr_summary_list


RR_traits_summ$tr_summary_list

## summary of the assemblages * species data frame:

fish_kelp <- as.matrix(fish_kelp)

fish_RR <- as.matrix(fish_RR)

kelp_traits_summ <- mFD::asb.sp.summary(asb_sp_w = fish_kelp) 

RR_traits_summ <- mFD::asb.sp.summary(asb_sp_w = fish_RR)

# retrieve occurrence matrix:

occ_matrix_kelp <- kelp_traits_summ$asb_sp_occ

occ_matrix_RR <- RR_traits_summ$asb_sp_occ

# total biomass of each species:

kelp_traits_summ$sp_tot_w

RR_traits_summ$sp_tot_w

# total biomass for each transect:

kelp_traits_summ$asb_tot_w

RR_traits_summ$asb_tot_w

# species richness per transect:

kelp_traits_summ$asb_sp_richn

RR_traits_summ$asb_sp_richn


### Compute functional distance:

# compute functional distance based on traits:
sp_dist_kelp <- mFD::funct.dist(
  sp_tr  = traits_kelp,
  tr_cat = traits_type,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)

sp_dist_RR <- mFD::funct.dist(
  sp_tr  = traits_RR,
  tr_cat = traits_type,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)

### Compute functional spaces and their quality:

# compute the quality of functional spaces:
fspaces_quality_kelp <- mFD::quality.fspaces(
  sp_dist             = sp_dist_kelp, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")

fspaces_quality_RR <- mFD::quality.fspaces(
  sp_dist             = sp_dist_RR, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")

# retrieve table of quality metric:
fspaces_quality_kelp$quality_fspaces

# We can see that the functional space with the lowest quality metric is the one with three dimensions (pcoa_3d=0.43)
# Will use 3D to build the functional space

fspaces_quality_RR$quality_fspaces

# We can see that the functional space with the lowest quality metric is the one with three dimensions (pcoa_3d=0.40)
# Will use 3D to build the functional space


### Test correlation between traits and functional axes 

# retrieve coordinates of species:
sp_faxes_coord_kelp <- fspaces_quality_kelp$"details_fspaces"$"sp_pc_coord"

sp_faxes_coord_RR <- fspaces_quality_RR$"details_fspaces"$"sp_pc_coord"

# test correlation between traits and axes:
cor_tr_faxes_kelp <- mFD::traits.faxes.cor(
  sp_tr          = traits_kelp, 
  sp_faxes_coord = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

cor_tr_faxes_RR <- mFD::traits.faxes.cor(
  sp_tr          = traits_RR, 
  sp_faxes_coord = sp_faxes_coord_RR[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# get the table of correlation:
cor_tr_faxes_kelp$tr_faxes_stat

cor_tr_faxes_RR$tr_faxes_stat


# get the plot:
traits_axes_kelp <- cor_tr_faxes_kelp$tr_faxes_plot

traits_axes_RR <- cor_tr_faxes_RR$tr_faxes_plot


################################################################################

#### 7 - Compute functional alpha diversity indices (for
#### kelp forest and RR separetely) ####

### Compute indices:

  fd_indices_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = fish_kelp,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

## Need to eliminate some transects as they have less than 3 sp

fish_kelp_2 <- fish_kelp[-c(20), ]

fish_kelp_2 <- as.matrix(fish_kelp_2)

fd_indices_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = fish_kelp_2,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_kelp <- fd_indices_kelp$functional_diversity_indices

write.csv(fd_kelp, file="data/fd_kelp.csv")

fd_indices_RR <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_RR[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = fish_RR,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

## Need to eliminate some transects as they have less than 3 sp

fish_RR_2 <- fish_RR[-c(33), ]

fish_RR_2 <- as.matrix(fish_RR_2)

fd_indices_RR <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_RR[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = fish_RR_2,
  ind_vect         = c("fric", "feve", "fdiv"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_RR <- fd_indices_RR$functional_diversity_indices

write.csv(fd_RR, file="data/fd_RR.csv")

################## end of code ##################
