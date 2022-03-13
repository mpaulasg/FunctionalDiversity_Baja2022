################################################################################
##
## Script for FD analysis 
##
## Paper: Fish functional diversity is modulated by 
## small-scale habitat complexity in a temperate ecosystem 
## (to submit in Hydrobiologia)
##
## Authors: Sgarlatta, M. Paula, Ramírez-Valdez, Arturo,
## Ladah, Lydia B., Calderon-Aguilera, Luis E.
## 
## Code by Paula Sgarlatta 
##
##
##  Scripts to run:
##
## 1/ Figures from statistical analysis (script 2)
##
###############################################################################


rm(list=ls()) # cleaning memory


# load packages

library(tidyverse)
library(ggplot2)
library(here)

#Loading data

fd_both <- read.csv(here::here("data", "fd_values_replicates.csv"))

metadata <- read.csv(here::here("data", "kelp&rocky_metadata_replicate.csv"))


#Prepare data 

fd_all <- fd_both %>% 
  left_join(metadata, by=c("X"="Code")) %>% 
  select(-X, -HabitatID, -Replicate) %>% 
  relocate(Site, Habitat, .before = "sp_richn")


# color code for the 2 habitats

hab_colors <- c(Kelp = "mediumseagreen", Rocky ="sienna")

#Species richness (S)

sp_rich <- ggplot(fd_all, aes(x=Site, y=sp_richn, fill = Habitat)) +
  geom_boxplot()+
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0))+
  facet_wrap(~Habitat, scales="free_x", labeller = labeller(Habitat = 
c("Kelp" = "Kelp forests", 
  "Rocky" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  labs(x="Sites", y="Species richness") 

  mytheme <- theme(panel.background=element_rect(fill="white"), 
  panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
  panel.grid.major = element_blank(),
  axis.line = element_line(size = 1, colour = "black"),
  axis.text = element_text(size = (14),colour = "black"), 
  axis.title = element_text(size= (16)),
  legend.position = "none" , 
  strip.text.x = element_text(size = 14))

  
plot_sp_rich <- print(sp_rich + mytheme + 
scale_x_discrete(labels=c("CK_K"="CK", "LB_K"="LB", "LR_K" = "LR",
"RB_K"="RB", "LR_R"="LR", "PM_R"="PM", "RB_R"="RB", "Z_R"="Z")))


#Functional richness (Fric)

fric <- ggplot(fd_all, aes(x=Site, y=fric, fill = Habitat)) +
  geom_boxplot() +
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0))+
  facet_wrap(~Habitat, scales="free_x", labeller = labeller(Habitat = 
   c("Kelp" = "Kelp forests", "Rocky" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2))+
  labs(x="Sites", y="Functional richness") 

mytheme <- theme(panel.background=element_rect(fill="white"), 
                 panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                 panel.grid.major = element_blank(),
                 axis.line = element_line(size = 1, colour = "black"),
                 axis.text = element_text(size = (14),colour = "black"), 
                 axis.title = element_text(size= (16)),
                 legend.position = "none" , 
                 strip.text.x = element_text(size = 14))


plot_fric <- print(fric + mytheme + 
                        scale_x_discrete(labels=c("CK_K"="CK", "LB_K"="LB", "LR_K" = "LR",
                                                  "RB_K"="RB", "LR_R"="LR", "PM_R"="PM", "RB_R"="RB", "Z_R"="Z")))

#Functional evenness (Feve)

feve <- ggplot(fd_all, aes(x=Site, y=feve, fill = Habitat)) +
  geom_boxplot() +
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0))+
  facet_wrap(~Habitat, scales="free_x", labeller = labeller(Habitat = 
  c("Kelp" = "Kelp forests",  "Rocky" = "Rocky reefs"))) + 
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Sites", y="Functional evenness") 

mytheme <- theme(panel.background=element_rect(fill="white"), 
                 panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                 panel.grid.major = element_blank(),
                 axis.line = element_line(size = 1, colour = "black"),
                 axis.text = element_text(size = (14),colour = "black"), 
                 axis.title = element_text(size= (16)),
                 legend.position = "none" , 
                 strip.text.x = element_text(size = 14))


plot_feve <- print(feve + mytheme + 
scale_x_discrete(labels=c("CK_K"="CK", "LB_K"="LB", "LR_K" = "LR",
"RB_K"="RB", "LR_R"="LR", "PM_R"="PM", "RB_R"="RB", "Z_R"="Z")))

#Functional divergence (Fdiv)

fdiv <- ggplot(fd_all, aes(x=Site, y=fdiv, fill = Habitat)) +
  geom_boxplot() +
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0))+
  facet_wrap(~Habitat, scales="free_x", labeller = labeller(Habitat = 
  c("Kelp" = "Kelp forests", "Rocky" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Sites", y="Functional divergence") 

mytheme <- theme(panel.background=element_rect(fill="white"), 
                 panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                 panel.grid.major = element_blank(),
                 axis.line = element_line(size = 1, colour = "black"),
                 axis.text = element_text(size = (14),colour = "black"), 
                 axis.title = element_text(size= (16)),
                 legend.position = "none" , 
                 strip.text.x = element_text(size = 14))


plot_fdiv <- print(fdiv + mytheme + 
                     scale_x_discrete(labels=c("CK_K"="CK", "LB_K"="LB", "LR_K" = "LR",
                                               "RB_K"="RB", "LR_R"="LR", "PM_R"="PM", "RB_R"="RB", "Z_R"="Z")))


## merging all plot into a single figure and saving as png ####
figure2 <- ( plot_sp_rich + plot_fric ) / ( plot_feve +  plot_fdiv )
ggsave(figure2, file=here::here("figures/", "figure2_new.png"),
       height = 22, width = 25, unit = "cm" )
