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
## 1/ Figures from statistical analysis (script 2)
##
###############################################################################


rm(list=ls()) # cleaning memory


# load packages

library(tidyverse)
library(ggplot2)
library(here)

#Loading data

fd_both <- read.csv("data/fd_both.csv")

# mean and sd of diversity among each site for each habitat type
ind_toplot <- fd_both %>%
  group_by(habitat, site) %>%
  summarise( 
    n = n(),
    sp_mean = mean(sp_richn),
    sp_sd = sd(sp_richn),
    FRic_mean = mean(fric),
    FRic_sd = sd(fric),
    FEve_mean = mean(feve),
    FEve_sd = sd(feve),
    Fdiv_mean = mean(fdiv),
    Fdiv_sd = sd(fdiv),
  ) %>%
  mutate( sp_se = sp_sd/sqrt(n))  %>%
  mutate( FRic_se = FRic_sd/sqrt(n)) %>% 
  mutate( FEve_se = FEve_sd/sqrt(n))  %>%
  mutate( Fdiv_se = Fdiv_sd/sqrt(n))

ind_toplot

# color code for the 2 habitats

hab_colors <- c(Kelp_forests = "mediumseagreen", Rocky_reefs ="sienna")

#Species richness (S)

sp_rich <- ggplot(ind_toplot) +
  geom_bar( aes(x=site, y=sp_mean, color = habitat, fill = habitat), 
  stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=site, ymin=sp_mean-sp_se, ymax=sp_mean+sp_se), width=0.1, size=0.8, colour="black" ) +
facet_wrap(~habitat, scales="free_x", labeller = labeller(habitat = 
c("Kelp_forests" = "Kelp forests", 
  "Rocky_reefs" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,12), breaks = seq(from=0, to=12, by=5)  ) +
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

fric <- ggplot(ind_toplot) +
  geom_bar( aes(x=site, y=FRic_mean, color = habitat, fill = habitat), 
            stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=site, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.8, colour="black" ) +
  facet_wrap(~habitat, scales="free_x", labeller = labeller(habitat = 
                                                            c("Kelp_forests" = "Kelp forests", 
                                                              "Rocky_reefs" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.6), breaks = seq(from=0, to=1, by=0.2)  ) +
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

feve <- ggplot(ind_toplot) +
  geom_bar( aes(x=site, y=FEve_mean, color = habitat, fill = habitat), 
            stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=site, ymin=FEve_mean-FEve_se, ymax=FEve_mean+FEve_se), width=0.1, size=0.8, colour="black" ) +
  facet_wrap(~habitat, scales="free_x", labeller = labeller(habitat = 
                                                              c("Kelp_forests" = "Kelp forests", 
                                                                "Rocky_reefs" = "Rocky reefs")))+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.8), breaks = seq(from=0, to=1, by=0.2)  ) +
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

fdiv <- ggplot(ind_toplot) +
  geom_bar( aes(x=site, y=Fdiv_mean, color = habitat, fill = habitat), 
            stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=site, ymin=Fdiv_mean-Fdiv_se, ymax=Fdiv_mean+Fdiv_se), width=0.1, size=0.8, colour="black" ) +
  facet_wrap(~habitat, scales="free_x", labeller = labeller(habitat = 
                                                              c("Kelp_forests" = "Kelp forests", 
                                                                "Rocky_reefs" = "Rocky reefs")))+
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
ggsave(figure2, file=here::here("figures/", "figure2.png"),
       height = 22, width = 25, unit = "cm" )
