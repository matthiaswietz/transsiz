
############################################################################################
           ###  RAS2016-17 - AMPLICON ANALYSIS  ###
############################################################################################

# This script: calculate rarefaction and diversity indices 

setwd("/AWI_MPI/collaborations/TRANSSIZ/Rstats")
load("TRANSSIZ.Rdata")

# Load packages
library(iNEXT)
library(phyloseq)
library(olsrr)
library(ggplot2)
library(cowplot)
library(dplyr)


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

## CALCULATION ##

iNEXT <- otu_table(
  ASV, taxa_are_rows=F)

iNEXT <- iNEXT(
  as.data.frame(otu_table(iNEXT)), q=c(0),
  datatype="abundance", conf=0.95, nboot=100)

###################################

## RAREFACTION ##

rarefac <- fortify(iNEXT, type=1) %>%
  left_join(ENV, by=c("site"="sample_title"))

rarefac.point <- rarefac[which(
  rarefac$method == "observed"),]
rarefac.line <- rarefac[which(
  rarefac$method != "observed"),]
rarefac.line$method <- factor(rarefac.line$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

###################################

## COVERAGE ##

cover <- fortify(iNEXT, type=2) %>%
  left_join(ENV, by=c("site"="sample_title"))

cover.point <- cover [which(
  cover$method == "observed"),]
cover.line <- cover [which(
  cover$method != "observed"),]
cover.line$method <- factor(cover.line$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

coverage <- ggplot(cover, 
  aes(x=x, y=y, colour=site))+ 
geom_line(aes(linetype = method), 
  lwd = 0.5, data = cover.line) +
#scale_colour_discrete(guide="none") +
  scale_color_manual(
    values = scico(34, palette='lajolla')) +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
labs(x="Sample size", y="Sample coverage") +
theme_bw(base_size=12) + 
theme(legend.position="none",
      axis.ticks = element_blank())

###################################

## RICHNESS ##

richness <- ggplot(rarefac, 
  aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = rarefac.line) +
 # scale_colour_discrete(guide="none") +
 # scale_colour_manual(values = rev(
  #  pnw_palette("Sunset", 34))) +
  scale_color_manual(
    values = scico(34, palette='lajolla')) +
  scale_x_continuous(limits = c(0,1e+5)) +
  labs(x="Sample size", y="Species richness") +
  theme_bw(base_size=12) + 
  theme(legend.position="none",
        axis.ticks = element_blank())

# Plot curves
plot_grid(
  richness, 
  coverage,
  ncol=2,
  align="tblr")

###################################

## CREATE SUMMARIES ##

richness <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site) 
simpson <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site) 
shannon <- estimate_richness(
  pseq.abs, measures=c("Shannon")) %>%
  rownames_to_column("sample_title") %>%
  arrange(sample_title)

###################################

# Merge all

AlphaDiv <- data.frame(
  sample_title = ENV$sample_title,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Shannon) %>%
  mutate(evenness = shannon/log(richness)) %>%
  left_join(ENV)

#################

# remove temp-data
rm(richness, simpson, shannon,
   rarefac, rarefac.point, rarefac.line,
   cover, cover.point, cover.line)



############################################################################################
   ###  Distance-Decay  ### 
############################################################################################

## Distances for PERMANOVA + BetaDiv ##

###################################

# Merge *div*
div.bac <- merge(d.bac, div.bac, by="ID_month") 
div.euk <- merge(d.euk, div.euk, by="ID_month") 

# Combine everything
div.all <- rbind(div.bac, div.euk) 

# Append ENV
div.all <- merge(div.all, ENV.euk[, c(
   "ID_month","ID_plot","month","temp","nutrients",
   "daylight","ice_act","ice_past","AW","PW","pH",
   "CO2","O2_conc","O2_sat","strat","season")], 
   by="ID_month", all.y = F) 

###################################
