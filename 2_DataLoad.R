
############################################################################################
   ###  TRANSSIZ - BACTERIAL COMMUNITIES  ###
############################################################################################

# Published in Gros et al. (doi 10.5194/bg-2022-150) 
# In total 38 samples were analyzed by DADA 
# Four omitted here: original lat/lon info not documented
# Deposited at ENA nonetheless, for fully reproducible ASV generation
# Here, the metadata file only lists the relevant samples
# Full ASV table is hence subsetted to 34 samples

#######################

# Set working directory - adjust for your own system
# Place all input files in this directory
setwd("")

# Load packages and colors
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(gtools)
library(stringr)
library(ggpubr)
library(scales)
library(phyloseq)
library(vegan)
library(iNEXT)
library(olsrr)
library(cowplot)
library(psych)
library(scico)

############################################################################################
   ###  RAW DATA  ###
############################################################################################

# Read ASVs
ASV <- read.table(
  "seqtab.txt",
  h = T,
  sep = "\t",
  row.names = 1, 
  check.names = F)

TAX <- read.table(
  "tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

# Rename NAs with last known taxrank + "uc"
k <- ncol(TAX)-1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) >1) {
    temp <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], "_uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- temp}
  if (sum(is.na(TAX[, i]))==1) {
    temp <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], "_uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX[is.na(TAX[, i]),] <- temp
  }
}
TAX[is.na(TAX[, (k+1)]), (k+1)] <- paste(
  TAX[is.na(TAX[, (k+1)]), k], "_uc", sep="")

# shorten/modify names
TAX <- TAX %>%
  mutate(across(everything(),~gsub("Clade","SAR11_Clade", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_clade_NAC11-7_lineage","Roseobacter_NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406_clade)","SAR406", ., fixed=T)))


############################################################################################
   ###  METADATA + VOCs ###
############################################################################################

ENV <- read.table(
  "metadata.txt", h=T, sep="\t", 
  stringsAsFactors=F) %>%
  mutate(lat_short=lat) %>%
  mutate(across(lat_short, round, 1))

# VOCs
VOCs <- read.table(
  "VOCs.txt", h=T, sep="\t") %>%
  na_if("") 

# watermasses
watermass <- read.csv(
  "watermasses.txt", h=T, 
  sep="\t", check.name =F) %>%
  na_if("") %>%
 mutate(lat_short=lat) %>%
 mutate(across(lat_short, round, 3)) %>%
 group_by(lat_short, watermass) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup

# MeSH-DMS details 
MeSH_DMS <- read.csv(
  "Appendix_MeSH-DMS.txt", h=T, 
  sep="\t", check.name =F) %>%
  na_if("") %>%
  mutate(lat_short=lat) %>%
  mutate(across(lat_short, round, 2)) %>%
  group_by(lat_short) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup 

# Chlorophyll a
Chl <- read.csv(
  "Appendix_chla.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

# Chlorophyll-DMS-Iso >80°N
# recalculate ng to ug Chla
Chl_80N <- read.csv(
  "Appendix_chla_80N.txt", h=T, 
  sep="\t", check.name=F) %>%
  drop_na() %>%
  mutate(Chla=Chla/1000)

# Henry's Law
Henry <- read.csv(
  "Appendix_HenryLaw.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

# Acetone
Acetone <- read.csv(
  "Appendix_Acetone.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

############################

# reformat  
ASV <- ASV[,mixedsort(names(ASV))]
ASV <- ASV[, ENV$clip_id]

# rename
colnames(ASV) = ENV$sample_title
row.names(ENV) = ENV$sample_title


############################################################################################
   ###  PHYLOSEQ LOAD  ###
############################################################################################

asv = otu_table(ASV, taxa_are_rows=T)
tax = tax_table(TAX)
rownames(tax) <- rownames(asv)

pseq.abs <- phyloseq(otu_table(
  asv, taxa_are_rows=F), 
  sample_data(ENV), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(pseq.abs@tax_table)<- c(
  "Kingdom","Phylum","Class",
  "Order","Family","Genus")

# calculate rel. abundances
pseq.rel = transform_sample_counts(
  pseq.abs, function(x) x / sum(x) * 100) 


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

iNEXT <- otu_table(
  ASV, taxa_are_rows=F)

iNEXT <- iNEXT(
  as.data.frame(otu_table(iNEXT)), q=c(0),
  datatype="abundance", conf=0.95, nboot=100)

###################################

## RAREFACTION ##

rarefac <- fortify(iNEXT, type=1) %>%
  left_join(ENV, by=c("Assemblage"="sample_title"))

rarefac.point <- rarefac[which(
  rarefac$Method == "Observed"),]
rarefac.line <- rarefac[which(
  rarefac$Method != "Observed"),]
rarefac.line$Method <- factor(
  rarefac.line$Method, 
  c("Rarefaction","Extrapolation"),
  c("Rarefaction","Extrapolation"))

richness <- ggplot(
  rarefac, aes(x=x, y=y, colour=Assemblage)) +
geom_line(aes(
  linetype = Method), 
  linewidth = 0.5, data = rarefac.line) +
scale_color_manual(
  values = scico(34, palette='lajolla')) +
scale_x_continuous(limits = c(0,1e+5)) +
labs(
  x="Sample size", 
  y="Species richness") +
theme_bw(base_size=12) + 
theme(legend.position="none",
      axis.ticks = element_blank())

###################################

## COVERAGE ##

cover <- fortify(iNEXT, type=2) %>%
  left_join(ENV, by=c("Assemblage"="sample_title"))

cover.point <- cover [which(
  cover$Method == "Observed"),]
cover.line <- cover [which(
  cover$Method != "Observed"),]
cover.line$Method <- factor(
  cover.line$Method,
  c("Rarefaction","Extrapolation"),
  c("Rarefaction","Extrapolation"))

coverage <- ggplot(
  cover, aes(x=x, y=y, colour=Assemblage)) + 
geom_line(
  aes(linetype = Method), 
  linewidth = 0.5, data = cover.line) +
scale_color_manual(
  values = scico(34, palette='lajolla')) +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
labs(
  x="Sample size", 
  y="Sample coverage") +
theme_bw(base_size=12) + 
theme(legend.position="none",
      axis.ticks = element_blank())

###################################

# Plot all
plot_grid(
  richness, 
  coverage,
  ncol=2,
  align="tblr")
