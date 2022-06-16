
############################################################################################
##  TRANSSIZ - bacterial communities  ##
############################################################################################

# In total 38 samples were analyzed by DADA
# Four omitted here: original lat/lon info not documented
# deposited at ENA nonetheless, for fully reproducible ASV generation
# in this workflow, the metadata file only lists the relevant samples
# Full ASV table will be hence subsetted to 34 samples

#######################

## This script: format dada2-amplicons and metadata 

# Set working directory
setwd("/AWI_MPI/collaborations/TRANSSIZ/Rstats")
load("TRANSSIZ.Rdata")

# Load packages and colors
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(tidyr)
library(psych)
library(stringr)
library(scico)
library(scales)
source("Colors.R")


############################################################################################
   ###  RAW DATA -- LOAD AND FORMATTING ###
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

# Rename BAC-NAs with last known taxrank + "uc"
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

# Full latitudinal data
VOC.full <- read.table(
  "VOCs.txt", h=T, sep="\t") %>%
  na_if("") 

# mean by "lat_short" for ASV correlations
VOC.short <- read.table(
  "VOCs.txt", h=T, sep="\t", 
  stringsAsFactors=F) %>%
  na_if("") %>%
  mutate(lat_short=lat) %>%
  mutate(across(lat_short, round, 1)) %>%
  group_by(lat_short) %>%
  summarise(across(everything(), mean)) %>%
  dplyr::select(-c("lat")) %>%
  ungroup

# MeSH-DMS details 
MeSH_DMS <- read.csv(
  "Appendix_MeSH-DMS.txt", h=T, 
  sep="\t", check.name =F) %>%
  na_if("") 

# Figure SX (Chl a)
Chl <- read.csv(
  "Appendix_chla.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

# Figure SX (Chl-DMS-Iso >80°N)
# recalculate ng to ug Chla
Chl_80N <- read.csv(
  "Appendix_chla_80N.txt", h=T, 
  sep="\t", check.name=F) %>%
  drop_na() %>%
  mutate(`Chla [ug/L]`=`Chla [ng/L]`/1000)

# Figure SX (Henry's Law)
Henry <- read.csv(
  "Appendix_HenryLaw.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

# Figure SX (Acetone)
Acetone <- read.csv(
  "Appendix_Acetone.txt", h=T, 
  sep="\t", check.name=F) %>%
  na_if("") 

############################

# merge
ENV <- ENV %>% 
  left_join(VOC.short) 

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

# remove temp-data
rm(asv,tax,i,j,k)

save.image("TRANSSIZ.Rdata")
