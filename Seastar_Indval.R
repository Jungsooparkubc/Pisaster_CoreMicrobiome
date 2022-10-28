library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the data.table package for better metadata manipulation
library(data.table)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)
#load the vegan library
library(vegan)
library(dplyr)
library(vegan)

#### environment settings ####
#set working directory
theme_set(theme_bw())
setwd("~/desktop/seastar_project/R_analysis/")

#########################
#########Indival#########

### INDVAL Indicator Species Analysis ###

library(dplyr)
library(vegan)
library(labdsv)
library(indicspecies)
library(phyloseq)


#write metadata from RDS file
#metadata <- as.data.frame(unclass(sample_data(project_data)))
#write.csv(metadata, file="metadata.csv", row.names=T)

### remove what I want###
#metadata <- subset(metadata, Sample_type !="water")

### from HERE after revising metadata####
metadata <- read.csv(file="metadata.csv")  ### after adding columns manually
View (metadata)
otu_table<- read.csv(file="genus_level_final.otu.csv")
View (otu_table)

### subsamples from metadata ### Pisastar vs Env in Bamfield
metadata <- subset(metadata, Sample_location == "Bamfield")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))

### subsamples from metadata ### Pisastar oral vs Env in PortMoody
metadata <- subset(metadata, Sample_location == "Port_Moody")
metadata <- subset(metadata, body_site %in% c("oral", "Env"))

### subsamples from metadata ### Pisastar vs Env in Hakai
metadata <- subset(metadata, Sample_location == "Hakai")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")

### subsamples from metadata ### Dermasterias vs Env in Hakai
metadata <- subset(metadata, Sample_location == "Hakai")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))
metadata <- subset(metadata, Type !="Pisaster")
metadata <- subset(metadata, Type !="Henrecia")

### subsamples from metadata ### Henrecia vs Env in Hakai
metadata <- subset(metadata, Sample_location == "Hakai")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))
metadata <- subset(metadata, Type !="Pisaster")
metadata <- subset(metadata, Type !="Dermasterias")

### subsamples from metadata ### Pisastar vs Env in Port_Moody
metadata <- subset(metadata, Sample_location == "Port_Moody")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))

### subsamples from metadata ### Pisastar vs Env in Vancouver
metadata <- subset(metadata, Sample_location !="Monterey_Bay")
metadata <- subset(metadata, Type =="Pisaster")
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))

### subsamples from metadata ### Pisastar*oral vs Env in Vancouver
metadata <- subset(metadata, Sample_location !="Monterey_Bay")
metadata <- subset(metadata, body_site %in% c("oral", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")

### subsamples from metadata ### Pisastar*cecum vs Env in Vancouver
metadata <- subset(metadata, Sample_location !="Monterey_Bay")
metadata <- subset(metadata, body_site %in% c("cecum", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")

### subsamples from metadata ### Pisastar*cecum vs Env in Hakai
metadata <- subset(metadata, Sample_location == "Hakai")
metadata <- subset(metadata, body_site %in% c("cecum", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")

### subsamples from metadata ### Pisastar*cecum vs Env in Bamfield
metadata <- subset(metadata, Sample_location == "Bamfield")
metadata <- subset(metadata, body_site %in% c("cecum", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")

### subsamples from metadata ### Pisastar*cecum vs Env in PortMoody
metadata <- subset(metadata, Sample_location == "Port_Moody")
metadata <- subset(metadata, body_site %in% c("cecum", "Env"))
metadata <- subset(metadata, Type !="Dermasterias")
metadata <- subset(metadata, Type !="Henrecia")
metadata <- subset(metadata, Type !="Pycnopodia")


### subsamples from metadata ### Pisastar (Monterey) vs Env (Vancouver)
metadata <- subset(metadata, group1 %in% c("Monterey_Bay", "Env"))

### subsamples from metadata ### Pisastar (Monterey) vs Env (All)
metadata <- subset(metadata, body_site %in% c("biopsy_punch", "Env", "aboral"))
metadata <- subset(metadata, Type != "Dermasterias")
metadata <- subset(metadata, Type != "Henrecia")

### subsamples from metadata ### 	Symptomatic vs Asymptomatic (Hakai)
metadata <- subset(metadata, Sample_location =="Hakai")
metadata <- subset(metadata, disease_status %in% c("Symptomatic", "Asymptomatic"))
metadata <- subset(metadata, body_site == "aboral")

### subsamples from metadata ### 	Asymptomatic vs Env (Hakai)
metadata <- subset(metadata, Sample_location =="Hakai")
metadata <- subset(metadata, disease_status %in% c("Asymptomatic", "Env"))
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))

### subsamples from metadata ### 	Symptomatic vs Env (Hakai)
metadata <- subset(metadata, Sample_location =="Hakai")
metadata <- subset(metadata, disease_status %in% c("Symptomatic", "Env"))
metadata <- subset(metadata, body_site %in% c("aboral", "Env"))

### subsamples from metadata ### 	Symptomatic vs Asymptomatic (Montery)
metadata <- subset(metadata, Sample_location =="Monterey_Bay")

### subsamples from metadata ### 	Symptomatic vs Asymptomatic (All)
metadata <- subset(metadata, Type =="Pisaster")
metadata <- subset(metadata, body_site != "oral")

### subsamples from metadata ### 	Fucus vs Env (Bamfield)
metadata <- subset(metadata, Sample_location == "Bamfield")
metadata <- subset(metadata, Type %in% c("Fucus", "Rock", "seawater"))

levels(metadata$body_site)
metadata$body_site = factor(metadata$body_site)
levels(metadata$body_site)
table(metadata$body_site)

master_table <- left_join(metadata, otu_table, by = "Informative_sampleID")
View(as.data.frame(master_table))
## remove NA##
#na.omit(master_table)
#na.omit(master_table)
## replace NA with 0##
#master_table[is.na(master_table)] <- 0


### Loading your data (this should contain metadata + taxa abundances)
#fucus <- read.table("your_sample(rows)_taxa(columns)_table.txt", header = T)
seastar <- master_table
#head(fucus)

### Set factors you are interested in i.e. want to get core for fucus but using seawater as control, then use the column with info on sample_type

## use one of the interesting info
sample_type <-seastar$body_site

class(sample_type)
levels(sample_type)

### Creating an object to store abundances only so you can run the analysis (i.e. remove 7 first columns of metadata with dplyr)
seastar_abund <- seastar %>% dplyr::select(-(1:17)) 
#head(seastar_abund)

#### Multipatt analysis: indval with fucus and water ####
multipatt.seastar <- multipatt(seastar_abund, sample_type, control = how(nperm=999))
summary(multipatt.seastar)

## get output and save it
indaval_output <- capture.output(summary(multipatt.seastar, indvalcomp=TRUE))
write.table(as.data.frame(indaval_output), file = "~/desktop/seastar_project/R_analysis/IndVal_Result/Oral_Env_PortMoody.txt", quote=F, row.names=F, col.names=T, sep="n")

#########################
#########################
