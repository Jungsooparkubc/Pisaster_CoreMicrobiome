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
# filter and reformat data frames
library(tibble)      
# Needed for converting column to row names

#### environment settings ####
#set working directory
theme_set(theme_bw())
setwd("~/desktop/seastar_project/R_analysis/")

####Read Phyloseq data to work
project_data <- readRDS("Seastar_project.RDS") 

#### rarefy data ####
set.seed(24) 
project_data.rarefied <- rarefy_even_depth(project_data, sample.size = 1500)
### Save rarifying data###
saveRDS(project_data.rarefied, "Seastar_project_1500.RDS")

project_data <- readRDS("Seastar_project_1500.RDS") # after rarefaction to 1500


#remove a specific taxa (e.g., contaminant)
project_data <- project_data %>%
  subset_taxa(family != "Alcaligenaceae")

#remove a specific sample
project_data <- project_data %>%
  subset_samples(Informative_sampleID != "B32_HAK")

## subset data for data analyses

project_data <- project_data %>%
  subset_samples(Type != "Pycnopodia") %>%
  subset_samples(Type != "Dermasterias") %>%
  subset_samples(Type != "Henrecia") %>%
  subset_samples(Category != "Env") %>%
  subset_samples(Category != "Other_hosts") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  #subset_samples(sample_type != "Pisaster_oral") %>%
  subset_samples(body_site2 != "lesion") %>%
  #subset_samples(sample_type != c("Dermasterias_cecum", "Henrecia_cecum")) %>%
  subset_samples(disease_status != "Symptomatic") %>%
  #subset_samples(Sample_location != "Hakai") %>%
  subset_samples(Sample_location != "Monterey_Bay")

project_data <- project_data %>%
  subset_samples(Type != "Pycnopodia") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  #subset_samples(body_site2 != "lesion") %>%
  subset_samples(sample_type %in% c("Pisaster_aboral", "Pisaster_cecum", "Pisaster_lesion")) %>%
  subset_samples(Sample_location == "Hakai")

project_data <- project_data %>%
  subset_samples(Sample_location %in% "Monterey_Bay") %>%
  subset_samples(disease_status != "Dead")


project_data <- project_data %>%
  subset_samples(body_site %in% c("aboral", "oral", "cecum", "Rock_surface", "Seawater"))

##Pisaster Only and oral/aboral
project_data <- project_data %>%
  subset_samples(sample_type %in% c("Pisaster_aboral", "Pisaster_oral" )) %>%
  subset_samples(disease_status != "Symptomatic")
  #subset_samples(body_site == "aboral")

project_data <- project_data %>%
  subset_samples(sample_type %in% c("Pisaster_aboral", "Pisaster_oral"))

project_data <- project_data %>%
  subset_samples(sample_type %in% c("Pisaster_aboral", "Pisaster_oral", "Pisaster_biopsy_punch")) %>%
  subset_samples(disease_status != "Dead") %>%
subset_samples(disease_status != "Symptomatic")

project_data <- project_data %>%
  subset_samples(sample_type %in% c("Pisaster_cecum")) %>%
  subset_samples(disease_status != "Dead") %>%
  subset_samples(disease_status != "Symptomatic")
  

#for viewing
View(as.data.frame(unclass(sample_data(project_data))))
View(as.data.frame(unclass(otu_table(project_data))))
otu <- as.data.frame(unclass(otu_table(project_data)))
asv <- colSums(otu)
View(asv)
View(as.data.frame(unclass(tax_table(project_data))))

##count taxa function##
ntaxa(project_data)

which(is.na(tax_table(project_data)) == TRUE)
tax_table(project_data)[is.na(tax_table(project_data))] <- 0


#### Create Phyloseq Data### 
##Manually revise otu/taxonomic tables to match the format at https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

##Reverse otu table
trans<- read.csv(file="genus_level_final.otu.csv")
trans <- t(trans)
View (trans)
write.csv(trans, file="otu_rev.csv", row.names=T)
###################
otu_mat<- read.csv(file="otu_rev.csv")
tax_mat<- read.csv(file="genus_level_final.tax.csv")
samples_df <- read.csv(file="metadata_new.csv")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample_ID") 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

project_data <- phyloseq(OTU, TAX, samples)
project_data
saveRDS(project_data, "Seastar_project_asv_added.RDS")



#########################
#####Rarefaction Curve  PLOT#####

options(scipen=10000) # number format

pdf("Hakai_disease_rarefaction.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.

p<- rarecurve(t(otu_table(project_data)), step=100, col=sample_data(project_data)$disease_status, lwd=2, ylab="ASVs", label=F)

#disease_status
#, xlim=c(0,5000)
dev.off()

#########################
#####TAXONOMIC PLOT#####

plot(sort(sample_sums(project_data)), xaxt = "n", xlab='Sample Name')
axis(1, at=1:length(sample_names(project_data)), labels=sort(sample_names(project_data)))

taxonomy_plot_obj <- project_data %>%
  tax_glom(taxrank = "asv") # agglomerate at your rank of interest (Rank5 is roughly equal to "family" in our example)

# recommend roughly 20 taxa maximum, since it becomes more difficult to distinguish colours with more taxa than that
topOTUs <- names(sort(taxa_sums(taxonomy_plot_obj), TRUE)[1:200]) #where N is the number of taxa you want to retain for plotting purposes

# REQUIRED: transform to relative abundance, melt to long format for plotting
taxonomy_plot_obj <- taxonomy_plot_obj %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
  psmelt() %>%                                          # Melt to long format
  arrange(asv)     

# identify rows (taxa) in plotting object that belong to topOTUs, assign them the taxonomy string "other"
taxonomy_plot_obj$asv <- as.character(taxonomy_plot_obj$asv) #must assign this as character in order to introduce a new string "other"
taxonomy_plot_obj$asv[-which(taxonomy_plot_obj$OTU %in% topOTUs)] <- "other"
taxonomy_plot_obj$asv <- factor(taxonomy_plot_obj$asv) #convert back to a factor #necessary for plotting properly

# intense- 18 colours
cbPalette <- c("#441f00","#9243f3","#00ad00","#002bc9","#45e16b","#f073ff","#a3d649","#361258","#dab600","#0289d2","#fa150d","#008460","#bf003f","#d0c96f","#ff80c3","#606300","#ff9365","#a75400")


# # 3. [Optional] order levels you are interested in the way you would like them plotted. many ways to do this.
taxonomy_plot_obj$asv <- factor(taxonomy_plot_obj$asv, levels = c("others", "ASV118", "ASV39", "ASV10", "ASV7", "ASV4", "ASV1", "ASV17", "ASV3", "ASV2")) #method 1

taxonomy_plot_obj$asv <- factor(taxonomy_plot_obj$asv, levels = c("others", "ASV15", "ASV13", "ASV7", "ASV1", "ASV3", "ASV2")) #method 1

#### Make Plots ####

pdf("taxonomiplot_core_ceca.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 16, # define plot width and height (this is in Inches). completely up to user.
    height = 9
)

ggplot(taxonomy_plot_obj, aes(x = Informative_sampleID, y = Abundance, fill = asv)) + 
  facet_wrap(~Sample_location, ncol=4, strip.position = "top", drop=TRUE, scales="free") + #OPTIONAL LINE: facet_wrap is the function for determining how plots are grouped within a multi-plot space
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  theme(strip.text = element_text(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_bar(stat = "identity", width = 1.0) + #geom_bar controls the bar plots themselves
  scale_y_continuous(expand = c(0.005,0.005)) + 
  #scale_fill_manual(values=c("#4fc083","#33d4d1","#b94a73","#bd7f37","#6777cf","#acb140","#699643","#b94c3f","#9750a1")) +
  scale_fill_manual(values=c("#5a667e","#a4aac1","#bd7f37","#acb140","#b94c3f","#9750a1")) +#this controls the y axis scale, for bigger margins above and below, increase the numbers provided
  #scale_fill_manual(values = cbPalette) + #here's where to use the raw colour palate derived above
  # for colourblind: scale_fill_manual(values = cbPalette) + #example with colourblind palette
  #NOTE: to preserve colours between plots, you must use the line below, not the ones above
  # for custom palette: scale_fill_manual(values = myCustomPalette) + #example with custom palette
  ylab("Relative Abundance") + #x and y axis labels
  xlab("Samples") +
  ggtitle("Taxonomy Summary") + #plot title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")) #another "theme" command. a lot of extra control can be established here. this line ensures that there is no padding in the multi-plot grid

dev.off()

#########################
#### Basic alpha div plot ####
#chao1
pdf("Chao1_bodysite_location_final.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.

#### basic alpha div plot [Evenness: Pielou's measure] #####

alpha_diversity <- estimate_richness(project_data, measure = c("Chao1", "Shannon", "Observed"))
H <- alpha_diversity$Shannon
S1 <- alpha_diversity$Observed
S <- log(S1)
evenness <- H/S
evenness
alpha_diversity$Evenness = evenness
alpha_diversity

abundance <- sample_sums(project_data)

### merge with metadata####
sample_data(project_data)$evenness <- alpha_diversity$Evenness #add to metadata (the rows are in the same order already)
sample_data(project_data)$richness <- alpha_diversity$Chao1
sample_data(project_data)$evenness <- as.numeric(sample_data(project_data)$evenness)
sample_data(project_data)$evenness <- as.numeric(sample_data(project_data)$richness)

sample_data(project_data)$total_reads <- sample_sums(project_data)
### dataframe for ANOVA ###
df <- data.frame(sample_data(project_data))
anova_result <- aov(richness ~ body_site, df)
anova_result2 <- aov(evenness ~ body_site, df)
summary(anova_result)
## ANOVA (effects on total_reads)
anova_result <- aov(total_reads ~ body_site*disease_status, df)
summary(anova_result)

### dataframe for ANOVA ###
df <- data.frame(sample_data(project_data))
anova_result <- aov(richness ~ body_site, df)
anova_result2 <- aov(evenness ~ body_site, df)
summary(anova_result2)

## load csv data for core abundance for each core ANOVA 
setwd("~/desktop/seastar_project/R_analysis/ANOVA_core")
#df<- read.csv(file="monterey_core.csv", header = TRUE)
df<- read.csv(file="monterey_core_noGoo.csv", header = TRUE)
#df<- read.csv(file="hakai_core.csv", header = TRUE)
#df<- read.csv(file="Hakai_ceca_core.csv", header = TRUE)

df<- filter(df, taxa == "ASV39.Gracilibacteria")
anova <- aov(abundance ~ disease_stages2, data = df)
summary(anova)

anova <- aov(abundance ~ body_site2, data = df)
summary(anova)


library(agricolae)
tukey_result <- HSD.test(anova_result, "body_site", group = TRUE)
print(tukey_result)

tukey_result <- HSD.test(anova_result2, "body_site", group = TRUE)
print(tukey_result)

pdf("Richness_bodysite_final.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.
#### step 2: use calculated alpha diversity to make basic plot with ggplot ####
#this plot lets you customize things a bit more than the plot_richness function, if desired
group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]

p <- ggplot(sample_data(project_data), aes(x = factor(body_site, level = c('aboral', 'oral', 'cecum', 'Other_hosts', 'Rock_surface', 'Seawater')), y=richness))
p + geom_boxplot(outlier.size = -11) +
  #facet_grid(~ FACTOR_2, drop=TRUE, scales="free", space="free") + #drop, scales, and space are used to ensure that each boxplot bar is equally sized, and that empty levels created by the facet grid are dropped from the final plot
  labs(title="Alpha Diversity (Evenness)", x="body_site", y="Evenness") + 
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.text=element_text(size=15)) +
  #scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()


pdf("Simpson_bodysite_final.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.
#### basic alpha div plot ####
#### step1: calculate alpha diversity and add it as a column in the metadata ####
project_data.simpson <- estimate_richness(project_data, split = TRUE, measures = c("Simpson")) #estimate richness
sample_data(project_data)$simpson <- project_data.simpson$Chao1 #add to metadata (the rows are in the same order already)
sample_data(project_data)$simpson <- as.numeric(sample_data(project_data)$simpson)
#### step 2: use calculated alpha diversity to make basic plot with ggplot ####
#this plot lets you customize things a bit more than the plot_richness function, if desired
p <- ggplot(sample_data(project_data), aes(x = factor(body_site, level = c('aboral', 'oral', 'cecum', 'Other_hosts', 'Rock_surface', 'Seawater')), y=simpson, color=Sample_location))
p + geom_boxplot(outlier.size = -11) +
  #facet_grid(~ FACTOR_2, drop=TRUE, scales="free", space="free") + #drop, scales, and space are used to ensure that each boxplot bar is equally sized, and that empty levels created by the facet grid are dropped from the final plot
  labs(title="Alpha Diversity (Chao1)", x="body_site", y="Chao1") + 
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.text=element_text(size=15)) +
  #scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()




pdf("Simpson_bodysite_final.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 12, height = 9)# define plot width and height. completely up to user.

#### basic alpha div plot ####
#### step1: calculate alpha diversity and add it as a column in the metadata ####
project_data.simpson <- estimate_richness(project_data, split = TRUE, measures = c("Simpson")) #estimate richness
sample_data(project_data)$simpson <- project_data.simpson$Simpson #add to metadata (the rows are in the same order already)
sample_data(project_data)$simpson <- as.numeric(sample_data(project_data)$simpson)
#### step 2: use calculated alpha diversity to make basic plot with ggplot ####
#this plot lets you customize things a bit more than the plot_richness function, if desired
p <- ggplot(sample_data(project_data), aes(x = factor(body_site, level = c('aboral', 'oral', 'cecum', 'Other_hosts', 'Rock_surface', 'Seawater')), y=simpson))
p + geom_boxplot(outlier.size = -11) +
  #facet_grid(~ FACTOR_2, drop=TRUE, scales="free", space="free") + #drop, scales, and space are used to ensure that each boxplot bar is equally sized, and that empty levels created by the facet grid are dropped from the final plot
  labs(title="Alpha Diversity (Chao1)", x="body_site", y="Chao1") + 
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.text=element_text(size=15)) +
  #scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()


which(is.na(tax_table(project_data)) == TRUE)
project_data[is.na(project_data)] <- 0
###########Beta Diversity##############
set.seed(24)
NMDS.bray <- ordinate(
  physeq = project_data,
  method = "NMDS",
  distance = "bray"
) # you can choose different methods and distance metrics, see the ordinate help page for details. this function works with "phyloseq" class objects.

#### making beta div plots ####
#we get more plotting control if we don't use the phyloseq plotting functions for ordination plots, and instead add the results of the ordination to our existing metadata
NMDS <- as.data.frame(unclass(sample_data(project_data)))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$type2, NMDS$body_site),]

## manually designate colours 15
cbPalette <- c("#ca2700", "#972a3a", "#ff5e37", "#d1ff9b", "#95f926", "#58e72a", "#b4943e", "#4500aa", "#7835d3", "#ffb3cf", "#d9a8d1", "#8b30bb", "#933b8c", "#131212", "#7fdcf6")
## manually designate colours 55
cbPalette <- c("#c1873d", "#a069be", "#6da452", "#cc5262", "#2baad3")
cbPalette <- c("#dac4ff","#dee799","#ed98c7","#7bdaaf","#f3938e","#7dfff4", "#ffcb9b","#3db8e2","#aadb98","#b5e1ff")

###############NMDS PLOT2: Sample Type###############
pdf("NMDS_seastar_Hakai_symptom_asymptom_lesion.pdf"
    , width = 11 # Default is 7
    , height = 10 # Change to 10; make it taller
)
p <- ggplot(NMDS.sort, aes(x=NMDS.bray1, y=NMDS.bray2, color = Sample_location, shape=body_site)) # change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=5, alpha=0.85) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values = c(16, 17, 8, 15)) +
  geom_text() +
  labs(title="NMDS_seastar_all_sampletypes_coloured") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#xlim(0.8, -0.8) + ylim(-0.4, 0.5)
dev.off()


###############NMDS PLOT2: Location###############
## manually designate colours 4
cbPalette <- c("#ca2700", "#00c250", "#006f99", "#8b30bb")
NMDS.sort <- NMDS[order(NMDS$Sample_location, NMDS$body_site),]

pdf("NMDS_Figure1B.pdf"
    , width = 11 # Default is 7
    , height = 10 # Change to 10; make it taller
)
p <- ggplot(NMDS.sort, aes(x=NMDS.bray1, y=NMDS.bray2, color = Sample_location, shape = body_site)) # change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=5, alpha=0.85) + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  labs(title="NMDS_seastar_all_sampletypes_coloured") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#xlim(-3, 3) + ylim(-3, 3)
dev.off()

############Make dot plot [BUBBLE PLOT]#############
####################################################
theme_set(theme_bw())
setwd("~/desktop/seastar_project/R_analysis/Bubble plot/")
###upload data
New_plots<- read.csv("Averages_for_Dot_plot_seastarcore3_final.csv")   #Pisaster overall Fig.2A
New_plots<- read.csv("Averages_for_Dot_plot_seastarcore3_final_noMont.csv") #Pisaster overall Fig.2A w/o Mont
New_plots<- read.csv("Averages_for_Dot_plot_cecacore3_final.csv") #Pisaster ceca Fig.2B
New_plots<- read.csv("Averages_for_Dot_plot_seastarcore_monterey.csv")
head(New_plots)

#### Go to wide/fat format
Long_np<-gather(New_plots,Taxa,Relative_Abundance, ASV2.Spirochete, ASV3.Spirochete, ASV17.Spirochete, ASV1.Salinispira,	ASV4.Salinispira,	ASV7.Reichenbachiella, ASV10.Reichenbachiella, ASV39.Gracilibacteria, ASV118.Amoebophilaceae, All.core.ASVs, factor_key=TRUE) #surface

Long_np<-gather(New_plots,Taxa,Relative_Abundance, ASV2.Spirochete, ASV3.Spirochete, ASV1.Salinispira, ASV7.Reichenbachiella, ASV13.Hepatoplasma, ASV15.Hepatoplasma, All.core.ASVs, factor_key=TRUE) #ceca

#Long_np<-gather(New_plots,Taxa,Relative_Abundance, ASV2.Spirochete, ASV3.Spirochete, ASV17.Spirochete, ASV1.Salinispira,	ASV4.Salinispira,	ASV7.Reichenbachiella, ASV10.Reichenbachiella, ASV39.Unidentified, ASV118.Unidentified, factor_key=TRUE)
head(Long_np)
#remove Rows having NAs
Long_np <- Long_np[rowSums(is.na(Long_np)) == 0,]

#Order your factors for plot 
#Long_np$Samples <- factor(Long_np$Samples, level = c("All", "Bamfield", "PortMoody", "Hakai", "Hakai_s", "Monterey", "Monterey_s"))
Long_np$Samples <- factor(Long_np$Samples, level = c("All Pisaster", "Pisaster Bamfield", "Pisaster PortMoody", "Pisaster Hakai", "Pisaster Hakai(S)","Pisaster Hakai(S-lesion)", "Pisaster Monterey", "Pisaster Monterey(S)", "Dermasterias Hakai", "Henricia Hakai", "Other hosts", "Rock", "Seawater"))

Long_np$Status <- factor(Long_np$Status, level = c("zero","one","two","three","four","five"))
#Long_np$Samples <- factor(Long_np$Samples, level = c("Pisaster", "Dermasterias", "Henrecia"))

###Make dot plot 
pdf("buble_plot_seastar_disease_status_core_final_OCT.22_yellow.pdf"
    , width = 7 # Default is 7
    , height = 5 # Change to 10; make it taller
)
Long_np$Relative_Abundance <- as.numeric(as.character(Long_np$Relative_Abundance))
Dot_plot_Survey_core<- ggplot(Long_np, aes(y= Taxa, x = Status)) + geom_point(aes(size = Relative_Abundance, color = "#7f64b9")) + #, color = Status))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #scale_color_manual(values = c("#bb7438","#7f64b9", "#72ac5c", "#b94b75")) +
  scale_color_manual(values = c("#bb7438")) +
  scale_size(range=c(-0.5,7)) +  theme_bw() + 
  #breaks = c("1","3","6","9") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))  +
  theme(axis.text.x=element_text(size=8, angle= 45, hjust= 1), axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold")) + ylab("Core Taxa") + xlab("Location") + ggtitle("")

Dot_plot_Survey_core

dev.off()

######################CECA##########
New_plots<- read.csv("Averages_for_Dot_plot_cecacore3_final.csv")
head(New_plots)

#### Go to wide/fat format
Long_np<-gather(New_plots,Taxa,Relative_Abundance, ASV2.Spirochete,	ASV3.Spirochete,	ASV1.Salinispira,	ASV7.Reichenbachiella, ASV13.Hepatoplasma,	ASV15.Hepatoplasma, All.core.ASVs, factor_key=TRUE)
head(Long_np)
#remove Rows having NAs
Long_np <- Long_np[rowSums(is.na(Long_np)) == 0,]


#Order your factors for plot 
#Long_np$Samples <- factor(Long_np$Samples, level = c("All", "Bamfield", "PortMoody", "Hakai", "Hakai_s"))
Long_np$Samples <- factor(Long_np$Samples, level = c("All Pisaster", "Pisaster Bamfield", "Pisaster PortMoody", "Pisaster Hakai", "Pisaster Hakai(S)", "Dermasterias Hakai", "Henricia Hakai", "Other hosts", "Rock", "Seawater"))
#Long_np$Samples <- factor(Long_np$Samples, level = c("Pisaster", "Dermasterias", "Henricia"))

###Make dot plot 
pdf("buble_plot_seastar_ceca_core3_final.pdf"
    , width = 7 # Default is 7
    , height = 5 # Change to 10; make it taller
)
Long_np$Relative_Abundance <- as.numeric(as.character(Long_np$Relative_Abundance))
Dot_plot_Survey_core<- ggplot(Long_np, aes(y= Taxa, x = Samples)) +geom_point(aes(size = Relative_Abundance, color = Status))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_size(range=c(-0.5,9)) +  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))  +
  theme(axis.text.x=element_text(size=8, angle= 45, hjust= 1), axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold")) + ylab("Core Taxa") + xlab("Location") + ggtitle("")

Dot_plot_Survey_core

dev.off()


####data subsetting#####
project_data <- readRDS("Seastar_project_1500_2.RDS")
project_data <- readRDS("Seastar_project_1500_ceca.RDS")

#Hakai healthy aboral vs lesion 
project_data <- project_data %>%
  subset_samples(Type == "Pisaster") %>%
  subset_samples(Sample_location == "Hakai") %>%
  subset_samples(body_site %in% c("aboral", "lesion"))
  #subset_samples(body_site %in% c("cecum"))


project_data <- project_data %>%
  subset_samples(Type != "Dermasterias") %>%
  subset_samples(Type != "Henrecia") %>%
  subset_samples(body_site2 != "lesion") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  subset_samples(body_site2 != "cecum") %>%
  subset_samples(disease_status != "Symptomatic") %>%
  subset_samples(Sample_location != "Monterey_Bay")


#only asymptomatic pisaster at hakai, bamfield and Port Moody
project_data <- project_data %>%
  subset_samples(Type == "Pisaster") %>%
  subset_samples(body_site2 != "lesion") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  #subset_samples(body_site2 != "cecum") %>%
  subset_samples(disease_status != "Symptomatic") %>%
  subset_samples(Sample_location != "Monterey_Bay")

project_data <- project_data %>%
  subset_samples(Category == "Star") %>%
  subset_samples(Type != "Pycnopodia") %>%
  subset_samples(body_site2 != "lesion") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  subset_samples(body_site2 != "cecum") %>%
  #subset_samples(Sample_location %in% "Monterey_Bay") %>%
  subset_samples(disease_status != "Dead")

project_data <- project_data %>%
  #subset_samples(Type != "Dermasterias") %>%
  #subset_samples(Type != "Henrecia") %>%
  #subset_samples(Type != "Pycnopodia") %>%
  #subset_samples(body_site2 != "lesion") %>%
  subset_samples(Type == "Pisaster") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  subset_samples(body_site2 != "cecum") %>%
  subset_samples(disease_status != "Symptomatic") %>%
  subset_samples(Sample_location == "Hakai") %>%
  subset_samples(body_site != "oral")
  #subset_samples(body_site %in% c("aboral", "oral"))

project_data <- project_data %>%
  subset_samples(Type %in% c("Pisaster", "abiotic")) %>%
  subset_samples(Sample_location != "Monterey_Bay")


View(as.data.frame(unclass(sample_data(project_data))))

project_bray <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))

adonis2(project_bray ~ body_site * Sample_location,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_bray ~ body_site,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_bray ~ disease_status,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_bray ~ Sample_location,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_bray ~ body_site + Sample_location + body_site:Sample_location,
        data= sample_df , permutations=10000, by_margin = "T")
adonis2(project_bray ~ disease_status + body_site + disease_status:body_site,
        data= sample_df , permutations=10000, by_margin = "T")

########################################################
########################################################
########################################################
########################################################
####### add phylogenetic tree to phyloseq class######
random_tree = rtree(ntaxa(project_data), rooted=TRUE, tip.label=taxa_names(project_data))
#plot(random_tree)
project_data <- merge_phyloseq(project_data, random_tree)


####### START PERMANOVA ANALAYSIS######

project_bray <- phyloseq::distance(project_data, method = "bray")
project_jaccard <- phyloseq::distance(project_data, method = "jaccard")
project_wunifrac <- phyloseq::distance(project_data, method = "wunifrac")

sample_df <- data.frame(sample_data(project_data))


adonis(project_bray ~ sample_type, data=sample_df, method="bray")
adonis(project_jaccard ~ Type, data=sample_df, method="jaccard")
adonis(project_wunifrac ~ Type, data=sample_df, method="wunifrac")

adonis2(project_bray ~ body_site,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_jaccard ~ Type,
        data= sample_df , permutations=10000, by = "margin")
adonis2(project_wunifrac ~ Type,
        data= sample_df , permutations=10000, by = "margin")

res.adonis.rarefied <- adonis(project_bray ~ Type, data=sample_df, method="bray")
#res.adonis.rarefied <- adonis(project_bray.rarefied ~ Month, data=sample_df, method="bray")
res.adonis.rarefied #run for full description of results
summary(res.adonis.rarefied) #run for results summary (this is less informative if I remember correctly)

#install adnois.pair
#devtools::install_github("hadley/devtools")
#install.packages("devtools")
library(devtools)
library(EcolUtils)


## running 10,000 permutations instead of the default 999###
adonis2(project_bray ~ sample_type2,
        data= sample_df , permutations=999, by = "margin")
options(scipen=999)

pair_type <- adonis.pair(project_bray, sample_df$body_site,
            nper = 10000, corr.method = "BH")
#adonis.pair(project_jaccard, sample_df$Phylum,
          #  nper = 1000, corr.method = "BH")

capture.output(file="adnois_pair_type_2.txt", pair_type)

#### analysis: beta dispersion test ####
beta.Type <- betadisper(project_bray, sample_df$Type) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
permutest(beta.Type) #do permutation test of above #this is more customizable, see documentation


beta.Month <- betadisper(project_bray.rarefied, sample_df$Month)
permutest(beta.Month)

#### saving output of these tests ####
# you can use the "capture.output()" function to save the output of any of the testing and summary commands above, like so
capture.output(file="beta_disp.output.txt", permutest(beta.Phylum))
capture.output(file="beta_disp.output.txt", permutest(beta.Order))
# where the form of the command is always:
capture.output(file="name_of_output_file", command)



#######PERMANOVA spearate pair-wise for each site
project_data <- readRDS("Seastar_project_1500.RDS")

#only asymptomatic pisaster at hakai, bamfield and Port Moody
project_data <- project_data %>%
  subset_samples(Type == "Pisaster") %>%
  subset_samples(body_site2 != "lesion") %>%
  subset_samples(body_site2 != "coelomic fluid") %>%
  #subset_samples(Sample_location == "Hakai") %>%
  subset_samples(Sample_location == "Port_Moody") %>%
  #subset_samples(Sample_location == "Bamfield") %>%
  subset_samples(disease_status != "Symptomatic") %>%
  subset_samples(Sample_location != "Monterey_Bay")

View(as.data.frame(unclass(sample_data(project_data))))

project_bray <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))
pair_type <- adonis.pair(project_bray, sample_df$body_site,
                         nper = 10000, corr.method = "BH")
pair_type

capture.output(file="adnois_pair_bodysite_PortMoody.txt", pair_type)
