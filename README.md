# R code for ecological data analysis

by Ling Wang

- Developed to process QIIME 2 objects using Phyloseq and qiime2R 
- Phyloseq: https://joey711.github.io/phyloseq/index.html
- Subjects: data handling, analysis and visulisation of amplicon sequencing data

## Data 

- Amplicon sequencing data (bacteria in mouse gut)
- Compare results with previous qiime graphs

## Study design

- WT_mucus
- WT_lumen 
- WT_colon
- KO_mucus
- KO_lumen
- KO_colon 

## 1. Check the environment

```
rm(list=ls()) # Remove all the variables to clear the workspace at start of job 
getwd() # get current working directory

# Change the working directory
filepath <- "~/Google Drive/R Programming" # If on Mac
#filepath <- "C:/Users/lwang18/Google Drive/R Programming" # If on PC

# Set working directory into a new path
setwd(filepath)
dir(filepath)
```
See a short list of the most useful R commands
https://www.personality-project.org/r/r.commands.html
http://stackoverflow.com/questions/2615128/where-does-r-store-packages

## 2. Load phyloseq, qiime2R and other packages

```
#BiocManager::install("phyloseq")
library(phyloseq) 

library(qiime2R)
library(readr) 
library(ape)
library(ggplot2)
library(Biostrings) 
library(biomformat) 
library(Hmisc) 
library(yaml) 
library(tidyr)
library(dplyr) 
library(stats) 
library(utils) 

#see where your libraries are stored
.libPaths()
```

## 3. Import data

In my filepath "~/Google Drive/R Programming", I have the following files:
- pyy_wtko_16s_table-with-phyla-no-mitochondria-no-chloroplast.qza
- pyy_wtko_sample-metadata-fixeds.tsv
- pyy_wtko_16s_taxonomy.qza
- pyy_wtko_16s_rooted-tree.qza

```
# Read Metadata
metadata<-read_q2metadata("pyy_wtko_sample-metadata-fixeds.tsv")
  head(metadata)

# Read OTU table
OTUs<-read_qza("pyy_wtko_16s_table-with-phyla-no-mitochondria-no-chloroplast.qza")
  OTUs$data[1:10,1:15] #show first 10 taxa within first 15 samples

# Read Taxonomy
taxonomy<-read_qza("pyy_wtko_16s_taxonomy.qza")
head(taxonomy$data)

#Creating a Phyloseq Object
phylo<-qza_to_phyloseq(
  features="pyy_wtko_16s_table-with-phyla-no-mitochondria-no-chloroplast.qza",
  tree="pyy_wtko_16s_rooted-tree.qza",
  "pyy_wtko_16s_taxonomy.qza",
  metadata = "pyy_wtko_sample-metadata-fixeds.tsv")
phylo

rank_names(phylo) 
sample_variables(phylo)
sample_names(phylo)
get_taxa_unique(phylo, "Family")
tax_table(phylo)[1:40] 
levels(sample_data(phylo)$geno)
levels(sample_data(phylo)$location)
```

## 4. Examine the data
```
# Examine rarefaction cureve
rarecurve(t(otu_table(phylo)), step=50, cex=0.2,label = F)

# Total reads per sample

readsumsdf = data.frame(nreads = sort(taxa_sums(phylo), TRUE), 
                        sorted = 1:ntaxa(phylo), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(phylo), 
                                                        TRUE), 
                                          sorted = 1:nsamples(phylo), 
                                          type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + 
  facet_wrap(~type, 1, scales = "free")+
  ylab("Numbers of reads") + 
  xlab(" ") 
```

## 5. Data normalization

```
# random subsampling
set.seed(28132)
phylo_rarefy = rarefy_even_depth(phylo, sample.size = 10000) 


# sanity-check: for phylo_rarefy, all sum up to 10000 reads
par(mfrow = c(1, 2))

plot(sort(sample_sums(phylo), TRUE), type = "h", main = "Before rarefaction", 
     ylab = "Numbers of reads", xlab = "samples", ylim = c(0, 30000))

plot(sort(sample_sums(phylo_rarefy), TRUE), type = "h", main = "After rarefaction", 
     ylab = "Numbers of reads", xlab = "samples", ylim = c(0, 12000))
```

```
# See what's removed after rarefaction
require(sqldf)
require(dplyr) 

a1 <- data.frame(a1=sample_names(phylo))
a2 <- data.frame(a2=sample_names(phylo_rarefy))

a1NotIna2 <- sqldf('SELECT * FROM a1 EXCEPT SELECT * FROM a2')
a1Ina2 <- sqldf('SELECT * FROM a1 INTERSECT SELECT * FROM a2')

a1NotIna2 # what's removed after rarefaction
a1Ina2 # what's left after rarefaction
```

## 6. Alpha diversity

- Observed OTU counts: simple richness of the community
- Shannon's index reflects ecological richness and eveness of the community

```
Observed <-
  plot_richness(phylo_rarefy, x = "location",color = "geno",
                       measures=c("Observed"), nrow=1) +geom_boxplot()+
  #ggtitle("Observed OTUs") +
  scale_color_manual(values = c("cornflowerblue", viridis::viridis(6)))+
  #theme_q2r() +
  theme(axis.title.y = element_text(size=8,color="Black"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=8,color="Black"),
      axis.text.x = element_text(size=8,color="Black"),
      axis.ticks = element_blank(),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"),
      #legend.position="bottom",
      
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "solid",size=0.2,colour = "grey"),
      panel.background = element_rect(fill = "white", colour = "black", size=1.2),
      
      strip.text = element_blank())
# Remove first laywer 
# see https://github.com/joey711/phyloseq/issues/836
# https://stackoverflow.com/questions/10493084/ggplot2-jitter-and-position-dodge-together
Observed$layers <- Observed$layers[-1]
Observed
```

```
shannon=
  plot_richness(phylo_rarefy, x = "location",color = "geno",
                      measures=c("shannon"), nrow=1) +geom_boxplot()  + 
  #ggtitle("Shannon's index") +
  scale_color_manual(values = c("cornflowerblue", viridis::viridis(6)))+
  #theme_q2r() +
  theme(axis.title.y = element_text(size=8,color="Black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8,color="Black"),
        axis.text.x = element_text(size=8,color="Black"),
        axis.ticks = element_blank(),
        
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"),
        #legend.position="bottom",
        
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "solid",size=0.2,colour = "grey"),
        panel.background = element_rect(fill = "white", colour = "black", size=1.2),
        
        strip.text = element_blank()) 
  
#+scale_y_continuous(limits = c(0, 7))

shannon$layers <- shannon$layers[-1]
shannon
```

Plot both graphs side by side 
```
require(gridExtra)
grid.arrange(Observed,  shannon, nrow=2)
```

## 7. Beta diversity

Ordination methods
- Non-metric multidimensional scaling (NMDS): species presence/absence
- Principal coordinates analysis (PCoA): 

Distance metrics
- Bray-Curtis similarity index
- Weighted Unifrac
- Unweighted Unifrac

```
# NMDS bray 
phylo_rarefy_nmds <- ordinate(phylo_rarefy, "NMDS", "bray")

NMDS=
  plot_ordination(phylo_rarefy, phylo_rarefy_nmds, 
                  color = "geno", shape="location") +
  geom_point(size = 4, alpha = 1) + 
  #ggtitle("PYY") +
  #scale_color_manual(values=taxa_palette_mg)  +
  scale_colour_viridis_d()+
  theme_q2r() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank()) +
  facet_wrap(~location)+stat_ellipse()
```

```
# PCoA_unifrac_W
phylo_rarefy_unifrac_W = ordinate(phylo_rarefy, method="PCoA", 
                                  distance="unifrac", weighted=T)


PCoA_unifrac_W=
  plot_ordination(phylo_rarefy, phylo_rarefy_unifrac_W, 
                  color =  "geno", shape="location") +
  geom_point(size = 4, alpha = 1) + 
  #ggtitle("PYY") +
  #scale_color_manual(values=taxa_palette_mg)  +
  scale_colour_viridis_d()+
  theme_q2r() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank()) 

PCoA_unifrac_W_wrap=
  plot_ordination(phylo_rarefy, phylo_rarefy_unifrac_W, 
                  color =  "geno", shape="location") +
  geom_point(size = 4, alpha = 1) + 
  #ggtitle("PYY") +
  #scale_color_manual(values=taxa_palette_mg)  +
  scale_colour_viridis_d()+
  theme_q2r() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank()) +
  facet_wrap(~location)+stat_ellipse()
```

```
grid.arrange(PCoA_unifrac_W_wrap, PCoA_unifrac_W, ncol=1)
```








