#### Import the workspace image ####

save.image(file = "euk_R_26-08.RData")

#### Import packages ####

#### base ####
library(BiocManager)
library(tidyverse)
library(openxlsx)
library(reshape2)

#### analysis ####
library(phyloseq)
library(vegan)
library(rstatix)
library(ggpmisc)
library(NetCoMi)
library(limma)
library(factoextra)
#library(lme4)
#library(lmerTest)
library(dtw)
library(dtwclust)

#### visual ####
library(ggConvexHull)
library(magrittr)
library(ggtext)
library(tidytext)
library(ggh4x)
library(patchwork)

# Since two of NetCoMi's dependencies are only available on GitHub, 
# it is recommended to install them first:
#FOR MAC: needed to install gfortran  mac first
#https://mac.r-project.org/tools/gfortran-12.2-universal.pkg
#devtools::install_github("zdk123/SpiecEasi")
#devtools::install_github("GraceYoon/SPRING")
#FOR LINUX:needed to instal gdl with
#sudo  apt-get install gsl-bin libgsl0-dev
#devtools::install_github("stefpeschel/NetCoMi",repos = c("https://cloud.r-project.org/",BiocManager::repositories()))

####| ####

#### Functions ####

#create phyloseq
create_phyloseq_SILVA <- function(x,y,z) {
  #Import ASV_table
  ASV_table <- read_tsv(x,show_col_types = F,num_threads = 5,skip = 1) %>%
    as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "#OTU ID")
  #Import Tax_table
  tax_table <- read_tsv(y,show_col_types = F,num_threads = 5) %>%
    as.data.frame() %>%
    select(-"Confidence") %>%
    dplyr::mutate(.,
                  Taxon = gsub("; ","%",Taxon),
                  Taxon = gsub("d__","",Taxon),
                  Taxon = gsub("%.__","%",Taxon)) %>%
    separate_wider_delim(names = c("Domain","Phylum","Class","Order",
                                   "Family","Genus","Species"),
                         Taxon,delim = "%", too_few = "align_start") %>%
    remove_rownames() %>% column_to_rownames(var = "Feature ID") %>%
    subset(rownames(.) %in% rownames(ASV_table))
  #Import metadata
  metadata <- read_tsv(z,show_col_types = F,num_threads = 5) %>%
    as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "sample-id")
  #physeq object
  return(phyloseq(OTU = otu_table(as.matrix(ASV_table), taxa_are_rows = T),
                  TAX = tax_table(as.matrix(tax_table)),
                  samples = sample_data(metadata)))
}

create_phyloseq_PR2 <- function(x,y,z) {
  #Import ASV_table
  ASV_table <- read_tsv(x,show_col_types = F,num_threads = 5,skip = 1) %>%
    as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "#OTU ID")
  #Import Tax_table
  tax_table <- read_tsv(y,show_col_types = F,num_threads = 5) %>%
    as.data.frame() %>%
    #select(-"Confidence") %>%
    # dplyr::mutate(.,Taxon = gsub(";","%",Taxon), Taxon = gsub("d__","",Taxon)) %>%
    separate_wider_delim(names = c("Domain","Supergroup","Division","Subdivision",
                                   "Class","Order","Family","Genus","Species","Sub_Species"),
                         Taxon,delim = ";", too_few = "align_start") %>%
    remove_rownames() %>% column_to_rownames(var = "Feature ID") %>%
    subset(rownames(.) %in% rownames(ASV_table))
  #Import metadata
  metadata <- read_tsv(z,show_col_types = F,num_threads = 5) %>%
    as.data.frame() %>% remove_rownames() %>% column_to_rownames(var = "sample-id")
  #physeq object
  return(phyloseq(OTU = otu_table(as.matrix(ASV_table), taxa_are_rows = T),
                  TAX = tax_table(as.matrix(tax_table)),
                  samples = sample_data(metadata)))
}

#inverse %in%
`%out%` <- Negate(`%in%`)

#get MDS output
get_MDS_output <- function(x,y,z) {
  #perform the MDS
  MDS <- cmdscale(y,k = z,eig = T)
  #extract the valuable outputs
  MDS_points <- as.data.frame(MDS$points) %>% rename("Axis1"="V1", "Axis2"="V2")
  MDS_species <- as.data.frame(wascores(MDS_points, t(x),expand = FALSE))
  Axis1 <- round(MDS$eig[[1]]/sum(MDS$eig),3)*100
  Axis2 <- round(MDS$eig[[2]]/sum(MDS$eig),3)*100
  #return a list
  return(list(MDS_points, MDS_species, Axis1, Axis2))
}

#Get ASVs acounting for >= 1% rel. abund. at least once
get_asv_above_threshold <- function(physeq_object, threshold = 1) {
  #get the ASV table
  otu_table <- otu_table(physeq_object)
  #transform into rel. abund.
  otu_table_percent <- apply(otu_table, 2, function(x) {
    x / sum(x) * 100
  })
  #get ASVs at least once >= 1%
  asv_above_threshold <- which(apply(otu_table_percent, 1, function(x) {
    any(x >= threshold)
  }), arr.ind = TRUE)
  #output
  asv_names <- rownames(otu_table)[asv_above_threshold]
  return(asv_names)
}

#Get ASVs in X% of samples
get_asv_min_samples <- function(physeq_object, percentage) {
  #get the ASV table
  otu_table <- otu_table(physeq_object)
  #nb. of samples for each ASVs
  presence_matrix <- otu_table > 0
  asv_presence <- rowSums(presence_matrix)
  #min. nb. of samples
  min_samples <- round(ncol(otu_table)*percentage)
  #Get ASVs in X% of samples
  asv_above_threshold <- which(asv_presence >= min_samples)
  #output
  asv_names <- rownames(otu_table)[asv_above_threshold]
  return(asv_names)
}

####| ####

#### Formulas ####
lm_formula_linear<-y~x
lm_formula_poly2<-y~x+I(x^2)
lm_formula_poly3<-y~x+I(x^2)+I(x^3)
lm_formula_poly4<-y~x+I(x^2)+I(x^3)+I(x^4)
lm_formula_poly5<-y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)

####| ####

#### Color palettes ####

palette_season<-c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B")

palette_Class <- c("#2C5234", #gyrista non-diatoms
                   '#5b6f0d','#809c13', #gyrista diatoms
                   '#71ab68','#9bc495','#c2ff94',#Chlorophyta
                   '#f3efa0','#dfd96d',#Cryptophyta and Dinoflagellata
                   '#F6BC88','#F3A057','#D8BB8C', #Ciliophora
                   '#d2af79','#bd9568', #Ciliophora
                   '#a67b58',#Kathablepharida
                   "#CF8A4F","#C27933",#fungi
                   '#d15c14', #Cercozoa
                   '#8B3D0D',"#74330B",#Bigyra
                   '#60372c', #Choanoflagellata
                   '#49231e', #Chrompodellids
                   "#321010", #Telonemia
                   '#E40000',#Fungi parasites
                   '#ff4343','#ff807d',#Perkinsea
                   'lightgrey') #not top 20 class
                         
pie(rep(1,length(palette_Class)), col=(palette_Class))

palette_Phyto_Class <- c('#2C5234','#35633F','#3E744A', #gyrista non-diatoms
                         '#5b6f0d','#809c13','#a5c919', #gyrista diatoms
                         '#71ab68','#c2ff94', #Chlorophyta
                         '#f3efa0','#dfd96d', #Cryptophyta and Dinoflagellata
                         '#439798', #Euglenozoa
                         '#43dd98','#A2EECC', #Streptophyta
                         '#CFF2FB') #Cyanobacteria

pie(rep(1,length(palette_Phyto_Class)), col=(palette_Phyto_Class))

trophic_state<- c("#28A448","#94D2A4","#2F6B9D","#97B5CD")

palette_lake_chla<-c("#094f29","#247549","#537B2F","#8DA750","#BBBE69","#E4EB9C",
                     "#8CB8E2","#5C9BD6","#005696")

pie(rep(1,length(palette_lake_chla)), col=(palette_lake_chla))

####________________________####