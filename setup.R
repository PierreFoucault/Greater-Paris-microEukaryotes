#### Import packages ####

library(BiocManager)
library(phyloseq)
library(vegan) 
library(tidyverse)
library(tidytext)
library(ggh4x)
library(patchwork)
library(rstatix)
library(reshape2)
library(lmerTest)
library(lme4)
library(dtw)
library(dtwclust)
library(ggConvexHull)
library(magrittr)
library(ggtext)

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

####| ####

#### Formulas ####
lm_formula_linear<-y~x
lm_formula_poly2<-y~x+I(x^2)
lm_formula_poly3<-y~x+I(x^2)+I(x^3)
lm_formula_poly4<-y~x+I(x^2)+I(x^3)+I(x^4)
lm_formula_poly5<-y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)

####| ####

#### Color palettes ####

# palette_easter<-c("#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC",
#                   "#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#88F393","#2E6D74","#DCBCF3","#EFE0DC","#F18B6C",
#                   "#B52239","#F5D2DB","#B7F7E6","#517BE6","#E0DBC1","#E19B8E","#BDCEC4","#99F8F3","#E45632","#D22A82",
#                   "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
#                   "#DA56EF","#EC97B3","#E0789F","#8CB8E2","#E89E34","#8D37A3","#71992A","#D8C5B2","#E75FA6","#EBE49A",
#                   "#ADAC5A","#BDF6D3","#3BDCA8","#47D1C8","#F8EBD4","#F4F6D8","#B6A7BD","#D6F0F2","#3EA34E","#B8E2EC",
#                   "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73","#A899F0")

palette_season<-c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B")

palette_Subdivision <- c('#325d3b','#71ab68','#c2ff94',
                         '#f3efa0','#dfd96d','#d2af79',
                         '#bd9568','#a67b58','#8f6349',
                         '#784c3a','#60372c','#49231e','#321010','#ff4343','#ff807d', 
                         '#ffb3b3','#4A4A4A')

trophic_state<- c("#28A448","#94D2A4","#2F6B9D","#97B5CD")

palette_lake_chla<-c("#28A448",
                     "#537B2F","#8DA750","#BBBE69","#E4EB9C","#94D2A4",
                     "#8CB8E2","#5C9BD6","#2F6B9D")

#Foucault 2025
#palette_lake_3T<-c("#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D")

pie(rep(1,length(palette_lake_chla)), col=(palette_lake_chla))


####|####

####________________________####
#### Old Amaury ####
####________________________####

#### Amaury colors ####
Groups <-rev(c("#28A448","#6F6F6F","#28A448",
               "#28A448","#28A448","#6F6F6F",
               "#5C9BD6","#5C9BD6","#5C9BD6"))

palette_lake_3T<-rev(c("#28A448","#1E7A36","#980000",
                       "#E40000","#E45000","#EB7947",
                       "#F2AB8C","#5C9BD6","#2F6B9D"))

palette_lake_3T <- c('#CCCCCC', '#B3B3B3', '#999999', '#808080', '#666666', '#4D4D4D', '#333333', '#1A1A1A', '#000000')

palette_lake_3T<-rev(c("#262626","#404040","#595959",
                       "#737373","#8C8C8C","#A6A6A6",
                       "#BFBFBF","#D9D9D9","#F2F2F2"))




palake <- rev(c('#BDEFA8','#81C384','#EFA8A8',
                '#C38B81','#EFD4A8','#C3BA81',
                '#EFEFA8','#A8EFEE','#8194C3'))

palette_season <- c('honeydew4', 'honeydew3', 'honeydew2', 'honeydew1' )
palette_season <- c("#BBB0BB", "#cfc3ab", "#E0EEED", "#ECFFED")
palette_season <- c("#808080", "#BEB050", "#F0C374", "#BF7A50")
palette_season <- c("#ffe4c4", "#FFDEC4", "#f0f8ff", "#f5fffa")

#((("#228B22","hotpink","#228B22","#228B22","deeppink3","#228B22","#005B96","#005B96","#005B96")))

palette_season <- c('#C5B993','#B7A99A','#A3B5C6','#A8BBA0')
TMPal <- c('#005B96','#8FDC8F','#228B22','#87CEEB')

Text <- c('#444444', '#DADADA','#ffffff', '#0e0e0e' )
Text <- c('black', 'black', 'black', 'black')
Text <- c('white', 'white', 'white', 'white')

palette_Group <- c('#17b92d','#2a57f0','#403b3f')
strip_color_Group<- strip_themed(
  background_x = elem_list_rect(fill =(TMPal),
                                color = "black"),
  text_x = elem_list_text(colour = Text,
                          face = "bold",size=10))

#STRIPS
strip_color_lake<- strip_themed(
  background_x = elem_list_rect(fill =(palette_lake_3T),
                                color = "black"),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=12))

v2strip_color_lake<- strip_themed(
  background_x = elem_list_rect(fill = 'black',
                                color = "black"),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=10))

TMPstrip<- strip_themed(
  background_x = elem_list_rect(fill =(code_saisons),
                                color = "black"),
  text_x = elem_list_text(color = 'white',face = 'bold', size=12))

TMPstrip<- strip_themed(
  background_x = elem_list_rect(fill =(palette_season),
                                color = "black"),
  text_x = elem_list_text(color = 'black',face = 'bold', size=12))

####
month_code <- c('Jun21', 'Jul21', 'Aug21', 'Sep21', 'Oct21', 'Nov21', 'Jan22', 'Fev22', 'Mar22', 
                'Apr22', 'May22', 'Jun22', 'Jul22', 'Aug22', 'Sep22', 'Oct22', 'Nov22', 'Dec22')




