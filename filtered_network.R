#### Filtered network ####
#### VSM  ####

####__VSM 1% data load ####
VSM_euk <- rar_euk %>% subset_samples(lakeID == "VSM")
asv_list <- get_asv_above_threshold(VSM_euk, threshold = 1)
VSM_euk <- prune_taxa(asv_list,VSM_euk )
asv_list <- get_asv_min_samples(VSM_euk, 0.4)
VSM_euk <- prune_taxa(asv_list,VSM_euk)
n_ASV = ntaxa(VSM_euk)
VSM_euk <- tax_glom(VSM_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(VSM_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
VSM_euk_renamed <- renameTaxa(VSM_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
VSM_pears <- netConstruct(VSM_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

####____network properties ####
VSM_props <- netAnalyze(VSM_pears, clustMethod = "cluster_fast_greedy")

VSM_net_prop.df <- cbind(lakeID= "VSM",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(VSM_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(VSM_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(VSM_props)[[2]][2,1],
                         modularity = summary(VSM_props)[[2]][3,1],
                         positive_edge_percentage = summary(VSM_props)[[2]][4,1],
                         natural_connectivity = summary(VSM_props)[[2]][6,1],
                         average_path_length = summary(VSM_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(VSM_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(VSM_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(VSM_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))

palette_mode <- c("#325d3b","#71ab68","#c2ff94",
                  "#f3efa0","#dfd96d","#d2af79",
                  "#bd9568","#a67b58","#8f6349",
                  "#784c3a","#60372c","#49231e",
                  "#321010","#ff4343","#ff807d",
                  "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa
palette_mode <- c("#71ab68","#49231e","#d2af79",
                  "#f3efa0","#dfd96d","#a67b58",
                  "#325d3b","#bd9568")

# palette_mode <- c("#784c3a","#bababa","#8f6349",
#                   "#71ab68","#60372c","#49231e",
#                   "#d2af79","#f3efa0","#dfd96d",
#                   "#ffb3b3","#a67b58","#325d3b",
#                   "#c2ff94","#ff807d","#bd9568",
#                   "#321010","#ff4343")

VSM_network <-plot(VSM_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   #title1 = paste0("VSM ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)
#### | ####

#### JAB  ####

####__JAB 1% data load ####
JAB_euk <- rar_euk %>% subset_samples(lakeID == "JAB")
asv_list <- get_asv_above_threshold(JAB_euk, threshold = 1)
JAB_euk <- prune_taxa(asv_list,JAB_euk )
asv_list <- get_asv_above_threshold(JAB_euk, threshold = 1)
JAB_euk <- prune_taxa(asv_list,JAB_euk )
n_ASV=ntaxa(JAB_euk)
JAB_euk <- tax_glom(JAB_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(JAB_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
JAB_euk_renamed <- renameTaxa(JAB_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
JAB_pears <- netConstruct(JAB_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

####____network properties ####
JAB_props <- netAnalyze(JAB_pears, clustMethod = "cluster_fast_greedy")

JAB_net_prop.df <- cbind(lakeID= "JAB",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(JAB_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(JAB_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(JAB_props)[[2]][2,1],
                         modularity = summary(JAB_props)[[2]][3,1],
                         positive_edge_percentage = summary(JAB_props)[[2]][4,1],
                         natural_connectivity = summary(JAB_props)[[2]][6,1],
                         average_path_length = summary(JAB_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(JAB_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(JAB_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(JAB_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))

# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa

palette_mode <- c("#784c3a","#8f6349","#71ab68",
                  "#49231e","#d2af79","#f3efa0",
                  "#dfd96d","#a67b58","#325d3b",
                  "#c2ff94","#ff807d","#bd9568",
                  "#321010","#ff4343")

# palette_mode <- c("#784c3a","#8f6349","#71ab68",
#                   "#49231e","#d2af79","#f3efa0",
#                   "#dfd96d","#a67b58","#325d3b",
#                   "#c2ff94","#ff807d","#bd9568",
#                   "#321010","#ff4343","grey",
#                   "#8f6349")

JAB_network <-plot(JAB_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("JAB ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)

#### | ####

#### CER-L  ####

####__CER-L 1% data load ####
CER_L_euk <- rar_euk %>% subset_samples(lakeID == "CER-L")
asv_list <- get_asv_above_threshold(CER_L_euk, threshold = 1)
CER_L_euk <- prune_taxa(asv_list,CER_L_euk )
asv_list <- get_asv_above_threshold(CER_L_euk, threshold = 1)
CER_L_euk <- prune_taxa(asv_list,CER_L_euk )
n_ASV=ntaxa(CER_L_euk)
CER_L_euk <- tax_glom(CER_L_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(CER_L_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
CER_L_euk_renamed <- renameTaxa(CER_L_euk, 
                                pat = "<name>", 
                                substPat = "<name>_<subst_name>(<subst_R>)",
                                numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
CER_L_pears <- netConstruct(CER_L_euk_renamed,taxRank = "ASV_ID",
                            zeroMethod = "multRepl", normMethod = "clr",
                            measure = "sparcc", sparsMethod = "threshold",
                            thresh = threshhold,verbose = 3)

####____network properties ####
CER_L_props <- netAnalyze(CER_L_pears, clustMethod = "cluster_fast_greedy")

CER_L_net_prop.df <- cbind(lakeID= "CER-L",
                           n_ASV=n_ASV,
                           nb_positive_edge = nrow(CER_L_pears$edgelist1 %>% subset(asso >0)),
                           nb_negative_edge = nrow(CER_L_pears$edgelist1 %>% subset(asso <0)),
                           clustering_coefficient = summary(CER_L_props)[[2]][2,1],
                           modularity = summary(CER_L_props)[[2]][3,1],
                           positive_edge_percentage = summary(CER_L_props)[[2]][4,1],
                           natural_connectivity = summary(CER_L_props)[[2]][6,1],
                           average_path_length = summary(CER_L_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(CER_L_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(CER_L_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(CER_L_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))
# 
# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa
#Bigyra Centroplasthelida Cercozoa Chlorophyta Choanoflagellata Chrompodellids Ciliophora Cryptophyta Dinoflagellata Fungi Gyrista Haptophyta Ichthyosporea Kathablepharida Perkinsea Telonemia
palette_mode <- c("#784c3a","#bababa","#8f6349",
                  "#71ab68","#60372c","#49231e",
                  "#d2af79","#f3efa0","#dfd96d",
                  "#a67b58","#325d3b","#c2ff94",
                  "#ff807d","#bd9568","#321010",
                  "#ff4343")

CER_L_network <-plot(CER_L_props,layout = "spring",
                     repulsion = 0.8,edgeInvisPar = threshhold,
                     shortenLabels = "none",charToRm = "none",
                     labelScale = FALSE,rmSingles = TRUE,
                     nodeSize = "clr",nodeSizeSpread = 4,
                     nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                     featVecCol = mode, colorVec =  palette_mode,
                     posCol = "#5e3c99", negCol = "#e66101",
                     edgeTranspLow = 0, edgeTranspHigh = 40,
                     cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                     title1 = paste0("CER-L ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                     showTitle = TRUE,cexTitle = 1)

#### | ####

#### CER-S  ####

####__CER-L 1% data load ####
CER_S_euk <- rar_euk %>% subset_samples(lakeID == "CER-S")
asv_list <- get_asv_above_threshold(CER_S_euk, threshold = 1)
CER_S_euk <- prune_taxa(asv_list,CER_S_euk )
asv_list <- get_asv_above_threshold(CER_S_euk, threshold = 1)
CER_S_euk <- prune_taxa(asv_list,CER_S_euk )
n_ASV=ntaxa(CER_S_euk)

CER_S_euk <- tax_glom(CER_S_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(CER_S_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
CER_S_euk_renamed <- renameTaxa(CER_S_euk, 
                                pat = "<name>", 
                                substPat = "<name>_<subst_name>(<subst_R>)",
                                numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
CER_S_pears <- netConstruct(CER_S_euk_renamed,taxRank = "ASV_ID",
                            zeroMethod = "multRepl", normMethod = "clr",
                            measure = "sparcc", sparsMethod = "threshold",
                            thresh = threshhold,verbose = 3)

####____network properties ####
CER_S_props <- netAnalyze(CER_S_pears, clustMethod = "cluster_fast_greedy")

CER_S_net_prop.df <- cbind(lakeID= "CER-S",
                           n_ASV=n_ASV,
                           nb_positive_edge = nrow(CER_S_pears$edgelist1 %>% subset(asso >0)),
                           nb_negative_edge = nrow(CER_S_pears$edgelist1 %>% subset(asso <0)),
                           clustering_coefficient = summary(CER_S_props)[[2]][2,1],
                           modularity = summary(CER_S_props)[[2]][3,1],
                           positive_edge_percentage = summary(CER_S_props)[[2]][4,1],
                           natural_connectivity = summary(CER_S_props)[[2]][6,1],
                           average_path_length = summary(CER_S_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(CER_S_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(CER_S_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(CER_S_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))
# 
# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa

palette_mode <- c("#784c3a","#8f6349","#71ab68",
                  "#60372c","#49231e","#d2af79",
                  "#f3efa0","#dfd96d","#a67b58",
                  "#325d3b","#c2ff94","#ff807d",
                  "#bd9568","#321010")

CER_S_network <-plot(CER_S_props,layout = "spring",
                     repulsion = 0.8,edgeInvisPar = threshhold,
                     shortenLabels = "none",charToRm = "none",
                     labelScale = FALSE,rmSingles = TRUE,
                     nodeSize = "clr",nodeSizeSpread = 4,
                     nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                     featVecCol = mode, colorVec =  palette_mode,
                     posCol = "#5e3c99", negCol = "#e66101",
                     edgeTranspLow = 0, edgeTranspHigh = 40,
                     cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                     title1 = paste0("CER-S ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                     showTitle = TRUE,cexTitle = 1)

#### | ####

#### CRE  ####

####__CRE 1% data load ####
CRE_euk <- rar_euk %>% subset_samples(lakeID == "CRE")
asv_list <- get_asv_above_threshold(CRE_euk, threshold = 1)
CRE_euk <- prune_taxa(asv_list,CRE_euk )
asv_list <- get_asv_above_threshold(CRE_euk, threshold = 1)
CRE_euk <- prune_taxa(asv_list,CRE_euk )
n_ASV=ntaxa(CRE_euk)
CRE_euk <- tax_glom(CRE_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(CRE_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
CRE_euk_renamed <- renameTaxa(CRE_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
CRE_pears <- netConstruct(CRE_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

####____network properties ####
CRE_props <- netAnalyze(CRE_pears, clustMethod = "cluster_fast_greedy")

CRE_net_prop.df <- cbind(lakeID= "CRE",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(CRE_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(CRE_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(CRE_props)[[2]][2,1],
                         modularity = summary(CRE_props)[[2]][3,1],
                         positive_edge_percentage = summary(CRE_props)[[2]][4,1],
                         natural_connectivity = summary(CRE_props)[[2]][6,1],
                         average_path_length = summary(CRE_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(CRE_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(CRE_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(CRE_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))

# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# # #Discoba ffb3b3
# # #Centroplasthelida bababa
palette_mode <- c("#784c3a","#8f6349","#71ab68",
                  "#60372c","#49231e","#d2af79",
                  "#f3efa0","#dfd96d","#a67b58",
                  "#325d3b","#bd9568","#FFFFFF",
                  "#321010")

CRE_network <-plot(CRE_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("CRE ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)

#### | ####

#### BLR  ####

####__BLR 1% data load ####
BLR_euk <- rar_euk %>% subset_samples(lakeID == "BLR")
asv_list <- get_asv_above_threshold(BLR_euk, threshold = 1)
BLR_euk <- prune_taxa(asv_list,BLR_euk)
asv_list <- get_asv_min_samples(BLR_euk, 0.4)
BLR_euk <- prune_taxa(asv_list,BLR_euk)
n_ASV = ntaxa(BLR_euk)
BLR_euk <- tax_glom(BLR_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(BLR_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
BLR_euk_renamed <- renameTaxa(BLR_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshold = 0.6
BLR_pears <- netConstruct(BLR_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc",
                          sparsMethod = "threshold",
                          thresh = threshold,verbose = 3)

####____network properties ####
BLR_props <- netAnalyze(BLR_pears, clustMethod = "cluster_fast_greedy")

BLR_net_prop.df <- cbind(lakeID= "BLR",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(BLR_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(BLR_pears$edgelist1 %>% subset(asso <0)),
                         positive_edge_percentage = summary(BLR_props)[[2]][4,1],
                         clustering_coefficient = summary(BLR_props)[[2]][2,1],
                         modularity = summary(BLR_props)[[2]][3,1],
                         natural_connectivity = summary(BLR_props)[[2]][6,1],
                         average_path_length = summary(BLR_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(BLR_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(BLR_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(BLR_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))
# 
# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# # #Discoba ffb3b3
# # #Centroplasthelida bababa
##Alveolata BF3232

palette_mode <- c("#8f6349","#71ab68",
                  "#d2af79","#f3efa0","#dfd96d",
                  "#a67b58","#325d3b","#bd9568")

# palette_mode <- c("#784c3a","#bababa","#8f6349",
#                   "#71ab68","#60372c","#49231e",
#                   "#d2af79","#f3efa0","#dfd96d",
#                   "#a67b58","#325d3b","#bd9568",
#                   "#FFFFFF","#321010","#ff4343")
# 
# palette_mode <- c("#BF3232","#784c3a","#bababa","#8f6349",
#                   "#71ab68","#60372c","#49231e",
#                   "#d2af79","#f3efa0","#dfd96d",
#                   "#a67b58","#325d3b","#c2ff94",
#                   "#ff807d","#bd9568","#FFFFFF","#321010")

BLR_network <-plot(BLR_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("BLR ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)

#### | ####

#### LGP  ####

####__LGP 1% data load ####
LGP_euk <- rar_euk %>% subset_samples(lakeID == "LGP")
asv_list <- get_asv_above_threshold(LGP_euk, threshold = 1)
LGP_euk <- prune_taxa(asv_list,LGP_euk )
asv_list <- get_asv_min_samples(LGP_euk, 0.4)
LGP_euk <- prune_taxa(asv_list,LGP_euk)
n_ASV = ntaxa(LGP_euk)

LGP_euk <- tax_glom(LGP_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(LGP_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
LGP_euk_renamed <- renameTaxa(LGP_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
LGP_pears <- netConstruct(LGP_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

####____network properties ####
LGP_props <- netAnalyze(LGP_pears, clustMethod = "cluster_fast_greedy")

LGP_net_prop.df <- cbind(lakeID= "LGP",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(LGP_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(LGP_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(LGP_props)[[2]][2,1],
                         modularity = summary(LGP_props)[[2]][3,1],
                         positive_edge_percentage = summary(LGP_props)[[2]][4,1],
                         natural_connectivity = summary(LGP_props)[[2]][6,1],
                         average_path_length = summary(LGP_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(LGP_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(LGP_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(LGP_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))
#
# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa

palette_mode <- c("#71ab68","#49231e","#d2af79",
                  "#f3efa0","#325d3b","#bd9568")

LGP_network <-plot(LGP_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("LGP ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)


#### CSM  ####

####__CSM 1% data load ####
CSM_euk <- rar_euk %>% subset_samples(lakeID == "CSM")
asv_list <- get_asv_above_threshold(CSM_euk, threshold = 1)
CSM_euk <- prune_taxa(asv_list,CSM_euk )
asv_list <- get_asv_min_samples(CSM_euk, 0.4)
CSM_euk <- prune_taxa(asv_list,CSM_euk)
n_ASV = ntaxa(CSM_euk)

CSM_euk <- tax_glom(CSM_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(CSM_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
CSM_euk_renamed <- renameTaxa(CSM_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshhold = 0.6
CSM_pears <- netConstruct(CSM_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

####____network properties ####
CSM_props <- netAnalyze(CSM_pears, clustMethod = "cluster_fast_greedy")

CSM_net_prop.df <- cbind(lakeID= "CSM",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(CSM_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(CSM_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(CSM_props)[[2]][2,1],
                         modularity = summary(CSM_props)[[2]][3,1],
                         positive_edge_percentage = summary(CSM_props)[[2]][4,1],
                         natural_connectivity = summary(CSM_props)[[2]][6,1],
                         average_path_length = summary(CSM_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(CSM_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(CSM_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(CSM_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))
#
# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# #Discoba ffb3b3
# #Centroplasthelida bababa
palette_mode <- c("#71ab68",
                  "#49231e","#d2af79","#f3efa0","#325d3b",
                  "#bd9568")

CSM_network <-plot(CSM_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("CSM ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)

#### | ####

#### VSS  ####

####__VSS 1% data load ####
VSS_euk <- rar_euk %>% subset_samples(lakeID == "VSS")
asv_list <- get_asv_above_threshold(VSS_euk, threshold = 1)
VSS_euk <- prune_taxa(asv_list,VSS_euk )
asv_list <- get_asv_min_samples(VSS_euk, 0.4)
VSS_euk <- prune_taxa(asv_list,VSS_euk)
n_ASV = ntaxa(VSS_euk)
VSS_euk <- tax_glom(VSS_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(VSS_euk), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
VSS_euk_renamed <- renameTaxa(VSS_euk, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "ASV_ID")

####____create network ####
threshold = 0.6
VSS_pears <- netConstruct(VSS_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshold,verbose = 3)

####____network properties ####
VSS_props <- netAnalyze(VSS_pears, clustMethod = "cluster_fast_greedy")

VSS_net_prop.df <- cbind(lakeID= "VSS",
                         n_ASV=n_ASV,
                         nb_positive_edge = nrow(VSS_pears$edgelist1 %>% subset(asso >0)),
                         nb_negative_edge = nrow(VSS_pears$edgelist1 %>% subset(asso <0)),
                         clustering_coefficient = summary(VSS_props)[[2]][2,1],
                         modularity = summary(VSS_props)[[2]][3,1],
                         positive_edge_percentage = summary(VSS_props)[[2]][4,1],
                         natural_connectivity = summary(VSS_props)[[2]][6,1],
                         average_path_length = summary(VSS_props)[[2]][10,1])

####____plot network ####
graph3 <- igraph::graph_from_adjacency_matrix(VSS_pears$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(VSS_pears$adjaMat1)

####____Trophic mode color
taxtab <- as(tax_table(VSS_euk_renamed), "matrix")
mode <- as.factor(taxtab[, "trophic_subdivision"])
names(mode) <- taxtab[, "ASV_ID"]

unique(mode)
# levels=c("Gyrista","Chlorophyta",'Haptophyta',
#          "Cryptophyta","Dinoflagellata","Ciliophora",
#          "Kathablepharida","Fungi","Cercozoa",
#          "Bigyra","Choanoflagellata", "Chrompodellids",
#          "Perkinsea","Telonemia",'Ichthyosporea',
#          'NA'))


# palette_mode <- c("#325d3b","#71ab68","#c2ff94",
#                   "#f3efa0","#dfd96d","#d2af79",
#                   "#bd9568","#a67b58","#8f6349",
#                   "#784c3a","#60372c","#49231e",
#                   "#321010","#ff4343","#ff807d",
#                   "#000000")
# # #Discoba ffb3b3
# # #Centroplasthelida bababa

palette_mode <- c("#71ab68","#49231e","#d2af79",
                  "#f3efa0",
                  "#325d3b","#bd9568")

VSS_network <-plot(VSS_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   posCol = "#5e3c99", negCol = "#e66101",
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("VSS ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = TRUE,cexTitle = 1)

#### | ####

#### Properties ####

####__prop df ####
network_prop.df <-
  rbind(VSM_net_prop.df,JAB_net_prop.df,CER_L_net_prop.df,
        CER_S_net_prop.df,CRE_net_prop.df,BLR_net_prop.df,
        LGP_net_prop.df,CSM_net_prop.df,VSS_net_prop.df) %>%
  as.data.frame()

write.csv(network_prop.df,"/home/jesus/Pierre/PhD/ms_Euk/network_prop.csv",row.names = F)

network_props.df <- read_delim("SparCC_network_prop.csv",delim = ";",show_col_types = F)

####__PCA ####
network.pca <- network_props.df %>%
  remove_rownames() %>%
  column_to_rownames("lakeID") %>% scale(.,center = T, scale = T)

network.dist <-network.pca %>% vegdist(., method="euclidean") 

PCA_network.mds<- cmdscale(network.dist,eig=TRUE, k=2)

PCA_network.df <- PCA_network.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,.before=1,
                lakeID=factor(rownames(.),levels=c("VSM", "JAB", "CER-L",
                                                   "CER-S", "CRE", "BLR","LGP", "CSM", 
                                                   "VSS")))
network_prop <- envfit(PCA_network.mds$points ~
                         n_ASV+
                         nb_edge+
                         clustering_coefficient+modularity+
                         positive_edge_percentage+
                         negative_edge_percentage+natural_connectivity+
                         average_path_length,
                       data = as.data.frame(network.pca), perm = 999)

network_vectors.df <- as.data.frame(scores(network_prop, display = "vectors")) 
network_vectors.df <- cbind(network_vectors.df, prop = rownames(network_vectors.df)) 

PCA_network <- PCA_network.df %>%
  ggplot(.,aes(dim1,dim2,color=lakeID,label=lakeID)) +
  #geom_point(size=2,shape=21,show.legend = F, color= "black")+
  geom_segment(data = network_vectors.df,
               mapping=aes(x = 0, xend = Dim1*2.5, y = 0, yend = Dim2*2.5),
               color = "gray40",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  geom_text(data = network_vectors.df,
            mapping=aes(x = Dim1*2.62, y = Dim2*2.62, label = prop),
            size = 3.5,color="gray40",fontface = "bold",inherit.aes = F)+
  geom_text(size=5, check_overlap = T,show.legend = F,fontface="bold")+
  #geom_label(size=4,show.legend = F,fontface="bold", fill = "lightgrey")+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10, color="black"),
        axis.ticks = element_line(color="black"),
        plot.title = element_text(size=11,face = "italic"),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  #scale_fill_manual(values=rev(palette_lake_chla))+
  scale_color_manual(values=rev(palette_lake_chla))+
  #scale_x_reverse()+
  #scale_y_reverse()+
  labs(x=paste0("PC1 [",round(PCA_network.mds$eig[1]*100/sum(PCA_network.mds$eig),1),"%]"),
       y=paste0("PC2 [",round(PCA_network.mds$eig[2]*100/sum(PCA_network.mds$eig),1),"%]"))
PCA_network   

ggsave("/Users/piefouca/Desktop/µEuk/Figures/PCA_network_prop.svg",
       units = "in",dpi = 300,width = 13.4,height = 9.8)

####________________________####