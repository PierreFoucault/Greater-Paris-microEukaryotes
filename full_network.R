#### Full network ####
#plot 626 * 595 pixels
####__VSM 1% data load ####
VSM_euk <- rar_euk %>% subset_samples(lakeID == "VSM")
asv_list <- get_asv_above_threshold(VSM_euk, threshold = 1)
VSM_euk <- prune_taxa(asv_list,VSM_euk)
n_ASV = ntaxa(VSM_euk)
VSM_euk <- tax_glom(VSM_euk, taxrank = "ASV_ID")

# Taxonomic table
taxtab <- as(tax_table(VSM_euk), "matrix")

# Rename taxonomic table and make ASV rank unique
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

VSM_props <- netAnalyze(VSM_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  "#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  "#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  '#E1FFCA', #Haptophyta
                  "#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                 #"#FFFFFF", #NA
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343', #Perkinsea
                  "#321010") #Telonemia

VSM_network <-plot(VSM_props,layout = "spring",
                   repulsion = 0.85,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 0,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("VSM ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
VSM_trophic_network.df <-
  data.frame(ASV=rownames(VSM_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(VSM_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(VSM_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(VSM_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(VSM_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  VSM_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(VSM_trophic_network.df))

VSM_trophic_cluster.df <-
  as.data.frame(VSM_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="VSM_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(VSM_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(VSM_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(VSM_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(VSM_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(VSM_trophic_cluster.df %>% filter(trophic_mode == 1))
trophic_mode_2 <- nrow(VSM_trophic_cluster.df %>% filter(trophic_mode == 2))
trophic_mode_3 <- nrow(VSM_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(VSM_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(VSM_trophic_cluster.df$trophic_mode)

VSM_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

VSM_trophic_edge.df <-
  as.data.frame(VSM_pears$edgelist1 %>%
       subset(v1 %in% VSM_trophic_network.df[,1]
              & v2 %in% VSM_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(VSM_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(VSM_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(VSM_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  VSM_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()
  
VSM_net_prop.df <- data.frame(lakeID= "VSM",
                         nb_ASV=n_ASV,
                         nb_components = summary(VSM_props)[[1]][1,1],
                         nb_clusters = length(unique(VSM_trophic_cluster.df$cluster_ID)),
                         nb_connected_node = nrow(VSM_trophic_network.df),
                         nb_edge = (nrow(VSM_pears$edgelist1 %>%
                                           subset(asso >0 & v1 %in% VSM_trophic_network.df[,1]
                                                          & v2 %in% VSM_trophic_network.df[,1]))
                                    + nrow(VSM_pears$edgelist1 %>%
                                             subset(asso <0 & v1 %in% VSM_trophic_network.df[,1]
                                                            & v2 %in% VSM_trophic_network.df[,1]))),
                         nb_positive_edge = nrow(VSM_pears$edgelist1 %>%
                                                   subset(asso >0 & v1 %in% VSM_trophic_network.df[,1]
                                                                  & v2 %in% VSM_trophic_network.df[,1])),
                         nb_negative_edge = nrow(VSM_pears$edgelist1 %>%
                                                   subset(asso <0 & v1 %in% VSM_trophic_network.df[,1]
                                                                  & v2 %in% VSM_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(VSM_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(VSM_props)[[1]][2,1],
                modularity = summary(VSM_props)[[1]][3,1],
                edge_density = summary(VSM_props)[[1]][5,1],
                natural_connectivity = summary(VSM_props)[[1]][6,1],
                percentage_nodes_LCC = summary(VSM_props)[[2]][1,1],
                LCC_average_path_length = summary(VSM_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)

#### | ####

####__JAB 1% data load ####
JAB_euk <- rar_euk %>% subset_samples(lakeID == "JAB")
asv_list <- get_asv_above_threshold(JAB_euk, threshold = 1)
JAB_euk <- prune_taxa(asv_list,JAB_euk)
n_ASV = ntaxa(JAB_euk)
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

JAB_props <- netAnalyze(JAB_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  #'#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  '#E1FFCA', #Haptophyta
                  "#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  #"#FFFFFF", #NA
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343', #Perkinsea
                  "#321010") #Telonemia

JAB_network <-plot(JAB_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("JAB ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
JAB_trophic_network.df <-
  data.frame(ASV=rownames(JAB_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(JAB_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(JAB_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(JAB_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(JAB_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  JAB_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(JAB_trophic_network.df))

JAB_trophic_cluster.df <-
  as.data.frame(JAB_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="JAB_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(JAB_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(JAB_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(JAB_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(JAB_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(JAB_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(JAB_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(JAB_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(JAB_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(JAB_trophic_cluster.df$trophic_mode)

JAB_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

JAB_trophic_edge.df <-
  as.data.frame(JAB_pears$edgelist1 %>%
                  subset(v1 %in% JAB_trophic_network.df[,1]
                         & v2 %in% JAB_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(JAB_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(JAB_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(JAB_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>%   dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  JAB_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

JAB_net_prop.df <- data.frame(lakeID= "JAB",
                              nb_ASV=n_ASV,
                              nb_components = summary(JAB_props)[[1]][1,1],
                              nb_clusters = length(unique(JAB_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(JAB_trophic_network.df),
                              nb_edge = (nrow(JAB_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% JAB_trophic_network.df[,1]
                                                       & v2 %in% JAB_trophic_network.df[,1]))
                                         + nrow(JAB_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% JAB_trophic_network.df[,1]
                                                         & v2 %in% JAB_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(JAB_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% JAB_trophic_network.df[,1]
                                                               & v2 %in% JAB_trophic_network.df[,1])),
                              nb_negative_edge = nrow(JAB_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% JAB_trophic_network.df[,1]
                                                               & v2 %in% JAB_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(JAB_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(JAB_props)[[1]][2,1],
                modularity = summary(JAB_props)[[1]][3,1],
                edge_density = summary(JAB_props)[[1]][5,1],
                natural_connectivity = summary(JAB_props)[[1]][6,1],
                percentage_nodes_LCC = summary(JAB_props)[[2]][1,1],
                LCC_average_path_length = summary(JAB_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)

#### | ####

####__CER-L 1% data load ####
CER_L_euk <- rar_euk %>% subset_samples(lakeID == "CER-L")
asv_list <- get_asv_above_threshold(CER_L_euk, threshold = 1)
CER_L_euk <- prune_taxa(asv_list,CER_L_euk)
n_ASV = ntaxa(CER_L_euk)
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

palette_mode <- c("#8B3D0D", #Bigyra
                  "#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  '#E1FFCA', #Haptophyta
                  "#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  #"#FFFFFF", #NA
                  #"#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343', #Perkinsea
                  "#321010") #Telonemia

CER_L_network <-plot(CER_L_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("CER-L ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
CER_L_trophic_network.df <-
  data.frame(ASV=rownames(CER_L_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CER_L_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CER_L_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CER_L_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CER_L_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  CER_L_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(CER_L_trophic_network.df))

CER_L_trophic_cluster.df <-
  as.data.frame(CER_L_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="CER_L_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CER_L_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CER_L_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CER_L_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CER_L_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(CER_L_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(CER_L_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(CER_L_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(CER_L_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(CER_L_trophic_cluster.df$trophic_mode)

CER_L_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

CER_L_trophic_edge.df <-
  as.data.frame(CER_L_pears$edgelist1 %>%
                  subset(v1 %in% CER_L_trophic_network.df[,1]
                         & v2 %in% CER_L_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(CER_L_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(CER_L_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(CER_L_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  CER_L_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

CER_L_net_prop.df <- data.frame(lakeID= "CER_L",
                              nb_ASV=n_ASV,
                              nb_components = summary(CER_L_props)[[1]][1,1],
                              nb_clusters = length(unique(CER_L_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(CER_L_trophic_network.df),
                              nb_edge = (nrow(CER_L_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% CER_L_trophic_network.df[,1]
                                                       & v2 %in% CER_L_trophic_network.df[,1]))
                                         + nrow(CER_L_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% CER_L_trophic_network.df[,1]
                                                         & v2 %in% CER_L_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(CER_L_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% CER_L_trophic_network.df[,1]
                                                               & v2 %in% CER_L_trophic_network.df[,1])),
                              nb_negative_edge = nrow(CER_L_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% CER_L_trophic_network.df[,1]
                                                               & v2 %in% CER_L_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(CER_L_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(CER_L_props)[[1]][2,1],
                modularity = summary(CER_L_props)[[1]][3,1],
                edge_density = summary(CER_L_props)[[1]][5,1],
                natural_connectivity = summary(CER_L_props)[[1]][6,1],
                percentage_nodes_LCC = summary(CER_L_props)[[2]][1,1],
                LCC_average_path_length = summary(CER_L_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)

#### | ####

####__CER-S 1% data load ####
CER_S_euk <- rar_euk %>% subset_samples(lakeID == "CER-S")
asv_list <- get_asv_above_threshold(CER_S_euk, threshold = 1)
CER_S_euk <- prune_taxa(asv_list,CER_S_euk)
n_ASV = ntaxa(CER_S_euk)
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

CER_S_props <- netAnalyze(CER_S_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  '#E1FFCA', #Haptophyta
                  "#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  #"#FFFFFF", #NA
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343' #Perkinsea
                  #"#321010"
                  ) #Telonemia

CER_S_network <-plot(CER_S_props,layout = "spring",
                     repulsion = 0.75,edgeInvisPar = threshhold,
                     shortenLabels = "none",charToRm = "none",
                     labelScale = FALSE,rmSingles = TRUE,
                     nodeSize = "clr",nodeSizeSpread = 4,
                     nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                     featVecCol = mode, colorVec =  palette_mode,
                     negCol = "#09264A", posCol = "#f99820",
                     edgeWidth = 1.2,
                     edgeTranspLow = 0, edgeTranspHigh = 40,
                     cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                     title1 = paste0("CER-S ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                     showTitle = F,cexTitle = 1)

####____network properties ####
CER_S_trophic_network.df <-
  data.frame(ASV=rownames(CER_S_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CER_S_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CER_S_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CER_S_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CER_S_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  CER_S_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(CER_S_trophic_network.df))

CER_S_trophic_cluster.df <-
  as.data.frame(CER_S_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="CER_S_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CER_S_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CER_S_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CER_S_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CER_S_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(CER_S_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(CER_S_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(CER_S_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(CER_S_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(CER_S_trophic_cluster.df$trophic_mode)

CER_S_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

CER_S_trophic_edge.df <-
  as.data.frame(CER_S_pears$edgelist1 %>%
                  subset(v1 %in% CER_S_trophic_network.df[,1]
                         & v2 %in% CER_S_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(CER_S_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(CER_S_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(CER_S_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  CER_S_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

CER_S_net_prop.df <- data.frame(lakeID= "CER_S",
                              nb_ASV=n_ASV,
                              nb_components = summary(CER_S_props)[[1]][1,1],
                              nb_clusters = length(unique(CER_S_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(CER_S_trophic_network.df),
                              nb_edge = (nrow(CER_S_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% CER_S_trophic_network.df[,1]
                                                       & v2 %in% CER_S_trophic_network.df[,1]))
                                         + nrow(CER_S_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% CER_S_trophic_network.df[,1]
                                                         & v2 %in% CER_S_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(CER_S_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% CER_S_trophic_network.df[,1]
                                                               & v2 %in% CER_S_trophic_network.df[,1])),
                              nb_negative_edge = nrow(CER_S_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% CER_S_trophic_network.df[,1]
                                                               & v2 %in% CER_S_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(CER_S_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(CER_S_props)[[1]][2,1],
                modularity = summary(CER_S_props)[[1]][3,1],
                edge_density = summary(CER_S_props)[[1]][5,1],
                natural_connectivity = summary(CER_S_props)[[1]][6,1],
                percentage_nodes_LCC = summary(CER_S_props)[[2]][1,1],
                LCC_average_path_length = summary(CER_S_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

####__CRE 1% data load ####
CRE_euk <- rar_euk %>% subset_samples(lakeID == "CRE")
asv_list <- get_asv_above_threshold(CRE_euk, threshold = 1)
CRE_euk <- prune_taxa(asv_list,CRE_euk)
n_ASV = ntaxa(CRE_euk)
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

CRE_props <- netAnalyze(CRE_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  #'#E1FFCA', #Haptophyta
                  #"#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  "#FFFFFF", #not_annotated
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343' #Perkinsea
                  #"#321010"
) #Telonemia
CRE_network <-plot(CRE_props,layout = "spring",
                   repulsion = 0.75,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("CRE ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
CRE_trophic_network.df <-
  data.frame(ASV=rownames(CRE_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CRE_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CRE_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CRE_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CRE_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  CRE_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(CRE_trophic_network.df))

CRE_trophic_cluster.df <-
  as.data.frame(CRE_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="CRE_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CRE_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CRE_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CRE_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CRE_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(CRE_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(CRE_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(CRE_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(CRE_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(CRE_trophic_cluster.df$trophic_mode)

CRE_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

CRE_trophic_edge.df <-
  as.data.frame(CRE_pears$edgelist1 %>%
                  subset(v1 %in% CRE_trophic_network.df[,1]
                         & v2 %in% CRE_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(CRE_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(CRE_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(CRE_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  CRE_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

CRE_net_prop.df <- data.frame(lakeID= "CRE",
                              nb_ASV=n_ASV,
                              nb_components = summary(CRE_props)[[1]][1,1],
                              nb_clusters = length(unique(CRE_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(CRE_trophic_network.df),
                              nb_edge = (nrow(CRE_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% CRE_trophic_network.df[,1]
                                                       & v2 %in% CRE_trophic_network.df[,1]))
                                         + nrow(CRE_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% CRE_trophic_network.df[,1]
                                                         & v2 %in% CRE_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(CRE_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% CRE_trophic_network.df[,1]
                                                               & v2 %in% CRE_trophic_network.df[,1])),
                              nb_negative_edge = nrow(CRE_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% CRE_trophic_network.df[,1]
                                                               & v2 %in% CRE_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(CRE_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(CRE_props)[[1]][2,1],
                modularity = summary(CRE_props)[[1]][3,1],
                edge_density = summary(CRE_props)[[1]][5,1],
                natural_connectivity = summary(CRE_props)[[1]][6,1],
                percentage_nodes_LCC = summary(CRE_props)[[2]][1,1],
                LCC_average_path_length = summary(CRE_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

####__BLR 1% data load ####
BLR_euk <- rar_euk %>% subset_samples(lakeID == "BLR")
asv_list <- get_asv_above_threshold(BLR_euk, threshold = 1)
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
threshhold = 0.6
BLR_pears <- netConstruct(BLR_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

BLR_props <- netAnalyze(BLR_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  "#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  #"#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  #'#E1FFCA', #Haptophyta
                  #"#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  "#FFFFFF", #not_annotated
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343', #Perkinsea
                  "#321010") #Telonemia

BLR_network <-plot(BLR_props,layout = "spring",
                   repulsion = 0.9,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("BLR ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
BLR_trophic_network.df <-
  data.frame(ASV=rownames(BLR_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(BLR_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(BLR_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(BLR_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(BLR_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  BLR_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(BLR_trophic_network.df))

BLR_trophic_cluster.df <-
  as.data.frame(BLR_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="BLR_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(BLR_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(BLR_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(BLR_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(BLR_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(BLR_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(BLR_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(BLR_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(BLR_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(BLR_trophic_cluster.df$trophic_mode)

BLR_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

BLR_trophic_edge.df <-
  as.data.frame(BLR_pears$edgelist1 %>%
                  subset(v1 %in% BLR_trophic_network.df[,1]
                         & v2 %in% BLR_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(BLR_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(BLR_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(BLR_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))
  

nb_different_trophic_edge <-
  BLR_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

BLR_net_prop.df <- data.frame(lakeID= "BLR",
                              nb_ASV=n_ASV,
                              nb_components = summary(BLR_props)[[1]][1,1],
                              nb_clusters = length(unique(BLR_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(BLR_trophic_network.df),
                              nb_edge = (nrow(BLR_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% BLR_trophic_network.df[,1]
                                                       & v2 %in% BLR_trophic_network.df[,1]))
                                         + nrow(BLR_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% BLR_trophic_network.df[,1]
                                                         & v2 %in% BLR_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(BLR_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% BLR_trophic_network.df[,1]
                                                               & v2 %in% BLR_trophic_network.df[,1])),
                              nb_negative_edge = nrow(BLR_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% BLR_trophic_network.df[,1]
                                                               & v2 %in% BLR_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(BLR_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(BLR_props)[[1]][2,1],
                modularity = summary(BLR_props)[[1]][3,1],
                edge_density = summary(BLR_props)[[1]][5,1],
                natural_connectivity = summary(BLR_props)[[1]][6,1],
                percentage_nodes_LCC = summary(BLR_props)[[2]][1,1],
                LCC_average_path_length = summary(BLR_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

####__LGP 1% data load ####
LGP_euk <- rar_euk %>% subset_samples(lakeID == "LGP")
asv_list <- get_asv_above_threshold(LGP_euk, threshold = 1)
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

LGP_props <- netAnalyze(LGP_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  #'#E1FFCA', #Haptophyta
                  #"#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  #"#FFFFFF", #not_annotated
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343' #Perkinsea
                  #"#321010"
) #Telonemia

LGP_network <-plot(LGP_props,layout = "spring",
                   repulsion = 0.75,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("LGP ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
LGP_trophic_network.df <-
  data.frame(ASV=rownames(LGP_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(LGP_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(LGP_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(LGP_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(LGP_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))
  

percent_trophic_mode <-
  LGP_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(LGP_trophic_network.df))

LGP_trophic_cluster.df <-
  as.data.frame(LGP_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="LGP_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(LGP_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(LGP_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(LGP_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(LGP_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(LGP_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(LGP_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(LGP_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(LGP_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(LGP_trophic_cluster.df$trophic_mode)

LGP_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

LGP_trophic_edge.df <-
  as.data.frame(LGP_pears$edgelist1 %>%
                  subset(v1 %in% LGP_trophic_network.df[,1]
                         & v2 %in% LGP_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(LGP_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(LGP_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(LGP_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  LGP_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

LGP_net_prop.df <- data.frame(lakeID= "LGP",
                              nb_ASV=n_ASV,
                              nb_components = summary(LGP_props)[[1]][1,1],
                              nb_clusters = length(unique(LGP_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(LGP_trophic_network.df),
                              nb_edge = (nrow(LGP_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% LGP_trophic_network.df[,1]
                                                       & v2 %in% LGP_trophic_network.df[,1]))
                                         + nrow(LGP_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% LGP_trophic_network.df[,1]
                                                         & v2 %in% LGP_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(LGP_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% LGP_trophic_network.df[,1]
                                                               & v2 %in% LGP_trophic_network.df[,1])),
                              nb_negative_edge = nrow(LGP_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% LGP_trophic_network.df[,1]
                                                               & v2 %in% LGP_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(LGP_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(LGP_props)[[1]][2,1],
                modularity = summary(LGP_props)[[1]][3,1],
                edge_density = summary(LGP_props)[[1]][5,1],
                natural_connectivity = summary(LGP_props)[[1]][6,1],
                percentage_nodes_LCC = summary(LGP_props)[[2]][1,1],
                LCC_average_path_length = summary(LGP_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

####__CSM 1% data load ####
CSM_euk <- rar_euk %>% subset_samples(lakeID == "CSM")
asv_list <- get_asv_above_threshold(CSM_euk, threshold = 1)
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

CSM_props <- netAnalyze(CSM_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  #'#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  #'#E1FFCA', #Haptophyta
                  "#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  "#FFFFFF", #not_annotated
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343' #Perkinsea
                  #"#321010"
) #Telonemia

CSM_network <-plot(CSM_props,layout = "spring",
                   repulsion = 0.75,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("CSM ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
CSM_trophic_network.df <-
  data.frame(ASV=rownames(CSM_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CSM_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CSM_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CSM_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CSM_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))


percent_trophic_mode <-
  CSM_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(CSM_trophic_network.df))

CSM_trophic_cluster.df <-
  as.data.frame(CSM_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="CSM_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(CSM_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(CSM_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(CSM_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(CSM_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(CSM_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(CSM_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(CSM_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(CSM_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(CSM_trophic_cluster.df$trophic_mode)

CSM_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

CSM_trophic_edge.df <-
  as.data.frame(CSM_pears$edgelist1 %>%
                  subset(v1 %in% CSM_trophic_network.df[,1]
                         & v2 %in% CSM_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(CSM_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(CSM_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(CSM_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  CSM_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

CSM_net_prop.df <- data.frame(lakeID= "CSM",
                              nb_ASV=n_ASV,
                              nb_components = summary(CSM_props)[[1]][1,1],
                              nb_clusters = length(unique(CSM_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(CSM_trophic_network.df),
                              nb_edge = (nrow(CSM_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% CSM_trophic_network.df[,1]
                                                       & v2 %in% CSM_trophic_network.df[,1]))
                                         + nrow(CSM_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% CSM_trophic_network.df[,1]
                                                         & v2 %in% CSM_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(CSM_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% CSM_trophic_network.df[,1]
                                                               & v2 %in% CSM_trophic_network.df[,1])),
                              nb_negative_edge = nrow(CSM_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% CSM_trophic_network.df[,1]
                                                               & v2 %in% CSM_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(CSM_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(CSM_props)[[1]][2,1],
                modularity = summary(CSM_props)[[1]][3,1],
                edge_density = summary(CSM_props)[[1]][5,1],
                natural_connectivity = summary(CSM_props)[[1]][6,1],
                percentage_nodes_LCC = summary(CSM_props)[[2]][1,1],
                LCC_average_path_length = summary(CSM_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

####__VSS 1% data load ####
VSS_euk <- rar_euk %>% subset_samples(lakeID == "VSS")
asv_list <- get_asv_above_threshold(VSS_euk, threshold = 1)
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
threshhold = 0.6
VSS_pears <- netConstruct(VSS_euk_renamed,taxRank = "ASV_ID",
                          zeroMethod = "multRepl", normMethod = "clr",
                          measure = "sparcc", sparsMethod = "threshold",
                          thresh = threshhold,verbose = 3)

VSS_props <- netAnalyze(VSS_pears, clustMethod = "cluster_fast_greedy")

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

palette_mode <- c("#8B3D0D", #Bigyra
                  #"#B89579", #Centroplasthelida
                  "#d15c14", #Cercozoa
                  "#71ab68", #Chlorophyta
                  '#60372c', #Choanoflagellata
                  '#49231e', #Chrompodellids
                  '#d2af79', #Ciliophora
                  '#f3efa0', #Cryptophyta
                  '#dfd96d', #Dinoflagellata
                  #"#ffb3b3", #Discoba
                  "#CF8A4F", #Fungi
                  "#809c13", #Gyrista · diatoms
                  '#2C5234', #Gyrista · non-diatoms
                  #'#E1FFCA', #Haptophyta
                  #"#B52239", #Ichthyosporea
                  '#a67b58', #Kathablepharida
                  #"#FFFFFF", #not_annotated
                  "#E40000", #Parasitic Fungi
                  #"#FF6565", #Parasitic Gyrista
                  '#ff4343' #Perkinsea
                  #"#321010"
) #Telonemia

VSS_network <-plot(VSS_props,layout = "spring",
                   repulsion = 0.8,edgeInvisPar = threshhold,
                   shortenLabels = "none",charToRm = "none",
                   labelScale = FALSE,rmSingles = TRUE,
                   nodeSize = "clr",nodeSizeSpread = 4,
                   nodeColor = "feature",borderCol = "black",nodeTransp = 0,
                   featVecCol = mode, colorVec =  palette_mode,
                   negCol = "#09264A", posCol = "#f99820",
                   edgeWidth = 1.2,
                   edgeTranspLow = 0, edgeTranspHigh = 40,
                   cexNodes = 2, cexLabels = 0, cexHubLabels = 0,
                   title1 = paste0("VSS ASVs Network\n(SparCC - min. corr. ",threshhold,")"),
                   showTitle = F,cexTitle = 1)

####____network properties ####
VSS_trophic_network.df <-
  data.frame(ASV=rownames(VSS_network$q1$layout)) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(VSS_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(VSS_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(VSS_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(VSS_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode))

percent_trophic_mode <-
  VSS_trophic_network.df %>%
  dplyr::mutate(trophic_mode = factor(trophic_mode,
                                      levels = c("phototrophs",
                                                 "mixotrophs",
                                                 "phagotrophs",
                                                 "parasites"))) %>%
  group_by(trophic_mode) %>%
  summarise(percent = (n()*100)/nrow(VSS_trophic_network.df))

VSS_trophic_cluster.df <-
  as.data.frame(VSS_props$clustering$clust1) %>%
  rownames_to_column(var = "ASV_ID") %>%
  rename("cluster_ID"="VSS_props$clustering$clust1") %>%
  filter(cluster_ID > 0) %>%
  dplyr::mutate(.,
                trophic_mode =
                  as.data.frame(VSS_euk@tax_table)$trophic_mode[match(.$ASV,
                                                                      as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                trophic_subdivision =
                  as.data.frame(VSS_euk@tax_table)$trophic_subdivision[match(.$ASV,
                                                                             as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                class =
                  as.data.frame(VSS_euk@tax_table)$Class[match(.$ASV,
                                                               as.data.frame(VSS_euk@tax_table)$ASV_ID)]) %>%
  dplyr::filter(!is.na(trophic_mode)) %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% summarise_all(n_distinct)

trophic_mode_1 <- nrow(VSS_trophic_cluster.df %>% filter(trophic_mode ==1))
trophic_mode_2 <- nrow(VSS_trophic_cluster.df %>% filter(trophic_mode ==2))
trophic_mode_3 <- nrow(VSS_trophic_cluster.df %>% filter(trophic_mode ==3))
trophic_mode_4 <- nrow(VSS_trophic_cluster.df %>% filter(trophic_mode ==4))

mean_trophic_richness <- mean(VSS_trophic_cluster.df$trophic_mode)

VSS_trophic_cluster.df %>%
  dplyr::select(c("cluster_ID","trophic_mode")) %>%
  group_by(cluster_ID) %>% 
  summarise(total_count=n(),
            .groups = 'drop')

VSS_trophic_edge.df <-
  as.data.frame(VSS_pears$edgelist1 %>%
                  subset(v1 %in% VSS_trophic_network.df[,1]
                         & v2 %in% VSS_trophic_network.df[,1])) %>%
  dplyr::mutate(.,
                trophic_mode_v1 =
                  as.data.frame(VSS_euk@tax_table)$trophic_mode[match(.$v1,
                                                                      as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                trophic_mode_v2 =
                  as.data.frame(VSS_euk@tax_table)$trophic_mode[match(.$v2,
                                                                      as.data.frame(VSS_euk@tax_table)$ASV_ID)],
                edge_same_trophic_mode = ifelse(trophic_mode_v1 == trophic_mode_v2,"YES","NO")) %>%
  dplyr::filter(!is.na(trophic_mode_v1)) %>% dplyr::filter(!is.na(trophic_mode_v2))

nb_different_trophic_edge <-
  VSS_trophic_edge.df %>%
  filter(edge_same_trophic_mode=="NO") %>%
  nrow()

VSS_net_prop.df <- data.frame(lakeID= "VSS",
                              nb_ASV=n_ASV,
                              nb_components = summary(VSS_props)[[1]][1,1],
                              nb_clusters = length(unique(VSS_trophic_cluster.df$cluster_ID)),
                              nb_connected_node = nrow(VSS_trophic_network.df),
                              nb_edge = (nrow(VSS_pears$edgelist1 %>%
                                                subset(asso >0 & v1 %in% VSS_trophic_network.df[,1]
                                                       & v2 %in% VSS_trophic_network.df[,1]))
                                         + nrow(VSS_pears$edgelist1 %>%
                                                  subset(asso <0 & v1 %in% VSS_trophic_network.df[,1]
                                                         & v2 %in% VSS_trophic_network.df[,1]))),
                              nb_positive_edge = nrow(VSS_pears$edgelist1 %>%
                                                        subset(asso >0 & v1 %in% VSS_trophic_network.df[,1]
                                                               & v2 %in% VSS_trophic_network.df[,1])),
                              nb_negative_edge = nrow(VSS_pears$edgelist1 %>%
                                                        subset(asso <0 & v1 %in% VSS_trophic_network.df[,1]
                                                               & v2 %in% VSS_trophic_network.df[,1]))) %>%
  dplyr::mutate(positive_edge_percentage = (nb_positive_edge*100)/nb_edge,
                negative_edge_percentage = (nb_negative_edge*100)/nb_edge,
                largest_hub_degree_norm = summary(VSS_props)[[5]]$degree1[1,1],
                clustering_coefficient = summary(VSS_props)[[1]][2,1],
                modularity = summary(VSS_props)[[1]][3,1],
                edge_density = summary(VSS_props)[[1]][5,1],
                natural_connectivity = summary(VSS_props)[[1]][6,1],
                percentage_nodes_LCC = summary(VSS_props)[[2]][1,1],
                LCC_average_path_length = summary(VSS_props)[[2]][10,1],
                phototroph_connected_node_percent = as.numeric(percent_trophic_mode[1,2]),
                mixotrophs_connected_node_percent = as.numeric(percent_trophic_mode[2,2]),
                hetero_phagotrophs_connected_node_percent = as.numeric(percent_trophic_mode[3,2]),
                hetero_parasites_connected_node_percent = as.numeric(percent_trophic_mode[4,2]),
                percent_different_trophic_edge = nb_different_trophic_edge*100/nb_edge,
                percent_cluster_1_trophic_mode = trophic_mode_1*100/nb_clusters,
                percent_cluster_2_trophic_mode = trophic_mode_2*100/nb_clusters,
                percent_cluster_3_trophic_mode = trophic_mode_3*100/nb_clusters,
                percent_cluster_4_trophic_mode = trophic_mode_4*100/nb_clusters,
                mean_cluster_trophic_richness = mean_trophic_richness)
#### | ####

#### Properties ####

####__prop df ####
full_SparCC_network_prop.df <-
  rbind(VSM_net_prop.df,JAB_net_prop.df,CER_L_net_prop.df,
        CER_S_net_prop.df,CRE_net_prop.df,BLR_net_prop.df,
        LGP_net_prop.df,CSM_net_prop.df,VSS_net_prop.df) %>%
  as.data.frame()

write.csv(full_SparCC_network_prop.df,"full_SparCC_network_prop.csv",row.names = F)

####__Heatmap ####
heatmap_network_prop <- full_SparCC_network_prop.df %>%
  dplyr::select(c("lakeID",
           "nb_clusters","nb_connected_node",
           "nb_edge","positive_edge_percentage",
           "largest_hub_degree_norm", "clustering_coefficient", "modularity",
           "edge_density", "natural_connectivity",
           "percentage_nodes_LCC","LCC_average_path_length",
           "phototroph_connected_node_percent","mixotrophs_connected_node_percent",
           "hetero_phagotrophs_connected_node_percent","hetero_parasites_connected_node_percent",
           "percent_cluster_1_trophic_mode","percent_cluster_2_trophic_mode",
           "percent_cluster_3_trophic_mode","percent_cluster_4_trophic_mode")) %>%
  remove_rownames() %>%
  column_to_rownames("lakeID") %>%
  mutate_at(1:ncol(.), as.numeric) %>%
  dplyr::mutate(hetero_parasites_connected_node_percent = ifelse(is.na(hetero_parasites_connected_node_percent),
                                                                 0,hetero_parasites_connected_node_percent))

heatmap_network_prop[1:ncol(heatmap_network_prop)] <-
  apply(heatmap_network_prop[1:ncol(heatmap_network_prop)], 2, standardize)

heatmap_network_prop <- heatmap_network_prop %>%
  rownames_to_column(.,var = "lakeID") %>%
  pivot_longer(!lakeID, names_to = "property", values_to = "value") %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels= rev(c("VSM", "JAB", "CER_L",
                                              "CER_S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS"))),
                property=factor(property,levels=c("nb_node","nb_connected_node",
                                                  "clustering_coefficient","modularity",
                                                  "positive_edge_percentage","nb_edge",
                                                  "largest_hub_degree_norm","edge_density","natural_connectivity",
                                                  "percentage_nodes_LCC","LCC_average_path_length",
                                                  "phototroph_connected_node_percent","mixotrophs_connected_node_percent",
                                                  "hetero_phagotrophs_connected_node_percent","hetero_parasites_connected_node_percent",
                                                  "percent_cluster_1_trophic_mode","percent_cluster_2_trophic_mode",
                                                  "percent_cluster_3_trophic_mode","percent_cluster_4_trophic_mode"))) %>%
  ggplot(.,aes(property,lakeID,))+
  geom_tile(aes(fill = value), color = "black") +
  theme_bw()+
  theme(aspect.ratio = 0.25,
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size=16, color="black",
                                   angle = 45, hjust =1),
        axis.text.y = element_text(size=16, color="black",
                                   face= "bold",),
        axis.ticks = element_line(color="black"),
        plot.title = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.ticks.length = unit(12.5,"pt"),
        legend.ticks = element_line(size=0.5, color = "#f7f7f7"),
        legend.key.height = unit(40,"pt"),
        legend.key.width = unit(25,"pt"),
        legend.text = element_text(size=12))+
  scale_y_discrete(labels = rev(c("VSM", "JAB", "CER-L",
                              "CER-S", "CRE", "BLR","LGP", "CSM",
                              "VSS")))+
  scale_x_discrete(labels = c("Nb. components","Nb. con. components",
                              "Clustering coef.","Modularity",
                              "Positive edges (%)",
                              "Nb. edges",
                              "Max. node degrees",
                              "Edge density", "Nat. connectivity",
                              "LCC nodes (%)","LCC avg. path length",
                              "Phototrophs con. nodes (%)","Mixotrophs con. nodes (%)",
                              "Phagotrophs con. nodes (%)","Parasites con. nodes (%)",
                              "Modules with 1 trophic mode (%)","Modules with 2 trophic modes (%)",
                              "Modules with 3 trophic modes (%)","Modules with 4 trophic modes (%)"))+
  scale_fill_gradient2(high = "#FF0000",
                       mid = "#f7f7f7",
                       low = "#0000FF"
                       )+
  coord_cartesian(expand = FALSE)

heatmap_network_prop

ggsave("/Users/piefouca/Desktop/µEuk/Figures/network_heatmap.png",
       units = "in",dpi = 300,width = 13.4,height = 9.8)

####__vs. Chla ####
full_SparCC_network_Chla <-
  full_SparCC_network_prop.df %>%
  dplyr::mutate(., Chla_range=chla_range.df$range[match(.$lakeID,chla_range.df$lakeID)])

cor_test(vars = as.numeric(positive_edge_percentage), vars2 = Chla_range,
         method = "spearman", data=full_SparCC_network_Chla)

full_SparCC_network_Chla %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS"))) %>%
  ggplot(.,
       aes(Chla_range,positive_edge_percentage,
           fill = lakeID, color = lakeID, group =1))+
  geom_point(shape=21, size= 5, color = "black")+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=14,face = "bold"),
        axis.text = element_text(size=12, color="black"),
        axis.ticks = element_line(color="black"),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = rev(palette_lake_chla))+
  guides(fill=guide_legend(override.aes=list(shape=22,size=8),nrow=1),
         color = FALSE)+
  ggplot2::annotate("text",x = 300, y= 55,hjust=0,size=5,
                    label= expression(paste(italic("p")," <0.05 \u03C1 0.77")))+
  coord_cartesian(xlim = c(0,400), ylim = c(50,100), expand = F)+
  labs(x=expression(paste(bold(" Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))),
       y="Network positive edge (%)", fill= "Lake")

ggsave("/Users/piefouca/Desktop/µEuk/Figures/positive_edge_chla_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__PCA ####
network.pca <- full_SparCC_network_prop.df %>%
  dplyr::select(c("lakeID",
                  "nb_clusters","nb_connected_node",
                  "nb_edge","positive_edge_percentage",
                  "largest_hub_degree_norm", "clustering_coefficient", "modularity",
                  "edge_density", "natural_connectivity",
                  "percentage_nodes_LCC","LCC_average_path_length",
                  "phototroph_connected_node_percent","mixotrophs_connected_node_percent",
                  "hetero_phagotrophs_connected_node_percent","hetero_parasites_connected_node_percent",
                  "percent_cluster_1_trophic_mode","percent_cluster_2_trophic_mode",
                  "percent_cluster_3_trophic_mode")) %>%
  rename("Nb. components"="nb_clusters",
         "Nb. con. components"="nb_connected_node",
         "Nb. edges"="nb_edge",
         "Positive edges (%)"="positive_edge_percentage",
         "Max. node degrees"="largest_hub_degree_norm",
         "Clustering coef."="clustering_coefficient",
         "Modularity"="modularity",
         "Edge density"="edge_density",
         "Nat. connectivity"="natural_connectivity",
         "LCC nodes (%)"="percentage_nodes_LCC",
         "LCC avg. path length"="LCC_average_path_length",
         "Phototrophs con. nodes (%)"="phototroph_connected_node_percent",
         "Mixotrophs con. nodes (%)"="mixotrophs_connected_node_percent",
         "Phagotrophs con. nodes (%)"="hetero_phagotrophs_connected_node_percent",
         "Parasites con. nodes (%)"="hetero_parasites_connected_node_percent",
         "Modules with 1 trophic mode (%)"="percent_cluster_1_trophic_mode",
         "Modules with 2 trophic modes (%)"="percent_cluster_2_trophic_mode",
         "Modules with 3 trophic modes (%)"="percent_cluster_3_trophic_mode") %>%
  dplyr::mutate(lakeID= ifelse(lakeID == "CER_L","CER-L",
                               ifelse(lakeID == "CER_S","CER-S",lakeID)),
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS"))) %>%
  remove_rownames() %>%
  column_to_rownames("lakeID") %>%
  mutate_at(1:ncol(.), as.numeric) %>%
  dplyr::mutate(`Parasites con. nodes (%)` = ifelse(is.na(`Parasites con. nodes (%)`),
                                                                 0,`Parasites con. nodes (%)`))

pca_var_trophic<-network.pca

pca_var_trophic[1:18] <-
  apply(pca_var_trophic[1:18], 2, standardize)

pca_var_trophic <- get_pca_var(prcomp(pca_var_trophic))

var_dim_trophic <- pca_var_trophic$cos2 %>% .[,c(1,2,3,4)]

network.pca[1:18] <-
  apply(network.pca[1:18], 2, standardize)

network.dist <-network.pca %>% vegdist(., method="euclidean")

PCA_network.mds<- cmdscale(network.dist,eig=TRUE, k=4)

PCA_network.df <- PCA_network.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,.before=1,
                lakeID=factor(rownames(.),levels=c("VSM", "JAB", "CER-L",
                                                   "CER-S", "CRE", "BLR","LGP", "CSM", 
                                                   "VSS")))
network_prop <- envfit(PCA_network.mds$points ~
                       #   `Nb. components`+`Nb. con. components`+
                       # `Clustering coef.`+`Modularity`+
                       # `Positive edges (%)`+
                       # `Nb. edges`+
                       # `Max. node degrees`+
                       # `Edge density`+ `Nat. connectivity`+
                       # `LCC nodes (%)`+`LCC avg. path length`+
                       # `Phototrophs con. nodes (%)`+`Mixotrophs con. nodes (%)`+
                       # `Phagotrophs con. nodes (%)`+`Parasites con. nodes (%)`+
                       # `Cluster with 1 trophic mode (%)`+`Cluster with 2 trophic modes (%)`+
                       # `Cluster with 3 trophic modes (%)`,
                         `Nb. components`+#`Nb. con. components`+
                         `Clustering coef.`+`Modularity`+
                         `Positive edges (%)`+
                         `Nb. edges`+
                         `Max. node degrees`+
                         `Edge density`+ `Nat. connectivity`+
                         `LCC nodes (%)`+`LCC avg. path length`+
                        # `Phototrophs con. nodes (%)`+`Mixotrophs con. nodes (%)`+
                         `Phagotrophs con. nodes (%)`+#`Parasites con. nodes (%)`+
                         #`Modules with 1 trophic mode (%)`+`Modules with 2 trophic modes (%)`+
                         `Modules with 3 trophic modes (%)`,
                       data = as.data.frame(network.pca), perm = 999)

network_vectors.df <- as.data.frame(scores(network_prop, display = "vectors")) 
network_vectors.df <- cbind(network_vectors.df, prop = rownames(network_vectors.df)) 

PCA_network <- PCA_network.df %>%
  ggplot(.,aes(dim1,dim2,color=lakeID,label=lakeID)) +
  geom_segment(data = network_vectors.df,
               mapping=aes(x = 0, xend = Dim1*2.5, y = 0, yend = Dim2*2.5),
               color = "gray40",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  geom_text(data = network_vectors.df,
            mapping=aes(x = Dim1*2.62, y = Dim2*2.62, label = prop),
            size = 5,color="gray40",fontface = "bold",inherit.aes = F)+
  geom_text(size=7, check_overlap = T,show.legend = F,fontface="bold")+
  #geom_point(size=2,show.legend = T)+
  theme_bw()+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=16,face = "bold"),
        axis.text = element_text(size=15, color="black"),
        axis.ticks = element_line(color="black"))+
  scale_color_manual(values=rev(palette_lake_chla))+
  labs(x=paste0("PC1 [",round(PCA_network.mds$eig[1]*100/sum(PCA_network.mds$eig),1),"%]"),
       y=paste0("PC2 [",round(PCA_network.mds$eig[2]*100/sum(PCA_network.mds$eig),1),"%]"))
PCA_network

ggsave("/Users/piefouca/Desktop/µEuk/Figures/PCA_network_prop.svg",
       units = "in",dpi = 300,width = 13.4,height = 9.8)

network_prop_corr <- network.pca %>%
  dplyr::mutate(.,chla_mean = chla_mean.df$chla_mean[match(rownames(.),
                                                           chla_mean.df$lakeID)])
var_tested=network_prop_corr[,19]
var_to_test=network_prop_corr[,1:18]

list_corr <- c("network_prop","pvalue","r2")

for (i in 1:length(var_to_test) ){
   tmp <- cor.test(var_tested, var_to_test[[i]], method="spearman")
   corr_results <- c(colnames(var_to_test)[i],tmp$p.value,tmp$estimate[[1]])
   list_corr<-rbind(list_corr,corr_results)
}
View(list_corr)


####________________________####