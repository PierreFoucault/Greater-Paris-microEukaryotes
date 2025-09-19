#### Import data ####
 
####__18S-PR2 (rarefied at 12,287 reads) ####
rar_euk <-
  create_phyloseq_PR2("18S_data/PR2_ASV_table.tsv",
                  "18S_data/PR2_ASV-tax.tsv",
                  "18S_data/metadata_18S.txt")


View(rar_euk@tax_table)
#VSM_B_W1 removed during the preprocess (rarefaction)

####__trophic mode ####
tmp<-
  read_delim("18S_data/pr2_5.1.0_ecological_function.csv",
             delim = ",",show_col_types = F) %>%
  filter(class %in%
           (as.data.frame(rar_euk@tax_table)$Class)) %>%
  #select(c("domain","subdivision","ecological_function")) %>%
  select(c("domain","subdivision","class","order","family","ecological_function")) %>%
  filter(domain == "Eukaryota") %>%
  distinct(.keep_all = TRUE)
Vew(tmp)
write.csv(tmp,"18S_data/trophic_mode_PR2_Class.csv",row.names = F)

####____Subdivision ####
# list_phototroph_subdivision <- c("Gyrista","Chlorophyta_X", "Haptophyta_X",
#                                  "Glaucophyta_X","Rhodophyta_X")
# 
# list_mixotroph_subdivision <- c("Cryptophyta_X", "Dinoflagellata")
# 
# list_heterotroph_subdivision <- c("Ciliophora","Kathablepharida","Fungi",
#                                 "Cercozoa","Bigyra","Choanoflagellata",
#                                 "Chrompodellids","Evosea_X","Tubulinea_X",
#                                 "Centroplasthelida_X","Filasterea","Rotosphaerida",
#                                 "Rigifilida_X","Hemimastigophora_X","Rhizaria_X",
#                                 "Pluriformea","Discosea_X","Colponemidia",
#                                 "Telonemia_X")
#    
# list_parasite_subdivision <- c("Perkinsea","Discoba_X",
#                                "Parabasalia",'Preaxostyla',
#                                "Apicomplexa","Alveolata_X")

#unknown trophic mode -> "Nebulidia_X" 'Stramenopiles_X' #"Ancyromonadida_X",
#"Cercozoa Novel-clade-10-12"
#651 ASVs not annotated at the Subdivision rank

# rar_euk@tax_table<-
#   tax_table(data.frame(tax_table(rar_euk)) %>%
#               dplyr::mutate(.,
#                             trophic_mode = ifelse(Subdivision %in% list_phototroph_subdivision,"Phototroph",
#                                                   ifelse(Subdivision %in% list_heterotroph_subdivision,"Heterotroph",
#                                                          ifelse(Subdivision %in% list_mixotroph_subdivision,"Mixotroph",
#                                                                 ifelse(Subdivision %in% list_parasite_subdivision,"Parasite",
#                                                                        "unknown")))),
#                             trophic_subdivision = ifelse(is.na(Subdivision),"NA",gsub("_X.","",Subdivision)),
#                             ASV_ID=paste("ASV", 1:ntaxa(rar_euk), sep = "_")) %>%
#               as.matrix(.))

####____Class ####
class_trophic_mode.df <- read_delim("18S_data/trophic_mode_PR2_Class.csv",";",
                                    show_col_types = F)
rar_euk <-
  create_phyloseq_PR2("18S_data/PR2_ASV_table.tsv",
                      "18S_data/PR2_ASV-tax.tsv",
                      "18S_data/metadata_18S.txt")

rar_euk@tax_table<-
  tax_table(data.frame(tax_table(rar_euk)) %>%
              dplyr::mutate(.,
                            trophic_mode =
                              class_trophic_mode.df$trophic_mode[match(.$Class,
                                                                       class_trophic_mode.df$class)],
                            trophic_subdivision = Subdivision,
                            trophic_subdivision = gsub("_X","",trophic_subdivision),
                            trophic_subdivision = ifelse(is.na(Subdivision),"not_annotated",trophic_subdivision),
            
                            trophic_subdivision = ifelse(trophic_subdivision == "Gyrista" &
                                                           Class %in% c(
                                                             "Bacillariophyceae",
                                                             "Mediophyceae",
                                                             "Coscinodiscophyceae"),
                                                         "Gyrista · diatoms",trophic_subdivision),
            
                            trophic_subdivision = ifelse(trophic_subdivision == "Gyrista" &
                                                           Class %out% c(
                                                             "Bacillariophyceae",
                                                             "Mediophyceae",
                                                             "Coscinodiscophyceae"),
                                                         "Gyrista · non-diatoms",trophic_subdivision),
                            trophic_subdivision = ifelse(trophic_subdivision == "Gyrista" &
                                                           trophic_mode == "heterotrophic parasites",
                                                         "Parasitic Gyrista",trophic_subdivision),
                            trophic_subdivision = ifelse(trophic_subdivision == "Fungi" &
                                                           trophic_mode == "heterotrophic parasites",
                                                         "Parasitic Fungi",trophic_subdivision),
                            trophic_subdivision = ifelse(Subdivision == "Fungi" & is.na(trophic_mode),
                                                         "Fungi",trophic_subdivision),
                            trophic_subdivision = ifelse(is.na(Subdivision) & is.na(trophic_mode),
                                                         "not_annotated",trophic_subdivision),
                            ASV_ID=paste("ASV", 1:ntaxa(rar_euk), sep = "_")) %>%
              as.matrix(.))


#write_tsv(as.data.frame(rar_euk@tax_table) %>% rownames_to_column(var = "QIIME_ID"),
          #"18S_data/18S_mode-tax.tsv")
View(rar_euk@tax_table)
#### | ####

#### Alpha-Diversity ####

####__Rar. Curve ####
Fig_rar_curve_18S <-
  rarecurve(otu_table(rar_euk) %>% as.data.frame() %>% t(), step=500, cex=0.2,tidy = T) %>%
  ggplot(.,aes(x=Sample,y=Species,group=Site),color="black")+
  geom_line(linewidth = 0.2)+theme_bw()+
  labs( y = "Eukaryota ASVs Richness", x = "Nb. of reads")+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        axis.ticks = element_blank())+
  scale_x_continuous(expand=c(0,0),
                     limits = c(0,12288),
                     breaks = c(0,2000,4000,6000,8000,10000,12000))+
  scale_y_continuous(expand=c(0,0),
                     limits = c(0,600),
                     breaks = c(0,200,400,600))
Fig_rar_curve_18S

####__Indices ####
alphadiv_18S <- rar_euk@otu_table %>% estimate_richness() %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_euk@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_euk@sam_data$lakeID[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                season_year=
                  (rar_euk@sam_data$season_year[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))
#View(alphadiv_18S)
write.csv(alphadiv_18S,"output/alphadiv_18S.csv")

View(alphadiv_18S %>%
  group_by(season_year) %>%
  summarise(Richness_mean=mean(Richness),
            Richness_sd=sd(Richness)))

####______Richness ####
richness_18S<- alphadiv_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1, show.legend = F)+
  geom_point(shape=21, color =  NA, size = 0, stroke = 0.2, show.legend = T, alpha = 0.6)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 600),
                     breaks = c(0, 200, 400, 600))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                     #labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  guides(fill=guide_legend(override.aes=list(shape=22,size=8, stroke =0.5, alpha = 1,color = "black"), nrow=1),
         color = FALSE)+
  labs(#y = 'ASV richness',
    title = 'ASV richness',
    
       fill = "Lake")
richness_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_18S.png",units = "in",dpi = 300,width = 13.4,height = 9.8)
 
####______Evenness ####
evenness_18S<- alphadiv_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Evenness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1)+
  #geom_point(shape=21, color = "black", size = 1, stroke = 0.2, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black"),
        axis.text = element_text(size=12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Eukaryota ASV Evenness')
evenness_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/evenness_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####______Shannon ####
shannon_18S<- alphadiv_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Shannon, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
       # axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0,5.5),
                     breaks = c(0, 1, 3, 5))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                     #labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(#y = 'ASV Shannon diversity',
       title = 'ASV Shannon diversity')
shannon_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/shannon_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####______ Rare 1% ####
asv_rare <- get_asv_above_threshold(rar_euk, threshold = 1)
length(asv_rare)
(ntaxa(rar_euk) - length(asv_rare)) *100 / ntaxa(rar_euk)

####______ Rare 0.1% ####
asv_rare <- get_asv_above_threshold(rar_euk, threshold = 0.1)
length(asv_rare)
(ntaxa(rar_euk) - length(asv_rare)) *100 / ntaxa(rar_euk)

####__Phototrophs ####
alphadiv_phototroph_18S <- subset_taxa(rar_euk, trophic_mode == "phototrophs")
ntaxa(alphadiv_phototroph_18S)*100/ntaxa(rar_euk)
sum(rowSums(alphadiv_phototroph_18S@otu_table))*100/sum(rowSums(rar_euk@otu_table))

alphadiv_phototroph_18S <- alphadiv_phototroph_18S@otu_table %>% estimate_richness() %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_euk@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_euk@sam_data$lakeID[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))

####______Richness ####
richness_phototroph_18S<- alphadiv_phototroph_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1, show.legend = F)+
  #geom_point(shape=21, color = "black", size = 1, stroke = 0.2, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 200),
                     breaks = c(0, 100, 200))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                     #labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(#y = 'Phototrophs ASV richness',
       title = 'Phototrophs ASV richness')
richness_phototroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_phototroph_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Het. phagotrophs ####
alphadiv_heterotroph_18S <- subset_taxa(rar_euk, trophic_mode == "heterotrophic phagotrophs")
ntaxa(alphadiv_heterotroph_18S)*100/ntaxa(rar_euk)
sum(rowSums(alphadiv_heterotroph_18S@otu_table))*100/sum(rowSums(rar_euk@otu_table))

alphadiv_heterotroph_18S <- alphadiv_heterotroph_18S@otu_table %>% estimate_richness() %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_euk@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_euk@sam_data$lakeID[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))

####______Richness ####
richness_heterotroph_18S<- alphadiv_heterotroph_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1, show.legend = F)+
  #geom_point(shape=21, color = "black", size = 1, stroke = 0.2, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 260),
                     breaks = c(0, 100, 200))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                     #labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(#y = 'Heterotrophic phagotrophs ASV richness',
       title = 'Heterotrophic phagotrophs ASV richness')
richness_heterotroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_heterotroph_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Mixotrophs ####
alphadiv_mixotroph_18S <- subset_taxa(rar_euk, trophic_mode == "mixotrophs")
ntaxa(alphadiv_mixotroph_18S)*100/ntaxa(rar_euk)
sum(rowSums(alphadiv_mixotroph_18S@otu_table))*100/sum(rowSums(rar_euk@otu_table))


alphadiv_mixotroph_18S <- alphadiv_mixotroph_18S@otu_table %>% estimate_richness() %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_euk@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_euk@sam_data$lakeID[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))

####______Richness ####
richness_mixotroph_18S<- alphadiv_mixotroph_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1,show.legend = F)+
  #geom_point(shape=21, color = "black", size = 1, stroke = 0.2, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100),
                     breaks = c(0, 50, 100))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                    # labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(#y = 'Mixotrophs ASV richness',
       title = 'Mixotrophs ASV richness')
richness_mixotroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_mixotroph_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Het. parasites ####
alphadiv_parasite_18S <- subset_taxa(rar_euk, trophic_mode == "heterotrophic parasites")
ntaxa(alphadiv_parasite_18S)*100/ntaxa(rar_euk)
sum(rowSums(alphadiv_parasite_18S@otu_table))*100/sum(rowSums(rar_euk@otu_table))

alphadiv_parasite_18S <- alphadiv_parasite_18S@otu_table %>% estimate_richness(measures = c("Observed","Shannon")) %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_euk@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_euk@sam_data$lakeID[match(.$sample_code,rar_euk@sam_data$sample_code)]),
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))

####______Richness ####
richness_parasite_18S<- alphadiv_parasite_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
  geom_smooth(method = "loess", alpha=.1, linewidth = 1,show.legend = F)+
  #geom_point(shape=21, color = "black", size = 1, stroke = 0.2, show.legend = F)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2, linetype = "dashed"),
        plot.title = element_text(size=12, color = "black", hjust = 0, face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_line(color="black", linewidth = 0.4),
        axis.ticks.length=unit(4, "pt"),
        axis.text = element_text(size=12, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 120),
                     breaks = c(0, 30, 60, 90, 120))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 90, 180, 360, 540))+
                     #labels = c("0", "30", " 60", "120", "180", "360", "540"))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(#y = 'Heterotrophic parasites ASV richness',
       title = 'Heterotrophic parasites ASV richness')
richness_parasite_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_parasite_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### | ####

#### Composition ####

####__season-Sub. ####

####____Subdivision ####
barplot_subdivision_season.df <-
  rar_euk %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_season,Subdivision) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_season) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_euk@sam_data$lakeID[match(.$lake_season,rar_euk@sam_data$lake_season)]),
        season_year=
          (rar_euk@sam_data$season_year[match(.$lake_season,rar_euk@sam_data$lake_season)]),
        .)

barplot_subdivision_season <- barplot_subdivision_season.df %>%
  dplyr::mutate(.,
                taxa_legend=ifelse(median_abundance<1,"Taxa < 1%",Subdivision)) %>%
  dplyr::mutate(taxa_legend = recode(taxa_legend, 'Cryptophyta_X' = 'Cryptophyta', 
                                     'Telonemia_X' = 'Telonemia', 
                                     'Chlorophyta_X' = 'Chlorophyta',
                                     'Haptophyta_X'='Haptophyta'),
                lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS")),
                taxa_legend=factor(taxa_legend,
                                   levels=c("Gyrista","Chlorophyta",'Haptophyta',
                                            "Cryptophyta","Dinoflagellata","Ciliophora",
                                            "Kathablepharida","Fungi","Cercozoa",
                                            "Bigyra","Choanoflagellata", "Chrompodellids",
                                            "Telonemia","Perkinsea",'Ichthyosporea',
                                            'Taxa < 1%')),
                season_year = factor(season_year, levels = c("summer_2021","fall_2021","winter_2021",
                                                             "spring_2022","summer_2022","fall_2022","winter_2022"))) %>%
  ggplot(.,aes(x=season_year,y=median_abundance,fill=taxa_legend))+
  theme_bw()+
  geom_col(width = .98,color = NA, linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank())+
  scale_x_discrete(expand = c(0,0))+ #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0),
                     labels = c('0%','25%','50%','75%','100%'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Subdivision)+ #set a color palette for the barplot
  labs(y = "", fill="Eukaryota\nSubdivision", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'white',
                #                               color = NA),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black") +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -5, ymax = -0.3, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black")+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,7.5), ylim = c(-5,NA),clip = "off")
barplot_subdivision_season

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_subdivision_season_barplot.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####____Class ####

barplot_class_season.df <-
  rar_euk %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_season,Class) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_season) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_euk@sam_data$lakeID[match(.$lake_season,rar_euk@sam_data$lake_season)]),
        season_year=
          (rar_euk@sam_data$season_year[match(.$lake_season,rar_euk@sam_data$lake_season)]),
        .)

top20_class <- barplot_class_season.df %>% group_by(Class) %>%
  summarise(sum_abund = sum(total_count)) %>%
  dplyr::mutate(.,
                Class=gsub("_XX","",Class),
                Class=gsub("_X","",Class)) %>%
  arrange(.,desc(sum_abund),Class) %>% .[1:25,]

#View(top20_class)
#unique(barplot_class_season$data$taxa_legend)

barplot_class_season <- barplot_class_season.df %>%
  dplyr::mutate(.,
                #taxa_legend=ifelse(median_abundance<5,"Taxa < 5%",Class)) %>%
                subdivision=
                  (as.data.frame(rar_euk@tax_table)$Subdivision[match(.$Class,as.data.frame(rar_euk@tax_table)$Class)]),
                subdivision=gsub("_X","",subdivision),
                Class=gsub("_XX","",Class),
                Class=gsub("_X","",Class),
                taxa_legend=ifelse(Class %in% top20_class$Class,paste0(subdivision," · ",Class),"Others")) %>%
  dplyr::mutate(
                lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM",
                                               "VSS")),
                taxa_legend=factor(taxa_legend,
                   levels=c("Gyrista · Chrysophyceae",
                            "Gyrista · Bacillariophyceae",
                            "Gyrista · Mediophyceae",
                            "Chlorophyta · Chlorophyceae",
                            "Chlorophyta · Chlorodendrophyceae",
                            "Chlorophyta · Trebouxiophyceae",
                            "Cryptophyta · Cryptophyceae",
                            "Dinoflagellata · Dinophyceae",
                            "Ciliophora · CONThreeP",
                            "Ciliophora · Oligohymenophorea",
                            "Ciliophora · Litostomatea",
                            "Ciliophora · Nassophorea",
                            "Ciliophora · Spirotrichea",
                            "Kathablepharida · Kathablepharidea",
                            "Fungi · Fungi",
                            "Fungi · NA",
                            "Cercozoa · Filosa-Imbricatea",
                            "Bigyra · Bicoecea",
                            "Bigyra · Opalozoa",
                            "Choanoflagellata · Choanoflagellatea",
                            "Chrompodellids · Colpodellidea",
                            "Telonemia · Telonemia",
                            "Fungi · Chytridiomycota",
                            "Perkinsea · Perkinsida",
                            "Perkinsea · Perkinsea",
                            "Others")),
                season_year = factor(season_year, levels = c("summer_2021","fall_2021","winter_2021",
                                                             "spring_2022","summer_2022","fall_2022","winter_2022"))) %>%
  ggplot(.,aes(x=season_year,y=median_abundance,
               fill= taxa_legend))+
  theme_bw()+
  geom_col(width = .98,color = NA,
           linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        #legend.position = "bottom",
       # legend.direction = "vertical",
        legend.title = element_text(face="bold", size= 14),
        legend.text = element_text(size=12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank())+
  scale_x_discrete(expand = c(0,0))+ #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0),
                     labels = c('0%','25%','50%','75%','100%'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Class)+ #set a color palette for the barplot
  #guides(fill=guide_legend(nrow=7))+
  guides(fill=guide_legend(ncol=1, color = "black"))+
  labs(y = "", fill="Microeukaryotes", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'white',
                #                               color = NA),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black") +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -5, ymax = -0.3, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black")+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,7.5), ylim = c(-5,NA),clip = "off")
barplot_class_season

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_class_season_barplot.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)
ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_class_season_barplot.svg",units = "in",dpi = "retina",width = 13.4,height = 9.8)

write.xlsx((barplot_class_season.df %>%
              dplyr::mutate(.,
                            #taxa_legend=ifelse(median_abundance<5,"Taxa < 5%",Class)) %>%
                            subdivision=
                              (as.data.frame(rar_euk@tax_table)$Subdivision[match(.$Class,as.data.frame(rar_euk@tax_table)$Class)]),
                            subdivision=gsub("_X","",subdivision),
                            Class=gsub("_XX","",Class),
                            Class=gsub("_X","",Class),
                            taxa_legend=ifelse(Class %in% top20_class$Class,paste0(subdivision," · ",Class),"Others")) %>%
              dplyr::mutate(
                lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM",
                                               "VSS")),
                taxa_legend=factor(taxa_legend,
                                   levels=c("Gyrista · Chrysophyceae",
                                            "Gyrista · Bacillariophyceae",
                                            "Gyrista · Mediophyceae",
                                            "Chlorophyta · Chlorophyceae",
                                            "Chlorophyta · Chlorodendrophyceae",
                                            "Chlorophyta · Trebouxiophyceae",
                                            "Cryptophyta · Cryptophyceae",
                                            "Dinoflagellata · Dinophyceae",
                                            "Ciliophora · CONThreeP",
                                            "Ciliophora · Oligohymenophorea",
                                            "Ciliophora · Litostomatea",
                                            "Ciliophora · Nassophorea",
                                            "Ciliophora · Spirotrichea",
                                            "Kathablepharida · Kathablepharidea",
                                            "Fungi · Fungi",
                                            "Fungi · NA",
                                            "Cercozoa · Filosa-Imbricatea",
                                            "Bigyra · Bicoecea",
                                            "Bigyra · Opalozoa",
                                            "Choanoflagellata · Choanoflagellatea",
                                            "Chrompodellids · Colpodellidea",
                                            "Telonemia · Telonemia",
                                            "Fungi · Chytridiomycota",
                                            "Perkinsea · Perkinsida",
                                            "Perkinsea · Perkinsea",
                                            "Others")),
                season_year = factor(season_year, levels = c("summer_2021","fall_2021","winter_2021",
                                                             "spring_2022","summer_2022","fall_2022","winter_2022")),
                taxa_legend=gsub(" · ","-",taxa_legend)) %>%
              group_by(lakeID,season_year,taxa_legend) %>%
              summarise(median_abundance=sum(median_abundance))),
           "output/euk_class_season_rel_abund.xlsx")
  

####__month-Sub. ####

####____Subdivision ####
barplot_subdivision.df <- rar_euk %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_month,Subdivision) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_euk@sam_data$lakeID[match(.$lake_month,rar_euk@sam_data$lake_month)]),
        month_code=
          (rar_euk@sam_data$month_code[match(.$lake_month,rar_euk@sam_data$lake_month)]),
        .)

# View(barplot_subdivision.df %>%
#   group_by(Subdivision) %>%
#   summarise(min=min(median_abundance),
#             max=max(median_abundance),
#             mean=mean(median_abundance),
#             sd=sd(median_abundance)))
# View(rar_euk@tax_table)
# head(sort(rowSums(rar_euk@otu_table),decreasing = T))

barplot_subdivision <- barplot_subdivision.df %>%
  dplyr::mutate(taxa_legend=ifelse(median_abundance<1,"Taxa < 1%",Subdivision),
                taxa_legend = recode(taxa_legend, 'Cryptophyta_X' = 'Cryptophyta', 
                               'Telonemia_X' = 'Telonemia', 
                               'Chlorophyta_X' = 'Chlorophyta',
                               'Haptophyta_X'='Haptophyta'),
                lakeID=factor(lakeID, levels =c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS")),
                taxa_legend=factor(taxa_legend,
                                    levels=c("Gyrista","Chlorophyta",'Haptophyta',
                                             "Cryptophyta","Dinoflagellata","Ciliophora",
                                             "Kathablepharida","Fungi","Cercozoa",
                                             "Bigyra","Choanoflagellata", "Chrompodellids",
                                             "Telonemia","Perkinsea",'Ichthyosporea',
                                             'Taxa < 1%')),
                month_code = recode(month_code,
                                      'A'=' 6','B'=' 7','C'=' 8','D'=' 9',
                                      'E'=' 10','F'=' 11','G'='1','H'='2',
                                      'I'='3','J'='4','K'='5','L'='6','M'='7',
                                      'N'='8','O'='9','P'='10','Q'='11','R'='12'),
                month_code=factor(month_code, levels=c(' 6',' 7',' 8',' 9',' 10',' 11','1',
                                             '2', '3','4','5','6','7','8','9','10',
                                             '11', '12'))) %>%
  ggplot(.,aes(x=month_code,y=median_abundance,fill=taxa_legend))+
  theme_bw()+
  geom_col(width = .98, color = NA, linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank())+
  scale_x_discrete(expand = c(0,0),
                   labels = c('J','J','A','S','O','N',
                              'J','F', 'M','A','M','J',
                              'J','A','S','O','N', 'D'))+ #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0),
                     labels = c('0%','25%','50%','75%','100%'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Subdivision)+ #set a color palette for the barplot
  labs(y = "", fill="Eukaryota\nSubdivision", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'lightgrey',
                #                               color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))+ #x axis header on the top (default=bottom)
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -5, ymax = -0.3, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(-5,NA),clip = "off")
barplot_subdivision

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_subdivision_month_barplot.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####____Class ####
barplot_class_month.df <-
  rar_euk %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_month,Class) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_euk@sam_data$lakeID[match(.$lake_month,rar_euk@sam_data$lake_month)]),
        month_code=
          (rar_euk@sam_data$month_code[match(.$lake_month,rar_euk@sam_data$lake_month)]),
        .)

View(barplot_class_month.df %>%
       group_by(Class) %>%
       summarise(min=min(median_abundance),
                 max=max(median_abundance),
                 mean=mean(median_abundance),
                 sd=sd(median_abundance)))

View(barplot_class_month$data)
barplot_class_month <- barplot_class_month.df %>%
  dplyr::mutate(.,
                subdivision=
                  (as.data.frame(rar_euk@tax_table)$Subdivision[match(.$Class,as.data.frame(rar_euk@tax_table)$Class)]),
                subdivision=gsub("_X","",subdivision),
                Class=gsub("_XX","",Class),
                Class=gsub("_X","",Class),
                taxa_legend=ifelse(Class %in% top20_class$Class,paste0(subdivision," · ",Class),"Others")) %>%
  dplyr::mutate(
    lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                   "CER-S", "CRE", "BLR","LGP", "CSM",
                                   "VSS")),
    taxa_legend=factor(taxa_legend,
                       levels=c("Gyrista · Chrysophyceae",
                                "Gyrista · Bacillariophyceae",
                                "Gyrista · Mediophyceae",
                                "Chlorophyta · Chlorophyceae",
                                "Chlorophyta · Chlorodendrophyceae",
                                "Chlorophyta · Trebouxiophyceae",
                                "Cryptophyta · Cryptophyceae",
                                "Dinoflagellata · Dinophyceae",
                                "Ciliophora · CONThreeP",
                                "Ciliophora · Oligohymenophorea",
                                "Ciliophora · Litostomatea",
                                "Ciliophora · Nassophorea",
                                "Ciliophora · Spirotrichea",
                                "Kathablepharida · Kathablepharidea",
                                "Fungi · Fungi",
                                "Fungi · NA",
                                "Cercozoa · Filosa-Imbricatea",
                                "Bigyra · Bicoecea",
                                "Bigyra · Opalozoa",
                                "Choanoflagellata · Choanoflagellatea",
                                "Chrompodellids · Colpodellidea",
                                "Telonemia · Telonemia",
                                "Fungi · Chytridiomycota",
                                "Perkinsea · Perkinsida",
                                "Perkinsea · Perkinsea",
                                "Others")),
    month_code = recode(month_code,
                        'A'=' 6','B'=' 7','C'=' 8','D'=' 9',
                        'E'=' 10','F'=' 11','G'='1','H'='2',
                        'I'='3','J'='4','K'='5','L'='6','M'='7',
                        'N'='8','O'='9','P'='10','Q'='11','R'='12'),
    month_code=factor(month_code, levels=c(' 6',' 7',' 8',' 9',' 10',' 11','1',
                                           '2', '3','4','5','6','7','8','9','10',
                                           '11', '12'))) %>%
  ggplot(.,aes(x=month_code,y=median_abundance,
               fill= taxa_legend))+
  theme_bw()+
  geom_col(width = .98, color = NA, linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank())+
  scale_x_discrete(expand = c(0,0),
                   labels = c('J','J','A','S','O','N',
                              'J','F', 'M','A','M','J',
                              'J','A','S','O','N', 'D'))+ #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0),
                     labels = c('0%','25%','50%','75%','100%'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Class)+ #set a color palette for the barplot
  #guides(fill=guide_legend(nrow=7))+
  guides(fill=guide_legend(ncol=1, color = "black"))+
  labs(y = "", fill="Microeukaryotic ", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'white',
                #                               color = NA),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -5, ymax = -0.3, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -5, ymax = -0.3, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = -5, ymax = -0.3, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = -5, ymax = -0.3, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(-5,NA),clip = "off")
barplot_class_month

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_class_month_barplot.svg",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### | ####

#### Beta-Diversity ####

BC_18S.dist <- vegdist(t(rar_euk@otu_table),method = "bray")
JC_18S.dist <- vegdist(t(rar_euk@otu_table),
                                method = "jaccard",binary = T)
mantel_stat_18S <-
  mantel(BC_18S.dist,JC_18S.dist,
         method = "spearman",permutations = 999)$statistic

####__Bray-Curtis ####
BC_18S_out <- get_MDS_output(rar_euk@otu_table,
                             BC_18S.dist,
                             2)

BC_18S_out[[1]] <- BC_18S_out[[1]] %>%
  cbind(.,
        lakeID = rar_euk@sam_data$lakeID[match(rownames(.),
                                               rownames(rar_euk@sam_data))],
        month_code = rar_euk@sam_data$month_code[match(rownames(.),
                                                       rownames(rar_euk@sam_data))],
        lake_season = rar_euk@sam_data$lake_season[match(rownames(.),
                                                         rownames(rar_euk@sam_data))],
        season_year = rar_euk@sam_data$season_year[match(rownames(.),
                                                         rownames(rar_euk@sam_data))],
        year = rar_euk@sam_data$year[match(rownames(.),
                                                         rownames(rar_euk@sam_data))])

BC_18S <- BC_18S_out[[1]] %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                            "CER-S", "CRE", "BLR","LGP", "CSM", 
                                            "VSS")),
                season_year=factor(season_year,levels = c("summer_2021","fall_2021","winter_2021",
                                                "spring_2022","summer_2022","fall_2022","winter_2022")),
                year = ifelse(season_year == "winter_2021", "2021", year),
                year=factor(year,levels = c("2021", "2022"))) %>%
  ggplot(.,aes(Axis1,Axis2,color=season_year,fill=season_year, shape = year, group = season_year))+
  geom_convexhull(alpha=0.6,show.legend = F,size=0.6,)+
  geom_point(show.legend = T,color="black",size=3)+
  theme_bw() +
  theme(aspect.ratio = 0.7,
        panel.grid = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10, color="black"),
        axis.ticks = element_line(color="black"),
        plot.title = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values = c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"),
                     labels = c("Summer 2021","Fall 2021","Winter 2021",
                                "Spring 2022","Summer 2022","Fall 2022","Winter 2022"))+
  scale_fill_manual(values = c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"))+
  scale_shape_manual(values = c(21,22))+
  #scale_x_reverse()+
 # scale_y_reverse()+
  scale_y_continuous(expand = c(0,0), trans = "reverse",
                     limits = c(0.36,-0.38),
                     breaks = c(0.2, 0, -0.2),
                     labels = c(0.2, 0, -0.2))+
  scale_x_continuous(expand = c(0,0),
                     limits = c(-0.36,0.5),
                     breaks = c(-0.2, 0, 0.2, 0.4),
                     labels = c(-0.2, 0, 0.2, 0.4))+
  guides(shape = F,fill= F,
         color=guide_legend(override.aes=list(alpha=1,size=6,
                                              fill = c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"),
                                              shape = c(21,21,21,22,22,22,22)),
                            color="black",nrow=2))+
  facet_wrap2(~lakeID,scale="fixed",nrow=3,remove_labels = F,
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'lightgrey',
                #                               color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  labs(color="Season",
       x=paste0("Axis 1 [",BC_18S_out[[3]],"%]"),
       y=paste0("Axis 2 [",BC_18S_out[[4]],"%]"))
BC_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/BC_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Adonis ####
Euk_factor<-  data.frame(BC_18S_out[[1]]) %>%
  dplyr::mutate(.,
                season=rar_euk@sam_data$season[match(rownames(.),rownames(rar_euk@sam_data))],
                trophic_status = ifelse(lakeID %in% c("VSM", "JAB", "CER-L"),
                                        "oligo-mesotrophic",
                                        ifelse(lakeID %in% c("CSM", "VSS"),
                                               "hypereutrophic","meso-eutrophic")))

permanova_euk<-adonis2(BC_18S.dist ~ lakeID*season ,data=Euk_factor,
                       permutations=999,method="bray",by = "terms")
permanova_euk$`Pr(>F)`[[1]]
permanova_euk

permanova_euk<-adonis2(BC_18S.dist ~ trophic_status ,data=Euk_factor,
                       permutations=999,method="bray",by = "terms")
permanova_euk$`Pr(>F)`[[1]]
permanova_euk

mantel(BC_ABC.dist,JC_ABC.dist,method = "spearman", permutations = 999, na.rm = TRUE)

####__Phototrophs ####
tmp<-subset_taxa(rar_euk, trophic_mode == "phototrophs")

BC_phototroph.dist <- vegdist(t(tmp@otu_table),method = "bray")
BC_phototroph_out <- get_MDS_output(tmp@otu_table,
                                    BC_phototroph.dist,
                             2)

BC_phototroph_out[[1]] <- BC_phototroph_out[[1]] %>%
  cbind(.,
        lakeID = tmp@sam_data$lakeID[match(rownames(.),
                                               rownames(tmp@sam_data))],
        month_code = tmp@sam_data$month_code[match(rownames(.),
                                                       rownames(tmp@sam_data))],
        lake_season = tmp@sam_data$lake_season[match(rownames(.),
                                                         rownames(tmp@sam_data))],
        season_year = tmp@sam_data$season_year[match(rownames(.),
                                                         rownames(tmp@sam_data))],
        year = tmp@sam_data$year[match(rownames(.),
                                           rownames(tmp@sam_data))])

phototroph_factor<-  data.frame(BC_phototroph_out[[1]]) %>%
  dplyr::mutate(.,
                season=tmp@sam_data$season[match(rownames(.),rownames(tmp@sam_data))],
                trophic_status = ifelse(lakeID %in% c("VSM", "JAB", "CER-L"),
                                        "oligo-meso","eutro-hyper"))

permanova_phototroph<-adonis2(BC_phototroph.dist ~ lakeID*season ,data=phototroph_factor,
                       permutations=999,method="bray",by = "terms")
permanova_phototroph$`Pr(>F)`[[1]]
permanova_phototroph

####__Het. phagotrophs ####
tmp<-subset_taxa(rar_euk, trophic_mode == "heterotrophic phagotrophs")

BC_heterotroph.dist <- vegdist(t(tmp@otu_table),method = "bray")
BC_heterotroph_out <- get_MDS_output(tmp@otu_table,
                                    BC_heterotroph.dist,
                                    2)

BC_heterotroph_out[[1]] <- BC_heterotroph_out[[1]] %>%
  cbind(.,
        lakeID = tmp@sam_data$lakeID[match(rownames(.),
                                           rownames(tmp@sam_data))],
        month_code = tmp@sam_data$month_code[match(rownames(.),
                                                   rownames(tmp@sam_data))],
        lake_season = tmp@sam_data$lake_season[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        season_year = tmp@sam_data$season_year[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        year = tmp@sam_data$year[match(rownames(.),
                                       rownames(tmp@sam_data))])

heterotroph_factor<-  data.frame(BC_heterotroph_out[[1]]) %>%
  dplyr::mutate(.,
                season=tmp@sam_data$season[match(rownames(.),rownames(tmp@sam_data))],
                trophic_status = ifelse(lakeID %in% c("VSM", "JAB", "CER-L"),
                                        "oligo-meso","eutro-hyper"))

permanova_heterotroph<-adonis2(BC_heterotroph.dist ~ lakeID*season ,data=heterotroph_factor,
                              permutations=999,method="bray",by = "terms")
permanova_heterotroph$`Pr(>F)`[[1]]
permanova_heterotroph

####__Mixotrophs ####
tmp<-subset_taxa(rar_euk, trophic_mode == "mixotrophs")

BC_mixotroph.dist <- vegdist(t(tmp@otu_table),method = "bray")
BC_mixotroph_out <- get_MDS_output(tmp@otu_table,
                                     BC_mixotroph.dist,
                                     2)

BC_mixotroph_out[[1]] <- BC_mixotroph_out[[1]] %>%
  cbind(.,
        lakeID = tmp@sam_data$lakeID[match(rownames(.),
                                           rownames(tmp@sam_data))],
        month_code = tmp@sam_data$month_code[match(rownames(.),
                                                   rownames(tmp@sam_data))],
        lake_season = tmp@sam_data$lake_season[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        season_year = tmp@sam_data$season_year[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        year = tmp@sam_data$year[match(rownames(.),
                                       rownames(tmp@sam_data))])

mixotroph_factor<-  data.frame(BC_mixotroph_out[[1]]) %>%
  dplyr::mutate(.,
                season=tmp@sam_data$season[match(rownames(.),rownames(tmp@sam_data))],
                trophic_status = ifelse(lakeID %in% c("VSM", "JAB", "CER-L"),
                                        "oligo-meso","eutro-hyper"))

permanova_mixotroph<-adonis2(BC_mixotroph.dist ~ lakeID*season ,data=mixotroph_factor,
                               permutations=999,method="bray",by = "terms")
permanova_mixotroph$`Pr(>F)`[[1]]
permanova_mixotroph

####__Het. parasites ####
tmp<-subset_taxa(rar_euk, trophic_mode == "heterotrophic parasites")

BC_parasites.dist <- vegdist(t(tmp@otu_table),method = "bray")
BC_parasites_out <- get_MDS_output(tmp@otu_table,
                                   BC_parasites.dist,
                                     2)
 
BC_parasites_out[[1]] <- BC_heterotroph_out[[1]] %>%
  cbind(.,
        lakeID = tmp@sam_data$lakeID[match(rownames(.),
                                           rownames(tmp@sam_data))],
        month_code = tmp@sam_data$month_code[match(rownames(.),
                                                   rownames(tmp@sam_data))],
        lake_season = tmp@sam_data$lake_season[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        season_year = tmp@sam_data$season_year[match(rownames(.),
                                                     rownames(tmp@sam_data))],
        year = tmp@sam_data$year[match(rownames(.),
                                       rownames(tmp@sam_data))])

heterotroph_factor<-  data.frame(BC_heterotroph_out[[1]]) %>%
  dplyr::mutate(.,
                season=tmp@sam_data$season[match(rownames(.),rownames(tmp@sam_data))],
                trophic_status = ifelse(lakeID %in% c("VSM", "JAB", "CER-L"),
                                        "oligo-meso","eutro-hyper"))

permanova_heterotroph<-adonis2(BC_heterotroph.dist ~ lakeID*season ,data=heterotroph_factor,
                               permutations=999,method="bray",by = "terms")
permanova_heterotroph$`Pr(>F)`[[1]]
permanova_heterotroph

####|####

#### Time Lag Analysis ####

####__BC by day_gap ####
TLA_euk.df<- as.matrix(phyloseq::distance(rar_euk,method = "bray")) %>%
  melt(.) %>% rename(G1=Var1,G2=Var2) %>%
  filter(as.character(G1) != as.character(G2)) %>%
  mutate_if(is.factor,as.character) %>%
  dplyr::rename("BC"="value") %>%
  cbind(.,
        G1_day=rar_euk@sam_data$day_number[match(.$G1,
                                                  rownames(rar_euk@sam_data))],
        G2_day=rar_euk@sam_data$day_number[match(.$G2,
                                                  rownames(rar_euk@sam_data))],
        G1_lakeID=rar_euk@sam_data$lakeID[match(.$G1,
                                                 rownames(rar_euk@sam_data))],
        G2_lakeID=rar_euk@sam_data$lakeID[match(.$G2,
                                                 rownames(rar_euk@sam_data))],
        G1_column=rar_euk@sam_data$column[match(.$G1,
                                                 rownames(rar_euk@sam_data))],
        G2_column=rar_euk@sam_data$column[match(.$G2,
                                                 rownames(rar_euk@sam_data))]) %>%
  dplyr::mutate(day_gap=abs(.$G1_day-.$G2_day),
                comp_type=paste0(G1_lakeID,"_",G1_lakeID),
                comp_col=paste0(G1_column,"_",G2_column)) %>%
  subset(.,G1_day-G2_day >0 & G1_lakeID==G2_lakeID & G1_column==G2_column)

####__linear vs polynomial ####
list_lake<- unique(TLA_euk.df$G1_lakeID)

comp_model<-c("lakeID","sig_linear","r_linear",
              "AIC_linear","AIC_diff_poly2","AIC_diff_poly3","AIC_diff_poly4","AIC_diff_poly5",
              "sig_poly2","r_poly2","sig_poly3","r_poly3","sig_poly4","r_poly4","sig_poly5","r_poly5")

for (i in 1:length(unique(TLA_euk.df$G1_lakeID))) {
  # select lake
  lm_lake<-TLA_euk.df %>% subset(G1_lakeID %in% list_lake[i])
  # linea model sig, r and AIC score
  lm_linear<-lm(lm_lake,formula=BC ~ day_gap) %>% summary()
  sig_linear<-round(as.numeric(pf(lm_linear$fstatistic[1],
                                  lm_linear$fstatistic[2],
                                  lm_linear$fstatistic[3],
                                  lower.tail = FALSE)),4)
  r_linear<-round(as.numeric(lm_linear$adj.r.squared),3)*100
  AIC_linear<-round(AIC(lm(lm_lake,formula=BC ~ day_gap)),1)
  # poly2 model sig, r and AIC score
  lm_poly2<-lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2)) %>% summary()
  sig_poly2<-round(as.numeric(pf(lm_poly2$fstatistic[1],
                                 lm_poly2$fstatistic[2],
                                 lm_poly2$fstatistic[3],
                                 lower.tail = FALSE)),4)
  r_poly2<-round(as.numeric(lm_poly2$adj.r.squared),3)*100
  AIC_poly2<-round(AIC(lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2))),1)
  # poly3 model sig, r and AIC score
  lm_poly3<-lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3)) %>% summary()
  sig_poly3<-round(as.numeric(pf(lm_poly3$fstatistic[1],
                                 lm_poly3$fstatistic[2],
                                 lm_poly3$fstatistic[3],
                                 lower.tail = FALSE)),4)
  r_poly3<-round(as.numeric(lm_poly3$adj.r.squared),3)*100
  AIC_poly3<-round(AIC(lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3))),1)
  # poly4 model sig, r and AIC score
  lm_poly4<-lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3)+I(day_gap^4)) %>% summary()
  sig_poly4<-round(as.numeric(pf(lm_poly4$fstatistic[1],
                                 lm_poly4$fstatistic[2],
                                 lm_poly4$fstatistic[3],
                                 lower.tail = FALSE)),4)
  r_poly4<-round(as.numeric(lm_poly4$adj.r.squared),3)*100
  AIC_poly4<-round(AIC(lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3)+I(day_gap^4))),1)
  # poly5 model sig, r and AIC score
  lm_poly5<-lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3)+I(day_gap^4)) %>% summary()
  sig_poly5<-round(as.numeric(pf(lm_poly5$fstatistic[1],
                                 lm_poly5$fstatistic[2],
                                 lm_poly5$fstatistic[3],
                                 lower.tail = FALSE)),4)
  r_poly5<-round(as.numeric(lm_poly5$adj.r.squared),3)*100
  AIC_poly5<-round(AIC(lm(lm_lake,formula=BC ~ day_gap + I(day_gap^2) +I(day_gap^3)+I(day_gap^4)+I(day_gap^5))),1)
  # diff AIC score in percent
  AIC_diff_poly2<-round(((AIC_poly2-AIC_linear)/AIC_poly2)*100,2)
  AIC_diff_poly3<-round(((AIC_poly3-AIC_linear)/AIC_poly3)*100,2)
  AIC_diff_poly4<-round(((AIC_poly4-AIC_linear)/AIC_poly4)*100,2)
  AIC_diff_poly5<-round(((AIC_poly5-AIC_linear)/AIC_poly5)*100,2)
  tmp<-c(list_lake[i],sig_linear,r_linear,
         AIC_linear,AIC_diff_poly2,AIC_diff_poly3,AIC_diff_poly4,AIC_diff_poly5,
         sig_poly2,r_poly2,sig_poly3,r_poly3,sig_poly4,r_poly4,sig_poly5,r_poly5)
  comp_model<-rbind(comp_model,tmp)
}

comp_model<-
  comp_model %>% as.data.frame() %>% row_to_names(row_number = 1) %>%
  remove_rownames()

write.csv(comp_model,"output/model_comp_18S.csv",row.names = F)

####__plot ####
TLA_euk.df <- TLA_euk.df %>%
  dplyr::mutate(.,
                poly=ifelse(G1_lakeID %in% c("CSM","CRE","JAB","CER-L"), "poly5",
                            ifelse(G1_lakeID %in% c("CER-S","VSM"), "poly3",
                                   ifelse(G1_lakeID %in% c("LGP","VSS","BLR"),"poly4","linear"))),
                G1_lakeID=factor(G1_lakeID,
                                 levels=c("VSM", "JAB", "CER-L",
                                          "CER-S", "CRE", "BLR","LGP", "CSM", 
                                          "VSS")))
TLA_euk <- TLA_euk.df %>%
  ggplot(.,aes(x=day_gap,y=BC,group=G1_lakeID,label=poly))+
  #timeframe
  geom_rect(xmin = 85, xmax=95, ymin = 0, ymax = 1,
            color = "#E5E4E2", fill = "#E5E4E2",
            show.legend = F)+
  geom_rect(xmin = 175, xmax=185, ymin = 0, ymax = 1,
            color = "#E5E4E2", fill = "#E5E4E2",
            show.legend = F)+
  geom_rect(xmin = 360, xmax=370, ymin = 0, ymax = 1,
            color = "#E5E4E2", fill = "#E5E4E2",
            show.legend = F)+
  #points
  geom_point(size=1,alpha=0.4,color="gray30",show.legend = F)+
  #linear model
  geom_text(x=250,y=0.12,size=3,check_overlap = T,fontface="italic",label="Linear",color="gray32")+
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,formula = lm_formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,formula = lm_formula_linear,
               label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #best poly model poly3
  geom_text(x=250,y=0.04,size=3,
            label= "Polynomial", check_overlap = T,
            fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  stat_poly_line(formula = lm_formula_poly3,show.legend = F,
                 color = "blue", fill = "lightblue",alpha=0.4,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.02, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  #best poly model poly4
  geom_text(x=250,y=0.04,size=3,
            label= "Polynomial", check_overlap = T,
            fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  stat_poly_line(formula = lm_formula_poly4,show.legend = F,
                 color = "blue", fill = "lightblue",alpha=0.4,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.02, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  #best poly model poly5
  geom_text(x=250,y=0.04,size=3,
            label= "Polynomial", check_overlap = T,
            fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  stat_poly_line(formula = lm_formula_poly5,show.legend = F,
                 color = "blue", fill = "lightblue",alpha=0.4,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly5,show.legend = F,
               label.y = 0.02, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.ticks = element_line(color="black"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     label = c("0","0.25","0.5","0.75","1"),
                     expand = c(0,0),limits = c(0,1))+
  scale_x_continuous(breaks = c(30, 90, 180, 365, 540),
                     expand = c(0,0),limits = c(-1,545))+
  labs(y="Bray-curtis dissimilarity",
       x="Time between sampling dates (days)")+
  facet_wrap2(~ G1_lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'white',
                #                               color = NA),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))
TLA_euk

ggsave("/Users/piefouca/Desktop/µEuk/Figures/TLA_euk.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__dist ####

comm.df<-
  as.matrix(phyloseq::distance(rar_euk,method = "bray")) %>%
  melt(.) %>% rename(G1=Var1,G2=Var2) %>%
  filter(as.character(G1) != as.character(G2)) %>%
  mutate_if(is.factor,as.character) %>%
  dplyr::rename("BC"="value") %>%
  cbind(.,
        G1_month=
          rar_euk@sam_data$month_code[match(.$G1,
                                            rownames(rar_euk@sam_data))],
        G2_month=
          rar_euk@sam_data$month_code[match(.$G2,
                                            rownames(rar_euk@sam_data))],
        G1_day=rar_euk@sam_data$day_number[match(.$G1,
                                                 rownames(rar_euk@sam_data))],
        G2_day=rar_euk@sam_data$day_number[match(.$G2,
                                                 rownames(rar_euk@sam_data))],
        G1_lakeID=rar_euk@sam_data$lakeID[match(.$G1,
                                                rownames(rar_euk@sam_data))],
        G2_lakeID=rar_euk@sam_data$lakeID[match(.$G2,
                                                rownames(rar_euk@sam_data))],
        G1_column=rar_euk@sam_data$column[match(.$G1,
                                                rownames(rar_euk@sam_data))],
        G2_column=rar_euk@sam_data$column[match(.$G2,
                                                rownames(rar_euk@sam_data))]) %>%
  dplyr::mutate(day_gap=abs(.$G1_day-.$G2_day),
                comp_type=paste0(G1_lakeID,"_",G1_lakeID),
                comp_col=paste0(G1_column,"_",G2_column)) %>%
  subset(.,G1_day-G2_day >0 & G1_lakeID==G2_lakeID & G1_column==G2_column)

  as.data.frame(TLA_euk.df %>%
  dplyr::mutate(lake_col=paste0(G1_lakeID,"_",G1_column)) %>%
  select(c("lake_col","day_gap","BC"))) %>%
  group_by(lake_col,day_gap) %>% summarise(BC=median(BC)) %>%
  pivot_wider(names_from = day_gap,values_from = BC) %>%
  select(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`) %>% remove_rownames() %>%
  column_to_rownames(var = "lake_ID")

clust.com <- tsclust(tmp_com, type="partitional", k=2L:5L, distance="dtw",centroid = "pam")

names(clust.com) <- paste0("k_", 2L:5L)
sapply(clust.com, cvi,type="DB")

clust.com.hier <- tsclust(tmp_com, type = "h", k =7L, distance = "dtw")
clust.com.hier_clust<-data.frame(cluster=cutree(clust.com.hier, k=4L))
plot(clust.com.hier, type="series")

dend <- clust.com.hier %>% as.dendrogram

com_dendro_data <- dendro_data(dend, type = "rectangle")
com_dendro_data$labels

Fig_com_dendro <- ggplot(com_dendro_data$segments) +
  annotate("rect", xmin = -1.6, xmax = -1.8, ymin = 0.8, ymax = 3.2, fill = "#2F6B9D",color="black")+
  annotate("rect", xmin = -1.6, xmax = 7.8, ymin = 0.8, ymax = 3.2,fill="#2F6B9D",alpha=0.2)+
  annotate("rect", xmin = -1.6, xmax = -1.8, ymin = 3.8, ymax = 8.2, fill = "#28A448",color="black")+
  annotate("rect", xmin = -1.6, xmax = 7.8, ymin = 3.8, ymax = 8.2,fill="#28A448",alpha=0.2)+
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend))+
  geom_text(data = com_dendro_data$labels, aes(-0.1,x, label = label),
            hjust =0, angle = 0, size = 4,fontface=2)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_x_reverse(limits=c(9,-2.5))+
  scale_y_continuous(limits = c(0,9),expand = c(0,0))

Fig_com_dendro

wilcox_test(chla ~ cluster, data = Fig_com_chla$data,
            exact = FALSE)

Fig_com_chla <- data.frame(lake_ID=clust.com.hier$labels) %>%
  cbind(.,chla=round(Chla_2022$chla_mean[match(.$lake_ID,Chla_2022$lake_ID)],1),
        subcluster=clust.com.hier_clust$cluster[match(.$lake_ID,rownames(clust.com.hier_clust))]) %>%
  dplyr::mutate(cluster=ifelse(subcluster ==1,"cluster 1","cluster 2")) %>%
  ggplot(.,aes(as.factor(cluster),chla,fill=cluster))+
  #geom_text(mapping = aes(x=50,y=1),label="*",size=3,inherit.aes = F)+
  geom_boxplot(show.legend = F,alpha=0.2,width=0.5)+
  geom_jitter(width = 0.1,show.legend = F,size=2)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.line.x.bottom = element_line(color = 'black'),
        panel.border = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = c("#2F6B9D","#28A448"))+
  scale_y_continuous(breaks = c(5,10,20,50,100,200),limits = c(2.6,NA),trans = 'log')+
  annotate("text", label = "*",x = 1.5, y = 50, size = 12, colour = "darkred")+
  labs(y=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))))+
  coord_flip()
Fig_com_chla

####|####

#### MOTA ####

####__Trajectory length ####
#56 Axis are needed to explain 90.17%
MOTA_lake_18S <- final_cum_distance_18S %>%
  melt(.,"time_label") %>% as.data.frame() %>%
  dplyr::rename("lakeID"="variable") %>%
  dplyr::mutate(lake_month=paste0(lakeID,"_",time_label),
                season_year = rar_euk@sam_data$season_year[match(.$time_label,
                                                       rar_euk@sam_data$month_number)],
                season_year=factor(season_year,levels = c("summer_2021","fall_2021","winter_2021",
                                                      "spring_2022","summer_2022","fall_2022","winter_2022")),
                lakeID=factor(lakeID,levels=rev(c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS")))) %>%
  group_by(lakeID) %>%
  dplyr::mutate(.,max_dist=max(value)) %>%
  ggplot(.,
         aes(x=value,y=lakeID,label=time_label,fill=season_year, group = lakeID))+
  geom_line(size=0.5,color="black",show.legend = F )+
  #color VSS
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$VSS[2]+((final_cum_distance_18S$VSS[3]-final_cum_distance_18S$VSS[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[2]+((final_cum_distance_18S$VSS[3]-final_cum_distance_18S$VSS[2])/2)),
           xmax=(final_cum_distance_18S$VSS[5]+((final_cum_distance_18S$VSS[6]-final_cum_distance_18S$VSS[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[5]+((final_cum_distance_18S$VSS[6]-final_cum_distance_18S$VSS[5])/2)),
           xmax=(final_cum_distance_18S$VSS[7]+((final_cum_distance_18S$VSS[8]-final_cum_distance_18S$VSS[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[7]+((final_cum_distance_18S$VSS[8]-final_cum_distance_18S$VSS[7])/2)),
           xmax=(final_cum_distance_18S$VSS[10]+((final_cum_distance_18S$VSS[11]-final_cum_distance_18S$VSS[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[10]+((final_cum_distance_18S$VSS[11]-final_cum_distance_18S$VSS[10])/2)),
           xmax=(final_cum_distance_18S$VSS[13]+((final_cum_distance_18S$VSS[14]-final_cum_distance_18S$VSS[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[13]+((final_cum_distance_18S$VSS[14]-final_cum_distance_18S$VSS[13])/2)),
           xmax=(final_cum_distance_18S$VSS[16]+((final_cum_distance_18S$VSS[17]-final_cum_distance_18S$VSS[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSS[16]+((final_cum_distance_18S$VSS[17]-final_cum_distance_18S$VSS[16])/2)),
           xmax=max(final_cum_distance_18S$VSS)+0.23)+
  #color CSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$CSM[2]+((final_cum_distance_18S$CSM[3]-final_cum_distance_18S$CSM[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[2]+((final_cum_distance_18S$CSM[3]-final_cum_distance_18S$CSM[2])/2)),
           xmax=(final_cum_distance_18S$CSM[5]+((final_cum_distance_18S$CSM[6]-final_cum_distance_18S$CSM[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[5]+((final_cum_distance_18S$CSM[6]-final_cum_distance_18S$CSM[5])/2)),
           xmax=(final_cum_distance_18S$CSM[7]+((final_cum_distance_18S$CSM[8]-final_cum_distance_18S$CSM[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[7]+((final_cum_distance_18S$CSM[8]-final_cum_distance_18S$CSM[7])/2)),
           xmax=(final_cum_distance_18S$CSM[10]+((final_cum_distance_18S$CSM[11]-final_cum_distance_18S$CSM[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[10]+((final_cum_distance_18S$CSM[11]-final_cum_distance_18S$CSM[10])/2)),
           xmax=(final_cum_distance_18S$CSM[13]+((final_cum_distance_18S$CSM[14]-final_cum_distance_18S$CSM[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[13]+((final_cum_distance_18S$CSM[14]-final_cum_distance_18S$CSM[13])/2)),
           xmax=(final_cum_distance_18S$CSM[16]+((final_cum_distance_18S$CSM[17]-final_cum_distance_18S$CSM[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CSM[16]+((final_cum_distance_18S$CSM[17]-final_cum_distance_18S$CSM[16])/2)),
           xmax=max(final_cum_distance_18S$CSM)+0.23)+
  #color LGP
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$LGP[2]+((final_cum_distance_18S$LGP[3]-final_cum_distance_18S$LGP[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[2]+((final_cum_distance_18S$LGP[3]-final_cum_distance_18S$LGP[2])/2)),
           xmax=(final_cum_distance_18S$LGP[5]+((final_cum_distance_18S$LGP[6]-final_cum_distance_18S$LGP[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[5]+((final_cum_distance_18S$LGP[6]-final_cum_distance_18S$LGP[5])/2)),
           xmax=(final_cum_distance_18S$LGP[7]+((final_cum_distance_18S$LGP[8]-final_cum_distance_18S$LGP[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[7]+((final_cum_distance_18S$LGP[8]-final_cum_distance_18S$LGP[7])/2)),
           xmax=(final_cum_distance_18S$LGP[10]+((final_cum_distance_18S$LGP[11]-final_cum_distance_18S$LGP[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[10]+((final_cum_distance_18S$LGP[11]-final_cum_distance_18S$LGP[10])/2)),
           xmax=(final_cum_distance_18S$LGP[13]+((final_cum_distance_18S$LGP[14]-final_cum_distance_18S$LGP[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[13]+((final_cum_distance_18S$LGP[14]-final_cum_distance_18S$LGP[13])/2)),
           xmax=(final_cum_distance_18S$LGP[16]+((final_cum_distance_18S$LGP[17]-final_cum_distance_18S$LGP[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_18S$LGP[16]+((final_cum_distance_18S$LGP[17]-final_cum_distance_18S$LGP[16])/2)),
           xmax=max(final_cum_distance_18S$LGP)+0.23)+
  #color BLR
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$BLR[2]+((final_cum_distance_18S$BLR[3]-final_cum_distance_18S$BLR[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[2]+((final_cum_distance_18S$BLR[3]-final_cum_distance_18S$BLR[2])/2)),
           xmax=(final_cum_distance_18S$BLR[5]+((final_cum_distance_18S$BLR[6]-final_cum_distance_18S$BLR[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[5]+((final_cum_distance_18S$BLR[6]-final_cum_distance_18S$BLR[5])/2)),
           xmax=(final_cum_distance_18S$BLR[7]+((final_cum_distance_18S$BLR[8]-final_cum_distance_18S$BLR[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[7]+((final_cum_distance_18S$BLR[8]-final_cum_distance_18S$BLR[7])/2)),
           xmax=(final_cum_distance_18S$BLR[10]+((final_cum_distance_18S$BLR[11]-final_cum_distance_18S$BLR[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[10]+((final_cum_distance_18S$BLR[11]-final_cum_distance_18S$BLR[10])/2)),
           xmax=(final_cum_distance_18S$BLR[13]+((final_cum_distance_18S$BLR[14]-final_cum_distance_18S$BLR[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[13]+((final_cum_distance_18S$BLR[14]-final_cum_distance_18S$BLR[13])/2)),
           xmax=(final_cum_distance_18S$BLR[16]+((final_cum_distance_18S$BLR[17]-final_cum_distance_18S$BLR[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_18S$BLR[16]+((final_cum_distance_18S$BLR[17]-final_cum_distance_18S$BLR[16])/2)),
           xmax=max(final_cum_distance_18S$BLR)+0.23)+
  #color CRE
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$CRE[2]+((final_cum_distance_18S$CRE[3]-final_cum_distance_18S$CRE[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[2]+((final_cum_distance_18S$CRE[3]-final_cum_distance_18S$CRE[2])/2)),
           xmax=(final_cum_distance_18S$CRE[5]+((final_cum_distance_18S$CRE[6]-final_cum_distance_18S$CRE[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[5]+((final_cum_distance_18S$CRE[6]-final_cum_distance_18S$CRE[5])/2)),
           xmax=(final_cum_distance_18S$CRE[7]+((final_cum_distance_18S$CRE[8]-final_cum_distance_18S$CRE[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[7]+((final_cum_distance_18S$CRE[8]-final_cum_distance_18S$CRE[7])/2)),
           xmax=(final_cum_distance_18S$CRE[10]+((final_cum_distance_18S$CRE[11]-final_cum_distance_18S$CRE[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[10]+((final_cum_distance_18S$CRE[11]-final_cum_distance_18S$CRE[10])/2)),
           xmax=(final_cum_distance_18S$CRE[13]+((final_cum_distance_18S$CRE[14]-final_cum_distance_18S$CRE[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[13]+((final_cum_distance_18S$CRE[14]-final_cum_distance_18S$CRE[13])/2)),
           xmax=(final_cum_distance_18S$CRE[16]+((final_cum_distance_18S$CRE[17]-final_cum_distance_18S$CRE[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_18S$CRE[16]+((final_cum_distance_18S$CRE[17]-final_cum_distance_18S$CRE[16])/2)),
           xmax=max(final_cum_distance_18S$CRE)+0.23)+
  #color `CER-S`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$`CER-S`[2]+((final_cum_distance_18S$`CER-S`[3]-final_cum_distance_18S$`CER-S`[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[2]+((final_cum_distance_18S$`CER-S`[3]-final_cum_distance_18S$`CER-S`[2])/2)),
           xmax=(final_cum_distance_18S$`CER-S`[5]+((final_cum_distance_18S$`CER-S`[6]-final_cum_distance_18S$`CER-S`[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[5]+((final_cum_distance_18S$`CER-S`[6]-final_cum_distance_18S$`CER-S`[5])/2)),
           xmax=(final_cum_distance_18S$`CER-S`[7]+((final_cum_distance_18S$`CER-S`[8]-final_cum_distance_18S$`CER-S`[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[7]+((final_cum_distance_18S$`CER-S`[8]-final_cum_distance_18S$`CER-S`[7])/2)),
           xmax=(final_cum_distance_18S$`CER-S`[10]+((final_cum_distance_18S$`CER-S`[11]-final_cum_distance_18S$`CER-S`[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[10]+((final_cum_distance_18S$`CER-S`[11]-final_cum_distance_18S$`CER-S`[10])/2)),
           xmax=(final_cum_distance_18S$`CER-S`[13]+((final_cum_distance_18S$`CER-S`[14]-final_cum_distance_18S$`CER-S`[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[13]+((final_cum_distance_18S$`CER-S`[14]-final_cum_distance_18S$`CER-S`[13])/2)),
           xmax=(final_cum_distance_18S$`CER-S`[16]+((final_cum_distance_18S$`CER-S`[17]-final_cum_distance_18S$`CER-S`[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-S`[16]+((final_cum_distance_18S$`CER-S`[17]-final_cum_distance_18S$`CER-S`[16])/2)),
           xmax=max(final_cum_distance_18S$`CER-S`)+0.23)+
  #color `CER-L`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$`CER-L`[2]+((final_cum_distance_18S$`CER-L`[3]-final_cum_distance_18S$`CER-L`[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[2]+((final_cum_distance_18S$`CER-L`[3]-final_cum_distance_18S$`CER-L`[2])/2)),
           xmax=(final_cum_distance_18S$`CER-L`[5]+((final_cum_distance_18S$`CER-L`[6]-final_cum_distance_18S$`CER-L`[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[5]+((final_cum_distance_18S$`CER-L`[6]-final_cum_distance_18S$`CER-L`[5])/2)),
           xmax=(final_cum_distance_18S$`CER-L`[7]+((final_cum_distance_18S$`CER-L`[8]-final_cum_distance_18S$`CER-L`[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[7]+((final_cum_distance_18S$`CER-L`[8]-final_cum_distance_18S$`CER-L`[7])/2)),
           xmax=(final_cum_distance_18S$`CER-L`[10]+((final_cum_distance_18S$`CER-L`[11]-final_cum_distance_18S$`CER-L`[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[10]+((final_cum_distance_18S$`CER-L`[11]-final_cum_distance_18S$`CER-L`[10])/2)),
           xmax=(final_cum_distance_18S$`CER-L`[13]+((final_cum_distance_18S$`CER-L`[14]-final_cum_distance_18S$`CER-L`[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[13]+((final_cum_distance_18S$`CER-L`[14]-final_cum_distance_18S$`CER-L`[13])/2)),
           xmax=(final_cum_distance_18S$`CER-L`[16]+((final_cum_distance_18S$`CER-L`[17]-final_cum_distance_18S$`CER-L`[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_18S$`CER-L`[16]+((final_cum_distance_18S$`CER-L`[17]-final_cum_distance_18S$`CER-L`[16])/2)),
           xmax=max(final_cum_distance_18S$`CER-L`)+0.23)+
  #color JAB
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$JAB[2]+((final_cum_distance_18S$JAB[3]-final_cum_distance_18S$JAB[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[2]+((final_cum_distance_18S$JAB[3]-final_cum_distance_18S$JAB[2])/2)),
           xmax=(final_cum_distance_18S$JAB[5]+((final_cum_distance_18S$JAB[6]-final_cum_distance_18S$JAB[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[5]+((final_cum_distance_18S$JAB[6]-final_cum_distance_18S$JAB[5])/2)),
           xmax=(final_cum_distance_18S$JAB[7]+((final_cum_distance_18S$JAB[8]-final_cum_distance_18S$JAB[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[7]+((final_cum_distance_18S$JAB[8]-final_cum_distance_18S$JAB[7])/2)),
           xmax=(final_cum_distance_18S$JAB[10]+((final_cum_distance_18S$JAB[11]-final_cum_distance_18S$JAB[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[10]+((final_cum_distance_18S$JAB[11]-final_cum_distance_18S$JAB[10])/2)),
           xmax=(final_cum_distance_18S$JAB[13]+((final_cum_distance_18S$JAB[14]-final_cum_distance_18S$JAB[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[13]+((final_cum_distance_18S$JAB[14]-final_cum_distance_18S$JAB[13])/2)),
           xmax=(final_cum_distance_18S$JAB[16]+((final_cum_distance_18S$JAB[17]-final_cum_distance_18S$JAB[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_18S$JAB[16]+((final_cum_distance_18S$JAB[17]-final_cum_distance_18S$JAB[16])/2)),
           xmax=max(final_cum_distance_18S$JAB)+0.23)+
  #color VSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=-0.16,
           xmax=(final_cum_distance_18S$VSM[2]+((final_cum_distance_18S$VSM[3]-final_cum_distance_18S$VSM[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[2]+((final_cum_distance_18S$VSM[3]-final_cum_distance_18S$VSM[2])/2)),
           xmax=(final_cum_distance_18S$VSM[5]+((final_cum_distance_18S$VSM[6]-final_cum_distance_18S$VSM[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[5]+((final_cum_distance_18S$VSM[6]-final_cum_distance_18S$VSM[5])/2)),
           xmax=(final_cum_distance_18S$VSM[7]+((final_cum_distance_18S$VSM[8]-final_cum_distance_18S$VSM[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[7]+((final_cum_distance_18S$VSM[8]-final_cum_distance_18S$VSM[7])/2)),
           xmax=(final_cum_distance_18S$VSM[10]+((final_cum_distance_18S$VSM[11]-final_cum_distance_18S$VSM[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[10]+((final_cum_distance_18S$VSM[11]-final_cum_distance_18S$VSM[10])/2)),
           xmax=(final_cum_distance_18S$VSM[13]+((final_cum_distance_18S$VSM[14]-final_cum_distance_18S$VSM[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[13]+((final_cum_distance_18S$VSM[14]-final_cum_distance_18S$VSM[13])/2)),
           xmax=(final_cum_distance_18S$VSM[16]+((final_cum_distance_18S$VSM[17]-final_cum_distance_18S$VSM[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_18S$VSM[16]+((final_cum_distance_18S$VSM[17]-final_cum_distance_18S$VSM[16])/2)),
           xmax=max(final_cum_distance_18S$VSM)+0.23)+
  geom_line(size=0.5,color="black",show.legend = F)+
  geom_point(show.legend = F)+
  geom_label(show.legend = F)+
  theme_bw()+
  theme(#aspect.ratio=1,
        panel.grid = element_blank(),
        plot.caption =element_text(face="italic",size=10,hjust=1),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_blank(),
        # axis.text.y=element_text(size=12,face="bold",
        #                          color=rev(c("#28A448","#E45000","#980000","#1E7A36","#EB7947","#F2AB8C","#E40000","#2F6B9D","#5C9BD6"))),
        axis.text.y=element_text(size=12,face="bold",color="black"),
        axis.text.x=element_text(size=12, color="black"),
        axis.ticks = element_line(color="black"),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"),
                    labels = c("Summer 2021","Fall 2021","Winter 2021",
                               "Spring 2022","Summer 2022","Fall 2022","Winter 2022"))+
  scale_x_continuous(expand = c(0,0),limits = c(-0.25,10.6),
                     breaks = c(0, 2, 4, 6, 8, 10),
                     labels = c("0", "2", "4", "6", "8", "10"))+
  labs(#caption=paste0("MOTA – Bray-Curtis – ",Threshold*100,"% of the total expl. var. (58 axis)"),
    x="Trajectory length",fill="Season")+
  guides(fill=guide_legend(override.aes=list(shape=21,size=8,color="black"),nrow=2))
MOTA_lake_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/MOTA_lake_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__ MOTA vs Chla ####
delta_chla <- read_delim("cum_delta_chla.csv",";",show_col_types = F)

MOTA_delta_chla$data %>%
  cbind(.,
                month_number=
                  delta_chla$month_number[match(.$lake_month,
                                                delta_chla$lake_month_number)])
corr_chla_MOTA <- MOTA_lake_18S$data %>%
  ungroup() %>%
  select(c(time_label,lakeID,value)) %>%
  pivot_wider(names_from = lakeID, 
              values_from = c(value)) %>%
  cbind(data.frame(delta_chla %>%
                    select(c(month_number,lakeID,delta_Chla)) %>%
                    dplyr::mutate(.,lakeID = paste0(lakeID,"_chla")) %>%
                    pivot_wider(names_from = lakeID, 
                                values_from = delta_Chla)))

MOTA_delta_chla <- 
  MOTA_lake_18S$data %>%
  ungroup() %>%
  select(c("lakeID","lake_month", "value")) %>%
  rename("MOTA"="value") %>%
  dplyr::mutate(.,
                chla_range=chla_range.df$range[match(.$lakeID,chla_range.df$lakeID)],
                chla_max=chla_range.df$chla_max[match(.$lakeID,chla_range.df$lakeID)],
                chla_mean=chla_mean.df$chla_mean[match(.$lakeID,chla_mean.df$lakeID)],
                delta_chla=delta_chla$delta_Chla[match(.$lake_month,delta_chla$lake_month_number)],
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS"))) %>%
  group_by(lakeID) %>% summarise(MOTA_max=max(MOTA),
                                 delta_chla_max=max(delta_chla),
                                 chla_max=max(chla_max),
                                 chla_mean=max(chla_mean),
                                 chla_range=max(chla_range)) %>%
  ggplot(.,aes(chla_mean, MOTA_max, group = lakeID, fill = lakeID, color = lakeID)) +
  geom_point(show.legend = F,shape=21,color="black", size= 3)+
  theme_bw()+
  theme(
   # aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size=14,face = "bold"),
    axis.text.x = element_text(size=12, color="black",hjust = 0.2),
    axis.text.y = element_text(size=12, color="black",vjust = 0.2),
    axis.ticks = element_line(color="black"),
    legend.title = element_text(size=12,face = "bold",hjust=0),
    legend.text = element_text(size=10),
    legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=rev(palette_lake_chla))+
  scale_color_manual(values=rev(palette_lake_chla))+
  scale_x_continuous(expand = c(0,0),limits = c(2.6,65),
                     trans = "log1p",
                     breaks = c(2.6,7.3,56),
                     labels = c("2.6", "7.3", "56"))+
  scale_y_continuous(expand = c(0,0),limits = c(6,10.6),
                     breaks = c(6,8,10),
                     labels = c("6", "8", "10"))+
  geom_text(label=expression(paste(italic("p"),"<0.005 ",rho," 0.9")),
            color="black",aes(x=15,y=6.2),size=4.5,
            hjust=0,check_overlap = T,
            inherit.aes = F)+
  labs(y="Trajectory length",
       x=expression(paste(bold("Chl"),bold(italic("a")),bold(" mean ["),bold(italic("µ")),bold(g.L^"-1"),bold("]"))))+
  guides(fill=guide_legend(override.aes=list(shape=22,size=8),nrow=1),
         color = FALSE)

MOTA_delta_chla

cor.test(MOTA_delta_chla$data$MOTA_max,MOTA_delta_chla$data$chla_mean,method = "spearman")
# delta_chla_max 0.01722 0.7833333
# chla_max 0.003075 0.8833333 S = 14
# chla_mean 0.002028 0.9
# chla_range 0.003075 0.8833333 S = 14

design_MOTA_chla<-"
AAAA#BB"

MOTA_lake_18S+MOTA_delta_chla+
  plot_layout(desig=design_MOTA_chla,guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold",size=15),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")

ggsave("/Users/piefouca/Desktop/µEuk/Figures/MOTA_delta_chla_18S.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####________________________####