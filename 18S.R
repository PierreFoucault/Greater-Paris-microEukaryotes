#### Import data ####
 
####__18S-PR2 (rarefied at 12,287 reads) ####

rar_euk <-
  create_phyloseq_PR2("18S_data/PR2_ASV_table.tsv",
                  "18S_data/PR2_ASV-tax.tsv",
                  "18S_data/metadata_18S.txt")
#VSM_B_W1 removed during the preprocess (rarefaction)

####__trophic mode ####
list_phototroph_subdivision <- c("Gyrista","Chlorophyta_X", "Haptophyta_X")


list_heterotroph_subdivision <- c("Cryptophyta_X", "Dinoflagellata")

list_mixotroph_subdivision <- c("Ciliophora","Kathablepharida","Fungi",
                                "Cercozoa","Bigyra","Choanoflagellata",
                                "Chrompodellids",
                                "Perkinsea","Telonemia_X")

list_parasite_subdivision <- c("Perkinsea","Telonemia_X")

rar_euk@tax_table<-
  tax_table(data.frame(tax_table(rar_euk)) %>%
              dplyr::mutate(.,
                            trophic_mode = ifelse(Subdivision %in% list_phototroph_subdivision,"Phototroph",
                                                  ifelse(Subdivision %in% list_heterotroph_subdivision,"Heterotroph",
                                                         ifelse(Subdivision %in% list_mixotroph_subdivision,"Mixotroph",
                                                                ifelse(Subdivision %in% list_parasite_subdivision,"Parasite",NA))))) %>%
                                                           
              as.matrix(.))

write_tsv(as.data.frame(rar_euk@tax_table) %>% rownames_to_column(var = "QIIME_ID"),
          "18S_data/18S_mode-tax.tsv")

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
                day_number=
                  (rar_euk@sam_data$day_number[match(.$sample_code,rar_euk@sam_data$sample_code)]))
#View(alphadiv_18S)
write.csv(alphadiv_18S,"output/alphadiv_18S.csv")

alphadiv_18S %>%
  group_by(lakeID) %>%
  summarise(Richness_mean=mean(Richness),
            Evenness_mean=mean(Evenness),
            Shannon_mean=mean(Shannon))

####______Richness ####
richness_18S<- alphadiv_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Richness, color= lakeID, fill= lakeID, group=lakeID))+
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 600),
                     breaks = c(0, 200, 400, 600))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Eukaryota ASV Richness')
richness_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

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

ggsave("/Users/piefouca/Desktop/µEuk/Figures/evenness_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####______Shannon ####
shannon_18S<- alphadiv_18S %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(., aes(x=day_number, y=Shannon, color= lakeID, fill= lakeID, group=lakeID))+
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
  scale_y_continuous(expand = c(0,0), limits = c(0,5.5),
                     breaks = c(0, 1, 3, 5))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Eukaryota ASV Shannon')
shannon_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/shannon_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Phototrophs ####
alphadiv_phototroph_18S <- subset_taxa(rar_euk, Subdivision %in% list_phototroph_subdivision)

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
  scale_y_continuous(expand = c(0,0), limits = c(0, 200),
                     breaks = c(0, 100, 200))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Phototrophs ASV Richness')
richness_phototroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_phototroph_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Heterotrophs ####
alphadiv_heterotroph_18S <- subset_taxa(rar_euk, Subdivision %in% list_heterotroph_subdivision)

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
  scale_y_continuous(expand = c(0,0), limits = c(0, 100),
                     breaks = c(0, 50, 100))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Heterotroph ASV Richness')
richness_heterotroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_heterotroph_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Mixotrophs ####
alphadiv_mixotroph_18S <- subset_taxa(rar_euk, Subdivision %in% list_mixotroph_subdivision)

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
  scale_y_continuous(expand = c(0,0), limits = c(0, 350),
                     breaks = c(0, 150, 300))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Mixotroph ASV Richness')
richness_mixotroph_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_mixotroph_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Parasite ####
alphadiv_parasite_18S <- subset_taxa(rar_euk, Subdivision %in% list_parasite_subdivision)

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
  scale_y_continuous(expand = c(0,0), limits = c(0, 20),
                     breaks = c(0, 5, 10))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Parasite ASV Richness')
richness_parasite_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_parasite_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### | ####

#### Composition ####

####__season-Sub. ####
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
                                            "Perkinsea","Telonemia",'Ichthyosporea',
                                            'Taxa < 1%')),
                season_year = factor(season_year, levels = c("summer_2021","fall_2021","winter_2021",
                                                             "spring_2022","summer_2022","fall_2022","winter_2022"))) %>%
  ggplot(.,aes(x=season_year,y=median_abundance,fill=taxa_legend))+
  theme_bw()+
  geom_col(width = .95,color = "black", linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid = element_blank())+
  scale_x_discrete(expand = c(0,0))+ #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0),
                     labels = c('0%','25%','50%','75%','100%'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Subdivision)+ #set a color palette for the barplot
  labs(y = "", fill="Eukaryota\nSubdivision", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+ #x axis header on the top (default=bottom)
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black") +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -4, ymax = -0.5, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black")+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,7.5), ylim = c(-4,NA),clip = "off")
barplot_subdivision_season

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_subdivision_season_barplot.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__month-Sub. ####
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
                                             "Perkinsea","Telonemia",'Ichthyosporea',
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
  geom_col(width = .95, color = "black", linewidth = 0.3)+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 8),
        axis.text.y =element_text(color = "black", size = 10),
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
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+ #x axis header on the top (default=bottom)
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -4, ymax = -0.5, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(-4,NA),clip = "off")
barplot_subdivision

ggsave("/Users/piefouca/Desktop/µEuk/Figures/euk_subdivision_month_barplot.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

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
                year=factor(year,levels = c("2021", "2022"))) %>%
  ggplot(.,aes(Axis1,Axis2,color=season_year,fill=season_year, shape = year, group = season_year))+
  geom_convexhull(alpha=0.6,show.legend = F,size=0.6,)+
  geom_point(show.legend = T,color="black",size=3)+
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10, color="black"),
        axis.ticks = element_line(color="black"),
        plot.title = element_text(size=11,face = "italic"),
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
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  labs(color="Season",
       x=paste0("Axis 1 [",BC_18S_out[[3]],"%]"),
       y=paste0("Axis 2 [",BC_18S_out[[4]],"%]"),
       title=paste0("Micro-eukaryotic Community · Bray-Curtis"))
BC_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/BC_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__Phototrophs ####
# Type trophique ~~~~~~~~~~~~~~~~~~~~
HET <- subset_taxa(rar_water, Order %in% c("Ciliophora",
                                           "Kathablepharida","Fungi","Cercozoa",
                                           "Bigyra","Choanoflagellata", "Chrompodellids"))
AUT <- subset_taxa(rar_water, Order %in% c("Cryptophyta_X", "Chryptophyta_X:nucl", "Gyrista","Chlorophyta_X",
                                           "Streptophyta_X"))
Mixotrophs <- subset_taxa(rar_water, Order %in% c("Dinoflagellata"))

Heterotrophs_parasites <- subset_taxa(rar_water, Order %in% c("Perkinsea","Telonemia_X"))

####__Heterotrophs ####
Heterotrophs_2022 <- subset_taxa(Subset_2022, Order %in% c("Ciliophora",
                                                           "Kathablepharida","Fungi","Cercozoa",
                                                           "Bigyra","Choanoflagellata", "Chrompodellids", "Perkinsea", "Telonemia_X"))
Autotrophs_2022 <- subset_taxa(Subset_2022, Order %in% c("Cryptophyta_X", "Chryptophyta_X:nucl", "Gyrista","Chlorophyta_X",
                                                         "Streptophyta_X"))
Mixotrophs <- subset_taxa(Subset_2022, Order %in% c("Dinoflagellata"))


####__Mixotrophs ####
Heterotrophs_parasites <- subset_taxa(Subset_2022, Order %in% c("Perkinsea","Telonemia_X"))

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
  #points
  geom_point(size=1,alpha=0.4,color="grey",show.legend = F)+
  #linear model
  geom_text(x=250,y=0.21,size=3,check_overlap = T,fontface="italic",label="linear",color="gray32")+
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,formula = lm_formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,formula = lm_formula_linear,
               label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #best poly model poly3
  geom_text(x=250,y=0.132,size=3,check_overlap = T,fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly3,show.legend = F,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly3"),])+
  #best poly model poly4
  geom_text(x=250,y=0.132,size=3,check_overlap = T,fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly4,show.legend = F,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly4"),])+
  #best poly model poly5
  geom_text(x=250,y=0.132,size=3,check_overlap = T,fontface="italic",color = "blue",
            data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly5,show.legend = F,
                 data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly5,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_euk.df[(TLA_euk.df$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=8),
        axis.ticks = element_line(color="black"),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0.1,1))+
  scale_x_continuous(breaks = c(30, 60, 120, 180, 360, 540),expand = c(0,0),limits = c(-1,545))+
  labs(y="Eukaryotic Bray-curtis dissimilarity",
       x="Time between sampling dates (days)")+
  facet_wrap2(~ G1_lakeID,nrow = 3,scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))
TLA_euk

ggsave("/Users/piefouca/Desktop/µEuk/Figures/TLA_euk.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### MOTA ####

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
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$VSS)+0.18)+
  #color CSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$CSM)+0.18)+
  #color LGP
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$LGP)+0.18)+
  #color BLR
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$BLR)+0.18)+
  #color CRE
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$CRE)+0.18)+
  #color `CER-S`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$`CER-S`)+0.18)+
  #color `CER-L`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$`CER-L`)+0.18)+
  #color JAB
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$JAB)+0.18)+
  #color VSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=-0.11,
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
           xmax=max(final_cum_distance_18S$VSM)+0.18)+
  geom_line(size=0.5,color="black",show.legend = F)+
  geom_point(show.legend = T)+
  geom_label(show.legend = F)+
  theme_bw()+
  theme(aspect.ratio=1,
        panel.grid = element_blank(),
        plot.caption =element_text(face="italic",size=10,hjust=1),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_blank(),
        # axis.text.y=element_text(size=12,face="bold",
        #                          color=rev(c("#28A448","#E45000","#980000","#1E7A36","#EB7947","#F2AB8C","#E40000","#2F6B9D","#5C9BD6"))),
        axis.text.y=element_text(size=12,face="bold",color="black"),
        axis.text.x=element_text(size=10, color="black"),
        axis.ticks = element_line(color="black"),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"),
                    labels = c("Summer 2021","Fall 2021","Winter 2021",
                               "Spring 2022","Summer 2022","Fall 2022","Winter 2022"))+
  scale_x_continuous(expand = c(0,0),limits = c(-0.25,10.6))+
  labs(#caption=paste0("MOTA – Bray-Curtis – ",Threshold*100,"% of the total expl. var. (58 axis)"),
    x="micro-Eukaryotic community trajectories",fill="Season")+
  guides(fill=guide_legend(override.aes=list(shape=21,size=8,color="black"),nrow=2))
MOTA_lake_18S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/MOTA_lake_18S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####________________________####
#### Old Amaury ####
####________________________####

# BC-month_gap-normalized --------------------------------------------------

sumBC <- phyloseq::distance(rar_water, method = 'bray') %>%
  as.matrix() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = rar_water@sam_data$month_number[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., month_2 = rar_water@sam_data$month_number[match(.$sample_2, rownames(rar_water@sam_data))]) %>%  
  dplyr::mutate(., Chla1 = rar_water@sam_data$Chla[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Chla2 = rar_water@sam_data$Chla[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., lake_ID = rar_water@sam_data$lake[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., lake_2 = rar_water@sam_data$lake[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., col_1 = rar_water@sam_data$column[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., col_2 = rar_water@sam_data$column[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_1 - month_2))%>%
  dplyr::mutate(., DeltaChla = abs(Chla1-Chla2))
rows_to_keep <- distance_sumBC$col_1 == distance_sumBC$col_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$lake_ID == distance_sumBC$lake_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]

distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., lake_ID=factor(lake_ID, levels=c("VSM", "JAB", "CER-L", 
                                                    "CER-S", "CRE", "BLR",
                                                    "LGP", "CSM", "VSS")))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser / nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month)
# filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  # filter(month_2 != 18)
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

## POLY TYPE COL
distance_normal_sumBC <- distance_normal_sumBC %>%
  mutate(poly = case_when(
    lake_ID %in% 'VSS' ~ 'poly3',
    lake_ID %in% c('BLR', 'CSM', 'CER-S', 'LGP', 'VSM') ~ 'poly4',
    lake_ID %in% c('CTL', 'CER-L', 'JAB') ~ 'poly5'))

v2palette_lake_3T = c('red','red','red','red','red','red','red','red',"red")
lightv2palette_lake_3T = c('violetred','violetred','violetred',
                           'violetred','violetred','violetred',
                           'violetred','violetred',"violetred")

### NORMAL_BC_TLA_plot
tmp <- distance_normal_sumBC %>%
  dplyr::mutate(., lake_ID=factor(lake_ID, levels=c("JAB", "VSM", "CER-L", 
                                              "CSM", "CER-S", "CRE",
                                              "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity, color = lake_ID, fill = lake_ID, label=poly))+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = 'none'
  )+
  geom_smooth(method = 'lm',
              formula =, size=1.2)+
  geom_smooth(inherit.aes = FALSE, method = 'lm',
              mapping = aes(x = normal_month_gap, y = BC_dissimilarity),
              color = 'black', size=1)+
  facet_wrap2(~ lake_ID, scales = "fixed", strip = v2strip_color_lake)+
  # scale_x_continuous(breaks = c(22, 79, 171, 386, 542),
  #                    labels = c('1', '3', '6', '12', '18'))+
  scale_color_manual(values = v2palette_lake_3T)+
  scale_fill_manual(values = lightv2palette_lake_3T)+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  #linear model
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,
                 formula = lm_formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,
               formula = lm_formula_linear,label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #polynomial model
  stat_poly_line(formula = lm_formula_poly4) +
  stat_poly_eq(use_label(c("R2","P")), size=3,
               formula = lm_formula_poly4,
               label.y = 0.01, label.x = 0.97,rr.digits = 2)+
  labs(x='', y='')
  
  
#  labs(title = 'Dissimilarité de Bray Curtis en fonction du nombre de mois séparant les échantillons',
       #x = "Nombre de jours d'écart entre",
       #y = "Dissimilarité de Bray-Curtis")
tmp
View(distance_normal_sumBC)

### NORMAL_BC_TLA_plot
distance_normal_sumBC %>%
  dplyr::mutate(., lake_ID=factor(lake_ID, levels=c("JAB", "VSM", "CER-L", 
                                                    "CER-S", "CRE", "BLR",
                                                    "LGP", "CSM", "VSS")))%>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity, group=lake_ID, color = lake_ID, fill = lake_ID))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks=element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=16),
        axis.text = element_text(size=15, face = 'bold'),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_fill_manual(values = palette_lake_3T)+
  geom_smooth(method = 'loess')+
  # geom_smooth(inherit.aes = FALSE, method = 'lm',
  #             mapping = aes(x = month_gap, y = BC_dissimilarity, group=lake_ID, color = lake_ID), size=.2)+
  scale_color_manual(values = palette_lake_3T)+
  scale_x_continuous(limit=c(0,17.5), expand=c(0,0), breaks = c(1,3,6,9,12,16), labels = c(1,3,6,9,12,16))+
  labs(title = 'Time Lag Analysis', fill = 'lac', colour = 'lac',
       x = "Nombre de mois d'écart",
       y = "Dissimilarité de Bray-Curtis")+
  coord_cartesian(ylim = c(0.5, 1), xlim = c(0.5, 16.5))
# caption = "L'axe y débute à 0.5 car aucune valeur inférieure n'a été obtenue"

###
# BC-T°C_gap-normalized --------------------------------------------------
sumBC <- phyloseq::distance(rar_water, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., Temp_1 = rar_water@sam_data$Temperature[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Temp_2 = rar_water@sam_data$Temperature[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., lake_1 = rar_water@sam_data$lake[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., lake_2 = rar_water@sam_data$lake[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Temp_compare = paste0(Temp_1, '_', Temp_2))
sumBC$Temp_1 <- as.numeric(as.character(sumBC$Temp_1))
sumBC$Temp_2 <- as.numeric(as.character(sumBC$Temp_2))
View(distance_sumBC)
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,Temp_gap = abs(Temp_2 - Temp_1))
rows_to_keep <- distance_sumBC$lake_1 == distance_sumBC$lake_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$Temp_2 > distance_sumBC$Temp_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]
distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., lake_1=factor(lake_1, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))


### 
distance_sumBC %>% na.omit() %>%
  ggplot(., aes(x = Temp_gap, y = BC_dissimilarity))+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  geom_smooth(inherit.aes = FALSE, method = 'loess', 
              mapping = aes(x = Temp_gap, y = BC_dissimilarity), 
              color = 'black', size=1.2)+
  facet_wrap2(~ lake_1, scales = "fixed", strip = strip_color_lake)+
  # scale_x_continuous(breaks = c(1, 65, 157, 276, 372, 521),
  #                    labels = c('1', '3', '6', '9', '12', '17'))+
  labs(title = 'Dissimilarité de Bray Curtis en fonction du degré de température que séparent les deux échantillons',
       x = "Nombre de degré Celcius d'écart",
       y = "Dissimilarité de Bray-Curtis")




distance_sumBC %>% na.omit() %>%
  dplyr::mutate(., lake_2=factor(lake_2, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))%>%
  ggplot(., aes(x = Temp_gap, y = BC_dissimilarity, group=lake_1, grp.label = lake_1))+
  geom_jitter(size=0, color='white')+
  theme_bw()+
  # stat_poly_line() +
  # stat_poly_eq(aes(label = paste(after_stat(grp.label), "*\": \"*",
  #                                after_stat(eq.label), "*\", \"*",
  #                                after_stat(rr.label), sep = ""))) +
  theme(
    panel.grid = element_blank()
  )+
  geom_smooth(inherit.aes = FALSE, method = 'lm', 
              mapping = aes(x = Temp_gap, y = BC_dissimilarity, group=lake_2, color = lake_2), size=1.2)+
  scale_color_manual(values = palette_lake_3T)+
  # scale_x_continuous(limit=c(0,17.5), expand=c(0,0), breaks = c(1,3,6,9,12,16), labels = c(1,3,6,9,12,16))+
  labs(colour = 'Lac', title = 'Dissimilarité de Bray Curtis en fonction de la température qui sépare les deux échantillons',
       x = "Nombre de degré Celcius d'écart",
       y = "Dissimilarité de Bray-Curtis",
       caption = "L'axe y débute à 0.4 car aucune valeur inférieure n'a été obtenue")+
  coord_cartesian(ylim = c(0.4, 1), xlim = c(-0.49, 24.3))

# Chla2by2 ----------------------------------------------------------------
Fig_month_BC <- final_sumBC %>%
  ggplot(., aes(x= month_1, y = ChlaComp, group= month_2))+
  # geom_boxplot()+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank()
  )+
  facet_wrap2(~ lake_1, scales = "free_y", strip = strip_color_lake)+
  geom_smooth(inherit.aes = FALSE, method = loess, 
              mapping = aes(x = month_1, y = BC_dissimilarity, group=lake_1), 
              color = 'black', size=1.2)+
  labs(title = 'Chlaf(lacs/Months)')
Fig_month_BC
# Chla --------------------------------------------------------------------
sam_dat <- as.data.frame(rar_water@sam_data)
Plottmp <- ggplot(sam_dat, aes(x=Chla, y=1, group='lake', color = lake))+
  geom_point(size = 4)+
  theme_bw()+
  ylim(1-0.05,1+0.05)+
  theme(
    axis.text.y = element_blank()
  )+
  scale_color_manual(values = rev(c("#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D"))) +
  scale_x_log10()+
  annotate("segment", x = 2.6, xend = 2.6, y = 0.95, yend = 1.05,
           colour = "#2F6B9D", size=1)+
  annotate("segment", x = 7.4, xend = 7.2, y = 0.95, yend = 1.05,
           colour = "#EB7947", size=1)+
  annotate("segment", x = 56, xend = 55, y = 0.95, yend = 1.05,
           colour = "#1E7A36", size=1)
Plottmp

#                               

# Chla concentrations -----------------------------------------------------------------
transform(rar_water@sam_data, Chla = as.numeric(Chla))
transform(rar_water@sam_data, Year = as.character(Year))

palette_Chla = c('lightsteelblue1','lightsteelblue4')
palette_Chla = c('azure1','azure3')
palette_Chla = c('#F4FFF4','#BFD6BF')

## -> Full years
CHLA_Year <- as.data.frame(as.tibble(rar_water@sam_data)) %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VSM", "CER-L", 
                                              "CSM", "CER-S", "CRE",
                                              "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(x=lake, y=Chla, fill=as.factor(Year))) + 
  geom_boxplot()+
  geom_hline(yintercept=2.6, linetype="dashed", color = "#B0D0E8", size=1.2)+
  annotate("text", x =6, y = 1.5, size = 5,
           label = "paste(bold(Oligotrophe))", parse = TRUE)+
  geom_hline(yintercept=7.2, linetype="dashed", color = "#C99594", size=1.2)+
  annotate("text", x =7, y = 3.2, size = 5,
           label = "paste(bold(Mésotrophe))", parse = TRUE)+
  geom_hline(yintercept=55, linetype="dashed", color = "#94C999", size=1.2)+
  annotate("text", x =1.5, y = 10, size = 5,
           label = "paste(bold(Eutrophe))", parse = TRUE)+
  annotate("text", x =1.5, y = 100, size = 5,
           label = "paste(bold(Hypereutrophe))", parse = TRUE)+
  # geom_jitter(size=.8, color='azure4')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.title = element_text(size=18),
        axis.text.x = element_text(size=15, vjust =.8, colour = palette_lake_3T),
        legend.title = element_text(size=17),
        aspect.ratio = 1,
        legend.text = element_text(size=14))+
  scale_fill_manual(values = palette_Chla)+
  theme(axis.text = element_text(size = 10, face = 4))+
  labs(y = 'Chlorophylle a (µg/L log10)', fill = 'Année', x = 'Lac')+
  scale_y_log10()
CHLA_Year

palette_lake_3T<-rev(c("purple","grey10","grey20",
                       "grey30","grey40","hotpink",
                       "grey50","grey60","grey70"))

## -> Shared months
CHLA_shared <- as.data.frame(as.tibble(ps_Shared@sam_data)) %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VSM", "CER-L", 
                                              "CSM", "CER-S", "CRE",
                                              "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(x=lake, y=Chla, fill=as.factor(Year))) + 
  geom_boxplot()+
  geom_hline(yintercept=2.6, linetype="dashed", color = "#B0D0E8", size=1.2)+
  annotate("text", x =6, y = 1,
           label = "paste(bold(Oligotroph))", parse = TRUE)+
  geom_hline(yintercept=7.2, linetype="dashed", color = "#9CCABA", size=1.2)+
  annotate("text", x =7, y = 3.4,
           label = "paste(bold(Mesotroph))", parse = TRUE)+
  geom_hline(yintercept=55, linetype="dashed", color = "#94C999", size=1.2)+
  annotate("text", x =1.5, y = 10,
           label = "paste(bold(Eutroph))", parse = TRUE)+
  annotate("text", x =1.5, y = 100,
           label = "paste(bold(Hypereutroph))", parse = TRUE)+
  geom_jitter(size=.8, color='ivory3')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=15, vjust =.8, colour = palette_lake_3T),
        legend.title = element_text(size=15),
        aspect.ratio = 1,
        legend.text = element_text(size=14))+
  scale_fill_manual(values = palette_Chla)+
  theme(axis.text = element_text(size = 10, face = 4))+
  labs(subtitle = 'Shared months between 2021 and 20222 in the dataset (June to November)', y = 'Chlorophyll a (log10)', fill = 'Year', x = 'Lake')+
  scale_y_log10()
CHLA_shared

## -> 4 summer monts
CHLA_Summer_shared <- as.data.frame(as.tibble(ps_4Summer@sam_data)) %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI",
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  ggplot(., aes(x=lake, y=Chla, fill=as.factor(Year))) + 
  geom_boxplot()+
  geom_hline(yintercept=2.6, linetype="dashed", color = "#B0D0E8", size=1.2)+
  annotate("text", x =6, y = 1,
           label = "paste(bold(Oligotrophe))", parse = TRUE)+
  geom_hline(yintercept=7.2, linetype="dashed", color = "#C99594", size=1.2)+
  annotate("text", x =6, y = 3.2,
           label = "paste(bold(Mésotrophe))", parse = TRUE)+
  geom_hline(yintercept=55, linetype="dashed", color = "#94C999", size=1.2)+
  annotate("text", x =1.5, y = 10,
           label = "paste(bold(Eutrophe))", parse = TRUE)+
  annotate("text", x =1.5, y = 100,
           label = "paste(bold(Hypereutrophe))", parse = TRUE)+
  geom_jitter(size=.8, color='ivory4')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=15, vjust =.8, colour = palette_lake_3T),
        legend.title = element_text(size=15),
        aspect.ratio = 1,
        legend.text = element_text(size=14))+
  scale_fill_manual(values = palette_Chla)+
  theme(axis.text = element_text(size = 10, face = 4))+
  labs(subtitle = "Mois d'été partagés (Juin à Août)", y = 'Chlorophylle a (log10)', fill = 'Année', x = 'Lac')+
  scale_y_log10()

CHLA_Year
CHLA_shared
CHLA_Summer_shared

##
tmp <- as.data.frame(rar_water@sam_data) %>%
  # group_by(lake, month_number) %>%
  # summarise_at(c("Chla"), mean) %>%
  ggplot(., aes(x = month_number, y = Chla)) +
  geom_point(size = .8)+
  theme_bw()+
  geom_hline(yintercept=2.6, linetype="dashed", color = "royalblue")+
  geom_hline(yintercept=7.4, linetype="dashed", color = "salmon")+
  geom_hline(yintercept=56, linetype="dashed", color = "seagreen")+
  facet_wrap(~lake, scale='free')+
  # scale_y_log10()+
  geom_smooth(method = loess)
tmp

tmp <- as.data.frame(rar_water@sam_data) %>%
  # dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI",
  #                                             "CRJ", "CTL", "BOI",
  #                                             "GDP", "CHA", "VER"))) %>%
  ggplot(., aes(x = month, y = Chla)) +
  geom_boxplot() +
  geom_jitter(size=.5, coor = 'grey')+
  theme_bw()+
  geom_hline(yintercept=2.6, linetype="dashed", color = "royalblue", size=.8)+
  geom_hline(yintercept=7.4, linetype="dashed", color = "salmon",size =.8)+
  geom_hline(yintercept=56, linetype="dashed", color = "seagreen", size=.8)+
  # scale_y_log10()+
  facet_wrap2(~ lake, scales = "free", strip = strip_color_lake)
tmp

View(tmp$data)

ggplot(sam_dat, aes(x=month_number, y=Chla))+
  geom_point()+
  geom_smooth(method=loess)+
  scale_y_log10()

tmp <- as.data.frame(rar_water@sam_data)
med_chla <- tmp %>%
  group_by(lake) %>%
  summarize(median_chla = median(Chla, na.rm = TRUE))
View(med_chla)
# MOTA - Chla variations ---------------------------------------------------------
sam_dat <- rar_water@sam_data
transform(sam_dat, Chla = as.numeric(Chla))

CHLA_raw <- sam_dat %>%
  group_by(lake, month_number) %>%
  summarise(mean_Chla = median(Chla))

# subset.data.frame(., lake %in% 'BOI')
CHLA_1 <- CHLA_raw  %>%
  mutate(diff = mean_Chla - lag(mean_Chla, default = first(mean_Chla)))
CHLA_2 <- na.omit(CHLA_1) 
create_new_column <- function(column) {
  new_column <- numeric(length(column))
  new_column[1] <- 0  # Mettre 0 comme première valeur
  for (i in 2:length(column)) {
    new_column[i] <- column[i - 1] + new_column[i - 1]
  }
  return(new_column)
}
coord <- create_new_column(CHLA_2$diff)
View(coord)
CHLA_3 <- CHLA_2
CHLA_3$diff <- abs(CHLA_3$diff)



DT <- data.table(CHLA_3, key = "lake")
DT[, csum := cumsum(diff), by = key(DT)]

View(DT)

CHLA_df <- melt(DT, value.name = DT$lake)

DT_1 <- reshape(DT, idvar = "month_number", timevar = "lake", direction = "wide")
View(DT_1)



palette_season <- c('#B7D7FF', '#FFFED5', '#FFD9C8', '#E9D2A4' )

MOTA_CHLA <- DT %>% as.data.frame() %>%
  # filter(!lake=='VER') %>%
  dplyr::mutate(lake_month=paste0(lake,"_",month_number),
                season=if_else(month_number %in% c("7","8", "18"),"hiver",
                               if_else(month_number %in% c("15","16","17","4","5","6"),"automne",
                                       if_else(month_number %in% c("9","10","11"),"printemps","été"))),
                season=factor(season,levels = c("hiver","printemps","été","automne"))) %>%
  group_by(lake) %>%
  dplyr::mutate(.,max_dist=max(csum)) %>%
  ggplot(.,
         aes(x=csum,y=(fct_reorder(lake,max_dist)),label=month_number,fill=season))+
  geom_line(size=0.5,color="black",show.legend = F )+
  geom_line(size=0.5,color="black",show.legend = F)+
  geom_point(show.legend = T)+
  geom_label(show.legend = F)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.caption =element_text(face="italic",size=10,hjust=1),
    axis.title.x=element_text(size=12,face="bold"),axis.title.y=element_blank(),
    axis.text.y=element_text(size=14,face="bold",
                             color=rev(c("#52C655","#52C655","#473A3B","#52C655","#52C655","#52AAC6","#473A3B","#52AAC6","#52AAC6"))),
    axis.text.x=element_text(size=10, face = 'bold'),
    axis.ticks = element_blank(),
    legend.title = element_text(size=12,face = "bold",hjust=0),
    legend.text = element_text(size=11),
    legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=palette_season)+
  labs(caption=paste0("MOTA – CHLA"),
       x="~ Distance parcourue, cumsum de médiane de Chla",fill="Saison")+
  guides(fill=guide_legend(override.aes=aes(shape=22,size=8),nrow=2))


MOTA_CHLA_noVER
MOTA_CHLA

Fig_MOTA_CHLA <- MOTA_CHLA + MOTA_CHLA_noVER  
Fig_MOTA_CHLA


# Params PCA --------------------------------------------------------------

# Table params numerical ~~~~~~~~~~~~
numerical <- rar_water@sam_data %>% subset(., select=c(Temperature, pH, NO3NO2, Salinity, PO4, NH4, Chla, TPC, TPN))
numerical$Temperature <- as.numeric(numerical$Temperature)
numerical$pH <- as.numeric(numerical$pH)
numerical$NH4 <- as.numeric(numerical$NH4)
numerical$PO4 <- as.numeric(numerical$PO4)
numerical$NO3NO2 <- as.numeric(numerical$NO3NO2)
numerical$TPN <- as.numeric(numerical$TPN)
numerical$TPC <- as.numeric(numerical$TPC)
numerical <- na.omit(numerical)
numerical <- numerical %>%
  scale(., center =T, scale =T) 


# Fig - corrélations params ~~~~~~~~~
pca_res <- prcomp(numerical)
var <- get_pca_var(pca_res)
var_dim1<-var$cos2 %>% .[,c(1,2,3,4)]
melted_corr <- melt(var_dim1)
corr_Fig<- melted_corr %>%  
  # subset(.,Var2=="Dim.1") %>%
  ggplot(.,aes(x = Var2, y = Var1, fill = value, size=value))+
  geom_point(color='black',shape=21) + theme_bw()+
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)), color = "black",size = 3.7)+
  #labs(fill = expression(Cos^"2"),size="",x="  PC1      -      PC2      -        PC3        -        PC4")+
  theme(panel.grid=element_blank(),
        axis.text.y = element_text(face="bold",size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",size = 10))+
  scale_fill_gradient2(midpoint = 0.6, low ="#F4FFF4", mid = "palegreen1", high = "skyblue4",
                       limits = c(0, 1)) +
  scale_size_continuous(range = c(9,18),breaks = c(0.2,0.4,0.8),labels = c("0.2","0.4","0.6"))
corr_Fig

# palette_Chla = c('#F4FFF4','#BFD6BF')
# mid ="hotpink3", low = "palegreen3", high = "grey30",

fviz_pca_biplot(pca_res, 
                # Individuals
                geom.ind = "point",
                fill.ind = PCdf@sam_data$lake, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Lake", color = "Contrib",
                                    alpha = "Contrib")
)

fviz_pca(pca_res)

# Cercle des corrélations ~~~~~~~~~~~
fviz_pca_var(comp)
fviz_pca_biplot(comp)
fviz_pca_var(comp, col.var="contrib")+
  scale_color_gradient2(low="tomato", mid="red",
                        high="black", midpoint=14, space ="Lab") +
  theme_minimal()

corrplot(BC$vectors, is.corr = FALSE)

# PCoA BC -----------------------------------------------------------------

# Ordination + lien avec params ~~~~~
BC<-ordinate(PCdf,method = 'PCoA' ,distance = "bray")
param_env <- envfit(BC$vectors, env=numerical, perm = 999)
param_env.df <- as.data.frame(scores(param_env,"vectors"))
param_env.df <- cbind(param_env.df, param = rownames(param_env.df))
coord.PCoA <- as.data.frame(BC$vectors)
coord.PCoA$lake = PCdf@sam_data$lake
coord.PCoA$Cal_month <- PCdf@sam_data$Cal_month


# Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig_mds_param <- coord.PCoA %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI", 
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  ggplot(.,aes(Axis.1,-Axis.2,color = lake, label=Cal_month)) +
  geom_convexhull(aes(color = lake,fill=factor(lake),group = lake),alpha=0.2,show.legend = F,size=0.4)+
  geom_segment(data = param_env.df,
               mapping=aes(x = 0, xend = Axis.1*1, y = 0, yend = -Axis.2*1),
               color = "black",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  geom_text(data = param_env.df, mapping=aes(x = Axis.1*1.1, y = -Axis.2*1.1, label = param),
            size = 3.5,fill="white",fontface = "bold",inherit.aes = F)+
  geom_text(size=6, key_glyph="point")+
  # geom_point(size=2,shape=16,show.legend = F) + 
  theme_bw() +
  theme(plot.title=element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.caption = element_text(size=12,face = "italic",color="black"),
        panel.grid = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=palette_lake_3T)+
  scale_fill_manual(values=palette_lake_3T)+
  #scale_x_reverse()+
  scale_y_reverse()+
  #scale_shape_manual(values = Summer_2021)+
  guides(color = guide_legend(override.aes=list(size=6)),
         fill = guide_legend(override.aes=list(size=6)))+
  labs(color="Lake",shape="Month",
       x=paste0("PC1 [",round(BC$values$Eigenvalues[1]*100/sum(BC$values$Eigenvalues),1),"%]"),
       y=paste0("PC2 [",round(BC$values$Eigenvalues[2]*100/sum(BC$values$Eigenvalues),1),"%]"))
# title=expression(paste("TPN, TPC, ",NH[4]^"-",", ",NO[3]^"-"," & ",NO[2]^"-",", ",PO[4]^"3-"))
# title="Enviromental parameters and lake morphometry PCA"
Fig_mds_param

# Bonne représentation ? ~~~~~~~~~~~~
comp <- prcomp(numerical, scale = TRUE)
comp$rotation <- -1*comp$rotation
comp$x <- -1*comp$x
fviz_pca_var(comp)
fviz_pca_biplot(comp)
fviz_pca_var(comp, col.var="contrib")+
  scale_color_gradient2(low="white", mid="red",
                        high="black", midpoint=6, space ="Lab") +
  theme_minimal()
Fig_mds_param

View(rar_water@sam_data)
# PC f(param) --------------------------------------------------------------
BC_SAMDAT <- merge(SamPCdf, coord.PCoA)
View(BC_SAMDAT)
BC_SAMDAT %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI", 
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  ggplot(., aes(Chla, Axis.1, color=lake, group=lake))+
  geom_point(size=3)+
  scale_color_manual(values=palette_lake_3T)
# + scale_x_log10()

BC_SAMDAT$Temperature <- as.character(BC_SAMDAT$Temperature)
BC_SAMDAT$pH <- as.character(BC_SAMDAT$pH)
BC_SAMDAT$Temperature <- as.numeric(BC_SAMDAT$Temperature)
BC_SAMDAT$pH <- as.numeric(BC_SAMDAT$pH)

BC_SAMDAT %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI",
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  ggplot(., aes(Axis.1, pH, color=lake))+
  geom_convexhull(aes(color = lake,fill=factor(lake),group = lake),alpha=0.2,show.legend = F,size=0.4)+
  geom_point(size=3)+
  scale_color_manual(values=palette_lake_3T)
# Winter_Summer - alphadiv ---------------------------------------------------------
# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
WIN_alpha <- WIN@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(WIN,measures = "Observed")$Observed,
                Shannon=estimate_richness(WIN,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
View(WIN_alpha)

# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_WIN_alpha <- WIN_alpha %>%
  group_by(lake) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
View(mean_WIN_alpha)

mean(mean_WIN_alpha$mean_Richness)
mean(mean_WIN_alpha$sd_Richness)

# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
SUM_alpha <- SUM@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(SUM,measures = "Observed")$Observed,
                Shannon=estimate_richness(SUM,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
View(SUM_alpha)

# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_SUM_alpha <- SUM_alpha %>%
  group_by(lake) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
View(mean_SUM_alpha)

mean(mean_SUM_alpha$mean_Richness)
mean(mean_SUM_alpha$sd_Richness)

# T-Student
lacs <- intersect(WIN_alpha$lake, SUM_alpha$lake)
results <- list()

# Boucle sur les lacs communs pour effectuer le test de Student
for (lac in lacs) {
  # Extraire les données correspondantes aux deux dataframes pour ce lac
  data_lac_df1 <- subset(WIN_alpha, lake == lac)$Richness
  data_lac_df2 <- subset(SUM_alpha, lake == lac)$Richness
  
  # Effectuer le test de Student
  test_result <- t.test(data_lac_df1, data_lac_df2)
  
  # Stocker les résultats dans la liste
  results[[lac]] <- test_result
}

for (lac in names(results)) {
  cat("Test de Student pour le lac", lac, ":\n")
  print(results[[lac]])
  cat("\n")
}
# HET_AUT - alphadiv ---------------------------------------------------------
# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
HET_alpha <- HET@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(HET,measures = "Observed")$Observed,
                Shannon=estimate_richness(HET,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
View(HET_alpha)

# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_HET_alpha <- HET_alpha %>%
  group_by(lake) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
View(mean_HET_alpha)

# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
AUT_alpha <- AUT@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(AUT,measures = "Observed")$Observed,
                Shannon=estimate_richness(AUT,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
View(AUT_alpha)

# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_AUT_alpha <- AUT_alpha %>%
  group_by(lake) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
View(mean_AUT_alpha)

# T-Student
lacs <- intersect(HET_alpha$lake, AUT_alpha$lake)
results <- list()

# Boucle sur les lacs communs pour effectuer le test de Student
for (lac in lacs) {
  # Extraire les données correspondantes aux deux dataframes pour ce lac
  data_lac_df1 <- subset(HET_alpha, lake == lac)$Richness
  data_lac_df2 <- subset(AUT_alpha, lake == lac)$Richness
  
  # Effectuer le test de Student
  test_result <- t.test(data_lac_df1, data_lac_df2)
  
  # Stocker les résultats dans la liste
  results[[lac]] <- test_result
}

for (lac in names(results)) {
  cat("Test de Student pour le lac", lac, ":\n")
  print(results[[lac]])
  cat("\n")
}

View(rar_water@otu_table)



sumBC <- phyloseq::distance(HET, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = HET@sam_data$month_number[match(.$sample_1, rownames(HET@sam_data))]) %>%
  dplyr::mutate(., month_2 = HET@sam_data$month_number[match(.$sample_2, rownames(HET@sam_data))]) %>%
  dplyr::mutate(., lake_1 = HET@sam_data$lake[match(.$sample_1, rownames(HET@sam_data))]) %>%
  dplyr::mutate(., lake_2 = HET@sam_data$lake[match(.$sample_2, rownames(HET@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_2 - month_1))
rows_to_keep <- distance_sumBC$lake_1 == distance_sumBC$lake_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]
distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., lake_1=factor(lake_1, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser par le nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month) %>%
  filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  filter(month_2 != 18) %>%
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

### NORMALIZED
distance_normal_sumBC %>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity))+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  geom_smooth(inherit.aes = FALSE, method = 'loess', 
              mapping = aes(x = normal_month_gap, y = BC_dissimilarity), 
              color = 'black', size=1.2)+
  facet_wrap2(~ lake_1, scales = "fixed", strip = strip_color_lake)+
  scale_x_continuous(breaks = c(1, 65, 157, 276, 372, 521),
                     labels = c('1', '3', '6', '9', '12', '17'))+
  labs(title = 'Dissimilarité de Bray Curtis en fonction du nombre de mois séparant les échantillons',
       x = "Nombre de mois d'écart (normalisé par le nombre de jours entre les échantillonnages)",
       y = "Dissimilarité de Bray-Curtis")



tmpHET <- distance_normal_sumBC %>%
  dplyr::mutate(., lake_2=factor(lake_2, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))%>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity, group=lake_1))+
  geom_jitter(size=3, color='black')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks=element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=16),
        axis.text = element_text(size=15, face = 'bold'),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_lake_3T)+
  scale_fill_manual(values = palette_lake_3T)+
  guides(color=guide_legend(override.aes = list(shape=16)))+
  geom_smooth(inherit.aes = FALSE, method = NULL, formula = y ~ poly(x,3), 
              mapping = aes(x = month_gap, y = BC_dissimilarity, group=lake_2, color = lake_2), size=1.2)+
  geom_smooth(inherit.aes = FALSE, method = 'lm',
              mapping = aes(x = month_gap, y = BC_dissimilarity, group=lake_2, color = lake_2), size=.2)+
  scale_color_manual(values = palette_lake_3T)+
  scale_x_continuous(limit=c(0,17.5), expand=c(0,0), breaks = c(1,3,6,9,12,16), labels = c(1,3,6,9,12,16))+
  labs(title = 'Heterotrophs TLA', colour = '      Lac',
       x = "Nombre de mois d'écart (normalisé par le nombre de jours entre les échantillonnages)",
       y = "Dissimilarité de Bray-Curtis", caption = "L'axe y débute à 0.5 car aucune valeur inférieure n'a été obtenue")+
  coord_cartesian(ylim = c(0.5, 1), xlim = c(0.5, 16.5))
tmpHET



sumBC <- phyloseq::distance(AUT, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = AUT@sam_data$month_number[match(.$sample_1, rownames(AUT@sam_data))]) %>%
  dplyr::mutate(., month_2 = AUT@sam_data$month_number[match(.$sample_2, rownames(AUT@sam_data))]) %>%
  dplyr::mutate(., lake_1 = AUT@sam_data$lake[match(.$sample_1, rownames(AUT@sam_data))]) %>%
  dplyr::mutate(., lake_2 = AUT@sam_data$lake[match(.$sample_2, rownames(AUT@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_2 - month_1))
rows_to_keep <- distance_sumBC$lake_1 == distance_sumBC$lake_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]
distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., lake_1=factor(lake_1, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser par le nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month) %>%
  filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  filter(month_2 != 18) %>%
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

### NORMALIZED
distance_normal_sumBC %>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity))+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  geom_smooth(inherit.aes = FALSE, method = 'loess', 
              mapping = aes(x = normal_month_gap, y = BC_dissimilarity), 
              color = 'black', size=1.2)+
  facet_wrap2(~ lake_1, scales = "fixed", strip = strip_color_lake)+
  scale_x_continuous(breaks = c(1, 65, 157, 276, 372, 521),
                     labels = c('1', '3', '6', '9', '12', '17'))+
  labs(title = 'Dissimilarité de Bray Curtis en fonction du nombre de mois séparant les échantillons',
       x = "Nombre de mois d'écart (normalisé par le nombre de jours entre les échantillonnages)",
       y = "Dissimilarité de Bray-Curtis")



tmpAUT <- distance_normal_sumBC %>%
  dplyr::mutate(., lake_2=factor(lake_2, levels=c("JAB", "VAI", "TRI", 
                                                  "CRJ", "CTL", "BOI",
                                                  "GDP", "CHA", "VER")))%>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity, group=lake_1))+
  geom_jitter(size=3, color='black')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks=element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=16),
        axis.text = element_text(size=15, face = 'bold'),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_lake_3T)+
  scale_fill_manual(values = palette_lake_3T)+
  guides(color=guide_legend(override.aes = list(shape=16)))+
  geom_smooth(inherit.aes = FALSE, method = NULL, formula = y ~ poly(x,3), 
              mapping = aes(x = month_gap, y = BC_dissimilarity, group=lake_2, color = lake_2), size=1.2)+
  geom_smooth(inherit.aes = FALSE, method = 'lm',
              mapping = aes(x = month_gap, y = BC_dissimilarity, group=lake_2, color = lake_2), size=.2)+
  scale_color_manual(values = palette_lake_3T)+
  scale_x_continuous(limit=c(0,17.5), expand=c(0,0), breaks = c(1,3,6,9,12,16), labels = c(1,3,6,9,12,16))+
  labs(title = 'Autotrophs TLA', colour = '      Lac',
       x = "Nombre de mois d'écart (normalisé par le nombre de jours entre les échantillonnages)",
       y = "Dissimilarité de Bray-Curtis", caption = "L'axe y débute à 0.5 car aucune valeur inférieure n'a été obtenue")+
  coord_cartesian(ylim = c(0.5, 1), xlim = c(0.5, 16.5))
tmpAUT

Fig <- ggarrange(tmp, tmpHET, tmpAUT, nrow = 1)
Fig


#### SIMPER ####

####__ 0.1% Filter ####

total_depth <- sum(taxa_sums(rar_water))
threshold <- 0.001 * total_depth #0.1% abundance
abundant.taxa <-(taxa_sums(rar_water) > threshold)
View(abundant.taxa)

# tax table simper
tax_simper<- rar_water@tax_table %>% t(.) %>% as.data.frame(.) %>% t(.)
tax_simper <- tax_simper %>% as.data.frame()

ps_summer_01_total<-prune_taxa(abundant.taxa,rar_water)
View(abundant.taxa)
# otu table simper
simper_01_total_com.df <- ps_summer_01_total@otu_table %>% as.data.frame() %>% t()
View(simper_01_total_com.df)
simper_total.df<-simper(simper_01_total_com.df,rar_water@sam_data$Group,permutations = 999)
simp<- summary(simper_total.df, ordered =T)[[1]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
View(simper_total.df)
# simpBOI<- summary(simper_BOI_vs_others.df, ordered =T)[[1]] %>% as.data.frame() %>%
#   subset(.,cumsum <= 0.7 & p <=0.001)
# list_simper_BOI<-rownames(simp)
View(simp %>% dplyr::mutate(.,Phylum=tax_simper$Phylum[match(rownames(.),rownames(tax_simper))]))
simp_Order <- simp %>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))])
# # BOI vs others
# simper_BOI_vs_others.df<-simper(simper_01_total_com.df,rar_water@sam_data$Boi,permutations = 999)
# View(simpBOI %>% dplyr::mutate(.,Phylum=tax_simper$Phylum[match(rownames(.),rownames(tax_simper))]))
# View(simpBOI %>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))
# View(simp)

### Trophic Status ###
simper_01_total.df<-simper(simper_01_total_com.df,A_D_rarW_freq@sam_data$lake_type,permutations = 999,parallel = 5)

###__ Contrast H _ M ###
simper_H_M <- summary(simper_01_total.df)[[3]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)

list_simper_H_M<-rownames(simper_H_M)

View(ps_summer_01_total@tax_table)
View(rar_water@otu_table)

ASV <- as_data_frame(rar_water@otu_table)
simptest <- simper(ASV, rar_water@sam_data$lake, permutations = 10)
View(rar_water@tax_table)

####

View(simp_Order)
simp_Order$ASVlist <- rownames(simp_Order)
simp_Order$average <- 100*simp_Order$average

ggplot(simp_Order, aes(x = ASVlist, y = average)) +
  geom_bar(stat = "identity", fill = "turquoise4") +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2) +
  labs(x = "ASV List", y = "Average", title = "Barplot with Error Bars") +
  theme_minimal()








# Contrib SIMPER ----------------------------------------------------------
df_SIMP <- data.frame(Taxa = c("Cryptophyta","Ciliophora","Chlorophyta","Gyrista","Kathablepharida"), Contrib = c('24.07','5.69','4.75','1.54','1.32'))
df_SIMP$Contrib <- as.numeric(df_SIMP$Contrib)
View(df_SIMP)
df_SIMP %>% 
  dplyr::mutate(., Taxa=factor(Taxa, levels=c("Cryptophyta","Ciliophora",
                                              "Chlorophyta","Gyrista","Kathablepharida"))) %>%
  ggplot(., aes(x=Taxa, y=Contrib, fill=Taxa))+
  geom_col(size = 2)+
  scale_fill_manual(values=SIMPalette)+
  labs(x = 'Taxa', y= 'Contribution', title = 'Part de la contribution à la variance entre tous les lacs pour les 10 ASV les plus explicatifs des différences de communautés')

SIMPalette = c('#ffd966','#f4b084','#a9d08e','#548235','#c65911')


num <- as.data.frame(numerical)

num <- as.data.frame(rar_water@sam_data)
num$pH <- as.numeric(num$pH) 
ggplot(num, aes(x=pH, y=Chla, fill=lake))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_fill_manual(values=c('grey80','grey70','grey60','grey50','grey40','grey30','grey20','grey10','black'))




# Test temature -----------------------------------------------------------
temperature_Summer_21vs22 <- subset_samples(rar_water, season %in% c('Summer_2021', 'Summer_2022')) 
temperature_Summer_21vs22 <- temperature_Summer_21vs22@sam_data %>% as.data.frame()
temperature_Summer_21vs22$Temperature <- as.numeric(temperature_Summer_21vs22$Temperature)

plot_temperature_summer_21vs22 <- temperature_Summer_21vs22 %>%
  ggplot(., aes(x=season, y=Temperature, fill=season))+
  geom_boxplot()+
  geom_signif(comparisons = list(c('Summer_2021', 'Summer_2022')), map_signif_level = TRUE )+
  theme_bw()+
  scale_fill_manual(values = c('#f5df4d', '#6667ab'))+
  labs(fill='Saison', y='Température (°C)', x='Saison', caption='p value < 0.001')

plot_temperature_summer_21vs22

#####
#####~Groups~#####

# BC ----------------------------------------------------------------------

BC<-ordinate(rar_water,method = "PCoA",distance = "bray")


## Barycentre ##
vectors <- as.data.frame(BC$vectors)
sam <- rar_water@sam_data
Coords <- data.frame(Axis.1 = vectors$Axis.1, Axis.2 = vectors$Axis.2, Group = sam$Group, month_number = sam$month_number)
BCarycentre <- Coords %>%
  group_by(Group) %>%
  dplyr::mutate(bary1 = median(Axis.1), bary2 = median(Axis.2))
#~#~#~#~#~#~#~#~#


Fig_rar_BC.df<-plot_ordination(rar_water,BC,color="Group", label="month_number")

Fig_rar_Group_BC <- as.data.frame(Fig_rar_BC.df$data) %>%
  dplyr::mutate(., bary1 = BCarycentre$bary1, bary2 = BCarycentre$bary2) %>%
  group_by(Group, month, month_number) %>%
  # summarise_at(c("Axis.1", "Axis.2"), mean) %>%
  dplyr::mutate(., Group=factor(Group, levels=c('Euplus','Mesomoins','Shift'))) %>%
  ggplot(., aes(Axis.1, Axis.2, color=Group, label=Group))+
  geom_convexhull(aes(fill=Group, group=Group), size=0.9, alpha=0.04, show.legend = FALSE)+
  geom_point(size=1)+
  geom_label(x=BCarycentre$bary1, y=BCarycentre$bary2, label=BCarycentre$Group, size = 5, show.legend = FALSE)+
  # geom_text(size=5, key_glyph="point")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks=element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=17),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_Group)+
  scale_fill_manual(values = palette_Group)+
  guides(color=guide_legend(override.aes = list(shape=16)))+
  labs(title="Bray-Curtis ~ PCoA",color="       Groupe",
       x=paste("Axis 1 [", round(BC$values$Relative_eig[1]*100, 2), "] %"),
       y=paste("Axis 2 [", round(BC$values$Relative_eig[2]*100, 2), "] %"))
# caption = "Ordination sur l'ensemble du jeu de données (9 lacs 18 mois), un polygone et barycentre par lac")
Fig_rar_Group_BC


###
FigBCs <- ggarrange( Fig_rar_BC,Fig_rar_Group_BC, nrow = 1)
FigBCs
###





# BC seasons --------------------------------------------------------------
BC<-ordinate(rar_water,method = "PCoA",distance = "bray")
TMPal <- c('seagreen','steelblue4','hotpink')
TMPal <- c('black', 'grey75', 'hotpink2')


tmp <- plot_ordination(rar_water, BC, color='Group')
Fig_BC_per_season <- tmp$data %>%
  group_by(season, Group, month_number)%>%
  dplyr::mutate(., group=factor(Group, levels=c('Eutrophe','Meso-EU','Mesotrophe','Eu-Meso'))) %>%
  # summarize_at(c('Axis.1', 'Axis.2'), mean) %>%
  ggplot(., aes(Axis.1, Axis.2, color=Group, label=month_number))+
  geom_point(size=1.5)+
  # geom_text(size=4, key_glyph="point")+
  theme_bw()+
  geom_convexhull(aes(fill=Group, group=Group), size=.8, alpha=0.1, show.legend = F)+
  scale_color_manual(values = TMPal)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face ='bold', size=13),
    axis.title.y = element_text(face ='italic', size=15),
    axis.title.x = element_blank(),
    legend.text = element_text(face ='bold', size=13),
    legend.title = element_text(face ='italic', size=15),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    aspect.ratio = 1
  )+
  labs(caption = 'Points separed by seasons keeping the same ordination', x='Axis.1')+
  guides(color=guide_legend(override.aes = list(shape=16)))+
  scale_fill_manual(values = TMPal)+
  facet_wrap2(~factor(season, levels=c('Summer_2021', 'Fall_2021',
                                       'Winter_2021', 'Spring_2022', 
                                       'Summer_2022', 'Fall_2022', 'Winter_2022')), nrow=1,
              labeller = as_labeller(c('Summer_2021' = 'Jun-Aug 2021', 'Fall_2021'='Sep-Nov 2021',
                                       'Winter_2021'='Jan-Fev 2022', 'Spring_2022'='Mar-May 2022', 
                                       'Summer_2022'='Jun-Aug 2022', 'Fall_2022'='Sep-Nov 2022', 
                                       'Winter_2022'='Dec 2022')),
              strip = TMPstrip)+
  labs(color='groupes')
Fig_BC_per_season

#

# SIMPER ------------------------------------------------------------------

total_depth <- sum(taxa_sums(FALL21))
threshold <- 0.001 * total_depth #0.1% abundance
abundant.taxa <-(taxa_sums(FALL21) > threshold)
View(abundant.taxa)

# tax table simper
tax_simper<- FALL21@tax_table %>% t(.) %>% as.data.frame(.) %>% t(.)
tax_simper <- tax_simper %>% as.data.frame()

ps_summer_01_total<-prune_taxa(abundant.taxa,FALL21)
View(abundant.taxa)
# otu table simper
simper_01_total_com.df <- ps_summer_01_total@otu_table %>% as.data.frame() %>% t()
View(simper_01_total_com.df)
simper_total.df<-simper(simper_01_total_com.df,FALL21@sam_data$Group,permutations = 999)


View(simper_total.df)

simp<- summary(simper_total.df, ordered =T)[[2]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
simp_Order <- simp %>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))])
View(simp_Order)
# Alpha -------------------------------------------------------------------
# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
rar_alphadiv <- rar_water@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(rar_water,measures = "Observed")$Observed,
                Shannon=estimate_richness(rar_water,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
View(rar_alphadiv)

# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_rar_alphadiv <- rar_alphadiv %>%
  group_by(Group, month) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
View(mean_rar_alphadiv)
# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rich
mean_GRichness <- mean_rar_alphadiv %>%
  dplyr::mutate(., Group=factor(Group, levels=c('Euplus','Mesomoins','Shift'))) %>%
  mutate(month = recode(month, 
                        'A'='Jun21','B'='Jul21','C'='Aug21','D'='Sep21',
                        'E'='Oct21','F'='Nov21','G'='Jan22','H'='Fev22','I'='Mar22',
                        'J'='Apr22','K'='May22','L'='Jun22','M'='Jul22','N'='Aug22','O'='Sep22',
                        'P'='Oct22','Q'='Nov22','R'='Dec22'))%>%
  dplyr::mutate(., month=factor(month, levels=month_code)) %>%
  ggplot(., aes(x=month, y=mean_Richness, color= Group, group=Group))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.text.y = element_text(face = 'bold', size =11),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=12, angle = 40, vjust =.8),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_Group)+
  #facet_wrap2(vars(Group))+
  geom_smooth(aes(ymin=mean_Richness-sd_Richness, ymax=mean_Richness+sd_Richness), alpha=.1, size = 1)+
  labs(x = '', y = 'Richesse observée moyenne', colour = 'Group')
mean_GRichness
# Even
mean_GEvenness <- mean_rar_alphadiv %>%
  dplyr::mutate(., Group=factor(Group, levels=c('Euplus','Mesomoins','Shift'))) %>%
  mutate(month = recode(month, 
                        'A'='Jun21','B'='Jul21','C'='Aug21','D'='Sep21',
                        'E'='Oct21','F'='Nov21','G'='Jan22','H'='Fev22','I'='Mar22',
                        'J'='Apr22','K'='May22','L'='Jun22','M'='Jul22','N'='Aug22','O'='Sep22',
                        'P'='Oct22','Q'='Nov22','R'='Dec22'))%>%
  dplyr::mutate(., month=factor(month, levels=month_code)) %>%
  ggplot(., aes(x=month, y=mean_Evenness, color= Group, group=Group))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.text.y = element_text(face = 'bold', size =11),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=12, angle = 40, vjust =.8),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_Group)+
  #facet_wrap2(vars(Group))+
  geom_smooth(aes(ymin=mean_Evenness-sd_Evenness, ymax=mean_Evenness+sd_Evenness), alpha=.1, size = 1)+
  labs(x = '', y = 'Evenness observée moyenne', colour = 'Group')
mean_GEvenness


# Alpha troph -------------------------------------------------------------

### HETEROTROPHS ###
taxa_to_keep <- c("Ciliophora",
                  "Kathablepharida","Fungi","Cercozoa",
                  "Bigyra","Choanoflagellata", "Chrompodellids",
                  "Perkinsea","Telonemia_X")
diversite_alpha_het <- subset_taxa(rar_water, Order %in% taxa_to_keep) %>%
  estimate_richness(., measures = "Observed")%>%
  dplyr::mutate(.,month_number = rar_water@sam_data$month_number,
                month = rar_water@sam_data$month,
                Group = rar_water@sam_data$Group, 
                column = rar_water@sam_data$column) %>%
  group_by(Group, month) %>%
  dplyr::mutate(., Group=factor(Group, levels=c("Euplus","Mesomoins","Shift"))) %>%
  summarise(mean_het_Observed = mean(Observed),
            sd_het_Observed = sd(Observed))
diversite_alpha_het
# View(diversite_alpha_het)

g <- diversite_alpha_het %>%
  mutate(month = recode(month, 
                        'A'='Jun21','B'='Jul21','C'='Aug21','D'='Sep21',
                        'E'='Oct21','F'='Nov21','G'='Jan22','H'='Fev22','I'='Mar22',
                        'J'='Apr22','K'='May22','L'='Jun22','M'='Jul22','N'='Aug22','O'='Sep22',
                        'P'='Oct22','Q'='Nov22','R'='Dec22'))%>%
  dplyr::mutate(., month=factor(month, levels=month_code)) %>%
  ggplot(., aes(y = mean_het_Observed, x = month, color = Group, group = Group))+
  geom_smooth(aes(ymin=mean_het_Observed-sd_het_Observed, 
                  ymax=mean_het_Observed+sd_het_Observed), alpha=.1)+
  scale_color_manual(values = palette_Group)+
  theme_bw()+
  labs(title = 'Richesse observée dans le temps chez les hétérotrophes', x = 'mois', y = 'Richesse observée',
       color = 'Groupe')
g

### AUTOTOPHS ###
taxa_to_keep <- c("Gyrista","Chlorophyta_X")
diversite_alpha_aut <- subset_taxa(rar_water, Order %in% taxa_to_keep) %>%
  estimate_richness(., measures = "Observed")%>%
  dplyr::mutate(.,month_number = rar_water@sam_data$month_number,
                month = rar_water@sam_data$month,
                Group = rar_water@sam_data$Group, 
                column = rar_water@sam_data$column) %>%
  group_by(Group, month) %>%
  dplyr::mutate(., Group=factor(Group, levels=c("Euplus","Mesomoins","Shift"))) %>%
  
  summarise(mean_het_Observed = mean(Observed),
            sd_het_Observed = sd(Observed))
g <- diversite_alpha_aut %>%
  mutate(month = recode(month, 
                        'A'='Jun21','B'='Jul21','C'='Aug21','D'='Sep21',
                        'E'='Oct21','F'='Nov21','G'='Jan22','H'='Fev22','I'='Mar22',
                        'J'='Apr22','K'='May22','L'='Jun22','M'='Jul22','N'='Aug22','O'='Sep22',
                        'P'='Oct22','Q'='Nov22','R'='Dec22'))%>%
  dplyr::mutate(., month=factor(month, levels=month_code)) %>%
  ggplot(., aes(y = mean_het_Observed, x = month, color = Group, group = Group))+
  geom_smooth(aes(ymin=mean_het_Observed-sd_het_Observed, 
                  ymax=mean_het_Observed+sd_het_Observed), alpha=.1)+
  scale_color_manual(values = palette_Group_3T)+
  theme_bw()+
  labs(title = 'Richesse observée dans le temps chez les autotrophes', x = 'mois', y = 'Richesse observée',
       color = 'Groupe')
g


# TLAAA -------------------------------------------------------------------

twelve_summer <- subset_samples(rar_water, month_number %in% c('1','2','3','4','5','6', '7','8','9','10','11','12'))
twelve_winter <- subset_samples(rar_water, month_number %in% c('6', '7','8','9','10','11','12','13','14','15','16','17'))


sumBC <- phyloseq::distance(twelve_winter, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = twelve_winter@sam_data$month_number[match(.$sample_1, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., month_2 = twelve_winter@sam_data$month_number[match(.$sample_2, rownames(twelve_winter@sam_data))]) %>%  
  dplyr::mutate(., Chla1 = twelve_winter@sam_data$Chla[match(.$sample_1, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., Chla2 = twelve_winter@sam_data$Chla[match(.$sample_2, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., Group_1 = twelve_winter@sam_data$Group_v2[match(.$sample_1, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., Group_2 = twelve_winter@sam_data$Group_v2[match(.$sample_2, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., col_1 = twelve_winter@sam_data$column[match(.$sample_1, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., col_2 = twelve_winter@sam_data$column[match(.$sample_2, rownames(twelve_winter@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_1 - month_2))%>%
  dplyr::mutate(., DeltaChla = abs(Chla1-Chla2))
rows_to_keep <- distance_sumBC$col_1 == distance_sumBC$col_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$Group_1 == distance_sumBC$Group_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]

distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., Group_1=factor(Group_1, levels=c('Euplus','Mesomoins','CHA', 'CRJ')))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser / nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month)
# filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  # filter(month_2 != 18)
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

## POLY TYPE COL
distance_normal_sumBC <- distance_normal_sumBC %>%
  mutate(poly = case_when(
    Group_1 %in% 'Euplus' ~ 'poly3',
    Group_1 %in% c('Euplus', 'CHA', 'CRJ', 'Mesomoins') ~ 'poly',
    Group_1 %in% c('CTL', 'TRI', 'JAB') ~ 'poly5'))






BC_21_22_time_all.df <- distance_normal_sumBC %>%
  dplyr::mutate(.,color=ifelse(Group_1=="Euplus","chartreuse2",
                               ifelse(Group_1=="Mesomoins","skyblue",
                                      ifelse(Group_1=="CHA","forestgreen",'royalblue'))))
BC_21_22_time_all.df<- BC_21_22_time_all.df %>% as.data.frame() %>%
  dplyr::mutate(., Group_1=factor(Group_1, levels=c('Euplus','Mesomoins','CHA','CRJ')))
Fig_BC_21_22_time_all <- BC_21_22_time_all.df %>%
  ggplot(.,aes(x=normal_month_gap,y=BC_dissimilarity,color=color,group=Group_1,fill=color))+
  #best poly model poly3
  stat_poly_line(alpha=0.2,formula = lm_formula_poly3,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly3"),])+
  #best poly model poly4
  stat_poly_line(alpha=0.2,formula = lm_formula_poly4,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly4"),])+
  #best poly model poly5
  stat_poly_line(alpha=0.2,formula = lm_formula_poly5,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0,1))+
  coord_cartesian(ylim = c(0.45, 1))+
  scale_color_identity()+
  scale_fill_identity()+
  labs(y="Dissimilarité de composition taxonomique microeucaryote",
       x="Jours séparant les échantillons")
Fig_BC_21_22_time_all
####______ lake 

## POLY TYPE COL

distance_normal_sumBC <- distance_normal_sumBC %>%
  dplyr::mutate(.,color=ifelse(Group_1=="Euplus","palegreen3",
                               ifelse(Group_1=="Mesomoins","skyblue",
                                      ifelse(Group_1=="CHA","forestgreen","royalblue"))))
Fig_BC_21_22_time_Group<- distance_normal_sumBC %>%
  ggplot(.,aes(x=normal_month_gap,y=BC_dissimilarity,color=color,group=Group_1,fill=Group_1, label='poly'))+
  #points
  geom_point(size=1,alpha=0.4,color="grey",show.legend = F)+
  #linear model
  geom_text(x=110,y=0.21,size=3,check_overlap = T,fontface="italic",label="linear",color="gray32")+
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,formula = lm_formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,formula = lm_formula_linear,
               label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #best poly model poly3
  geom_text(x=110,y=0.132,size=3,check_overlap = T,fontface="italic",
            data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly3,show.legend = F,
                 data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,
               data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=8),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0.1,1))+
  scale_color_identity()+
  scale_fill_manual(values = TMPal)+
  #scale_color_manual(values = rev(palette_Group_3T))+
  labs(y="",
       x="Jours séparant les échantillons")+
  facet_wrap2(~ Group_1,nrow = 1,strip=strip_color_Group,scales = "fixed")
Fig_BC_21_22_time_Group

ggarrange(Fig_BC_21_22_time_all, Fig_BC_21_22_time_Group, nrow = 2)
TLA <- Fig_BC_21_22_time_all+Fig_BC_21_22_time_Group+
  plot_layout(guides = 'collect')+
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")
#
TLA
# 

# TLa ---------------------------------------------------------------------


sumBC <- phyloseq::distance(rar_water, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = rar_water@sam_data$month_number[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., month_2 = rar_water@sam_data$month_number[match(.$sample_2, rownames(rar_water@sam_data))]) %>%  
  dplyr::mutate(., Chla1 = rar_water@sam_data$Chla[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Chla2 = rar_water@sam_data$Chla[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Group_ID = rar_water@sam_data$Group[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., Group_2 = rar_water@sam_data$Group[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., col_1 = rar_water@sam_data$column[match(.$sample_1, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., col_2 = rar_water@sam_data$column[match(.$sample_2, rownames(rar_water@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_1 - month_2))%>%
  dplyr::mutate(., DeltaChla = abs(Chla1-Chla2))
rows_to_keep <- distance_sumBC$col_1 == distance_sumBC$col_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$Group_ID == distance_sumBC$Group_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]

distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., Group_ID=factor(Group_ID, levels=c("Euplus","Mesomoins","Shift")))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser / nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month)
# filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  # filter(month_2 != 18)
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

# ## POLY TYPE COL
# distance_normal_sumBC <- distance_normal_sumBC %>%
#   mutate(poly = case_when(
#     Group_ID %in% 'VAI' ~ 'poly3',
#     Group_ID %in% c('BOI', 'CHA', 'CRJ', 'GDP', 'VER') ~ 'poly4',
#     Group_ID %in% c('CTL', 'TRI', 'JAB') ~ 'poly5'))

### NORMAL_BC_TLA_plot
tmp <- ggplot(data =distance_normal_sumBC, aes(x = normal_month_gap, y = BC_dissimilarity, color = Group_ID, fill = Group_ID, label=poly))+
  geom_jitter(size=01, alpha=0.4, color='grey')+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  geom_smooth(method = 'loess', size=1.2)+
  # geom_smooth(inherit.aes = FALSE, method = 'lm',
  #             mapping = aes(x = normal_month_gap, y = BC_dissimilarity),
  #             color = 'black', size=1)+
  facet_wrap2(~ Group_ID, scales = "fixed", strip = strip_color_Group)+
  # scale_x_continuous(breaks = c(22, 79, 171, 386, 542),
  #                    labels = c('1', '3', '6', '12', '18'))+
  scale_color_manual(values = palette_Group_3T)+
  scale_fill_manual(values = palette_Group_3T)+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  #linear model
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,
                 formula = formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,
               formula = formula_linear,label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #polynomial model
  stat_poly_line(formula = formula_poly4) +
  stat_poly_eq(use_label(c("R2","P")), size=3,
               formula = formula_poly4,
               label.y = 0.01, label.x = 0.97,rr.digits = 2)+
  
  
  labs(title = 'Dissimilarité de Bray Curtis en fonction du nombre de mois séparant les échantillons',
       x = "Nombre de jours d'écart entre",
       y = "Dissimilarité de Bray-Curtis")
tmp
View(distance_normal_sumBC)

### NORMAL_BC_TLA_plot
Allgroup <- distance_normal_sumBC %>%
  dplyr::mutate(., Group_ID=factor(Group_ID, levels=c("Euplus","Mesomoins","Shift")))%>%
  ggplot(., aes(x = normal_month_gap, y = BC_dissimilarity, group=Group_ID, color = Group_ID, fill = Group_ID))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks=element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=16),
        axis.text = element_text(size=15, face = 'bold'),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_fill_manual(values = palette_Group)+
  geom_smooth(formula = lm_formula_poly3)+
  # geom_smooth(inherit.aes = FALSE, method = 'lm',
  #             mapping = aes(x = month_gap, y = BC_dissimilarity, group=Group_ID, color = Group_ID), size=.2)+
  scale_color_manual(values = palette_Group)+
  # scale_x_continuous(limit=c(0,17.5), expand=c(0,0), breaks = c(1,3,6,9,12,16), labels = c(1,3,6,9,12,16))+
  labs(title = 'Time Lag Analysis', fill = 'groupe', colour = 'groupe',
       x = "Nombre de jours d'écart",
       y = "Dissimilarité de Bray-Curtis")+
  coord_cartesian(ylim = c(0.63, 0.9), xlim = c(22, 542))
# caption = "L'axe y débute à 0.5 car aucune valeur inférieure n'a été obtenue"

# TEst --------------------------------------------------------------------


palette_Order <- c('#17380a','#37a03b',
                   '#C3C04B','#F1FF7E','#d28951',
                   '#b87347','#9e5d3c','#844931',
                   '#6a3726','#50251b','red3', 
                   'tomato','grey','#4A4A4A')

Cleaned_orders <- subset_taxa(rar_water, Order %notin% 'Streptophyta_X')

Order_AR <-  Cleaned_orders %>% # Name of your final object and your current phyloseq object
  tax_glom(.,taxrank="Order") %>% #reducing the taxonomy (tax_glom) at the Order levelfor a barplot
  psmelt(.) %>% #melting from a wide (more columns than rows) to a long (more rows than columns,computer-friendly) dataframe)
  dplyr::group_by(lake, month, Order) %>% #using group_by to compute the Abundance of each Phyla BY each SampleType) 
  summarise(median_count=median(Abundance)) %>% #compute the median count for each group_by-designated group #another group_by to get the overall count for each SampleType
  mutate(median_abundance=as.numeric(
    paste0((round(median_count/sum(median_count),4))*100))) %>% #compute the relative abundance for each SampleType and keep only 4 digits
  dplyr::mutate(., Order_legend=ifelse(median_abundance<4,"Ordres < 4%",Order)) %>%
  mutate(Order_legend = recode(Order_legend, 'Cryptophyta_X' = 'Cryptophyta', 
                               'Telonemia_X' = 'Telonemia', 
                               'Chlorophyta_X' = 'Chlorophyta', 
                               'Cryptophyta_X:nucl' = 'Cryptophyta')) %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI", 
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  dplyr::mutate(., Order_legend=factor(Order_legend, levels=c(
    "Gyrista","Chlorophyta",
    "Cryptophyta","Dinoflagellata","Ciliophora",
    "Kathablepharida","Fungi","Cercozoa",
    "Bigyra","Choanoflagellata","Perkinsea",
    "Telonemia", "Chrompodellids", "Ordres < 4%")))%>%
  ggplot(.,aes(x="",y=median_abundance,fill=Order_legend))+
  theme_void()+
  geom_col()+
  coord_polar(theta = 'y')+
  geom_text(aes(label = 'y'),
            position = position_stack(vjust = 0.5))+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.x = element_text(size = 12),
        axis.ticks=element_blank(),
        panel.grid = element_blank())+
  expand_limits(y = c(0, 0)) + #remove extra space on the x axis
  scale_y_continuous(expand = c(0,0))+ #keep a little extra space on the y axis #set a color palette for the barplot
  scale_fill_manual(values=palette_Order)+ #set a color palette for the barplot
  labs(y = "Abondance relative (%)", fill="Ordre", x = "Mois (1=Janvier, 2=Février...)") #change your axis and legend title 
Order_AR
View(Cleaned_orders@tax_table)
View(Order_AR$data)
write.csv(Order_AR$data, 'ABD.csv')

palpie <- c('#37a03b','#C3C04B','#9e5d3c','tomato','grey50')

PieDF <- data.frame(group = c('Autotrophes','Mixotrophes','Hétérotrophes','Parasites','NA & <4%'), abdrel = c(406875.5,692282,521627.5,17073,158110.5))
PieDF <- data.frame(group = c('Autotrophes','Mixotrophes','Hétérotrophes','Parasites','<4%'), abdrel = c(22.65,38.55,29.86,0.95,7.99))

PieDF %>%
  dplyr::mutate(., group=factor(group, levels=c('Autotrophes','Mixotrophes',
                                                'Hétérotrophes','Parasites','<4%')))%>%
  ggplot(., aes(x="", y=abdrel, fill=group))+
  theme_void()+
  geom_col(color='black')+
  coord_polar(theta = 'y')+
  geom_label(aes(label = abdrel),
             position = position_stack(vjust=0.4),
             show.legend = FALSE)+
  scale_fill_manual(values=palpie)+
  theme(legend.text = element_text(size=11))+
  labs(fill = 'fonctions écologiques')
View(PieDF)

###__ 0.1% Filter ###

total_depth <- sum(taxa_sums(rar_water))
threshold <- 0.001 * total_depth #0.1% abundance
abundant.taxa <-(taxa_sums(rar_water) > threshold)
View(abundant.taxa)

# tax table simper
tax_simper<- rar_water@tax_table %>% t(.) %>% as.data.frame(.) %>% t(.)
tax_simper <- tax_simper %>% as.data.frame()

ps_summer_01_total<-prune_taxa(abundant.taxa,rar_water)

# otu table simper
simper_01_total_com.df <- ps_summer_01_total@otu_table %>% as.data.frame() %>% t()
View(simper_01_total_com.df)
simper_total.df<-simper(simper_01_total_com.df,rar_water@sam_data$Group,permutations = 999)
simp<- summary(simper_total.df, ordered =T)[[1]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
View(simper_total.df)
# simpBOI<- summary(simper_BOI_vs_others.df, ordered =T)[[1]] %>% as.data.frame() %>%
#   subset(.,cumsum <= 0.7 & p <=0.001)
# list_simper_BOI<-rownames(simp)
View(simp %>% dplyr::mutate(.,Phylum=tax_simper$Phylum[match(rownames(.),rownames(tax_simper))]))
View(simp %>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))
# # BOI vs others
# simper_BOI_vs_others.df<-simper(simper_01_total_com.df,rar_water@sam_data$Boi,permutations = 999)
# View(simpBOI %>% dplyr::mutate(.,Phylum=tax_simper$Phylum[match(rownames(.),rownames(tax_simper))]))
# View(simpBOI %>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))

View(simp)

View(simp$Euplus_Mesomoins %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))
View(simp$Euplus_Shift %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))
View(simp$Shift_Mesomoins %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))]))



a <- simp$Euplus_Mesomoins %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))])
View(a)
b <- simp$Euplus_Shift %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))])
View(b)
c <- simp$Shift_Mesomoins %>% subset(., p <= 0.001 & cumsum <= 0.7)%>% dplyr::mutate(.,Order=tax_simper$Order[match(rownames(.),rownames(tax_simper))])
View(c)







####
#####~Rare_ASV ----------------------------------------------------------------
# ASVs dépassant les 1% ABDreL dans au moins un échantillon
Rare <- read.csv('ASV_rare.csv', header=TRUE, sep = ';')
Rare[[482]] <- NULL
Rare[[481]] <- NULL
Rare <- na.omit(Rare)
View(Rare)
Rare_list <- list(Rare$X)
View(Rare_list)
RARE <- prune_taxa(x=rownames(rar_water@tax_table), taxa = Rare_list)

tmp_no_rare <-(rownames(rar_water@tax_table) %in% Rare$X)
tmp_rare <-(rownames(rar_water@tax_table) %notin% Rare$X)
View(tmp_rare)
RARE<-prune_taxa(tmp_rare,rar_water)
NO_RARE <- prune_taxa(tmp_no_rare, rar_water)



View(RARE@otu_table)
View(NO_RARE@otu_table)





####
# Rare ASV Alphadiv -------------------------------------------------------

# dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~
rar_alphadiv <- NO_RARE@sam_data %>% as.tbl() %>% as.data.frame() %>% 
  dplyr::mutate(., Richness=estimate_richness(NO_RARE,measures = "Observed")$Observed,
                Shannon=estimate_richness(NO_RARE,measures = "Shannon")$Shannon, 
                Evenness=Shannon/log(Richness))
# mean dataframe ~~~~~~~~~~~~~~~~~~~~
mean_rar_alphadiv <- rar_alphadiv %>%
  group_by(lake, month) %>%
  summarise(mean_Shannon = mean(Shannon),
            mean_Richness = mean(Richness),
            mean_Evenness = mean(Evenness),
            sd_Shannon = sd(Shannon),
            sd_Richness = sd(Richness),
            sd_Evenness = sd(Evenness))
# Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rich
mean_Richness_rare <- mean_rar_alphadiv %>%
  dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "TRI", 
                                              "CRJ", "CTL", "BOI",
                                              "GDP", "CHA", "VER"))) %>%
  mutate(month = recode(month, 
                        'A'='Jun21','B'='Jul21','C'='Aug21','D'='Sep21',
                        'E'='Oct21','F'='Nov21','G'='Jan22','H'='Fev22','I'='Mar22',
                        'J'='Apr22','K'='May22','L'='Jun22','M'='Jul22','N'='Aug22','O'='Sep22',
                        'P'='Oct22','Q'='Nov22','R'='Dec22'))%>%
  dplyr::mutate(., month=factor(month, levels=month)) %>%
  ggplot(., aes(x=month, y=mean_Richness, color= lake, group=lake))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=19,face="italic"),
        axis.ticks.y =element_blank(),
        axis.text.y = element_text(face = 'bold', size =11),
        axis.title = element_text(size=16),
        aspect.ratio=1,
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, angle = 40, vjust =.8),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_color_manual(values = palette_lake_3T)+
  #facet_wrap2(vars(lake))+
  geom_smooth(aes(ymin=mean_Richness-sd_Richness, ymax=mean_Richness+sd_Richness), alpha=.1, size = 1)+
  labs(x = 'ASV > 1% ABDrel dans au moins un échantillon', y = 'Richness moyenne', colour = 'Lac')+
  coord_cartesian(ylim = c(0, 450))

mean_Richness_rare

Fig_rich <- mean_Richness + mean_Richness_rare
View(mean_Richness$data)
View(mean_Richness_rare$data)

# TLA STATS RAREs ------------------------------------

sumBC <- phyloseq::distance(NO_RARE, method = 'bray') %>%
  as.matrix() %>%
  as.tibble() %>%
  as.data.frame()%>%
  dplyr::mutate(., rownames=colnames(.)) %>%
  melt(., value.name = 'BC_dissimilarity', id = 'rownames') %>%
  rename('sample_1' = 'rownames', 'sample_2' = 'variable') %>%
  dplyr::mutate(., month_1 = NO_RARE@sam_data$month_number[match(.$sample_1, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., month_2 = NO_RARE@sam_data$month_number[match(.$sample_2, rownames(NO_RARE@sam_data))]) %>%  
  dplyr::mutate(., Chla1 = NO_RARE@sam_data$Chla[match(.$sample_1, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., Chla2 = NO_RARE@sam_data$Chla[match(.$sample_2, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., lake_ID = NO_RARE@sam_data$lake[match(.$sample_1, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., lake_2 = NO_RARE@sam_data$lake[match(.$sample_2, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., col_1 = NO_RARE@sam_data$column[match(.$sample_1, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., col_2 = NO_RARE@sam_data$column[match(.$sample_2, rownames(NO_RARE@sam_data))]) %>%
  dplyr::mutate(., month_compare = paste0(month_1, '_', month_2))
sumBC$month_1 <- as.numeric(as.character(sumBC$month_1))
sumBC$month_2 <- as.numeric(as.character(sumBC$month_2))
#
distance_sumBC <- sumBC %>%
  dplyr::mutate(.,month_gap = abs(month_1 - month_2))%>%
  dplyr::mutate(., DeltaChla = abs(Chla1-Chla2))
rows_to_keep <- distance_sumBC$col_1 == distance_sumBC$col_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$lake_ID == distance_sumBC$lake_2
distance_sumBC <- distance_sumBC[rows_to_keep, ]
rows_to_keep <- distance_sumBC$month_2 > distance_sumBC$month_1
distance_sumBC <- distance_sumBC[rows_to_keep, ]

distance_sumBC <- distance_sumBC %>%
  dplyr::mutate(., lake_ID=factor(lake_ID, levels=c("JAB", "VAI", "TRI", 
                                                    "CRJ", "CTL", "BOI",
                                                    "GDP", "CHA", "VER")))
#
normal <- data.frame(
  month_1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
              11, 12, 13, 14, 15, 16, 17, 18),
  month_1_numeric = c(1, 30, 65, 100, 128, 157, 220, 248, 276, 311,
                      337, 372, 401, 429, 464, 492, 521, 543)
)
# Normaliser / nombre de jours
normal <- normal %>%
  mutate(diff_month = c(diff(month_1_numeric), NA)) %>%
  mutate(normalized_value = month_1_numeric / diff_month)
# filter(month_1_numeric != 18)
# 
distance_sumBC <- distance_sumBC %>%
  # filter(month_2 != 18)
  dplyr::mutate(., normalized_value = normal$month_1_numeric[match(.$month_1, normal$month_1)])%>%
  dplyr::mutate(., normalized_value2 = normal$month_1_numeric[match(.$month_2, normal$month_1)])
distance_normal_sumBC <- distance_sumBC %>%
  dplyr::mutate(.,normal_month_gap = abs(normalized_value2 - normalized_value))

## POLY TYPE COL
distance_normal_sumBC <- distance_normal_sumBC %>%
  mutate(poly = case_when(
    lake_ID %in% 'VAI' ~ 'poly5',
    lake_ID %in% c('BOI', 'CHA', 'CRJ', 'GDP', 'VER') ~ 'poly5',
    lake_ID %in% c('CTL', 'TRI', 'JAB') ~ 'poly5'))


BC_21_22_time_all.df <- distance_normal_sumBC %>%
  dplyr::mutate(.,color=ifelse(lake_ID=="JAB","#2F6B9D",
                               ifelse(lake_ID=="VAI","#5C9BD6",
                                      ifelse(lake_ID=="TRI","#F2AB8C",
                                             ifelse(lake_ID=="CRJ","#EB7947",
                                                    ifelse(lake_ID=="CTL","#E45000",
                                                           ifelse(lake_ID=="BOI","#E40000",
                                                                  ifelse(lake_ID=="GDP","#980000",
                                                                         ifelse(lake_ID=="CHA","#1E7A36","#28A448")))))))))
BC_21_22_time_all.df<- BC_21_22_time_all.df %>% as.data.frame() %>%
  dplyr::mutate(., lake_ID=factor(lake_ID, levels=c("JAB", "VAI", "TRI", 
                                                    "CRJ", "CTL", "BOI",
                                                    "GDP", "CHA", "VER")))
Fig_BC_21_22_time_all <- BC_21_22_time_all.df %>%
  ggplot(.,aes(x=normal_month_gap,y=BC_dissimilarity,color=color,group=lake_ID,fill=color))+
  #best poly model poly3
  stat_poly_line(alpha=0.2,formula = lm_formula_poly3,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly3"),])+
  #best poly model poly4
  stat_poly_line(alpha=0.2,formula = lm_formula_poly4,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly4"),])+
  #best poly model poly5
  stat_poly_line(alpha=0.2,formula = lm_formula_poly5,show.legend = F,
                 data =BC_21_22_time_all.df[(BC_21_22_time_all.df$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0,1.2))+
  coord_cartesian(ylim = c(0.35, 1))+
  scale_color_identity()+
  scale_fill_identity()+
  labs(y="Dissimilarité de composition taxonomique microeucaryote",
       x="Jours séparant les échantillons")
Fig_BC_21_22_time_all
####______ lake 
distance_normal_sumBC <- distance_normal_sumBC %>%
  dplyr::mutate(.,color=ifelse(lake_ID=="JAB","#2F6B9D",
                               ifelse(lake_ID=="VAI","#5C9BD6",
                                      ifelse(lake_ID=="TRI","#F2AB8C",
                                             ifelse(lake_ID=="CRJ","#EB7947",
                                                    ifelse(lake_ID=="CTL","#E45000",
                                                           ifelse(lake_ID=="BOI","#E40000",
                                                                  ifelse(lake_ID=="GDP","#980000",
                                                                         ifelse(lake_ID=="CHA","#1E7A36","#28A448")))))))))
Fig_BC_21_22_time_lake<- distance_normal_sumBC %>%
  ggplot(.,aes(x=normal_month_gap,y=BC_dissimilarity,color=color,group=lake_ID,fill=lake_ID,label=poly))+
  #points
  geom_point(size=1,alpha=0.4,color="grey",show.legend = F)+
  #linear model
  geom_text(x=110,y=0.21,size=3,check_overlap = T,fontface="italic",label="linear",color="gray32")+
  stat_poly_line(color="gray32",fill="gray32",alpha=0.4,formula = lm_formula_linear) +
  stat_poly_eq(use_label(c("R2","P")),color="gray32",size=3,formula = lm_formula_linear,
               label.y = 0.1, label.x = 0.97,rr.digits =2)+
  #best poly model poly3
  geom_text(x=110,y=0.132,size=3,check_overlap = T,fontface="italic",
            data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly3,show.legend = F,
                 data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,
               data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly3"),])+
  #best poly model poly4
  geom_text(x=110,y=0.132,size=3,check_overlap = T,fontface="italic",
            data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly4"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly4,show.legend = F,
                 data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly4"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,
               data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly4"),])+
  #best poly model poly5
  geom_text(x=110,y=0.132,size=3,check_overlap = T,fontface="italic",
            data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly5"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly5,show.legend = F,
                 data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly5"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly5,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,
               data =distance_normal_sumBC[(distance_normal_sumBC$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=8),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0.1,1.2))+
  scale_color_identity()+
  scale_fill_manual(values = palette_lake_3T)+
  #scale_color_manual(values = rev(palette_lake_3T))+
  labs(y="",
       x="Jours séparant les échantillons")+
  facet_wrap2(~ lake_ID,nrow = 3,strip=strip_color_lake,scales = "fixed")
Fig_BC_21_22_time_lake


TLA <- Fig_BC_21_22_time_all+Fig_BC_21_22_time_lake+
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")
#
TLA

View(rar_water@sam_data)
# PERMANOVA GROUP ####
# matrice de dist
SubNoTri <- subset_samples(rar_water, month_number %notin% '1')
features_table<- t(SubNoTri@otu_table) %>% as.data.frame()
dist_matrix <- features_table %>% vegdist(., method="bray")

perm_group<-adonis2(dist_matrix ~ SubNoTri@sam_data$Group_v2,method =
                      "bray",permutations = 999)
perm_group #<0.001 **
perm_group_betadisp<-as.data.frame(anova(betadisper(dist_matrix,SubNoTri@sam_data$Group,type="centroid")))
perm_group_betadisp #<0.0001 ***
pair_perm_group<-pairwise.adonis(dist_matrix,SubNoTri@sam_data$Group,p.adjust.m = "bonferroni",perm = 999)
#pair.perm_PER_WU_lake

View(perm_group)
View(perm_group_betadisp)
View(pair_perm_group)

# PERMANOVA SEASON####
# matrice de dist
SubNoTri <- subset_samples(rar_water, month_number %notin% '1')
features_table<- t(SubNoTri@otu_table) %>% as.data.frame()
dist_matrix <- features_table %>% vegdist(., method="bray")



perm_season<-adonis2(dist_matrix ~ SubNoTri@sam_data$season,method =
                       "bray",permutations = 999)
perm_season #<0.001 **
perm_season_betadisp<-as.data.frame(anova(betadisper(dist_matrix,SubNoTri@sam_data$season,type="centroid")))
perm_season_betadisp #<0.0001 ***
pair_perm_season<-pairwise.adonis(dist_matrix,SubNoTri@sam_data$season*SubNoTri@sam_data$Group,p.adjust.m = "bonferroni",perm = 999)
#pair.perm_PER_WU_lake

View(perm_season)
View(perm_season_betadisp)
View(pair_perm_season)









write.csv(distance_normal_sumBC, 'TLA$Data.csv')





# PERMANOVA ---------------------------------------------------------------

SAM <- rar_water@sam_data %>% as.data.frame()
OTU <- rar_water@otu_table %>% as.data.frame()
dist_matrix <- vegdist(rar_water@otu_table, method = "bray")

#  permanova avec vegan
permanova_ko<-adonis2(MDS_FUNC.hellinger~ Status + Month,data=dist_matrix,permutations=999,method="bray")
summary(perm)


View(BC$vectors)
perm<-adonis2(OTU ~ Group + Season,data = dist_matrix, permutations=999, method="bray", na.rm = TRUE)
summary(perm)

# lmm avec lme4 ET lmerTest

testlm<-lmer(turnover~ Status + month_number + (1|lake),data=subset(Fig_box_ASV_KO_turnover$data,features=="KOs"),REML = F)
summary(testlm)

# MOTA - Params centerscaled ----------------------------------------------

View(df_cumsum_params)


BASE <- rar_water@sam_data %>% as.data.frame
numerical <- rar_water@sam_data %>% subset(., select=c(Temperature, pH, NO3NO2, Salinity, PO4, NH4, Chla, TPC, TPN))
numerical$Temperature <- as.numeric(numerical$Temperature)
numerical$pH <- as.numeric(numerical$pH)
numerical <- numerical %>%
  scale(., center =T, scale =T)
numerical <- as.data.frame(abs(numerical))
numerical$lake = BASE$lake
numerical$month=BASE$month
numerical <- numerical %>%
  filter(month != 'A')
numerical_medians <- numerical %>% group_by(lake, month) %>% summarize(across(where(is.numeric), median, na.rm = TRUE))
df_cumsum_params <- numerical_medians %>%
  group_by(lake) %>%
  mutate(across(where(is.numeric), ~ ifelse(month == 'B', 0, .), .names = "cumsum_{.col}")) %>%
  mutate(across(starts_with("cumsum_"), cumsum))



BASE <- rar_water@sam_data
transform(BASE, Chla = as.numeric(Chla))
median_BASE <- BASE %>%
  group_by(lake, month) %>%
  summarise(Chla = median(Chla))
filt_BASE <- median_BASE %>%
  filter(month != 'A')
write.csv(filt_BASE, 'tmp_chla_values.csv')
df_cumsum_Chla <- read.csv('tmp_chla_values.csv', header = TRUE, sep=';')







View(df_cumsum_Chla)








# subset.data.frame(., lake %in% 'BOI')
CHLA_1 <- CHLA_raw  %>%
  mutate(diff = mean_Chla - lag(mean_Chla, default = first(mean_Chla)))
CHLA_2 <- na.omit(CHLA_1) 
create_new_column <- function(column) {
  new_column <- numeric(length(column))
  new_column[1] <- 0  # Mettre 0 comme première valeur
  for (i in 2:length(column)) {
    new_column[i] <- column[i - 1] + new_column[i - 1]
  }
  return(new_column)
}
coord <- create_new_column(CHLA_2$diff)
View(coord)
CHLA_3 <- CHLA_2
CHLA_3$diff <- abs(CHLA_3$diff)



DT <- data.table(CHLA_3, key = "lake")
DT[, csum := cumsum(diff), by = key(DT)]

View(DT)

CHLA_df <- melt(DT, value.name = DT$lake)

DT_1 <- reshape(DT, idvar = "month_number", timevar = "lake", direction = "wide")
View(DT_1)





palette_season <- c('#B7D7FF', '#FFFED5', '#FFD9C8', '#E9D2A4' )

MOTA_Params <- BASE %>% as.data.frame() %>%
  filter(!lake=='VER') %>%
  dplyr::mutate(lake_month=paste0(lake,"_",month),
                season=if_else(month %in% c("7","8", "18"),"hiver",
                               if_else(month %in% c("15","16","17","4","5","6"),"automne",
                                       if_else(month %in% c("9","10","11"),"printemps","été"))),
                season=factor(season,levels = c("hiver","printemps","été","automne"))) %>%
  group_by(lake) %>%
  dplyr::mutate(.,max_dist=max(csum)) %>%
  ggplot(.,aes(x=csum,y=(fct_reorder(lake,max_dist)),label=month,fill=season))+
  geom_line(size=0.5,color="black",show.legend = F )+
  geom_line(size=0.5,color="black",show.legend = F)+
  geom_point(show.legend = T)+
  geom_label(show.legend = F)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    plot.caption =element_text(face="italic",size=10,hjust=1),
    axis.title.x=element_text(size=12,face="bold"),axis.title.y=element_blank(),
    axis.text.y=element_text(size=14,face="bold",
                             color=rev(c("#52C655","#52C655","#473A3B","#52C655","#52AAC6","#52C655","#52AAC6","#473A3B","#52AAC6"))),
    axis.text.x=element_text(size=10, face = 'bold'),
    axis.ticks = element_blank(),
    legend.title = element_text(size=12,face = "bold",hjust=0),
    legend.text = element_text(size=11),
    legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=palette_season)+
  labs(caption=paste0("MOTA – Params"),
       x="~ Distance parcourue, cumsum des params centrés réduits",fill="Saison")+
  guides(fill=guide_legend(override.aes=aes(shape=22,size=8),nrow=2))

MOTA_Params

























# Correl MOTA - BC/PARAMS-------------------------------------------------------------


View(df_cumsum_params) # paramètres colonne VALUE
View(df_cumsum_BC) # paramètres colonne VALUE (month = time_)

final_corr_mot <- read.csv('C:/Users/Amaury/OneDrive/Bureau/Rpertory/motdata.csv', header=TRUE, sep=';')

write.csv(df_cumsum_BC, 'BC_mot.csv')
write.csv(df_cumsum_params, 'PARAMS_mot.csv')
df_cumsum_params <- read.csv('PARAMS_mot.csv', header = TRUE, sep = ';')


final_corr_mot %>% dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VAI", "CRJ2", 
                                                               "CRJ1", "CTL", "BOI",
                                                               "GDP", "CHA", "VER"))) %>%
  ggplot(., aes(x = mota, y=param, gr))+
  facet_wrap2(~factor(lake), strip = strip_color_lake)+
  stat_poly_line(color = 'black') +
  stat_poly_eq(use_label(c("R2", "p"))) +
  geom_point(size=3, aes(color=lake))+
  # geom_convexhull(aes(fill=lake, group=lake), size=.8, alpha=0.1, show.legend = FALSE)+
  # geom_label(size=5, aes(label = lake, color = lake), show.legend = FALSE)+
  scale_color_manual(values = palette_lake_3T)+
  scale_fill_manual(values = palette_lake_3T)+
  labs(title = 'sommes cumulées des variations en compositions taxonomiques (BC) et paramètres physicochimiques (centrés réduits) des mois 2 à 18',
       x = 'Variations de composition des communautés (Bray-Curtis)',
       y= 'Variations des paramètres physicochimiques (Centré-Réduit)')




# Correl MOTA - BC/CHLA + NORMALIZED x=y ---------------------------------------------------

final_corr_mot <- read.csv('C:/Users/Amaury/OneDrive/Bureau/Rpertory/motdatah.csv', header=TRUE, sep=';')
PAL<-c('#00FF00', '#33FF33', '#66FF66', '#33FF99', '#00FFCC', '#00FFFF', '#33CCFF', '#3399FF', '#0000FF')
PAL<- rev(c('#2E7D32', '#388E3C', '#43A047', '#66BB6A', '#4CAF50', '#39A1D0', '#0288D1', '#0277BD', '#01579B'))
PAL <-rev(c('#1B5E20', '#2E7D32', '#388E3C', '#4CAF50', '#66BB6A', '#42A5F5', '#1E88E5', '#1976D2', '#1565C0')
)
PAL<-c('darkblue', 'blue', '#1E88E5', 'lightblue', 'seagreen1', 'darkseagreen3', 'olivedrab1', 'seagreen', 'darkslategray')


final_corr_mot %>% dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VSM", "CER-L", 
                                                               "CSM", "CER-S", "CRE",
                                                               "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(x = Chla, y=mota, group=lake))+
  # facet_wrap2(~factor(lake), strip = strip_color_lake)+
  stat_poly_line(color='black') +
  #stat_poly_eq(use_label(c("R2", "p"))) +
  geom_point(size=3, aes(color=lake), show.legend = TRUE)+
  # geom_convexhull(aes(fill=lake, group=lake), size=.8, alpha=0.1, show.legend = FALSE)+
  # geom_label(size=5, aes(label = lake, color = lake), show.legend = FALSE)+
  scale_color_manual(values = palette_lake_3T)+
  scale_fill_manual(values = palette_lake_3T)
#labs(title = 'sommes cumulées des variations en compositions taxonomiques (BC) et concentration en Chla des mois 2 à 18',
       #x = " Variations de composition des communautés (sommes cumulées de l'indice Bray-Curtis)",
       #y= 'Variations de concentration en Chla (sommes cumulées des concentrations dans le temps)')
#normaliser les données


View(final_corr_mot)


normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
final_corr_mot_normalized <- final_corr_mot %>%
  group_by(lake) %>%
  mutate(
    mota_normalized = normalize(mota),
    Chla_normalized = normalize(Chla)
  ) %>%
  ungroup()

# plot normalized

final_corr_mot_normalized %>% dplyr::mutate(., lake=factor(lake, levels=c("JAB", "VSM", "CER-L", 
                                                                          "CSM", "CER-S", "CRE",
                                                                          "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(x = Chla_normalized, y=mota_normalized, group=lake))+
  facet_wrap2(~factor(lake), strip = v2strip_color_lake)+
  theme_bw()+
  stat_poly_line(color='black', size=.5) +
  geom_line()+
  #stat_poly_eq(use_label(c("R2", "p"))) +
  geom_point(size=3, aes(color=lake), show.legend = FALSE)+
  # geom_convexhull(aes(fill=lake, group=lake), size=.8, alpha=0.1, show.legend = FALSE)+
  # geom_label(size=5, aes(label = lake, color = lake), show.legend = FALSE)+
  scale_color_manual(values = palette_lake_3T)+
  scale_fill_manual(values = palette_lake_3T)+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=1)+
  coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  annotate("text", x = 0.2, y = 0.5, label = "x = y", color = "red", size = 5, fontface = "bold") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, 
               label.x = 0.05, label.y = 0.95,  
               parse = TRUE, size = 3.5)
#labs(title = 'sommes cumulées des variations en compositions taxonomiques (BC) et concentration en Chla des mois 2 à 18',
#x = " Variations de composition des communautés (sommes cumulées de l'indice Bray-Curtis)",
#y= 'Variations de concentration en Chla (sommes cumulées des concentrations dans le temps)')




# Chla-BCSeasons ----------------------------------------------------------



View(SubSHIFT)

Chla <- rar_water@sam_data %>% as.data.frame() %>%
  group_by(season, lake, month_number)%>%
  dplyr::mutate(., lake=factor(lake, levels=c( "CHA","CRJ"))) %>%
  summarize_at(c('Chla'), median) %>%
  na.omit() %>%
  ggplot(., aes(lake, Chla, color=lake, fill = lake))+
  theme_bw()+
  geom_col(linewidth=.5)+
  scale_fill_manual(values = c('#005B96','#228B22'))+
  
  scale_color_manual(values = c('black', 'black'))+
  facet_wrap2(~factor(season, levels=c('Summer_2021', 'Fall_2021',
                                       'Winter_2021', 'Spring_2022', 
                                       'Summer_2022', 'Fall_2022', 'Winter_2022')), nrow=1,
              labeller = as_labeller(c('Summer_2021' = 'Jun-Aug 2021', 'Fall_2021'='Sep-Nov 2021',
                                       'Winter_2021'='Jan-Fev 2022', 'Spring_2022'='Mar-May 2022', 
                                       'Summer_2022'='Jun-Aug 2022', 'Fall_2022'='Sep-Nov 2022', 
                                       'Winter_2022'='Dec 2022')),
              strip = TMPstrip)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face ='bold', size=13),
    axis.title.y = element_text(face ='italic', size=15),
    axis.title.x = element_blank(),
    legend.text = element_text(face ='bold', size=13),
    legend.title = element_text(face ='italic', size=15),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # aspect.ratio = 1,
    strip.text = element_text(face='bold'),
    legend.position = 'none'
  )


Chla
Fig_BC_per_season

ggarrange(Chla, Fig_BC_per_season, nrow = 2 )



  
mean_Richness
mean_Richness_rare
g
h
View(rar_water@tax_table)
ggarrange(mean_Richness, mean_Richness_rare, g, h)

SAM <- rar_water@sam_data %>% as.data.frame()
chloro_year <- SAM %>%
  group_by(Cal_month) %>%
  summarise_at(c("Chla"), mean)
View(tmp)
ggplot(tmp, aes(x = Cal_month, y=Chla))+
  geom_smooth(method = 'loess')+
  scale_y_log10()

chloro_year

View(SAM)
  
  

SAM <- rar_water@sam_data %>% as.data.frame()
SAM$Temperature <- as.numeric(SAM$Temperature)
tmp <- SAM %>% group_by(cat, Temperature) %>%
  summarise(test = mean(Temperature))
View(tmp)
write.csv(tmp, 'tmp.csv')





Fig_MAG_size_barplot<-MAG_abund.df %>% melt(id = "dMAG_ID", value.name = "rel.abund") %>%
  rename("sample"="variable") %>%
  dplyr::mutate(.,
                phylum=MAG_summary.df$phylum[match(.$dMAG_ID,MAG_summary.df$dMAG_ID)],
                class=MAG_summary.df$class[match(.$dMAG_ID,MAG_summary.df$dMAG_ID)],
                genus=MAG_summary.df$genus[match(.$dMAG_ID,MAG_summary.df$dMAG_ID)],
                month_number=metadata_metaG$month_number[match(.$sample,metadata_metaG$sample)],
                lake_ID=metadata_metaG$lake_ID[match(.$sample,metadata_metaG$sample)],
                month_ID=metadata_metaG$month_ID[match(.$sample,metadata_metaG$sample)],
                genome_size=MAG_summary.df$genome_type[match(.$dMAG_ID,MAG_summary.df$dMAG_ID)],
                quality=MAG_summary.df$MAG_quality[match(.$dMAG_ID,MAG_summary.df$dMAG_ID)]) %>%
  subset(.,quality=="High_Quality") %>%
  group_by(lake_ID,month_ID,month_number,genome_size) %>%
  summarise(rel.abund=sum(rel.abund)) %>%
  group_by(lake_ID,month_ID,month_number) %>%
  dplyr::mutate(.,
                rel_abund=rel.abund/sum(rel.abund)*100,
                lake_ID=factor(lake_ID,levels=rev(c("VSS","CRE","LGP","CER-S","CSM","CER-L","JAB","VSM")))) %>%
  ggplot(.,aes(month_number,rel_abund,fill=genome_size))+
  geom_col(width=0.8,color="black",size=0.2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        legend.title = element_text(face="bold",size=12),
        legend.text = element_text(size=10),
        axis.title = element_blank(),
        axis.text.x = element_text(size=7.5,hjust=0.5),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank())+
  scale_fill_manual(values=palette_genomz_size)+
  guides(fill=guide_legend(nrow= 1,reverse = TRUE))+
  scale_x_continuous(n.breaks = 12)+
  scale_y_continuous(expand = c(0.01,0.01),breaks = c(0,25,50,75,100))+
  annotate("rect", xmin = 0.57, xmax = 2.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black") +
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -4, ymax = -0.5, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black") +
  annotate("rect", xmin = 11.5, xmax = 12.43, ymin = -4, ymax = -0.5, fill = "#63849B",color="black") +
  labs(fill="Genome size category")+
  facet_wrap2(~ lake_ID,scale="free_x",nrow = 2,
              strip = strip_color_lake<- strip_themed(
                background_x = elem_list_rect(fill = "black",
                                              color = "black"),
                text_x = elem_list_text(colour = "white",face = "bold",size=10)))+
  coord_cartesian(xlim =c(1,12), ylim = c(-4,NA),clip = "off")

("#CCEBC5","#308238","#67B16B","lightgrey","#008080","#173C5B","#225A88","#245D8E","#2E78B6","#97BCDB","#CBDDE")



# PieChart ----------------------------------------------------------------

library(phyloseq)
library(ggplot2)
library(dplyr)

# Chargement de l'objet phyloseq
# phyloseq_object <- load("votre_fichier_phyloseq.Rdata")

genres_heterotrophes <- c("Ciliophora",
                  "Kathablepharida","Fungi","Cercozoa",
                  "Bigyra","Choanoflagellata", "Chrompodellids")
genres_autotrophes <- c("Gyrista","Chlorophyta_X", "Haptophyta_X")
genres_parasites <- c("Perkinsea","Telonemia_X")
genres_mixotrophes <- c("Cryptophyta_X", "Chryptophyta_X:nucl", "Dinoflagellata")

tax_data <- as.data.frame(tax_table(rar_water))
sample_data_df <- as.data.frame((sample_data(rar_water)))
otu_data <- as.data.frame(otu_table(rar_water))
View(sam_data_df)

sample_data_df <- sample_data_df %>%
  rename(SampleID = rownames(sample_data_df))


tax_data$Mode_Trophique <- ifelse(tax_data$Order %in% genres_heterotrophes, "Hétérotrophe",
                                  ifelse(tax_data$Order %in% genres_autotrophes, "Autotrophe",
                                         ifelse(tax_data$Order %in% genres_parasites, 'Parasite',
                                                ifelse(tax_data$Order %in% genres_mixotrophes, 'Mixotrophe', 'NA'))))

otu_tax_data <- merge(otu_data, tax_data, by = "row.names")
View(otu_tax_data)
otu_tax_data <- merge(otu_tax_data, sam_data_df, by = "colnams")
richesse_par_lac <- otu_tax_data %>%
  group_by(Lake, Mode_Trophique) %>%
  summarise(Richesse = n_distinct(OTU))

richesse_par_lac %>%
  ggplot(aes(x = "", y = Richesse, fill = Mode_Trophique)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ Lac) +
  theme_void() +
  labs(title = "Richesse des genres autotrophes et hétérotrophes par lac",
       fill = "Mode Trophique")





write.csv("/Users/Amaury/OneDrive/Bureau/Rpertory/Plot/pie.csv")
write.csv(x = otu_tax_data, file = "/Users/Amaury/OneDrive/Bureau/Rpertory/Plot/Pie.csv")
pie <- read.csv(file = "/Users/Amaury/OneDrive/Bureau/Rpertory/pie.txt", header = TRUE)
View(pie)

View(rar_water@sam_data)



# Pierre ------------------------------------------------------------------

tmp <- as.data.frame(rar_water@sam_data)
tmpp <- tmp %>%
  group_by(lake, month) %>%
  summarise(mean_chla = mean(Chla, na.rm = TRUE))
# Calculer la différence max - min de chla par lac
diff_chla <- tmpp %>%
  group_by(lake) %>%
  summarise(chla_range = max(mean_chla, na.rm = TRUE) - min(mean_chla, na.rm = TRUE))


View(diff_chla)


TEST <- data.table(Lake = c('BLR', 'CSM', 'CER-S', 'CER-L', 'CRE', 'LGP', 'JAB', 'VSM', 'VSS'), MOTA = c(7.8622944,6.8877732,7.8858792,9.5464144,9.7833574,7.78393354,8.5721360,7.6973117,10.3744089),
                   RangeChla = c(50.36667,108.93333,14.9,12.5333,87.32333,64.26667,19.26667,17.16667,270.5333))


normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
normal_test <- TEST %>%
  mutate(
    mota_normalized = normalize(MOTA),
    Chla_normalized = normalize(RangeChla)
  )
View(normal_test)


RangeChlaMOTA <- normal_test %>%
  dplyr::mutate(., Lake=factor(Lake, levels=c("JAB", "VSM", "CER-L", 
                                              "CSM", "CER-S", "CRE",
                                              "LGP", "BLR", "VSS"))) %>%
ggplot(., aes(y = mota_normalized, x = Chla_normalized, color = Lake)) +
  geom_point(size = 4) +     
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=1)+
  labs(y = "normal_mota", x = "normal_ChlaRange")+
  scale_color_manual(values=palette_lake_3T)+
  theme(
    aspect.ratio = 1
  )

RangeChlaMOTA

###

RangeChlaMOTA <- TEST %>%
  dplyr::mutate(., Lake=factor(Lake, levels=c("JAB", "VSM", "CER-L", 
                                              "CSM", "CER-S", "CRE",
                                              "LGP", "BLR", "VSS"))) %>%
  ggplot(., aes(y = MOTA, x = RangeChla, color = Lake)) +
  theme_bw()+
  geom_point(size = 6)+
  theme(aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  #geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=1)+
  #labs(x = "longueur trajectoire", y = "Range de CHLA du lac", title = "Graphique RangeChla vs MOTA")+
  scale_color_manual(values=palette_lake_3T)+
  scale_x_log10()

RangeChlaMOTA


####________________________####