#### Import data ####

####__16S-SILVA (rarefied at 8,283 reads) ####

rar_prok <-
  create_phyloseq_SILVA("16S_data/SILVA_ASV_table.tsv",
                        "16S_data/SILVA_ASV-tax.tsv",
                        "16S_data/metadata_16S.txt") %>%
  subset_samples(sample_code %out% "VSM_B_W1")

#### | ####

#### Alpha-Diversity ####

####__ Rar. Curve ####
Fig_rar_curve_16S <-
  rarecurve(otu_table(rar_prok) %>% as.data.frame() %>% t(), step=500, cex=0.3,tidy = T) %>%
  ggplot(.,aes(x=Sample,y=Species,group=Site),color="black")+
  geom_line(linewidth = 0.2)+theme_bw()+
  labs( y = "Prokaryota ASVs Richness", x = "Nb. of reads")+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        axis.ticks = element_blank())+
  scale_x_continuous(expand=c(0,0),
                     limits = c(0,8284),
                     breaks = c(0,2000,4000,6000,8000))+
  scale_y_continuous(expand=c(0,0),
                     limits = c(0,600),
                     breaks = c(0,200,400,600))
Fig_rar_curve_16S

####__Indices ####
alphadiv_16S <- rar_prok@otu_table %>% estimate_richness() %>%
  select(c("Observed","Shannon")) %>% rename("Richness"="Observed") %>%
  dplyr::mutate(.,
                .after=1,Evenness=(Shannon/log(Richness))) %>%
  cbind(.,sample_code = rar_prok@sam_data$sample_code) %>%
  dplyr::mutate(.,
                lakeID=
                  (rar_prok@sam_data$lakeID[match(.$sample_code,rar_prok@sam_data$sample_code)]),
                day_number=
                  (rar_prok@sam_data$day_number[match(.$sample_code,rar_prok@sam_data$sample_code)]))
write.csv(alphadiv_16S,"output/alphadiv_16S.csv")

alphadiv_16S %>%
  group_by(lakeID) %>%
  summarise(Richness_mean=mean(Richness),
            Evenness_mean=mean(Evenness),
            Shannon_mean=mean(Shannon))

####____Richness ####
richness_16S<- alphadiv_16S %>%
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 500),
                     breaks = c(0, 200, 400))+
  scale_x_continuous(expand = c(0,0), limits = c(-1, 545),
                     breaks = c(0, 30, 60, 120, 180, 360, 540))+
  scale_color_manual(values = rev(palette_lake_chla))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  labs(y = 'Prokaryota ASV Richness')
richness_16S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/richness_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####____Evenness ####
evenness_16S<- alphadiv_16S %>%
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
  labs(y = 'Prokaryota ASV Evenness')
evenness_16S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/evenness_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####____Shannon ####
shannon_16S<- alphadiv_16S %>%
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
  labs(y = 'Prokaryota ASV Shannon')
shannon_16S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/shannon_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### | ####

#### Composition ####

####__season-Sub. ####
barplot_phylum_season.df <-
  rar_prok %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_season,Phylum) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_season) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_prok@sam_data$lakeID[match(.$lake_season,rar_prok@sam_data$lake_season)]),
        season_year=
          (rar_prok@sam_data$season_year[match(.$lake_season,rar_prok@sam_data$lake_season)]),
        .)

palette_Phylum <- c("#B3CDE3","#FED9A6","#DECBE4","#CCEBC5","#E5D8BD",
                    "#FFFFCC","#862C56","gray38","#fddaec","#99B093",
                    "#FC8F36","#C03E7B","black","#357227","#F2F2F2",'#ffb3b3')



barplot_phylum_season <- barplot_phylum_season.df %>%
  dplyr::mutate(.,
                taxa_legend=ifelse(median_abundance<1,"Taxa < 1%",Phylum)) %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS")),
                season_year = factor(season_year, levels = c("summer_2021","fall_2021","winter_2021",
                                                             "spring_2022","summer_2022","fall_2022","winter_2022"))) %>%
  ggplot(.,aes(x=season_year,y=median_abundance,
               fill=fct_rev(fct_reorder(taxa_legend,median_abundance))))+
  theme_bw()+
  geom_col(width = .95, color = "black", linewidth = 0.3)+
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
  scale_fill_manual(values=c("#B3CDE3","#FED9A6","#DECBE4","#CCEBC5","#E5D8BD",
                             "#FFFFCC","#862C56","gray38","#fddaec","#99B093",
                             "#FC8F36","#C03E7B","black","#357227","#F2F2F2",'#ffb3b3'))+ #set a color palette for the barplot
  labs(y = "", fill="Prokaryota\nPhyla", x = "")+ #change your axis and legend title 
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
barplot_phylum_season

ggsave("/Users/piefouca/Desktop/µEuk/Figures/prok_phylum_season_barplot.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__month-Sub. ####
barplot_phylum.df <- rar_prok %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  group_by(lake_month,Phylum) %>%
  summarise(total_count=sum(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(total_count/sum(total_count),4))*100))) %>%
  cbind(lakeID=
          (rar_prok@sam_data$lakeID[match(.$lake_month,rar_prok@sam_data$lake_month)]),
        month_code=
          (rar_prok@sam_data$month_code[match(.$lake_month,rar_prok@sam_data$lake_month)]),
        .)

barplot_phylum <- barplot_phylum.df %>%
  dplyr::mutate(taxa_legend=ifelse(median_abundance<1,"Taxa < 1%",Phylum),
                lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS")),
                month_code = recode(month_code,
                                    'A'=' 6','B'=' 7','C'=' 8','D'=' 9',
                                    'E'=' 10','F'=' 11','G'='1','H'='2',
                                    'I'='3','J'='4','K'='5','L'='6','M'='7',
                                    'N'='8','O'='9','P'='10','Q'='11','R'='12'),
                month_code=factor(month_code, levels=c(' 6',' 7',' 8',' 9',' 10',' 11','1',
                                                       '2', '3','4','5','6','7','8','9','10',
                                                       '11', '12'))) %>%
  ggplot(.,aes(x=month_code,y=median_abundance,
               fill=fct_rev(fct_reorder(taxa_legend,median_abundance))))+
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
  scale_fill_manual(values=c("#B3CDE3","#FED9A6","#CCEBC5","#DECBE4","#E5D8BD",
                             "#FFFFCC","#357227","#C03E7B","#F2F2F2","#fddaec",
                             "#FC8F36","black","#99B093","gray38","#862C56",
                             "gray20","#f6b26b","blue",'#ffb3b3'))+#set a color palette for the barplot
  labs(y = "", fill="Prokaryota\nPhyla", x = "")+ #change your axis and legend title 
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
barplot_phylum

ggsave("/Users/piefouca/Desktop/µEuk/Figures/prok_phylum_month_barplot.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### | ####

#### Beta-Diversity ####

BC_16S.dist <- vegdist(t(rar_prok@otu_table),method = "bray")
JC_16S.dist <- vegdist(t(rar_prok@otu_table),
                       method = "jaccard",binary = T)
mantel_stat_16S <-
  mantel(BC_16S.dist,JC_16S.dist,
         method = "spearman",permutations = 999)$statistic

####__Bray-Curtis ####
BC_16S.dist<-ordinate(rar_prok,"PCoA","bray")

BC_16S_out<- plot_ordination(rar_prok,BC_16S.dist,label="season_year",color="lakeID")

BC_rar_prok.df<- as.data.frame(BC_16S_out$data)

BC_16S <- BC_rar_prok.df %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS")),
                season_year=factor(season_year,levels = c("summer_2021","fall_2021","winter_2021",
                                                          "spring_2022","summer_2022","fall_2022","winter_2022")),
                year=factor(year,levels = c("2021", "2022"))) %>%
  ggplot(.,aes(Axis.1,Axis.2,color=season_year,fill=season_year, shape = year, group = season_year))+
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
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.38,- -0.62),
                     breaks = c(-0.4, -0.2, 0, 0.2, 0.4),
                     labels = c(-0.4, -0.2, 0, 0.2, 0.4))+
  scale_x_continuous(expand = c(0,0), trans = "reverse",
                     limits = c(0.25,-0.62),
                     breaks = c(0.2, 0, -0.2, -0.4, -0.6),
                     labels = c(0.2, 0, -0.2, -0.4, -0.6))+
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
       x=paste0("Axis 1 [",round(BC_16S.dist$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 1 [",round(BC_16S.dist$values$Relative_eig[2],3)*100,"%]"),
       title=paste0("Prokaryotic Community · Bray-Curtis"))
BC_16S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/BC_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### Time Lag Analysis ####

####__BC by day_gap ####
TLA_prok.df<- as.matrix(phyloseq::distance(rar_prok,method = "bray")) %>%
  melt(.) %>% rename(G1=Var1,G2=Var2) %>%
  filter(as.character(G1) != as.character(G2)) %>%
  mutate_if(is.factor,as.character) %>%
  dplyr::rename("BC"="value") %>%
  cbind(.,
        G1_day=rar_prok@sam_data$day_number[match(.$G1,
                                                  rownames(rar_prok@sam_data))],
        G2_day=rar_prok@sam_data$day_number[match(.$G2,
                                                  rownames(rar_prok@sam_data))],
        G1_lakeID=rar_prok@sam_data$lakeID[match(.$G1,
                                                 rownames(rar_prok@sam_data))],
        G2_lakeID=rar_prok@sam_data$lakeID[match(.$G2,
                                                 rownames(rar_prok@sam_data))],
        G1_column=rar_prok@sam_data$column[match(.$G1,
                                                 rownames(rar_prok@sam_data))],
        G2_column=rar_prok@sam_data$column[match(.$G2,
                                                 rownames(rar_prok@sam_data))]) %>%
  dplyr::mutate(day_gap=abs(.$G1_day-.$G2_day),
                comp_type=paste0(G1_lakeID,"_",G1_lakeID),
                comp_col=paste0(G1_column,"_",G2_column)) %>%
  subset(.,G1_day-G2_day >0 & G1_lakeID==G2_lakeID & G1_column==G2_column)

####__linear vs polynomial ####
list_lake<- unique(TLA_prok.df$G1_lakeID)

comp_model<-c("lakeID","sig_linear","r_linear",
              "AIC_linear","AIC_diff_poly2","AIC_diff_poly3","AIC_diff_poly4","AIC_diff_poly5",
              "sig_poly2","r_poly2","sig_poly3","r_poly3","sig_poly4","r_poly4","sig_poly5","r_poly5")

for (i in 1:length(unique(TLA_prok.df$G1_lakeID))) {
  # select lake
  lm_lake<-TLA_prok.df %>% subset(G1_lakeID %in% list_lake[i])
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

write.csv(comp_model,"output/model_comp_16S.csv",row.names = F)

####__plot ####
TLA_prok.df <- TLA_prok.df %>%
  dplyr::mutate(.,
                poly=ifelse(G1_lakeID %in% c("BLR", "CSM","CRE"), "poly5",
                            ifelse(G1_lakeID %in% c("CER-S","JAB"), "poly3",
                                   ifelse(G1_lakeID %in% c("CER-L","LGP","VSM","VSS"),"poly4","linear"))),
                G1_lakeID=factor(G1_lakeID,
                                 levels=c("VSM", "JAB", "CER-L",
                                          "CER-S", "CRE", "BLR","LGP", "CSM", 
                                          "VSS")))
TLA_prok <- TLA_prok.df %>%
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
            data =TLA_prok.df[(TLA_prok.df$poly=="poly3"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly3,show.legend = F,
                 data =TLA_prok.df[(TLA_prok.df$poly=="poly3"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_prok.df[(TLA_prok.df$poly=="poly3"),])+
  #best poly model poly4
  geom_text(x=250,y=0.132,size=3,check_overlap = T,fontface="italic",color = "blue",
            data =TLA_prok.df[(TLA_prok.df$poly=="poly4"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly4,show.legend = F,
                 data =TLA_prok.df[(TLA_prok.df$poly=="poly4"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly3,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_prok.df[(TLA_prok.df$poly=="poly4"),])+
  #best poly model poly5
  geom_text(x=250,y=0.132,size=3,check_overlap = T,fontface="italic",color = "blue",
            data =TLA_prok.df[(TLA_prok.df$poly=="poly5"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly5,show.legend = F,
                 data =TLA_prok.df[(TLA_prok.df$poly=="poly5"),])+
  stat_poly_eq(use_label(c("R2","P")),size=3,formula = lm_formula_poly5,show.legend = F,
               label.y = 0.01, label.x = 0.97,rr.digits =2,color = "blue",
               data =TLA_prok.df[(TLA_prok.df$poly=="poly5"),])+
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
  labs(y="Prokaryotic Bray-curtis dissimilarity",
       x="Time between sampling dates (days)")+
  facet_wrap2(~ G1_lakeID,nrow = 3,scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))
TLA_prok

ggsave("/Users/piefouca/Desktop/µEuk/Figures/TLA_prok.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### MOTA ####
#61 Axis are needed to explain 90.07 %

MOTA_lake_16S <- final_cum_distance_16S %>%
  melt(.,"time_label") %>% as.data.frame() %>%
  dplyr::rename("lakeID"="variable") %>%
  dplyr::mutate(lake_month=paste0(lakeID,"_",time_label),
                season_year = rar_prok@sam_data$season_year[match(.$time_label,
                                                                  rar_prok@sam_data$month_number)],
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
           xmax=(final_cum_distance_16S$VSS[2]+((final_cum_distance_16S$VSS[3]-final_cum_distance_16S$VSS[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[2]+((final_cum_distance_16S$VSS[3]-final_cum_distance_16S$VSS[2])/2)),
           xmax=(final_cum_distance_16S$VSS[5]+((final_cum_distance_16S$VSS[6]-final_cum_distance_16S$VSS[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[5]+((final_cum_distance_16S$VSS[6]-final_cum_distance_16S$VSS[5])/2)),
           xmax=(final_cum_distance_16S$VSS[7]+((final_cum_distance_16S$VSS[8]-final_cum_distance_16S$VSS[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[7]+((final_cum_distance_16S$VSS[8]-final_cum_distance_16S$VSS[7])/2)),
           xmax=(final_cum_distance_16S$VSS[10]+((final_cum_distance_16S$VSS[11]-final_cum_distance_16S$VSS[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[10]+((final_cum_distance_16S$VSS[11]-final_cum_distance_16S$VSS[10])/2)),
           xmax=(final_cum_distance_16S$VSS[13]+((final_cum_distance_16S$VSS[14]-final_cum_distance_16S$VSS[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[13]+((final_cum_distance_16S$VSS[14]-final_cum_distance_16S$VSS[13])/2)),
           xmax=(final_cum_distance_16S$VSS[16]+((final_cum_distance_16S$VSS[17]-final_cum_distance_16S$VSS[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=0.84,ymax=1.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSS[16]+((final_cum_distance_16S$VSS[17]-final_cum_distance_16S$VSS[16])/2)),
           xmax=max(final_cum_distance_16S$VSS)+0.18)+
  #color CSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$CSM[2]+((final_cum_distance_16S$CSM[3]-final_cum_distance_16S$CSM[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[2]+((final_cum_distance_16S$CSM[3]-final_cum_distance_16S$CSM[2])/2)),
           xmax=(final_cum_distance_16S$CSM[5]+((final_cum_distance_16S$CSM[6]-final_cum_distance_16S$CSM[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[5]+((final_cum_distance_16S$CSM[6]-final_cum_distance_16S$CSM[5])/2)),
           xmax=(final_cum_distance_16S$CSM[7]+((final_cum_distance_16S$CSM[8]-final_cum_distance_16S$CSM[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[7]+((final_cum_distance_16S$CSM[8]-final_cum_distance_16S$CSM[7])/2)),
           xmax=(final_cum_distance_16S$CSM[10]+((final_cum_distance_16S$CSM[11]-final_cum_distance_16S$CSM[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[10]+((final_cum_distance_16S$CSM[11]-final_cum_distance_16S$CSM[10])/2)),
           xmax=(final_cum_distance_16S$CSM[13]+((final_cum_distance_16S$CSM[14]-final_cum_distance_16S$CSM[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[13]+((final_cum_distance_16S$CSM[14]-final_cum_distance_16S$CSM[13])/2)),
           xmax=(final_cum_distance_16S$CSM[16]+((final_cum_distance_16S$CSM[17]-final_cum_distance_16S$CSM[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=1.84,ymax=2.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CSM[16]+((final_cum_distance_16S$CSM[17]-final_cum_distance_16S$CSM[16])/2)),
           xmax=max(final_cum_distance_16S$CSM)+0.18)+
  #color LGP
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$LGP[2]+((final_cum_distance_16S$LGP[3]-final_cum_distance_16S$LGP[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[2]+((final_cum_distance_16S$LGP[3]-final_cum_distance_16S$LGP[2])/2)),
           xmax=(final_cum_distance_16S$LGP[5]+((final_cum_distance_16S$LGP[6]-final_cum_distance_16S$LGP[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[5]+((final_cum_distance_16S$LGP[6]-final_cum_distance_16S$LGP[5])/2)),
           xmax=(final_cum_distance_16S$LGP[7]+((final_cum_distance_16S$LGP[8]-final_cum_distance_16S$LGP[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[7]+((final_cum_distance_16S$LGP[8]-final_cum_distance_16S$LGP[7])/2)),
           xmax=(final_cum_distance_16S$LGP[10]+((final_cum_distance_16S$LGP[11]-final_cum_distance_16S$LGP[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[10]+((final_cum_distance_16S$LGP[11]-final_cum_distance_16S$LGP[10])/2)),
           xmax=(final_cum_distance_16S$LGP[13]+((final_cum_distance_16S$LGP[14]-final_cum_distance_16S$LGP[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[13]+((final_cum_distance_16S$LGP[14]-final_cum_distance_16S$LGP[13])/2)),
           xmax=(final_cum_distance_16S$LGP[16]+((final_cum_distance_16S$LGP[17]-final_cum_distance_16S$LGP[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=2.84,ymax=3.16,alpha=0.9,
           xmin=(final_cum_distance_16S$LGP[16]+((final_cum_distance_16S$LGP[17]-final_cum_distance_16S$LGP[16])/2)),
           xmax=max(final_cum_distance_16S$LGP)+0.18)+
  #color BLR
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$BLR[2]+((final_cum_distance_16S$BLR[3]-final_cum_distance_16S$BLR[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[2]+((final_cum_distance_16S$BLR[3]-final_cum_distance_16S$BLR[2])/2)),
           xmax=(final_cum_distance_16S$BLR[5]+((final_cum_distance_16S$BLR[6]-final_cum_distance_16S$BLR[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[5]+((final_cum_distance_16S$BLR[6]-final_cum_distance_16S$BLR[5])/2)),
           xmax=(final_cum_distance_16S$BLR[7]+((final_cum_distance_16S$BLR[8]-final_cum_distance_16S$BLR[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[7]+((final_cum_distance_16S$BLR[8]-final_cum_distance_16S$BLR[7])/2)),
           xmax=(final_cum_distance_16S$BLR[10]+((final_cum_distance_16S$BLR[11]-final_cum_distance_16S$BLR[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[10]+((final_cum_distance_16S$BLR[11]-final_cum_distance_16S$BLR[10])/2)),
           xmax=(final_cum_distance_16S$BLR[13]+((final_cum_distance_16S$BLR[14]-final_cum_distance_16S$BLR[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[13]+((final_cum_distance_16S$BLR[14]-final_cum_distance_16S$BLR[13])/2)),
           xmax=(final_cum_distance_16S$BLR[16]+((final_cum_distance_16S$BLR[17]-final_cum_distance_16S$BLR[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=3.84,ymax=4.16,alpha=0.9,
           xmin=(final_cum_distance_16S$BLR[16]+((final_cum_distance_16S$BLR[17]-final_cum_distance_16S$BLR[16])/2)),
           xmax=max(final_cum_distance_16S$BLR)+0.18)+
  #color CRE
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$CRE[2]+((final_cum_distance_16S$CRE[3]-final_cum_distance_16S$CRE[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[2]+((final_cum_distance_16S$CRE[3]-final_cum_distance_16S$CRE[2])/2)),
           xmax=(final_cum_distance_16S$CRE[5]+((final_cum_distance_16S$CRE[6]-final_cum_distance_16S$CRE[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[5]+((final_cum_distance_16S$CRE[6]-final_cum_distance_16S$CRE[5])/2)),
           xmax=(final_cum_distance_16S$CRE[7]+((final_cum_distance_16S$CRE[8]-final_cum_distance_16S$CRE[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[7]+((final_cum_distance_16S$CRE[8]-final_cum_distance_16S$CRE[7])/2)),
           xmax=(final_cum_distance_16S$CRE[10]+((final_cum_distance_16S$CRE[11]-final_cum_distance_16S$CRE[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[10]+((final_cum_distance_16S$CRE[11]-final_cum_distance_16S$CRE[10])/2)),
           xmax=(final_cum_distance_16S$CRE[13]+((final_cum_distance_16S$CRE[14]-final_cum_distance_16S$CRE[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[13]+((final_cum_distance_16S$CRE[14]-final_cum_distance_16S$CRE[13])/2)),
           xmax=(final_cum_distance_16S$CRE[16]+((final_cum_distance_16S$CRE[17]-final_cum_distance_16S$CRE[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=4.84,ymax=5.16,alpha=0.9,
           xmin=(final_cum_distance_16S$CRE[16]+((final_cum_distance_16S$CRE[17]-final_cum_distance_16S$CRE[16])/2)),
           xmax=max(final_cum_distance_16S$CRE)+0.18)+
  #color `CER-S`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$`CER-S`[2]+((final_cum_distance_16S$`CER-S`[3]-final_cum_distance_16S$`CER-S`[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[2]+((final_cum_distance_16S$`CER-S`[3]-final_cum_distance_16S$`CER-S`[2])/2)),
           xmax=(final_cum_distance_16S$`CER-S`[5]+((final_cum_distance_16S$`CER-S`[6]-final_cum_distance_16S$`CER-S`[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[5]+((final_cum_distance_16S$`CER-S`[6]-final_cum_distance_16S$`CER-S`[5])/2)),
           xmax=(final_cum_distance_16S$`CER-S`[7]+((final_cum_distance_16S$`CER-S`[8]-final_cum_distance_16S$`CER-S`[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[7]+((final_cum_distance_16S$`CER-S`[8]-final_cum_distance_16S$`CER-S`[7])/2)),
           xmax=(final_cum_distance_16S$`CER-S`[10]+((final_cum_distance_16S$`CER-S`[11]-final_cum_distance_16S$`CER-S`[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[10]+((final_cum_distance_16S$`CER-S`[11]-final_cum_distance_16S$`CER-S`[10])/2)),
           xmax=(final_cum_distance_16S$`CER-S`[13]+((final_cum_distance_16S$`CER-S`[14]-final_cum_distance_16S$`CER-S`[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[13]+((final_cum_distance_16S$`CER-S`[14]-final_cum_distance_16S$`CER-S`[13])/2)),
           xmax=(final_cum_distance_16S$`CER-S`[16]+((final_cum_distance_16S$`CER-S`[17]-final_cum_distance_16S$`CER-S`[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=5.84,ymax=6.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-S`[16]+((final_cum_distance_16S$`CER-S`[17]-final_cum_distance_16S$`CER-S`[16])/2)),
           xmax=max(final_cum_distance_16S$`CER-S`)+0.18)+
  #color `CER-L`
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$`CER-L`[2]+((final_cum_distance_16S$`CER-L`[3]-final_cum_distance_16S$`CER-L`[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[2]+((final_cum_distance_16S$`CER-L`[3]-final_cum_distance_16S$`CER-L`[2])/2)),
           xmax=(final_cum_distance_16S$`CER-L`[5]+((final_cum_distance_16S$`CER-L`[6]-final_cum_distance_16S$`CER-L`[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[5]+((final_cum_distance_16S$`CER-L`[6]-final_cum_distance_16S$`CER-L`[5])/2)),
           xmax=(final_cum_distance_16S$`CER-L`[7]+((final_cum_distance_16S$`CER-L`[8]-final_cum_distance_16S$`CER-L`[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[7]+((final_cum_distance_16S$`CER-L`[8]-final_cum_distance_16S$`CER-L`[7])/2)),
           xmax=(final_cum_distance_16S$`CER-L`[10]+((final_cum_distance_16S$`CER-L`[11]-final_cum_distance_16S$`CER-L`[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[10]+((final_cum_distance_16S$`CER-L`[11]-final_cum_distance_16S$`CER-L`[10])/2)),
           xmax=(final_cum_distance_16S$`CER-L`[13]+((final_cum_distance_16S$`CER-L`[14]-final_cum_distance_16S$`CER-L`[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[13]+((final_cum_distance_16S$`CER-L`[14]-final_cum_distance_16S$`CER-L`[13])/2)),
           xmax=(final_cum_distance_16S$`CER-L`[16]+((final_cum_distance_16S$`CER-L`[17]-final_cum_distance_16S$`CER-L`[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=6.84,ymax=7.16,alpha=0.9,
           xmin=(final_cum_distance_16S$`CER-L`[16]+((final_cum_distance_16S$`CER-L`[17]-final_cum_distance_16S$`CER-L`[16])/2)),
           xmax=max(final_cum_distance_16S$`CER-L`)+0.18)+
  #color JAB
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$JAB[2]+((final_cum_distance_16S$JAB[3]-final_cum_distance_16S$JAB[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[2]+((final_cum_distance_16S$JAB[3]-final_cum_distance_16S$JAB[2])/2)),
           xmax=(final_cum_distance_16S$JAB[5]+((final_cum_distance_16S$JAB[6]-final_cum_distance_16S$JAB[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[5]+((final_cum_distance_16S$JAB[6]-final_cum_distance_16S$JAB[5])/2)),
           xmax=(final_cum_distance_16S$JAB[7]+((final_cum_distance_16S$JAB[8]-final_cum_distance_16S$JAB[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[7]+((final_cum_distance_16S$JAB[8]-final_cum_distance_16S$JAB[7])/2)),
           xmax=(final_cum_distance_16S$JAB[10]+((final_cum_distance_16S$JAB[11]-final_cum_distance_16S$JAB[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[10]+((final_cum_distance_16S$JAB[11]-final_cum_distance_16S$JAB[10])/2)),
           xmax=(final_cum_distance_16S$JAB[13]+((final_cum_distance_16S$JAB[14]-final_cum_distance_16S$JAB[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[13]+((final_cum_distance_16S$JAB[14]-final_cum_distance_16S$JAB[13])/2)),
           xmax=(final_cum_distance_16S$JAB[16]+((final_cum_distance_16S$JAB[17]-final_cum_distance_16S$JAB[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=7.84,ymax=8.16,alpha=0.9,
           xmin=(final_cum_distance_16S$JAB[16]+((final_cum_distance_16S$JAB[17]-final_cum_distance_16S$JAB[16])/2)),
           xmax=max(final_cum_distance_16S$JAB)+0.18)+
  #color VSM
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=-0.11,
           xmax=(final_cum_distance_16S$VSM[2]+((final_cum_distance_16S$VSM[3]-final_cum_distance_16S$VSM[2])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[2]+((final_cum_distance_16S$VSM[3]-final_cum_distance_16S$VSM[2])/2)),
           xmax=(final_cum_distance_16S$VSM[5]+((final_cum_distance_16S$VSM[6]-final_cum_distance_16S$VSM[5])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[5]+((final_cum_distance_16S$VSM[6]-final_cum_distance_16S$VSM[5])/2)),
           xmax=(final_cum_distance_16S$VSM[7]+((final_cum_distance_16S$VSM[8]-final_cum_distance_16S$VSM[7])/2)))+
  annotate("rect",color="#EBAF47",fill="#EBAF47",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[7]+((final_cum_distance_16S$VSM[8]-final_cum_distance_16S$VSM[7])/2)),
           xmax=(final_cum_distance_16S$VSM[10]+((final_cum_distance_16S$VSM[11]-final_cum_distance_16S$VSM[10])/2)))+
  annotate("rect",color="#80A53F",fill="#80A53F",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[10]+((final_cum_distance_16S$VSM[11]-final_cum_distance_16S$VSM[10])/2)),
           xmax=(final_cum_distance_16S$VSM[13]+((final_cum_distance_16S$VSM[14]-final_cum_distance_16S$VSM[13])/2)))+
  annotate("rect",color="#E36414",fill="#E36414",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[13]+((final_cum_distance_16S$VSM[14]-final_cum_distance_16S$VSM[13])/2)),
           xmax=(final_cum_distance_16S$VSM[16]+((final_cum_distance_16S$VSM[17]-final_cum_distance_16S$VSM[16])/2)))+
  annotate("rect",color="#63849B",fill="#63849B",size=0.5,ymin=8.84,ymax=9.16,alpha=0.9,
           xmin=(final_cum_distance_16S$VSM[16]+((final_cum_distance_16S$VSM[17]-final_cum_distance_16S$VSM[16])/2)),
           xmax=max(final_cum_distance_16S$VSM)+0.18)+
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
  scale_x_continuous(expand = c(0,0),limits = c(-0.25,10))+
  labs(#caption=paste0("MOTA – Bray-Curtis – ",Threshold*100,"% of the total expl. var. (58 axis)"),
    x="Prokaryotic community trajectories",fill="Season")+
  guides(fill=guide_legend(override.aes=list(shape=21,size=8,color="black"),nrow=2))
MOTA_lake_16S

ggsave("/Users/piefouca/Desktop/µEuk/Figures/MOTA_lake_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####________________________####
#### Old Pierre ####
####________________________####

#### turnover #####
## melt dataframe

prok_codyn.df <- rar_water_prok %>%
  psmelt(.) %>%
  dplyr::group_by(lake_month,OTU,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>% dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(paste0((round(median_count/sum(median_count),4))*100))) %>%
  cbind(lake_ID=
          (rar_water_prok@sam_data$lake_ID[match(.$lake_month,rar_water_2122@sam_data$lake_month)]),
        month_full=
          (rar_water_2122@sam_data$month[match(.$lake_month,rar_water_2122@sam_data$lake_month)])) %>%
  dplyr::mutate(.,sample_codyn=paste0(lake,"_",month_full,"_",ASV_codyn))

ASV_codyn<-as.data.frame(cbind(
  sample_codyn=Water_2122_codyn.df$sample_codyn,
  median_abund=Water_2122_codyn.df$median_abundance)) %>%
  dplyr::mutate(.,median_abund=as.numeric(median_abund)) %>%
  dplyr::group_by(sample_codyn) %>%
  dplyr::summarise(median_rel_abund=sum(median_abund)) %>%
  cbind(lake=
          Water_2122_codyn.df$lake[match(.$sample_codyn,Water_2122_codyn.df$sample_codyn)],
        month_full=
          Water_2122_codyn.df$month_full[match(.$sample_codyn,Water_2122_codyn.df$sample_codyn)],
        ASV_codyn=
          Water_2122_codyn.df$ASV_codyn[match(.$sample_codyn,Water_2122_codyn.df$sample_codyn)],
        median_rel_abund=.$median_rel_abund) %>% .[,-c(1,2)] %>%
  dplyr::mutate(.,
                month_number=rar_water_2122@sam_data$month_numer[match(.$month_full,rar_water_2122@sam_data$month)])

ASV_codyn$month_number=as.numeric(ASV_codyn$month_number)

# codyn's turnover function
ASV_turnover<-turnover(ASV_codyn,time.var = "month_number",species.var = "ASV_codyn",
                       abundance.var = "median_rel_abund",replicate.var = "lake",metric = "total") %>%
  dplyr::rename("ASV_turnover"="total") %>%
  dplyr::mutate(.,
                lake_month=paste0(lake,"_",month_number),
                lake=factor(.$lake,levels=c("CSM","VSS","LGP","BLR","CRE","CER-S","CER-L","VSM","JAB")))

# ggplot
Fig_ASV_turnover<-
  ASV_turnover %>%
  dplyr::mutate(.,
                lake=factor(lake,levels=rev(c("VSS","CSM","LGP","BLR","CRE","CER-S","CER-L","VSM","JAB"))),
                month=rar_water_2122@sam_data$month[match(.$lake_month,
                                                          rar_water_2122@sam_data$lake_ID_month)]) %>%
  ggplot(.,aes(x=month,y=ASV_turnover,group=1))+
  geom_line(color="black")+
  geom_point(size=1,show.legend = F)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        axis.title = element_blank())+
  guides(fill=guide_legend(nrow=3,shape=22,size=4),strip=T)+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))+
  facet_wrap2(~ lake,scale="fixed",nrow = 3,
              strip = strip_color_lake<- strip_themed(
                background_x = elem_list_rect(fill = rev(palette_lake_3T),
                                              color = "black"),
                text_x = elem_list_text(colour = "white",face = "bold",size=9)))
Fig_ASV_turnover

turnover_2122.df<-Fig_ASV_turnover$data

write.csv(turnover_2122.df,"turnover_2122.csv")

#### | ####


####________________________####