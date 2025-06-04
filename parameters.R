#### Chla ####

chla.df<-read_delim("param_data/chla_all_values.txt","\t",show_col_types = F) %>%
  filter(!is.na(Chla)) %>%
  subset(sample_code %out% "VSM_B_W1")

chla.df %>%
  group_by(lakeID) %>%
  summarise(chla_min=min(Chla),
            chla_max=max(Chla))

chla.df %>%
  group_by(lakeID) %>%
  summarise(chla_mean=mean(Chla))

chla_timeseries<- chla.df %>%
  group_by(lakeID,month_code) %>%
  summarise(chla_median=median(Chla)) %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  ggplot(.,aes(month_code,chla_median, group=lakeID))+
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 0.1, ymax = 2.6),
            fill = "#97B5CD",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 2.6, ymax = 7.3),
            fill = "#2F6B9D",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 7.3, ymax = 56),
            fill = "#94D2A4",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 56, ymax = 400),
            fill = "#28A448",color=NA,alpha=0.4) +
  geom_line(linewidth = 0.5)+
  geom_point(size = 2)+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_blank(),
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
  scale_y_continuous(breaks = c(2.6,7.3,56,400),limits = c(0.1,400),
                     labels = c(2.6,7.3,56,400),
                     expand = c(0,0),trans = "log")+
                     #breaks = c('2.6','7.4','56','100','400'))+ #keep a little extra space on the y axis #set a color palette for the barplot
  #scale_fill_manual(values=palette_Subdivision)+ #set a color palette for the barplot
 # labs(y = "", fill="Eukarya\nSubdivision", x = "")+ #change your axis and legend title 
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  labs(y=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))))+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.1, ymax = 0.18, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0.1, ymax = 0.18, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 0.1, ymax = 0.18, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = 0.1, ymax = 0.18, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = 0.1, ymax = 0.18, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = 0.1, ymax = 0.18, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = 0.1, ymax = 0.18, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(0.1,NA),clip = "off")
chla_timeseries

ggsave("/Users/piefouca/Desktop/µEuk/Figures/chla_timeseries.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### Phytoplanckton ####

phyto_data.df<-read.csv("param_data/phyto_count.tsv", sep = "\t") %>% 
  subset(column == "W2") %>%
  merge(., read.csv("param_data/genera_biovolumes.tsv", sep = "\t"), by = "Taxa") %>%
  dplyr::mutate(.,rel_biovolume=nb*Biovolume)

phyto_domain.df<-phyto_data.df %>% .[,c(2,3,10,11,14)] %>%
  drop_na(rel_biovolume) %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                                "CER-S", "CRE", "BLR","LGP", "CSM", 
                                                "VSS")),
                domain=if_else(Phylum=="Cyanobacteria","Prokaryota (Cyanobacteria)","Eukaryota"),
                domain=factor(domain,levels=c("Eukaryota","Prokaryota (Cyanobacteria)"))) %>%
  group_by(lakeID,month_code,domain) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,rel_biovolume=(rel_biovolume/sum(rel_biovolume))*100)

phyto_domain_barplot<- phyto_domain.df %>%
  # dplyr::mutate(.,
  #               month_code = recode(month_code,
  #                                   'A'=' 6','B'=' 7','C'=' 8','D'=' 9',
  #                                   'E'=' 10','F'=' 11','G'='1','H'='2',
  #                                   'I'='3','J'='4','K'='5','L'='6','M'='7',
  #                                   'N'='8','O'='9','P'='10','Q'='11','R'='12'),
  #               month_code=factor(month_code, levels=c(' 6',' 7',' 8',' 9',' 10',' 11','1',
  #                                                      '2', '3','4','5','6','7','8','9','10',
  #                                                      '11', '12'))) %>%
  ggplot(., aes(fill=domain, y=rel_biovolume, x=month_code))+
  geom_col(width = 0.8,color="black")+theme_bw()+
  theme(axis.title = element_blank(),
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y =element_text(color = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  guides(fill=guide_legend(nrow = 1,size=2))+
  scale_fill_manual(values = c("#308238","#CCEBC5"))+
  #scale_x_continuou(n.breaks = 12)+
  scale_x_discrete(expand = c(0,0),
                   labels = c('J','J','A','S','O','N',
                              'J','F', 'M','A','M','J',
                              'J','A','S','O','N', 'D'))+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0,25,50,75,100),
                     labels = c("0%","25%","50%","75%","100%"))+
  labs(fill="Phytoplankton domain")+
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -4, ymax = -0.5, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = -4, ymax = -0.5, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = -4, ymax = -0.5, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = -4, ymax = -0.5, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(-4,NA),clip = "off")
phyto_domain_barplot

ggsave("/Users/piefouca/Desktop/µEuk/Figures/phyto_domain_barplot.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

cyano_biovolume.df <- phyto_domain.df %>%
  subset(domain == "Prokaryota (Cyanobacteria)") %>%
  dplyr::mutate(.,lake_month=paste0(lakeID,"_",month_code))

phyto_domain.df %>%
  drop_na(TotalBiovolume) %>%
  group_by(lake_ID,month_number,domain) %>%
  summarise(TotalBiovolume=sum(TotalBiovolume)) %>%
  group_by(lake_ID,month_number) %>%
  dplyr::mutate(.,biovol=TotalBiovolume/sum(TotalBiovolume)*100) %>%
  group_by(lake_ID,domain) %>%
  dplyr::summarise_at(vars("biovol"),list(mean=mean, sd=sd,max=max,min=min))

####________________________####
#### Old Pierre ####
####________________________####

read_delim("chla_2022.csv",delim = ";",show_col_types = F) %>%
  subset(lake_ID!="BLR") %>%
  #subset(lake_ID %in% c("VSM","VSS")) %>%
  group_by(lake_ID,month_number) %>%
  drop_na(chla) %>%
  dplyr::summarise(chla_mean = mean(chla)) %>%
  dplyr::mutate(.,
                trophic_status=ifelse(lake_ID %in% c("VSS","CRE","LGP","CER-S"),
                                      "eutrophic","oligo-mesotrophic"),
                show_lake=ifelse(lake_ID %in% c("VSM","VSS"),"YES","NO"),
                lake_month=paste0(lake_ID,"_",month_number),
                #lake_ID=factor(lake_ID,levels=rev(c("VSS","CRE","LGP","CER-S","CSM","CER-L","JAB","VSM"))),
                show_lake=factor(show_lake,levels=c("YES","NO")),
                trophic_status=factor(trophic_status,levels=c("oligo-mesotrophic","eutrophic"))) %>%
  #lake_ID=factor(lake_ID,levels=c("VSS","CRE","LGP","CER-S","CSM","CER-L","JAB","VSM"))) %>%
  ggplot(.,aes(month_number,chla_mean,group=lake_ID,color=show_lake,linetype=show_lake,linewidth=show_lake))+
  #Carlson status background
  geom_rect(show.legend = F,mapping = aes(xmin = 0.4, xmax = 12.6, ymin = 0.2, ymax = 2.6),
            fill = "#97B5CD",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.4, xmax = 12.6, ymin = 2.6, ymax = 7.3),
            fill = "#2F6B9D",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.4, xmax = 12.6, ymin = 7.3, ymax = 56),
            fill = "#94D2A4",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.4, xmax = 12.6, ymin = 56, ymax = 400),
            fill = "#28A448",color=NA,alpha=0.4) +
  #Carlson status line
  #geom_hline(yintercept=2.6,color="black",linetype="dashed")+
  #geom_hline(yintercept=7.3,color="black",linetype="dashed")+
  #geom_hline(yintercept=56,color="black",linetype="dashed")+
  #geom_hline(yintercept=200,color="black",linetype="dashed")+
  #chla line and point
  geom_line(#size=0.8,
    show.legend = F)+
  geom_point(size=1.5,
             show.legend = F)+
  #visualisation
  theme_bw()+
  theme(aspect.ratio = 1,
        #strip.background = element_blank(),
        #strip.text = element_blank(),
        plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text = element_text(size=15),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  scale_color_manual(values = c("black","gray38"))+
  scale_linewidth_manual(values = c(1.3,0.8))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),expand = c(0,0),limits = c(0,13))+
  scale_y_continuous(breaks = c(2.6,7.3,56,200),limits = c(0.2,400),
                     labels = c(2.6,7.3,56,200),
                     expand = c(0,0),trans = "log")+
  # geom_rect(fill ="#28A448" ,xmin = 0.4,xmax = 12.6,
  #           ymin = 56,ymax = 400, alpha = 1, color = "black") +
  #season
  annotate("rect", xmin = 0.4, xmax = 2.5, ymin = 0.2, ymax = 0.3, fill = "#63849B",color="black") +
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = 0.2, ymax = 0.3, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 0.2, ymax = 0.3, fill = "#80A53F",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = 0.2, ymax = 0.3, fill = "#E36414",color="black") +
  annotate("rect", xmin = 11.5, xmax = 12.6, ymin = 0.2, ymax = 0.3, fill = "#63849B",color="black") +
  facet_wrap2(~ trophic_status,ncol = 2,strip=strip_color_lake<- strip_themed(
    background_x = elem_list_rect(fill = 'black',
                                  color = "black"),
    text_x = elem_list_text(colour = "white",
                            face = "bold",size=10)),scales = "fixed")+
  guides(color=guide_legend(nrow=3))+
  coord_cartesian(xlim =c(0.4,12.6), ylim = c(NA,NA),clip = "off")

ggsave("/Users/pierre/Desktop/PhD/Congrès/AquaOmics_2025/Figures/chla.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)
  
####________________________####