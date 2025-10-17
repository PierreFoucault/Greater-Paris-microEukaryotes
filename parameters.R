#### Chla ####

View(chla.df)
chla.df<-read_delim("param_data/chla_all_values.txt","\t",show_col_types = F) %>%
  filter(!is.na(Chla)) %>%
  subset(sample_code %out% "VSM_B_W1")

chla_range.df <- chla.df %>%
  group_by(lakeID) %>%
  summarise(chla_min=min(Chla),
            chla_max=max(Chla),
            range = abs(chla_max-chla_min))

chla_mean.df <- chla.df %>%
  group_by(lakeID) %>%
  summarise(chla_mean=mean(Chla),
            chla_sd=sd(Chla))

View(chla.df %>%
  group_by(lakeID,month_code) %>%
  summarise(Chla_mean=mean(Chla),
            Chla_sd=sd(Chla)))
  
chla_timeseries<- chla.df %>%
  group_by(lake_month,lakeID,month_code) %>%
  summarise(Chla_mean=mean(Chla),
            Chla_sd=sd(Chla),
            Chla_w_sd=Chla_mean+Chla_sd,
            Chla_wo_sd=Chla_mean-Chla_sd) %>%
  dplyr::mutate(lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS")),
                Chla_wo_sd = ifelse(Chla_wo_sd<=0,0.1,Chla_wo_sd)) %>%
  ggplot(.,aes(month_code,Chla_mean, group=lakeID))+
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 0.1, ymax = 2.6),
            fill = "#97B5CD",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 2.6, ymax = 7.3),
            fill = "#2F6B9D",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 7.3, ymax = 56),
            fill = "#94D2A4",color=NA,alpha=0.4) +
  geom_rect(show.legend = F,mapping = aes(xmin = 0.5, xmax = 18.5, ymin = 56, ymax = 400),
            fill = "#28A448",color=NA,alpha=0.4) +
  geom_linerange(show.legend = F,
                 aes(#ymin = ifelse(Chla_mean-Chla_sd <=0,0,Chla_mean-Chla_sd),
                      ymin = Chla_w_sd,
                     ymax = Chla_wo_sd,
                     group =lake_month ))+
                                  
  geom_line(linewidth = 0.5)+
  geom_point(size = 2)+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_blank(),
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
  scale_y_continuous(breaks = c(2.6,7.3,56,400),
                     labels = c(2.6,7.3,56,400),
                     expand = c(0,0),trans = "log10",
                     limits = c(0.1,400))+
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                # background_x = elem_list_rect(fill = 'white',
                #                               color = NA),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=12)))+
  labs(y=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))))+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.1, ymax = 0.14, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0.1, ymax = 0.14, fill = "#E36414",color="black") +
  annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 0.1, ymax = 0.14, fill = "#63849B",color="black") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = 0.1, ymax = 0.14, fill = "#EBAF47",color="black") +
  annotate("rect", xmin = 11.5, xmax = 14.5, ymin = 0.1, ymax = 0.14, fill = "#80A53F",color="black")+
  annotate("rect", xmin = 14.5, xmax = 17.5, ymin = 0.1, ymax = 0.14, fill = "#E36414",color="black")+
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = 0.1, ymax = 0.14, fill = "#63849B",color="black")+
  coord_cartesian(xlim =c(0.5,18.5), ylim = c(0.1,400),clip = "off")
chla_timeseries

ggsave("/Users/piefouca/Desktop/µEuk/Figures/chla_timeseries.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### Env.parameters ####

####__Values ####
env.df<-read_delim("param_data/env_parameters.csv",";",show_col_types = F) %>%
  dplyr::select(-Chla,-salinity,-volume_TPCN,
                -Na_589_592,-P_178_287,-S_182_034,-Fe_238_204) %>%
  #filter(!is.na(Chla)) %>%
  subset(sample_code %out% "VSM_B_W1")

env.df %>%
  melt(id.vars = c("sample_code","lakeID","month_code","season",
                   "season_year","column","lake_month","year"),
       variable.name = "parameter") %>%
  group_by(lakeID,parameter) %>%
  dplyr::summarise_at(vars("value"),list(mean=mean,sd=sd,max=max,min=min))

openxlsx2::write_xlsx(phyto_class.df, "param_data/env_data.xlsx")

####__PCA ####

PCA_param.dist <- env.df %>%
  remove_rownames() %>% column_to_rownames(var = "sample_code") %>%
  dplyr::select(-lakeID,-month_code,-season,
                -season_year,-column,-lake_month,-year)

PCA_param.dist[1:7] <-
  apply(PCA_param.dist[1:7], 2, standardize)

PCA_param.dist<- PCA_param.dist %>%
  vegdist(., method="euclidean") 

PCA_param.mds<- cmdscale(PCA_param.dist,eig=TRUE, k=2)

max(PCA_param.df$dim2)
min(PCA_param.df$dim2)

PCA_param.df <- PCA_param.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,
                lakeID=env.df$lakeID[match(rownames(.),env.df$sample_code)],
                season_year=env.df$season_year[match(rownames(.),env.df$sample_code)],
                year=env.df$year[match(rownames(.),env.df$sample_code)],
                lake_month=env.df$lake_month[match(rownames(.),env.df$sample_code)],
                month_code=env.df$month_code[match(rownames(.),env.df$sample_code)])

Fig_PCA_param <- PCA_param.df %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS")),
                season_year=factor(season_year,levels = c("summer_2021","autumn_2021","winter_2021",
                                                          "spring_2022","summer_2022","autumn_2022","winter_2022")),
                year = ifelse(season_year == "winter_2021", "2021", year),
                year=factor(year,levels = c("2021", "2022"))) %>%
  ggplot(.,aes(dim1,dim2,color=season_year,fill=season_year, shape = year, group = season_year))+
  geom_convexhull(alpha=0.6,show.legend = F,size=0.6,)+
  geom_point(show.legend = T,color="black",size=3)+
  theme_bw() +
  theme(#aspect.ratio = 0.7,
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
  guides(shape = F,fill= F,
         color=guide_legend(override.aes=list(alpha=1,size=6,
                                              fill = c("#80A53F","#E36414","#63849B","#EBAF47","#80A53F","#E36414","#63849B"),
                                              shape = c(21,21,21,22,22,22,22)),
                            color="black",nrow=2))+
  facet_wrap2(~lakeID,scale="fixed",nrow=3,remove_labels = F,
              strip = strip_themed(
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  labs(color="Season",
       x=paste0("PC1 [",round(PCA_param.mds$eig[1]*100/sum(PCA_param.mds$eig),1),"%]"),
       y=paste0("PC2 [",round(PCA_param.mds$eig[2]*100/sum(PCA_param.mds$eig),1),"%]"))

Fig_PCA_param

####__Corr. plot ####

pca_var <- env.df %>%
  remove_rownames() %>% column_to_rownames(var = "sample_code") %>%
  dplyr::select(-lakeID,-month_code,-season,-season_year,-column,-lake_month,-year)

pca_var[1:7] <-
  apply(pca_var[1:7], 2, standardize)

pca_var <- get_pca_var(prcomp(pca_var))
                     
var_dim<-pca_var$cos2 %>% .[,c(1,2)]

melted_corr <- melt(var_dim)

melted_corr <- melted_corr %>%
  rename("param"="Var1","PC"="Var2") %>%
  dplyr::mutate(.,
                PC=str_replace(PC,'Dim.', 'PC'),
                var_explain=if_else(PC=="PC1",round(PCA_param.mds$eig[1]*100/sum(PCA_param.mds$eig),1),
                                    round(PCA_param.mds$eig[2]*100/sum(PCA_param.mds$eig),1)),
                PC_percent=paste0(PC,"\n[",var_explain,"%]"))

PCA_corr<- melted_corr %>%
  dplyr::mutate(param = factor(param, levels = rev(c("Temperature","pH",
                                                 "TPN","TPC","NH4","NO3_NO2","PO43")))) %>%
  ggplot(.,aes(x = PC_percent, y = param, fill = value, size=value))+
  geom_point(color='black',shape=21) + theme_bw()+
  geom_text(aes(x = PC_percent, y = param, label = round(value, 2)),
            color = "black",size = 5)+
  theme(
    panel.grid=element_blank(),
    axis.text.x = element_markdown(face="bold",size = 15),
    axis.text.y = element_text(size = 12, face = "bold",angle =90, hjust = 0.5),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 15),
    legend.frame = element_rect(color = "black"),
    legend.key.height = unit(20, "pt"),
    legend.key.width = unit(40, "pt"),
    legend.ticks.length = unit(10,"pt"),
    legend.ticks = element_line(color="black",linewidth = 0.25),
    legend.title = element_text(face="bold",size = 25, vjust = 0.8))+
  scale_fill_gradient2(midpoint = 0.65,
                       low = "#FCFBF4",
                       mid = "#24b19c",
                       high = "#1c8776",
                       labels = c(0,0.2,0.4,0.6,0.8,1),
                       breaks = c(0,0.2,0.4,0.6,0.8,1),
                       limits = c(0, +1))+
  scale_y_discrete(labels = c("T","pH",
                                "TPN","TPC","NH<sub>4</sub><sup>+</sup>",
                                " NO<sub>3</sub><sup>-</sup>+NO<sub>2</sub><sup>-</sup>",
                                "   PO<sub>4</sub><sup>3-</sup>"),
                   limits =(rev))+
  scale_size_continuous(range = c(12,20),
                        breaks = c(0.2,0.4,0.8),
                        labels = c("0.2","0.4","0.6"))+
  guides(size=F)+
  labs(fill = expression(bold(Cos^"2")))+
  coord_flip()
PCA_corr

####__Chla - PCs ####
chla_PC1.df<-data_frame(chla_timeseries$data[,c(1:4)]) %>%
  cbind(.,
        PC1=PCA_param.df$dim1[match(.$lake_month,PCA_param.df$lake_month)],
        PC2=PCA_param.df$dim2[match(.$lake_month,PCA_param.df$lake_month)])

cor.test(chla_PC1.df$Chla_mean,chla_PC1.df$PC1,method = "spearman")
cor.test(chla_PC1.df$Chla_mean,chla_PC1.df$PC2,method = "spearman")

Fig_PC1_chla<-chla_PC1.df %>%
  ggplot(.,aes(Chla_mean,PC1))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=12),
        axis.ticks = element_line(color="black"))+
  guides(color=guide_legend(nrow = 2))+
  scale_y_continuous(trans="reverse",limits=c(7,-7),
                     breaks = c(6,3,0,-3,-6),
                     labels = c(6,3,0,-3,-6))+
  scale_x_continuous(breaks = c(1,5,10,20,50,100,200),
                     limits = c(0.43,300),transform = "log1p")+
  annotate("text",x = 80,y = 6,
    label = latex2exp::TeX("$\\textit{p} < 0.001 \\rho  -0.40$", output = "character"),
    parse = TRUE,size = 4.5,color="black")+
  labs(y="env. parameters - PC1 [33.6%]",
       x=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))))
Fig_PC1_chla

Fig_PC2_chla<-chla_PC1.df %>%
  ggplot(.,aes(Chla_mean,PC2))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=12),
        axis.ticks = element_line(color="black"))+
  guides(color=guide_legend(nrow = 2))+
  scale_y_continuous(trans="reverse",limits=c(1.3,-6.2),
                     breaks = c(1,0,-1,-3,-6),
                     labels = c(1,0,-1,-3,-6))+
  scale_x_continuous(breaks = c(1,5,10,20,50,100,200),
                     limits = c(0.43,300),transform = "log1p")+
  annotate("text",x = 150,y = 0.9,
           label = latex2exp::TeX("$\\textit{p} < 0.001 \\rho  -0.69$", output = "character"),
           parse = TRUE,size = 4.5,color="black")+
  labs(y="PC2 [28.6%]",
       x=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))))
Fig_PC2_chla

####__Overall corr. plot ####
design_param<-"
AAAA
AAAA
AAAA
BBCC"

Fig_PCA_param+PCA_corr+Fig_PC2_chla+
  plot_layout(desig=design_param,guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold",size=15),
        legend.title = element_text(size=20,face = "bold",hjust=0,vjust = 1.2),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "horizontal")

ggsave("/Users/piefouca/Desktop/µEuk/Figures/param_corr.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### Phytoplanckton ####

phyto_data.df<-read.csv("param_data/phyto_count_PR2.csv", sep = ";") %>% 
  subset(column == "W2") %>%
  merge(., read.csv("param_data/genera_biovolumes_PR2.csv", sep = ";"), by = "Taxa") %>%
  dplyr::mutate(.,rel_biovolume=nb*Biovolume)

####__Domain ####
phyto_domain.df<-phyto_data.df %>% .[,c(2,3,10,11,14)] %>%
  drop_na(rel_biovolume) %>%
  dplyr::mutate(.,
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS")),
                domain=if_else(Eukaryota.Subdivision_Prokaryota.Phylum=="Cyanobacteria","Prokaryota (Cyanobacteria)","Eukaryota"),
                domain=factor(domain,levels=c("Eukaryota","Prokaryota (Cyanobacteria)"))) %>%
  group_by(lakeID,month_code,domain) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,rel_biovolume=(rel_biovolume/sum(rel_biovolume))*100)

phyto_domain.df %>%
  drop_na(rel_biovolume) %>%
  group_by(lakeID,month_code,domain) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,biovol=rel_biovolume/sum(rel_biovolume)*100) %>%
  group_by(domain) %>%
  dplyr::summarise_at(vars("biovol"),list(mean=mean, sd=sd,max=max,min=min))

phyto_domain.df %>%
  drop_na(rel_biovolume) %>%
  group_by(lakeID,month_code,domain) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,biovol=rel_biovolume/sum(rel_biovolume)*100) %>%
  group_by(lakeID,domain) %>%
  dplyr::summarise_at(vars("biovol"),list(mean=mean, sd=sd,max=max,min=min))



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
  geom_col(width = .98, linewidth = 0.3, color ="black")+
  theme(
    axis.title = element_blank(),
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
                     labels = c('0%','25%','50%','75%','100%'))+ 
  scale_fill_manual(values = c("#308238","#CCEBC5"))+
  #scale_x_continuou(n.breaks = 12)+
  labs(fill="Phytoplankton\ndomain")+
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
phyto_domain_barplot

ggsave("/Users/piefouca/Desktop/µEuk/Figures/phyto_domain_barplot.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)



####__Class ####
phyto_class.df<-phyto_data.df %>% .[,c(2,3,4,10,11,14)] %>%
  drop_na(rel_biovolume) %>%
  dplyr::mutate(.,.after=1,
                Eukaryota.Subdivision_Prokaryota.Phylum = 
                  ifelse(Eukaryota.Subdivision_Prokaryota.Phylum=="Cyanobacteria",
                         "Cyanobacteriota",Eukaryota.Subdivision_Prokaryota.Phylum),
                Class_legend=paste0(Eukaryota.Subdivision_Prokaryota.Phylum,
                                    " · ",Class),
                lakeID=factor(lakeID,levels=c("VSM", "JAB", "CER-L",
                                              "CER-S", "CRE", "BLR","LGP", "CSM", 
                                              "VSS")),
                Class_legend=ifelse(Class_legend=="Cyanobacteriota · Cyanophyceae",
                             "Cyanobacteriota",Class_legend),
                Class_legend=factor(Class_legend,levels = c("Gyrista · Chrysophyceae",
                                                            "Gyrista · Eustigmatophyceae",
                                                            "Gyrista · Xanthophyceae",
                                                            "Gyrista · Bacillariophyceae",
                                                            "Gyrista · Mediophyceae",
                                                            "Gyrista · Coscinodiscophyceae",
                                                            "Chlorophyta · Chlorophyceae",
                                                            "Chlorophyta · Trebouxiophyceae",
                                                            "Cryptophyta · Cryptophyceae",
                                                            "Dinoflagellata · Dinophyceae",
                                                            "Euglenozoa · Euglenophyceae",
                                                            "Streptophyta · Klebsormidiophyceae",
                                                            "Streptophyta · Zygnematophyceae",
                                                            "Cyanobacteriota"))) %>%
  select(-Eukaryota.Subdivision_Prokaryota.Phylum, -Class) %>%
  group_by(lakeID,month_code,Class_legend) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,rel_biovolume=(rel_biovolume/sum(rel_biovolume))*100) %>%
  cbind(season_year=rar_euk@sam_data$season_year[match(.$month_code,
                                               rar_euk@sam_data$month_code)])
                

openxlsx2::write_xlsx(phyto_class.df, "param_data/phyto_class_biovol.xlsx")
write.csv(phyto_class.df,"param_data/phyto_class_biovol.csv",row.names = F)

phyto_class_barplot<- phyto_class.df %>%
  ggplot(., aes(fill=Class_legend, y=rel_biovolume, x=month_code))+
  geom_col(width = .98, linewidth = 0.3, color ="black")+
  theme(
    axis.title = element_blank(),
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
                     labels = c('0%','25%','50%','75%','100%'))+ 
  scale_fill_manual(values = palette_Phyto_Class)+
  #scale_x_continuou(n.breaks = 12)+
  labs(fill="Phytoplankton")+
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
phyto_class_barplot

ggsave("/Users/piefouca/Desktop/µEuk/Figures/phyto_class_barplot.png",units = "in",dpi = "retina",width = 13.4,height = 9.8)

cyano_biovolume.df <- phyto_class.df %>%
  subset(Class_legend == "Cyanobacteria · Cyanophyceae") %>%
  dplyr::mutate(.,lake_month=paste0(lakeID,"_",month_code))

phyto_class.df %>%
  drop_na(rel_biovolume) %>%
  
  dplyr::mutate(Class_legend = as.character(Class_legend),
                Class_legend = ifelse(Class_legend %in% c("Gyrista · Bacillariophyceae",
                                                          "Gyrista · Mediophyceae",
                                                          "Gyrista · Coscinodiscophyceae"),
                                      "diatoms",Class_legend),
                Class_legend=factor(Class_legend,levels = c("Gyrista · Chrysophyceae",
                                                            "Gyrista · Eustigmatophyceae",
                                                            "Gyrista · Xanthophyceae",
                                                            "diatoms",
                                                            "Chlorophyta · Chlorophyceae",
                                                            "Chlorophyta · Trebouxiophyceae",
                                                            "Cryptophyta · Cryptophyceae",
                                                            "Dinoflagellata · Dinophyceae",
                                                            "Euglenozoa · Euglenophyceae",
                                                            "Streptophyta · Klebsormidiophyceae",
                                                            "Streptophyta · Zygnematophyceae",
                                                            "Cyanobacteria · Cyanophyceae"))) %>%
  group_by(lakeID,month_code,Class_legend) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,biovol=rel_biovolume/sum(rel_biovolume)*100) %>%
  group_by(Class_legend) %>%
  dplyr::summarise_at(vars("biovol"),list(mean=mean, sd=sd,max=max,min=min))

phyto_class.df %>%
  drop_na(rel_biovolume) %>%
  group_by(lakeID,month_code,Class_legend) %>%
  summarise(rel_biovolume=sum(rel_biovolume)) %>%
  group_by(lakeID,month_code) %>%
  dplyr::mutate(.,biovol=rel_biovolume/sum(rel_biovolume)*100) %>%
  group_by(Class_legend) %>%
  dplyr::summarise_at(vars("biovol"),list(mean=mean, sd=sd,max=max,min=min))

####________________________####