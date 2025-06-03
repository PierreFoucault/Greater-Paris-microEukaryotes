#### TLA joint ####

####__BC by day_gap ####
joint_day_gap.df<-
  TLA_prok.df %>% .[,c(1:3,6,10,13)] %>% dplyr::mutate(.,domain="Prokaryota") %>%
  rbind(.,(TLA_euk.df[,c(1:3,6,10,13)] %>% dplyr::mutate(.,domain="Eukaryota"))) %>%
  dplyr::mutate(.,
                G1_lakeID=factor(G1_lakeID,
                                 levels=c("VSM", "JAB", "CER-L",
                                          "CER-S", "CRE", "BLR","LGP", "CSM", 
                                          "VSS")),
                domain=factor(domain,levels = c("Eukaryota","Prokaryota")))

####__plot ####
ggplot(joint_day_gap.df,aes(x=day_gap,y=BC,group=domain,fill=domain,color=domain))+
  geom_point(size=0.5,alpha=0.4,show.legend = T)+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly3,show.legend = F,
                 data =joint_day_gap.df[(joint_day_gap.df$poly=="poly3"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly4,show.legend = F,
                 data =joint_day_gap.df[(joint_day_gap.df$poly=="poly4"),])+
  stat_poly_line(alpha=0.4,formula = lm_formula_poly5,show.legend = F,
                 data =joint_day_gap.df[(joint_day_gap.df$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=8),
        axis.ticks = element_line(color="black"),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom",
        legend.direction = "vertical")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),expand = c(0,0),limits = c(0.1,1))+
  scale_x_continuous(breaks = c(30, 60, 120, 180, 360, 540),expand = c(0,0),limits = c(-1,545))+
  scale_fill_manual(values = c("#F57D05","#0A82FA"))+
  scale_color_manual(values = c("#F57D05","#0A82FA"))+
  guides(fill=F,
         color=guide_legend(nrow=1,
                            override.aes = list(shape=22,size=5,alpha=1,fill=c("#F57D05","#0A82FA"),color="black")))+
  labs(y="Bray-Curtis dissimilarity",color="Domain",
       x="Time between sampling dates (days)")+
  facet_wrap2(~ G1_lakeID,nrow= 3,scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))

ggsave("/Users/piefouca/Desktop/µEuk/Figures/TLA_joint.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####|####

#### MOTA joint ####

###__corr. ####

corr_MOTA<-MOTA_lake_16S$data %>%
  group_by(lakeID) %>% summarise(max_dist_prok=median(max_dist)) %>%
  cbind(MOTA_lake_18S$data %>% summarise(max_dist_euk=median(max_dist))) %>%
  .[,-3]

cor_test(data=corr_MOTA,
         vars = max_dist_euk,
         vars2 = max_dist_prok,
         method = "spearman")

grDevices::cairo_pdf("/Users/piefouca/Desktop/µEuk/Figures/corr_MOTA_maxdist.pdf",
                     width = 13.4,height = 9.8,fallback_resolution = 300)


corr_MOTA %>%
  dplyr::mutate(lakeID= factor(lakeID,
                               levels = c("VSM", "JAB", "CER-L",
                                          "CER-S", "CRE", "BLR","LGP", "CSM", 
                                          "VSS"))) %>%
  ggplot(.,aes(max_dist_euk,max_dist_prok, fill=lakeID))+
  geom_smooth(method = "lm",fill="lightgrey",color="gray40")+
  geom_point(shape= 21, color= "black", show.legend = F, size= 2.5)+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(size=12, color = "black"),
        axis.title = element_text(face = 'bold', size =12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12, face = 'bold'))+
  scale_y_continuous(expand = c(0,0),limits = c(3.8,11),
                     breaks = c(4, 6, 8, 10))+
  scale_x_continuous(expand = c(0,0),limits = c(3.8,11),
                     breaks = c(4, 6, 8, 10))+
  scale_fill_manual(values = rev(palette_lake_chla))+
  annotate("text",x = 4.2, y= 10,hjust=0,
           expression(paste(italic("p"),">0.05")),
           label= expression(paste(italic("p"),"<0.01 \u03C1 0.87")))+
  labs(x='micro-Eukaryotic community trajectories',
       y="Prokaryotic community trajectories")
dev.off()


ggsave("/Users/piefouca/Desktop/µEuk/Figures/corr_MOTA_maxdist.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####__BC by day_gap ####

####________________________####
#### Figure PhD defense ####
####________________________####

phd_TLA<-joint_day_gap.df %>%
  # subset(lake %in% c("JAB")) %>%
  subset(lake %in% c("JAB","VSS")) %>%
  dplyr::mutate(.,community=factor(community,levels = c("Prokaryotes","µ-Eukaryotes")),
                lake_community=paste0(lake,"_",community)) %>% #subset(community =="Prokaryotes") %>%
  subset(month_gap<13)
# dplyr::mutate(.,community=factor(community,levels = c("Prokaryotes","µ-Eukaryotes")),
#               lake_community=paste0(lake,"_",community))

phd_TLA %>%
  ggplot(.,
         aes(x=month_gap,y=BC,group=community,fill=lake,color=lake,linetype = community))+
  geom_point(size=0.5,alpha=0.4,show.legend = F)+#color="grey")+
  stat_poly_line(alpha=0.3,formula = lm_formula_poly3,show.legend = F,size=1.2,
                 data =phd_TLA[(phd_TLA$poly=="poly3"),])+
  stat_poly_line(alpha=0.3,formula = lm_formula_poly4,show.legend = F,size=1.2,
                 data =phd_TLA[(phd_TLA$poly=="poly4"),])+
  stat_poly_line(alpha=0.3,formula = lm_formula_poly5,show.legend = F,size=1.2,
                 data =phd_TLA[(phd_TLA$poly=="poly5"),])+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(size=12),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "right")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),expand = c(0,0),limits = c(0,1))+
  scale_x_continuous(breaks = c(1,3,6,12),expand = c(0,0),limits = c(0.5,12.5))+
  scale_fill_manual(values = c("#2F6B9D","#28A448"))+
  scale_color_manual(values = c("#2F6B9D","#28A448"))+
  guides(fill=F,linetype=F,
         color=guide_legend(ncol=1,
                            override.aes = list(shape=22,size=5,alpha=1,fill=c("#2F6B9D","#28A448"),color="black")))+
  labs(y="Taxa-content dissimilarity",color="Domain",
       x="Month interval")
# facet_wrap2(~ lake,ncol = 3,scale="fixed",
#             #strip.position="left",
#             strip = strip_color_lake<- strip_themed(
#               background_x = elem_list_rect(fill = "black",
#                                             color = "black"),
#               text_x = elem_list_text(colour = "white",face = "bold",size=9)))

ggsave("/Users/pierre/Desktop/PhD/ms2/Figures/TLA_Prok_Euk.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)
