#### Rar.curve ####

Fig_rar_curve_18S+Fig_rar_curve_16S

ggsave("/Users/piefouca/Desktop/Figures/rar_curves.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### BC ####

#Fig_rar_euk_BC+Fig_rar_prok_BC+

design_F2<-"
AAAA
BBBB"

Fig_rar_euk_BC_linear+Fig_rar_prok_BC_linear+
  plot_layout(desig=design_F2,guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold",size=15),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")

#ggsave("/Users/pierre/Desktop/PhD/ms2/Figures/redun.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### Chla + Cyano biovolume ####

chla_joint.df <- chla.df%>%
  group_by(lakeID,month_code) %>%
  summarise(chla_median=median(Chla)) %>%
  dplyr::mutate(.,
               # chla_median = round(log(chla_median),2),
                lake_month=paste0(lakeID,"_",month_code),
                lakeID=factor(lakeID, levels=c("VSM", "JAB", "CER-L",
                                               "CER-S", "CRE", "BLR","LGP", "CSM", 
                                               "VSS"))) %>%
  cbind(.,
        cyano_biovolume=
          cyano_biovolume.df$rel_biovolume[match(.$lake_month,
                                                 cyano_biovolume.df$lake_month)]) %>%
  dplyr::mutate(.,
               cyano_biovolume=if_else(is.na(cyano_biovolume), 0, cyano_biovolume),
               cyano_biovolume= round(cyano_biovolume,2))

#10^log(max(chla_joint.df$chla_median))

#max(chla_joint.df$cyano_biovolume)/20

ggplot(chla_joint.df, aes(month_code)) +
  geom_col(aes(y = cyano_biovolume/20 + 1),
           width = 0.8,color="black",fill="grey") +
  geom_point(aes(y = log10(chla_median)),
                 size = 1,color = "black") +
  scale_y_continuous(
    labels = function(x) scales::comma(exp(x)),
   # breaks = c(2.6,7.3,56,400),
    minor_breaks = NULL,
    sec.axis = sec_axis(~(.-1)*20))+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 8),
        axis.text.y =element_text(color = "black", size = 10),
        panel.grid = element_blank())+
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))

#%>%
  #dplyr::mutate(.,
#View(chla_joint.df)
#ratio = with(chla_joint.df, max(cyano_biovolume)/max(chla_median))
#scaleRight <- max(chla_joint.df$cyano_biovolume)/log(100)

# ggplot(CountryData) +
#   geom_area(aes(x=DateReport1,y=TotalCases_log10*ratio), fill="grey80") +
#   geom_col(aes(x=DateReport1,y=NewCases), fill="#D86422",size=0.6) +
#   scale_y_continuous(labels=comma, name="Daily Cases",
#                      sec.axis = sec_axis(~ ./ratio, 
#                                          labels=function(x) comma(round(10^x), accuracy=1),
#                                          breaks=c(0:5, 0:5+0.5),
#                                          name="Cumulative Cases")) +
#   theme(axis.text.y=element_text(colour="#D86422"),
#         axis.text.y.right=element_text(colour="grey50"),
#         axis.title.y=element_text(colour="#D86422"),
#         axis.title.y.right=element_text(colour="grey50"))

#chla_cyanobacteria_timeseries <-
#chla_joint.df %>%
  ggplot(.,aes(month_code,cyano_biovolume))+
  geom_col(width = 0.8,color="black",fill="grey")+
  geom_line(linewidth = 0.5,
            mapping=aes(x=month_code, y=chla_median*ratio, group = lakeID))+
            #data = chla_joint.df, inherit.aes = F)+
  geom_point(size = 2,
             mapping=aes(x=month_code, y=chla_median*ratio))+
             #data = chla_joint.df, inherit.aes = F)+
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
                              'J','A','S','O','N', 'D'))+
  scale_y_continuous(#first axis
    # name = 'Cyanobacteria biovolume',
    expand = c(0,0),
    #limits = c(NA,400),
   # breaks = c(0,100,200,300,400),
    #labels = c("0%","25%","50%","75%","100%"),
    #second axis
    sec.axis=
      sec_axis(~ .,trans = "log",
               #labels=function(x) comma(round(10^x), accuracy=1),
               #labels=function(x) 10^(x),
             #breaks=c(0:5, 0:5+0.5),
               #breaks = c(log(2.6),log(7.3),log(56),log(400)),
               #labels = c("2.6","7.3","56","400"),
               #name=expression(paste(bold("Chl"),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))),
              ))+
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))+
  labs(y="Cyanobacteria biovolume")
  # annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.1, ymax = 0.18, fill = "#80A53F",color="black")+
  # annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0.1, ymax = 0.18, fill = "#E36414",color="black") +
  # annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 0.1, ymax = 0.18, fill = "#63849B",color="black") +
  # annotate("rect", xmin = 8.5, xmax = 11.5, ymin = 0.1, ymax = 0.18, fill = "#EBAF47",color="black") +
  # annotate("rect", xmin = 11.5, xmax = 14.5, ymin = 0.1, ymax = 0.18, fill = "#80A53F",color="black")+
  # annotate("rect", xmin = 14.5, xmax = 17.5, ymin = 0.1, ymax = 0.18, fill = "#E36414",color="black")+
  # annotate("rect", xmin = 17.5, xmax = 18.5, ymin = 0.1, ymax = 0.18, fill = "#63849B",color="black")+
 # coord_cartesian(xlim =c(0.5,18.5), ylim = c(0,400),clip = "off")
#chla_cyanobacteria_timeseries


#scaleFactor <- log(chla_joint.df$chla_median)

# ymax <- max(cyano_biovolume.df$rel_biovolume)
#

chla_cyanobacteria_timeseries<- chla.df %>%
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
  geom_col(width = 0.8,color="black",fill="grey",
           mapping=aes(x=month_code, y=rel_biovolume),
           data = cyano_biovolume.df, inherit.aes = F)+
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
#chla_cyanobacteria_timeseries




set.seed(2)
CountryData = tibble(
  DateReport1 = seq(as.Date("2020-03-01"), as.Date("2020-05-01"), "1 day"),
  NewCases = sample(100:10000, length(DateReport1)),
  TotalCases_log10 = log10(cumsum(NewCases))
)

# Ratio of largest NewCases value to largest TotalCases_log10 value
# Also multiplied by 1.5 to make the grey area plot extend above the bars
ratio = 1.5 * with(CountryData, max(NewCases)/max(TotalCases_log10))
View(CountryData)
ggplot(CountryData) +
  geom_area(aes(x=DateReport1,y=TotalCases_log10*ratio), fill="grey80") +
  geom_col(aes(x=DateReport1,y=NewCases), fill="#D86422",size=0.6) +
  scale_y_continuous(labels=comma, name="Daily Cases",
                     sec.axis = sec_axis(~ ./ratio, 
                                         labels=function(x) comma(round(10^x), accuracy=1),
                                         breaks=c(0:5, 0:5+0.5),
                                         name="Cumulative Cases")) +
  theme(axis.text.y=element_text(colour="#D86422"),
        axis.text.y.right=element_text(colour="grey50"),
        axis.title.y=element_text(colour="#D86422"),
        axis.title.y.right=element_text(colour="grey50"))


#Voici un exemple de code en R utilisant ggplot2 pour générer un graphique avec deux axes Y, l'un en échelle logarithmique et l'autre en échelle linéaire. Les données pour l'axe linéaire varient entre 0 et 100, et pour l'axe logarithmique, elles varient entre 0.2 et 400.

# Charger la bibliothèque ggplot2
#library(ggplot2)

# Créer un data frame avec des données d'exemple
data <- data.frame(
  x = 1:10,
  y_linear = seq(0, 100, length.out = 10),
  y_log = seq(0.1, 400, length.out = 10)
)

# Générer le graphique avec deux axes Y
ggplot(data) +
  # Ajouter les données pour l'axe Y linéaire
  geom_line(aes(x = x, y = y_linear*4), color = "blue") +
  # Ajouter les données pour l'axe Y logarithmique
  geom_line(aes(x = x, y = y_log+1), color = "red") +
  # Définir l'échelle de l'axe Y de gauche (linéaire)
  scale_y_continuous(
    name = "Axe Y linéaire",
    breaks = c(0, 100, 200, 300, 400),
    labels = c("0%","25%","50%","75%","100%"),
    sec.axis = sec_axis(~ .,
                        name = "Axe Y logarithmique",
                        trans = "log",
                        breaks = c(log10(0.1+1),log(2.6+1),log(56+1),log(400+1)))) +
  # Ajouter un titre et des labels
  labs(title = "Graphique avec deux axes Y", x = "Axe X") +
  # Ajouter une légende
  theme_bw()


ggplot(chla_joint.df, aes(month_code)) +
  geom_col(aes(y = cyano_biovolume/20 + 1),
           width = 0.8,color="black",fill="grey") +
  geom_point(aes(y = log10(chla_median)),
             size = 1,color = "black") +
  scale_y_continuous(
    labels = function(x) scales::comma(exp(x)),
    # breaks = c(2.6,7.3,56,400),
    minor_breaks = NULL,
    sec.axis = sec_axis(~(.-1)*20))+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 8),
        axis.text.y =element_text(color = "black", size = 10),
        panel.grid = element_blank())+
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))


data <- data.frame(
  x = 1:10,
  y_linear = seq(0, 100, length.out = 10),
  y_log = seq(0.1, 400, length.out = 10)
)

ggplot(chla_joint.df) +
  geom_line(aes(x = month_code, y = log(chla_median), group = lakeID), color = "black") +
  geom_point(aes(x = month_code, y = log(chla_median)), color = "black") +
  geom_col(aes(x = month_code, y = cyano_biovolume), color = "black", fill="lightgrey") +
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(face="bold"),
        axis.ticks.x =element_blank(),
        axis.ticks.y =element_line(color = "black"),
        axis.text.x =element_text(color = "black", size = 8),
        axis.text.y =element_text(color = "black", size = 10),
        panel.grid = element_blank())+
  scale_y_continuous(name = "Cyanobacteria biovolume",
                     limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%","25%","50%","75%","100%"),
                     sec.axis = sec_axis(~ ., name = "Chla",
                                         trans = "log",
                                        breaks = c(0.1, 2.6, 56, 100)
                                        )) +
  facet_wrap2(~ lakeID, scales = "fixed",remove_labels = 'x',
              strip = strip_themed(
                background_x = elem_list_rect(fill = 'lightgrey',
                                              color = "black"),
                text_x = elem_list_text(colour = "black",
                                        face = "bold",size=10)))


#### F1 ####

#map chla+cyano_biovol
#trends 18S barplot season
#trends 18S barplot season

#### F2 ####

#richess 18S PCoA
#ratio et autre data mode trophique

#### F3 ####

#TLA_joint+MOTA_joint

#### F4 ####
#network lake et PCA des characteristics (1 point = 1 lake-saison)

####________________________####