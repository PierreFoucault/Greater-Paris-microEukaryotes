


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

