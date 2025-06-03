library(phyloseq)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(mixOmics)
library(readr)
library(stringr)
library(vegan)

#### Omic data ####
features_table<- t((rar_euk %>% subset_samples(month_code != "A"))@otu_table) %>% as.data.frame()
metadata <- (rar_euk %>% subset_samples(month_code != "A"))@sam_data %>% as.data.frame()

features_table_2<- t((rar_prok %>% subset_samples(month_code != "A"))@otu_table) %>% as.data.frame()
metadata_2 <- (rar_prok %>% subset_samples(month_code != "A"))@sam_data %>% as.data.frame()

#### parameters ####
Threshold= 0.9 #in ]0;1[
distance_metrics="median" #"median" or "mean"
if (distance_metrics=="median") {
  distance_metrics=median
} else distance_metrics=mean

#### Omics1 ####
omic_1.dist <- features_table %>% vegdist(., method="bray")
omic_1.ord <- cmdscale(omic_1.dist,eig = T,k=2)

TotExplVar_omic_1=((omic_1.ord$eig[1])+(omic_1.ord$eig[2]))/sum(omic_1.ord$eig)
#  as.data.frame(omic_1.ord$cum.var)[nrow(as.data.frame(omic_1.ord$cum.var)),1] #View(TotExplVar_omic_1)

if (TotExplVar_omic_1 < Threshold) {
  TotExplVar_omic_1=0
  nb_axis_omic_1=1
  while (TotExplVar_omic_1<Threshold) {
    nb_axis_omic_1=nb_axis_omic_1+1
    omic_1.ord<-cmdscale(omic_1.dist,eig = T,k=nb_axis_omic_1)
    for (i in 1:nb_axis_omic_1) {
      ExplVar_omic_1=(omic_1.ord$eig[i])
      TotExplVar_omic_1=ExplVar_omic_1+TotExplVar_omic_1
    }
    TotExplVar_omic_1=TotExplVar_omic_1/sum(omic_1.ord$eig)
    Axis_needed_omic_1<-nb_axis_omic_1
    Axis_needed_end_omic_1<-nb_axis_omic_1+2
  }
  print(paste(nb_axis_omic_1,
              "Axis are needed to explain",round(TotExplVar_omic_1,4)*100,"% of the total variance"))
  TotExplVar_omic_1=0
  nb_axis_omic_1=0
}

omic_1_df<-
  data.frame(group=metadata$lake_month,
             object=metadata$lakeID,
             time=metadata$month_code,
             time_group=metadata$season_year,
             omic_1.ord$points) %>%
  dplyr::mutate(
    lake_time_group=paste0(object,"_",time_group),.before=1) %>%
  dplyr::mutate(
    time_group=factor(time_group,levels = c("summer_2021","fall_2021","winter_2021",
                                            "spring_2022","summer_2022","fall_2022","winter_2022")))
omic_1_df<-omic_1_df[order(omic_1_df$time_group),]
omic_1_df_split<-split(omic_1_df,f = omic_1_df$object) #as many df as treatment group

for (i in 1:length(omic_1_df_split)) {
  #for each lake
  tmp_df<-omic_1_df_split[[i]] %>%
    #season centroids
    group_by(time_group,lake_time_group) %>%
    dplyr::summarise(across(starts_with("X"),.fns=distance_metrics))
  #starting by the first row and celaring a storing variable
  p=1
  tmp_distance<-as.numeric(0)
  # now compute the distance from each point to the season centroids
  for (r in 1:nrow(omic_1_df_split[[i]])) {
    if (r > p*9) { #one centroids for each season = 9 samples
      p=p+1
    }
    # for all needed axis
    for (c in 6:ncol(omic_1_df_split[[i]])) {
      distance_centroids = sum((tmp_df[p,c-3] - omic_1_df_split[[i]][r,c])^2)
      distance_centroids = sqrt(distance_centroids)
    }
    # add the distance to a vector
    tmp_distance<-c(tmp_distance,as.numeric(distance_centroids))
  }
  tmp_d.df<-data.frame(samples=rownames(omic_1_df_split[[i]]),
                       lake_month=omic_1_df_split[[i]]$group,
                       lake=omic_1_df_split[[i]]$group)
  tmp_d.df$distance_to_centroid<-tmp_distance[-1] #remove the first 0 created with the storing variable
  assign(paste0("df_ord_omic_1_",i),tmp_df)
  assign(paste0("distance_centroids_omic_1_",i,".df"),tmp_d.df)
}

for (i in 1:length(omic_1_df_split)) {
  tmp_df<- omic_1_df_split[[i]] %>% group_by(omic_1_df_split[[i]]$group) %>%
    dplyr::summarise(across(starts_with("X"),.fns=distance_metrics))
  #tmp_df[nrow(tmp_df) + 1,1]<-paste0(tmp_df[1,1],"_bis")
  #tmp_df[nrow(tmp_df),2:44] <- tmp_df[1,2:44]
  names(tmp_df)[1] <-"lake_month"
  assign(paste0("df_ord_omic_1_",i),tmp_df)
}


#### df_omic_1_BLR (1_1)
# add distance_label
Nlabel_omic_1=1
label_omic_1_1<-as.character(df_ord_omic_1_1[1,1]) #View(df_ord_omic_1_1)

while (Nlabel_omic_1<nrow(df_ord_omic_1_1)) {
  label_next_omic_1_1<-as.character(paste(df_ord_omic_1_1[Nlabel_omic_1,1],"_",df_ord_omic_1_1[Nlabel_omic_1+1,1]))
  label_omic_1_1<-c(label_omic_1_1,label_next_omic_1_1)
  Nlabel_omic_1=Nlabel_omic_1+1
}
df_ord_omic_1_1$label<-label_omic_1_1 #View(df_ord_omic_1_1)

# compute euclidean norm
eucl_norm_omic_1_1<-as.numeric(0) #View(eucl_norm_omic_1_1)
for (i in 2:nrow(df_ord_omic_1_1)) {
  sum_distance=0
  for (n in 2:(Axis_needed_omic_1+1)) {
    distance_0<-sum((df_ord_omic_1_1[i,n]-df_ord_omic_1_1[i-1,n])^2)
    sum_distance=sum_distance+distance_0
  }
  eucl_norm_omic_1_1<-c(eucl_norm_omic_1_1,sqrt(sum_distance))
}
print("Distances:")
print(eucl_norm_omic_1_1)
df_ord_omic_1_1$eucl_norm_omic_1_1<-eucl_norm_omic_1_1 #View(df_ord_omic_1_1)

# compute relative cumulative euclidean norm
rel_norm_omic_1_1<-as.numeric(0)
for (i in 2:nrow(df_ord_omic_1_1)) {
  distance<-sum(((df_ord_omic_1_1$eucl_norm_omic_1_1[i]/sum(eucl_norm_omic_1_1))*100))
  rel_norm_omic_1_1<-c(rel_norm_omic_1_1,distance)
}
print("Relative distances:")
print(rel_norm_omic_1_1)

# compute cumulative euclidean norm
cum_norm_omic_1_1<-as.numeric(0)

for (i in 2:nrow(df_ord_omic_1_1)) {
  distance<-sum(df_ord_omic_1_1$eucl_norm_omic_1_1[i]+cum_norm_omic_1_1[i-1])
  cum_norm_omic_1_1<-c(cum_norm_omic_1_1,distance)
}
print("Cumulative distances:")
print(cum_norm_omic_1_1)

# compute relative cumulative euclidean norm
cum_rel_norm_omic_1_1<-as.numeric(0)
for (i in 2:nrow(df_ord_omic_1_1)) {
  distance<-sum(((df_ord_omic_1_1$eucl_norm_omic_1_1[i]/sum(eucl_norm_omic_1_1))*100)+cum_rel_norm_omic_1_1[i-1])
  cum_rel_norm_omic_1_1<-c(cum_rel_norm_omic_1_1,distance)
}
print("Relative cumulative distances:")
print(cum_rel_norm_omic_1_1)

df_ord_omic_1_1$cum_norm_omic_1_1<-cum_norm_omic_1_1
df_ord_omic_1_1$rel_norm_omic_1_1<-rel_norm_omic_1_1
df_ord_omic_1_1$cum_rel_norm_omic_1_1<-cum_rel_norm_omic_1_1 #head(df_ord_omic_1_1)
