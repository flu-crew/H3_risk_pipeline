#! /usr/bin/env Rscript
# Last Update: 2019/09/20

# ==== User Variables
# Set working Directory
getwd()
setwd("/Users/Carine.Souza/Desktop/Human\ sera_R2")

# Input files
post_exposure.xlsx = "Copy\ of\ Human\ sera_H3\ results_post\ exposure\ cohort_v9.xlsx"
post_exposure.sheet = "Data_Mod"
post_exposure.sheetFC = "FoldChange"

post_vac.xlsx = "Copy_of_Human_sera_H3_results_Post_Vacc_cohort_v9.xlsx"
post_vac.sheet = "Data_Mod"
post_vac.sheetFC = "FoldChange"

# ==== Libraries
library(ggplot2)      # Plots
library(tidyverse)    # Clean data
library(magrittr)     # Pipe style coding
library(readxl)       # Loads Excell File
library(reshape2)     # Rearrange cols and rows
library(dplyr)        # generates summaries across groups of cols or rows

# ==== Function

# Geometric Mean Titers
# X = LOG2(Titer/10)
HI_to_GMT <- function(hi_titers){
  gmt <- log(hi_titers/10, 2)
  return(gmt)
}

# === Read in Data
data <- read_excel(post_exposure.xlsx, 
                   sheet = post_exposure.sheet)

dataFC <- read_excel(post_exposure.xlsx, 
                   sheet = post_exposure.sheetFC)

# === Plot a comparison of antigens against each other
cdata <- data %>% t(.) %>% data.frame(.)
#names(cdata)<-paste(data$header,data$Clade, sep=":")
names(cdata)<-paste(data$Clade, sep=":")
plot(cdata,lower.panel=NULL,
     main="Comparing antigens based on human sera")

# === Correlation of Antigens
# Tutorial: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
cor_input <- data[,4:11] %>% t(.)  %>% HI_to_GMT(.)
#colnames(cor_input) <- paste(data$header, data$Clade, sep = ".")
colnames(cor_input) <- paste(data$Clade, sep = ".")
cdata_corr <- cor(cor_input, method="spearman")

# === draft
mcdata_corr <- melt(cdata_corr) %>%
  mutate(
    HuVac=Var1,
    Strain=Var2,
    Spearman=value,
    Var1=NULL,
    Var2=NULL,
    value=NULL
  ) %>%{
    subset(., grepl("HuVac",.$HuVac))
  } %>% {
    subset(., !grepl("HuVac",.$Strain))
  } %>%
  mutate(
    HuVac=factor(HuVac, levels=c("VI75.HuVac","BE92.HuVac","WU95.HuVac","HK14.HuVac"))
  )

(p<-mcdata_corr %>% 
    ggplot(. , aes(x=HuVac, y=Strain, fill=Spearman))+
    geom_tile()+
    geom_text(aes(label=round(Spearman,digits=1)),size=4)+
    theme_classic()+
    scale_fill_gradient2(low="#0000FF",mid = "#FFFFFF",midpoint=0, high = "#FF0000")+
    labs(title="Correlation of HI titers",x="Human vaccine strains",y="Swine H3N2 strains")
)

ggsave("Post_exposure_correlation_heatmap2_v9.png",plot=p,width=6, height=4)


 
# === end draft

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

#Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
 cormat[lower.tri(cormat)] <- NA
return(cormat)
 }

#ant_levels=c("NY11.3.1990.4a","IN16.3.1990.4b","MX15.3.1990.2","AB14.3.1990.4","OH17.3.2010.1",
#             "OK17.3.2010.2","VI75.HuVac","BE92.HuVac","HK14.HuVac")

ant_levels=c("MO18.1990.1","NY11.1990.4a","IN16.1990.4b","MX15.1990.2","AB14.1990.4","OH17.2010.1",
            "OK17.2010.2","VI75.HuVac","BE92.HuVac","HK14.HuVac")
#ant_levels=c("3.1990.1","3.1990.4a","3.1990.4b","3.1990.2","3.1990.4","3.2010.1",
            # "3.2010.2","VI75.HuVac","BE92.HuVac","WU95.HuVac","HK14.HuVac")

mcdata_corr<-cdata_corr[ant_levels,ant_levels] %>% get_lower_tri(.) %>% melt(.,na.rm=TRUE) %>%
  mutate(
    Var1=factor(Var1, levels=ant_levels),
    Var2=factor(Var2, levels=ant_levels)
  )

#%>% get_lower_tri(.)

(p<-mcdata_corr %>% ggplot(.,aes(x=Var1,y=Var2)) +
  geom_tile(aes(fill=value),color="white") +
  geom_text(aes(label=round(value,digits = 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1),name="Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(title="Correlation of Antigens based on human sera",y="",x="") +
  coord_fixed())

ggsave("Post_exposure_heatmap_v9.png",plot=p,width=10.4, height = 7.36)
# mdata <- data %>% melt(.,id=c("Antigen","header"))

# 1) Plot the HI titers for each swine virus against the HI titers 
#    for the HK/14 and Beth/15 strains for the post-exposure cohort
#    to see if they are correlated.
#data <- read_excel("Human\ sera_H3\ results_post\ exposure\ cohort.xlsx",sheet="Data_Mod")
tdata <- data[,4:12] %>% t(.) %>% data.frame(.)
#names(tdata)<-paste(data$header, data$Clade, sep=".")
names(tdata)<-paste(data$Clade, sep=".")
tdata$person=rownames(tdata)

VI75 <- tdata %>% 
  melt(.,id=c("person","VI75.HuVac"))  %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  subset(.,!grepl("Circ",.$variable)) %>%
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="VI75.HuVac",
    Ref_Titor=VI75.HuVac
  ) %>%
  select(person,Swine_Strain,Swine_Titor,Ref_Strain,Ref_Titor)

BE92 <- tdata %>% 
  melt(.,id=c("person","BE92.HuVac"))  %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  subset(.,!grepl("Circ",.$variable)) %>%
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="BE92.HuVac",
    Ref_Titor=BE92.HuVac
  ) %>%
  select(person,Swine_Strain,Swine_Titor,Ref_Strain,Ref_Titor)

WU95 <- tdata %>% 
  melt(.,id=c("person","WU95.HuVac"))  %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  subset(.,!grepl("Circ",.$variable)) %>%
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="WU95.HuVac",
    Ref_Titor=WU95.HuVac
  ) %>%
  select(person,Swine_Strain,Swine_Titor,Ref_Strain,Ref_Titor)

HK14 <- tdata %>% 
  melt(.,id=c("person","HK14.HuVac"))  %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  subset(.,!grepl("Circ",.$variable)) %>%
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="HK14.HuVac",
    Ref_Titor=HK14.HuVac
  ) %>%
  select(person,Swine_Strain,Swine_Titor,Ref_Strain,Ref_Titor)



mdata <- rbind(VI75, BE92, WU95, HK14) %>%
  mutate(
    Ref_Strain=factor(Ref_Strain,levels=c("VI75.HuVac", "BE92.HuVac", "WU95.HuVac","HK14.HuVac")),
    Swine_Strain=factor(Swine_Strain, levels = c("MO18.1990.1","NY11.1990.4a","IN16.1990.4b","MX15.1990.2","AB14.1990.4","OH17.2010.1",
                                                 "OK17.2010.2",
                                                 "VI75.HuVac", "BE92.HuVac", "HK14.HuVac"))
  )

# ==== Correlation scatterplot with line and text
mcdata_corr<-cdata_corr  %>% melt(.,na.rm=TRUE) 
ss<-mdata$Swine_Strain %>% unique(.) %>% as.character(.)
rs<-mdata$Ref_Strain %>% unique(.) %>% as.character(.)

temp <- mcdata_corr %>% 
  mutate(
    Swine_Strain = Var1 %>% as.character(.),
    Ref_Strain = Var2 %>% as.character(.),
    Corr=value,
    Var1=NULL,
    Var2=NULL,
    value=NULL,
    x=2,
    y=7
  ) %>%
  subset(., Swine_Strain %in% ss) %>%
  subset(., Ref_Strain %in% rs) %>%
  mutate(
    Ref_Strain=factor(Ref_Strain,levels=c("VI75.HuVac", "BE92.HuVac","WU95.HuVac","HK14.HuVac")),
    Swine_Strain=factor(Swine_Strain, levels = c("MO18.1990.1","NY11.1990.4a","IN16.1990.4b","MX15.1990.2","AB14.1990.4","OH17.2010.1",
                                                 "OK17.2010.2",
                                                 "VI75.HuVac", "BE92.HuVac", "WU95.HuVac","HK14.HuVac")
                        )
  )


(p<-mdata %>%
  ggplot(., aes(x = log(Ref_Titor/10, 2),
                y = log(Swine_Titor/10, 2))) +
#  geom_smooth(method="auto",se=FALSE) +
  geom_point(size=2)+
  geom_text(data = temp, mapping=aes(x=x,y=y,label=round(Corr,digits=2)))+
  facet_grid(Swine_Strain~Ref_Strain) +
  labs(#title="Swine titers against reference titers", 
       x="GMT HuVac ",y="GMT swine")+
  theme_minimal()+
  theme(strip.text.y = element_text(angle = 360),
        panel.border = element_rect(colour = "black",fill=NA)
        )+
  expand_limits(x=0, y=0)+
  xlim(0,13)+
  coord_fixed())

ggsave("Post_exposure_correlation_bw_v9.png",plot=p,width=5.82, height=7)

(p<-mdata %>%
    ggplot(., aes(x = log(Ref_Titor/10, 2),
                  y = log(Swine_Titor/10, 2))) +
    #geom_smooth(method="lm",se=FALSE) +
    geom_jitter(size=2,alpha=0.8)+ #
    geom_text(data = temp, mapping=aes(x=x,y=y,label=round(Corr,digits=2)))+
    facet_grid(Swine_Strain~Ref_Strain) +
    labs(#title="Swine titers against reference titers", 
      x="GMT HuVac",y="GMT swine")+
    theme_minimal()+
    theme(strip.text.y = element_text(angle = 360),
          panel.border = element_rect(colour = "black",fill=NA)
    )+
    expand_limits(x=0, y=0)+
    xlim(0,13)+
    coord_fixed())
ggsave("Post_exposure_correlation_color_v9.png",plot=p,width=5.82, height=7)

#mdata %>%
#  ggplot(.,aes(person, log(Swine_Titor)))+
#  geom_col()+
#  geom_line(data=mdata, aes(x=Ref_Strain,y=log(Ref_Titor),group=person))+
#  facet_grid(~Swine_Strain)+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#cdata<-data %>% melt(., id=c("Antigen","header","Clade")) %>% 
cdata<-data %>% melt(., id=c("Antigen","Clade")) %>% 
  mutate(
    #Strain=paste(header,Clade,sep="."),
    Strain=paste(Clade,sep="."),
    Strain=factor(Strain,levels=c("VI75.HuVac","BE92.HuVac", "WU95.HuVac", "HK14.HuVac",
                                  "MO18.1990.1","NY11.1990.4a","IN16.1990.4b","MX15.1990.2","AB14.1990.4","OH17.2010.1",
            "OK17.2010.2"))
  )
#Overall_post-exposure Violin plots
str(cdata)
cdata$vac = grepl("HuVac" , cdata$Clade)
annot_df <- cdata %>%
  group_by(Strain) %>%
  summarise(
    average = round(mean (value)),
    pos = round(sum(HI_to_GMT(value) >= 2)/10*100)
  ) 
cdata %>% ggplot(., aes(x=Strain,y=HI_to_GMT(value)#,
                        #color=Strain,
                        #fill=Strain
                        ))+
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  
  geom_violin(width=1,trim=FALSE)+
  geom_boxplot(color="black",
               aes(fill=vac) , ##remove fill
               width=0.25) +
  scale_fill_manual(values=c("TRUE"="#969696","FALSE"="#FFFFFF"))+ # remove fill
  #stat_summary(fun.y="mean")+
  geom_jitter(color="black",width=0.05,height=0, size=1.2)+
  stat_summary(fun.y="mean", color="red", size=0.3)+
  geom_text(data = annot_df, aes(label = average, 
                              y=16))+
  geom_text(data = annot_df, aes(label = paste("(",
                                               pos,
                                               "%)",
                                               sep=""
                                                ),
                                 y=15)
            )+                         
  theme_minimal()+
  #geom_violin(trim=FALSE, fill='white', color="white")+
  facet_grid(.~Strain,switch="both",scale="free")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill="NA"),
        legend.position = "none",
        axis.title.x=element_blank(),
       # axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
       text=element_text(size = 14)
  ) +
  scale_y_continuous(breaks=c(0:8)*2,limits=c(0,16))+
  scale_x_discrete(breaks=NULL)+
  #  ylim(0,14)+
  labs(title="Post-exposure sera cohort - HI Titers",y="Log2 GMT",x="Strain")

  
ggsave("Post_exposure_cohort_HI violin plots_v9.png",dpi=600,width=14,height=4.88)

#Overall_post-exposure Box plots
cdata %>% ggplot(., aes(x=Strain,y=HI_to_GMT(value),
                        color=Strain,fill=Strain))+
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  
  #geom_violin(width=1,trim=FALSE)+
  geom_boxplot(color="black",fill= "NA", width=0.25)+
  geom_jitter(color="black",width=0.05,height=0)+
  theme_minimal()+
  #geom_violin(trim=FALSE, fill='white', color="white")+
  facet_grid(.~Strain,switch="both",scale="free")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill="NA"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  scale_y_continuous(breaks=c(0:7)*2,limits=c(0,12))+
  scale_x_discrete(breaks=NULL)+
#  ylim(0,14)+
  labs(title="Post-exposure sera cohort - HI Titers",y="GMT",x="Strain")

ggsave("Post_exposure_cohort_HI_boxplots_v8.png",dpi=600,width=12.3,height=4.88)

cdata %>% 
  mutate(
    responder=case_when(grepl("High",variable)~"High responder (HI > 320)",
                        grepl("Med", variable)~"Med responder (HI = 80 - 320)"),
    responder=factor(responder, levels=c("Med responder (HI = 80 - 320)","High responder (HI > 320)"))
  ) %>%
  ggplot(., aes(x=responder,y=HI_to_GMT(value),
                        color=responder,fill=responder))+
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(width=1,trim=FALSE)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(color="black",width=0.05,height=0)+
  theme_minimal()+
  facet_grid(.~Strain,switch="both",scale="free")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  scale_y_continuous(breaks=c(0:7)*2,limits=c(0,13))+
  scale_x_discrete(breaks=NULL)+
  #  ylim(0,14)+
  labs(title="Post Exposure Cohort HI Titers",y="GMT",x="Strain")

ggsave("Post_exposure_cohort_HI_responders.png",dpi=300,width=12.3,height=4.88)

# Attempt at Fold Change (Did not use this part)-------------------------------------
dataFC <- read_excel(post_exposure.xlsx, 
                     sheet = post_exposure.sheetFC)
names(dataFC)
library(scales)
#mdataFC <- dataFC %>% melt(id=c("Antigen","header","Clade"))
mdataFC <- dataFC %>% melt(id=c("Antigen", "Clade"))
mdataFC %>% 
  # mutate(
  #   header=paste(header,Clade, sep="."),
  #   header=factor(header, levels=c("HK14.HuVac","VI75.HuVac","BE92.HuVac","BT15.Circ.H3",
  #                                  "MO18.1990.1","NY11.3.1990.4a","IN16.3.1990.4b","MX15.3.1990.2","AB14.3.1990.4","OH17.3.2010.1","OK17.3.2010.2"))
  # )  %>% 
  mutate(
    header=paste(Clade, sep="."),
    header=factor(levels=c("HK14.HuVac","VI75.HuVac","BE92.HuVac","BT15.Circ.H3",
                                   "3.1990.1","3.1990.4a","3.1990.4b","3.1990.2","3.1990.4","3.2010.1","3.2010.2"))
  )  %>% 
  ggplot(.,aes(x=header,y=value,fill=header))+
  geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
  geom_violin(aes(color=header))+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(height = 0, width=0.2)+
  #scale_y_continuous(breaks=c(-10:4),limits=c(-10,4))+
  #scale_y_continuous(breaks=c(-40:4),labels=paste("2^",c(-40:4),sep=""))+
  scale_y_continuous(breaks=c(-9:9),labels=c(-9:9))+
  labs(y="Fold Change",title="Post-exposure cohort - Fold change to HK14.HuVac")+
  theme_minimal()+
  facet_grid(.~header,switch="both",scale="free")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  scale_x_discrete(breaks=NULL)
#+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 2^x),
                #labels = trans_format("log10", math_format(2^.(2^x)))) 

ggsave("Post_exposure_cohort_foldchange.png",dpi=300,width=12.3,height=4.88)
# ----------------------------------------------------
# 2)	Group the male and female together,
# and plot the HI titers for each swine virus against the HI titers
# for the 4 human vaccine strains for the post-vacc cohort to see if they are correlated.
  
# === Read in Data
data2 =  readxl::read_excel(post_vac.xlsx, sheet=post_vac.sheet)
  
#tdata <- data[,9:17] %>% t(.) %>% data.frame(.)
 
VI75 <- data2 %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                             "pre_vaccine","post_vaccine2","Fold_change",
                             "VI75.HuVac")) %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="VI75.HuVac",
    Ref_Titor=VI75.HuVac,
    VI75.HuVac=NULL,
    variable=NULL,
    value=NULL
  ) 

BE92 <- data2 %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                             "pre_vaccine","post_vaccine2","Fold_change",
                             "BE92.HuVac")) %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="BE92.HuVac",
    Ref_Titor=BE92.HuVac,
    BE92.HuVac=NULL,
    variable=NULL,
    value=NULL
  ) 
WU95 <- data2 %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                             "pre_vaccine","post_vaccine2","Fold_change",
                             "WU95.HuVac")) %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="WU95.HuVac",
    Ref_Titor=WU95.HuVac,
    WU95.HuVac=NULL,
    variable=NULL,
    value=NULL
  ) 

HK14 <- data2 %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                             "pre_vaccine","post_vaccine2","Fold_change",
                             "HK14.HuVac")) %>% 
  subset(.,!grepl("HuVac",.$variable)) %>% 
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value,
    Ref_Strain="HK14.HuVac",
    Ref_Titor=HK14.HuVac,
    HK14.HuVac=NULL,
    variable=NULL,
    value=NULL
  ) 

mdata <- rbind(VI75,BE92,WU95,HK14)

# == correlation 
cor_input<-data2[,-c(1:8)] %>% HI_to_GMT(.)
cdata_corr <- cor(cor_input, method="spearman")

# === draft
str(cdata_corr)
mcdata_corr <- melt(cdata_corr) %>%
  mutate(
    HuVac=Var1,
    Strain=Var2,
    Spearman=value,
    Var1=NULL,
    Var2=NULL,
    value=NULL
  ) %>% {
    subset(., grepl("HuVac", .$HuVac))
  } %>% {
    subset(., !grepl("HuVac", .$Strain))
  }

# Heatmap of correlation
(p<-mcdata_corr %>% 
    ggplot(. , aes(x=HuVac, y=Strain, fill=Spearman))+
    geom_tile()+
    geom_text(aes(label=round(Spearman,digits=1)),size=4)+
    theme_classic()+
    scale_fill_gradient2(low="#0000FF",mid = "#FFFFFF",midpoint=0, high = "#FF0000")+
    labs(title="Correlation of HI titers",x="Human vaccine strains",y="swine H3N2 strains")
)

ggsave("Post_vac_correlation_heatmap_v9.png",plot=p,width=6, height=4)

# === end draft
ss=mdata$Swine_Strain %>% unique(.) %>% as.character(.)
rs=mdata$Ref_Strain %>% unique(.) %>% as.character(.)

text_corr <- cdata_corr %>% melt() %>%
  mutate(
    Swine_Strain=Var1,
    Ref_Strain=Var2,
    Corr=value,
    Var1=NULL,
    Var2=NULL,
    value=NULL,
    Decade=40
  ) %>%
  subset(.,Swine_Strain %in% ss) %>%
  subset(.,Ref_Strain %in% rs) %>%
  mutate(
    Ref_Strain=factor(Ref_Strain, levels=c("VI75.HuVac","BE92.HuVac","HK14.HuVac"))
  )

mdata<-mdata %>%
  mutate(
    Ref_Strain=factor(Ref_Strain, levels=c("VI75.HuVac","BE92.HuVac","HK14.HuVac"))
  )

(p<-mdata %>% 
  ggplot(., aes(x=log(Ref_Titor/10,2),y=log(Swine_Titor/10,2),color=Decade)) +
  #geom_smooth(method="lm",se=FALSE,color="darkgray",size=0.5,fullrange=TRUE)+
  geom_jitter() +
  #facet_grid(Swine_Strain~Ref_Strain, switch="both") +
  geom_text(data=text_corr,aes(x=2,y=7,label=round(Corr,digits = 2)),color="black")+
  facet_grid(Swine_Strain~Ref_Strain) +
  theme_minimal()+
  labs(title="Swine HI titers against HuVac HI titers",x="GMT HuVac",y="GMT swine")+
  theme(strip.text.y = element_text(angle = 360),
        panel.border = element_rect(colour = "black",fill=NA)
        )+
  coord_fixed())

ggsave("Post_vac_correlation.png",plot=p,width=6.5, height=9)

(p<-mdata %>% 
#    subset(Ref_Titor>200) %>%
    ggplot(., aes(x=log(Ref_Titor/10,2),y=log(Swine_Titor/10,2),color=Decade)) +
    #geom_smooth(method="lm",se=FALSE,color="darkgray",size=0.5,fullrange=TRUE)+
    geom_jitter() +
    #facet_grid(Swine_Strain~Ref_Strain, switch="both") +
    geom_text(data=text_corr,aes(x=2,y=7,label=round(Corr,digits = 2)),color="black")+
    facet_grid(Swine_Strain~Ref_Strain) +
    theme_minimal()+
    labs(title="Swine HI titers against HuVac HI titers",x="GMT HuVac",y="GMT swine")+
    theme(strip.text.y = element_text(angle = 360),
          panel.border = element_rect(colour = "black",fill=NA)
    ))
  #+
    #coord_fixed())

ggsave("Post_vac_correlation.png",plot=p,width=6.5, height=9)

# 3) Group the male and female together in a similar plot distributed by age by decade 
# for the post-vacc cohort.

input <- data2 %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                     "pre_vaccine","post_vaccine2","Fold_change")) %>%
  mutate(
    Swine_Strain=variable,
    Swine_Titor=value
  )
#%>% 
#  subset(.,!grepl("HuVac",.$variable)) %>% 
#  mutate(
#    Swine_Strain=variable,
#    Swine_Titor=value,
#    Ref_Strain="HK14.HuVac",
#    Ref_Titor=HK14.HuVac,
#    HK14.HuVac=NULL,
#    variable=NULL,
#    value=NULL
#  ) 

input %>%
  ggplot(.,aes(x=factor(Swine_Strain),y=log(Swine_Titor/10,2),color=factor(Swine_Strain),fill=factor(Swine_Strain))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=FALSE,width=1)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by Antigen",x="",y="GMT",fill="Decade",color="Decade")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none"
  ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))

input$vac = grepl("HuVac",input$Swine_Strain)

(p_all<-input %>%
  ggplot(.,aes(x=factor(Swine_Strain),y=log(Swine_Titor/10,2)))+ #,color=factor(Swine_Strain),fill=factor(Swine_Strain))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=FALSE,width=1,color="black",fill="white")+
  geom_boxplot(color="black",width=0.25,aes(fill=vac))+
  scale_fill_manual(values=c("TRUE"="#969696","FALSE"="#FFFFFF"))+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by strain",x="",y="GMT",fill="Decade",color="Decade")+
  theme(strip.text.y = element_text(angle = 180),
        theme(axis.text.x = element_text(size=5),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none"
  ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))
))
ggsave("Post_vac_titersall_huvac_violin plots.png",dpi=300,width=16,height=4.88)

input$vac2 <- "Strain"
input$vac2[input$vac] <- "HuVac"
annot_df <- input %>%
  group_by(Swine_Strain) %>%
  summarise(
    average = round (mean(Swine_Titor)),
    pos = sum(HI_to_GMT(Swine_Titor) >= 2)/40*100
  )
(p_all2<-input %>%
    ggplot(.,aes(x=factor(Swine_Strain),y=log(Swine_Titor/10,2)))+ #,color=factor(Swine_Strain),fill=factor(Swine_Strain))) +
    geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
    geom_violin(trim=FALSE,width=1)+
    geom_boxplot(color="black",aes(fill=vac), width=0.25)+
    scale_fill_manual(values=c("TRUE"="#969696","FALSE"="#FFFFFF"))+
    scale_color_manual(values=c("TRUE"="#969696","FALSE"="#000000"))+
    geom_jitter(color="black",width=0.1,height=0)+
    stat_summary(fun.y="mean", color="red", size=0.3)+ #red dots (mean)
    geom_text(data = annot_df, aes(label = average, 
                                   y=11))+
    geom_text(data = annot_df, aes(label = paste("(",
                                                 pos,
                                                 "%)",
                                                 sep=""
    ),
    y=10)
    )+ 
    facet_grid(.~Swine_Strain,switch="both",scale="free")+
    theme_minimal()+
    labs(title="HI Titers by strain",x="",y="Log2 GMT",fill="Decade",color="Decade")+
    theme(strip.text.y = element_text(angle = 180),
          panel.border = element_rect(colour = "black",fill=NA),
          legend.position = "none",
          text=element_text(size = 12)
    ) +
    scale_x_discrete(breaks=c(0:5)*2)+ylim(c(0,10))+
    scale_y_continuous(breaks=c(0:6)*2,limits=c(0,11))
)

(p_all_box<-input %>%
    ggplot(.,aes(x=factor(Swine_Strain),y=log(Swine_Titor/10,2)))+ #,color=factor(Swine_Strain),fill=factor(Swine_Strain))) +
    geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
#    geom_violin(trim=FALSE,width=1,aes(fill=vac,color=vac))+
    geom_boxplot(color="black",
                 width=0.25,aes(fill=vac))+
    scale_fill_manual(values=c("TRUE"="#969696","FALSE"="#FFFFFF"))+
    scale_color_manual(values=c("TRUE"="#969696","FALSE"="#000000"))+
    geom_jitter(color="black",width=0.1,height=0)+
    facet_grid(.~Swine_Strain,switch="both",scale="free")+
    theme_minimal()+
    labs(title="HI Titers by strain",x="",y="GMT",fill="Decade",color="Decade")+
    theme(strip.text.y = element_text(angle = 180),
          panel.border = element_rect(colour = "black",fill=NA),
          legend.position = "none"
    ) +
    scale_x_discrete(breaks=NULL)+ylim(c(0,8))
)

#Overall_PostVacc_ boxplots
input %>%
  ggplot(.,aes(x=factor(Swine_Strain),y=log(Swine_Titor/10,2),color=factor(Swine_Strain),fill=factor(Swine_Strain))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  #geom_violin(trim=FALSE,width=1)+
  geom_boxplot(color="black",fill= NA, width=0.25)+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by Antigen",x="",y="GMT",fill="Decade",color="Decade")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none"
  ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))
ggsave("Post_vac_titersall_huvac_boxplots.png",dpi=300,width=12,height=4.88)

input %>% subset(.,grepl("HuVac",Swine_Strain)) %>%
  ggplot(.,aes(x=factor(Decade),y=log(Swine_Titor/10,2),color=factor(Decade),fill=factor(Decade))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=FALSE,width=1)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by decade",x="",y="GMT",fill="Decade of birth",color="Decade of birth")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom"
        ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))
ggsave("Post_vac_titersbydecade_huvac.png",dpi=300,width=8.15,height=4.88)

input %>% subset(.,!grepl("HuVac",Swine_Strain)) %>%
  ggplot(.,aes(x=factor(Decade),y=log(Swine_Titor/10,2),color=factor(Decade),fill=factor(Decade))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=FALSE,width=1)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by year of birth",x="",y="GMT",fill="Decade of birth",color="Decade of birth")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom"
  ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))

ggsave("Post_vac_titersbydecade_swine.png",dpi=300,width=16,height=4.88)

# ==== This one

(p_birth =input %>% 
  ggplot(.,aes(x=factor(Decade),y=log(Swine_Titor/10,2)))+
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=FALSE,width=1, color="white", 
              aes(fill=factor(Decade))
              )+
  geom_jitter(width=0.1,height=0,alpha=0.9, size=2,color="black")+
  geom_boxplot( width=0.25,alpha=0,size=0.8,color="black")+
  facet_wrap(.~Swine_Strain,nrow=3)+
#  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  scale_color_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
  scale_fill_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
  labs(title="HI titers by decade of birth",x="",y="Log2 GMT",fill="Decade of birth",color="Decade of birth")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom", 
        text=element_text(size=16)
  ) +
#  guides(boxplot=FALSE)+
  scale_x_discrete(breaks=NULL)+ylim(c(0,8)))
#par("din")
# Boxplot
(p_decade_box<-input %>% 
  ggplot(.,aes(x=factor(Decade),y=log(Swine_Titor/10,2)))+
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
#  geom_violin(trim=FALSE,width=1, color="white", 
#              aes(fill=factor(Decade))
#  )+
  geom_jitter(width=0.1,height=0,alpha=0.9, size=2,color="black")+
  geom_boxplot( width=0.25,alpha=0,size=0.8,aes(color=factor(Decade)))+
  facet_wrap(.~Swine_Strain,nrow=3)+
  #  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  scale_color_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
  scale_fill_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
  labs(title="HI titers by decade of birth",x="",y="Log2 GMT",fill="Decade of birth",color="Decade of birth")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom"
  ) +
  #  guides(boxplot=FALSE)+
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))
)

# Violin_box
(p_decade_violinbox<-input %>% 
    ggplot(.,aes(x=factor(Decade),y=log(Swine_Titor/10,2)))+
    geom_hline(yintercept=2,linetype="dashed",color="#222222")+
    geom_violin(trim=FALSE, width=1, aes(color=factor(Decade),fill=factor(Decade)))+ 
    geom_jitter(width=0.1,height=0,alpha=0.9, size=2,color="black")+
    geom_boxplot( width=0.25,alpha=0,size=0.8,#aes(color=factor(Decade)),
                  color="black")+
    facet_wrap(.~Swine_Strain,nrow=3)+
    #  facet_grid(.~Swine_Strain,switch="both",scale="free")+
    theme_minimal()+
    scale_color_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
    scale_fill_manual(values=c("1940s and 50s"="#c6dbef","1960s"="#9ecae1","1970s"="#6baed6","1980s"="#3182bd","1990s"="#08519c"))+
    labs(title="HI titers by decade of birth",x="",y="Log2 GMT",fill="Decade of birth",color="Decade of birth")+
    theme(strip.text.y = element_text(angle = 180),
          panel.border = element_rect(colour = "black",fill=NA),
          legend.position = "bottom"
    ) +
    #  guides(boxplot=FALSE)+
    scale_x_discrete(breaks=NULL)+ylim(c(0,8))
)

library(cowplot)

#plot_grid(p_all, p_decade_box, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(1.5,3))
#ggsave("Post_vac_cowplot_boxplot.png",dpi=300,width=8.4,height=8)

plot_grid(p_all_box, p_decade_box, labels = c('A', 'B'), label_size = 10, ncol=1, rel_heights = c(1.5,3))
ggsave("Post_vac_cowplot_boxplot_v9.png",dpi=300,width=9.0,height=12)

#plot_grid(p_all, p_decade_violinbox, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(1.5,3))
#ggsave("Post_vac_cowplot_violinbox.png",dpi=300,width=8.4,height=8)

plot_grid(p_all2, p_birth, labels = c('A', 'B'), label_size = 14, ncol=1, rel_heights = c(4,9))
ggsave("Post_vac_cowplot_violinbox_v9.png",dpi=300,width=11.29,height=12)
#?plot_grid
par("din")
#ggsave("Post_vac_titersbydecade_both.png",dpi=300)
#ggsave("Post_vac_titersbydecade_both.png",dpi=300,width=7,height=9)

low_responders= input %>% subset(Swine_Strain=="HK14.HuVac" & Swine_Titor<=160) %>% {.$Study_ID}
#med_responders= input %>% subset(Swine_Strain=="HK14.HuVac" & Swine_Titor>=160 & Swine_Titor <=320) %>% {.$Study_ID}
high_responders = input %>% subset(Swine_Strain=="HK14.HuVac" & Swine_Titor>160) %>% {.$Study_ID}

input %>% mutate(
  responder=case_when(Study_ID %in% low_responders ~ "low responder (HI = 20 - 160)",
                      #Study_ID %in% med_responders ~ "med responder (HI = 160 - 320)",
                      Study_ID %in% high_responders ~ "high responder (HI = 320 - 640)"),
  responder=factor(responder, levels=c("low responder (HI = 20 - 160)",
                                       #"med responder (HI = 160 - 320)", 
                                       "high responder (HI = 320 - 640)"))
) %>%
  ggplot(.,aes(x=factor(responder),y=log(Swine_Titor/10,2),color=factor(responder),fill=factor(responder))) +
  geom_hline(yintercept=2,linetype="dashed",color="darkgray")+
  geom_violin(trim=TRUE,width=1)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(color="black",width=0.1,height=0)+
  facet_grid(.~Swine_Strain,switch="both",scale="free")+
  theme_minimal()+
  labs(title="HI Titers by Antigen",x="",y="GMT",fill="",color="")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "bottom"
  ) +
  scale_x_discrete(breaks=NULL)+ylim(c(0,8))
ggsave("Post_vac_titersresponder_huvac.png",dpi=300,width=12,height=4.88)

##----- Fold Change post vaccine
data2FC <- read_excel(post_vac.xlsx, sheet=post_vac.sheetFC)

#tdata <- data[,9:17] %>% t(.) %>% data.frame(.)

mdata2FC <- data2FC %>% melt(.,id=c("Study_ID","Age","Decade","Gender","Study_ID2",
                              "pre_vaccine","post_vaccine2","Fold_change"))
mdata2FC %>%
   mutate(
     variable=factor(variable, levels=c("HK14.HuVac","VI75.HuVac","BE92.HuVac",
                                      "NY11.1990.4a","IN16.1990.4b","MX15.1990.2",
                                      "AB14.1990.4","OH17.2010.1","OK17.2010.2"))
   ) %>%
  ggplot(.,aes(x=variable,y=value,color=variable,fill=variable)) +
  geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
  geom_violin(trim=TRUE,width=1)+
  geom_boxplot(color="black",width=0.25)+
  geom_jitter(height = 0, width=0.2,color="black")+
  #scale_y_continuous(breaks=c(-10:6),labels=c(-10:6))+
  scale_y_continuous(breaks=c(-5:5),limits=c(-5,5))+
  labs(y="Fold Change",title="Fold change against HK14.HuVac")+
  theme_minimal()+
  facet_grid(.~variable,switch="both",scale="free")+
  theme(strip.text.y = element_text(angle = 180),
        panel.border = element_rect(colour = "black",fill=NA),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  scale_x_discrete(breaks=NULL)
  
ggsave("Post_vaccine_cohort_foldchange.png",dpi=300,width=12.3,height=4.88)
  
