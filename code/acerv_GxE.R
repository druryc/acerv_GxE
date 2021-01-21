library(tidyverse);library(readxl);library(janitor)
library(lubridate)
library(scales)
library(caret)
library(ggrepel)
library(matrixStats)
library(doMC)
library(pdp)
library(cowplot)
library(zoo)
library(gbm)
library(factoextra)
library(RccpCNPy)
library(vegan)
setwd("~/CRD_GBS/3FD_Acerv/data")
############## MAP ############################################################################################ #####
library(rgdal);library(ggsn);library(tidyverse);library(sf);library(ggrepel)
florida<-st_read("../map/cb_2018_us_state_500k.shp")
sites<-read_tsv("../map/sites.txt")%>%select(-tag)%>%mutate(site=case_when(site=="MB"~"Miami Beach",
                                                                           site=="CVFD"~"Cheetos",
                                                                           site=="StruggleBus"~"Struggle Bus",
                                                                           site=="Site211"~"Site 211",
                                                                           site=="GovtCut"~"Govt Cut",
                                                                           TRUE ~ as.character(site)))

quartz(w=2.5,h=3.5)
c<-ggplot()+
     geom_sf(data=florida)+
     coord_sf(xlim=c(-80.5,-79.95),ylim=c(25.2,25.9))+
     geom_point(aes(long,lat,color=Type,size=Type),data=sites)+
     geom_label_repel(aes(long,lat,label=site),vjust=1,size=2,nudge_x = 0.2,nudge_y=0,data=sites%>%filter(Type=="Collection"))+
     scale_size_manual(values=c(5,3,2))+
     scale_color_manual(values=c("blue","lightgreen","orange"))+
     theme_minimal(base_size=8)+
     scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks=seq(24,27,0.1))+
     scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks=seq(-83.5,-79.8,0.1))+
     theme(legend.position=c(0.2,0.88),
           legend.title=element_blank(),
           legend.key.width = unit(0.25,"cm"),
           legend.key.height = unit(.6,"cm"),
           panel.border = element_rect(fill = "NA",colour = "gray35",size=.9))+
     xlab("Longitude (W)")+
     ylab("Latitude (N)")+
     #annotate("text",x=-80.4,y=25.9,label="Florida")+
     scalebar(x.min=-80.4,x.max=-80.3,y.min=25.7,y.max=25.8,transform=TRUE,model="WGS84",dist=10,dist_unit="km",
              height = 0.1, st.dist = 0.2,
              box.fill = c("black", "white"),
              box.color = "gray35", border.size = .5,st.size=2);c



############## TEMPERATURE TIMELINE ########################################################################### #####
rawdata<-read_xlsx("raw_data.xlsx",sheet="Temp")
daily<-rawdata%>%
     mutate(Time=ymd_hms(Time))%>%
     mutate(mean=rowMeans(.[,2:9],na.rm=TRUE))%>%
     gather(Site,Temp,-Time)%>%filter(Temp!="NA")%>%
     filter(Time >=as_datetime("2015-05-1"))%>%
     filter(Time <=as_datetime("2015-09-30"))%>%
     group_by(Site,Time=floor_date(Time,"day"))%>%summarise(Temp=mean(Temp))%>%
     mutate(group=case_when(Site=="mean" ~ "Mean",
                            Site!="mean" ~ "Outplant Sites"))


temp<-ggplot(daily)+
     geom_hline(yintercept=30.5,color="black",linetype="dotted")+
     geom_line(aes(Time,Temp,group=Site,color=group,size=group))+
     geom_line(data = subset(daily, Site == 'mean'),aes(Time,Temp,group=Site,color=group,size=group))+
     scale_color_manual(values=c("orange","darkgray"))+
     scale_size_manual(values=c(0.5,0.1))+
     theme_classic(base_size=7)+
     scale_y_continuous(limits=c(24,33),breaks=seq(24,33,1))+
     theme(axis.title.x=element_blank(),
           legend.position="bottom",
           legend.title=element_blank(),
           legend.key.size = unit(0.2,"line"))+
     ylab("Temperature (°C)")+
     scale_x_datetime(breaks = date_breaks("month"),labels = date_format("%b\n%Y"));temp

detail<-rawdata%>%
     mutate(Time=ymd_hms(Time))%>%
     mutate(mean=rowMeans(.[,2:9],na.rm=TRUE))%>%
     gather(Site,Temp,-Time)%>%filter(Temp!="NA")%>%
     filter(Time >=as_datetime("2015-05-1"))%>%
     filter(Time <=as_datetime("2015-09-30"))%>%
     filter(Site!="mean")

DHW<-detail%>%mutate(DHH=Temp-29.7)%>%filter(DHH>0)%>%group_by(Site)%>%summarise(DHH=sum(DHH))%>%mutate(DHH/(24*7))%>%rename(DHW=3)%>%
     mutate(Site=case_when(Site=="CVFD"~"Cheetos",Site=="MB"~"Miami Beach",Site=="StruggleBus"~ "Struggle Bus", TRUE ~ as.character(Site)))
DHWmean<-daily%>%filter(group=="Mean")%>%mutate(DHD=Temp-29.7)%>%mutate(DHW=DHD/7)%>%filter(DHD>0)%>%summarise(sum=sum(DHW))

############## TEMPERATURE TABLE ############################################################################## #####

#temp table
a<-detail%>%group_by(Site,Time=floor_date(Time,"hour"))%>%summarise(Temp=mean(Temp))%>%ungroup()%>%group_by(Site)%>%
     summarise(max=max(Temp),min=min(Temp),sd=sd(Temp),Average=mean(Temp),hr30_5=count(Temp>30.5),hr31=count(Temp>31),hr32=count(Temp>32),hr33=count(Temp>33))
b<-detail%>%mutate(date=date(Time))%>%group_by(Site,date)%>%summarise(min=min(Temp),max=max(Temp),range=abs(min-max))%>%group_by(Site)%>%summarise(mean=mean(range))
depths<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(site,Depth)%>%
     mutate(site=case_when(site=="CVFD"~"Cheetos",site=="MB"~"Miami Beach",site=="StruggleBus"~ "Struggle Bus", TRUE ~ as.character(site)))%>%rename(Site=1)
table1<-bind_cols(a,b$mean)%>%rename(Range=10,Minimum=min,Maximum=max)%>%
     select(Site,Average,sd,Range,Maximum,Minimum,Average,hr30_5,hr31,hr32,hr33)%>%
     mutate(Site=case_when(Site=="CVFD"~"Cheetos",Site=="MB"~"Miami Beach",Site=="StruggleBus"~ "Struggle Bus",TRUE ~ as.character(Site)))%>%
     left_join(.,depths,by="Site")%>%left_join(.,DHW,by="Site")%>%select(-DHH)
     
#write.table(table1,"../T1.txt",sep="\t",quote=FALSE,row.names=FALSE)

#supplemental temp table
a<-detail%>%filter(Time<"2015-6-15 00:00:00")%>%group_by(Site,Time=floor_date(Time,"hour"))%>%summarise(Temp=mean(Temp))%>%ungroup()%>%group_by(Site)%>%
     summarise(max=max(Temp),min=min(Temp),sd=sd(Temp),Average=mean(Temp),hr30_5=count(Temp>30.5),hr31=count(Temp>31),hr32=count(Temp>32),hr33=count(Temp>33))
b<-detail%>%mutate(date=date(Time))%>%group_by(Site,date)%>%summarise(min=min(Temp),max=max(Temp),range=abs(min-max))%>%group_by(Site)%>%summarise(mean=mean(range))
depths<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(site,Depth)%>%
     mutate(site=case_when(site=="CVFD"~"Cheetos",site=="MB"~"Miami Beach",site=="StruggleBus"~ "Struggle Bus", TRUE ~ as.character(site)))%>%rename(Site=1)
tables1<-bind_cols(a,b$mean)%>%rename(Range=10,Minimum=min,Maximum=max)%>%
     select(Site,Average,sd,Range,Maximum,Minimum,Average,hr30_5,hr31,hr32,hr33)%>%
     mutate(Site=case_when(Site=="CVFD"~"Cheetos",Site=="MB"~"Miami Beach",Site=="StruggleBus"~ "Struggle Bus",TRUE ~ as.character(Site)))%>%
     left_join(.,depths,by="Site")%>%left_join(.,DHW,by="Site")%>%select(-DHH,-DHW,-Depth,-contains('hr'))
#write.table(tables1,"../ST1.txt",sep="\t",quote=FALSE,row.names=FALSE)


############## TEMPERATURE PCA ################################################################################ #####


x<-as.data.frame(t(daily%>%filter(group!="Mean")%>%filter(Time<=as_datetime("2015-08-15"))%>%select(-group)%>%spread(Site,Temp))%>%row_to_names(1))
y<-mutate_if(x,is.factor,as.character)%>%mutate_if(is.character,as.numeric)

pca<-prcomp(y)
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(rownames(x),axes$data)%>%rename(site=1)

raw_bl<-read_xlsx("raw_data.xlsx",sheet="Condition")
bl_data<-raw_bl%>%
     filter(Tag!="NA")%>%
     select(Tag,Site,Genet,T1,T2,T3)%>%
     mutate(T1=case_when(T1=="P"~1,T1=="PB"~2,T1=="3"~3,is.na(T1)~0,TRUE ~ as.numeric(T1)))%>%
     mutate(T2=case_when(T2=="P"~1,T2=="PB"~2,T2=="3"~3,is.na(T2)~0,TRUE ~ as.numeric(T2)))%>%
     mutate(T3=case_when(T3=="P"~1,T3=="PB"~2,T3=="3"~3,is.na(T3)~0,TRUE ~ as.numeric(T3)))%>%
     filter(Genet!="NA")%>%filter(Tag!="N/A")%>%filter(T1!=3)%>%
     mutate(score=rowMeans(.[,4:6]))%>%select(-T1,-T2,-T3)%>%
     clean_names()%>%
     filter(genet!="Control1")%>%filter(genet!="Control2")%>%
     group_by(genet,site)%>%summarise(mean=mean(score))%>%ungroup()%>%
     group_by(site)%>%mutate(site_mean=mean(mean))%>%select(site,site_mean)%>%distinct()

final_plotdata<-left_join(plotdata,bl_data,by="site")%>%
     mutate(site=case_when(site=="MB"~"Miami Beach",site=="CVFD"~"Cheetos",site=="StruggleBus"~"Struggle Bus",
                           TRUE ~as.character(site)))%>%left_join(.,(DHW%>%rename(site=Site)),by="site")

quartz(w=3.5,h=3)
a<-ggplot(final_plotdata)+
     geom_hline(yintercept=0,linetype="dashed",color="gray")+
     geom_vline(xintercept=0,linetype="dashed",color="gray")+
     geom_point(aes(x,y,fill=site_mean,size=DHW),pch=21,)+
     geom_text_repel(aes(x,y,label=site),hjust=1, vjust=1,size=3,  nudge_x = -0.5,nudge_y=0.5)+
     scale_fill_gradient(low="blue",high = "red",name="Mean\nBleaching\nScore")+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-4,3))+
     scale_x_continuous(limits=c(-7,5))+
     ylab("PC2 (20.5%)")+xlab("PC1 (54.3%)")+
     theme(legend.key.size = unit(0.25,"cm"));a




############## PHENOTYPE  DATA ################################################################################ #####
raw_bl<-read_xlsx("raw_data.xlsx",sheet="Condition")
bl_data<-raw_bl%>%
     filter(Tag!="NA")%>%
     select(Tag,Site,Genet,T1,T2,T3)%>%
     mutate(T1=case_when(T1=="P"~1,T1=="PB"~2,T1=="3"~3,is.na(T1)~0,TRUE ~ as.numeric(T1)))%>%
     mutate(T2=case_when(T2=="P"~1,T2=="PB"~2,T2=="3"~3,is.na(T2)~0,TRUE ~ as.numeric(T2)))%>%
     mutate(T3=case_when(T3=="P"~1,T3=="PB"~2,T3=="3"~3,is.na(T3)~0,TRUE ~ as.numeric(T3)))%>%
     filter(Genet!="NA")%>%filter(Tag!="N/A")%>%filter(T1!=3)%>%
     mutate(score=rowMeans(.[,4:6]))%>%select(-T1,-T2,-T3)%>%
     clean_names()%>%
     filter(genet!="Control1")%>%filter(genet!="Control2")%>%
     group_by(genet,site)%>%summarise(mean=mean(score))%>%ungroup()%>%
     group_by(site)%>%mutate(site_mean=mean(mean))%>%ungroup()%>%
     mutate(residual=mean-site_mean)
res<-bl_data%>%group_by(genet)%>%summarise(mean=mean(residual))%>%arrange(-mean)
code<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(tag,site)%>%rename(genet=2)
residual<-inner_join(res,code,by="genet")%>%select(tag,genet,mean)%>%rename(residual=mean)%>%rename(sample=tag)
#saveRDS(residual,"residual")

bl_data$genet<-factor(bl_data$genet,levels=c("Coopers","GovtCut","MB","CVFD","Stephs","Jons","Site211","StruggleBus","Grounding","Inshore"),labels=c("Coopers","Govt Cut","Miami Beach","Cheetos","Stephs","Jons","Site 211","Struggle Bus","Grounding","Inshore"))
residplot<-ggplot(bl_data)+geom_hline(yintercept=0,color="gray",linetype="dashed")+
     geom_boxplot(aes(reorder(genet,-residual),residual,fill=genet),outlier.color=NA)+
     geom_point(aes(reorder(genet,residual),residual),color="black",size=0.5)+
     scale_fill_viridis_d()+
     theme_classic(base_size=7)+
     theme(legend.position="none",
           axis.text.x= element_text(angle = 90,size=6,hjust=1))+
     ylab("Bleaching Residual")+
     xlab("Genotype")+
     annotate("text", x = 7.4, y = -0.6, label = "Less Bleaching",hjust=0,size=2,fontface="italic")+
     annotate("text", x = 7.4, y = 1, label = "More Bleaching",hjust=0,size=2,fontface="italic")+
     annotate("text", x= 0.75, y=-0.6,label="Genotype p=0.012",hjust=0,size=2,fontface="italic")

bl_data$site<-factor(bl_data$site,levels=c("Coopers","MB","CVFD","Stephs","Jons","StruggleBus","Grounding","Inshore"),labels=c("Coopers","Miami Beach","Cheetos","Stephs","Jons","Struggle Bus","Grounding","Inshore"))
grid<-ggplot(bl_data)+geom_tile(aes(genet,site,fill=mean))+theme_classic(base_size=7)+
     theme(legend.key.size = unit(0.4,"line"),
           axis.text.x= element_text(angle = 90,size=6,hjust=1,vjust=1))+
     scale_fill_gradient(low="blue",high = "red",name="Mean\nScore")+
     annotate("text",x=1,y=5,label="ND",size=2)+
     annotate("point",x=1,y=1,size=0.5)+
     annotate("point",x=3,y=2,size=0.5)+
     annotate("point",x=4,y=3,size=0.5)+
     annotate("point",x=5,y=4,size=0.5)+
     annotate("point",x=6,y=5,size=0.5)+
     annotate("point",x=8,y=6,size=0.5)+
     annotate("point",x=9,y=7,size=0.5)+
     annotate("point",x=10,y=8,size=0.5)+
     xlab("Genotype")+ylab("Site")
     
quartz(w=7.2,h=2.25)
plot_grid(temp,grid,residplot,nrow=1,align="h", axis="b",labels=c("A","B","C"),label_size=8,rel_widths=c(0.8,1.1,0.8),label_x=c(0,0,-0.05))
############## SURVIVORSHIP ################################################################################### #####
raw_bl<-read_xlsx("raw_data.xlsx",sheet="Mortality")
bl_data<-raw_bl%>%
     filter(Tag!="NA")%>%
     select(Tag,Site,Genet,T1,T2,T3,T4)%>%
     mutate(T1=case_when(T1=="P"~1,T1=="PB"~2,T1=="3"~3,is.na(T1)~0,TRUE ~ as.numeric(T1)))%>%
     mutate(T2=case_when(T2=="P"~1,T2=="PB"~2,T2=="3"~3,is.na(T2)~0,TRUE ~ as.numeric(T2)))%>%
     mutate(T3=case_when(T3=="P"~1,T3=="PB"~2,T3=="3"~3,is.na(T3)~0,TRUE ~ as.numeric(T3)))%>%
     mutate(T4=case_when(T4=="H"~"1",T4==0~"0"))%>%
     mutate(T4=as.numeric(T4))%>%
     filter(Genet!="NA")%>%filter(Tag!="N/A")%>%
     filter(T1!=3)%>%
     mutate(score=rowMeans(.[,4:6]))%>%select(-T1,-T2,-T3)%>%
     clean_names()%>%
     filter(genet!="Control1")%>%filter(genet!="Control2")

points<-bl_data%>%mutate(t4=case_when(t4==1~"Alive",t4=="0"~"Dead"))
a<-ggplot(points)+geom_boxplot(aes(t4,score,fill=t4))+
     geom_jitter(aes(t4,score),height=0.01,width=0.25,alpha=0.2)+
     theme_classic(base_size=8)+
     theme(legend.key.size = unit(0.35,"cm"),legend.position="none")+
     scale_fill_manual(values=c("gray","orange"))+
     ylab("Visual Bleaching Score")+xlab("Final Survivorship")+
     annotate('text',x=0.5,y=2,label=paste("italic(Wilcox~p==0.001)"), size=2,hjust=0,parse=TRUE);a

wilcox.test(score~t4,data=points)

averages<-bl_data%>%group_by(site,genet)%>%summarise(score=mean(score),t4=mean(t4))
b<-ggplot(averages)+geom_point(aes(score,t4))+theme_classic(base_size=8)+
     xlab("Visual Bleaching Score")+ylab("Proportion Survivorship")+
     scale_x_continuous(limits=c(0,2),breaks=seq(0,2,0.5))+
     scale_y_continuous(limits=c(-0.1,1),breaks=seq(0,1,0.25));b

summary(lm(score~t4,data=averages%>%filter(t4>0.01)))

c<-ggplot(averages%>%filter(t4>0.01))+geom_point(aes(score,t4))+geom_smooth(aes(score,t4),method="lm")+theme_classic(base_size=8)+
     xlab("Visual Bleaching Score")+ylab("Proportion Survivorship")+
     scale_x_continuous(limits=c(0,2),breaks=seq(0,2,0.5))+
     scale_y_continuous(limits=c(-0.1,1),breaks=seq(0,1,0.25))+
     annotate('text',x=0,y=-0.1,label=paste("italic(lm~p<0.001)","~italic(R)^2==0.35"), size=2,hjust=0,parse=TRUE);c


quartz(w=5.5,h=2)
plot_grid(a,b,c,nrow=1,labels=c("A","B","C"),label_size=8)



############## BLEACHING STATS ################################################################################ #####
raw_bl<-read_xlsx("raw_data.xlsx",sheet="Condition")
int<-raw_bl%>%
     filter(Tag!="NA")%>%
     select(Tag,Site,Genet,T1,T2,T3)%>%
     mutate(T1=case_when(T1=="P"~1,T1=="PB"~2,T1=="dead"~3,is.na(T1)~0,TRUE ~ as.numeric(T1)))%>%
     mutate(T2=case_when(T2=="P"~1,T2=="PB"~2,T2=="dead"~3,is.na(T2)~0,TRUE ~ as.numeric(T2)))%>%
     mutate(T3=case_when(T3=="P"~1,T3=="PB"~2,T3=="dead"~3,is.na(T3)~0,TRUE ~ as.numeric(T3)))%>%
     filter(Genet!="NA")%>%filter(Tag!="N/A")%>%filter(T1!=3)%>%
     mutate(score=rowMeans(.[,4:6]))%>%select(-T1,-T2,-T3)%>%
     clean_names()%>%
     filter(genet!="Control1")%>%filter(genet!="Control2")%>%
     mutate(site=as.factor(site),genet=as.factor(genet))
site_mean<-int%>%group_by(site)%>%summarise(site_mean=mean(score))
data<-inner_join(int,site_mean,by="site")%>%mutate(residual=score-site_mean)

library(car)
library(MASS)
### MODEL1: score by site*genotype
model1<-glm(sqrt(score)~genet*site,data=data) 
qqPlot(residuals(model1))
leveneTest(residuals(model1)~data$genet*data$site) #homogeneity of variance
Anova(model1,test.statistic=c("F"))

### MODEL2: residual by genotype
model2<-glm((residual)~genet,data=data) 
qqPlot(residuals(model2))
leveneTest(residuals(model2)~data$genet) #homogeneity of variance
Anova(model2,test.statistic=c("F"))

detach(package:MASS);detach(package:car)
############## SYMBIONTS ###################################################################################### #####
code<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(tag,site)%>%rename(genet=2)
list<-as.data.frame(rep(c("A","B","C","D"),10))%>%rename(symbiont=1)
genets<-as.data.frame(rep(c("C","D","E","F","G","H","I","J","K","L"),4))%>%rename(tag=1)%>%arrange(tag)
metadata<-left_join(bind_cols(genets,list),code,by="tag")
counts<-read_tsv("./symbionts/symbiont_summary.txt",col_names = FALSE,n_max=40)%>%rename(count=1)

data<-bind_cols(metadata,counts)%>%select(-tag)%>%
     group_by(genet)%>%
     mutate(n=sum(count))%>%
     mutate(prop=count/n)%>%
     select(-count,-n)

data$symbiont <- factor(data$symbiont,levels=c("A","B","C","D"),labels = c("Symbiodinium","Breviolum","Cladocopium","Durusdinium"))
data$genet<-factor(data$genet,levels=c("Coopers","GovtCut","MB","CVFD","Stephs","Jons","Site211","StruggleBus","Grounding","Inshore"),labels=c("Coopers","Govt Cut","Miami Beach","Cheetos","Stephs","Jons","Site 211","Struggle Bus","Grounding","Inshore"))
quartz(w=3,h=2.5)
symb<-ggplot(data)+geom_bar(aes(genet,prop,fill=symbiont),stat="identity",width=0.7)+
     theme_classic(base_size=7)+
     theme(legend.position ="bottom",
           legend.title=element_blank(),
           legend.text = element_text(face="italic"),
           legend.key.size = unit(0.5,"line"),
           axis.text.x= element_text(angle = 90,size=6,hjust=1))+
     scale_fill_manual(values=c("orange","blue","red","gray"))+
     ylab("Proportion Reads Aligned" )+xlab("Genotype")+
     guides(fill=guide_legend(nrow=1));symb
     
############## IBS ############################################################################################ #####
matrix<-as.matrix(read_tsv("./angsd/genolike.ibsMat",col_names = FALSE)%>%dplyr::select(-X11))
diag(matrix)<-NA
#matrix[upper.tri(matrix)] <-NA
rownames(matrix)<-c("C","D","E","F","G","H","I","J","K","L")
colnames(matrix)<-c("C","D","E","F","G","H","I","J","K","L")
code<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(tag,site)%>%rename(genet_1=1)%>%mutate(genet_2=genet_1)
code[1,2]<-"Miami Beach"
code[7,2]<-"Cheetos"

pwdgrid<-as.data.frame(matrix)%>%rownames_to_column(var="genet_1")%>%gather(genet_2,PWD,-genet_1)%>%left_join(.,code,by="genet_1")%>%select(-genet_2.y)%>%
     rename(genet_2=2)%>%left_join(.,code,by="genet_2")
pwdgrid$site.x<-factor(pwdgrid$site.x,levels=c("Coopers","GovtCut","Miami Beach","Cheetos","Stephs","Jons","Site211","StruggleBus","Grounding","Inshore"),labels=c("Coopers","Govt Cut","Miami Beach","Cheetos","Stephs","Jons","Site 211","Struggle Bus","Grounding","Inshore"))
pwdgrid$site.y<-factor(pwdgrid$site.y,levels=c("Coopers","GovtCut","Miami Beach","Cheetos","Stephs","Jons","Site211","StruggleBus","Grounding","Inshore"),labels=c("Coopers","Govt Cut","Miami Beach","Cheetos","Stephs","Jons","Site 211","Struggle Bus","Grounding","Inshore"))


quartz()
a<-ggplot(pwdgrid)+geom_tile(aes(site.x,site.y,fill=PWD))+theme_classic(base_size=8)+
     scale_fill_gradient(na.value = 'red',low="white",high="lightgray")+
     theme(axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.text.x= element_text(angle = 90,hjust=1,vjust=1),
           legend.key.size = unit(0.25,"cm"),
           legend.position="left");a

############## PCANGSD  ####################################################################################### #####
data<-as.matrix(read.table("./angsd/pcangsd_out.cov"))
labels<-read_tsv("./angsd/bamlist.txt",col_names=FALSE)%>%rename(file=1)%>%separate(file,into=c("file","garbo"),sep="_")%>%select(-garbo)

e <- eigen(data)
axes<-as.data.frame(e$vectors[,1:2])%>%rename(PC1=1,PC2=2)
code<-read_excel("raw_data.xlsx",sheet="Genets")%>%select(tag,site)%>%filter(tag!="NA")
residual<-readRDS("residual")%>%select(-genet)%>%rename(tag=1)
plotdata<-bind_cols(axes,left_join(code,residual,by="tag"))%>%rename(genet=site)%>%
     mutate(genet=case_when(genet=="GovtCut"~"Govt Cut",
                            genet=="StruggleBus"~"Struggle Bus",
                            genet=="Site211"~"Site 211",
                            genet=="MB" ~ "Miami Beach",
                            genet=="CVFD" ~ "Cheetos",
                            TRUE~as.character(genet)))
summary(prcomp(data))
saveRDS(plotdata,"PCA_results")

quartz()
b<-ggplot(plotdata)+geom_point(aes(PC1,PC2,fill=residual),pch=21,size=3)+
     geom_text_repel(aes(PC1,PC2,label=genet),size=2,force=20)+
     theme_classic(base_size=8)+
     #scale_y_continuous(limits=c(-50,25))+
     #scale_x_continuous(limits=c(-50,50))+
     ylab("PC2 (15.9%)")+xlab("PC1 (15.9%)")+
     scale_fill_viridis_c(direction = -1,name="Bleaching\nResidual")+
     theme(legend.position="bottom",
           legend.key.height = unit(0.25,"cm"))+
     scale_size(guide = 'none')+
     annotate("text", x = -0.5, y = -0.9, label = "NGSadmix K=1",hjust=0,size=2,fontface="italic")

quartz(w=5.5,h=2.5)
plot_grid(a,b,align="h",axis="tb",rel_widths=c(2.3,2),labels=c("A","B"),label_size=8)


############## RESIDUAL REGRESSION ############################################################################ #####
residual<-readRDS('residual')
depth<-read_tsv("./angsd/genolike_depth.counts")%>%select(-X11)%>%mutate(mean=rowMeans(.[,1:10]))
loci<-read_tsv("./angsd/genolike_depth.beagle")%>%select(marker,allele1,allele2)
cols<-colnames(read_tsv("./angsd/genolike_depth.beagle")%>%select(marker,contains('Ind')))

rawdat<-read_tsv("./angsd/genolike.geno",col_names=FALSE)%>%select(-X33)%>%unite(marker,X1,X2,sep="_")
colnames(rawdat)<-cols

ex<-rawdat%>%select(marker,contains('_'))%>%
     rowwise()%>%
     mutate(C=Ind0_1+(2*Ind0_2))%>%select(-Ind0_1,-Ind0_2)%>%
     mutate(D=Ind1_1+(2*Ind1_2))%>%select(-Ind1_1,-Ind1_2)%>%
     mutate(E=Ind2_1+(2*Ind2_2))%>%select(-Ind2_1,-Ind2_2)%>%
     mutate(F=Ind3_1+(2*Ind3_2))%>%select(-Ind3_1,-Ind3_2)%>%
     mutate(G=Ind4_1+(2*Ind4_2))%>%select(-Ind4_1,-Ind4_2)%>%
     mutate(H=Ind5_1+(2*Ind5_2))%>%select(-Ind5_1,-Ind5_2)%>%
     mutate(I=Ind6_1+(2*Ind6_2))%>%select(-Ind6_1,-Ind6_2)%>%
     mutate(J=Ind7_1+(2*Ind7_2))%>%select(-Ind7_1,-Ind7_2)%>%
     mutate(K=Ind8_1+(2*Ind8_2))%>%select(-Ind8_1,-Ind8_2)%>%
     mutate(L=Ind9_1+(2*Ind9_2))%>%select(-Ind9_1,-Ind9_2)%>%
     gather(sample,GL,-marker)

list<-loci%>%select(-allele1,-allele2)%>%
     separate(marker,into=c("CHR","POS"),sep="_")%>%arrange(CHR,POS)%>%
     mutate(diff=as.numeric(POS)-lag(as.numeric(POS),n=1L))%>%unite(POS,CHR,POS,sep="_")%>%
     replace_na(list(diff=-1))%>%
     rename(marker=1)%>%filter(diff<0|diff>0)%>%select(-diff)

complete<-inner_join(ex,list%>%select(marker),by="marker")%>%
     spread(marker,GL)%>%
     inner_join(.,residual,by="sample")%>%
     select(sample,residual,everything(),-genet)

data<-complete[,-nearZeroVar(complete)]
input<-data[,2:ncol(data)]
col <- names(input)[-1]

lm.loop <- vector("list", length(col))
for(i in seq_along(col)){
     lm.loop[[i]] <- lm(reformulate(col[i],"residual"), data = input)
}
smry <- lapply(lm.loop, summary)
results <- data.frame(p_value=numeric(length(col)),r2=numeric(length(col)))
for (i in 1:(ncol(input)-1)){
     results$p_value[i]<-smry[[i]]$coefficients[2,4]
     results$r2[i]<-smry[[i]]$adj.r.squared
}

output<-bind_cols(as.data.frame(col),results)%>%rename(marker=1)%>%
     filter(p_value<0.01)%>%
     separate(marker,into=c("CHR","POS"),sep="_")%>%arrange(CHR,POS)%>%
     mutate(diff=as.numeric(POS)-lag(as.numeric(POS),n=1L))%>%unite(POS,CHR,POS,sep="_")%>%
     replace_na(list(diff=-1))%>%
     rename(marker=1)%>%filter(diff<0|diff>100)%>%select(-diff)
saveRDS(output,"top_individual_markers")

heats<-bind_cols(as.data.frame(col),results)%>%rename(marker=1)%>%
     separate(marker,into=c("CHR","POS"),sep="_")%>%arrange(CHR,POS)%>%
     mutate(diff=as.numeric(POS)-lag(as.numeric(POS),n=1L))%>%unite(POS,CHR,POS,sep="_")%>%
     replace_na(list(diff=-1))%>%
     rename(marker=1)%>%filter(diff<0|diff>0)%>%select(-diff)
saveRDS(heats,"raw_lm_output")

rf_test_data<-semi_join(ex,output,by="marker")%>%arrange(marker)%>%spread(marker,GL)%>%
     inner_join(.,residual,by="sample")%>%select(sample,residual,everything(),-genet)
saveRDS(rf_test_data,"rf_test_data")



############## IMPORTANT LOCI VS STRUCTURE #################################################################### #####
pca<-readRDS("PCA_results")%>%rename(sample=tag)
output<-readRDS("top_individual_markers")
cov_test<-semi_join(ex,output,by="marker")%>%arrange(marker)%>%inner_join(.,pca,by="sample")%>%select(marker,sample,GL,PC1)%>%
     spread(marker,GL)%>%select(-sample)

col<-names(cov_test)[-1]
lm.loop<-vector("list",length(names(cov_test))-1)
for(i in seq_along(col)){
     lm.loop[[i]] <- lm(reformulate(col[i],"PC1"), data = cov_test)
}
smry <- lapply(lm.loop, summary)
pca_results <- data.frame(p_value=numeric(length(col)),r2=numeric(length(col)))%>%rownames_to_column()%>%mutate(rowname=as.numeric(rowname)+4)
for (i in 1:(ncol(cov_test)-1)){
     pca_results$p_value[i]<-smry[[i]]$coefficients[2,4]
     pca_results$r2[i]<-smry[[i]]$adj.r.squared
}

quartz()
a<-ggplot(pca_results)+geom_density(aes(r2))+geom_rug(aes(r2))+
     scale_x_continuous(limits=c(-0.5,0.5),breaks=seq(-0.5,0.5,0.1))+
     theme_classic(base_size=8)+
     xlab("r2 of predicted derived allele count against PC1\nfor 58 loci highly associated with bleaching residual")
max(abs(pca_results$r2))

b<-ggplot(pca)+geom_point(aes(PC1,residual))+
     theme_classic(base_size=8)+
     ylab("Bleaching Residual")+
     annotate("text",y=0.2,x=-0.2,label="r2=0.0001, p=0.97",size=3)+
     geom_smooth(aes(PC1,residual),method="lm")
summary(lm(PC1~residual,data=pca))
plot_grid(a,b,align="h",axis="b",rel_widths=c(1.5,1))

############## TOP ANNOTATIONS - BLAST ######################################################################## #####
size<-1000
list<-readRDS("top_individual_markers")%>%
     select(marker)%>%rename(POS=1)%>%
     separate(POS,into=c("CHR","POS"),sep="_")%>%
     mutate(start=as.numeric(POS)-(size/2))%>%
     mutate(start=case_when(start<0 ~ 0,TRUE ~ as.numeric(start)))%>%
     mutate(end=as.numeric(POS)+(size/2))%>%
     unite(name,CHR,POS,sep="_",remove=FALSE)%>%select(CHR,start, end, name)
write.table(list,"./annotation_working/blast_list.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
system(paste("/anaconda3/envs/vcf/bin/bedtools getfasta -bed ./annotation_working/blast_list.bed -fi ~/CRD_GBS/genome/Amil_v2.02/Amil.v2.02.chrs.fasta -nameOnly -fo ./annotation_working/blast_list.fasta"))

#ncbi against taxid 6073, save as hit-table.csv
#https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
accessions<-read_csv("./annotation_working/ZZ3WEHG3016-Alignment-HitTable.csv",col_names=FALSE)%>%select(X1,X2,X11)%>%rename(POS=1,accession=2,e=3)%>%
     group_by(POS)%>%arrange(e)%>%slice(1)%>%filter(e<1e-4)

library(magicfor)
library(rentrez)
magic_for(print,silent=FALSE)
for (i in unique(accessions$accession)){
     print(entrez_summary(db="nucleotide", id=(entrez_search(db="nucleotide", term=i,use_history=TRUE)$ids))$title)
}

raw_lm_output<-readRDS("raw_lm_output")
blast_out<-magic_result_as_dataframe()%>%rename(accession=1,name=2)%>%select(accession,name)%>%left_join(.,accessions,by="accession")%>%
     select(POS,accession,name,e)%>%left_join(.,raw_lm_output%>%rename(POS=1),by="POS")%>%select(POS,p_value,r2,accession,e,name)%>%
     separate(name,into=c("name","trash"),sep="\\(")%>%select(-trash)%>%separate(name,into=c("garbo","name"),sep=":")%>%select(-garbo)%>%
     separate(name,into=c("garbo","species"),sep=9)%>%select(-garbo)%>%mutate(genus="A.")%>%unite(name,genus,species,sep=" ")

saveRDS(blast_out,"annotations")
detach("package:magicfor",unload=TRUE)
############## FIT ############################################################################################ #####
residual<-readRDS("residual")
rf_test_data<-readRDS("rf_test_data")%>%select(-sample)
set.seed(3839)
seeds <- vector(mode = "list", length = 61) #nrepeats * nfolds +1
for(i in 1:61) seeds[[i]]<- sample.int(n=100, 100) #last value is at least mtry
seeds[[61]]<-sample.int(1000, 1)
mtry <- sqrt(ncol(rf_test_data));tunegrid <- expand.grid(.mtry=seq(1,25,1))
modelFit <- train(residual~., 
                  data=rf_test_data,
                  method="rf",
                  preProcess=c("center","scale"),
                  trControl=trainControl(method="repeatedcv",repeats=20,number=2,seeds=seeds,savePredictions = TRUE),
                  metric="Rsquared",
                  nTree=1000,
                  tuneGrid=tunegrid, 
                  #tuneLength=100,
                  importance=TRUE,
                  verbose=FALSE);modelFit

saveRDS(modelFit,"modelFit")
#modelFit<-readRDS("modelFit")
plot(modelFit)
modelFit$bestTune
modelFit$resample$Rsquared
mean(modelFit$resample$Rsquared,na.rm=TRUE)
sd(modelFit$resample$Rsquared,na.rm=TRUE)
out<-cbind(residual,as.data.frame(predict(modelFit, newdata = rf_test_data)))%>%rename(prediction=4)
prediction<-predict(modelFit,newdata=rf_test_data,type=c("raw"))
rsq<-as.data.frame(modelFit$resample$Rsquared)%>%rename(r2=1)%>%filter(r2<1)%>%filter(r2!="NA")

#histogram of resampling accuracies
hist<-ggplot(rsq)+
     geom_histogram(aes(r2),color="NA",fill="orange",binwidth=0.01,boundary=0)+
     geom_density(aes(r2))+
     theme_classic(base_size=7)+ylab("Resamplings (n)")+
     xlab(expression("Variance explained (R"^2~") in CV"))+
     scale_y_continuous(breaks=seq(0,30,5))+
     scale_x_continuous(limits=c(0.75,1),breaks=seq(0.75,1,0.05));hist

#predictions and figures
out<-cbind(residual%>%arrange(sample),as.data.frame(predict(modelFit, newdata = rf_test_data)))%>%rename(prediction=4)%>%mutate(prediction=prediction)
summary(lm(residual~prediction,data=out))
anno1=expression("R"^2~"=0.964; p<0.001")
reg<-ggplot(out)+
     geom_abline(slope=1)+
     geom_hline(yintercept=0,linetype="dotted",color="gray")+
     geom_vline(xintercept=0,linetype="dotted",color="gray")+
     geom_smooth(aes(residual,prediction),method="lm",color="orange",fill="orange")+
     geom_point(aes(residual,prediction),size=0.5)+
     theme_classic(base_size=7)+
     ylab("Prediction")+xlab("Observed Bleaching Residual")+
     #scale_x_continuous(limits=c(-0.6,1.05))+
     #scale_y_continuous(limits=c(-0.6,1.05))+
     annotate("text",x=0.1,y=-0.15,label=anno1,size=2,fontface = 'italic');reg

quartz(w=4,h=1.75)
plot_grid(hist,reg,nrow=1,align="h",axis="b",labels=c("A","B"),label_size=8)


x<-test_data%>%gather(POS,freq,-residual)
ggplot(x[401:630,])+geom_point(aes(residual,freq))+
     geom_smooth(aes(residual,freq),method="lm")+
     theme_classic(base_size=4)+
     facet_wrap(~POS,nrow=1)



############## EFFECTS OF MISSING VARIANTS #################################################################### #####
residual<-readRDS('residual')%>%arrange(sample)
rf_test_data<-readRDS('rf_test_data')%>%select(-sample)
modelFit<-readRDS("modelFit")
count_na <- function(x) sum(is.na(x))

getValue <- function() {
     missing_data<-as.data.frame(lapply(rf_test_data[,2:67], function(cc) cc[ sample(c(TRUE, NA), prob = c(0.5, 0.5), size = length(cc), replace = TRUE) ]))
     missingness<-missing_data%>%mutate(count_na = apply(., 1, count_na))%>%select(count_na);mean(missingness$count_na)/65
     missing_data[is.na(missing_data)]<-1
     prediction<-predict(modelFit,newdata=missing_data,type=c("raw"))
     out<-cbind(residual,prediction)%>%rename(prediction=4)%>%mutate(prediction=2*prediction)
     ggplot(out)+geom_point(aes(prediction,residual))+geom_smooth(aes(prediction,residual),method="lm")
     report<-summary(lm(residual~prediction,data=out))
     return(report$r.squared)
}
missing_accuracy<-replicate(100, getValue())
mean(missing_accuracy)
sd(missing_accuracy)

############## VARIABLE IMPORTANCE ############################################################################ #####
lm_output<-readRDS("raw_lm_output")
blast_out<-readRDS("annotations")%>%rename(marker=POS)
modelFit<-readRDS("modelFit")
VI<-varImp(modelFit,scale=TRUE)$importance%>%rename(VI=1)%>%rownames_to_column()%>%rename(marker=1)
stable2<-left_join(VI,blast_out,by="marker")%>%select(-p_value,-r2)%>%left_join(lm_output,by="marker")%>%
     replace_na(list(accession="Unknown"))%>%separate(p_value,into=c("p_value","garbo"),sep=6)%>%select(-garbo)%>%
     mutate(R=sqrt(r2))%>%separate(R,into=c("R","garbo"),sep=6)%>%select(-r2,-garbo)%>%select(marker,p_value,R,VI,accession,name,e)%>%
     rename(Position=1,Accession=5,Gene=6,e_score=7)
write.table(stable2,"../revision/ST2.txt",sep="\t",quote=FALSE,row.names=FALSE)


############## WRITE SITE LIST FOR ADDITIONAL SAMPLES ######################################################### #####
rf_test_data<-readRDS("rf_test_data")%>%select(-residual,-sample)
size<-0
list<-as.data.frame(colnames(rf_test_data))%>%rename(POS=1)%>%
     separate(POS,into=c("CHR","POS"),sep="_")%>%
     mutate(start=as.numeric(POS)-(size/2))%>%
     mutate(start=case_when(start<0 ~ 0,TRUE ~ as.numeric(start)))%>%
     mutate(end=as.numeric(POS)+(size/2))%>%
     unite(name,CHR,POS,sep="_",remove=FALSE)%>%select(CHR,start, end, name)
write.table(list,"short_list_additional.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


############## SEQDATA - ADDITIONAL SAMPLES ################################################################### #####
cols<-colnames(read_tsv("./angsd/acerv_extra.beagle")%>%select(marker,contains('Ind')))
toplist<-readRDS("top_individual_markers")
samples<-read_tsv("./angsd/extra_bamlist.txt",col_names=FALSE)%>%separate(X1,into=c("sample","garbo"),sep=-4)%>%select(-garbo)%>%rowid_to_column()%>%
     mutate(rowid=as.numeric(rowid)-1)%>%mutate(Ind="Ind")%>%unite(name,Ind,rowid,sep="")%>%rename(sample=1,fullname=2)

rawdat<-read_tsv("./angsd/acerv_extra.geno",col_names=FALSE)%>%unite(marker,X1,X2,sep="_")%>%select(-X1833)%>%semi_join(.,toplist,by="marker")
colnames(rawdat)<-cols

match<-read_tsv("./angsd/acerv_extra.beagle")%>%select(marker,allele1,allele2)%>%inner_join(.,read_tsv("./angsd/genolike_depth.beagle")%>%select(marker,allele1,allele2),by="marker")%>%
     distinct()%>%rename(other1=2,other2=3,main1=4,main2=5)%>%mutate(match=case_when(other1==main1&other2==main2 ~"match",
                                                                                     other1==main2&other2==main1 ~"swap",
                                                                                     other1==main1&other2!=main2 ~"other"))

depth_data<-bind_cols(read_tsv("./angsd/acerv_extra.pos")%>%unite(pos,chr,pos,sep="_")%>%select(pos),read_tsv("./angsd/acerv_extra.counts"))%>%select(pos,everything(),-X611)%>%
     rename(marker=1)%>%semi_join(.,match,by="marker")%>%
     gather(ind,depth,-marker)%>%
     separate(ind,into=c("ind","trash"),sep=-8)%>%select(-trash)%>%
     mutate(sample = str_to_title(ind))

rf_data_extra<-rawdat%>%select(marker,contains('Ind'))%>%
     gather(sample,prob,-marker)%>%
     left_join(.,match,by="marker")%>%
     separate(sample,into=c("sample","rep"),sep="_")%>%
     replace_na(list(rep="0"))%>%
     mutate(rep=case_when(match=="swap"&rep==0 ~ "B",
                          match=="swap"&rep==2 ~ "A",
                          TRUE ~ as.character(rep)))%>%
     mutate(rep=case_when(rep=="B" ~ "2",
                          rep=="A" ~ "0",
                          TRUE ~ as.character(rep)))%>%
     filter(rep!=0)%>%select(marker,sample,rep,prob)%>%
     mutate(n=as.numeric(rep))%>%select(-rep)%>%
     mutate(prob=n*prob)%>%group_by(sample,marker)%>%
     summarise(sum=sum(prob))%>%ungroup()%>%
     left_join(.,depth_data,by=c("marker","sample"))%>%
     left_join(.,samples,by="sample")%>%
     select(ind,fullname,marker,sum,depth)%>%
     group_by(fullname)%>%
     mutate(count=length(fullname[depth>0]))%>%
     filter(count>=29)%>%
     select(marker,sum,fullname)%>%
     spread(marker,sum)%>%
     mutate(chr10_4753179=as.numeric("0.666"))

saveRDS(rf_data_extra,"rf_data_extra")

############## PREDICTIONS #################################################################################### #####
modelFit<-readRDS("modelFit")
fitdata<-readRDS("rf_data_extra")
residuals<-readRDS("residual")
sites<-read_excel("../map/coords.xlsx",sheet="working")%>%clean_names()

predictions<-as.data.frame(predict(modelFit,newdata=fitdata))%>%bind_cols(fitdata%>%select(fullname))%>%rename(predicted=1,filename=2)%>%mutate(dataset="new")%>%mutate(predicted=5.448*predicted)
sd(residuals$residual)/sd(predictions$predicted)


############## VALIDATION ##################################################################################### #####
mote<-read_xlsx("raw_data.xlsx",sheet="Validation_Phenotypes")%>%clean_names()%>%select(-trial,-date)%>%group_by(filename,genotype,status)%>%
     summarise(mean=mean(y_value))%>%
     spread(status,mean)%>%clean_names()%>%
     mutate(decline=(pre_bleaching-post_bleaching)/pre_bleaching)%>%
     left_join(.,predictions,by="filename")

summary(lm(decline~predicted,data=mote))

anno1=expression(italic("R"^2~"= 0.274,"~"p=0.098"))

a<-ggplot(mote)+
     geom_point(aes(decline,predicted))+
     geom_smooth(aes(decline,predicted),method="lm",color="orange")+
     #geom_label(aes(decline,predicted,label=genotype))+
     theme_classic(base_size=8)+
     ylab("Bleaching Residual (Predicted)")+
     xlab("Validation Bleaching\n(relative decline fv/fm)")+
     annotate("text", x = 0.8, y = 0.12, label = "More Bleaching",hjust=1,size=2,fontface="italic")+
     annotate("text", x = 0.5, y = -0.2, label = "Less Bleaching",hjust=0,size=2,fontface="italic")+
     annotate("text",x=0.8,y=-0.2,label=anno1,size=2,hjust=1)


############## PREDICTION STATS ############################################################################### #####
joint<-bind_rows(readRDS("residual")%>%select(-sample)%>%rename(sample=1)%>%mutate(dataset="core"),predictions%>%select(filename,predicted,dataset)%>%rename(sample=1,residual=2))%>%
     full_join(.,sites,by="sample")%>%
     mutate(region=case_when(dataset=="core"~"Experimental", TRUE ~as.character(region)))%>%
     mutate(keep=case_when(dataset=="core"~1,TRUE~as.numeric(keep)))%>%filter(keep==1)
write.table(joint,"~/Desktop/list.txt")
wilcox.test(residual~dataset,data=joint)

library(car);library(emmeans);library(lsmeans);library(MASS);library(multcomp)

### MODEL1: score by region
data<-joint%>%filter(dataset=="new")
model1<-glm((residual)~region,data=data) 
qqPlot(residuals(model1))
leveneTest(residuals(model1)~data$region) #homogeneity of variance
Anova(model1,test.statistic=c("F"))
emm<-emmeans(model1,~region)
cld(emm)
detach(package:multcomp);detach(package:car);detach(package:lsmeans);detach(package:TH.data);detach(package:MASS);detach(package:emmeans)

joint$region<-factor(joint$region,levels=c("Broward","Miami-Dade","Upper Keys","Middle Keys","Lower Keys","Dry Tortugas","Experimental"))
b<-ggplot(joint%>%filter(region!="Experimental"))+
     geom_hline(yintercept=0,color="gray",linetype="dashed")+
     geom_boxplot(aes(region,residual),fill="gray",outlier.color=NA)+
     geom_jitter(aes(region,residual,color=residual),size=1,width=0.05)+
     scale_color_viridis_c(direction=-1,limits=c(-0.3,0.4),name="Bleaching\nResidual")+
     theme_classic(base_size=8)+
     ylab("Bleaching Residual")+
     theme(axis.title.x=element_blank(),
           legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
     scale_y_continuous(limits=c(-0.3,.4))

c<-ggplot(joint%>%filter(region=="Experimental"))+
     geom_hline(yintercept=0,color="gray",linetype="dashed")+
     geom_boxplot(aes(region,residual),fill="gray",outlier.color=NA)+
     geom_jitter(aes(region,residual,color=residual),size=1,width=0.05)+
     scale_color_viridis_c(direction=-1,name="Bleaching\nResidual",limits=c(-0.3,0.4))+
     theme_classic(base_size=8)+
     ylab("Bleaching Residual")+
     theme(axis.title.x=element_blank(),
           legend.key.size = unit(0.25,"cm"),
           axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
           axis.text.y=element_blank(),
           axis.ticks.y=element_blank(),
           axis.line.y=element_blank(),
           axis.title.y=element_blank(),
           legend.position="none")+
     scale_y_continuous(limits=c(-0.3,.4))
     



############## PREDICTIONS MAP ################################################################################ #####
data<-joint%>%filter(dataset=="new")
library(rgdal);library(ggsn);library(sf)
florida<-st_read("../map/cb_2018_us_state_500k.shp")

d<-ggplot()+
     geom_sf(data=florida)+
     coord_sf(xlim=c(-83,-80),ylim=c(24.3,26.2))+
     geom_jitter(aes(long,lat,color=residual),data=data,width=0.02,height=0.02,size=1)+
     scale_color_viridis_c(direction=-1,breaks=seq(-0.2,0.4,0.1),name="Bleaching\nResidual",limits = c(min(data$residual),max(data$residual)))+
     theme_minimal(base_size=8)+
     scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks=seq(24,27,0.5))+
     scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks=seq(-83.5,-79.8,0.5))+
     theme(legend.position=c(0.1,0.6),
           legend.key.width = unit(0.25,"cm"),
           legend.key.height = unit(.6,"cm"),
           panel.border = element_rect(fill = "NA",colour = "gray35",size=.9))+
     xlab("Longitude (W)")+
     ylab("Latitude (N)")+
     annotate("text",x=-80.8,y=25.8,label="Florida")+
     scalebar(x.min=-81,x.max=-80,y.min=24.4,y.max=24.6,transform=TRUE,model="WGS84",dist=25,dist_unit="km",
              height = .25, st.dist = 0.4,
              box.fill = c("black", "white"),
              box.color = "black", border.size = .1,st.size=2)



quartz(w=7.2,h=3)
plot_grid(a,b,c,d,nrow=1,align="h",axis="tb",rel_widths=c(2.2,1.2,0.3,4),labels=c("A","B","","C"),label_size=8)




############## PREP GOMWU  #################################################################################### #####
loci<-readRDS("raw_lm_output")%>%select(marker)%>%rename(POS=1)

size<-5000 #insert size
lengths<-read_tsv("~/CRD_GBS/genome/Amil_v2.02/lengths.txt",col_names = FALSE)%>%rename(CHR=1,length=2) #produced with awk code see GO_protocol.docx

list<-left_join(loci%>%select(POS)%>%separate(POS,into=c("CHR","POS"),sep="_")%>%arrange(CHR,POS)%>%
                     mutate(POS=as.numeric(POS)),lengths,by="CHR")%>%rename(chr=CHR)%>%
     mutate(start=POS-(size/2))%>%
     mutate(start=case_when(start<0 ~ 0,TRUE ~ as.numeric(start)))%>%
     mutate(end=POS+(size/2))%>%
     mutate(end=case_when(end>length ~ length-1,TRUE ~as.numeric(end)))%>%
     select(chr,start,end,POS)%>%mutate(y=chr)%>%unite(name,y,POS,sep="_")
write.table(list,"./annotation_working/go_list_5k.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
system(paste("/anaconda3/envs/vcf/bin/bedtools getfasta -bed ./annotation_working/go_list_5k.bed -fi ~/CRD_GBS/genome/Amil_v2.02/Amil.v2.02.chrs.fasta -nameOnly -fo ./annotation_working/go_list_5k.fasta"))
system(paste("split -l 2000 ./annotation_working/go_list_5k.fasta ./annotation_working/list-"))
       

############## GO MWU ######################################################################################### #####
library(doMC)
registerDoMC(cores=4)
heats<-readRDS("raw_lm_output")%>%select(marker,r2)%>%rename(col=marker)
blast<-read_tsv("./annotation_working/blast_output.txt",col_names=FALSE)%>%
     select(X1,X2,X11)%>%rename(POS=X1,gene=X2,e_value=X11)%>%
     separate(gene,into=c("trash","gene"),sep="\\|",extra='drop')%>%select(-trash)%>%
     group_by(POS)%>%arrange(e_value,.by_group = TRUE)%>%group_by(POS)%>%slice_min(e_value,n=1,with_ties=T)

GO<-read_tsv("./annotation_working/GO_output.txt",col_names=FALSE)%>%gather(gene,GO,-X1)%>%select(-gene)%>%rename(gene=X1)%>%filter(GO!="NA")

table<-inner_join(GO,blast,by="gene")%>%select(POS,GO,-gene,-e_value)%>%distinct()#%>%semi_join(.,positions,by="POS")
write.table(heats,"../GO_MWU/0_heats.txt",sep=",",quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(table,"../GO_MWU/GO_rawtable.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
system(paste(perlPath="perl","../GO_MWU/nrify_GOtable.pl ../GO_MWU/GO_rawtable.txt > ../GO_MWU/0_table.txt"))
system(paste("rm ../GO_MWU/GO_rawtable.txt"))
setwd("~/CRD_GBS/3FD_Acerv/GO_MWU/")

#heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="MF";source("gomwu.functions.R")
#gomwuStats(heat_file, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.1, smallest=6, clusterCutHeight=0.25,Alternative="g")

heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="BP";source("gomwu.functions.R")
gomwuStats(heat_file, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.2, smallest=5, clusterCutHeight=0.05,Alternative="g")

#heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="CC";source("gomwu.functions.R")
#gomwuStats(heat_file, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.1, smallest=6, clusterCutHeight=0.2,Alternative="g")

################## PLOT BP #################################################################################### ######
setwd("~/CRD_GBS/3FD_Acerv/GO_MWU/")

absValue=0.2
level1=0.1
level2=0.05
level3=0.01
txtsize=1.2
treeHeight=0.5
adjusted=TRUE
require(ape)
library(ggdendro)

heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="BP";source("gomwu.functions.R")
input=heat_file
in.mwu=paste("MWU",goDivision,input,sep="_")
in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

cutoff=-log(level1,10)
pv=read.table(in.mwu,header=T)
row.names(pv)=pv$term
in.raw=paste(goDivision,input,sep="_")
rsq=read.table(in.raw,sep="\t",header=T)
rsq$term=as.factor(rsq$term)

if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval+1e-15,10)
heat$direction=0
heat$direction[pv$delta.rank>0]=1
if (cutoff>0) { 
     goods=subset(heat,pval>=cutoff) 
} else {
     goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
     goods=heat[row.names(heat) %in% goods.names,]
}

if (is.null(colors) | length(colors)<4 ) {
     colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
     if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
          colors=c("black","black","grey50","grey50")
     }
}
goods.names=row.names(goods)
GO_terms<-as.data.frame(goods.names)
# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods=diss[goods.names,goods.names]

# how many genes out of what we started with we account for with our best categories?
good.len=c();good.genes=c()
for (g in goods.names) {
     sel=rsq[rsq$term==g,]	
     pass=abs(sel$value)>=absValue
     sel=sel[pass,]
     good.genes=append(good.genes,as.character(sel$seq))
     good.len=append(good.len,nrow(sel))
}
ngenes=length(unique(good.genes))

#hist(rsq$value)
totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")

# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
goods=goods[labs,]
labs=sub(" activity","",labs)
labels<-bind_cols(as.data.frame(labs),GO_terms)%>%unite(label,labs,goods.names,sep=" - ")%>%
     separate(label,into=c("label","term"),sep=" - ")%>%separate(term,into=c("term","garbo"),sep=32)%>%
     select(-garbo)%>%unite(label,label,term,sep=" - ")
bp_print<-as.vector(labels$label)
bp_p_value<-semi_join(pv,GO_terms%>%rename(term=1),by="term")%>%mutate(font_display=case_when((p.adj<0.1&p.adj>0.05)~"p<0.1",
                                                                                              (p.adj<0.05&p.adj>0.01)~"p<0.05",
                                                                                              (p.adj<0.01)~"p<0.01"))

bp_hcdata<-dendro_data(cl.goods)
bp_hcdata$segments<-bp_hcdata$segments%>%mutate(yend=-yend/3,y=-y/3)


################## PLOT MF #################################################################################### ######
heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="MF";source("gomwu.functions.R")

input=heat_file
in.mwu=paste("MWU",goDivision,input,sep="_")
in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

cutoff=-log(level1,10)
pv=read.table(in.mwu,header=T)
row.names(pv)=pv$term
in.raw=paste(goDivision,input,sep="_")
rsq=read.table(in.raw,sep="\t",header=T)
rsq$term=as.factor(rsq$term)

if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval+1e-15,10)
heat$direction=0
heat$direction[pv$delta.rank>0]=1
if (cutoff>0) { 
     goods=subset(heat,pval>=cutoff) 
} else {
     goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
     goods=heat[row.names(heat) %in% goods.names,]
}

if (is.null(colors) | length(colors)<4 ) {
     colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
     if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
          colors=c("black","black","grey50","grey50")
     }
}
goods.names=row.names(goods)
GO_terms<-as.data.frame(goods.names)
# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods=diss[goods.names,goods.names]

# how many genes out of what we started with we account for with our best categories?
good.len=c();good.genes=c()
for (g in goods.names) {
     sel=rsq[rsq$term==g,]	
     pass=abs(sel$value)>=absValue
     sel=sel[pass,]
     good.genes=append(good.genes,as.character(sel$seq))
     good.len=append(good.len,nrow(sel))
}
ngenes=length(unique(good.genes))

#hist(rsq$value)
totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")

# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
goods=goods[labs,]
labs=sub(" activity","",labs)
labels<-bind_cols(as.data.frame(labs),GO_terms)%>%unite(label,labs,goods.names,sep=" - ")%>%
     separate(label,into=c("label","term"),sep=" - ")%>%separate(term,into=c("term","garbo"),sep=32)%>%
     select(-garbo)%>%unite(label,label,term,sep=" - ")
mf_print<-as.vector(labels$label)

mf_p_value<-semi_join(pv,GO_terms%>%rename(term=1),by="term")%>%mutate(font_display=case_when((p.adj<0.1&p.adj>0.05)~"p<0.1",
                                                                                              (p.adj<0.05&p.adj>0.01)~"p<0.05",
                                                                                              (p.adj<0.01)~"p<0.01"))%>%
     separate(term,into=c("term","garbo"),sep=32)%>%select(-garbo)


mf_hcdata<-dendro_data(cl.goods)
mf_hcdata$segments<-mf_hcdata$segments%>%mutate(yend=-yend/3,y=-y/3)


################## PLOT CC #################################################################################### ######
heat_file="0_heats.txt"; goAnnotations="0_table.txt";goDatabase="go.obo"; goDivision="CC";source("gomwu.functions.R")

input=heat_file
in.mwu=paste("MWU",goDivision,input,sep="_")
in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")

cutoff=-log(level1,10)
pv=read.table(in.mwu,header=T)
row.names(pv)=pv$term
in.raw=paste(goDivision,input,sep="_")
rsq=read.table(in.raw,sep="\t",header=T)
rsq$term=as.factor(rsq$term)

if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval+1e-15,10)
heat$direction=0
heat$direction[pv$delta.rank>0]=1
if (cutoff>0) { 
     goods=subset(heat,pval>=cutoff) 
} else {
     goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
     goods=heat[row.names(heat) %in% goods.names,]
}

if (is.null(colors) | length(colors)<4 ) {
     colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
     if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
          colors=c("black","black","grey50","grey50")
     }
}
goods.names=row.names(goods)
GO_terms<-as.data.frame(goods.names)
# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods=diss[goods.names,goods.names]

# how many genes out of what we started with we account for with our best categories?
good.len=c();good.genes=c()
for (g in goods.names) {
     sel=rsq[rsq$term==g,]	
     pass=abs(sel$value)>=absValue
     sel=sel[pass,]
     good.genes=append(good.genes,as.character(sel$seq))
     good.len=append(good.len,nrow(sel))
}
ngenes=length(unique(good.genes))

#hist(rsq$value)
totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")

# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
goods=goods[labs,]
labs=sub(" activity","",labs)
labels<-bind_cols(as.data.frame(labs),GO_terms)%>%unite(label,labs,goods.names,sep=" - ")%>%
     separate(label,into=c("label","term"),sep=" - ")%>%separate(term,into=c("term","garbo"),sep=32)%>%
     select(-garbo)%>%unite(label,label,term,sep=" - ")
cc_print<-as.vector(labels$label)

cc_p_value<-semi_join(pv,GO_terms%>%rename(term=1),by="term")%>%mutate(font_display=case_when((p.adj<0.1&p.adj>0.05)~"p<0.1",
                                                                                              (p.adj<0.05&p.adj>0.01)~"p<0.05",
                                                                                              (p.adj<0.01)~"p<0.01"))

cc_hcdata<-dendro_data(cl.goods)
cc_hcdata$segments<-cc_hcdata$segments%>%mutate(yend=-yend/3,y=-y/3)



################## FINAL GO FIGURE ############################################################################ #####

fontsize <- c("p<0.1"=2.5,"p<0.05"=2.5,"p<0.01"=3)
fontcolor <- c("p<0.1"="gray","p<0.05"="black","p<0.01"="red")
bp<-ggplot() + 
     geom_segment(data=segment(bp_hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
     geom_text(data=label(bp_hcdata), aes(x=x, y=y+0.1, label=bp_print, hjust=0,color=bp_p_value$font_display,size=bp_p_value$font_display),angle=0)+
     scale_color_manual(values=fontcolor)+
     scale_size_manual(values=fontsize)+
     scale_y_continuous(limits=c(min(bp_hcdata$segments$y)-0.1,5))+
     scale_x_continuous(limits=c(0.25,nrow(bp_p_value)+0.25))+
     theme_classic(base_size=8)+
     theme(axis.ticks.x=element_blank(),
           axis.text.x=element_blank(),
           axis.line.x=element_blank(),
           axis.title.x=element_blank(),
           axis.ticks.y=element_blank(),
           axis.line.y=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
           legend.position="none",
           legend.title=element_blank())+coord_flip()+
     xlab("BP");bp

mf<-ggplot() + 
     geom_segment(data=segment(mf_hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
     geom_text(data=label(mf_hcdata), aes(x=x, y=y+0.1, label=mf_print, hjust=0,color=mf_p_value$font_display,size=mf_p_value$font_display),angle=0)+
     scale_color_manual(values=fontcolor)+
     scale_size_manual(values=fontsize)+
     scale_y_continuous(limits=c(min(mf_hcdata$segments$y)-0.1,5))+
     scale_x_continuous(limits=c(0.25,nrow(mf_p_value)+0.25))+
     theme_classic(base_size=8)+
     theme(axis.ticks.x=element_blank(),
           axis.text.x=element_blank(),
           axis.line.x=element_blank(),
           axis.title.x=element_blank(),
           axis.ticks.y=element_blank(),
           axis.line.y=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
           legend.position="none",
           legend.title=element_blank())+
     xlab("MF")+
     coord_flip();mf

cc<-ggplot() + 
     geom_segment(data=segment(cc_hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
     geom_text(data=label(cc_hcdata), aes(x=x, y=y+0.1, label=cc_print, hjust=0,color=cc_p_value$font_display,size=cc_p_value$font_display),angle=0)+
     scale_color_manual(values=fontcolor)+
     scale_size_manual(values=fontsize)+
     scale_y_continuous(limits=c(min(cc_hcdata$segments$y)-0.1,5))+
     scale_x_continuous(limits=c(0.25,nrow(cc_p_value)+0.25))+
     theme_classic(base_size=8)+
     theme(axis.ticks.x=element_blank(),
           axis.text.x=element_blank(),
           axis.line.x=element_blank(),
           axis.title.x=element_blank(),
           axis.ticks.y=element_blank(),
           axis.line.y=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
           legend.position="none",
           plot.background = element_rect(fill = "transparent",colour = NA))+
     annotate("text",x=2.5,y=4.8,label="p<0.1",color="gray",size=2,hjust=0)+
     annotate("text",x=1,y=4.8,label="p<0.05",color="black",size=2,hjust=0)+
     coord_flip()+
     xlab("CC");cc

quartz(h=3,w=6)
plot_grid(bp,mf,ncol=1,rel_heights=c(19,4),labels=c("A","B","C"),label_size=8,align="v",axis="lr",label_y=c(1,1.1,1.1))




