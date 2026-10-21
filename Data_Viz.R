library(readxl)
library(stringr)
library(dplyr)
library(base)
library(stats)
require(gridExtra)
library(graphics)
library(ggpubr)
library(MASS)
library(grid)
library(RColorBrewer)


#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="") # QRJM(WBC, ANC, PLT) QRJM(WBC-ANC, ANC, PLT)
data_org=readRDS(paste(path,"/Transformed_DATA.rds",sep=""))



#Multinornal QQ plot
x<-data_org[,c('ANC','PLT')]

mqqnorm(x, main = "Multi-normal Q-Q Plot (Y_a,Y_b)")
d_sub=data_org%>%
  group_by(Child.Id)%>%
  summarise(Mean_ANC=mean(ANC),Mean_PLT=mean(PLT),
            Relapsed=ifelse((Relapsed[1])==0,'No','Yes'))



#Inital check
lnty=c("dotted","21","4452","41","solid")[c(1,5)]
p0<-NULL
colour_not=c("black","black")
p0<-ggplot(data = d_sub, aes(x = Mean_ANC, y = Mean_PLT)) +
  theme_bw()+theme(aspect.ratio = 1,panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())+
  scale_y_continuous(labels=c(expression(Q[b](0.1)),expression(Q[b](0.25)),expression(Q[b](0.5)),expression(Q[b](0.75)),expression(Q[b](0.9))),
                     breaks=c(quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  scale_x_continuous(labels=c(expression(Q[a](0.1)),expression(Q[a](0.25)),expression(Q[a](0.5)),expression(Q[a](0.75)),expression(Q[a](0.9))),
                     breaks=c(quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  #geom_point(aes(color=Relapsed))+ scale_color_manual(values=c("red","blue"))+
  geom_density2d(aes(linetype=Relapsed))+
  scale_linetype_manual(values=lnty)+
  labs(x='Y_a', y = 'Y_b')+#+xlim(c(6.5,8.25))+ylim(c(11,13.5))+
  guides(color= guide_legend("Density"),linetype= guide_legend("Relapsed", override.aes = list(color="black")), size=2)+
  theme(axis.title.x = element_text(size=12, face="bold", color=colour_not[1]),
        axis.title.y = element_text( size=12, face="bold",color=colour_not[2]),
        axis.text.x=element_text(size=11 , face= "bold", color=colour_not[1], angle=-22.5),
        axis.text.y=element_text(size=11 , face= "bold", color=colour_not[2]),
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"))
p0                 
                 
#Contour                 
                 

lnty=c("dotted","21","4452","41","dashed","solid")[c(2,6)]
colour_not=c("brown","#22771c")

vect=c(0.05,0.1,0.15,0.2,0.225,0.3)
p4<-NULL
colour_not=c("black","black")
p4<-ggplot(data = d_sub, aes(x = Mean_ANC, y = Mean_PLT)) +
  theme_bw()+theme(aspect.ratio = 1,panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())+
  scale_y_continuous(labels=c(expression(Q[b](0.1)),expression(Q[b](0.25)),expression(Q[b](0.5)),expression(Q[b](0.75)),expression(Q[b](0.9))),
                     breaks=c(quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  scale_x_continuous(labels=c(expression(Q[a](0.1)),expression(Q[a](0.25)),expression(Q[a](0.5)),expression(Q[a](0.75)),expression(Q[a](0.9))),
                     breaks=c(quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))))+
  #geom_point(aes(color=Relapsed))+ scale_color_manual(values=c("red","blue"))+
  geom_density2d(aes(linetype=Relapsed, color=as.factor(..level..)), breaks=vect, size=1.2, alpha=1 )+ 
  scale_color_manual(values=c("red","yellow2","blue","darkgrey","maroon1","purple"))+ 
  scale_linetype_manual(values=lnty)+
  labs(x='Y_a', y = 'Y_b')+#+xlim(c(6.5,8.25))+ylim(c(11,13.5))+
  guides(color= guide_legend("Density"),linetype= guide_legend("Relapsed", override.aes = list(color="black")), size=2)+
  theme(axis.title.x = element_text(size=12, face="bold", color=colour_not[1]),
        axis.title.y = element_text( size=12, face="bold",color=colour_not[2]),
        axis.text.x=element_text(size=11 , face= "bold", color=colour_not[1], angle=-22.5),
        axis.text.y=element_text(size=11 , face= "bold", color=colour_not[2]),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        #legend.key.size = unit(0.5, 'cm'),
        legend.text=element_text(size=11),legend.key.size = unit(0.5, 'cm'),
        legend.key.height= unit(0.8, 'cm'),
        legend.key.width= unit(1, 'cm'))

q_anc=quantile(d_sub$Mean_ANC, probs=c(0.1,0.25,0.5,0.75,0.9))
q_plt=quantile(d_sub$Mean_PLT, probs=c(0.1,0.25,0.5,0.75,0.9))

p41<-NULL
p41<-p4+geom_segment(aes(x = mean(q_anc[c(2,3)]), y = mean(q_plt[c(4,5)]), xend = q_anc[1]-0.15, yend = mean(q_plt[c(4,5)])),
                     arrow = arrow(length = unit(0.3, "cm")), color="brown", size=1.1)

p42<-p41+geom_segment(aes(x = mean(q_anc[c(4,5)])-0.5, y = mean(q_plt[c(2,3)])+0.3, xend = q_anc[5]+0.7, yend = mean(q_plt[c(2,3)])+0.3),
                      arrow = arrow(length = unit(0.3, "cm")), color="green", size=1.1)

p43<-p42+geom_segment(aes(x = mean(q_anc[c(3,4)]), y = mean(q_plt[c(2,3)]), xend = mean(q_anc[c(3,4)])+0.5, yend = mean(q_plt[c(2,3)])-1),
                      arrow = arrow(length = unit(0.3, "cm")), color="orange", size=1.1)


p44<-p43+geom_segment(aes(x = mean(q_anc[c(3,4)]), y = mean(q_plt[c(2,3)])+0.1, xend = mean(q_anc[c(3,4)]), yend = mean(q_plt[c(2,3)])+1.2),
                      arrow = arrow(length = unit(0.3, "cm")), color="black", size=1.1)
p44


