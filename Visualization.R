library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(ggplot2)
library(gridExtra)
library(grid)


#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="") # QRJM(WBC, ANC, PLT) QRJM(WBC-ANC, ANC, PLT)
data_org=readRDS(paste(path,"/Transformed_DATA.rds",sep=""))


T_set=matrix(c(0.25,0.25,0.25,
               0.25,0.75,0.75,
               0.50,0.50,0.50,
               0.75,0.25,0.25,
               0.75,0.75,0.75),ncol=3, byrow=T)

overlap=function(x1,x2,y1,y2)
{
  intersects=ifelse((max(x1,x2)<=min(y1,y2))|(max(y1,y2)<=min(x1,x2)),"No","Yes")
  significance_x=ifelse(x1*x2<0,"No","Yes")
  significance_y=ifelse(y1*y2<0,"No","Yes")
  z=NULL
  if(intersects=="No")
  {
    z=0
  }
  else{
    
    if(significance_x==significance_y)
    {
      
      z=(min(max(x1,x2),max(y1,y2))-max(min(x1,x2),min(y1,y2)))/(max(x1,x2,y1,y2)-min(x1,x2,y1,y2))
    }
    else
    {
      z=0
    }
  }
  
  return(z)
}

overlap1=function(x1,x2,y1,y2)
{
  z=NULL
  if((max(x1,x2)<=min(y1,y2))|(max(y1,y2)<=min(x1,x2)))
  {
    z=0
  }
  else{
    
    
    z=(min(max(x1,x2),max(y1,y2))-max(min(x1,x2),min(y1,y2)))/(max(x1,x2,y1,y2)-min(x1,x2,y1,y2))
  }
  
  return(z)
}

Ov=function(x,y)
{
  qr1=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  qr2=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[y,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[y,]*100, collapse="_"),".rds",sep=""))
  params=names(qr1$mean)
  params=params[!params%in%c("u","deviance","Process","y",
                             "tau.gam.h0","tau_beta","tau_theta")]
  V=NULL
  for(i in 1:length(params))
  {
    prm=params[i]
    l1=qr1$q2.5[[prm]]
    u1=qr1$q97.5[[prm]]
    l2=qr2$q2.5[[prm]]
    u2=qr2$q97.5[[prm]]
    z=NA
    if(prm=="D")
    {
      z1=as.vector(sapply(1:nrow(l1),function(r) sapply(1:ncol(l1), function(c)as.numeric(ifelse(r<c,NA,overlap(l1[r,c],u1[r,c],l2[r,c],u2[r,c]))) )) )
      z1=z1[!is.na(z1)]
      z=c(median(z1),mean(z1),length(z1),prm)
      
    }
    if(prm!="D")
    {
      if(length(l1)<=1)
      {
        z1=overlap(l1,u1,l2,u2)
        z= c(median(z1),mean(z1),1,prm)
      }
      if(length(l1)>1)
      {
        if(is.na(dim(l1)[2]))
        {
          z1=sapply(1:length(l1),function(r) overlap(l1[r],u1[r],l2[r],u2[r]))
          z=c(median(z1),mean(z1), length(l1),prm)
        }
        if(!is.na(dim(l1)[2]))
        {
          z1=as.vector(sapply(1:nrow(l1),function(r) sapply(1:ncol(l1), function(c) overlap(l1[r,c],u1[r,c],l2[r,c],u2[r,c]))))
          z=c(median(z1),mean(z1),prod(dim(l1)),prm)
        }
      }
    }
    V=rbind(V,z)
    i=i+1
  }
  V=data.frame(V)
  colnames(V)=c('Median','Mean','length','Parameters')
  V$Median=as.numeric(V$Median)
  V$Mean=as.numeric(V$Mean)
  V$length=as.numeric(V$length)
  V$Overlapp=V$Median
  #Median Overlapp of parameters
  out=as.numeric(median(V$Overlapp))
  #out=as.numeric(mean(V$Mean))
  #out=sum(V$Mean*V$length)/sum(V$length)
  return(out)
  #return(V)
}




gamma=matrix(rep(NA,(nrow(T_set)^2)), ncol=nrow(T_set))
colnames(gamma)=sapply(1:nrow(T_set), function(x)paste0(T_set[x,]*100, collapse="_") )
rownames(gamma)=colnames(gamma)


gamma[1:nrow(T_set),1:nrow(T_set)]=sapply(1:nrow(T_set), function(r)sapply(1:nrow(T_set), function(c)Ov(r,c)))
M=nrow(T_set)


saveRDS(gamma,paste(path_out,"/QRJM_comp.rds", sep=""))



library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(ggplot2)
library(gridExtra)
library(grid)
gamma=readRDS(paste(path_out,"/QRJM_comp.rds", sep=""))


library(lattice)


M=ncol(gamma)
plot(x=c(.5,M+.5), y=c(.5,M+.5), type='n', xlab='', ylab='', 
     asp=1, xaxt='n', yaxt='n',	main=paste(), cex.main=1.0,
     cex.lab=1.0)

for(x in 1:M){ for(y in 1:M){
  points(x,y, pch=22, cex=3.9*sqrt(gamma[x,y]), bg='grey')}}
for(x in 1:M){
  points(x,x, pch=22, cex=3.9, bg='grey')}
axis(1 , at=c(1:5), labels=c(colnames(gamma)) , las=2)
axis(2 , at=c(1:5), labels=c(colnames(gamma)), las=2)



library(reshape2)
library(ggplot2)
library("RColorBrewer")
#display.brewer.pal(n = 9, name = 'YlGnBu')
cols=brewer.pal(n = 9, name = 'YlGnBu')[c(2,5,8)]


G=melt(gamma)
G$Var1=as.factor(G$Var1)
G$Var2=as.factor(G$Var2)
G$value=round(as.numeric(G$value)*100,1)

ggplot(data = G, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = cols[1], mid=cols[2],high = cols[3], 
                       midpoint=(100-0)/2,limit = c(0,100), space = "Lab", 
                       name="Overlap") +
  theme_minimal()+ coord_fixed()+
  theme(axis.text.x=element_text(size=10 , face="bold", angle = 90),
        axis.text.y=element_text(size=10 , face="bold"))+
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 4, fontface = "bold") +ylab("")+xlab("")



#Parameter graphs across taus
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(ggplot2)
library(gridExtra)
library(grid)

#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="") # QRJM(WBC, ANC, PLT) QRJM(WBC-ANC, ANC, PLT)
data_org=readRDS(paste(path,"/Transformed_DATA.rds",sep=""))


T_set=matrix(c(0.25,0.25,0.25,
               0.25,0.75,0.75,
               0.50,0.50,0.50,
               0.75,0.25,0.25,
               0.75,0.75,0.75),ncol=3, byrow=T)



#psi

prm=c("psi")
Par_WBC_ANC=NULL
Par_ANC=NULL
Par_PLT=NULL
for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  z_w_a=c(qr$q2.5[[prm]][1],qr$mean[[prm]][1],qr$q97.5[[prm]][1],paste0(T_set[x,]*100, collapse="_"))
  z_a=c(qr$q2.5[[prm]][2],qr$mean[[prm]][2],qr$q97.5[[prm]][2],paste0(T_set[x,]*100, collapse="_"))
  z_p=c(qr$q2.5[[prm]][3],qr$mean[[prm]][3],qr$q97.5[[prm]][3],paste0(T_set[x,]*100, collapse="_"))
  
  Par_WBC_ANC=rbind(Par_WBC_ANC,z_w_a)
  Par_ANC=rbind(Par_ANC,z_a)
  Par_PLT=rbind(Par_PLT,z_p)
  x=x+1
  
  
}
Par_WBC_ANC=data.frame(Par_WBC_ANC)
Par_ANC=data.frame(Par_ANC)
Par_PLT=data.frame(Par_PLT)
colnames(Par_WBC_ANC)=c("lci","mean","uci","tau")
colnames(Par_ANC)=c("lci","mean","uci","tau")
colnames(Par_PLT)=c("lci","mean","uci","tau")

Par_WBC_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_WBC_ANC[,x]))
Par_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_ANC[,x]))
Par_PLT[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_PLT[,x]))

rownames(Par_WBC_ANC)=1:nrow(T_set)
rownames(Par_ANC)=1:nrow(T_set)
rownames(Par_PLT)=1:nrow(T_set)




Par_WBC_ANC$Response="LYMPH"
Par_ANC$Response="ANC"
Par_PLT$Response="PLT"
Par=data.frame(rbind(Par_WBC_ANC,Par_ANC,Par_PLT))



Par$Tau_num=ifelse(Par$tau=="25_25_25",1,
                   ifelse(Par$tau=="25_75_75",2,
                          ifelse(Par$tau=="50_50_50",3,
                                 ifelse(Par$tau=="75_25_25",4,5))))
p<- ggplot(Par, aes(x=Tau_num, y=mean ,color=Response )) + scale_color_manual(values=c("red","blue","green"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Association Estimates")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)
p




#Trace plots of psi
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="") # QRJM(WBC, ANC, PLT) QRJM(WBC-ANC, ANC, PLT)
data_org=readRDS(paste(path,"/Transformed_DATA.rds",sep=""))


T_set=matrix(c(0.25,0.25,0.25,
               0.25,0.75,0.75,
               0.50,0.50,0.50,
               0.75,0.25,0.25,
               0.75,0.75,0.75),ncol=3, byrow=T)


summary_qrjm=function(x)
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  v=NULL
  for(i in 1:3)
  {
    d=data.frame(cbind(rep(i,dim(qr$sims.list$psi)[1]), qr$sims.list$psi[,i],1:dim(qr$sims.list$psi)[1]))
    v=rbind(d,v)
    i=i+1
  }
  
  v=data.frame(v)
  rownames(v)=NULL
  colnames(v)=c("Response", "psi","iterations")
  v$Response=ifelse(v$Response==1,"LYMPH", ifelse(v$Response==2,"ANC","PLT"))
  
  p1<-ggplot(data = v, aes(x=iterations, y=psi, group=Response))+
    geom_line(size=0.5,aes(color=Response))+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
    geom_hline(yintercept = qr$mean$psi, color=c("blue","red","green"), size=1.1)+
    xlab("Iterations") + ylab("Association Trace")
  p1<-p1+theme(legend.position="none")
  # leg<- get_legend(p1)
  
  p2<-ggplot(data = v, aes(x=psi, color=Response))+
    geom_density()+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
    geom_vline(xintercept = qr$mean$psi, color=c("blue","red","green"), size=1.1)+
    xlab("Posterior Samples") + ylab("Density")
  p2<-p2+theme(legend.position="none")
  p3<-ggarrange(p1,p2, nrow=1, common.legend = F)
  annotate_figure(p3, top = text_grob(paste("QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""), 
                                      color = "black", face = "bold", size = 10))
  
}



plt<-vector(mode = "list", length = nrow(T_set))
for(i in 1:nrow(T_set))
{
  plt[[i]]<-summary_qrjm(i)
  i=i+1
  
}
ind_leg=function(x)
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  v=NULL
  for(i in 1:3)
  {
    d=data.frame(cbind(rep(i,dim(qr$sims.list$psi)[1]), qr$sims.list$psi[,i],1:dim(qr$sims.list$psi)[1]))
    v=rbind(d,v)
    i=i+1
  }
  
  v=data.frame(v)
  rownames(v)=NULL
  colnames(v)=c("Response", "psi","iterations")
  v$Response=ifelse(v$Response==1,"LYMPH", ifelse(v$Response==2,"ANC","PLT"))
  
  p1<-ggplot(data = v, aes(x=iterations, y=psi, group=Response))+
    geom_line(size=0.5,aes(color=Response))+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
    geom_hline(yintercept = qr$mean$psi, color=c("blue","red","green"), size=1.1)+
    xlab("Iterations") + ylab("Association Trace")
  leg<-get_legend(p1)
  return(leg)
}

leg<-ind_leg(1)

ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],plt[[5]],leg, nrow=3, ncol=2, common.legend = F)







#beta medicine


prm=c("b")
Par_WBC_ANC=NULL
Par_ANC=NULL
Par_PLT=NULL


medicine=c("6MP","MTX")

for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  z_w_a_mp=c(qr$q2.5[[prm]][1,4],qr$mean[[prm]][1,4],qr$q97.5[[prm]][1,4],paste0(T_set[x,]*100, collapse="_"),"6MP")
  z_a_mp=c(qr$q2.5[[prm]][2,4],qr$mean[[prm]][2,4],qr$q97.5[[prm]][2,4],paste0(T_set[x,]*100, collapse="_"),"6MP")
  z_p_mp=c(qr$q2.5[[prm]][3,4],qr$mean[[prm]][3,4],qr$q97.5[[prm]][3,4],paste0(T_set[x,]*100, collapse="_"),"6MP")
  
  z_w_a_mtx=c(qr$q2.5[[prm]][1,5],qr$mean[[prm]][1,5],qr$q97.5[[prm]][1,5],paste0(T_set[x,]*100, collapse="_"),"MTX")
  z_a_mtx=c(qr$q2.5[[prm]][2,5],qr$mean[[prm]][2,5],qr$q97.5[[prm]][2,5],paste0(T_set[x,]*100, collapse="_"),"MTX")
  z_p_mtx=c(qr$q2.5[[prm]][3,5],qr$mean[[prm]][3,5],qr$q97.5[[prm]][3,5],paste0(T_set[x,]*100, collapse="_"),"MTX")
  
  
  Par_WBC_ANC=rbind(Par_WBC_ANC,z_w_a_mp,z_w_a_mtx)
  Par_ANC=rbind(Par_ANC,z_a_mp,z_a_mtx)
  Par_PLT=rbind(Par_PLT,z_p_mp,z_p_mtx)
  x=x+1
  
  
}

Par_WBC_ANC=data.frame(Par_WBC_ANC)
Par_ANC=data.frame(Par_ANC)
Par_PLT=data.frame(Par_PLT)

colnames(Par_WBC_ANC)=c("lci","mean","uci","tau","med")
colnames(Par_ANC)=c("lci","mean","uci","tau","med")
colnames(Par_PLT)=c("lci","mean","uci","tau","med")

Par_WBC_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_WBC_ANC[,x]))
Par_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_ANC[,x]))
Par_PLT[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_PLT[,x]))

rownames(Par_WBC_ANC)=NULL
rownames(Par_ANC)=NULL
rownames(Par_PLT)=NULL




Par_WBC_ANC$Response="LYMPH"
Par_ANC$Response="ANC"
Par_PLT$Response="PLT"
Par=data.frame(rbind(Par_WBC_ANC,Par_ANC,Par_PLT))

Par$Tau_num=ifelse(Par$tau=="25_25_25",1,
                   ifelse(Par$tau=="25_75_75",2,
                          ifelse(Par$tau=="50_50_50",3,
                                 ifelse(Par$tau=="75_25_25",4,5))))

p1<- ggplot(Par[Par$med=="6MP",], aes(x=Tau_num, y=mean,  color=Response)) + scale_color_manual(values=c("red","blue","green"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("6MP")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)

#p1

p2<- ggplot(Par[Par$med=="MTX",], aes(x=Tau_num, y=mean,  color=Response)) + scale_color_manual(values=c("red","blue","green"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("MTX")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)
#p2


ggarrange(p1,p2, ncol=1, common.legend = TRUE, legend="right")



p3<- ggplot(Par[Par$Response=="ANC",], aes(x=Tau_num, y=mean, color=med)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=-22.5))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("ANC Effect")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)+
  guides(color=guide_legend(title="Medicine"))

#p3

p4<- ggplot(Par[Par$Response=="PLT",], aes(x=Tau_num, y=mean, color=med)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=-22.5))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("PLT Effect")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)+
  guides(color=guide_legend(title="Medicine"))

#p4

p5<- ggplot(Par[Par$Response=="LYMPH",], aes(x=Tau_num, y=mean, color=med)) + scale_color_manual(values=c("red","blue"))+
  geom_line(size=1.2) +
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.5,size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=-22.5))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("LYMPH Effect")+xlab("tau (LYMPH, ANC, PLT)")+
  scale_x_continuous(labels=c("25_25_25","25_75_75","50_50_50","75_25_25","75_75_75"),
                     breaks=c(1,2,3,4,5))+
  geom_hline(yintercept = 0)+
  guides(color=guide_legend(title="Medicine"))

#p5
#ggarrange(p5,p3, p4, ncol=3, common.legend = TRUE, legend="bottom")
legend<-get_legend(p3)
p3<-p3+theme(legend.position="none")
p4<-p4+theme(legend.position="none")
p5<-p5+theme(legend.position="none")

grid.arrange(p5,p3,p4,legend, ncol=2)





#Fixed variable significance
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(ggplot2)
library("RColorBrewer")
library(ggpubr)


#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="") # QRJM(WBC, ANC, PLT) QRJM(WBC-ANC, ANC, PLT)
data_org=readRDS(paste(path,"/Transformed_DATA.rds",sep=""))


T_set=matrix(c(0.25,0.25,0.25,
               0.25,0.75,0.75,
               0.50,0.50,0.50,
               0.75,0.25,0.25,
               0.75,0.75,0.75),ncol=3, byrow=T)


si=function(x)
{
  if(min(x)<0 & max(x)>0)
  {
    z="No"
  }
  else{z="Yes"}
  return(z)
}
si(c(1,2))

sig=function(x)
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  
  l=qr$q2.5
  u=qr$q97.5
  
  z1=sapply(1:3, function(i) sapply(6:19, function(j) si(c(l$b[i,j],u$b[i,j])) ) )
  colnames(z1)=c("WBC-ANC","ANC","PLT")
  z2=sapply(1:14 , function(i) si(c(l$theta[i],u$theta[i])))
  z1=as.matrix(cbind(z1,z2))
  colnames(z1)=c("WBC-ANC","ANC","PLT","Surv")
  z1=data.frame(z1)
  z1$tau=rep(paste0(T_set[x,]*100, collapse="_"),14)
  z1$var=c("Age","Gender","Lineage","Wbc_p","NCI","Bulky","CNS","Risk","Risk_p","Day8","M_CR1","MRD","Day35","Status")
  return(z1)
  
}

df=NULL
for(x in 1:5)
{
  df=rbind(df,sig(x))
  x=x+1
}
df=data.frame(df)
rownames(df)=NULL

#Remove the Status and Risk covariates
df=subset(df, !df$var%in%c("Risk","Status"))
p1<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+ #coord_fixed()+
  geom_tile(aes(fill = ANC), colour = "white")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("ANC")+scale_fill_discrete(name = "Significant")

p2<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+# coord_fixed()+
  geom_tile(aes(fill = PLT), colour = "white")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("PLT")+scale_fill_discrete(name = "Significant")

p3<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+# coord_fixed()+
  geom_tile(aes(fill = WBC.ANC), colour = "white")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("LYMPH")+scale_fill_discrete(name = "Significant")

p4<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+# coord_fixed()+
  geom_tile(aes(fill = Surv), colour = "white")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("Survival")+scale_fill_discrete(name = "Significant")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


ggarrange(p3, p1,p2,p4,ncol=2,nrow=2,common.legend = TRUE, legend="bottom")



#Significance Polarity


si=function(x)
{
  if(min(x)<0 & max(x)>0)
  {
    z="Not"
  }
  if(max(x)<=0)
  {
    z="Negatively"
  }
  if(min(x)>=0)
  {
    z="Positively"
  }
  
  return(z)
}
si(c(1,1))

sig=function(x)
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  
  l=qr$q2.5
  u=qr$q97.5
  
  z1=sapply(1:3, function(i) sapply(6:19, function(j) si(c(l$b[i,j],u$b[i,j])) ) )
  colnames(z1)=c("WBC-ANC","ANC","PLT")
  z2=sapply(1:14 , function(i) si(c(l$theta[i],u$theta[i])))
  z1=as.matrix(cbind(z1,z2))
  colnames(z1)=c("WBC-ANC","ANC","PLT","Surv")
  z1=data.frame(z1)
  z1$tau=rep(paste0(T_set[x,]*100, collapse="_"),14)
  z1$var=c("Age","Gender","Lineage","Wbc_p","NCI","Bulky","CNS","Risk","Risk_p","Day8","M_CR1","MRD","Day35","Status")
  return(z1)
  
}

df=NULL
for(x in 1:5)
{
  df=rbind(df,sig(x))
  x=x+1
}
df=data.frame(df)
rownames(df)=NULL

#Remove the Status and Risk covariates
colrs=c("#6666FF","#FFFF33","#FF3333")
df=subset(df, !df$var%in%c("Risk","Status"))
p1<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=ANC),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("ANC")

p2<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=PLT),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("PLT")

p3<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=WBC.ANC),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("LYMPH")

p4<-ggplot(df, aes(x=tau, y=var)) + theme_minimal()+
  geom_tile(aes(fill=Surv),color="white", size=0.75)+
  scale_fill_manual(drop=FALSE, values=colrs,name = "Significant")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=10 ,color="black", face="bold", angle = -22.5),
        axis.text.y=element_text(size=10 ,color="black", face="bold"))+ylab("")+xlab("")+ggtitle("Relapse-Time")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


ggarrange(p3, p1,p2,p4,ncol=2,nrow=2,common.legend = TRUE, legend="bottom")



#Rhat
z1=NULL
for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  names_Rhat=names(qr$Rhat)
  z=quantile(sapply(names_Rhat,function(x) quantile(qr$Rhat[[x]], probs=c(0.75))), probs=c(0.75), names=F)
  z0=c( paste0(T_set[x,]*100, collapse="_"),round(z,3))
  z1=rbind(z1,z0)
  x=x+1
  
  
}
z1=data.frame(z1)
colnames(z1)=c("Model","Rhat")
z1$Rhat=as.numeric(z1$Rhat)
rownames(z1)=NULL
z1


#D


prm=c("D")
Par_ANC=NULL
Par_PLT=NULL
Par_Corr=NULL
Par_Cov=NULL
for ( x in 1:nrow(T_set))
{
  qr=readRDS(paste(paste(path_out,"/QRJM_",paste0(T_set[x,]*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(T_set[x,]*100, collapse="_"),".rds",sep=""))
  
  z_a=c(qr$q2.5[[prm]][1,1],qr$mean[[prm]][1,1],qr$q97.5[[prm]][1,1],paste0(T_set[x,]*100, collapse="_"))
  z_p=c(qr$q2.5[[prm]][2,2],qr$mean[[prm]][2,2],qr$q97.5[[prm]][2,2],paste0(T_set[x,]*100, collapse="_"))
  z_cov=c(qr$q2.5[[prm]][1,2],qr$mean[[prm]][1,2],qr$q97.5[[prm]][1,2],paste0(T_set[x,]*100, collapse="_"))
  
  z_c=c(qr$q50[[prm]][1,2]/sqrt(qr$q50[[prm]][1,1]*qr$q50[[prm]][2,2]))
  
  Par_ANC=rbind(Par_ANC,z_a)
  Par_PLT=rbind(Par_PLT,z_p)
  Par_Cov=rbind(Par_Cov,z_cov)
  
  Par_Corr=c(Par_Corr,z_c)
  x=x+1
  
  
}

Par_ANC=data.frame(Par_ANC)
Par_PLT=data.frame(Par_PLT)
Par_Cov=data.frame(Par_Cov)
colnames(Par_ANC)=c("lci","mean","uci","tau")
colnames(Par_PLT)=c("lci","mean","uci","tau")
colnames(Par_Cov)=c("lci","mean","uci","tau")

Par_ANC[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_ANC[,x]))
Par_PLT[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_PLT[,x]))
Par_Cov[,c("lci","mean","uci")]=sapply(c("lci","mean","uci"), function(x) as.numeric(Par_Cov[,x]))

rownames(Par_ANC)=1:nrow(T_set)
rownames(Par_PLT)=1:nrow(T_set)
rownames(Par_Cov)=1:nrow(T_set)


p11<-ggplot(Par_ANC, aes(x=1:nrow(T_set), y=mean)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1) +
  geom_line() +
  theme_bw()+
  geom_hline(yintercept = 1)+
  geom_point(size=2)+ scale_x_discrete(limits=as.character(c(1:nrow(Par_ANC))),
                                       labels=Par_ANC[,"tau"])+xlab('')+ylab('D_ANC')+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))



p12<-ggplot(Par_PLT, aes(x=1:nrow(T_set), y=mean)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1) +
  geom_line() +
  theme_bw()+
  geom_hline(yintercept = 1)+
  geom_point(size=2)+ scale_x_discrete(limits=as.character(c(1:nrow(Par_PLT))),
                                       labels=Par_PLT[,"tau"])+xlab('')+ylab('D_PLT')+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))


grid.arrange(p11,p12,nrow=2,heights=c(2.5,2.5))

p13<-ggplot(Par_Cov, aes(x=1:nrow(T_set), y=mean)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.1) +
  geom_line() +
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_point(size=2)+ scale_x_discrete(limits=as.character(c(1:nrow(Par_PLT))),
                                       labels=Par_PLT[,"tau"])+xlab('')+ylab('Cov(ANC,PLT)')+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))

Par_PLT$corr=Par_Corr
p14<-ggplot(Par_PLT,aes(x=1:nrow(T_set), y=Par_Corr))+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_point(size=2, color="blue")+
  scale_x_discrete(limits=as.character(c(1:length(Par_Corr))),
                   labels=Par_PLT[,"tau"])+xlab('')+ylab('Corr(ANC,PLT)')+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9, angle=90))
grid.arrange(p13,p14,nrow=2,heights=c(2.5,2.5))


grid.arrange(p11,p12,p13,p14,nrow=4,heights=c(2.5,2.5,2.5,2.5))

