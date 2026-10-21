#The essential parameter
tau=c(0.10,0.10,0.10)



#Library
library(rjags)
library(R2jags)
library(jagsUI)
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)
library(lqmm)


#-------------------------------------------------------------------------------------------#
#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/QRJM/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/QRJM(WBC-ANC, ANC, PLT)",sep="")
#dir.create(path_out)


#Complete Data
data_org=readRDS(paste(path,"/data_0.rds",sep=""))
data_org$id_new=as.numeric(data_org$id_new)
data_org$time_id=as.numeric(data_org$time_id)
data_org=data_org[
  with(data_org, order(id_new, time_id)),
]


#Making every id to the same starting point.
min_w=data_org%>%
  group_by(Child.Id)%>%
  summarise(min_week=min(Week)*277)
data_org=merge(x=data_org, y=min_w,by=c("Child.Id"),all.x=T)
#data_org$Week=data_org$Week*277-data_org$min_week+1
#data_org$Last_follow_up_date=data_org$Last_follow_up_date*277-data_org$min_week+1
#data_org$Time_to_event=data_org$Time_to_event*277-data_org$min_week+1

#Scaling the time (max week is now 271 weeks instead of previous 277)
max_week=max(max(data_org$Week),max(data_org$Last_follow_up_date),
             max(data_org$Time_to_event))
data_org$Week=data_org$Week/max_week
data_org$Time_to_event=data_org$Time_to_event/max_week
data_org$Last_follow_up_date=data_org$Last_follow_up_date/max_week

#Removing columns
data_org=data_org[,!colnames(data_org)%in%c("Height","Weight","BSA","min_week")]



#Removing id's with < 5 times points and updating the id names
id_set=unique(data_org$Child.Id)
dat=NULL
for(i in 1:length(id_set))
{
  d=subset(data_org, data_org$Child.Id==id_set[i])
  dat=rbind(dat, c(id_set[i],nrow(d)))
  i=i+1
}
dat=data.frame(dat)
rownames(dat)=1:nrow(dat)
dat$X2=as.numeric(dat$X2)
dat1=subset(dat,dat$X2<5)
data_org1=subset(data_org, !data_org$Child.Id%in%unique(dat1$X1))
unique(data_org1$id_new)

dat1=data_org1[,c("Child.Id","id_new")]
dat1=dat1[!duplicated(dat1), ]
dat1$id_new=sapply(dat1$Child.Id,function(x)as.numeric(strsplit(x,"_")[[1]][2]))
dat1=dat1[with(dat1, order(id_new)), ]
dat1$id_new=1:nrow(dat1)
data_org1=merge(x=data_org1[,!colnames(data_org1)%in%c("id_new")], y=dat1, by=c("Child.Id"), all.x=T)
data_org1=data_org1[,c(1,26,2:25)]
data_org=data_org1
#----------------------------------------------------------------------------------------------#
#                                 Covariate transformation
#----------------------------------------------------------------------------------------------#

#Medicine (transformed --> untransformed )
#-------------------------------------------#
data_org$X6_mp=data_org$X6_mp*275
data_org$MTX=data_org$MTX*12.5

#Normalizing medicine
#---------------------#
data_org$X6_mp= (data_org$X6_mp)/median(data_org$X6_mp)
data_org$MTX= (data_org$MTX)/median(data_org$MTX)


#Normalizing fixed covariates
#-----------------------------#
data_org$WBC_high=log(exp(data_org$WBC_high)*21250)
data_org$WBC_high=data_org$WBC_high-median(unique(data_org$WBC_high))


data_org$Age=data_org$Age*4.280008 + 6.589453
data_org$Age=(data_org$Age-mean(unique(data_org$Age)))/sd(unique(data_org$Age))

#----------------------------------------------------------------------------------------------#
#                                 Response transformation
#----------------------------------------------------------------------------------------------#


#Center response y to 0 (untransform --> center to 0)
#----------------------------------------------------#
data_org$WBC=data_org$WBC +log(3200)
data_org$ANC=data_org$ANC+log(1564)

data_org$WBC_new=exp(data_org$WBC)-exp(data_org$ANC)
data_org$WBC_new=log(data_org$WBC_new)
data_org$WBC=0
data_org$WBC=data_org$WBC_new

data_org=data_org[,!colnames(data_org)%in%c("WBC_new")]
data_org$Platelet=data_org$Platelet+log(192000)

#Centering to 0
#---------------#
data_org$WBC=data_org$WBC-median(data_org$WBC[!is.na(data_org$WBC)])
data_org$ANC=data_org$ANC-median(data_org$ANC[!is.na(data_org$ANC)])
data_org$Platelet=data_org$Platelet-median(data_org$Platelet[!is.na(data_org$Platelet)])

saveRDS(data_org,paste(path,"/Transformed_DATA.rds",sep=""))

#----------------------------------------------------------------------------------------------#


#File path
path_out1=paste(path_out,"/QRJM_",paste0(tau*100, collapse="_"), sep="")
dir.create(path_out1)
setwd(path_out1)

#number of ids 
#-----------------------------------------------------------------------------------------------------
n <- length(unique(data_org$id_new))

#Response
#-----------------------------------------------------------------------------------------------------
y <- as.matrix(data_org[,c("id_new","time_id","WBC","ANC","Platelet")])


#Data with only unduplicated values of fixed variables for UPN's
dataL=data_org
#This is a fake variable (to be discarded later), ie, ANC_0
dataL$ANC_0=ifelse(!is.na(dataL$ANC),dataL$ANC,0)
id <- as.integer(transform(dataL , id_new=as.numeric(factor(id_new)))$id_new)
#offset
#-----------------------------------------------------------------------------------------------------
offset <- as.vector(c(1, 1+cumsum(tapply(id,id,length))))
timeVar <- "Week"
times <-as.vector(dataL[[timeVar]])


dat1=dataL[,!colnames(dataL)%in%c("Week","WBC","ANC","Platelet","MTX","X6_mp","ANC_0")]
dataL.id=dat1[!duplicated(dat1$Child.Id),]

survObject <- coxph(Surv( time =Time_to_event ,as.numeric(Relapsed)) ~id_new+Age+Gender+Lineage+WBC_high+NCI+Bulky+CNS+Risk+Risk_present+Day8+Morph_CR1+MRD+Day35+Status , data = dataL.id , x=TRUE , model=TRUE)
#W
#----------------------------------------------------------------------------------------------------------
W <-survObject$x #Baseline Covariates
W=as.matrix(W)
rownames(W)=1:nrow(W)

SurvInf <- survObject$y
typeSurvInf <-attr(SurvInf, "type")

#for right censored
if(typeSurvInf=="right")
{
  Time <- SurvInf[,"time"]
  Time[Time < 1e-04] <-1e-04
  nT <- length(Time)
  event <-SurvInf[,"status"]
  LongFormat <- FALSE
}

#zeros
#-----------------------------------------------------------------------------------------------------
zeros <- rep(0,nrow(dataL.id))


#c
#-----------------------------------------------------------------------------------------------------
C <- 10^7



#Information regarding K control points used in Gauss-Konrod
#quadrature rule to approx the integral  in the hazard function ,
#which does not have a closed form solution 
#-----------------------------------------------------------------------------------------------------

#G-K points
Gk_points=as.matrix(cbind(c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                            0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                            -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                            0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329),
                          c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                            0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                            0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                            0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)))

colnames(Gk_points)=c("cuts","weights")
Gk_points=data.frame(Gk_points)
Gk_points= Gk_points[with(Gk_points, order(cuts)),]
rownames(Gk_points)=1:nrow(Gk_points)
#G-K cuts
#sk
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
sk <- Gk_points$cuts
#G-K weights
#wk
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
wk <- Gk_points$weights
#K_pt , P
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
K_pt <- length(sk)   #Points in G-K quadrature
P <-0.5*(Time-dataL$Week[offset[c(-length(offset))]])
st <- outer(P , (sk+1))+ matrix(rep(dataL$Week[offset[c(-length(offset))]],each=length(sk)),byrow=T,ncol=length(sk))
id.GK <-rep(seq_len(nT) , each = K_pt)


#medicine covariate
medi=function(id, t)
{
  d=subset(data_org, data_org$id_new==id)
  d=d[with(d, order( time_id)),]
  z=NULL
  if(t <=max(d$Week) & t>=min(d$Week))
  {
    for(j in 1:(nrow(d)-1))
    {
      if(t>=d[j,c("Week")] & t<d[(j+1),c("Week")])
      {
        z = d[j,c("X6_mp","MTX")]
      }
    }
    if(t==max(d$Week))
    {
      z=d[nrow(d),c("X6_mp","MTX")]
    }
  }
  else{
    z=c(0,0)
  }
  return(as.numeric(z))
}



#X
#----------------------------------------------------------------------------------------------------------

formYx <-ANC_0~id_new+Week+I(Week^2)+X6_mp+MTX+Age+Gender+Lineage+WBC_high+NCI+Bulky+CNS+Risk+Risk_present+Day8+Morph_CR1+MRD+Day35+Status+Time_to_event+I(Time_to_event^2)
mfX <-model.frame(terms(formYx) , data=dataL)
TermsX <-attr(mfX , "terms")
X <-model.matrix(formYx, mfX) 
X <- as.matrix(X)
X <- X[,c(2,1,3:ncol(X))]
X1=data.frame(X)
X1=X1[with(X1, order( id_new,Week)),]
colnames(X1)=colnames(X)

X=as.matrix(X1[,!colnames(X1)%in%c("Time_to_event","I(Time_to_event^2)")]) #Design matrix of fixed effects

#XT
#-------------------------------------------------------------------------------------------------------
XT=X1[,!colnames(X1)%in%c("X6_mp","MTX","Week","I(Week^2)")]
XT=XT[!duplicated(X1$id_new),]
rownames(XT)=1:nrow(XT)
XT$X6_mp=NA
XT$MTX=NA
for( i in 1:nrow(XT))
{
  v=as.numeric(medi(as.numeric(XT[i,"id_new"]),as.numeric(XT[i,"Time_to_event"])))
  XT[i,c("X6_mp","MTX")]=v
  i=i+1
}
XT=XT[,colnames(XT)[c(1:2,(ncol(XT)+c(-3:0)), (2+1:14))]] #Design matrix for Part 1 survival.
colnames(XT)=c(colnames(XT)[1:2],c("Week","I(Week^2)"),colnames(XT)[5:ncol(XT)])

XT=as.matrix(XT)
rownames(XT)=1:nrow(XT)


#Xs
#-----------------------------------------------------------------------------------------------------
dataL.id2 <- dataL.id[id.GK,]
dataL.id2 <- data.frame(dataL.id2)
dataL.id2$ANC_0=0
dataL.id2$X6_mp=0
dataL.id2$MTX=0
dataL.id2[[timeVar]] <-pmax(c(t(st)),0)
mfX.id2 <-model.frame(TermsX , data=dataL.id2)
Xs <- model.matrix(formYx,mfX.id2)
Xs <-Xs[,!colnames(Xs)%in%c("Time_to_event","I(Time_to_event^2)")]
for(j in 1:nrow(Xs))
{
  vs=as.numeric(medi(as.numeric(Xs[j,"id_new"]),as.numeric(Xs[j,"Week"])))
  Xs[j,c("X6_mp","MTX")]=vs
  j=j+1
}
Xs <- Xs [,colnames(Xs)[c(2,1,3:ncol(Xs))]]
Xs <- as.matrix(Xs)
rownames(Xs)=1:nrow(Xs)
#W2T , W2s
#------------------------------------------------------------------------------------------------------
#Obtain Bspline basis for computation
tpower <- function(x,t,p) ((x-t)^p)*as.numeric(x>t)

bbase <- function(x,xl,xr,ndx,deg)
{
  dx <- (xr-xl)/ndx
  kn=tail(head(seq(xl, xr, length.out = (ndx + 1)), -1), -1)
  kn <- kn[kn < max(Time)]
  knots<- sort(c(rep(range(Time, st), (deg + 1)), kn)) #Extended knots
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n) , diff = deg +1)/gamma((deg + 1)*dx^deg)
  B <- (-1)^(deg+1) *P%*%t(D)
  B
}
#P splines of degree = deg with (ndx +1) equispaced knots on [xl,xr]
W2T <-bbase(x=Time ,xl=quantile(Time, probs =c(0.05),names=F) , xr=quantile(Time, probs =c(0.95),names=F) , ndx=14 , deg=3)

W2s <- bbase(x=c(t(st)) ,xl=quantile(Time, probs =c(0.05),names=F) , xr=quantile(Time, probs =c(0.95),names=F) , ndx=14 , deg=3)


DD <-diag(ncol(W2T))
M.gam.h0 <-crossprod(diff(DD,differences =2))+1e-06*DD




#Longitudinal indices according to Process grid points
ProL=NULL
ProS=NULL
ids=sort(unique(data_org$id_new))
for(i in 1:length(ids))
{
  d=subset(data_org, data_org$id_new==ids[i])
  tpts=sort(as.vector(c(min(d$Week),(d$Time_to_event[1]-min(d$Week))*(sk+1)/2 + min(d$Week),d$Time_to_event[1])))
  d=d[with(d, order(Week)),]
  
  d$index=(i-1)*(length(sk)+1)+sapply(1:nrow(d),  function(j) ifelse(d$Week[j]>=max(tpts),(length(tpts)-1),sum(sapply(2:length(tpts), function(i) ifelse(d$Week[j]>=tpts[i-1]&d$Week[j]<tpts[i], (i-1),0)))))
  tpts1=sqrt(c(min(tpts),diff(tpts[1:(length(tpts)-1)])))
  
  ProL=rbind(ProL,d[,c("id_new","index")])
  ProS =rbind(ProS,cbind(rep(d$id_new[1],length(tpts1)),c(1:length(tpts1))+(i-1)*(length(sk)+1),tpts1))
  colnames(ProS)=c("id_new","Process_time_id","Process_time")
  i=i+1
}
ProL=as.matrix(ProL)
ProS=as.matrix(ProS)

#Process offset
offset_Pr=(0:(n))*(length(wk)+1)+1


#-------------------------------------------------------------------------------------------#
#                                  QR Model for initialization
#-------------------------------------------------------------------------------------------#

sink(paste(path_out1,"QR.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 #BM Process
 #----------
 u[offset_Pr[i],1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[offset_Pr[i],1:ncy] = u[offset_Pr[i],1:ncy] * ProS[offset_Pr[i],3]
 
 for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
 {
 u[w,1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[w,1:ncy] = Process[(w-1),1:ncy] + u[w,1:ncy] * ProS[w,3]
 }
 
 
 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
 
 for(k in 1:ncy)
    {
    #Exponential variable
    Ex[j,k] ~ dexp(lambda[k])
    
    #mean
    mu[j,k] =  inprod(b[k,1:(ncX-1)], X[j,2:ncX]) +  Process[ProL[j,2],k]  + t1[k]*Ex[j,k]
    tau[j,k] = lambda[k]/(t2[k]*Ex[j,k])
    y[j,(2+k)] ~ dnorm(mu[j,k], tau[j,k])

    } #loop of k
  } #loop of j

 
}#loop of i

#Prior
#-------------------------------------
  #Exponential parameter
  #-----------------------------------
  for(k in 1:ncy)
  {
  lambda[k] ~ dgamma(0.001,0.001)
  sig[k] = 1/lambda[k]
  }
  
  
   
  #Cov matrix for BM Processes
  #-----------------------------------
 
  D_inv[1:ncy,1:ncy] ~ dwish(Omega[,], (ncy+1))
  D[1:ncy,1:ncy] <-inverse(D_inv[1:ncy,1:ncy])
  
  
  
  #Polynomial and medicine and fixed variables
  #-----------------------
  for(k in 1:ncy)
  {
  for (j in 1:(ncX-1)) 
  {
   b[k,j] ~ dnorm(0,0.001)
  }
  }
  
 

  
}
",fill = TRUE)
sink()

#-------------------------------------------------------------------------------------------#
# Initialize chain
#-------------------------------------------------------------------------------------------#

#Data list for QR
#------------------------#
dat_list_qr =list('n'=n,'y'=y ,
                  'offset'=offset,'offset_Pr'=offset_Pr ,'ProL'=ProL,'ProS'=ProS,
                  'X'=X  , 'ncX'=ncol(X) , 'ncy'=(ncol(y)-2) ,
                  'ze'=rep(0,(ncol(y)-2)),'Omega' = diag(rep(1,(ncol(y)-2))) ,
                  't1'=as.vector(sapply(1:length(tau), function(t) (1-2*tau[t])/(tau[t]*(1-tau[t])))), #ALD decomposition constant 1
                  't2'=as.vector(sapply(1:length(tau), function(t) 2/(tau[t]*(1-tau[t]))))) #ALD decomposition constant 2) 

#Initialization for QR (substitute responses by 0: not to be used later)
#------------------------#
data_subs_response=dataL
data_subs_response$WBC=ifelse(is.na(data_subs_response$WBC),0,data_subs_response$WBC)
data_subs_response$ANC=ifelse(is.na(data_subs_response$ANC),0,data_subs_response$ANC)
data_subs_response$Platelet=ifelse(is.na(data_subs_response$Platelet),0,data_subs_response$Platelet)

fit_y1<- suppressWarnings(lqmm(fixed = WBC ~ Week+I(Week^2)+X6_mp+MTX+Age+Gender+Lineage+WBC_high+NCI+Bulky+CNS+Risk+Risk_present+Day8+Morph_CR1+MRD+Day35+Status,
                               random = ~ sqrt(Week), group = id_new,	control =list(LP_max_iter = 2000, UP_max_iter = 2500),
                               data = data_subs_response, tau = tau[1], nK = 15, type = "normal"))


fit_y2<- suppressWarnings(lqmm(fixed = ANC ~ Week+I(Week^2)+X6_mp+MTX+Age+Gender+Lineage+WBC_high+NCI+Bulky+CNS+Risk+Risk_present+Day8+Morph_CR1+MRD+Day35+Status,
                               random = ~ sqrt(Week), group = id_new,	control =list(LP_max_iter = 2000, UP_max_iter = 2500),
                               data = data_subs_response , tau = tau[2], nK = 15, type = "normal"))


fit_y3<- suppressWarnings(lqmm(fixed = Platelet ~ Week+I(Week^2)+X6_mp+MTX+Age+Gender+Lineage+WBC_high+NCI+Bulky+CNS+Risk+Risk_present+Day8+Morph_CR1+MRD+Day35+Status,
                               random = ~ sqrt(Week), group = id_new,	control =list(LP_max_iter = 2000, UP_max_iter = 2500),
                               data = data_subs_response, tau = tau[3], nK = 15, type = "normal"))



initsList_qr = function(){list( D_inv = solve(diag(as.vector(c(VarCorr.lqmm(fit_y1)[2],VarCorr.lqmm(fit_y2)[2],VarCorr.lqmm(fit_y3)[2])))) ,
                                lambda=c(1/fit_y1$scale,1/fit_y2$scale,1/fit_y3$scale),
                                b = rbind(coef(fit_y1),coef(fit_y2),coef(fit_y3)))}

#Run QR
#-------#
jagsfit_qr <- jags(data=dat_list_qr, 
                   inits=initsList_qr, 
                   parameters.to.save=c('u','b','D','sig'),
                   model.file=paste(path_out1,'QR.model.txt', sep="/"),parallel = TRUE,
                   n.chains=2, n.iter=15000, n.burnin=10000,
                   n.thin=100,n.adapt = 1000,
                   DIC=TRUE)



#Save QR
saveRDS(jagsfit_qr,paste(path_out1,"/init_QR_fit_",paste0(tau*100, collapse="_"),".rds",sep=""))

#Read QR
qr=readRDS(paste(path_out1,"/init_QR_fit_",paste0(tau*100, collapse="_"),".rds",sep=""))

u=as.matrix(qr$mean$u)
Pr=matrix(rep(NA,(prod(dim(u)))), nrow=nrow(u))
for(i in 1:n) 
{
  #BM Process
  #----------
  
  Pr[offset_Pr[i],] = u[offset_Pr[i],] * ProS[offset_Pr[i],3]
  
  for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
  {
    
    Pr[w,] = Pr[(w-1),] + u[w,] * ProS[w,3]
    w=w+1
  }
  i=i+1
}

#Calculation longitudinal QR lines
long1=X[,c(-1)]%*% qr$mean$b[1,] +Pr[ProL[,c("index")],1]
long2=X[,c(-1)]%*% qr$mean$b[2,] +Pr[ProL[,c("index")],2]
long3=X[,c(-1)]%*% qr$mean$b[3,] +Pr[ProL[,c("index")],3]

#Fitting time-dependent Cox model
#-------------------------------------------------------------------------------------------#
DF=data.frame(id=id,Time=Time[id],event=event[id])
Wdat=as.data.frame(W)
DF=cbind(DF,Wdat[id,c(-1)])
DF <- cbind(DF, long1,long2,long3)
DF$start <- times
splitID <- split(DF[c("start","Time")] , DF$id)
DF$stop <- unlist(lapply(splitID , function(d) c(d$start[-1],d$Time[1])))
DF$event <- with( DF , ave(event, id , FUN = function(x) c(rep(0,length(x)-1),x[1])))
DF <- DF[!names(DF)%in%c("Time","id")]

fit_s =coxph(Surv(start, stop, event)~. , data = DF[DF$stop>DF$start,])


#--------------------------------------------------------------------------------------------------------#
#                         Final initialization for QRJM
#--------------------------------------------------------------------------------------------------------#

initsList_qrjm = function(){list( u = qr$mean$u ,
                                  D_inv = solve(qr$mean$D) ,
                                  gam.h0 = rnorm(n = dat_list_qrjm$ncW2T , mean = 0 , sd = 1),
                                  lambda=1/qr$mean$sig,
                                  b = qr$mean$b ,
                                  psi=unname(coef(fit_s)[c("long1","long2","long3")]) ,
                                  theta = unname(coef(fit_s)[!names(coef(fit_s))%in%c("long1","long2","long3")]))}

#Data list for QRJM
dat_list_qrjm =list('n'=n,'y'=y,'zeros'=zeros ,'C'=10^7,
                    'K'=length(wk) , 'P'=P , 'wk'=wk , 'offset'=offset,'offset_Pr'=offset_Pr ,'ProL'=ProL,'ProS'=ProS,
                    'X'=X , 'XT'=XT , 'Xs'=Xs , 'ncX'=ncol(X) ,'ncXT'=ncol(XT) ,'ncXs'=ncol(Xs) , 'ncy'=(ncol(y)-2) ,
                    'W'=W , 'W2T'=W2T , 'W2s'=W2s , 'ncW'=ncol(W) , 'ncW2T'=ncol(W2T) , 'ncW2s'=ncol(W2s) ,
                    'event'=event ,'ze'=rep(0,(ncol(y)-2)),'Omega' = diag(rep(1,(ncol(y)-2))) ,'priorMean.gam.h0'=seq(-8.5,0,len =ncol(W2s)) , 'M.gam.h0'=M.gam.h0 ,
                    'a.gam.h0'=1 , 'b.gam.h0'=0.05,
                    't1'=as.vector(sapply(1:length(tau), function(t) (1-2*tau[t])/(tau[t]*(1-tau[t])))), #ALD decomposition constant 1
                    't2'=as.vector(sapply(1:length(tau), function(t) 2/(tau[t]*(1-tau[t]))))) #ALD decomposition constant 2) 



#-------------------------------------------------------------------------------------------#
#                                       QRJM Model
#-------------------------------------------------------------------------------------------#
sink(paste(path_out1,"QRJM.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 #BM Process
 #----------
 u[offset_Pr[i],1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[offset_Pr[i],1:ncy] = u[offset_Pr[i],1:ncy] * ProS[offset_Pr[i],3]
 
 for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
 {
 u[w,1:ncy] ~ dmnorm(ze[],D_inv[,])
 Process[w,1:ncy] = Process[(w-1),1:ncy] + u[w,1:ncy] * ProS[w,3]
 }
 
 
 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
 
 for(k in 1:ncy)
    {
    #Exponential variable
    Ex[j,k] ~ dexp(lambda[k])
    
    #mean
    mu[j,k] =  inprod(b[k,1:(ncX-1)], X[j,2:ncX]) +  Process[ProL[j,2],k]  + t1[k]*Ex[j,k]
    tau[j,k] = lambda[k]/(t2[k]*Ex[j,k])
    y[j,(2+k)] ~ dnorm(mu[j,k], tau[j,k])

    } #loop of k
  } #loop of j

 
 
 #Survival part 1
 #---------------
 
 log.h0T[i] = inprod(gam.h0[1:ncW2T],W2T[i,1:ncW2T]) #Baseline log hazard
 theta.W[i] = inprod(theta[1:(ncW-1)],W[i,2:ncW]) #Baseline Covariate 
 for(k in 1:ncy)
 {
 muT[i,k] = inprod(b[k,1:(ncX-1)],XT[i,2:ncXT]) +  Process[(offset_Pr[i+1]-1),k]  #Fixed part of Longitudinal
 } 
 

 
 
 #Survival part 2
 #---------------
 
 for(l in 1:K)
  {
  log.h0s[i,l] = inprod(gam.h0[1:ncW2s],W2s[(K*(i-1)+l),1:ncW2s]) #Baseline log hazard
  for(k in 1:ncy)
 {
 mus[i,k,l] = inprod(b[k,1:(ncX-1)],Xs[(K*(i-1)+l),2:ncXs]) +  Process[((i-1)*(K+1) +1+l),k]  #Fixed part of Longitudinal
 } 
  SurvLong[i,l] = (wk[l]* exp(log.h0s[i,l]  + inprod(psi[1:ncy], mus[i,1:ncy,l])))

  }

  
  
  #Zeros trick in JAGS for JM
  #---------------------------------
  phi[i] = C- ((event[i]*(log.h0T[i] + theta.W[i] + inprod(psi[1:ncy],muT[i,1:ncy]))) -exp(theta.W[i]) * P[i]* sum(SurvLong[i,1:K]))
  zeros[i] ~ dpois(phi[i])
}#loop of i

#Prior
#-------------------------------------
  #Exponential parameter
  #-----------------------------------
  for(k in 1:ncy)
  {
  lambda[k] ~ dgamma(0.001,0.001)
  sig[k] = 1/lambda[k]
  }
  
  
   
  #Cov matrix for BM Processes
  #-----------------------------------
 
  D_inv[1:ncy,1:ncy] ~ dwish(Omega[,], (ncy+1))
  D[1:ncy,1:ncy] <-inverse(D_inv[1:ncy,1:ncy])
  
  #Beta
  #-----
  #Polynomial and medicine and fixed variables
  #-------------------------------------------
  for(k in 1:ncy)
  {
  for (j in 1:(ncX-1)) {
   b[k,j] ~ dnorm(0,0.001)
  }
  }
  
  

  #psi
  #-----
  for(k in 1:ncy)
  {
  psi[k] ~ dnorm(0,0.001)
  }

  #theta
  #-----
  for(h in 1:(ncW-1))
  {
  theta[h] ~  dnorm(0,0.001)
  }
  
  
  
  #Prior for penalized B-spline basis coefficients
  #------------------------------------------------
  gam.h0[1:ncW2T] ~  dmnorm(priorMean.gam.h0[],tau.gam.h0*M.gam.h0[,])
  
  #Smoothing parameter
  #------------------------------------------------
  tau.gam.h0 ~ dgamma(a.gam.h0,b.gam.h0)
  
  
}
",fill = TRUE)
sink()

#Run QRJM
jagsfit_qrjm <- jags(data=dat_list_qrjm, 
                     inits=initsList_qrjm, 
                     parameters.to.save=c('u','b','D','sig','psi','theta','gam.h0','tau.gam.h0','y'),
                     model.file=paste(path_out1,'QRJM.model.txt', sep="/"),parallel = TRUE,
                     n.chains=2, n.iter=20000, n.burnin=10000,
                     n.thin=50,n.adapt = 2000,
                     DIC=TRUE)


u=NULL
u=as.matrix(jagsfit_qrjm$mean$u)
Process=matrix(rep(NA,(prod(dim(u)))), nrow=nrow(u))
for(i in 1:n) 
{
  #BM Process
  #----------
  
  Process[offset_Pr[i],] = u[offset_Pr[i],] * ProS[offset_Pr[i],3]
  
  for( w in (offset_Pr[i]+1):(offset_Pr[i+1]-1))
  {
    
    Process[w,] = Process[(w-1),] + u[w,] * ProS[w,3]
    w=w+1
  }
  i=i+1
}
jagsfit_qrjm$mean$Process=Process

#Save QRJM

saveRDS(jagsfit_qrjm,paste(path_out1,"/QRJM_fit_",paste0(tau*100, collapse="_"),".rds",sep=""))


#pD
jagsfit_qrjm$pD

#DIC
jagsfit_qrjm$DIC
#-------------------------------------------------------------------------------------------#
#                                       Summarize QRJM
#-------------------------------------------------------------------------------------------#

Stats_pos=function(tau,statis)
{
  
  qrjm=readRDS(paste(paste(path_out,"/QRJM_",paste0(tau*100, collapse="_"), sep=""),"/QRJM_fit_",paste0(tau*100, collapse="_"),".rds",sep=""))
  model_names=paste0("QRJM(",paste0(tau*100, collapse=","),")",sep="")
  
  
  
  beta_names=c( "(Intercept)","Week","I(Week^2)",
                "X6_mp","MTX",
                "Age","Gender","Lineage","WBC_high","NCI",
                "Bulky","CNS","Risk","Risk_present","Day8",
                "Morph_CR1","MRD","Day35","Status")
  
  gamma_names=beta_names[-c(1:5)]
  D_names=c("WBC","ANC","PLT")
  psi_names=D_names
  
  betas=NULL
  Ds=NULL
  sigs=NULL
  psis=NULL
  gammas=NULL
  tau_beta=NULL
  tau_theta=NULL
  ci=function(a,b)
  {
    z=paste0("(",round(a,3),",",round(b,3),")",collapse="")
    return(z)
  }
  
  if(statis%in%c("mean","sd","q50","Rhat"))
  {
    stat=qrjm[[statis]]
    #beta
    betas=stat$b
    colnames(betas)=beta_names
    rownames(betas)=D_names
    
    #D
    Ds=stat$D
    colnames(Ds)=D_names
    rownames(Ds)=D_names
    
    #sig
    sigs=stat$sig
    names(sigs)=D_names
    
    #psi
    psis=stat$psi
    names(psis)=D_names
    
    #gamma
    gammas=stat$theta
    names(gammas)=gamma_names
    
    #tau_beta
    tau_beta=NA
    
    #tau_theta
    tau_theta=NA
    
    df=data.frame(round(rbind(t(betas),c(NA,NA,NA),
                              sigs,c(NA,NA,NA),
                              Ds,c(NA,NA,NA),
                              psis,c(NA,NA,NA),
                              cbind(gammas,rep(NA, length(gammas)),rep(NA, length(gammas))),c(NA,NA,NA),
                              c(tau_beta,NA,NA),c(NA,NA,NA),
                              c(tau_theta,NA,NA)),3)) 
    
  }
  if(statis%in%c("CI"))
  {
    stat1=qrjm[["q2.5"]]
    stat2=qrjm[["q97.5"]]
    
    #beta
    betas=sapply(1:ncol(stat1$b), function(c) sapply(1:nrow(stat1$b),function(r) ci(stat1$b[r,c],stat2$b[r,c])))
    colnames(betas)=beta_names
    rownames(betas)=D_names
    
    #D
    Ds=sapply(1:ncol(stat1$D), function(c) sapply(1:nrow(stat1$D),function(r) ci(stat1$D[r,c],stat2$D[r,c])))
    colnames(Ds)=D_names
    rownames(Ds)=D_names
    
    #sig
    sigs=sapply(1:length(stat1$sig), function(r) ci(stat1$sig[r],stat2$sig[r]))
    names(sigs)=D_names
    
    #psi
    psis=sapply(1:length(stat1$psi), function(r) ci(stat1$psi[r],stat2$psi[r]))
    names(psis)=D_names
    
    #gamma
    gammas=sapply(1:length(stat1$theta), function(r) ci(stat1$theta[r],stat2$theta[r]))
    names(gammas)=gamma_names
    
    #tau_beta
    tau_beta="(NA,NA)"
    
    #tau_theta
    tau_theta="(NA,NA)"
    
    df=data.frame(rbind(t(betas),c(NA,NA,NA),
                        sigs,c(NA,NA,NA),
                        Ds,c(NA,NA,NA),
                        psis,c(NA,NA,NA),
                        cbind(gammas,rep(NA, length(gammas)),rep(NA, length(gammas))),c(NA,NA,NA),
                        c(tau_beta,NA,NA),c(NA,NA,NA),
                        c(tau_theta,NA,NA))) 
  }
  
  
  df$points=c(beta_names,"--",
              "<-->","--",
              D_names,"--",
              "<-->","--",
              gamma_names,"--",
              "tau_beta","--",
              "tau_theta")
  df$Parameters=c(rep("beta",length(beta_names)),c("--"),
                  "sigma",c("--"),
                  rep("D",length(D_names)),c("--"),
                  "psi",c("--"),
                  rep("gamma",length(gamma_names)),c("--"),
                  c("tau_beta"),c("--"),
                  c("tau_theta"))
  df$statistics=statis
  
  df$Model=rep(model_names,nrow(df))
  rownames(df)=NULL
  df=df[,c("Model","statistics","Parameters","points","WBC","ANC","PLT")]
  return(df)
}

pos=data.frame(rbind(Stats_pos(tau,"mean"),
                     Stats_pos(tau,"sd"),
                     Stats_pos(tau,"CI")))
write.csv(pos,paste(path_out1,"/",pos$Model[1],".csv",sep=""))


#-------------------------------------------------------------------------------------------#
#Checks
#-------------------------------------------------------------------------------------------#
jagsfit_qrjm$q2.5$b[1,1:5]
jagsfit_qrjm$q97.5$b[1,1:5]
tau
jagsfit_qrjm$q2.5$b[2,]
jagsfit_qrjm$q97.5$b[2,]

jagsfit_qrjm$q2.5$psi
jagsfit_qrjm$q97.5$psi
paste(round(jagsfit_qrjm$q2.5$b[1,],3),round(jagsfit_qrjm$q97.5$b[1,],3),sep=" , ")
paste(round(jagsfit_qrjm$q2.5$b[2,],3),round(jagsfit_qrjm$q97.5$b[2,],3),sep=" , ")
#paste(round(jagsfit_qrjm$q2.5$b[3,],3),round(jagsfit_qrjm$q97.5$b[3,],3),sep=" , ")


v=sapply(1:ncol(jagsfit_qrjm$mean$b),
         function(j) ifelse(jagsfit_qrjm$q2.5$b[2,j]*jagsfit_qrjm$q97.5$b[2,j]<0,"Ins","Sig"))
table(v)
v=sapply(1:ncol(jagsfit_qrjm$mean$b),
         function(j) ifelse(jagsfit_qrjm$q2.5$b[1,j]*jagsfit_qrjm$q97.5$b[1,j]<0,"Ins","Sig"))
table(v)


v=sapply(1:length(jagsfit_qrjm$mean$theta),
         function(j) ifelse(jagsfit_qrjm$q2.5$theta[j]*jagsfit_qrjm$q97.5$theta[j]<0,"Ins","Sig"))
table(v)
