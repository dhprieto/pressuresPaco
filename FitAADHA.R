library(bbmle)
library(deSolve)
library(tidyr)
library(nlme)
library(ggplot2)

Time<-seq(0,100,1)
# dAA  = -k1*AA + k2*DHAA
# dDHA =  k1*AA - k2*DHAA - k3 * DHA
# dDKG =                  + k3 * DHA
aadhaa<-function(t,y,parms){
  aa<-y[1]
  dhaa<-y[2]
  k1<-parms["k1"]
  k2<-parms["k2"]
  k3<-parms["k3"]
  daa<- -k1*aa + k2*dhaa
  ddhaa<- k1*aa - k2*dhaa - k3*dhaa
  dkg <- k3*dhaa
  return(list(c(daa,ddhaa,dkg)))
}
C0<-c(aa=2,dhaa=1,dkg=0)
pars<-c(k1=.2,k2=.05,k3=.1)
y<-ode(y=C0,times=Time,func=aadhaa,parms=pars)
plot(y)
##  System of linear differential equations:
aadha.lde<-function(Time,k1,k2,k3,aa0,dhaa0,dkg0=0){
M <- matrix(c(-pars["k1"],pars["k2"],0, pars["k1"],-(pars["k2"]+pars["k3"]),0,0,pars["k3"],0), 3, 3, byrow=TRUE)
M
yt <-  lin.diff.eqn(A=M,initial=C0,t=Time)
yt<-as.data.frame(cbind(Time,yt))
names(yt)<-c("Time","aa","dhaa","dkg")
return(yt)
}
yt<-aadha.lde(Time,.2,.05,.1,2,1,0)
plot(y,obs=yt)

presiones<-read.csv("presiones1.csv",sep=";",dec=",")
VitC<-subset(presiones,select=c("X","rep","processing","sweetener","tiempo","Ac.Ascorbico","Ac.Dehidroascorbico"))
AA<-VitC
AA$Concentracion<-AA$Ac.Ascorbico
AA$Species<-"Ac.Ascorbico"
DHA<-VitC
DHA$Concentracion<-DHA$Ac.Dehidroascorbico
DHA$Species<-"Ac.Dehidroascorbico"
VitC.l<-rbind(AA,DHA)
VitC.l$Ac.Ascorbico<-NULL
VitC.l$Ac.Dehidroascorbico<-NULL
summary(VitC.l)
# create temperature
VitC.l$Temp<-factor(substr(VitC.l$X,0,3))
levels(VitC.l$Temp)<-c(4,20,4,20,4,20)
VitC.l$Temp<-as.numeric(as.character(VitC.l$Temp))
#create factors for covariate modelling
VitC.l$sweetener<-factor(VitC.l$sweetener)
VitC.l$processing<-factor(VitC.l$processing)
VitC.l$rep<-factor(VitC.l$rep)
VitC.l$Species<-factor(VitC.l$Species,ordered=F)

ggplot(data=VitC.l,aes(x=tiempo,y=Concentracion,col=Species))+geom_point()+facet_grid(processing~sweetener)

aadha.ldef<-function(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0=0){
  yhat<-rep(NA,length(Time))
  iPred<-data.frame(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0)
  for (i in 1:length(Time)){
    M <- matrix(c(-iPred$k1[i],iPred$k2[i],0,iPred$k1[i],-(iPred$k2[i]+iPred$k3[i]),0,0,iPred$k3[i],0), 3, 3, byrow=TRUE)
    C0<-c(iPred$aa0[i],iPred$dhaa0[i],iPred$dkg0[i])
    yt <-  lin.diff.eqn(A=M,initial=C0,t=Time[i])
    yhat[i]<-switch(as.character(iPred$Species[i]),
                    Ac.Ascorbico=yt[1],
                    Ac.Dehidroascorbico=yt[2])
  }
  return(yhat)
}

VitC.l$yhat<-with(VitC.l,aadha.ldef(tiempo,Species,.02,.05,.1,12,2.5,0))
ggplot(data=VitC.l,aes(x=tiempo,y=Concentracion,col=Species))+
  geom_point()+
  geom_line(aes(y=yhat))+
  facet_grid(processing~sweetener)


VitC.fm0<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                        k1=exp(lk1-Ea1/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea2/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea3/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=aa0,dhaa0=dhaa0),
               data=VitC.l,
               param=list(lk1~1,Ea1~1,lk2~1,Ea2~1,lk3~1,Ea3~1,aa0~1,dhaa0~1),
               start=c(lk1=0.02,Ea1=10,lk2=0.05,Ea2=10,lk3=0.1,Ea3=10,aa0=15,dhaa0=2.5),
               verbose=T)

library(minpack.lm)

VitC.fm0<-nlsLM(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                        k1=exp(lk1-Ea1/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea2/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea3/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=aa0,dhaa0=dhaa0),
               data=VitC.l,
               control=nls.lm.control(maxiter=200,maxfev=2000),
               start=c(lk1=0.02,Ea1=10,lk2=0.05,Ea2=10,lk3=0.1,Ea3=10,aa0=15,dhaa0=2.5))

VitC.fm00<-nlsLM(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                         k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         aa0=15,dhaa0=2),
                data=VitC.l,
                control=nls.lm.control(maxiter=200,maxfev=2000),
                start=c(lk1=0.02,Ea=10,lk2=0.05,lk3=0.1))


VitC.fm000<-nls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                          k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                          k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                          k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                          aa0=15,dhaa0=2),
                 data=VitC.l,
                 control=nls.control(maxiter=200),
                 start=c(lk1=-2.5,Ea=32,lk2=-3,lk3=-1.5))


VitC.fm000<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                         k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                         aa0=15,dhaa0=2),
                data=VitC.l,
                control=nls.control(maxiter=200),
                params=lk1+Ea+lk2+lk3~1,
                start=c(lk1=-2.5,Ea=32,lk2=-3,lk3=-1.5))


summary(VitC.fm0,cor=T)
summary(VitC.fm00,cor=T)
summary(VitC.fm000,cor=T)


VitC.pred<-expand.grid(tiempo=seq(0,90,length=50),
                       Species=levels(factor(VitC.l$Species)),
                       grouping=levels(factor(VitC.l$X))
)
VitC.pred$Temp<-factor(substr(VitC.pred$grouping,0,3))
levels(VitC.pred$Temp)<-c(4,20,4,20,4,20)
VitC.pred$Temp<-as.numeric(as.character(VitC.pred$Temp))
#create factors for covariate modelling
VitC.pred$sweetener<-factor(with(VitC.pred,substr(grouping,0,2)))
VitC.pred$processing<-factor(with(VitC.pred,substr(grouping,3,3)))
VitC.pred$sweetener<-factor(VitC.pred$sweetener)
VitC.pred$processing<-factor(VitC.pred$processing)
VitC.pred$Species<-factor(VitC.pred$Species,ordered=F)

VitC.pred$Concentration<-predict(VitC.fm000,newdata = VitC.pred)
yhat<-predict_gnls(VitC.fm000,newdata=VitC.pred,interval="prediction")
VitC.pred$Q2.5<-yhat[,"Q2.5"]
VitC.pred$Q97.5<-yhat[,"Q97.5"]

ggplot(VitC.l, 
       aes(x = tiempo, y = Concentracion, col =Species)) +
  facet_grid(processing~sweetener)+
  geom_line(data=VitC.pred,aes(x=tiempo,y=Concentration,col=Species))+
  geom_ribbon(data=VitC.pred,aes(ymax=Q97.5,ymin=Q2.5,y=Concentration,col=Species,fill=Species),alpha=0.2)+
  geom_point()+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")
ggsave(filename="Figure1AA.pdf")


