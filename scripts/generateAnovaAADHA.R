# libraries ----

library(deSolve)
library(tidyr)
library(nlme)
library(ggplot2)
library(nlraa)
library(rmutil)
library(ggpubr)
library(rstatix)

# differential equations for reaction mechanism

# dAA  = -k1*AA + k2*DHAA
# dDHA =  k1*AA - k2*DHAA - k3 * DHA
# dDKG =                  + k3 * DHA

# function to solve by general solver

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

# define parameters

Time<-seq(0,100,1)
C0<-c(aa=2,dhaa=1,dkg=0)
pars<-c(k1=.2,k2=.05,k3=.1)

# solve ordinary differential equations with fixed parameters

y<-ode(y=C0,times=Time,func=aadhaa,parms=pars)
plot(y)

#  System of linear differential equations function: ----
# It defines the matrix with the coefficients of 3 reactions an solves it by matrix exponentiation.

aadha.lde<-function(Time,k1,k2,k3,aa0,dhaa0,dkg0=0){
  M <- matrix(c(-pars["k1"],pars["k2"],0, pars["k1"],-(pars["k2"]+pars["k3"]),0,0,pars["k3"],0), 3, 3, byrow=TRUE)
  print(M)
  yt <- lin.diff.eqn(A=M,initial=C0,t=Time)
  yt<-as.data.frame(cbind(Time,yt))
  names(yt)<-c("Time","aa","dhaa","dkg")
  return(yt)
}

yt<-aadha.lde(Time,k1=.2,k2=.05,k3=.1,aa0=2,dhaa=1,dkg0=0)
plot(y,obs=yt)

# Experimental data ----

presiones<-read.csv("data/presiones1.csv",sep=";",dec=",")

# reading and sortening experimental data, it can be optimized.

VitC<-subset(presiones,select=c("X","rep","processing","sweetener","tiempo","Ac.Ascorbico","Ac.Dehidroascorbico"))
AA<-VitC
AA$Concentracion<-AA$Ac.Ascorbico
AA$Species<-"Ascorbic"

DHA<-VitC
DHA$Concentracion<-DHA$Ac.Dehidroascorbico
DHA$Species<-"Dehydroascorbic"

VitC.l<-rbind(AA,DHA)
VitC.l$Ac.Ascorbico<-NULL
VitC.l$Ac.Dehidroascorbico<-NULL
VitC.l$Species<-factor(VitC.l$Species)

VitC.l$Temp<-factor(substr(VitC.l$X,0,3))
levels(VitC.l$Temp)<-c(4,20,4,20,4,20)
VitC.l$Temp<-as.numeric(as.character(VitC.l$Temp))

# create factors for covariate modelling

VitC.l$sweetener<-factor(VitC.l$sweetener)
VitC.l$processing<-factor(VitC.l$processing)
VitC.l$rep<-factor(VitC.l$rep)
VitC.l$Species<-factor(VitC.l$Species,ordered=F)

aadha.ldef<-function(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0=0){
  yhat<-rep(NA,length(Time))
  # define the components of the matrix as a data frame
  iPred<-data.frame(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0)
  
  for (i in 1:length(Time)){
    # build the coefficients matrix from iPred for every row. 
    # Is Time selected for loop arbitrary?
    M <- matrix(c(-iPred$k1[i],iPred$k2[i],0,iPred$k1[i],-(iPred$k2[i]+iPred$k3[i]),0,0,iPred$k3[i],0), 3, 3, byrow=TRUE)
    # define initial concentrations of compounds and runs solver function
    C0<- c(iPred$aa0[i],iPred$dhaa0[i],iPred$dkg0[i])
    yt <- lin.diff.eqn(A=M,initial=C0,t=Time[i])
    # add the value to the species evaluated for every step
    yhat[i]<-switch(as.character(iPred$Species[i]),
                    Ascorbic=yt[1],
                    Dehydroascorbic=yt[2])
  }
  return(yhat)
}

# add the result of prediction to the experimental data
VitC.l$yhat<-with(VitC.l,aadha.ldef(tiempo,Species,.02,.05,.1,12,2.5,0))

# model the parameters by non lineal generalized least squares

VitC.fm0<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                        k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=15,dhaa0=2),
               data=VitC.l,
               control=nls.control(maxiter=200),
               params=lk1+Ea+lk2+lk3~1,
               start=c(lk1=-2.5,Ea=32,lk2=-3,lk3=-1.5),
               weights=varIdent(~1|Species),
               verbose=T)


summary(VitC.fm0,cor=T)

# adding dependance of both lk3 and ldhaa0 on sweetener 

VitC.fm4<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                        k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=exp(laa0),dhaa0=exp(ldhaa0)),
               data=VitC.l,
               weights=varIdent(~1|Species),
               control=nls.control(maxiter=200),
               params=list(lk1~1,Ea~1,lk2~1,lk3~sweetener,
                           laa0~processing+sweetener,
                           ldhaa0~sweetener),
               start=c(lk1=-2.5,
                       Ea=32,
                       lk2=-3,
                       lk3=-1.5,0.001,
                       laa=log(15),rep(0.001,3),
                       ldhaa0=log(2),rep(0.001,1)),
               verbose=T)
summary(VitC.fm4)


VitC.fm8<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                                        k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=exp(laa0),dhaa0=exp(ldhaa0)),
               data=VitC.l,
               weights=varIdent(~1|Species),
               control=nls.control(maxiter=200),
               params=list(lk1~1,Ea~1,lk2~processing+sweetener,lk3~processing+sweetener,
                           laa0~1,
                           ldhaa0~1),
               start=c(lk1=-2.5,
                       Ea=32,
                       lk2=-3,rep(-0.001,4),
                       lk3=-1.5,rep(-0.001,2),
                       laa=15,
                       ldhaa0=2),
               verbose=T)

summary(VitC.fm4)
summary(VitC.fm8)

anova(VitC.fm0,VitC.fm4,VitC.fm8)
VitC.fmF<-VitC.fm8

# simulation ----

VitC.SIM <- expand.grid(tiempo=seq(0,90,length=50),
                        Species=levels(factor(VitC.l$Species)),
                        Temp=c(4,20),
                        processing=levels(factor(VitC.l$processing)),
                        sweetener=levels(factor(VitC.l$sweetener))
)


yhatSIM<-predict_gnls(VitC.fmF,newdata=VitC.SIM,interval="confidence")

VitC.SIM$concentration<-yhatSIM[,"Estimate"]
VitC.SIM$Q2.5<-yhatSIM[,"Q2.5"]
VitC.SIM$Q97.5<-yhatSIM[,"Q97.5"]

# anova ----

VitC.SIM.w <- VitC.SIM %>% pivot_wider(names_from = "Species", 
                                     values_from =c("concentration","Q2.5", "Q97.5"), names_vary = "slowest")
VitC.SIM.w$id <- as.character(with(data = VitC.SIM.w, factor(Temp):processing:sweetener))

VitC.w <- subset(presiones,select=c("X","rep","processing","sweetener","tiempo","Ac.Ascorbico","Ac.Dehidroascorbico"))

VitC.w$id <- as.character(with(data = VitC.w, factor(rep):factor(processing):factor(sweetener)))

anova_test(concentration_Delphinidin.3.O.sambubioside.5.O.glucoside~tiempo*processing*sweetener, data=VitC.SIM.w)

# rawExperimentalData

write.anova.table <- function(datos, compuesto){
  get_anova_table(anova_test(data = datos, dv=i, wid=id, 
                             between = c(processing, sweetener), within= c(tiempo)), correction = "auto" )
}

VitC.raw.Anova <- list()
for (i in c("Ac.Ascorbico", "Ac.Dehidroascorbico")){
  VitC.raw.Anova <- append(VitC.raw.Anova, c(i, write.anova.table(VitC.w,i)))
}

write.csv(x=VitC.raw.Anova, file="results/VitCRawAnova.csv")

# simulated data

VitC.SIM.Anova.paired <- list()
for (i in c("concentration_Ascorbic", "concentration_Dehydroascorbic")){
  
  VitC.SIM.Anova.paired <- append(VitC.SIM.Anova.paired, c(i,write.anova.table(VitC.SIM.w, i)))
  
}

write.csv(x=VitC.SIM.Anova.paired, file="results/VitCSimAnova_paired.csv")



formulas_anova <- lapply(list("concentration_Ascorbic", "concentration_Dehydroascorbic"),
                         function(x) formula(paste0(x,"~tiempo*sweetener*processing"),env=globalenv()))       


VitC.SIM.Anova.unpaired <- list()
for (i in seq(1:2)){
  
  VitC.SIM.Anova.unpaired <- append(VitC.SIM.Anova.unpaired, 
                                   c(as.character(formulas_anova[[i]][[2]]),
                                     get_anova_table(anova_test(formulas_anova[[i]], data=VitC.SIM.w), 
                                                     correction = "auto" )))
  
  
}


write.csv(x=VitC.SIM.Anova.unpaired, file="results/VitCSimAnova_unpaired.csv")


# plots ----

ggplot(VitC.SIM, aes(x=tiempo, y= concentration, col = processing)) +
  facet_grid(Species~factor(Temp)+sweetener, scales ="free") +
  geom_line()+
  geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = concentration, fill = processing), alpha= 0.1) +
  geom_point(data=VitC.l, aes(x=tiempo, y = Concentracion, fill =processing))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of Vitamin C by processing", 
          subtitle = "Dots are experimental points")


ggplot(VitC.SIM, aes(x=tiempo, y= concentration, col = sweetener)) +
  facet_grid(Species~factor(Temp)+processing, scales ="free") +
  geom_line()+
  geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = concentration, fill = sweetener), alpha= 0.1) +
  geom_point(data=VitC.l, aes(x=tiempo, y = Concentracion, fill =sweetener))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of Vitamin C by sweetener", 
          subtitle = "Dots are experimental points")







