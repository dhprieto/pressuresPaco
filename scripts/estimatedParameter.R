library(texreg)
library(nlraa)
library(nlme)
library(dplyr)
r2<-function(obj){
  rss<-sum(resid(obj)^2)
  dep<-fitted(obj)+resid(obj)
  tss<-sum((mean(dep)-dep)^2)
  r2<-1-rss/tss
  n<-length(dep)
  p<-length(coef(obj))
  r2adj<-1-(n-1)/(n-p)*(1-r2)
  return(c(r2,r2adj))
}
recode_factor(data$TrialType, congruent = "Con", 
              incongruent = "InCon")

anthocyanins$compound <- recode_factor(anthocyanins$compound, Delphinidin.3.O.sambubioside.5.O.glucoside = "DpSG",
                                           Delphinidin.3.5.O.diglucoside="DpGG",
                                           Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside = "CySG+CyGG",	
                                           Delphinidin.3.O.sambubioside="DpS" , Delphinidin.3.O.glucoside="DpG",	
                                           Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside = "CyS+CyG")

list.of.compound.fit<-list()
list.of.r2adj<-list()
for (i in levels(factor(anthocyanins$compound))){
  print(i)
  list.of.compound.fit[[i]]<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
                                  data=subset(anthocyanins, compound == i),
                                  param=list(C0~compound+compound:sweetener+compound:processing,
                                             Cinf~compound+compound:sweetener+compound:processing,
                                             lk~compound+compound:sweetener+compound:processing,
                                             Ea~compound),
                                  start=c(C0=c(coef(ant.fm0)[1:6],rep(0.001,18)),
                                          Cinf=c(coef(ant.fm0)[7:12],rep(0.001,18)),
                                          lk=c(coef(ant.fm0)[13:18],rep(0.001,18)),
                                          Ea=c(coef(ant.fm0)[19:24],rep(0.001,0))),
  )
  
  list.of.r2adj[[i]]<-round(r2(list.of.compound.fit[[i]])[2],3)
}


screenreg(list.of.compound.fit,single.row=T,ci.force=T)

htmlreg(list.of.compound.fit,single.row=T,ci.force=T,
        groups=list("C0 [mg/100mL]"=1:24,"Cinf [mg/100mL]"=25:48,"ln(k) [1/min]"=49:72,"Ea [KJ/Mol/K]"=73:78),
        custom.gof.rows=list("R2adj"=list.of.r2adj),
        file="anthocyanins_params.html")
        #caption.above=T,
        #caption="Table X. Comparison of the apparent first order kinetic parameters and effects of processing and sweetener used",)


fit1 <- gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
             data=anthocyanins,
             param=list(C0~compound+compound:sweetener+compound:processing,
                        Cinf~compound+compound:sweetener+compound:processing,
                        lk~compound+compound:sweetener+compound:processing,
                        Ea~compound),
             start=c(C0=c(coef(ant.fm0)[1:6],rep(0.001,18)),
                     Cinf=c(coef(ant.fm0)[7:12],rep(0.001,18)),
                     lk=c(coef(ant.fm0)[13:18],rep(0.001,18)),
                     Ea=c(coef(ant.fm0)[19:24],rep(0.001,0))))

fit2 <- gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
             data=anthocyanins,
             param=list(C0~compound+sweetener+processing,
                        Cinf~compound+sweetener+processing,
                        lk~compound+sweetener+compound:processing,
                        Ea~compound),
             start=c(C0=c(coef(ant.fm0)[1:6],rep(0.001,3)),
                     Cinf=c(coef(ant.fm0)[7:12],rep(0.001,3)),
                     lk=c(coef(ant.fm0)[13:18],rep(0.001,13)),
                     Ea=c(coef(ant.fm0)[19],rep(0.001,5))
             )
)

screenreg(fit1,groups=list("C0 [mg/100mL]"=1:24,"Cinf [mg/100mL]"=25:48,"ln(k) [1/min]"=49:72,"Ea [KJ/Mol/K]"=73:78),
          single.row=T,ci.force=T)


wordreg(fit1,single.row=T,ci.force=T,
        groups=list("C0 [mg/100mL]"=1:24,"Cinf [mg/100mL]"=25:48,"ln(k) [1/min]"=49:72,"Ea [KJ/Mol/K]"=73:78),
        custom.gof.rows=list("R2adj"=round(r2(fit1)[2],3)),
        file="anthocyanins_paramsOnly1.doc")


wordreg(list(fit2, fit1),single.row=T,ci.force=T,
        groups=list("C0 [mg/100mL]"=1:9,"Cinf [mg/100mL]"=10:18,"ln(k) [1/min]"=19:37,"Ea [KJ/Mol/K]"=38:43,
                    "C0 [mg/100mL]" = 44:61, "Cinf [mg/100mL]"= 62:79, "ln(k) [1/min]"=80:85),
                custom.gof.rows=list("R2adj"=c(round(r2(fit2)[2],3),round(r2(fit1)[2],3))),
                file="anthocyanins_params2.doc")




# Vitamin C ----

library(deSolve)
library(rmutil)

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


list.of.compound.fit<-list()
list.of.r2adj<-list()
for (i in levels(factor(VitC.l$Species))){
  print(i)
  list.of.compound.fit[[i]]<- gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
                              k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                              k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                              k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                              aa0=exp(laa0),dhaa0=exp(ldhaa0)),
     data=subset(VitC.l, Species==i),
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
             ldhaa0=2))
  
  list.of.r2adj[[i]]<-round(r2(list.of.compound.fit[[i]])[2],3)
}

subset(VitC.l, Species=="Ascorbic" & rep == "rep1")

screenreg(list.of.compound.fit,single.row=T,ci.force=T)

htmlreg(list.of.compound.fit,single.row=T,ci.force=T,
        groups=list("C0 [mg/100mL]"=1:24,"Cinf [mg/100mL]"=25:48,"ln(k) [1/min]"=49:72,"Ea [KJ/Mol/K]"=73:78),
        custom.gof.rows=list("R2adj"=list.of.r2adj),
        file="anthocyanins_params.html")  
  
wordreg(list.of.compound.fit,single.row=T,ci.force=T,
        groups=list("ln(k1) [1/min]"=1,"Ea [KJ/Mol/K]"=2,"ln(k2) [1/min]"=3:6, "ln(k3) [1/min]"=7:10, 
                    "ln([AA]_0) [100ml/mg]"=11, "ln([DHAA]_0) [100ml/mg]" = 12),
        custom.gof.rows=list("R2adj"=list.of.r2adj),
        file="vitC_params.doc")
  