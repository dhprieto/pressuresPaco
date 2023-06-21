# Fit Flavanones ----
# libraries ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(texreg)
library(nlme)
library(nlraa)
source("scripts/r2adj.R")

# read the data ----

presiones1.1 <- read.csv("data/presiones1.csv", sep = ";", dec = ",")
presiones.l <- gather(presiones1.1, compound, concentration,
                      Delphinidin.3.O.sambubioside.5.O.glucoside:Ac.Dehidroascorbico,factor_key=TRUE)

# generate temperature 

presiones.l$Temp<-factor(substr(presiones.l$X,0,3))
levels(presiones.l$Temp)<-c(4,20,4,20,4,20)
presiones.l$Temp<-as.numeric(as.character(presiones.l$Temp))

#create factors for covariate modelling

presiones.l$sweetener<-factor(presiones.l$sweetener)
presiones.l$processing<-factor(presiones.l$processing)
presiones.l$rep<-factor(presiones.l$rep)
presiones.l$compound<-factor(presiones.l$compound,ordered=F)
  
# select flavanones 
flavanones <-  filter(presiones.l, compound == c("O.trygycosil.naringenin", "Eriodyctiol.7.O.rutinoside", 
                                      "Naringenin..7.O.rutinoside","Hespertin..7.O.rutinoside"))


# first order apparent plot to study the importance of Temperature

ggplot(flavanones, 
       aes(x = tiempo, y = log(concentration), group =sweetener:processing:factor(Temp),color=factor(Temp))) +
  facet_wrap(compound~.,scales="free")+
  geom_point()+ geom_smooth(method="loess",se = F)+
  xlab("Storage Time [Days]")+ylab("Concentration mg/100mL)")

# 1. in this case, opposite to anthocyanins, the storage temperature does not seem to be the most important factor. 
# The differences between lines of the same colour are relevant. So, there must be another factor more important
# 2. Only the Herpetin-rutinoside could reasonably be modelled with an apparent first order kinetic

# the baseline model. Temperature dependence but processing or sweetener don't affect

flav.fm0<-gnls(concentration~Cinf+(C0-Cinf)*
                 exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=flavanones,
              param=list(C0~compound,
                         Cinf~compound,
                         lk~compound,
                         Ea~compound),
              start=c(C0=c(1,rep(0.001,3)),
                      Cinf=c(0,rep(0.001,3)),
                      lk=c(0.01,rep(0.001,3)),
                      Ea=c(10,rep(0.001,3))),
              control=gnlsControl(nlsTol = 50)
)

summary(flav.fm0,cor=F)





# until here it runs ----


# With the experimental design we have we can estimate only main effects
# of the sweetener and processing variables.
# this model sweetener and processing are general effects (same for all compounds)

flav.fm1<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=flavanones,
              params=list(C0~compound+sweetener+processing,
                         Cinf~compound+sweetener+processing,
                         lk~compound+sweetener+processing,
                         Ea~1),
              start=c(C0=c(coef(flav.fm0)[1:4],rep(0.001,3)),
                      Cinf=c(coef(flav.fm0)[5:8],rep(0.001,3)),
                      lk=c(coef(flav.fm0)[9:12],rep(0.001,3)),
                      Ea=c(coef(flav.fm0)[13])
              )
)
summary(flav.fm1)


# let's be more corageous and try to estimate different reaction rate constflavs between different compounds and processing conditions

flav.fm2<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=flavanones,
              param=list(C0~compound+sweetener+processing,
                         Cinf~compound+sweetener+processing,
                         lk~compound+sweetener+compound:processing,
                         Ea~compound),
              start=c(C0=c(coef(flav.fm0)[1:4],rep(0.001,3)),
                      Cinf=c(coef(flav.fm0)[5:8],rep(0.001,3)),
                      lk=c(coef(flav.fm0)[9:12],rep(0.001,7)),
                      Ea=c(coef(flav.fm0)[13],rep(0.001,5))
              )
)
summary(flav.fm2)

# 4 or so significflav parameters. Ea is only different in D.3.5.0.diglucoside

anova(flav.fm0,flav.fm1,flav.fm2)

# improvements in models are very significflav and above noise.

flav.fm3<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=flavanones,
              param=list(C0~compound+compound:sweetener+compound:processing,
                         Cinf~compound+compound:sweetener+compound:processing,
                         lk~compound+compound:sweetener+compound:processing,
                         Ea~compound),
              start=c(C0=c(coef(flav.fm0)[1:6],rep(0.001,18)),
                      Cinf=c(coef(flav.fm0)[7:12],rep(0.001,18)),
                      lk=c(coef(flav.fm0)[13:18],rep(0.001,18)),
                      Ea=c(coef(flav.fm0)[13:18])
              )
)
screenreg(flav.fm3,single.row=T,ci.force=T)

anova(flav.fm0,flav.fm1,flav.fm2,flav.fm3)
r2(flav.fm3)
# diagnostic plots
plot(flav.fm3)
qqnorm(flav.fm3)
plot(flav.fm3,resid(.)~Temp)
plot(flav.fm3,compound~resid(.))
plot(flav.fm3,processing~resid(.))
plot(flav.fm3,sweetener~resid(.))

# a nice table, maybe
list.of.compound.fit<-list()
list.of.r2adj<-list()
for (i in levels(factor(flavanones$compound))){
  print(i)
  list.of.compound.fit[[i]]<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
                                  data=subset(flavanones,compound==i),
                                  param=list(C0~sweetener+processing,
                                             Cinf~sweetener+processing,
                                             lk~sweetener+processing,
                                             Ea~1),
                                  start=c(C0=c(5.05,rep(0.001,3)),
                                          Cinf=c(0,rep(0.001,3)),
                                          lk=c(-4.38,rep(0.001,3)),
                                          Ea=c(67)
                                  )
  )
  list.of.r2adj[[i]]<-round(r2(list.of.compound.fit[[i]])[2],3)
}

screenreg(list.of.compound.fit,single.row=T,ci.force=T)

htmlreg(list.of.compound.fit,single.row=T,ci.force=T,
        groups=list("C0 [mg/100mL]"=1:4,"Cinf [mg/100mL]"=5:8,"ln(k) [1/min]"=9:12,"Ea [KJ/Mol/K]"=13),
        custom.gof.rows=list("R2adj"=list.of.r2adj),
        file="flavanones.html",
        caption.above=T,
        caption="Table X. Comparison of the apparent first order kinetic parameters and effects of processing and sweetener used",)


