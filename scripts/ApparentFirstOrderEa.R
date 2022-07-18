# libraries
library(tidyr)
library(ggplot2)
library(texreg)
source("scripts/r2adj.R")
#read the data
presiones1.1 <- read.csv("data/presiones1.csv", sep = ";", dec = ",")
presiones.l <- gather(presiones1.1, compound, concentration,Delphinidin.3.O.sambubioside.5.O.glucoside:Ac.Dehidroascorbico,factor_key=TRUE)

# generate temperature 
presiones.l$Temp<-factor(substr(presiones.l$X,0,3))
levels(presiones.l$Temp)<-c(4,20,4,20,4,20)
presiones.l$Temp<-as.numeric(as.character(presiones.l$Temp))

#create factors for covariate modelling
presiones.l$sweetener<-factor(presiones.l$sweetener)
presiones.l$processing<-factor(presiones.l$processing)
presiones.l$rep<-factor(presiones.l$rep)
presiones.l$compound<-factor(presiones.l$compound,ordered=F)

# select anthocyanins
anthocyanins <- presiones.l[which(presiones.l$compound %in% c("Delphinidin.3.O.sambubioside.5.O.glucoside",	
                                                            "Delphinidin.3.5.O.diglucoside",
                                                            "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",	
                                                            "Delphinidin.3.O.sambubioside",	"Delphinidin.3.O.glucoside",	
                                                            "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside")),]

# first order apparent plot  to study the importance of Temperature
ggplot(anthocyanins, 
       aes(x = tiempo, y = log(concentration), group =sweetener:processing:factor(Temp),color=factor(Temp))) +
  facet_wrap(compound~.,scales="free")+
  geom_point()+ geom_smooth(method="lm",se = F)+
  xlab("Storage Time [Days]")+ylab("Concentration mg/100mL)")

# 1. the storage temperature is the most important factor in the study. 
# The differences between colours is much bigger than the differences between the lines of the same colour.
# 2. all the anthocyanins, except for the Delphinidin.3.5.O.diglucoside could reasonably be modelled with an apparent first order kinetic
library(nlme)
# the baseline model. Temperature dependence but processing or sweetener don't affect
ant.fm0<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=anthocyanins,
              param=list(C0~compound,
                         Cinf~compound,
                         lk~compound,
                         Ea~compound),
              start=c(C0=c(1,rep(0.001,5)),
                      Cinf=c(0,rep(0.001,5)),
                      lk=c(0.01,rep(0.001,5)),
                      Ea=c(10,rep(0.001,5)))
)
summary(ant.fm0,cor=F)
# With the experimental design we have we can estimate only main effects
# of the sweetener and processing variables.
# this model sweetener and processing are general effects (same for all compounds)
ant.fm1<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=anthocyanins,
              param=list(C0~compound+sweetener+processing,
                         Cinf~compound+sweetener+processing,
                         lk~compound+sweetener+processing,
                         Ea~1),
              start=c(C0=c(coef(ant.fm0)[1:6],rep(0.001,3)),
                      Cinf=c(coef(ant.fm0)[7:12],rep(0.001,3)),
                      lk=c(coef(ant.fm0)[13:18],rep(0.001,3)),
                      Ea=c(coef(ant.fm0)[19])
                      )
)
summary(ant.fm1)
# let's be more corageous and try to estimate different reaction rate constants between different compounds and processing conditions
ant.fm2<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
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
summary(ant.fm2)
# 4 or so significant parameters. Ea is only different in D.3.5.0.diglucoside
anova(ant.fm0,ant.fm1,ant.fm2)
# improvements in models are very significant and above noise.


ant.fm3<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
               data=anthocyanins,
               param=list(C0~compound+compound:sweetener+compound:processing,
                          Cinf~compound+compound:sweetener+compound:processing,
                          lk~compound+compound:sweetener+compound:processing,
                          Ea~compound),
               start=c(C0=c(coef(ant.fm0)[1:6],rep(0.001,18)),
                       Cinf=c(coef(ant.fm0)[7:12],rep(0.001,18)),
                       lk=c(coef(ant.fm0)[13:18],rep(0.001,18)),
                       Ea=c(coef(ant.fm0)[13:18])
               )
)
screenreg(ant.fm3,single.row=T,ci.force=T)

anova(ant.fm0,ant.fm1,ant.fm2,ant.fm3)
r2(ant.fm3)
# diagnostic plots
plot(ant.fm3)
qqnorm(ant.fm3)
plot(ant.fm3,resid(.)~Temp)
plot(ant.fm3,compound~resid(.))
plot(ant.fm3,processing~resid(.))
plot(ant.fm3,sweetener~resid(.))

# a nice table, maybe
list.of.compound.fit<-list()
list.of.r2adj<-list()
for (i in levels(factor(anthocyanins$compound))){
  print(i)
  list.of.compound.fit[[i]]<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
                                  data=subset(anthocyanins,compound==i),
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
        file="anthocyanins.html",
        caption.above=T,
        caption="Table X. Comparison of the apparent first order kinetic parameters and effects of processing and sweetener used",)

# Figure of Predictions
# first order apparent plot  to study the importance of Temperature
anthocyanins$grouping<-with(anthocyanins,sweetener:processing:factor(Temp))
ant.pred<-expand.grid(tiempo=seq(0,90,length=50),
                       compound=levels(factor(anthocyanins$compound)),
                       grouping=levels(factor(anthocyanins$grouping))
                       )
ant.pred$sweetener<-factor(with(ant.pred,substr(grouping,0,2)))
ant.pred$processing<-factor(with(ant.pred,substr(grouping,4,4)))
ant.pred$Temp<-as.numeric(with(ant.pred,substr(grouping,6,8)))
ant.pred$concentration<-predict(ant.fm3,newdata = ant.pred)
ggplot(anthocyanins, 
       aes(x = tiempo, y = concentration, col =grouping)) +
  facet_wrap(compound~.,scales="free")+
  geom_point()+geom_line(data=ant.pred,aes(x=tiempo,y=concentration,col=grouping))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")
ggsave(filename="Figure1.pdf")
