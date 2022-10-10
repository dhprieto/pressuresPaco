#create factors for covariate modelling
# libraries
library(tidyverse)
library(ggplot2)
library(texreg)
library(nlme)
source("scripts/r2adj.R")
library(tidyr)
library(nlraa)

# read the data
presiones1.1 <- read.csv("data/presiones1.csv", sep = ";", dec = ",")
presiones.l <- gather(presiones1.1, compound, concentration,Delphinidin.3.O.sambubioside.5.O.glucoside:Ac.Dehidroascorbico,factor_key=TRUE)

# generate temperature 
presiones.l$Temp<-factor(substr(presiones.l$X,0,3))
levels(presiones.l$Temp)<-c(4,20,4,20,4,20)
presiones.l$Temp<-as.numeric(as.character(presiones.l$Temp))

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
#model 0

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


#model 3

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

anthocyanins$grouping<-with(anthocyanins,sweetener:processing:factor(Temp))

ant.simul2<-expand.grid(tiempo=seq(0,90,length=30),
                        Temp = c(20,4),
                        compound=levels(factor(anthocyanins$compound)),
                        grouping=levels(factor(anthocyanins$grouping))
)



ant.simul2$sweetener<-factor(with(ant.simul2,substr(grouping,0,2)))
ant.simul2$processing<-factor(with(ant.simul2,substr(grouping,4,4)))
ant.simul2$grouping<-with(ant.simul2,sweetener:processing:factor(Temp))
ant.simul2$prediction<-predict(ant.fm3,newdata = ant.simul2)

# añadimos confidence bands

fm3.dm <- predict_gnls(ant.fm3, newdata = ant.simul2, interval = "conf")
ant.simul2.dm <- cbind(ant.simul2, fm3.dm)


ggplot(ant.simul2, 
       aes(x = tiempo, y = prediction, col = factor(Temp))) +
  facet_grid(processing+sweetener~factor(compound),scales="free")+
  geom_point()+geom_line()+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")
# ggsave(filename="Figure1.pdf")




# ST:1

ggplot(data = subset(ant.simul2.dm,grouping == "ST:1:20"), aes(x=tiempo, y = prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+
  geom_point(data=subset(anthocyanins,grouping == "ST:1:4"), aes(x=tiempo, y=concentration, col = factor(Temp)))+
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("ST:1 20º modelled | 4º experimental") + xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")

# ST:2
ggplot(data = ant.simul2.dm %>% filter(grouping == "ST:2:4"), aes(x=tiempo, y = prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+
  geom_point(data=anthocyanins %>% filter(grouping == "ST:2:20"), aes(x=tiempo, y=concentration, col = factor(Temp)))+
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("SU:1 20º experimental | 4º modelled")+ xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")

# ST:P
ggplot(data = ant.simul2.dm %>% filter(grouping == "ST:P:20"), aes(x=tiempo, y = prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+geom_point(data=anthocyanins %>% filter(grouping == "ST:P:4"), aes(x=tiempo, y=concentration, col = factor(Temp)))+
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("ST:P 20º modelled | 4º experimental")+ xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")

# SU:1

ggplot(data = ant.simul2.dm %>% filter(grouping == "SU:1:4"), aes(x=tiempo, y = prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+
  geom_point(data=anthocyanins %>% filter(grouping == "SU:1:20"), aes(x=tiempo, y=concentration, col = factor(Temp))) +
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("SU:1 20º experimental | 4º modelled")+ xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")

# SU:2

ggplot(data = ant.simul2.dm %>% filter(grouping == "SU:2:20"), aes(x=tiempo, y =prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+
  geom_point(data=anthocyanins %>% filter(grouping == "SU:2:4"), aes(x=tiempo, y=concentration, col = factor(Temp))) +
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("SU:2 20º modelled | 4º experimental")+ xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")

# SU:P

ggplot(data = ant.simul2.dm %>% filter(grouping == "SU:P:4"), aes(x=tiempo, y = prediction, col = factor(Temp)))+
  facet_wrap(compound~.,scales="free")+
  geom_point()+
  geom_point(data=anthocyanins %>% filter(grouping == "SU:P:20"), aes(x=tiempo, y=concentration, col = factor(Temp))) +
  geom_errorbar(aes(x = tiempo, ymin= Q2.5, ymax = Q97.5), alpha=0.2)+
  ggtitle("SU:P 20º experimental | 4º modelled")+ xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")



# pruebas


data(Orange)

fmoL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)

fmo <- nlme(fmoL, random = pdDiag(Asym + xmid + scal ~ 1))

sdat10 <- simulate_nlme(fmo, nsim = 100, psim = 1, level = 0, value = "data.frame")

ggplot(data = sdat10) + 
  geom_line(aes(x = age, y = sim.y, group = ii), color = "gray", alpha = 0.5) + 
  geom_point(aes(x = age, y = circumference)) + 
  ggtitle("psim = 1, level = 0") + 
  ylab("circumference")

set.seed(123)

fmoL <- nlsList(prediction ~ compound+Temp, data = ant.simul2)


ant_sim_unc <- simulate_gnls(ant.fm3, newdata = ant.simul2, nsim = 100, psim = 1, level = 1, value = "data.frame")
ant.simul2$sim <- ant_sim_unc

#

