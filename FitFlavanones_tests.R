# Fit Flavanones ----
# libraries ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(texreg)
library(nlme)
library(nlraa)
library(nls2)
library(aomisc)

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

presiones.l$tiempo[presiones.l$tiempo == 0] <- 0.1 
presiones.l$tiempo[presiones.l$tiempo == -1] <- 0.01

# select flavanones 
flavanones <-  filter(presiones.l, compound == c("O.trygycosil.naringenin", "Eriodyctiol.7.O.rutinoside", 
                                                 "Naringenin..7.O.rutinoside","Hespertin..7.O.rutinoside"))

flavanones$tiempoPLus <- flavanones$tiempo+1.001

# first order apparent plot to study the importance of Temperature

ggplot(flavanones, 
       aes(x = tiempo, y = log(concentration), group =sweetener:processing:factor(Temp),color=factor(Temp))) +
  facet_wrap(compound~.,scales="free")+
  geom_point()+ geom_smooth(method="loess",se = F)+
  xlab("Storage Time [Days]")+ylab("Concentration mg/100mL)")


getInitial(concentration ~ SSweibull(tiempoPLus, Asym, Drop, lrc, pwr), data = flavanones)


flav.weibull <- nls2(concentration ~ SSweibull(tiempoPLus, Asym, Drop, lrc, pwr), data = flavanones)


fit1 <- nls(concentration ~ SSweibull(tiempo+1.1, Asym, Drop, lcr, pwr), data = flavanones)


# ajuste gnls ----


model <- drm(concentration ~ tiempoPLus, fct = W2.4(), data = flavanones)

coef(model)

W2.4()


model.LL3<- drm(concentration~tiempoPLus, data=flavanones, fct=LL.3(fixed=c(NA, 3, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(model.LL3, fctList = list(W1.3(fixed=c(NA, 3, NA)),W1.4(), W2.3(fixed=c(NA, 3, NA)), W2.4(),  LL.4()),linreg=TRUE) 

model.W2.3<- drm(concentration~tiempoPLus, data=flavanones, fct=W2.3(fixed=c(NA, 3, NA), names = c("Slope", "Upper Limit", "ED50")))


ED(model.W2.3, c(1,90), interval = "delta")


stats::getInitial(concentration ~ SSlogis5(tiempoPLus, asym1,asym2,xmid,iscal,theta), data=flavanones)

flav.fm15 <- gnls(concentration~W2.4(fixed=c(1.1586e-02,NA, 1.4334, NA),names = c("b", "c", "d", "e")),
                  data=flavanones,
                  start=c(b=3.4928144,c=3.7674222,
                          d=53.5916814,e=-145.0211324),
                  control = gnlsControl(nlsTol = 100, maxIter = 200),verbose=T)

model2 <- W2.4(fixed=c(1.1586e-02,NA, 1.4334, NA), 
               names = c("b", "c", "d", "e"))


summary(model)
plot(model)
coef(model)


W2.4()

flav.fm0<-gnls(concentration~weibull2(tiempoPLus, c, e),
               data=flavanones,
               param=list(c~compound,e~1),
               start=c(c=c(4, rep(0.001,2)),
                       e=c(2.3e+71,0.001)),
               control = gnlsControl(nlsTol = 5000, maxIter = 200),
      
               weights=varIdent(~1|compound),
               verbose=T#,rep(0.001,3)),
                       )#,rep(0.001,3)))

# test hill ----

library(minpack.lm)

flav.nls <- nlsLM(concentration ~ SShill3(tiempoPLus, Ka, n,a), data = flavanones,
                  control = c(maxiter = 500))

summary(flav.nls)

flav.fm0.Hill<-gnls(concentration~SShill3(tiempoPLus, Ka, n, a), data = flavanones,
                    param=list(Ka~1,
                               n~1,
                               a~1),
                    start=coef(flav.nls),
                    control = gnlsControl(nlsTol = 500, maxIter = 200))


flav.fm1.Hill<-gnls(concentration~SShill3(tiempoPLus, Ka, n, a), data = flavanones,
                    param=list(Ka~compound,
                               n~compound,
                               a~compound),
                    start= c(Ka=c(coef(flav.nls)[1],rep(0.001,3)),
                             n=c(coef(flav.nls)[2],rep(0.001,3)),#,rep(0.001,3)),
                             a=c(coef(flav.nls)[3],rep(0.001,3))),
                    control=gnlsControl(nlsTol = 50, maxIter = 200))

summary(flav.fm0.Hill)




# tests SSweibull ----

flav_weibull <- flavanones %>% filter(compound == c("O.trygycosil.naringenin","Eriodyctiol.7.O.rutinoside"))

## NaNs produced

fit <- nls(concentration ~ SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
                          data = flavanones)

## NaNs produced


flav.fm0<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
               data=flavanones,
               start=c(Asym=c(160),#,rep(0.001,3)),
                       Drop=c(160),#,rep(0.001,3)),
                       lrc=c(-15.5),#,rep(0.001,3)),
                       pwr=c(5)),
               control = gnlsControl(nlsTol = 500, maxIter = 200)#,rep(0.001,3)))
)

anova(flav.fm0)
coef(flav.fm0)

flav.fm1<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
               data=flavanones,
               param=list(Asym~1,
                          Drop~1,
                          lrc~1,
                          pwr~1),
               start=c(Asym=c(3.5759363  , rep(0.001,0)),#,rep(0.001,3)),
                       Drop=c(-0.1944744),#,rep(0.001,3)),
                       lrc=c(-15.5),#,rep(0.001,3)),
                       pwr=c(5)),
               control = gnlsControl(nlsTol = 500, maxIter = 200)#,rep(0.001,3)))
)

anova(flav.fm1)
anova(flav.fm0,flav.fm1)
coef(flav.fm0)

flav.fm2<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
               data=flavanones,
               param=list(Asym~compound+processing,
                          Drop~Temp,
                          lrc~1,
                          pwr~1),
               start=c(Asym=c(300.5),rep(0.001,5),
                       Drop=c(-200),rep(0.1,1),
                       lrc=c(-15.5),rep(0.1,0),#,rep(0.001,3)),
                       pwr=c(3),rep(0.1,0)),
               control = gnlsControl(nlsTol = 500, maxIter = 200),
               verbose=T#,rep(0.001,3)))
)

anova(flav.fm2)



flav.fm3<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
               data=flavanones,
               param=list(Asym~compound+processing+sweetener,
                          Drop~Temp,
                          lrc~1,
                          pwr~1),
               start=c(Asym=c(300,rep(0.001,6)),
                       Drop=c(-200, rep(0.001,1)),
                       lrc=c(-15.5, rep(0.001,0)),
                       pwr=c(5, rep(0.001,0))),#rep(0.001,3)),
               control = gnlsControl(nlsTol = 500, maxIter = 2000),
               weights=varIdent(~1|compound),
               verbose=T)

flav.fm4<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc, pwr),
               data=flavanones,
               param=list(Asym~compound+processing+sweetener,
                          Drop~Temp,
                          lrc~1,
                          pwr~1),
               start=c(Asym=c(300,rep(0.001,6)),
                       Drop=c(-300, rep(0.001,1)),
                       lrc=c(-20, rep(0.001,0)),
                       pwr=c(7, rep(0.001,0))),#rep(0.001,3)),
               control = gnlsControl(nlsTol = 500, maxIter = 2000),
               weights=varIdent(~1|compound),
               verbose=T)


flav.fm5<-gnls(concentration~SSweibull(tiempoPLus, Asym, Drop, lrc=-19.84387002, pwr=6.94345611),
               data=flavanones,
               param=list(Asym~compound+processing,
                          Drop~compound:Temp),
               start=c(Asym=c(1.24463356,1.68524659,1.47584541,5.62738908, -0.11354167,
                              0.38479167, rep(0.001,0)),
                       Drop=c(0.15378739, -0.03367208, rep(0.001,3))),
               control = gnlsControl(nlsTol = 100, maxIter = 2000),
               weights=varIdent(~1|compound),
               verbose=T)

coef(flav.fm5)
anova(flav.fm5)
anova(flav.fm0,flav.fm1,flav.fm2,flav.fm3,flav.fm4, flav.fm5)

coef(flav.fm4)

flav.sim<-rbind(flavanones,flavanones,flavanones,flavanones)
flav.sim$tiempoPLus<-runif(144*4,0.001,91.001)

flav.sim$yhat<-predict(flav.fm4,newdata=flav.sim)

ggplot(flavanones, aes(x=tiempoPLus, y= concentration, col = processing)) +
  facet_grid(compound~factor(Temp)+sweetener, scales ="free") +
  geom_point()+
  geom_line(data=flav.sim,aes(y=yhat))
  

#geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = yhat, fill = processing), alpha= 0.1) +
  geom_point(data=flavanones, aes(x=tiempo, y = concentration, fill =processing))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of flavanones by processing", 
          subtitle = "Dots are experimental points")
# test ----

flav.SIM <- expand.grid(tiempoPLus=seq(0,90,length=50)+1.1,
                        compound=levels(factor(flavanones$compound)),
                        Temp=c(4,20),
                        processing=levels(factor(flav_weibull$processing)),
                        sweetener=levels(factor(flav_weibull$sweetener))
)





yhatSIM<-predict_gnls(flav.fm5,newdata=flav.SIM,interval="prediction")

flav.SIM$yhat<-yhatSIM[,"Estimate"]
flav.SIM$Q2.5<-yhatSIM[,"Q2.5"]
flav.SIM$Q97.5<-yhatSIM[,"Q97.5"]

predict2<-predict(flav.fm5, flav.SIM)

ggplot(flav.SIM, aes(x=tiempoPLus, y= predict2, col = processing)) +
  facet_grid(compound~factor(Temp)+sweetener, scales ="free") +
  geom_line()+
  #geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = yhat, fill = processing), alpha= 0.1) +
  geom_point(data=flavanones, aes(x=tiempo, y = concentration, fill =processing))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of flavanones by processing", 
          subtitle = "Dots are experimental points")

ggplot(flav.SIM, aes(x=tiempoPLus, y= predict2, col = sweetener)) +
  facet_grid(compound~factor(Temp)+processing, scales ="free") +
  geom_line()+
  #geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = yhat, fill = processing), alpha= 0.1) +
  geom_point(data=flavanones, aes(x=tiempo, y = concentration, fill =sweetener))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of flavanones by processing", 
          subtitle = "Dots are experimental points")


library(nlme)
flavanones$ExpInd<-with(flavanones,factor(compound:processing:sweetener:factor(Temp)))

flav.g<-groupedData(concentration~tiempoPLus|ExpInd,data=flavanones)
plot(flav.g)
flav.fml<-nlsList(concentration~SSbg4rp(tiempoPLus,w.max,lt.e,ldtm,ldtb),
                  data=flav.g)

flav.fml
plot(flav.fml)
augPred(flav.fml)


flav_nls_fit_func <- function(flav_expind) {
  nls.multstart::nls_multstart(concentration~SSbg4rp(tiempoPLus,w.max,lt.e,ldtm,ldtb),
                               data = flav_expind,
                               lower=c(w=0, w.max=0, lt.e=0,ldtm=0,ldtb=0),
                               upper=c(w=Inf, w.max=Inf, lt.e=Inf,ldtm=Inf,ldtb=Inf),
                               iter=500,
                               supp_errors = "Y")
}

flavdata<-flavanones
library(purrr)
library(nls.multstart)
flavdata <- flavdata %>% 
  mutate(flav_nls_fit = map(ExpInd, ~flav_nls_fit_func(.x)))

