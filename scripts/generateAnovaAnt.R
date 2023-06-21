# libraries ----

library(tidyr)
library(ggplot2)
library(nlraa)
library(nlme)
library(rstatix)

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


anova(ant.fm0, ant.fm3)


# simulation ----

ant.SIM <- expand.grid(tiempo=seq(-1,90,length=50),
                        compound=levels(factor(anthocyanins$compound)),
                        Temp=c(4,20),
                        processing=levels(factor(anthocyanins$processing)),
                        sweetener=levels(factor(anthocyanins$sweetener))
)

yhatSIM <- predict_gnls(ant.fmf, newdata = ant.SIM, interval = "confidence")

ant.SIM$yhat<-yhatSIM[,"Estimate"]
ant.SIM$Q2.5<-yhatSIM[,"Q2.5"]
ant.SIM$Q97.5<-yhatSIM[,"Q97.5"]

# anova ----

ant.SIM.w <- ant.SIM %>% pivot_wider(names_from = "compound", 
                                     values_from =c("concentration","Q2.5", "Q97.5"), names_vary = "slowest")

ant.SIM.w$id <- as.character(with(data = ant.SIM.w, factor(Temp):processing:sweetener))

anthocyanins.w <- anthocyanins %>% pivot_wider(names_from="compound", values_from = "concentration")
anthocyanins.w$id <- as.character(with(data = anthocyanins.w, rep:factor(Temp):processing:sweetener))

# rawExperimentalData

write.anova.table <- function(datos, compuesto){
  get_anova_table(anova_test(data = datos, dv=compuesto, wid=id, 
                             between = c(processing, sweetener), within= c(tiempo)), correction = "auto" )
}

# simulated data

write.anova.table.unpaired <- function(datos, compuesto){
  get_anova_table(anova_test(datos[,compuesto]~tiempo*processing*sweetener, data=datos), correction = "auto" )
}

write.anova.table.unpaired(and.SIM.W, )
ant.SIM.Anova <- list()
for (i in c("concentration_Delphinidin.3.O.sambubioside.5.O.glucoside",	
            "concentration_Delphinidin.3.5.O.diglucoside",
            "concentration_Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",	
            "concentration_Delphinidin.3.O.sambubioside",	"concentration_Delphinidin.3.O.glucoside",	
            "concentration_Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside")){
  
  ant.SIM.Anova <- append(ant.SIM.Anova, c(i,write.anova.table.unpaired(ant.SIM.w, i)))
  
}


formulas_anova <- lapply(list("concentration_Delphinidin.3.O.sambubioside.5.O.glucoside",	
                           "concentration_Delphinidin.3.5.O.diglucoside",
                           "concentration_Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",	
                           "concentration_Delphinidin.3.O.sambubioside",	"concentration_Delphinidin.3.O.glucoside",	
                           "concentration_Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside"),
      function(x) formula(paste0(x,"~tiempo*sweetener*processing"),env=globalenv()))       


ant.SIM.Anova.unpaired <- list()
for (i in seq(1:6)){
  
  ant.SIM.Anova.unpaired <- append(ant.SIM.Anova.unpaired, 
                                   c(as.character(formulas_anova[[i]][[2]]),
                                     get_anova_table(anova_test(formulas_anova[[i]], data=ant.SIM.w), 
                                                     correction = "auto" )))
  
  
  }


write.csv(x=ant.SIM.Anova.unpaired, file="results/simAntAnova_unpaired.csv")

# graphics ----

plot(x=ant.SIM$tiempo[which(ant.SIM$compound=="Delphinidin.3.O.sambubioside.5.O.glucoside")], 
     y=ant.SIM$yhat[which(ant.SIM$compound=="Delphinidin.3.O.sambubioside.5.O.glucoside")], type = "l")

ggplot(dplyr::filter(ant.SIM, compound=="Delphinidin.3.O.sambubioside.5.O.glucoside"), aes(x=tiempo, y=yhat))+
  geom_line()

View(dplyr::filter(ant.SIM, compound=="Delphinidin.3.O.sambubioside.5.O.glucoside"))

ggplot(ant.SIM, aes(x=tiempo, y= yhat, col = processing)) +
  facet_wrap(compound[2]~factor(Temp)+sweetener) +
  geom_line()+
  geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = yhat, fill = processing)) +
  geom_point(data=anthocyanins, aes(x=tiempo, y = concentration, fill =processing))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of anthocyanins by processing", 
          subtitle = "Dots are experimental points")

ggsave("antxprocessing.png", plot = p1 )

p1
ggplot(ant.SIM, aes(x=tiempo, y= yhat, col = sweetener)) +
  facet_grid(compound~factor(Temp)+processing, scales ="free") +
  geom_line()+
  geom_ribbon(aes(ymax=Q97.5, ymin=Q2.5, y = yhat, fill = sweetener), alpha= 0.1) +
  geom_point(data=anthocyanins, aes(x=tiempo, y = concentration, fill =sweetener))+
  xlab("Storage Time [Days]")+ylab("Concentration [mg/100mL]")+
  ggtitle("Degradation of anthocyanins by sweetener", 
          subtitle = "Dots are experimental points")


