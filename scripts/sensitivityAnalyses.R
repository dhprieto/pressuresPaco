# sensitivity analysis

library(tidyr)
library(ggplot2)
library(nlraa)
library(nlme)
library(dplyr)
library(data.table)
library(deSolve)
library(rmutil)
#library(sensobol)
#library(sensemakr)
library(sensitivity)


# Anthocyanins ----

anthocyanins.norm <- anthocyanins
anthocyanins.norm$concentration <- scales::rescale(anthocyanins.norm$concentration) 

ant.fm4<- gnls(model = concentration ~ Cinf + (C0 - Cinf) * exp(-exp(lk - 
                                                                       Ea/0.008314 * (1/(Temp + 273) - 1/(16 + 273))) * tiempo), 
               data = anthocyanins, params = list(C0 ~ compound + compound:sweetener + 
                                                    compound:processing, Cinf ~ compound + compound:sweetener + 
                                                    compound:processing, lk ~ compound + compound:sweetener + 
                                                    compound:processing, Ea ~ 1), start = c(C0 = c(coef(ant.fm0)[1:6], 
                                                                                                   rep(0.001, 18)), Cinf = c(coef(ant.fm0)[7:12], rep(0.001, 
                                                                                                                                                      18)), lk = c(coef(ant.fm0)[13:18], rep(0.001, 18)), Ea = c(coef(ant.fm0)[19], 
                                                                                                                                                                                                               rep(0.001, 0))))
anova(ant.fm4)  

tiempo.Ant <- rep(unique(anthocyanins$tiempo),125)
names.compounds.Ant <- names(list.of.compound.fit)
names.compounds <- c("Delphinidin.3.O.sambubioside.5.O.glucoside","Delphinidin.3.5.O.diglucoside",
                     "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",
                     "Delphinidin.3.O.sambubioside","Delphinidin.3.O.glucoside",
                     "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside")



n <- 1000
X1 <- data.frame(matrix(c(runif(4*n), rep(c(4,20), 500),tiempo.Ant,
                          c(rep(names.compounds,166),names.compounds[1:4]), 
                          rep(c("SU","ST"), 500),c(rep(c("1", "2", "P"), 333), "1")),nrow=n))


colnames(X1) <- c("C0", "Cinf", "lk", "Ea", "Temp","tiempo","compound", "sweetener", "processing")

X1 <- X1 %>% mutate_at(c("C0", "Cinf", "lk", "Ea", "Temp","tiempo"), as.numeric)

X2 <- data.frame(matrix(c(runif(4*n), sample(rep(c(4,20), 500)),sample(tiempo.Ant),
                          sample(c(rep(names(list.of.compound.fit),166),names(list.of.compound.fit)[1:4])), 
                          sample(rep(c("SU","ST"), 500)),sample(c(rep(c("1", "2", "P"), 333), "1"))),nrow=n))

colnames(X2) <- c("C0", "Cinf", "lk", "Ea", "Temp","tiempo","compound", "sweetener", "processing")

X2 <- X2 %>% mutate_at(c("C0", "Cinf", "lk", "Ea", "Temp","tiempo"), as.numeric)


ant.fm3.sobol.martinez <- sobolmartinez(model = ant.fm4, X1,X2)
ant.fm3.sobol.mara <- sobolmara(model = ant.fm4, X1)

ggplot(ant.fm3.sobol.martinez) + ggtitle("Sobol sensitivity test for Anthocyanins model")

# pruebas plot ----

y <- ant.fm3.sobol.martinez$S[-c(6,7),]

x <- colnames(t(y))

ggplot(y, aes(x,y[,1]))+geom_point()



# Vitamin C ----
VitC.l$tiempo <- VitC.l$tiempo+1.0001 

VitC.fmF<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = Species,
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
                       laa0=15,
                       ldhaa0=2))


tiempo.VitC <- rep(unique(VitC.l$tiempo),125)
names.compounds <- rep(unique(VitC.l$Species), 500)



n <- 1000
X1 <- data.frame(matrix(c(runif(6*n), rep(c(4,20), 500),tiempo.VitC,
                          names.compounds, 
                          rep(c("SU","ST"), 500),
                          c(rep(c("1", "2", "P"), 333), "1")),nrow=n))


colnames(X1) <- c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo","Species", "sweetener", "processing")

X1 <- X1 %>% mutate_at(c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo"), as.numeric)

X1$Species[X1$Species == "1"] <-"Ascorbic" 
X1$Species[X1$Species == "2"] <-"Dehydroascorbic"

X2 <- data.frame(matrix(c(runif(6*n), sample(rep(c(4,20), 500)),sample(tiempo.VitC),
                          sample(names.compounds), 
                          sample(rep(c("SU","ST"), 500)),sample(c(rep(c("1", "2", "P"), 333), "1"))),nrow=n))

colnames(X2) <- c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo","Species", "sweetener", "processing")

X2 <- X2 %>% mutate_at(c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo"), as.numeric)

X2$Species[X2$Species == "1"] <-"Ascorbic" 
X2$Species[X2$Species == "2"] <-"Dehydroascorbic"


VitC.fmF.sobol.martinez <- sobolmartinez(model = VitC.fmF, X1,X2)
VitC.fm3.sobol.mara <- sobolmara(model = VitC.fm4, X1)

ggplot(VitC.fmF.sobol.martinez) + ggtitle("Sobol sensitivity analysis for Vitamin C model")


# pruebas plot ----

y <- VitC.fm3.sobol.martinez$S[-c(6,7),]
x <- colnames(t(y))

ggplot(y, aes(x,y[,1]))+geom_point()


# salvando resultados

saveRDS(VitC.fmF.sobol.martinez, file ="results/VitCSobol.RDS")
saveRDS(ant.fm3.sobol.martinez, file ="results/AntSobol.RDS")

