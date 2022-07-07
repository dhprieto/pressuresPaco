# Visualización resultados experimento HPP Paco ----

# Carga de librerías ----

library(tidyverse)

# Pruebas ----

presiones1 <- read.csv("../presiones//presiones1.csv", sep = ";", dec = ",")

plot(presiones1)

presiones1 <- presiones1 %>% select(X, tiempo, mg.Ac..Ascórbico..AA..100ml.de.zumo, mg.Ac..Dehidroascórbico.100ml.de.zumo)

presiones1 <- na.omit(presiones1)

plot(y = presiones1$mg.Ac..Ascórbico..AA..100ml.de.zumo, x = presiones1$tiempo)



presiones1 <- read.csv("../presiones//presiones1.csv", sep = ";", dec = ",")

for (i in seq(1, nrow(presiones1))){

  if (grepl("ST1", presiones1$X[i])){
    presiones1$tratamiento[i] <- "presion_1"
  }
  
  if (grepl("ST2", presiones1$X[i])){
    presiones1$tratamiento[i] <- "presion_2"
  }
  
  if (grepl("STP", presiones1$X[i])){
    presiones1$tratamiento[i] <- "pasteurizacion"
  }
  
  if (grepl("ST", presiones1$X[i])){
    presiones1$endulzante[i] <- "ST"
  }
  
  if (grepl("SU", presiones1$X[i])){
    presiones1$endulzante[i] <- "SU"
  }
  
}


# Comienzo del proceso ----

# Lectura fichero

presiones1.1 <- read.csv("../presiones/presiones1.csv", sep = ";", dec = ",")

# Loop para añadir temperatura

for (i in seq(1, nrow(presiones1.1))){
  
  if (grepl("ST1", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "4"
  }
  
  if (grepl("ST2", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "20"
  }
  
  if (grepl("STP", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "4"
  }
  
  if (grepl("SU1", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "20"
  }
  
  if (grepl("SU2", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "4"
  }
  
  if (grepl("SUP", presiones1$X[i])){
    presiones1.1$temperatura[i] <- "20"
  }
  
}

# from wide to long : http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/

presiones.l <- gather(presiones1.1, compound, concentration,Delphinidin.3.O.sambubioside.5.O.glucoside:Ac.Dehidroascorbico, 
                      factor_key=TRUE)

# Gráficas preliminares sin temperatura
ggplot(subset(presiones.l, compound != "Ac.Ascorbico" & compound != "Ac.Dehidroascorbico"), 
  aes(x = tiempo, y = concentration, colour = compound)) +
  facet_grid(sweetener~processing)+
  geom_jitter()+ geom_smooth(se = F)

ggplot(subset(presiones.l, compound != "Ac.Ascorbico" & compound != "Ac.Dehidroascorbico"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_wrap(processing~compound, scales = "free", nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)


ggplot(subset(presiones.l, compound != "Ac.Ascorbico" & compound != "Ac.Dehidroascorbico"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_wrap(sweetener~compound,scales = "free", nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)


# Gráficas siguientes, añadiendo temperatura ----

# Agrupación en variable de las familias de compuestos

flavanones <- c("Delphinidin.3.O.sambubioside.5.O.glucoside",	"Delphinidin.3.5.O.diglucoside",
                "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",	
                "Delphinidin.3.O.sambubioside",	"Delphinidin.3.O.glucoside",	
                "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside")

anthocyanins <- c("O.trygycosil.naringenin","Eriodyctiol.7.O.rutinoside",	
                  "Naringenin..7.O.rutinoside","Hespertin..7.O.rutinoside")

ascorbic <- c("Ac.Ascorbico", "Ac.Dehidroascorbico")

## Gráficas flavanonas 20º ----

# Seleccion familia
presiones.flav <- presiones1.1 %>% select(-anthocyanins, -all_of(ascorbic))

# from wide to long 

presiones.flav.l <- gather(presiones.flav, compound, concentration, 
                           Delphinidin.3.O.sambubioside.5.O.glucoside:Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside,
                           factor_key=TRUE)
# Graficado 
ggplot(subset(presiones.flav.l, temperatura == "20"), 
  aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Flavanones Sweetener ~ compound 20º")

ggplot(subset(presiones.flav.l, temperatura == "20"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Flavanones Processing ~ compound 20º")

## Gráficas flavanones 4º ----

ggplot(subset(presiones.flav.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Flavavones Processing ~ compound 4º")

ggplot(subset(presiones.flav.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Flavanones Sweetener ~ compound 4º")


## Gráficas anthocyanins 20º ----

# Selección familia
presiones.ant <- presiones1.1 %>% select(-flavanones, -all_of(ascorbic))

# from wide to long 

presiones.ant.l <- gather(presiones.ant, compound, concentration, 
                             O.trygycosil.naringenin:Hespertin..7.O.rutinoside,
                             factor_key=TRUE)

# Graficado

ggplot(subset(presiones.ant.l, temperatura == "20"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Anthocyanins Sweetener ~ compound 20º")

ggplot(subset(presiones.ant.l, temperatura == "20"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Anthocyanins Processing ~ compound 20º")

## antocianos 4º ----

ggplot(subset(presiones.ant.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Anthocyanins Processing ~ compound 4º")

ggplot(subset(presiones.ant.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Anthocyanins Sweetener ~ compound 4º")


## ascórbicos 20 º ----

# Selección familia
presiones.asc <- presiones1.1 %>% select(-flavanones, -all_of(anthocyanins))

# from wide to long
presiones.asc.l <- gather(presiones.asc, compound, concentration, 
                          Ac.Ascorbico:Ac.Dehidroascorbico,
                          factor_key=TRUE)

# Graficado

ggplot(subset(presiones.asc.l, temperatura == "20"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Ascorbic Sweetener ~ compound 20º")

ggplot(subset(presiones.asc.l, temperatura == "20"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Ascorbic Processing ~ compound 20º")

## ascórbicos 4º ----


# Graficado 
ggplot(subset(presiones.asc.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = sweetener)) +
  facet_grid(processing~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Ascorbic Processing ~ compound 4º")

ggplot(subset(presiones.asc.l, temperatura == "4"), 
       aes(x = tiempo, y = concentration, colour = processing)) +
  facet_grid(sweetener~compound,scales = "free")+ #nrow = 3, ncol = 10)+
  geom_jitter()+ geom_smooth(se = F)+ ggtitle("Ascorbic Sweetener ~ compound 4º")

