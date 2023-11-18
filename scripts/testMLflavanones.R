# libraries ----

library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)
library(xgboost)

# functions ----

normalizingNumeric <- function(tableComplete) {
  
  for (i in colnames(tableComplete)){
    
    if (is.numeric(tableComplete[,i]) && i != "numVol"){
      
      tableComplete[,i] <- scales::rescale(tableComplete[,i])
      
    }
    
  }
  return (tableComplete)    
}



# read the data ----

presiones1.1 <- read.csv("data/presiones1.csv", sep = ";", dec = ",")
flavanones <-  select(presiones1.1, -c("Delphinidin.3.O.sambubioside.5.O.glucoside",                          
                                       "Delphinidin.3.5.O.diglucoside",                                     
                                       "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",
                                       "Delphinidin.3.O.sambubioside",                                        
                                       "Delphinidin.3.O.glucoside",                                           
                                       "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside",
                                       "Ac.Ascorbico",                                                        
                                       "Ac.Dehidroascorbico")) 
                      
flavanones$Temp<-factor(substr(flavanones$X,0,3))
levels(flavanones$Temp)<-c(4,20,4,20,4,20)
flavanones$Temp<-as.character(flavanones$Temp)
flavanones$tiempo <- as.character(flavanones$tiempo)
flavanonesNorm <- normalizingNumeric(flavanones)

flavanonesNorm.t <- flavanonesNorm %>% filter(tiempo != "-1")

# predictors

set.seed(123)

ctrl.reg <- trainControl(method="repeatedcv", repeats = 3)

tune_grid <- expand.grid(nrounds = 200,
                         max_depth = 5,
                         eta = 0.05,
                         gamma = 0.01,
                         colsample_bytree = 0.75,
                         min_child_weight = 0,
                         subsample = 0.5)

xgbFit.3gN <- train(O.trygycosil.naringenin ~ rep+processing+sweetener+tiempo+Temp, data = flavanonesNorm.t, method = "xgbTree",
                    trControl=ctrl.reg,
                    tuneGrid = tune_grid,
                    tuneLength = 20)

xgbFit.E7O <- train(Eriodyctiol.7.O.rutinoside ~rep+processing+sweetener+tiempo+Temp, data = flavanonesNorm.t, method = "xgbTree",
                    trControl=ctrl.reg,
                    tuneGrid = tune_grid,
                    tuneLength = 20)

xgbFit.N7O <- train(Naringenin..7.O.rutinoside ~rep+processing+sweetener+tiempo+Temp, data = flavanonesNorm.t, method = "xgbTree",
                    trControl=ctrl.reg,
                    tuneGrid = tune_grid,
                    tuneLength = 20)

xgbFit.H7O <- train(Hespertin..7.O.rutinoside ~rep+processing+sweetener+tiempo+Temp, data = flavanonesNorm.t, method = "xgbTree",
                    trControl=ctrl.reg,
                    tuneGrid = tune_grid,
                    tuneLength = 20)


# testing ----



flavPred <- flavanonesNorm.t
flavPred$Temp<-factor(substr(flavPred$X,0,3))
levels(flavPred$Temp)<-c(20,4,20,4,20,4)

flavPred$O.trygycosil.naringenin <- predict(xgbFit.3gN,flavPred)
flavPred$Eriodyctiol.7.O.rutinoside <- predict(xgbFit.E7O,flavPred)
flavPred$Naringenin..7.O.rutinoside <- predict(xgbFit.N7O,flavPred)
flavPred$Hespertin..7.O.rutinoside <- predict(xgbFit.H7O,flavPred)

flavPred.m1 <- flavanonesNorm %>% filter (tiempo == "-1")
flavPred.m1$Temp<-factor(substr(flavPred.m1$X,0,3))
levels(flavPred.m1$Temp)<-c(20,4,20,4,20,4)

flavaPredDef <- rbind(flavPred.m1, flavPred)

flavanonesNorm.l <- melt(flavanonesNorm)
flavanonesNorm.l$tiempo <- as.numeric(flavanonesNorm.l$tiempo)


flavPred.l <- melt(flavaPredDef)
flavPred.l$tiempo <- as.numeric(flavPred.l$tiempo) 

ggplot(flavPred.l, aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+geom_smooth(method="loess", se=F)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")

flavPred.l$processing <- factor(flavPred.l$processing)
levels(flavPred.l$processing) <- c("high pressures 1", "high pressures 2", "thermal")

flavanonesNorm.l$processing <- factor(flavanonesNorm.l$processing)
levels(flavanonesNorm.l$processing) <- c("high pressures 1", "high pressures 2", "thermal")

ggplot(flavPred.l %>% filter(Temp == "20" & variable == "Eriodyctiol.7.O.rutinoside"), aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+geom_smooth(method="loess", se=F)+
  geom_point(data = flavanonesNorm.l %>% filter(Temp == "4" & variable == "Eriodyctiol.7.O.rutinoside"), aes(x= tiempo, y = value), shape=4)+
  geom_smooth(method="loess", se=F)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")+
  ggtitle("Predictions on flavanones degradation at 20ºC | Crosses = experimental values at 4ºC")+
  theme(title = element_text(size=16, face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size=14))


ggplot(flavPred.l %>% filter(Temp == "4"), aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+geom_smooth(method="loess", se=F)+
  geom_point(data = flavanonesNorm.l %>% filter(Temp == "20"), aes(x= tiempo, y = value), shape=4)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")+
  ggtitle("Predictions on flavanones degradation at 4ºC")


# for t in unique(flavanones$tiempo):
# 
#   if X == ST1_rep1 &
#   predict(xgbFit.3gN, newdata =  )







