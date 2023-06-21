# libraries ----

library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)

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

flavPred.l <- melt(flavaPredDef)


ggplot(flavPred.l, aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+geom_smooth(method="lm", se=F)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")


ggplot(flavPred.l %>% filter(Temp == "20"), aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+
  geom_point(data = flavanonesNorm.l %>% filter(Temp == "4"), aes(x= tiempo, y = value), shape=4)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")+
  ggtitle("Predictions on flavanones degradation at 20ºC | Crosses = experimental values at 4ºC")


ggplot(flavPredN.l %>% filter(Temp == "4"), aes(x = tiempo, y = value, col= processing)) +
  facet_wrap(variable~sweetener,scales="free")+
  geom_point()+geom_smooth(method="lm", se=F)+
  xlab("Storage Time [Days]")+ylab("Concentration [normalized value]")+
  ggtitle("Predictions on flavanones degradation at 4ºC")


# for t in unique(flavanones$tiempo):
# 
#   if X == ST1_rep1 &
#   predict(xgbFit.3gN, newdata =  )







