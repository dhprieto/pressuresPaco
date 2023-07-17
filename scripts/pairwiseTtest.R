library(dplyr)
library(rstatix)


presiones1.1 <- read.csv("data/presiones1.csv", sep = ";", dec = ",")
anthocyanins <- presiones %>% dplyr::select(c("processing", "sweetener","tiempo", "Delphinidin.3.O.sambubioside.5.O.glucoside","Delphinidin.3.5.O.diglucoside",
                                              "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",
                                              "Delphinidin.3.O.sambubioside","Delphinidin.3.O.glucoside",
                                              "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside"))



anthocyanins <- anthocyanins %>% rename("Time" = "tiempo", "Sweetener" = "sweetener") 
anthocyanins$Time <- anthocyanins$Time + 1.001

pairwise_t_test(data = anthocyanins, formula = Delphinidin.3.O.glucoside ~Time,
                paired = F, p.adjust.method = "bonferroni") 


pwResult <- list()

pwResult <- pairwiseTTest(anthocyanins, varsRemoved = c("Time","Delphinidin.3.O.glucoside", 
                                            "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside"))

write.csv(pwResult, file = "pairwiseAnt.csv")

pairwiseTTest <- function(tabla, varsRemoved){

  for (i in colnames(tabla)){
    
    if (is.numeric(tabla[,i]) & !(i %in% varsRemoved)){
      
      message(paste("Analized variable: ", i))
      
      message(paste("Time comparisons ", i)) 
      
      timeComp <- pairwise_t_test(data = tabla , formula = as.formula(paste(sym(i),"~ Time")),
                            paired = F, p.adjust.method = "bonferroni") 
      
      tablaGr <- group_by(tabla, Sweetener, processing)
      
      message(paste("Time-Sweetener-Processing comparisons", i))
      
      timeSwetProcs <- pairwise_t_test(data = tablaGr, formula = as.formula(paste(sym(i),"~ Time")),
                            paired = F, p.adjust.method = "bonferroni") 
      
      
      tablaGr <- group_by(tabla, Sweetener)
      
      message(paste("Time-Sweetener comparisons", i))
      
      timeSweet <- pairwise_t_test(data = tablaGr , formula = as.formula(paste(sym(i),"~ Time")),
                            paired = F, p.adjust.method = "bonferroni") 
      
      
      tablaGr <- group_by(tabla, processing)
      
      message(paste("Time-Processing comparisons", i))
      
      timeProc <- pairwise_t_test(data = tablaGr , formula = as.formula(paste(sym(i),"~ Time")),
                            paired = F, p.adjust.method = "bonferroni") 
      
      pwResult <- append(pwResult, c(paste("Analized variable: ", i), paste("Time comparisons ", i), as.data.frame(timeComp), 
                                     paste("Time-Sweetener-Processing comparisons", i), as.data.frame(timeSwetProcs),
                                     paste("Time-Sweetener comparisons", i), as.data.frame(timeSweet), 
                                     paste("Time-Processing comparisons", i), as.data.frame(timeProc)))
      }
      
  }
return (pwResult)
}
