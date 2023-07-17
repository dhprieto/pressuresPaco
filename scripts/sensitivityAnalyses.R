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

# funciones: ----

aadha.ldef<-function(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0=0){
  yhat<-rep(NA,length(Time))
  # define the components of the matrix as a data frame
  iPred<-data.frame(Time,Species,k1,k2,k3,aa0,dhaa0,dkg0)
  
  for (i in 1:length(Time)){
    # build the coefficients matrix from iPred for every row. 
    # Is Time selected for loop arbitrary?
    M <- matrix(c(-iPred$k1[i],iPred$k2[i],0,iPred$k1[i],-(iPred$k2[i]+iPred$k3[i]),0,0,iPred$k3[i],0), 3, 3, byrow=TRUE)
    # define initial concentrations of compounds and runs solver function
    C0<- c(iPred$aa0[i],iPred$dhaa0[i],iPred$dkg0[i])
    yt <- lin.diff.eqn(A=M,initial=C0,t=Time[i])
    # add the value to the species evaluated for every step
    yhat[i]<-switch(as.character(iPred$Species[i]),
                    Ascorbic=yt[1],
                    Dehydroascorbic=yt[2])
  }
  return(yhat)
}

ggplot.sobolmartinez.custom <- function(data, mapping = aes(), ylim = c(0, 1), 
                                        y_col = NULL, y_dim3 = NULL, ..., environment = parent.frame()) {
  x <- data
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    if(is(x$y,"numeric")){
      nodeggplot(listx = list(x$S[row.names(x$S)%in% c("Temp","tiempo", "sweetener", "processing", "compound"),],
                              x$T[row.names(x$T)%in% c("Temp","tiempo", "sweetener", "processing", "compound"),]), 
                 xname = c("Main effect","Total effect"), ylim = ylim, pch = pch)
    } else if(is(x$y,"matrix") | is(x$y,"array")){
      if(is.null(y_col)) y_col <- 1
      if(is(x$y,"matrix") && !is.null(y_dim3)){
        y_dim3 <- NULL
        warning("Argument \"y_dim3\" is ignored since the model output is ",
                "a matrix")
      }
      if(is(x$y,"array") && !is(x$y,"matrix") && is.null(y_dim3)) y_dim3 <- 1
      nodeggplot(listx = list(x$S,x$T), xname = c("Main effect","Total effect"), ylim = ylim, pch = pch, y_col = y_col, y_dim3 = y_dim3)
    }
  }
}


nodeggplot <- function(listx, xname, xlim = NULL, ylim = NULL, labels = TRUE, title = NULL, 
                       col = par("col"), pch = 21, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  
  ngraphs <- length(listx)
  x <- unlist(listx)
  n <- nrow(listx[[1]])
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (n<=10){
    angle <- 0
    hjust <- 0
  }
  if (n>10 & n <=20){
    angle <- 45
    hjust <- 1
  }
  if (n>20){
    angle <- 90
    hjust <- 1
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (is.null(title)){
    title <- title
  }
  
  #  if (class(labels) == "logical") {
  if (inherits(labels, "logical")){
    if (labels) {
      l <- rownames(listx[[1]])
    } else {
      l <- NULL
    }
    #  } else if (class(labels) == "character") {
  } else if (inherits(labels, "character")){
    l <- labels
  }
  
  # bias
  
  d <- NULL
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if ("bias" %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]] - x[["bias"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, "original", y_col] - x[, "bias", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
      }
    } else {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, y_col, y_dim3]
      }
    }
    d <- rbind(d,data.frame(x=at,y=xx,id=xname[i]))
  }
  
  # confidence intervals
  d2 <- NULL
  n2 <- 0
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[["min. c.i."]]
        max_ci <- x[["max. c.i."]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col]
        max_ci <- x[, "max. c.i.", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col, y_dim3]
        max_ci <- x[, "max. c.i.", y_col, y_dim3]
      }
      n2 <- n2 +1
    }else{
      min_ci <- rep(NA,n)
      max_ci <- rep(NA,n)
    }
    d2 <- rbind(d2,data.frame(min_ci=min_ci,max_ci=max_ci))
  }
  
  d <- cbind(d,d2)
  
  if (ngraphs>1){
    pd <- position_dodge(0.3)
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, position = pd) + 
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }else{
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }
  }else{
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, colour=col) + 
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15))
    }else{
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15)) 
    }
  }
  return(g)
}




# Anthocyanins ----

# reading data
anthocyanins <- read.csv("data/anthocyanins_longFormat.csv")
coef_ant.fm0 <- read.csv("data/coef_ant.fm0.csv")

anthocyanins <- anthocyanins %>% mutate_at(c("processing", "sweetener", "compound"), factor)

# fit model ant

ant.fm3<-gnls(concentration~Cinf+(C0-Cinf)*exp(-exp(lk-Ea/8.314e-3*(1/(Temp+273)-1/(16+273)))*tiempo),
              data=anthocyanins,
              param=list(C0~compound+compound:sweetener+compound:processing,
                         Cinf~compound+compound:sweetener+compound:processing,
                         lk~compound+compound:sweetener+compound:processing,
                         Ea~compound),
              start=c(C0=c(coef_ant.fm0[1:6,2],rep(0.001,18)),
                      Cinf=c(coef_ant.fm0[7:12,2],rep(0.001,18)),
                      lk=c(coef_ant.fm0[13:18,2],rep(0.001,18)),
                      Ea=c(coef_ant.fm0[19:24,2],rep(0.001,0)))
)


# generating variables for Sobol-Martinez

tiempo.Ant <- rep(unique(anthocyanins$tiempo),125)
names.compounds.ant <- c("Delphinidin.3.O.sambubioside.5.O.glucoside","Delphinidin.3.5.O.diglucoside",
                     "Cyanidin.3.O.sambubioside.5.O.glucoside...Cyanidin.3.5.O.diglucoside",
                     "Delphinidin.3.O.sambubioside","Delphinidin.3.O.glucoside",
                     "Cyanidin.3.O.sambubioside..Cyanidin.3.O.glucoside")




n <- 1000

X1.ant <- data.frame(matrix(c(runif(4*n), rep(c(4,20), 500),tiempo.Ant,
                          c(rep(names.compounds.ant,166),names.compounds.ant[1:4]), 
                          rep(c("SU","ST"), 500),c(rep(c("1", "2", "P"), 333), "1")), nrow=n))


colnames(X1.ant) <- c("C0", "Cinf", "lk", "Ea", "Temp","tiempo","compound", "sweetener", "processing")

X1.ant <- X1.ant %>% mutate_at(c("C0", "Cinf", "lk", "Ea", "Temp","tiempo"), as.numeric)
X1.ant <- X1.ant %>% mutate_at(c("compound", "sweetener", "processing"), factor)

X2.ant <- data.frame(matrix(c(runif(4*n), sample(rep(c(4,20), 500)),sample(tiempo.Ant),
                          sample(c(rep(names.compounds.ant,166),names.compounds.ant[1:4])), 
                          sample(rep(c("SU","ST"), 500)),sample(c(rep(c("1", "2", "P"), 333), "1"))),nrow=n))

colnames(X2.ant) <- c("C0", "Cinf", "lk", "Ea", "Temp","tiempo","compound", "sweetener", "processing")
X2.ant <- X2.ant %>% mutate_at(c("C0", "Cinf", "lk", "Ea", "Temp","tiempo"), as.numeric)
X2.ant <- X2.ant %>% mutate_at(c("compound", "sweetener", "processing"), factor)

# running Sobol-Martinez

ant.fm3.sobol.martinez <- sobolmartinez(model = ant.fm3, X1.ant,X2.ant)

ggplot.sobolmartinez.custom(ant.fm3.sobol.martinez)+ggtitle("Sobol sensitivity analysi for Anthocyanin Model")

saveRDS(ant.fm3.sobol.martinez, file ="results/AntSobol.RDS")

capture.output(ant.fm3.sobol.martinez, file="results/AntindicesSobol.csv")


# Vitamin C ----

# reading data

VitC.l <- read.csv("data/vitC_longFormat.csv")

VitC.l$tiempo <- VitC.l$tiempo+1.0001 
VitC.l<- select(VitC.l, -Species)
VitC.l<- VitC.l %>% mutate_at(c("processing", "sweetener", "compound"), factor)

# fit model
VitC.fmF<-gnls(Concentracion~aadha.ldef(Time=tiempo,Species = compound,
                                        k1=exp(lk1-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k2=exp(lk2-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        k3=exp(lk3-Ea/8.314e-3*(1/(Temp+273)-1/(16+273))),
                                        aa0=exp(laa0),dhaa0=exp(ldhaa0)),
               data=VitC.l,
               weights=varIdent(~1|compound),
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

# generating vars for Sobol-Martinez
tiempo.VitC <- rep(unique(VitC.l$tiempo),125)
names.compounds.VitC <- rep(unique(VitC.l$compound), 500)

n <- 1000
X1.VitC <- data.frame(matrix(c(runif(6*n), rep(c(4,20), 500),tiempo.VitC,
                          names.compounds.VitC, 
                          rep(c("SU","ST"), 500),
                          c(rep(c("1", "2", "P"), 333), "1")),nrow=n))


colnames(X1.VitC) <- c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo","compound", "sweetener", "processing")

X1.VitC <- X1.VitC %>% mutate_at(c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo"), as.numeric)

X1.VitC$compound[X1.VitC$compound == "1"] <-"Ascorbic" 
X1.VitC$compound[X1.VitC$compound == "2"] <-"Dehydroascorbic"

X2.VitC <- data.frame(matrix(c(runif(6*n), sample(rep(c(4,20), 500)),sample(tiempo.VitC),
                          sample(names.compounds.VitC), 
                          sample(rep(c("SU","ST"), 500)),sample(c(rep(c("1", "2", "P"), 333), "1"))),nrow=n))

colnames(X2.VitC) <- c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo","compound", "sweetener", "processing")

X2.VitC <- X2.VitC %>% mutate_at(c("lk1", "lk2", "lk3", "Ea","laa0", "ldhaa0", "Temp","tiempo"), as.numeric)

X2.VitC$compound[X2.VitC$compound == "1"] <-"Ascorbic" 
X2.VitC$compound[X2.VitC$compound == "2"] <-"Dehydroascorbic"

# running Sobol-Martinez

VitC.fmF.sobol.martinez <- sobolmartinez(model = VitC.fmF, X1.VitC,X2.VitC)

ggplot.sobolmartinez.custom(VitC.fmF.sobol.martinez)+ggtitle("Sobol sensitivity analysi for Vitamin C Model")

saveRDS(VitC.fmF.sobol.martinez, file ="results/VitCSobol.RDS")

capture.output(VitC.fmF.sobol.martinez, file="results/VitCindicesSobol.csv")

# test ----

# coef_C0=c(5.1476533,-0.6894156, -2.9251113, -3.6165833, -1.2016035, -3.7776822)
# coef_Cinf=c(1.030496, 0.1705629, -0.6456348, -0.9549181, -0.7362231, -0.8933354)
# coef_lk=c(-3.973978,0.6668254,0.1156777,0.1142587,0.3928976,0.3214167)
# coef_Ea=c(72.08446,-5.107041,-3.558397,-11.21674,-6.684252,-8.100521)
# 
# c0 <- c(rep(coef_C0, 166),5.1476533,-0.6894156, -2.9251113, -3.6165833)
# cinf <- c(rep(coef_Cinf, 166),1.030496, 0.1705629, -0.6456348, -0.9549181)
# lk <- c(rep(coef_lk, 166), -3.973978,0.6668254,0.1156777,0.1142587)
# Ea <- c(rep(coef_Ea, 166), 72.08446,-5.107041,-3.558397,-11.21674)
