require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(forecast)

path <- "C:/Users/ducvi/OneDrive/Documents/Mathématiques/Séries temporelles/Projet"
setwd(path) #definit l'espace de travail (working directory ou "wd")

datafile <- "valeurs_mensuelles.csv" #definit le fichier de donnees

data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame
xm.source <- zoo(data[[2]]) #convertit le deuxieme element de data en serie temporelle de type "zoo"
T <- length(xm.source)
xm <- xm.source[4:T] #supprime les 3 premieres valeurs qui sont des titres et pas des donnees ; on obtient une serie de janvier 2023 à janvier 1990
xm <- rev(xm) #passage d'un ordre antechronologique a un ordre chronologique ; on obtient une serie de janvier 1990 à janvier 2023
xm <- xm[(1+24):(T-24-3)] #supprime les 2 premières années i.e. les 24 valeurs de 4 à 27 ou encore de janvier 1990 à décembre 1991 et les 2 dernieres années i.e. les 24 valeurs de 400 à 377 ou encore de janvier 2023 à février 2021
#serie de 28 à 376 i.e. de janvier 1992 à janvier 2021 

plot(xm, ylim=c(56,124), xaxt="n") #chronogramme, ylim nécessaire sinon error, supprime les labels et graduations de l'axe des abcisses
axis(side=1, ylim=c(56,124), at=seq(0,376,12)) #crée les labels et graduations de l'axe des abcisses

xm <- as.numeric(xm) #data[[2]] est de classe character ; on convertit en classe numeric pour acf
acf(xm)

dates <- as.yearmon(seq(from=1992, to=2021, by=1/12)) 
summary(lm(xm ~ dates))

require(fUnitRoots) #tests de racine unitaire plus modulables
adf <- adfTest(xm, lag=0, type="ct") #
adf

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

series <- xm; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(xm,24,adftype="ct")

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

adf

dxm <- diff(xm,1)
summary(lm(dxm ~ dates[-1]))
adf <- adfTest_valid(dxm,24,"nc")
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf

plot(dxm,type="l",xaxt="n")
axis(side=1, ylim=c(56,124), at=seq(0,376,12)) 

par(mfrow=c(1,2))
pacf(dxm,23);acf(dxm,23) #on regarde jusqu'a deux ans de retard


pmax=2;qmax=0

y <- dxm - mean(dxm) #
arima200 <- arima(y,c(2,0,0)) 
arima200
Box.test(arima200$residuals, lag=24, type="Ljung-Box", fitdf=2) 
Qtests(arima200$residuals, 24, 2) 
round(Qtests(arima200$residuals,24,fitdf=2),3)

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}
signif(arima200) 


mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #matrice vide `a remplir
rownames(mat) <- paste0("p=",0:pmax) #renomme les lignes
colnames(mat) <- paste0("q=",0:qmax) #renomme les colonnes
AICs <- mat #matrice des AIC non remplie
BICs <- mat #matrice des BIC non remplie
pqs <- expand.grid(0:pmax,0:qmax) #toutes les combinaisons possibles de p et q
for (row in 1:dim(pqs)[1]){ #boucle pour chaque (p,q)
  p <- pqs[row,1] #r´ecup`ere p
  q <- pqs[row,2] #r´ecup`ere q
  estim <- try(arima(y,c(p,0,q),include.mean = F)) #tente d’estimer l’ARIMA
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigne l’AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigne le BIC
}
AICs
BICs

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}

estim <- arima(y,c(1,0,0)); arimafit(estim)

Mod(polyroot(sort(arima200$coef[c('intercept','ar1', 'ar2')])))

models <- c("arima200")
preds <- zoo(matrix(NA,ncol=1,nrow=2),order.by=tail(index(xm.source),2))
colnames(preds) <- models
desaisonp <- preds #series vierges pour les pr´ediction de desaison
xmp <- preds #series vierges pour les pr´ediction de xm
##prediction de desaison et xm par chaque mod`ele
for (m in models){
  pred1 <- mean(dxm) + zoo(predict(get(m),2)$pred, order.by=tail(index(xm.source),2))
  pred2 <- as.numeric(tail(xm,12))[1:2] + pred1
  desaisonp[,m] <- pred1
  xmp[,m] <- pred2
}
xmp


arima_xm <- arima(xm, c(2, 1, 0))
prev <- forecast(arima_xm, h=2)
par(mfrow=c(1,1))
plot(prev, xlim=c(200, 350))

