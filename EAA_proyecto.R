##Required Packages
library(ggplot2)
library(dplyr)
library(purrr)
library(coda)
library(rjags)
library(corrplot)
library(tibble)
library(matlib)
#data load
datosMerk <- read.csv("femcare.csv")
str(datosMerk)
#fecha means date in Spanish
colnames(datosMerk)[1] <- "fecha"
datosMerk <- datosMerk %>% mutate(fecha = as.Date(fecha, "%d/%m/%Y")) %>% arrange(fecha) %>% 
  mutate(logit_B = logit(share_B))
summary(datosMerk)
##Separar por pais
datosMerk_CR <- datosMerk %>% filter(pais == "CR")
datosMerk_EL <- datosMerk %>% filter(pais == "EL")
datosMerk_GT <- datosMerk %>% filter(pais == "GT")
datosMerk_HD <- datosMerk %>% filter(pais == "HD")
datosMerk_NIC <- datosMerk %>% filter(pais == "NIC")
datosMerk_PTY <- datosMerk %>% filter(pais == "PTY")
nomPais <- c("Total", "CR", "EL", "GT", "HD", "NIC", "PTY")
names(paises) <- nomPais
paises <- list(datosMerk, datosMerk_CR, datosMerk_EL, datosMerk_GT, datosMerk_HD, datosMerk_NIC, datosMerk_PTY)
for(i in 1:length(paises)){
  paises[[i]] <- paises[[i]][, c(1, 2, 8, 3, 4, 5, 6, 7, 9, 10, 16, 11, 12, 13,14,15, 17)]
}
for(i in 1:length(paises)){
  #par(mfrow= c(1,2))
  corrplot(cor(paises[[i]][,3:16]), tl.col="black", type = "upper")
  title(nomPais[i])
  corrplot(cor(paises[[i]][,11:16]), tl.col="black", type = "upper")
  title(nomPais[i])
  #par(mfrow= c(1,1))
}
Ajuste <- paises
Test <- Ajuste
names(Ajuste) <- nomPais
names(Test) <- nomPais

for(i in 1:length(paises)){
  Ajuste[[i]] <- paises[[i]] %>% filter(fecha <= as.Date("31/10/2019", "%d/%m/%Y"))
  Test[[i]] <- paises[[i]] %>% filter(fecha > as.Date("31/10/2019", "%d/%m/%Y"))
  
}
str(Ajuste[[5]])
hist(datosMerk$share_A, breaks = 15)
plot(x = c(1:length(datosMerk_CR$share_B)),y= datosMerk_CR$share_B, type = "l")


datMat <- list()
datMat_Test <- list()
#Numero de coeficientes
p <- 15

##Standarize covariates (different scales)
for(j in 1: 6){
  datMat[[j]] <- matrix(rep(1, p*length(Ajuste[[j+1]]$share_A)), ncol = p, nrow = length(Ajuste[[j+1]]$share_A))
  for(i in 2:15){
    datMat[[j]][, i] <- (Ajuste[[j+1]][,i+2]-mean(Ajuste[[j+1]][,i+2]))/sd(Ajuste[[j+1]][,i+2])
  }
  colnames(datMat[[j]] ) <- c("int", "A1", "A2", "A3", "A4", "A5", "A6", "A7","B1", "B2", "B3", "B4", "B5", "B6", "B7")
}
for(j in 1: 6){
  datMat_Test[[j]] <- matrix(rep(1, p*length(Test[[j+1]]$share_A)), ncol = p, nrow = length(Test[[j+1]]$share_A))
  for(i in 2:15){
    if(sd(Test[[j+1]][,i+2])>0){
    datMat_Test[[j]][, i] <- (Test[[j+1]][,i+2]-mean(Test[[j+1]][,i+2]))/sd(Test[[j+1]][,i+2])
    } else {
      datMat_Test[[j]][, i] <- Test[[j+1]][,i+2]
    }
  }
  colnames(datMat_Test[[j]] ) <- c("int", "A1", "A2", "A3", "A4", "A5", "A6", "A7","B1", "B2", "B3", "B4", "B5", "B6", "B7")
}

##
##Model Fitting

#Whic covariates are included in the model
coefs <- c(1, 4, 5, 12)
p2 <- length(coefs)
##Country to adjust (CR=1, EL =2, GT=3, HD=4, NIC=5, PTY=6)
pais <- 1
##Prior Selection (normal conjugates, double exponential conjugate, flat ), RJAGS
model_string <- textConnection("model{
                      #Verosimilitud
                      for(i in 1:n){
                        Y[i] ~ dnorm(inprod(X[i,], beta[]), taue)
                               }
                       #Priors
                      beta[1] ~ dnorm(0,0.01)
                       for ( j in 2:p){
                          beta[j] ~ dnorm(0, taue*taub )
                       }
                      taue ~ dgamma( 1, 0.1)
                      taub ~ dgamma( 1, 0.1)
                      sigma <- 1/sqrt(taue)
                               }")
##Data for RJAGS simulation
data_Jags <- list(Y = logit(Ajuste[[pais+1]]$share_A), X = datMat[[pais]][, coefs],
                  n = length(Ajuste[[pais+1]]$share_A), p = p2)
##Initial Values
inits <- list(beta = rnorm(p2), taue = rgamma(1,1) ,taub = rgamma(1,1) 
              )
##Fiited model
modelo_1 <- jags.model(model_string, data = data_Jags, inits = inits, n.chains = 2)
update(modelo_1, 10000, progres.bar="none")
samples <- coda.samples(modelo_1, variable.names = c("beta", "sigma"),
                        n.iter = 10000, progess.bar="none")
summary(samples)[[1]][1:p2,1]
min(abs(summary(samples)[[1]][1:p2,1]))
#plot(samples)
#View(samples)
##DIC deviance criteia
DIC <- dic.samples(modelo_1, n.iter = 2500, progress.bar="none")
DIC

#Point Estimate
#Square errors
#Posterior predictive simulation
logit_Yp <- matrix(rep(0,length(Test[[pais+1]]$share_A)*10000), nrow=10000, ncol =6)
Yp <- logit_Yp
for(s in 1:10000){
  logit_Yp[s, ] <- datMat_Test[[pais]][, coefs] %*% samples[[1]][s, 1:p2] + rnorm(nrow(datMat_Test[[pais]]), mean = 0, 
                                                                        sd = samples[[2]][s, p2+1])
  Yp[s, ] <- exp(logit_Yp[s, ])/(1+exp(logit_Yp[s, ]))
  
}
er <- c(1:6)
for(j in 1:6){
  print(quantile(Yp[,j], probs = c(0.0275, 0.50, 0.975)))
  print(paste("Valor Real: ", Test[[pais+1]][j,3]))
  er[j] <- mean(Yp[,j]) - Test[[pais+1]][j,3]
  print(paste("Error: ", er[j]))
}
print(paste("MSE: ", mean(er^2)))


##Jeffreys Prior
beta_LS <- solve(t(datMat[[pais]][, coefs]) %*% datMat[[pais]][, coefs]) %*% 
  t(datMat[[pais]][, coefs])%*% logit(Ajuste[[pais+1]]$share_A)
sigma2 <- mean((logit(Ajuste[[pais+1]]$share_A) - datMat[[pais]][, coefs] %*% beta_LS)^2)
beta_cov <- sigma2 *  solve(t(datMat[[pais]][, coefs]) %*% datMat[[pais]][, coefs])
beta_sd <- sqrt(diag(beta_cov))
df <- length( logit(Ajuste[[pais+1]]$share_A))
beta_025 <- beta_LS + beta_sd*qt(0.025, df = df)
beta_975 <- beta_LS + beta_sd*qt(0.975, df = df)
for(i in 1:p2){
  print(c(beta_025[i], beta_LS[i], beta_975[i]))
}


##Zellner prior model (Correlated covariates)
modeloZellner <- gibbsZellner(10000, betaLS = beta_LS, XtX = t(datMat[[pais]][, coefs]) %*% datMat[[pais]][, coefs], 
             c=1, a=1, b=1, p=p2, n=length(Ajuste[[pais+1]]$share_A))
modeloZellner <- cbind(modeloZellner, 1/sqrt(modeloZellner[, p2+1]))
summary(modeloZellner)

logit_YpZ <- matrix(rep(0,length(Test[[pais+1]]$share_A)*10000), nrow=10000, ncol =6)
YpZ <- logit_Yp
for(s in 1:10000){
  logit_YpZ[s, ] <- datMat_Test[[pais]][, coefs] %*% modeloZellner[s, 1:p2] + rnorm(nrow(datMat_Test[[pais]]), mean = 0, 
                                                                                  sd =  modeloZellner[s, p2+2])
  YpZ[s, ] <- exp(logit_YpZ[s, ])/(1+exp(logit_YpZ[s, ]))
  
}

erZ <- c(1:6)
for(j in 1:6){
  
  erZ[j] <- mean(YpZ[,j]) - Test[[pais+1]][j,3]

}
print(paste("MSE: ", mean(erZ^2), "Zellner"))

if(mean(erZ^2)/mean(er^2) > 1){
  print("a-priori without covariate correlation has a better performance")
} else {
  print("Zellner's model has a better performance")
}
summary(lm(logit(share_A) ~ ventas_unidades_A + ppu_A + ventas_unidades_B, data = Ajuste[[pais+1]]))
colnames(Ajuste[[pais+1]])[3]
##Final fitted models
## CR: logit(share_A) ~ b0 + b1*A3 + b2*A4 + b3*B4, prior normal flat; bj ~ N(0, sigma2*100^2) 
## EL: logit(sahre_A) ~ b0 + b1*A3 + b2*A4 +b3*B4 + b4*B7
## GT:logit(share_A) ~ b0 + b1*A3 + b2*B7; prior doble expoencial
## HD: logit(share_A) ~ b0 + b1*A3 + b2*B4 + b3*B7, prior normales con precision comun entre coefs
## NIC: logit(share_A) ~ b0 + b1*A3 + b2*A4 + b3*B5 + b4*B6 + b5*B7 ; prior normales con precision comun entre coefs
## PTY: logit(share_A) ~ b0 + b1*A3 + b2*B7; prior doble exponencial con precision comun entre coefs



##Auxilary functions
logit <- function(x){
  return(log(x/(1-x)))
}

#Multivariate normal sample simulation
rmvnorm <- function(n, mu, Sigma){
  Y <- matrix( rep(0, n * length(mu)), ncol = length(mu), nrow = n)
  colN <- rep("s", length(mu))
  for(i in 1:n){
     z <- rnorm(length(mu), mean = 0, sd = 1)
     U <- chol(Sigma)
     Y[i, ] <- t(U) %*% z + mu
  }
  for(j in 1:length(mu)){
  colN[j] <- paste("Y", j) 
  }
  colnames(Y) <- colN
  return(Y)
}
rmvnorm(10, c(1, 1, -1), matrix(c(2,1, -1,1,2, -.5, -1, -.5, 1), nrow=3, ncol=3))

#Gibbs sampler to simulate Zellner's model posterior distribution
gibbsZellner <- function(S, betaLS, XtX, c, a, b, p, n){
  XtXin <- solve(XtX)
  chain <- array(dim = c(2*S+1, p+1))
  chain[1, ] <- c(rmvnorm(1, mu=rep(0,p), Sigma=XtXin), rgamma(1, a,b))
  for(s in 1:(2*S)){
    chain[s+1, 1:p] <- sqrt(c/chain[s, p+1])*rmvnorm(1, mu = c*betaLS, Sigma = XtXin)
    by <- 1/c* (t(chain[s+1, 1:p] - betaLS) %*% XtX %*% (chain[s+1, 1:p] - betaLS))
    chain[s+1, p+1] <- rgamma(1, shape = a+(n+p)/2, rate = by)
  }
  return(chain[-(1:S+1), ])
  
}



