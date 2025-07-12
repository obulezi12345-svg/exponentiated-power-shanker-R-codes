
setwd("C:\\Users\\1012 G2\\Documents\\R files\\EPS")


require(survival)
require(AdequacyModel)
# Function to calculate information criteria (AIC, BIC, CAIC, HQIC)
criterios <- function(value, num_params, num_data_points) {
  # value - the value from optimization
  # num_params - the number of parameters
  # num_data_points - the number of data points
  
  # Calculate log-likelihood
  l <- 2 * value
  
  # Calculate AIC
  AIC <- l + 2 * num_params
  
  # Calculate BIC
  BIC <- l + num_params * log(num_data_points)
  
  # Calculate CAIC
  CAIC <- AIC + (2 * (num_params + 2) * (num_params + 3)) / (num_data_points - num_params - 3)
  
  # Calculate HQIC
  HQIC <- l + 2 * log(log(num_data_points)) * num_params
  
  # Combine results into a matrix
  result <- cbind(AIC, CAIC, BIC, HQIC)
  
  # Return results
  return(result)
}



# Function to calculate confidence intervals and p-values of parameters
IC <- function(parametros, hessiana, n) {
  # Calculate variance-covariance matrix

  hessiana <- hessiana + diag(rep(1e-6, length(parametros)))
 
  # Calculate variance-covariance matrix
  Var <- solve(hessiana, tol = 1e-15)
  # Calculate standard errors of parameters
  SE <- sqrt(diag(Var))
  
  # Calculate t-values of parameters
  tvalue <- parametros / SE
  
  # Calculate p-values of parameters
  pvalue <- 2 * pt(abs(tvalue), n, lower.tail = FALSE)
  
  # Calculate confidence intervals of parameters
  LI <- parametros - qt(0.975, n) * SE
  LS <- parametros + qt(0.975, n) * SE
  
  # Combine results into a matrix
  resul <- cbind(parametros, SE, tvalue, pvalue, LI, LS)
  colnames(resul) <- c("Parameter", "SE", "t-value", "p-value", "Lower Bound", "Upper Bound")
  
  # Return results
  return(resul)
}



dados = read.table("data.txt", header = T)
t = dados$times
hist(t,border = "black",col = "red", xlab = "age ( in years)",main = "")
y=log(t)
censur = dados$censor





# Log-exponentiated power Shanker

EPS<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z <- (y - mu) / sigma
  E1 = exp(z)
  E2 = (exp(-(2*mu)/sigma) + 1)
  E3 = exp(-exp(z))
  E4 = exp((y-2*mu)/sigma)
  E5 = (exp(-(mu)/sigma)+exp(y/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.43775507,  2.41063213,  0.33183901, 0.11822759, 0.80031389, 0.81063213,  2.05183901, 0.01822759, 0.88031389)

#LEPS<-optim(valorin,EPS,method="BFGS",hessian=T);
#LEPS <- optim(valorin, EPS, method="Nelder-Mead", hessian=T, control=list(maxit=10000))
valorin <- valorin + runif(length(valorin), -0.1, 0.1)  # Small random perturbation
LEPS <- optim(valorin, EPS, method="BFGS", hessian=T)

LEPS
n=nrow(dados)
criterios(LEPS$value,length(LEPS$par),n)
parametros1=LEPS$par
hessiana1=LEPS$hessian
IC(parametros1,hessiana1,n)



# Log-exponentiated power chris-Jerry

EPCJ<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 2)
  E2 = (exp(-(mu)/sigma) + 2)
  E3 = exp(-exp(z))
  E4 = exp((y-2*mu)/sigma)
  E5 = (1 + exp(((2*y)-mu)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.08775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LEPCJ<-optim(valorin,EPCJ,method="BFGS",hessian=T);
LEPCJ
n=nrow(dados)
criterios(LEPCJ$value,length(LEPCJ$par),n)
parametros1=LEPCJ$par
hessiana1=LEPCJ$hessian
IC(parametros1,hessiana1,n)


# Log-exponentiated power Akash

EPA<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 2)
  E2 = (exp(-(2*mu)/sigma) + 2)
  E3 = exp(-exp(z))
  E4 = exp((y-3*mu)/sigma)
  E5 = (1 + exp((2*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.88775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LEPA<-optim(valorin,EPA,method="BFGS",hessian=T);
LEPA
n=nrow(dados)
criterios(LEPA$value,length(LEPA$par),n)
parametros1=LEPA$par
hessiana1=LEPA$hessian
IC(parametros1,hessiana1,n)

#summary(dados)





# Log-power Prakaamy 

PP<-function(par){
  c     = 1
  sigma = par[1]
  b0    = par[2]
  b1    = par[3]
  b2    = par[4]
  b3    = par[5]
  b4    = par[6]
  b5    = par[7]
  b6    = par[8]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  z <- (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 6)
  E2 = (exp(-(4*mu)/sigma) + 6)
  E3 = exp(-exp(z))
  E4 = exp((y-4*mu)/sigma)
  E5 = (exp(-mu/sigma) + exp((3*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.88775507,  0.01063213,  0.55183901, 2.21822759, 0.28031389, 2.81063213,  2.55183901, 0.01822759)
LPP<-optim(valorin,PP,method="BFGS",hessian=T);
LPP
n=nrow(dados)
criterios(LPP$value,length(LPP$par),n)
parametros2=LPP$par
hessiana2=LPP$hessian
IC(parametros2,hessiana2,n)

# Log-exponentiated power Lindley

EPL<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y - mu) / sigma
  E1 = exp(z)
  E2 = (exp(-(mu)/sigma) + 1)
  E3 = exp(-exp(z))
  E4 = exp((y-2*mu)/sigma)
  E5 = (1 + exp(y/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.08775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LEPL<-optim(valorin,EPL,method="BFGS",hessian=T);
LEPL
n=nrow(dados)
criterios(LEPL$value,length(LEPL$par),n)
parametros1=LEPL$par
hessiana1=LEPL$hessian
IC(parametros1,hessiana1,n)



# Log-power Rama 

PR<-function(par){
  c     = 1
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  z <- (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 6)
  E2 = (exp(-(3*mu)/sigma) + 6)
  E3 = exp(-exp(z))
  E4 = exp((y-4*mu)/sigma)
  E5 = (1 + exp((3*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.08775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LPR<-optim(valorin,PR,method="BFGS",hessian=T);
LPR
n=nrow(dados)
criterios(LPR$value,length(LPR$par),n)
parametros2=LPR$par
hessiana2=LPR$hessian
IC(parametros2,hessiana2,n)




# Log- power Lomax

PLO<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 2)
  E2 = (exp(-(2*mu)/sigma) + 2)
  E3 = exp(-(c*mu)/sigma)
  E4 = exp((y-c*mu)/sigma)
  E5 = (exp(-(mu)/sigma) + exp(y/sigma))
  g  = (E5^(-c-1))
  G  =  1 - E3*E5^(-c)
  pdf = c*E4*g/sigma
  sf = 1 - G
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.08775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LPLO<-optim(valorin,PLO,method="BFGS",hessian=T);
LPLO
n=nrow(dados)
criterios(LPLO$value,length(LPLO$par),n)
parametros1=LPLO$par
hessiana1=LPLO$hessian
IC(parametros1,hessiana1,n)


# Log-power Zeghdoudi


PZ<-function(par){
  c     = 1
  sigma = par[1]
  b0    = par[2]
  b1    = par[3]
  b2    = par[4]
  b3    = par[5]
  b4    = par[6]
  b5    = par[7]
  b6    = par[8]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z <- (y - mu) / sigma
  E1 = exp(z)*(exp(z) +2+exp(-(mu)/sigma))
  E2 = (exp(-(mu)/sigma)+2)
  E3 = exp(-exp(z))
  E4 = exp(((2*y)-(3*mu))/sigma)
  E5 = (1+ exp(y/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.88775507,  0.01063213,  0.55183901, 2.21822759, 0.28031389, 2.81063213,  2.55183901, 0.01822759)
LPZ<-optim(valorin,PZ,method="BFGS",hessian=T);
LPZ
n=nrow(dados)
criterios(LPZ$value,length(LPZ$par),n)
parametros2=LPZ$par
hessiana2=LPZ$hessian
IC(parametros2,hessiana2,n)


# Log-power Suja


PSU<-function(par){
  c     = 1
  sigma = par[1]
  b0    = par[2]
  b1    = par[3]
  b2    = par[4]
  b3    = par[5]
  b4    = par[6]
  b5    = par[7]
  b6    = par[8]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z <- (y - mu) / sigma
  E1 = exp(4*z)+4*exp(3*z)+12*exp(2*z)+ 24*exp(z)
  E2 = (exp(-(4*mu)/sigma)+24)
  E3 = exp(-exp(z))
  E4 = exp((y-(5*mu))/sigma)
  E5 = (1+ exp((4*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.88775507,  0.01063213,  0.55183901, 2.21822759, 0.28031389, 2.81063213,  2.55183901, 0.01822759)
LPSU<-optim(valorin,PSU,method="BFGS",hessian=T);
LPSU
n=nrow(dados)
criterios(LPSU$value,length(LPSU$par),n)
parametros2=LPSU$par
hessiana2=LPSU$hessian
IC(parametros2,hessiana2,n)

# Log-exponentiated power Ishita

EPI<-function(par){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  z <- (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 2)
  E2 = (exp(-(3*mu)/sigma) + 2)
  E3 = exp(-exp(z))
  E4 = exp((y-3*mu)/sigma)
  E5 = (exp(-mu/sigma) + exp((2*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.43775507,  2.41063213,  0.33183901, 0.11822759, 0.80031389, 0.81063213,  2.05183901, 0.01822759, 0.88031389)

LEPI<-optim(valorin,EPI,method="BFGS",hessian=T);
LEPI
n=nrow(dados)
criterios(LEPI$value,length(LEPI$par),n)
parametros1=LEPI$par
hessiana1=LEPI$hessian
IC(parametros1,hessiana1,n)

# Log-power Ishita 

PI<-function(par){
  c     = 1
  sigma = par[1]
  b0    = par[2]
  b1    = par[3]
  b2    = par[4]
  b3    = par[5]
  b4    = par[6]
  b5    = par[7]
  b6    = par[8]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  z <- (y - mu) / sigma
  E1 = exp(z)*(exp(z) + 2)
  E2 = (exp(-(3*mu)/sigma) + 2)
  E3 = exp(-exp(z))
  E4 = exp((y-3*mu)/sigma)
  E5 = (exp(-mu/sigma) + exp((2*y)/sigma))
  g  = (E4*E5*E3)/(sigma*E2)
  G  =  1 - (1 + E1/E2)*E3
  pdf = c*g*G^(c-1)
  sf = 1 - G^c
  lv<-censur*(log(pdf))+(1-censur)*(log(sf))
  sum(-lv)
}
valorin<-c(0.88775507,  0.01063213,  0.55183901, 2.21822759, 0.28031389, 2.81063213,  2.55183901, 0.01822759)
LPI<-optim(valorin,PI,method="BFGS",hessian=T);
LPI
n=nrow(dados)
criterios(LPI$value,length(LPI$par),n)
parametros2=LPI$par
hessiana2=LPI$hessian
IC(parametros2,hessiana2,n)

# Log-exponentiated WEibull

WE<-function(par){
  a     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  f = (a/sigma)*exp(z)*exp(-exp(z))*(1 - exp(-exp(z)))^(a-1)
  h = 1-(1-exp(-exp(z)))^(a)
  lv = censur*(log(f))+(1-censur)*(log(h))
  sum(-lv)
}
valorin<-c(0.08775507,  0.81063213,  2.55183901, 0.01822759, -0.28031389, 0.81063213,  2.55183901, 0.01822759, -0.28031389)
LWE<-optim(valorin,WE,method="BFGS",hessian=T);
LWE
n=nrow(dados)
criterios(LWE$value,length(LWE$par),n)
parametros4=LWE$par
hessiana4=LWE$hessian
IC(parametros4,hessiana4,n)

# quantile residuals



FEPS<-function(par,y){
  c     = par[1]
  sigma = par[2]
  b0    = par[3]
  b1    = par[4]
  b2    = par[5]
  b3    = par[6]
  b4    = par[7]
  b5    = par[8]
  b6    = par[9]
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$age + b2*dados$diab+b3*dados$heart+b4*dados$asthma+b5*dados$neuro+b6*dados$obesity
  z = (y-mu)/(sigma)
  v = (sigma*z -mu)/sigma
  z1= ((c*exp(v))/(sigma*(exp(-(2*mu)/sigma) + 1)))
  z2 = exp(-((sigma^2)*z+sigma*mu+mu)/sigma)+exp(-((sigma^2)*z+sigma*mu+sigma*z+mu)/sigma)
  z3 = exp(-exp(z))
  z4 = ((exp(z))/(exp(-(2*mu)/sigma) + 1))
  Ft = (1 - (1 + z4)*z3)^(c)
  return(Ft)
}

# Index plot

plot(qnorm(FEPS(LEPS$par,y)),pch=16,ylab = "Quantile residuals",xlab = "Index")
abline(h = c(-3,3), lty=2)

# Normal probability plot

q = rnorm(221)
qqplot(q,qnorm(FEPS(LEPS$par, y)), ylab="Sample Quantiles",
       xlab = "Theoretical Quantiles", pch=16)
qqline((qnorm(FEPS(LEPS$par, y))),lwd = 2)

