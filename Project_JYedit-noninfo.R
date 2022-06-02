#Added comment###

rm(list=ls())

library(MASS)
library(truncnorm)
library(mvtnorm)

PGA<- as.data.frame(read.csv("PGASample.csv", header = TRUE,stringsAsFactors = T))

## Treat as Factor
PGA$PlayerName<-as.factor(PGA$PlayerName)
PGA$Round<-as.factor(PGA$Round)
PGA$ParValue<-as.factor(PGA$ParValue)
PGA$RoundCat<-ifelse(PGA$Round==4,"Round4","Round123")
PGA$TotalPutts[PGA$TotalPutts>=3]<-3
#PGA$Distance<-PGA$Distance/12

## Baseline 
PGA$Slope<-relevel(as.factor(PGA$Slope), ref ="Level")
PGA$ParValue<-relevel(as.factor(PGA$ParValue), ref ="3")
PGA$TotalPutts<-relevel(as.factor(PGA$TotalPutts), ref ="1")
PGA$Conversion<-relevel(as.factor(PGA$Conversion), ref ="Birdie")
PGA$PlayerName<-relevel(as.factor(PGA$PlayerName), ref ="Others")

## Subsetted Data
PGA2017<-PGA[PGA$Year==2017,]
PGA2017<- droplevels(PGA2017)

PGA1516<-PGA[(PGA$Year==2015 | PGA$Year==2016),]
PGA1516<-PGA1516[PGA1516$PlayerName != "Tiger Woods" & PGA1516$PlayerName != "TigerWoods",] # to remove "Tiger Woods" and "TigerWoods"
PGA1516<- droplevels(PGA1516)


################### OLR 2015-2016
m1.2 <- polr(TotalPutts ~ Distance+Conversion+PlayerName, 
             data = PGA1516, Hess=TRUE,method = c("logistic"))

## Estimate of mean vector and var-cov matrix

mu0 <- m1.2$coefficients # dimension = 34
sigma0 <- vcov(m1.2)[1:(dim(vcov(m1.2))[2]-2), 1:(dim(vcov(m1.2))[2]-2)] # remove two intercepts
sd1 <- sqrt(vcov(m1.2)[(dim(vcov(m1.2))[2]-1), (dim(vcov(m1.2))[2]-1)]) ; sd2 <- sqrt(vcov(m1.2)[(dim(vcov(m1.2))[2]),(dim(vcov(m1.2))[2])]); 
# standard error of intercepts results from m1.2

## Response y following multinominal

N <- length(PGA2017$TotalPutts) # the number of observations

y <- matrix(nrow = N, ncol = 3)
for(i in 1:N){
  if(PGA2017$TotalPutts[i]==1){
    y[i,] <- c(1,0,0) 
  }
  if(PGA2017$TotalPutts[i]==2){
    y[i,] <- c(0,1,0) 
  }
  if(PGA2017$TotalPutts[i]==3){
    y[i,] <- c(0,0,1) 
  }
}
head(y)

## Design matrix X

m1.3 <- polr(TotalPutts ~ Distance+Conversion+PlayerName, 
             data = PGA2017, Hess=TRUE,method = c("logistic"))
design <- model.matrix(m1.3)
design <- design[,-1]
dim(design)[2] == length(m1.2$coefficients) # checking dimension

## combine response and design matrix for MCMC
yx <- cbind(y, design)

########################################################
#NonInformative


################### Constructing full conditional density

## alpha1

cond.alpha1 <- function(yx, alpha1, alpha2, beta){
  temp <- -(yx[,c(-1:-3)])%*%beta
  temp2<-(yx[,1]*(alpha1 + temp - log(1+exp(alpha1 + temp))) + yx[,2]*log((exp(alpha2 + temp)/(1+exp(alpha2 + temp)))-(exp(alpha1+temp)/(1+exp(alpha1 + temp)))))  
  temp3<-sum(temp2)}

## alpha2

cond.alpha2 <- function(yx, alpha1, alpha2, beta){
  
  temp <- -(yx[,c(-1:-3)])%*%beta
  temp2<-yx[,2]*log((exp(alpha2 + temp)/(1+exp(alpha2 + temp)))-(exp(alpha1+temp)/(1+exp(alpha1 + temp)))) - yx[,3]*log(1+exp(alpha2 + temp))
  temp3<-sum(temp2)}


## beta(vector)
cond.beta <- function(yx, alpha1, alpha2, beta){
  
  temp <- -yx[,c(-1:-3)]%*%beta
  temp2<-yx[,1]*(alpha1 + temp - log(1+exp(alpha1 + temp)))+ yx[,2]*log((exp(alpha2 + temp)/(1+exp(alpha2 + temp)))-(exp(alpha1+temp)/(1+exp(alpha1 + temp)))) - yx[,3]*log(1+exp(alpha2 + temp))
  temp3<-sum(temp2)}


################### Implement MH and RWMH with the proposal 
################### N(alpha1; alpha1^(t),10^2)
################### N(alpha2; alpha2^(t),10^2, lower=alpha1) 
################### MVN(beta; beta^(t), sigma0)

## Initial values

alpha1 <- m1.2$zeta[1] # from OLR with 2015 2016 data
alpha2 <- m1.2$zeta[2] # from OLR with 2015 2016 data
beta <- m1.2$coefficients # from OLR with 2015 2016 data # dim = 34

T <- 100; ### number of outer iterations of the chain
T1 <- 100; ### number of inner iterations of the alpha1 chain
T2 <- 10; ### number of inner iterations of the alpha2 chain
T3 <- 100; ### number of inner iterations of the beta chain
B <- T/2; ### number of burn-in iterations
alpha1_store <- rep(NA, T); ### storage for the values of the chain
alpha2_store <- rep(NA, T); ### storage for the values of the chain
beta_store <- matrix(rep(NA, length(beta)*T), ncol = length(beta)); ### storage for the values of the chain
accept1 <- 0 ; accept2 <- 0 ; accept3 <- 0 ; ### acceptance count

go <- proc.time()
set.seed(12341234)
for (t in 1:T) {
  
  for(t1 in 1:T1)
  {
    # draw new sample from proposal distribution given alpha1^(t)
    alpha1_star <- rnorm(1, alpha1, sd1)
    # symmertic proposal so random walk M-H
    logMH <- cond.alpha1(yx, alpha1_star, alpha2, beta) - cond.alpha1(yx, alpha1, alpha2, beta);
    MH <- exp(logMH);
    u <- runif(1); ### uniform variable to determine acceptance
    if (u < MH) {
      alpha1 <- alpha1_star
    }
  }
  if (u < MH) {
    accept1 <- accept1 + 1;
  }
  alpha1_store[t] <- alpha1;
  
  for(t1 in 1:T2)
  {
    
    # draw new sample from proposal distribution given alpha2^(t)
    alpha2_star <- rtruncnorm(1, a= alpha1, mean = alpha2, sd = sd2)
    # assymmertic proposal so plain M-H
    logMH1 <- cond.alpha2(yx, alpha1, alpha2_star, beta) - cond.alpha2(yx, alpha1, alpha2, beta) + log(dtruncnorm(alpha2, a = alpha1, mean = alpha2 , sd = sd2)) - log(dtruncnorm(alpha2_star, a = alpha1, mean = alpha2 , sd = sd2));
    MH1 <- exp(logMH1);
    u <- runif(1); ### uniform variable to determine acceptance
    if (u < MH1) {
      alpha2 <- alpha2_star
    }
  }
  if (u < MH1) {
    accept2 <- accept2 + 1;
  }
  alpha2_store[t] <- alpha2;
  
  for(t1 in 1:T3)
  {
    
    # draw new sample from proposal distribution given beta^(t)
    beta_star <- rmvnorm(1, mean = beta, sigma = sigma0)
    beta_star <- as.vector(beta_star)
    # symmertic proposal so random walk M-H
    logMH2 <- cond.beta(yx, alpha1, alpha2, beta = beta_star) - cond.beta(yx, alpha1, alpha2, beta = beta)
    MH2 <- exp(logMH2);
    u <- runif(1); ### uniform variable to determine acceptance
    if (u < MH2) {
      beta <- beta_star
    }
  }
  if (u < MH2) {
    accept3 <- accept3 + 1;
  }
  
  beta_store[t,] <- beta;
  
  print(t)
}

proc.time()-go

accept1/T;accept2/T;accept3/T;

############ Comparison with frequentist

result_bayes <- c(mean(alpha1_store[(B + 1):T]), mean(alpha2_store[(B + 1):T]), apply(beta_store[(B + 1):T,],2,mean))
result_freq <- c(m1.3$zeta[1], m1.3$zeta[2], m1.3$coefficients)
comp<-cbind(result_bayes, result_freq) 


############ Plot of the result

###### Graphics for alpha1

plot(1:T, alpha1_store, type = "l", xlab = "Iteration", ylab = expression(alpha1), main = paste("alpha1", "acceptance rate=", round(accept1/T, 3)));
abline(h = m1.2$zeta[1], col = 2); # initial value
hist(alpha1_store[(B + 1):T], main = "Histogram of the sampled values(alpha1)", xlab = expression(alpha_1));
abline(v = alpha1, col = 2); # initial value
abline(v = mean(alpha1_store[(B + 1):T]), col = 4, lty = 2);
abline(v = quantile(alpha1_store[(B + 1):T], 0.025), col = 3, lty = 2);
abline(v = quantile(alpha1_store[(B + 1):T], 0.975), col = 3, lty = 2);

###### Graphics for alpha2

plot(1:T, alpha2_store, type = "l", xlab = "Iteration", ylab = expression(alpha2), main = paste("alpha2", "acceptance rate=", round(accept2/T, 3)));
abline(h = m1.2$zeta[2], col = 2); # initial value
hist(alpha2_store[(B + 1):T], main = "Histogram of the sampled values(alpha1)", xlab = expression(alpha_1));
abline(v = alpha2, col = 2); # initial value
abline(v = mean(alpha2_store[(B + 1):T]), col = 4, lty = 2);
abline(v = quantile(alpha2_store[(B + 1):T], 0.025), col = 3, lty = 2);
abline(v = quantile(alpha2_store[(B + 1):T], 0.975), col = 3, lty = 2);

###### Graphics for beta

par(mfrow = c(2,1))
b <- 1
plot(1:T, beta_store[,b], type = "l", xlab = "Iteration", ylab = expression(beta1), main = paste("beta1", "acceptance rate=", round(accept3/T, 3)));
abline(h = m1.2$coefficients[b], col = 2); # initial value
hist(beta_store[(B + 1):T, b], main = "Histogram of the sampled values(beta1)", xlab = expression(beta1));
abline(v = beta[b], col = 2); # initial value
abline(v = mean(beta_store[(B + 1):T,b]), col = 4, lty = 2);
abline(v = quantile(beta_store[(B + 1):T,b], 0.025), col = 3, lty = 2);
abline(v = quantile(beta_store[(B + 1):T,b], 0.975), col = 3, lty = 2);


## convergence diagnosis
par(mfrow= c(2,1))
acf(alpha1_store[(B + 1):T], main = "Series alpha1")
acf(alpha2_store[(B + 1):T], main = "Series alpha2")
acf(beta_store[(B + 1):T,1], main = "Series beta1")
acf(beta_store[(B + 1):T,2], main = "Series beta2")
acf(beta_store[(B + 1):T,8], main = "Series beta4")


#############DIC calculation#########

ll<- function(yx, alpha1, alpha2, beta){
  temp <- -yx[,c(-1:-3)]%*%beta
  temp2<-yx[,1]*(alpha1 + temp - log(1+exp(alpha1 + temp)))+ yx[,2]*log((exp(alpha2 + temp)/(1+exp(alpha2 + temp)))-(exp(alpha1+temp)/(1+exp(alpha1 + temp)))) - yx[,3]*log(1+exp(alpha2 + temp))
  return(sum(temp2))
}

llist<-rep(0,T-B)
for (i in (B+1):T){
  llist[i]<-ll(yx,alpha1_store[i],alpha2_store[i],beta_store[i,])
}

DIC<-mean(-2*llist)-ll(yx,mean(alpha1_store),mean(alpha2_store),colMeans(beta_store))

##########

est.b=colMeans(data.frame(alpha1_store,alpha2_star,beta_store)[B:T,])
names(est.b)=c(names(m1.2$zeta),names(m1.2$coefficients))

est.2017=c(m1.3$zeta,m1.3$coefficients)


cbind(est.b,est.2017)

m1.3$fitted.values

a1 = mean(alpha1_store[B:T])
a2 = mean(alpha2_store[B:T])
beta1 = colMeans(beta_store[B:T,]) 

XB1 = -(yx[,c(-1:-3)])%*%beta1

pi.1 = exp(a1+XB1)/(1+exp(a1+XB1))
pi.2 = exp(a2+XB1)/(1+exp(a2+XB1))-exp(a1+XB1)/(1+exp(a1+XB1))
pi.3 = 1/(1+exp(a2+XB1))

Y.pred = c()

for(i1 in 1:length(pi.1))
{
  V1 = c(pi.1[i1],pi.2[i1],pi.3[i1])
  imax = which.max(V1)
  Y.pred[i1] = imax
  
}

Y.pred.freq = c()

fitted=m1.3$fitted.values

for(i1 in 1:length(fitted[,1]))
{
  V2 = c(fitted[i1,1],fitted[i1,2],fitted[i1,3])
  imax = which.max(V2)
  Y.pred.freq[i1] = imax
  
}



pred.pi = data.frame(pi.1,pi.2,pi.3,m1.3$fitted.values)

putts2017=PGA2017$TotalPutts

diff.pred.bayesian = Y.pred - as.numeric(putts2017)


diff.pred.bayesian = Y.pred - as.numeric(putts2017)
diff.pred.freq = Y.pred.freq - as.numeric(putts2017)

match.pred=data.frame(Y.pred,Exact=PGA2017$TotalPutts)

sum(diff.pred.bayesian==0)/(length(diff.pred.bayesian))
sum(diff.pred.freq==0)/(length(diff.pred.freq))




