mr <- function(df,selectsnp,exposureName,outcomeName,mr_model="linear",init="median",n.iter=500){
  options(mc.cores = parallel::detectCores())
  library(rstan)
  rstan_options(auto_write = TRUE)
  library(MendelianRandomization)
  library(ivreg)
  
  stanmodelcode <-'
    /* lg_t.stan */
    functions {
    // Vector square root
    vector vecsqrt(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1]){
            res[m] = sqrt(x[m]);
        }
    return res; }
    }
    data {
    int<lower=0> N;
    int<lower=0> J;
    matrix[N,J] Z;
    vector[N] X;
    vector[N] Y;
    }
    parameters {
    real <lower=0> sigmax;
    real <lower=0> sigmay;
    real <lower=0> sigmaalpha;
    real <lower=0> r1_global;
    real <lower=0> r2_global;
    real mualpha;
    real omegax;
    real omegay;
    real deltax;
    real deltay;
    real theta;
    vector[N] u;
    vector[J] z;
    vector<lower=0>[J] r1_local;
    vector<lower=0>[J] r2_local;
    vector[J] alpha;
    }
    transformed parameters {
    real<lower=0> tau;
    vector<lower=0> [J] lambda;
    vector[J] beta;
    tau      = r1_global * sqrt(r2_global);
    lambda	  = r1_local .* vecsqrt(r2_local);
    beta	    =  z .* lambda*tau;
    }
    model {
    X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
    Y   ~ normal(omegay+Z*beta+X*theta+u*deltay, sigmay);
    u 	~ normal(0,1);

    for(k in 1:J){
        alpha[k] ~ normal(mualpha, sigmaalpha);
    }
    // Constructing the prior for the lambda vector
        z ~ normal(0, 1);
        r1_local ~ normal(0.0, 1.0);
        r2_local ~ inv_gamma(0.5, 0.5);
    // Constructing the prior for tau
        r1_global ~ normal(0.0, 1.0);
        r2_global ~ inv_gamma(0.5, 0.5);
        }
    '

    stanmodelcodeLogit <-'
        /* lg_t.stan */
        functions {
        // Vector square root
        vector vecsqrt(vector x) {
            vector[dims(x)[1]] res;
            for (m in 1:dims(x)[1]){
                res[m] = sqrt(x[m]);
            }
        return res; }
        }
        data {
        int<lower=0> N;
        int<lower=0> J;
        matrix[N,J] Z;
        vector[N] X;
        int Y[N];
        }
        parameters {
        real <lower=0> sigmax;
        real <lower=0> sigmay;
        real <lower=0> sigmaalpha;
        real <lower=0> r1_global;
        real <lower=0> r2_global;
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real theta;
        vector[N] u;
        vector[J] z;
        vector<lower=0>[J] r1_local;
        vector<lower=0>[J] r2_local;
        vector[J] alpha;
        }
        transformed parameters {
        real<lower=0> tau;
        vector<lower=0> [J] lambda;
        vector[J] beta;
        vector<lower=0>[N] odds;
        vector<lower=0, upper=1>[N] prob;
        tau      = r1_global * sqrt(r2_global);
        lambda	  = r1_local .* vecsqrt(r2_local);
        beta	    =  z .* lambda*tau;
        for (i in 1:N){
            odds[i] = exp(omegay+Z[i]*beta+X[i]*theta+u[i]*deltay);
            prob[i] = odds[i] / (odds[i] + 1);
        }
        }
        model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ bernoulli(prob);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
        }
        // Constructing the prior for the lambda vector
            z ~ normal(0, 1);
            r1_local ~ normal(0.0, 1.0);
            r2_local ~ inv_gamma(0.5, 0.5);
        // Constructing the prior for tau
            r1_global ~ normal(0.0, 1.0);
            r2_global ~ inv_gamma(0.5, 0.5);
            }
        '

    df <- as.data.frame(df) #something will be wrong in converting array for tibble format
    exposure <- df[,exposureName]
    outcome <- df[,outcomeName]
    s <- selectsnp
    N <- nrow(df)
    J <- length(s)
    X <- array(exposure,dim=N)
    Y <- array(outcome,dim=N)
    Z <- as.matrix(df[,s],dim=c(N,J))
    mydata <- list(N=N,J=J,X=X,Y=Y,Z=Z)

    betaX <- array(NA, dim=J)
    betaY <- array(NA, dim=J)
    sebetaY <- array(NA, dim=J)
    sebetaX <- array(NA, dim=J)
    for(isnp in 1:J){
    regX <- lm(X ~ Z[,isnp])
    regY <- lm(Y ~ Z[,isnp])
    betaX[isnp] <- summary(regX)$coefficients[2,1]
    sebetaX[isnp] <- summary(regX)$coefficients[2,2]
    betaY[isnp] <- summary(regY)$coefficients[2,1]
    sebetaY[isnp] <- summary(regY)$coefficients[2,2]
    }

    oggetto = mr_input(bx = as.numeric(betaX),
                    bxse = as.numeric(sebetaX),
                    by = as.numeric(betaY),
                    byse = as.numeric(sebetaY),
                    correlation = cor(Z),
                    exposure = "X ", outcome = "Y",
                    snps = colnames(Z))
    if(init=="TSLS"){
    lm = ivreg(outcome~exposure|Z,method="OLS")
    thetamedianestimate = lm$coefficients[2]
    }else if(init=="median"){
    if(J<3){
        return(message("median initial needs at least 3 IVs"))
    }
    risultato = mr_allmethods(oggetto, method = "median")
    thetamedianestimate = risultato$Values[3,2]
    }else if(init=="egger"){
    if(J<3){
        return(message("egger initial needs at least 3 IVs"))
    }
    risultato = mr_allmethods(oggetto, method = "egger")
    thetamedianestimate = risultato$Values[7,2]
    }else if(init=="ivw"){
    risultato = mr_allmethods(oggetto, method = "ivw")
    thetamedianestimate = risultato$Values[4,2]
    }else if(class(init)=="numeric"){
    thetamedianestimate = init[1]
    }else{
    return(message("no this method to get initial value!"))
    }

    init_list = list(c1=list(theta=thetamedianestimate,
                            beta=rep(0,J),alpha=betaX,deltax=0,
                            deltay=0,u=rep(0,N)))
    if(mr_model=="linear"){
    fit <- stan(model_code=stanmodelcode, init=init_list, iter=n.iter,
                chains=1, verbose=F,data=mydata)
    }else if(mr_model=="logit"){
    fit <- stan(model_code=stanmodelcodeLogit, init=init_list, iter=n.iter,
                chains=1, cores="mc.cores",verbose=F,data=mydata)
    }else{
    return(message("no this MR model!"))
    }

    theta = extract(fit,pars='theta',permuted=FALSE)
    motheta = monitor(theta,digits_summary=5)
    re = list(thetaList=theta,mean=motheta$mean,se=motheta$se_mean,sd=motheta$sd,
            lower=motheta$`2.5%`,upper=motheta$`97.5%`,Rhat=motheta$Rhat)
    return(re)
}

library("tidyverse")
set.seed(111)
library(MendelianRandomization)
df <- read.csv("simData.csv",check.names=F)
N <- nrow(df)
epsilon <- runif(N,0,10)
pai <- runif(1,0.1,0.4)
GX <- sample(c(0,1,2),N,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
#V <- rnorm(n,30,5)
F1 <- df[,"F1"]
F2 <- df[,"F2"]
Z1 <- df[,"1GB"]
Z2 <- df[,"1GA"]
D1 <- 1.2*F1+epsilon
D2 <- 3*F2+epsilon
D3 <- 1.2*F1+3*F2+epsilon
D4 <- 1.2*F1+10*GX+epsilon
D5 <- 1.2*F1+30*Z1+epsilon
D6 <- 1.2*F1+30*Z1+20*Z2+epsilon
D7 <- 1.2*F1+3*F2+10*GX+30*Z1+20*Z2+epsilon
IV <- c("1GB","2GB","3GB","4GB","5GB","1GA","2GA","3GA")

outList <- paste("D",seq(1,7),sep="")

library(ivreg)
FUN <- function(i){
    library(MendelianRandomization)
    library(ivreg)
    df$outcome <- get(i)
    exposure <- df[,"F1"]
    outcome <- df[,"outcome"]
    s <- IV
    N <- nrow(df)
    J <- length(s)
    X <- array(exposure,dim=N)
    Y <- array(outcome,dim=N)
    Z <- as.matrix(df[,s],dim=c(N,J))
    mydata <- list(N=N,J=J,X=X,Y=Y,Z=Z)

    betaX <- array(NA, dim=J)
    betaY <- array(NA, dim=J)
    sebetaY <- array(NA, dim=J)
    sebetaX <- array(NA, dim=J)
    for(isnp in 1:J){
      regX <- lm(X ~ Z[,isnp])
      regY <- lm(Y ~ Z[,isnp])
      betaX[isnp] <- summary(regX)$coefficients[2,1]
      sebetaX[isnp] <- summary(regX)$coefficients[2,2]
      betaY[isnp] <- summary(regY)$coefficients[2,1]
      sebetaY[isnp] <- summary(regY)$coefficients[2,2]
    }

    oggetto = mr_input(bx = as.numeric(betaX),
                   bxse = as.numeric(sebetaX),
                   by = as.numeric(betaY),
                   byse = as.numeric(sebetaY),
                   correlation = cor(Z),
                   exposure = "X ", outcome = "Y",
                   snps = colnames(Z))
    risultato = mr_allmethods(oggetto, method = "median")
    median = risultato$Values[3,2]
    risultato = mr_allmethods(oggetto, method = "egger")
    egger = risultato$Values[7,2]
    risultato = mr_allmethods(oggetto, method = "ivw")
    ivw = risultato$Values[4,2]
    lm = ivreg(outcome~exposure|Z,method="OLS")
    TSLS = lm$coefficients[2]
    re <- mr(df,IV,"F1","outcome",n.iter=1000)
    R <- re$Rhat
    while(R>1.05){
       re <- mr(df,IV,"F1","outcome",n.iter=1000)
       R <- re$Rhat
    }
    bnmr = re$mean
    L <- re$lower
    U <- re$upper
    res <- c(median,egger,ivw,TSLS,bnmr,L,U,R)
    return(res)
}

library(parallel)
cl <- makeCluster(8)
clusterExport(cl,varlist=c("df",outList,"mr","N","IV"))
resL <- parLapply(cl, outList, FUN)
stopCluster(cl)
#resL <- lapply(outList,FUN)
dfre = data.frame(t(matrix(unlist(resL),ncol=length(outList))))
colnames(dfre) <- c("median","egger","ivw","TSLS","bnmr","L","U","R")
write_csv(dfre,"mr_re.csv")
