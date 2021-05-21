bnmr <- function(dfsnp,pheno,outcome,bn_method,mr_model,cutoff,repeats,nsam,n.iter){
  library("bnlearn")
  library(rstan)
  library(MendelianRandomization)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  learnBNpc <- function(df,nsam){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    model <- pc.stable(dfSam)
    dfarc <- data.frame(model$arcs)
    return(dfarc)
  }
  
  learnBNhc <- function(df,nsam){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    model <- hc(dfSam)
    dfarc <- data.frame(model$arcs)
    return(dfarc)
  }
  
  BNbootstrapPC <- function(df,repeats,nsam){
    arcsL <- replicate(repeats,learnBNpc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }
  
  BNbootstrapHC <- function(df,repeats,nsam){
    arcsL <- replicate(repeats,learnBNhc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }
  
  BNbootstrapBi <- function(df,repeats,nsam){
    arcsL <- replicate(repeats/2,learnBNpc(df,nsam),simplify = FALSE)
    arcsL2 <- replicate(repeats/2,learnBNhc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL2 <- do.call(rbind.fill,arcsL2)
    arcsL <- rbind(arcsL,arcsL2)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }
  
  getscore <- function(dfre,pheno,snp,repeats){
    #pheno is a str, snp is a vector of str.
    score <- rep(0,length(snp))
    for(i in 1:length(snp)){
      sn <- snp[i]
      count1 <- dfre[which(dfre$from==sn&dfre$to==pheno),"count"]
      count2 <- dfre[which(dfre$from==pheno&dfre$to==sn),"count"]
      if(length(count1)==0){
        count1 <- 0
      }
      if(length(count2)==0){
        count2 <- 0
      }
      score[i] <- (count1+count2)/repeats
    }
    dfscore <- data.frame(snp=snp,score=score)
    return(dfscore)
  }
  
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
  
  df1 <- cbind(dfsnp,pheno)
  
  if(bn_method=="pc"){
    dfre <- BNbootstrapPC(df1,repeats,nsam)
  }else if(bn_method=="hc"){
    dfre <- BNbootstrapHC(df1,repeats,nsam)
  }else if(bn_method=="mix"){
    dfre <- BNbootstrapBi(df1,repeats,nsam)
  }else{
    print("no this method")
    break
  }
  
  snp <- colnames(dfsnp)
  dfscore <- getscore(dfre,pheno,snp,repeats)
  selectsnp <- dfscore[which(dfscore$score>=cutoff),"snp"]
  
  s <- selectsnp
  N <- nrow(dfsnp)
  J <- length(s)
  X <- array(pheno,dim=N)
  Y <- array(outcome,dim=N)
  Z <- as.matrix(dfsnp[,s],dim=c(N,J))
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
  risultato = mr_allmethods(oggetto, method = "ivw")
  thetamedianestimate = risultato$Values[2,2]
  
  init_list = list(c1=list(theta=thetamedianestimate,
                           beta=rep(0,J),alpha=betaX,deltax=0,
                           deltay=0,u=rep(0,N)))
  if(mr_model=="linear"){
    fit <- stan(model_code=stanmodelcode, init=init_list, iter=n.iter,
                chains=1, verbose=F,data=mydata)
  }else if(mr_model=="logit"){
    fit <- stan(model_code=stanmodelcodeLogit, init=init_list, iter=n.iter,
                chains=1, verbose=F,data=mydata)
  }
  
  theta = extract(fit,pars='theta',permuted=FALSE)
  motheta = monitor(theta,digits_summary=5)
  return(theta)
}