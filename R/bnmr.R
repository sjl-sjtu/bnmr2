#' Causal inference between traits using the Bayesian Network-based Mendelian randomization
#'
#' @param df a data frame which contains data of SNPs and specified exposure. The values of snps in the data frame should be either numeric or factors (not integers) for BN learning.
#' @param snp a vector of string belonging to column names of df, which is the name of SNPs included in BN structure learning.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param outcomeName a string which is a column name of df corresponding to the outcome studied.
#' @param bn_method method for BN structure learning. Possible values are the function name of structure learning algorithm implemented in bnlearn. Default is "hc".
#' @param repeats an integer standing for the number of subsamples or bootstraps. Default is 1000.
#' @param selectNum the number of instrument to select. Default is NA.
#' @param alpha a number between 0 and 1 to specify the threshold for IV selection. We will use a threshold for variant selection as alpha*psam/length(snp). If selectNum is specified, the parameter will not be used. Default is 0.9.
#' @param nsam the size of individuals in each subsample of random graph forest. Default is 1000.
#' @param psam the size of variants in each subsample of random graph forest. Default is 100.
#' @param sample_replace is a boolean value to determine the sampling methods for individuals. TRUE with replacement and FALSE without replacement. Default is TRUE.
#' @param mr_model model for MR. Possible values are "linear" or "logit". Default is "linear".
#' @param prior a string represented shrinkage prior used in estimation. It can be "horseshoe", "lasso", "hyperlasso", "spikeslabBernoulli", "spikeslabUniform". Default is "horseshoe".
#' @param init the init value of theta for MCMC estimation. It can be a specific numeric or a string of "TSLS", "median", "egger" and "ivw", which means the initial value of the iteration will be calculated automatically by the above method. Default is "median".
#' @param n.iter an integer standing for the number of iterations. Default is 5000.
#' @param n.chain the number of chains in MCMC sampling. Default is 4.
#'
#' @return a list containing:
#'   \item{selectsnp}{a vector of string containing the colnames of df corresponding to selected SNPs.}
#'   \item{dfscore}{a data frame containing the score calculated for each SNP.}
#'   \item{betaList}{a vector cantaining the posterior of the causal parameter of interest using MCMC sampling.}
#'   \item{mean}{the mean estimate of the causal parameter.}
#'   \item{se}{the standard error of the estimation.}
#'   \item{lower}{the lower boundary of the 95\% CI of the causal estimation.}
#'   \item{upper}{the upper boundary of the 95\% CI of the causal estimation.}
#'   \item{Rhat}{an indicator to measure the convergence (at convergence, Rhat <= 1.1).}
#'   \item{fit_detail}{an S4 class stanfit object containing the details of Bayesian MR estimation}
#'
#' @export
#'
#' @examples
#' n <- 2000
#' p <- 200
#' snps <- replicate(p,sample(1:3,n,replace = TRUE))
#' snps <- apply(snps,2,as.numeric)
#' snpname <- paste0("g",1:p)
#' df <- as.data.frame(snps)
#' colnames(df) <- snpname
#' truesnp <- paste0("g",sample(1:p,50))
#' df$x <- as.matrix(df[,truesnp])%*%rnorm(50,0.05,0.05)+rnorm(n,0,1)
#' df$y <- 0.5*df$x+rnorm(n,0,1)
#' model <- bnmr(df,snpname,"x","y")
#'
bnmr <- function(df,snp,exposureName,outcomeName,bn_method="hc",repeats=1000,selectNum=NA,alpha=0.9,nsam=1000,psam=100,sample_replace=TRUE,mr_model="linear",prior="horseshoe",init="median",n.iter=5000,n.chain=4){
  library(bnlearn)
  library(plyr)
  library(dplyr)
  library(parallel)
  library(rstan)
  library(MendelianRandomization)
  library(AER)
  library(ivmodel)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- 2L
  } else {
    # use all cores in devtools::test()
    cores <- parallel::detectCores(logical = FALSE)
  }

  options(mc.cores = cores)
  rstan_options(auto_write = TRUE)

  learnBN <- function(iter,df,snp,exposureName,nsam,psam,bn_method,sample_replace){
    n <- nrow(df)
    if(sample_replace==TRUE){
      iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    }else{
      if(nsam>n){
        stop("subsample size is larger than the original sample size")
      }else{
        iSam <- sample(seq(1,n),size = nsam,replace=FALSE)
      }
    }
    if(psam>length(snp)){
      stop("subsample features is more than the number of features")
    }else{
      jSam <- sample(snp,size = psam, replace = FALSE)
    }
    dfSam <- df[iSam,c(jSam,exposureName)]
    rmFlag <- 0
    if(bn_method=="pc.stable"){
      model <- pc.stable(dfSam,undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="gs"){
      model <- gs(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="iamb"){
      model <- iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="fast.iamb"){
      model <- fast.iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="inter.iamb"){
      model <- inter.iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="iamb.fdr"){
      model <- iamb.fdr(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="hc"){
      model <- hc(dfSam)
    }else if(bn_method=="tabu"){
      model <- tabu(dfSam)
    }else if(bn_method=="mmhc"){
      model <- mmhc(dfSam)
    }else if(bn_method=="rsmax2"){
      model <- rsmax2(dfSam)
    }else if(bn_method=="h2pc"){
      model <- h2pc(dfSam)
    }else if(bn_method=="mmpc"){
      model <- mmpc(dfSam, undirected = TRUE)
    }else if(bn_method=="si.hiton.pc"){
      model <- si.hiton.pc(dfSam, undirected = TRUE)
    }else if(bn_method=="hpc"){
      model <- hpc(dfSam, undirected = TRUE)
    }else if(bn_method=="chow.liu"){
      model <- chow.liu(dfSam)
    }else if(bn_method=="aracne"){
      model <- aracne(dfSam)
    }else{
      stop("no this bn learning method")
    }
    dfarc <- data.frame(model$arcs)
    if(rmFlag==1){
      dfarc <- rmBidire(dfarc)
    }
    return(dfarc)
  }

  rmBidire <- function(df){
    df <- arrange(df,from,to)
    for(i in 1:nrow(df)){
      p <- which(df$from==df$from[i]|df$from==df$to[i])
      for(j in setdiff(p,i)){
        if(all(df[j,]%in%df[i,])){
          df <- df[-j,]
        }
      }
    }
    return(df)
  }

  getscore <- function(dfre,exposureName,snp,repeats){
    #exposureName is a str, snp is a vector of str.
    calc <- function(sn) {
      count1 <- dfre[which(dfre$from==sn&dfre$to==exposureName),"count"]
      count2 <- dfre[which(dfre$from==exposureName&dfre$to==sn),"count"]
      if(length(count1)==0){
        count1 <- 0
      }
      if(length(count2)==0){
        count2 <- 0
      }
      return((count1+count2)/repeats)
    }
    score <- sapply(snp, calc)
    dfscore <- data.frame(snp=snp,score=score)
    dfscore$snp <- as.character(dfscore$snp)
    return(dfscore)
  }

  BNbootstrap <- function(df,snp,exposureName,repeats,nsam,psam,bn_method,sample_replace){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {library("bnlearn")
      library("plyr")
      library("dplyr")
      library("parallel")
    })
    clusterExport(cl,deparse(substitute(learnBN)),envir=environment())
    arcsL <- parLapply(cl,seq(1,repeats),learnBN,df,snp,exposureName,nsam,psam,bn_method,sample_replace)
    stopCluster(cl)
    # arcsL <- replicate(repeats,learnBN(df,nsam,bn_method),simplify = FALSE)
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }

  stanmodelcodeHorseshoe <-'
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
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real beta;
        real<lower=0> tau;
        vector[N] u;
        vector[J] z;
        vector[J] alpha;
        vector[J] gamma;
        vector<lower=0>[J] phi;
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            phi[k] ~ cauchy(0, 1);
            gamma[k] ~ normal(0, phi[k]*tau);
        }
        // constructing the prior for tau
            tau ~ cauchy(0, 1);
        }
    '

  stanmodelcodeSpikeSlabUniform <-'
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
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real beta;
        vector<lower=0>[J] tau;
        vector[N] u;
        vector[J] z;
        vector[J] alpha;
        vector[J] gamma;
        vector<lower=0,upper=1>[J] lambda;
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            lambda[k] ~ uniform(0,1);
            gamma[k] ~ normal(0, lambda[k]*tau[k]);
            //prior for tau
            tau[k]  ~  inv_gamma(0.5,0.5);
        }
    }
    '

  stanmodelcodeSpikeSlabBernoulli <-'
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
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real beta;
        vector<lower=0>[J] tau;
        real<lower=0,upper=1> pi;
        vector[N] u;
        vector[J] z;
        vector[J] alpha;
        vector[J] gamma;
        vector<lower=0,upper=1>[J] lambda;
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            target += log_sum_exp(bernoulli_lpmf(0 | pi) + normal_lpdf(gamma[k] | 0, 0.001),
 						  bernoulli_lpmf(1 | pi) + normal_lpdf(gamma[k] | 0, tau[k]));
 						//prior for tau
            tau[k]  ~  inv_gamma(0.5,0.5);
        }
    }
    '

  stanmodelcodeLasso <-'
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
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real beta;
        vector[N] u;
        vector[J] z;
        vector[J] alpha;
        vector[J] gamma;
        real <lower=0> lambda;
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            gamma[k] ~ double_exponential(0,1/lambda);
        }
        // prior for lambda
        lambda ~ cauchy(0,1);
    }
    '

  stanmodelcodeHyperlasso <-'
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
        real mualpha;
        real omegax;
        real omegay;
        real deltax;
        real deltay;
        real beta;
        vector[N] u;
        vector[J] z;
        vector[J] alpha;
        vector[J] gamma;
        vector<lower=0>[J] tau;
        real<lower=0> lambda;
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            gamma[k] ~ double_exponential(0,sqrt(2*tau[k]));
            tau[k] ~ gamma(0.5,1/lambda^2);
        }
        // prior for lambda
        lambda ~ cauchy(0,1);
    }
    '


  stanmodelcodeLogitHorseshoe <-'
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
      real <lower=0> sigmaalpha;
      real mualpha;
      real omegax;
      real omegay;
      real deltax;
      real deltay;
      real beta;
      real<lower=0> tau;
      vector[N] u;
      vector[J] z;
      vector[J] alpha;
      vector[J] gamma;
      vector<lower=0>[J] phi;
    }
    model {
      X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
      Y   ~ bernoulli_logit(omegay+Z*gamma+X*beta+u*deltay);
      u 	~ normal(0,1);
      for(k in 1:J){
          alpha[k] ~ normal(mualpha, sigmaalpha);
          phi[k] ~ cauchy(0, 1);
          gamma[k] ~ normal(0, phi[k]*tau);
      }
      // constructing the prior for tau
          tau ~ cauchy(0, 1);
    }
    '

  stanmodelcodeLogitSpikeSlabUniform <-'
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
      real <lower=0> sigmaalpha;
      real mualpha;
      real omegax;
      real omegay;
      real deltax;
      real deltay;
      real beta;
      vector<lower=0>[J] tau;
      vector[N] u;
      vector[J] z;
      vector[J] alpha;
      vector[J] gamma;
      vector<lower=0,upper=1>[J] lambda;
    }
    model {
      X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
      Y   ~ bernoulli_logit(omegay+Z*gamma+X*beta+u*deltay);
      u 	~ normal(0,1);
      for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            lambda[k] ~ uniform(0,1);
            gamma[k] ~ normal(0, lambda[k]*tau[k]);
            //prior for tau
            tau[k]  ~  inv_gamma(0.5,0.5);
        }
    }
    '

  stanmodelcodeLogitSpikeSlabBernoulli <-'
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
      real <lower=0> sigmaalpha;
      real mualpha;
      real omegax;
      real omegay;
      real deltax;
      real deltay;
      real beta;
      vector<lower=0>[J] tau;
      real<lower=0,upper=1> pi;
      vector[N] u;
      vector[J] z;
      vector[J] alpha;
      vector[J] gamma;
      vector<lower=0,upper=1>[J] lambda;
    }
    model {
      X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
      Y   ~ bernoulli_logit(omegay+Z*gamma+X*beta+u*deltay);
      u 	~ normal(0,1);
      for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            target += log_sum_exp(bernoulli_lpmf(0 | pi) + normal_lpdf(gamma[k] | 0, 0.001),
 						  bernoulli_lpmf(1 | pi) + normal_lpdf(gamma[k] | 0, tau[k]));
 						//prior for tau
            tau[k]  ~  inv_gamma(0.5,0.5);
        }
    }
    '

  stanmodelcodeLogitLasso <-'
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
      real <lower=0> sigmaalpha;
      real mualpha;
      real omegax;
      real omegay;
      real deltax;
      real deltay;
      real beta;
      vector[N] u;
      vector[J] z;
      vector[J] alpha;
      vector[J] gamma;
      real <lower=0> lambda;
    }
    model {
      X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
      Y   ~ bernoulli_logit(omegay+Z*gamma+X*beta+u*deltay);
      u 	~ normal(0,1);
      for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            gamma[k] ~ double_exponential(0,1/lambda);
        }
        // prior for lambda
        lambda ~ cauchy(0,1);
    }
    '

  stanmodelcodeLogitHyperlasso <-'
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
      real <lower=0> sigmaalpha;
      real mualpha;
      real omegax;
      real omegay;
      real deltax;
      real deltay;
      real beta;
      vector[N] u;
      vector[J] z;
      vector[J] alpha;
      vector[J] gamma;
      vector<lower=0>[J] tau;
      real<lower=0> lambda;
    }
    model {
      X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
      Y   ~ bernoulli_logit(omegay+Z*gamma+X*beta+u*deltay);
      u 	~ normal(0,1);
      for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            gamma[k] ~ double_exponential(0,sqrt(2*tau[k]));
            tau[k] ~ gamma(0.5,1/lambda^2);
        }
        // prior for lambda
        lambda ~ cauchy(0,1);
    }
    '

  df <- as.data.frame(df)
  df <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]
  dfre <- BNbootstrap(df,snp,exposureName,repeats,nsam,psam,bn_method,sample_replace)

  dfscore <- getscore(dfre,exposureName,snp,repeats)
  dfscore <- arrange(dfscore,desc(score))

  if(!is.na(selectNum)){
    selectsnp <- dfscore[1:selectNum,"snp"]
  }else{
    selectsnp <- dfscore%>%filter(score>=alpha*psam/length(snp))%>%pull(snp)
  }

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

  if(mr_model=="linear"){
    for(isnp in 1:J){
      regX <- lm(X ~ Z[,isnp])
      regY <- lm(Y ~ Z[,isnp])
      betaX[isnp] <- summary(regX)$coefficients[2,1]
      sebetaX[isnp] <- summary(regX)$coefficients[2,2]
      betaY[isnp] <- summary(regY)$coefficients[2,1]
      sebetaY[isnp] <- summary(regY)$coefficients[2,2]
    }
    oggetto <- mr_input(bx = as.numeric(betaX),
                        bxse = as.numeric(sebetaX),
                        by = as.numeric(betaY),
                        byse = as.numeric(sebetaY),
                        correlation = cor(Z),
                        exposure = "X ", outcome = "Y",
                        snps = colnames(Z))
    if(init=="TSLS"){
      lm <- summary(ivreg(outcome~exposure|Z))
      betainitestimate <- lm$coefficients[2,1]
    }else if(init=="LIML"){
      m <- ivmodel(outcome,exposure,Z)
      risultato <- LIML(m)
      betainitestimate <- risultato$point.est
    }else if(init=="median"){
      if(J<3){
        stop("median initial needs at least 3 IVs")
      }
      risultato <- mr_allmethods(oggetto, method = "median")
      betainitestimate <- risultato$Values[3,2]
    }else if(init=="mode"){
      if(J<3){
        stop("mode initial needs at least 3 IVs")
      }
      risultato <- mr_mbe(oggetto,weighting = "weighted")
      betainitestimate <- risultato@Estimate
    }else if(init=="egger"){
      if(J<3){
        stop("egger initial needs at least 3 IVs")
      }
      risultato <- mr_allmethods(oggetto, method = "egger")
      betainitestimate <- risultato$Values[7,2]
    }else if(init=="ivw"){
      risultato <- mr_allmethods(oggetto, method = "ivw")
      betainitestimate <- risultato$Values[4,2]
    }else if(is.numeric(init)){
      betainitestimate <- init[1]
    }else{
      stop("no this method to get initial value!")
    }
  }else if(mr_model=="logit"){
    for(isnp in 1:J){
      regX <- lm(X ~ Z[,isnp])
      regY <- glm(Y ~ Z[,isnp],family=binomial(link='logit'))
      betaX[isnp] <- summary(regX)$coefficients[2,1]
      sebetaX[isnp] <- summary(regX)$coefficients[2,2]
      betaY[isnp] <- summary(regY)$coefficients[2,1]
      sebetaY[isnp] <- summary(regY)$coefficients[2,2]
    }
    oggetto <- mr_input(bx = as.numeric(betaX),
                        bxse = as.numeric(sebetaX),
                        by = as.numeric(betaY),
                        byse = as.numeric(sebetaY),
                        correlation = cor(Z),
                        exposure = "X ", outcome = "Y",
                        snps = colnames(Z))
    if(init=="median"){
      if(J<3){
        stop("median initial needs at least 3 IVs")
      }
      risultato <- mr_allmethods(oggetto, method = "median")
      betainitestimate <- risultato$Values[3,2]
    }else if(init=="mode"){
      if(J<3){
        stop("mode initial needs at least 3 IVs")
      }
      risultato <- mr_mbe(oggetto,weighting = "weighted")
      betainitestimate <- risultato@Estimate
    }else if(init=="egger"){
      if(J<3){
        stop("egger initial needs at least 3 IVs")
      }
      risultato <- mr_allmethods(oggetto, method = "egger")
      betainitestimate <- risultato$Values[7,2]
    }else if(init=="ivw"){
      risultato <- mr_allmethods(oggetto, method = "ivw")
      betainitestimate <- risultato$Values[4,2]
    }else if(is.numeric(init)){
      betainitestimate <- init[1]
    }else{
      stop("no this method to get initial value!")
    }
  }else{
    stop("no this MR model!")
  }


  init_list <- list(c1=list(beta=betainitestimate,
                            gamma=rep(0,J),alpha=betaX,deltax=0,
                            deltay=0,u=rep(0,N)))
  init_list <- rep(init_list,n.chain)
  if(mr_model=="linear"){
    if(prior=="horseshoe"){
      fit <- stan(model_code=stanmodelcodeHorseshoe, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="spikeslabUniform"){
      fit <- stan(model_code=stanmodelcodeSpikeSlabUniform, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="spikeslabBernoulli"){
      fit <- stan(model_code=stanmodelcodeSpikeSlabBernoulli, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="lasso"){
      fit <- stan(model_code=stanmodelcodeLasso, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="hyperlsso"){
      fit <- stan(model_code=stanmodelcodeHyperlasso, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else{
      stop("no this prior!")
    }
  }else if(mr_model=="logit"){
    if(prior=="horseshoe"){
      fit <- stan(model_code=stanmodelcodeLogitHorseshoe, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="spikeslabUniform"){
      fit <- stan(model_code=stanmodelcodeLogitSpikeSlabUniform, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="spikeslabBernoulli"){
      fit <- stan(model_code=stanmodelcodeLogitSpikeSlabBernoulli, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="lasso"){
      fit <- stan(model_code=stanmodelcodeLogitLasso, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else if(prior=="hyperlsso"){
      fit <- stan(model_code=stanmodelcodeLogitHyperlasso, init=init_list, iter=n.iter,
                  chains=n.chain, verbose=F,data=mydata,control=list(adapt_delta=0.85))
    }else{
      stop("no this prior!")
    }
  }else{
    stop("no this MR model!")
  }


  beta <- rstan::extract(fit,pars='beta',permuted=FALSE)
  mobeta <- monitor(beta,digits_summary=5)
  re <- list(IV=selectsnp,score=dfscore,
             betaList=beta,mean=mobeta$mean,se=mobeta$sd,
             lower=mobeta$`2.5%`,upper=mobeta$`97.5%`,
             Rhat=mobeta$Rhat,fit_detail=fit)
  return(re)
}
