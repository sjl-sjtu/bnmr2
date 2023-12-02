#' Causal estimation by Bayesian Mendelian randomization with shrinkage prior to cope with pleiotropy
#'
#'
#' @param df a data frame which contains data of IVs, specified exposure and outcome.
#' @param selectsnp a vector of string containing the column names of df corresponding to the IV used in MR.
#' @param exposureName a string which is a column name of df corresponding to the exposure studied.
#' @param outcomeName a string which is a column name of df corresponding to the outcome studied.
#' @param mr_model model for MR. Possible values are "linear" or "logit". Default is "linear".
#' @param prior a string represented shrinkage prior used in estimation. It can be "horseshoe", "lasso", "hyperlasso", "spikeslabBernoulli", "spikeslabUniform". Default is "horseshoe".
#' @param init the init value of theta for MCMC estimation. It can be a specific numeric or a string of "TSLS", "median", "egger" and "ivw", which means the initial value of the iteration will be calculated automatically by the above method. Default is "median".
#' @param n.iter an integer standing for the number of iterations. Default is 5000.
#' @param n.chain the number of chains in MCMC sampling. Default is 4.
#'
#' @return a list containing:
#'   \item{betaList}{a vector cantaining the posterior of the causal parameter of interest using MCMC sampling.}
#'   \item{mean}{the mean estimate of the causal parameter.}
#'   \item{se}{the standard error of the estimation.}
#'   \item{lower}{the lower boundary of the 95\% CI of the causal estimation.}
#'   \item{upper}{the upper boundary of the 95\% CI of the causal estimation.}
#'   \item{Rhat}{a indicator to measure the convergence (at convergence, Rhat <= 1.1).}
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
#' model <- mr(df,truesnp,"x","y")

#'

mr <- function(df,selectsnp,exposureName,outcomeName,mr_model="linear",prior="horseshoe",init="median",n.iter=5000,n.chain=4){
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
  re <- list(betaList=beta,mean=mobeta$mean,se=mobeta$sd,
             lower=mobeta$`2.5%`,upper=mobeta$`97.5%`,
             Rhat=mobeta$Rhat,fit_detail=fit)
  return(re)
}

