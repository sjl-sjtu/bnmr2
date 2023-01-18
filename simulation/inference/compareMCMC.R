

library(rstan)
library(MendelianRandomization)
library(AER)
library(ivmodel)

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

stanmodelcodeLaplace <-'
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
    }
    model {
        X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
        Y   ~ normal(omegay+Z*gamma+X*beta+u*deltay, sigmay);
        u 	~ normal(0,1);    
        for(k in 1:J){
            alpha[k] ~ normal(mualpha, sigmaalpha);
            gamma[k] ~ double_exponential(0, 2*tau^2);
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
        real<lower=0> c;
        real<lower=0,upper=1> pi;
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
        real<lower=0> c;
        real<lower=0,upper=1> pi;
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

library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
df <- read_csv("inference/fixedIV_2_n12.csv")
df_iv <- read_csv("inference/snp_iv.csv")
selectsnp <- df_iv %>% arrange(desc(score)) %>% slice_head(n=20) %>% pull(snp)
snp <- grep("^g",colnames(df))
exposureName <- "X"
outcomeName <- "Y7"

# df <- df[sample(nrow(df),1000),]

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

oggetto <- mr_input(bx = as.numeric(betaX),
                    bxse = as.numeric(sebetaX),
                    by = as.numeric(betaY),
                    byse = as.numeric(sebetaY),
                    correlation = cor(Z),
                    exposure = "X ", outcome = "Y",
                    snps = colnames(Z))

init="median"
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
}else if(class(init)=="numeric"){
  betainitestimate <- init[1]
}else{
  stop("no this method to get initial value!")
}

betainitestimate



n.iter <- 5000
n.chain <- 4
init_list <- list(c1=list(beta=betainitestimate,
                          gamma=rep(0,J),alpha=betaX,deltax=0,
                          deltay=0,u=rep(0,N)))
init_list <- rep(init_list,n.chain)

set.seed(10)

# timestart <- Sys.time()
# fit1 <- stan(model_code=stanmodelcodeHorseshoe, init=init_list, iter=5000,
#             chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t1 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit1,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit1,pars="beta",inc_warmup=TRUE)+ggtitle("Horseshoe")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC5000_1.png", width=5.5, height=4, dpi = 300)
# b1 <- mobeta$mean
# s1 <- mobeta$sd
# r1 <- mobeta$Rhat
# l1 <- mobeta$`2.5%`
# u1 <- mobeta$`97.5%`
# bess1 <- mobeta$Bulk_ESS
# tess1 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit2 <- stan(model_code=stanmodelcodeHyperlasso, init=init_list, iter=5000,
#             chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t2 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit2,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit2,pars="beta",inc_warmup=TRUE)+ggtitle("Hyperlasso")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC5000_2.png", width=5.5, height=4, dpi = 300)
# b2 <- mobeta$mean
# s2 <- mobeta$sd
# r2 <- mobeta$Rhat
# l2 <- mobeta$`2.5%`
# u2 <- mobeta$`97.5%`
# bess2 <- mobeta$Bulk_ESS
# tess2 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit3 <- stan(model_code=stanmodelcodeSpikeSlabUniform, init=init_list, iter=5000,
#             chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t3 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit3,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit3,pars="beta",inc_warmup=TRUE)+ggtitle("Uniform Spike & Slab")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC5000_3.png", width=5.5, height=4, dpi = 300)
# b3 <- mobeta$mean
# s3 <- mobeta$sd
# r3 <- mobeta$Rhat
# l3 <- mobeta$`2.5%`
# u3 <- mobeta$`97.5%`
# bess3 <- mobeta$Bulk_ESS
# tess3 <- mobeta$Tail_ESS

timestart <- Sys.time()
fit4 <- stan(model_code=stanmodelcodeLasso, init=init_list, iter=5000,
              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
timeend <- Sys.time()
t4 <- difftime(timeend,timestart,units="min")
beta <- rstan::extract(fit4,pars='beta',permuted=FALSE)
mobeta <- monitor(beta,digits_summary=5)
traceplot(fit4,pars="beta",inc_warmup=TRUE)+ggtitle("Bayesian Lasso")+theme(plot.title = element_text(hjust = 0.5))
ggsave("MCMC5000_4.png", width=5.5, height=4, dpi = 300)
b4 <- mobeta$mean
s4 <- mobeta$sd
r4 <- mobeta$Rhat
l4 <- mobeta$`2.5%`
u4 <- mobeta$`97.5%`
bess4 <- mobeta$Bulk_ESS
tess4 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit5 <- stan(model_code=stanmodelcodeSpikeSlabBernoulli, init=init_list, iter=5000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t5 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit5,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit5,pars="beta",inc_warmup=TRUE)+ggtitle("Bernoulli Spike & Slab")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC5000_5.png", width=5.5, height=4, dpi = 300)
# b5 <- mobeta$mean
# s5 <- mobeta$sd
# r5 <- mobeta$Rhat
# l5 <- mobeta$`2.5%`
# u5 <- mobeta$`97.5%`
# bess5 <- mobeta$Bulk_ESS
# tess5 <- mobeta$Tail_ESS

# re <- tibble(prior=c("horseshoe","hyperlasso","spikeslabuniform","lasso","spikeslabbernoulli"),
#              t=c(t1,t2,t3,t4,t5),bias=c(b1,b2,b3,b4,b5)-2,se=c(s1,s2,s3,s4,s5),
#              rhat=c(r1,r2,r3,r4,r5),low=c(l1,l2,l3,l4,l5),up=c(u1,u2,u3,u4,u5),bulk_ess=c(bess1,bess2,bess3,bess4,bess5),
#              tail_ess=c(tess1,tess2,tess3,tess4,tess5)) %>% write_csv("compareMCMC_5000.csv")


tibble(prior="lasso",t=t4,b=b4-2,se=s4,rhat=r4,low=l4,up=u4,bulk_ess=bess4,taol_ess=tess4) %>% write_csv("compareMCMC_5000.csv",append=TRUE)

# n.iter <- 10000
# n.chain <- 4
# init_list <- list(c1=list(beta=betainitestimate,
#                           gamma=rep(0,J),alpha=betaX,deltax=0,
#                           deltay=0,u=rep(0,N)))
# init_list <- rep(init_list,n.chain)

# timestart <- Sys.time()
# fit1 <- stan(model_code=stanmodelcodeHorseshoe, init=init_list, iter=10000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t1 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit1,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit1,pars="beta",inc_warmup=TRUE)+ggtitle("Horseshoe")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC10000_1.png", width=5.5, height=5, dpi = 300)
# b1 <- mobeta$mean
# s1 <- mobeta$sd
# r1 <- mobeta$Rhat
# l1 <- mobeta$`2.5%`
# u1 <- mobeta$`97.5%`
# bess1 <- mobeta$Bulk_ESS
# tess1 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit2 <- stan(model_code=stanmodelcodeHyperlasso, init=init_list, iter=10000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t2 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit2,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit2,pars="beta",inc_warmup=TRUE)+ggtitle("Hyperlasso")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC10000_2.png", width=5.5, height=5, dpi = 300)
# b2 <- mobeta$mean
# s2 <- mobeta$sd
# r2 <- mobeta$Rhat
# l2 <- mobeta$`2.5%`
# u2 <- mobeta$`97.5%`
# bess2 <- mobeta$Bulk_ESS
# tess2 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit3 <- stan(model_code=stanmodelcodeSpikeSlabUniform, init=init_list, iter=10000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t3 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit3,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit3,pars="beta",inc_warmup=TRUE)+ggtitle("Uniform Spike & Slab")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC10000_3.png", width=5.5, height=5, dpi = 300)
# b3 <- mobeta$mean
# s3 <- mobeta$sd
# r3 <- mobeta$Rhat
# l3 <- mobeta$`2.5%`
# u3 <- mobeta$`97.5%`
# bess3 <- mobeta$Bulk_ESS
# tess3 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit4 <- stan(model_code=stanmodelcodeLasso, init=init_list, iter=10000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t4 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit4,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit4,pars="beta",inc_warmup=TRUE)+ggtitle("Bayesian Lasso")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC10000_4.png", width=5.5, height=5, dpi = 300)
# b4 <- mobeta$mean
# s4 <- mobeta$sd
# r4 <- mobeta$Rhat
# l4 <- mobeta$`2.5%`
# u4 <- mobeta$`97.5%`
# bess4 <- mobeta$Bulk_ESS
# tess4 <- mobeta$Tail_ESS

# timestart <- Sys.time()
# fit5 <- stan(model_code=stanmodelcodeSpikeSlabBernoulli, init=init_list, iter=10000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# timeend <- Sys.time()
# t5 <- difftime(timeend,timestart,units="min")
# beta <- rstan::extract(fit5,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit5,pars="beta",inc_warmup=TRUE)+ggtitle("Bernoulli Spike & Slab")+theme(plot.title = element_text(hjust = 0.5))
# ggsave("MCMC5_4.png", width=5.5, height=5, dpi = 300)
# b5 <- mobeta$mean
# s5 <- mobeta$sd
# r5 <- mobeta$Rhat
# l5 <- mobeta$`2.5%`
# u5 <- mobeta$`97.5%`
# bess5 <- mobeta$Bulk_ESS
# tess5 <- mobeta$Tail_ESS

# re <- tibble(prior=c("horseshoe","hyperlasso","spikeslabuniform","lasso","spikeslabbernoulli"),
#              t=c(t1,t2,t3,t4,t5),bias=c(b1,b2,b3,b4,b5)-0.5,se=c(s1,s2,s3,s4,s5),
#              rhat=c(r1,r2,r3,r4,r5),low=c(l1,l2,l3,l4,l5),up=c(u1,u2,u3,u4,u5),bulk_ess=c(bess1,bess2,bess3,bess4,bess5),
#              tail_ess=c(tess1,tess2,tess3,tess4,tess5)) %>% write_csv("compareMCMC_4.csv")


# 
# fit <- stan(model_code=stanmodelcodeLasso, init=init_list, iter=2000,
#              chains=4, verbose=F,data=mydata,control=list(adapt_delta=0.8))
# beta <- rstan::extract(fit,pars='beta',permuted=FALSE)
# mobeta <- monitor(beta,digits_summary=5)
# traceplot(fit,pars="beta",inc_warmup=TRUE)
# re <- list(betaList=beta,mean=mobeta$mean,se=mobeta$se_mean,sd=mobeta$sd,
#            lower=mobeta$`2.5%`,upper=mobeta$`97.5%`,Rhat=mobeta$Rhat,fit_detail=fit)