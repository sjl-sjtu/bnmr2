library(tidyverse)
# library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# source("../bn.R")
# source("mr.R")
library(AER)
library(ivmodel)
library(sisVIVE)
library(parallel)

set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")
dfgwas <- df_iv
truesnp <- df_iv %>% filter(label==1) %>%pull(snp)

# snp <- df_iv%>%filter(p<1e-2)%>%pull(snp)
snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)



FUN2 <- function(i){

  library(tidyverse)
  library(AER)
  library(ivmodel)
  library(sisVIVE)
  library(MendelianRandomization)

  set.seed(0)
  df <- read_csv("flexibleIV_n.csv")
  df_iv <- read_csv("snp_iv.csv")
  dfgwas <- df_iv
  truesnp <- df_iv %>% filter(label==1) %>%pull(snp)

  snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)

  # generate
  n <- nrow(df)
  # balanced
  sigma_y <- 1 #
  delta_y <- 0.1
  epsilon_y <- rnorm(n,0,sigma_y)
  beta0 <- 5
  beta <- 0.5

  snp <- grep("^g",colnames(df),value = TRUE)
  snp2 <- dfgwas%>%filter(p<1e-1)%>%pull(snp)
  pm <- 100
  snpy <- sample(reduce(list(snp,truesnp,snp2),setdiff),pm)
  gamma1 <- rnorm(pm,0.1,0.05)

  #
  df$Y1 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+0
  pt <- 50 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- sample(truesnp,pt)
  df$Y2 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 25+25
  pt <- 50 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y3 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+50
  pt <- 50 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y4 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 100+0
  pt <- 100 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- sample(truesnp,pt)
  df$Y5 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+50
  pt <- 100 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y6 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+100
  pt <- 100 #
  gamma <- rnorm(pt,0,0.05)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y7 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # directional
  # 50+0
  pt <- 50 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- sample(truesnp,pt)
  df$Y8 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 25+25
  pt <- 50 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y9 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+50
  pt <- 50 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y10 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 100+0
  pt <- 100 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- sample(truesnp,pt)
  df$Y11 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+50
  pt <- 100 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y12 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+100
  pt <- 100 #
  gamma <- rnorm(pt,0.05,0.05)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y13 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)
    
  # IV <- df_iv%>%filter(p<1e-2)%>%pull(snp)
  # IV <- df_iv%>%dplyr::arrange(desc(score))%>%slice_head(n=20)%>%pull(snp)
  IV <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise

  exposure <- df$X
  outcome <- df%>%pull(paste0("Y",i))
  s <- IV
  N <- nrow(df)
  J <- length(s)
  X <- array(exposure,dim=N)
  Y <- array(outcome,dim=N)
  Z <- as.matrix(df[,s],dim=c(N,J))
  
  risultato <-  summary(ivreg(outcome~exposure|Z))
  TSLS <-  risultato$coefficients[2,1]-beta
  TSLS_se <-  risultato$coefficients[2,2]
  
  m <- ivmodel(outcome,exposure,Z)
  risultato <- LIML(m)
  liml <- risultato$point.est-beta
  liml_se <- risultato$std.err
  
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
  
  risultato <- mr_allmethods(oggetto, method = "median")
  median <- risultato$Values[3,2]-beta
  median_se <- risultato$Values[3,3]
  
  risultato <- mr_allmethods(oggetto, method = "egger")
  egger <- risultato$Values[7,2]-beta
  egger_se <- risultato$Values[7,3]
  
  risultato <- mr_allmethods(oggetto, method = "ivw")
  ivw <- risultato$Values[4,2]-beta
  ivw_se <- risultato$Values[4,3]
  
  risultato <- mr_mbe(oggetto,weighting = "weighted")
  mode <- risultato@Estimate-beta
  mode_se <- risultato@StdError
  
  
  # sisVIVE
  snp <- df_iv%>%filter(p<1e-2)%>%pull(snp)
  exposure <- df$X
  outcome <- df%>%pull(paste0("Y",i))
  N <- nrow(df)
  J <- length(snp)
  Z <- as.matrix(df[,snp],dim=c(N,J))
  sis <- cv.sisVIVE(outcome,exposure,Z,K=10)
  sisvive_est <- sis$beta-beta
  
  boot <- function(){
    dfn <- df %>% sample_n(n(),replace = TRUE)
    Z <- as.matrix(df[,snp],dim=c(N,J))
    sis <- cv.sisVIVE(outcome,exposure,Z,K=10)
    return(sis$beta-beta)
  }
  samplers <- replicate(100,boot())
  sisvive_se <- sd(samplers)
  
  
  IV <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise
  D <- as.numeric(df$X)
  Y <- as.numeric(df%>%pull(paste0("Y",i)))
  s <- IV
  N <- nrow(df)
  J <- length(s)
  # X <- array(exposure,dim=N)
  # Y <- array(outcome,dim=N)
  Z <- as.matrix(df[,s],dim=c(N,J))
  # Z <- apply(Z,1,as.numeric)
  
  source("../CIIV_Functions.R")
  ciiv <- CIIV(Y,D,Z)
  ci_est <- ciiv$Coefficients_CIM[1]-beta
  cigmm_est <- ciiv$Coefficients_CIM_GMM[1]-beta
  ci_se <- ciiv$sd_CIM[1]
  cigmm_se <- ciiv$sd_CIM_GMM[1]
  ci_ci <- unname(ciiv$ci_CIM)
  cigmm_ci <- unname(ciiv$ci_CIM_GMM)

  re <- c(ivw,ivw_se,egger,egger_se,median,
          median_se,mode,mode_se,ci_est,ci_se,ci_ci,cigmm_est,cigmm_se,cigmm_ci,
          sisvive_est-beta,sisvive_se,TSLS-beta,TSLS_se,liml-beta,liml_se)
  return(re)
}

FUN3 <- function(i){
  rep <- 50
  res <- replicate(rep,FUN2(i))
  sumre <- rowMeans(res)
  # sumre <- mean(res)
  # sere <- sd(res)
  re <- c(paste0("Y",i),sumre)
  return(re)
}

k <- 13
items <- seq(1,13)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
resL <- parLapply(cl, items, FUN2)
stopCluster(cl)


dfre <- data.frame(t(matrix(unlist(resL),ncol=length(items))))

