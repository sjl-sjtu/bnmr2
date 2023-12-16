library(doParallel)
library(foreach)
# n_cores <- parallel::detectCores()
n_cores <- 60
registerDoParallel(n_cores)
getDoParWorkers()

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("../bn.R")
library(AER)
library(ivmodel)
library(sisVIVE)
library(parallel)
library(tidyverse)
library(AER)
library(ivmodel)
library(sisVIVE)
library(MendelianRandomization)
library(cause)
library(R2BGLiMS)

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

# set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")
dfgwas <- df_iv
truesnp <- df_iv %>% filter(label==1) %>%pull(snp)

snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)
IV <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise
IV1 <- df_iv%>%filter(IV==1)%>%pull(snp) #bn

FUN2 <- function(r,i){

  # library(tidyverse)
  # library(AER)
  # library(ivmodel)
  # library(sisVIVE)
  # library(MendelianRandomization)
  # library(cause)
  # library(R2BGLiMS)
  
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
  # IV <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise

  #LD IV
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
  TSLS_cr <- as.numeric(TSLS-qnorm(0.975)*TSLS_se<=0 & TSLS+qnorm(0.975)*TSLS_se>=0)
  
  
  m <- ivmodel(outcome,exposure,Z)
  risultato <- LIML(m)
  liml <- risultato$point.est-beta
  liml_se <- risultato$std.err
  liml_ci <- risultato$ci
  liml_cr <- as.numeric(as.numeric(liml_ci[2])>=beta & as.numeric(liml_ci[1])<=beta)
  
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
  median_cr <- as.numeric(risultato$Values[3,5]>=beta & risultato$Values[3,4]<=beta)
  
  risultato <- mr_allmethods(oggetto, method = "egger")
  egger <- risultato$Values[7,2]-beta
  egger_se <- risultato$Values[7,3]
  egger_cr <- as.numeric(risultato$Values[7,5]>=beta & risultato$Values[7,4]<=beta)
  
  risultato <- mr_allmethods(oggetto, method = "ivw")
  ivw <- risultato$Values[4,2]-beta
  ivw_se <- risultato$Values[4,3]
  ivw_cr <- as.numeric(risultato$Values[4,5]>=beta & risultato$Values[4,4]<=beta)
  
  risultato <- mr_mbe(oggetto,weighting = "weighted")
  mode <- risultato@Estimate-beta
  mode_se <- risultato@StdError
  mode_cr <- as.numeric(risultato@CIUpper>=beta & risultato@CILower<=beta)


  #BN IV
  exposure <- df$X
  outcome <- df%>%pull(paste0("Y",i))
  s <- IV1
  N <- nrow(df)
  J <- length(s)
  X <- array(exposure,dim=N)
  Y <- array(outcome,dim=N)
  Z <- as.matrix(df[,s],dim=c(N,J))
  
  risultato <-  summary(ivreg(outcome~exposure|Z))
  TSLS1 <-  risultato$coefficients[2,1]-beta
  TSLS1_se <-  risultato$coefficients[2,2]
  TSLS1_cr <- as.numeric(TSLS-qnorm(0.975)*TSLS_se<=0 & TSLS+qnorm(0.975)*TSLS_se>=0)
  
  
  m <- ivmodel(outcome,exposure,Z)
  risultato <- LIML(m)
  liml1 <- risultato$point.est-beta
  liml1_se <- risultato$std.err
  liml1_ci <- risultato$ci
  liml1_cr <- as.numeric(as.numeric(liml_ci[2])>=beta & as.numeric(liml_ci[1])<=beta)
  
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
  median1 <- risultato$Values[3,2]-beta
  median1_se <- risultato$Values[3,3]
  median1_cr <- as.numeric(risultato$Values[3,5]>=beta & risultato$Values[3,4]<=beta)
  
  risultato <- mr_allmethods(oggetto, method = "egger")
  egger1 <- risultato$Values[7,2]-beta
  egger1_se <- risultato$Values[7,3]
  egger1_cr <- as.numeric(risultato$Values[7,5]>=beta & risultato$Values[7,4]<=beta)
  
  risultato <- mr_allmethods(oggetto, method = "ivw")
  ivw1 <- risultato$Values[4,2]-beta
  ivw1_se <- risultato$Values[4,3]
  ivw1_cr <- as.numeric(risultato$Values[4,5]>=beta & risultato$Values[4,4]<=beta)
  
  risultato <- mr_mbe(oggetto,weighting = "weighted")
  mode1 <- risultato@Estimate-beta
  mode1_se <- risultato@StdError
  mode1_cr <- as.numeric(risultato@CIUpper>=beta & risultato@CILower<=beta)


  re <- c(i,r,TSLS,TSLS_se,TSLS_cr,liml,liml_se,liml_cr,
          ivw,ivw_se,ivw_cr,egger,egger_se,egger_cr,median,
          median_se,median_cr,mode,mode_se,mode_cr,
          TSLS1,TSLS1_se,TSLS1_cr,liml1,liml1_se,liml1_cr,
          ivw1,ivw1_se,ivw1_cr,egger1,egger1_se,egger1_cr,median1,
          median1_se,median1_cr,mode1,mode1_se,mode1_cr) 
  names(re) <- NULL
  return(re)
}

# FUN3 <- function(i){
#   rep <- 2
#   res <- foreach(r = seq(1,rep), .combine = 'rbind',.packages = loaded_packages) %dopar% {
#     FUN2(r, i)
#   }
#   return(res)
# }

# k <- 2
# re <- foreach(i = seq(1,k),.combine = 'rbind',
#               .packages = loaded_packages) %dopar% {
#   FUN3(i)
# }

# set.seed(0)
k <- 13
reps <- 100
re <- foreach(i = seq(1,k),.combine = 'rbind',
              .packages = loaded_packages) %:%
              foreach(r = seq(1,reps),
                 .combine = 'rbind',.packages = loaded_packages) %dopar% {
                  FUN2(r, i)
                 }
stopImplicitCluster()


# k <- 13
# items <- seq(1,13)
# cores <- detectCores(logical = FALSE)
# cl <- makeCluster(cores)
# resL <- parLapply(cl, items, FUN2)
# stopCluster(cl)

saveRDS(re,"sim_supp_iv.rds")

re <- as.matrix(re)
colnames(re) <- c("outcome","repeats","TSLS_est","TSLS_se","TSLS_cr","liml_est","liml_se","liml_cr",
                  "ivw_est","ivw_se","ivw_cr","egger_est","egger_se","egger_cr",
                  "median_est","median_se","median_cr","mode_est","mode_se","mode_cr",
                  "TSLS1_est","TSLS1_se","TSLS1_cr","liml1_est","liml1_se","liml1_cr",
                  "ivw1_est","ivw1_se","ivw1_cr","egger1_est","egger1_se","egger1_cr",
                  "median1_est","median1_se","median1_cr","mode1_est","mode1_se","mode1_cr") 
re %>% as_tibble() %>% write_csv("mr_re_type_supp_IV3.csv")
stopImplicitCluster()


# dfre <- data.frame(t(matrix(unlist(resL),ncol=length(items))))

