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
# source("../bn.R")
# source("mr.R")
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
source("../JMCode.R")
source("../CIIV_Functions.R")
source("bnmr.R")

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

# set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")
dfgwas <- df_iv
truesnp <- df_iv %>% filter(label==1) %>%pull(snp)

# snp_all <- df_iv%>%pull(snp)
snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)

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
  gamma1 <- rnorm(pm,0.1,0.01)

  #
  df$Y1 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+0
  pt <- 50 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- sample(truesnp,pt)
  df$Y2 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 25+25
  pt <- 50 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y3 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+50
  pt <- 50 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y4 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 100+0
  pt <- 100 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- sample(truesnp,pt)
  df$Y5 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+50
  pt <- 100 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y6 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+100
  pt <- 100 #
  gamma <- rnorm(pt,0,0.01)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y7 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # directional
  # 50+0
  pt <- 50 #
  gamma <- rnorm(pt,0.01,0.01)
  pleiosnp <- sample(truesnp,pt)
  df$Y8 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 25+25
  pt <- 50 #
  gamma <- rnorm(pt,0.01,0.01)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y9 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+50
  pt <- 50 #
  gamma <- rnorm(pt,0.01,0.01)
  pleiosnp <- sample(setdiff(snp1,truesnp),pt)
  df$Y10 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 100+0
  pt <- 100 #
  gamma <- rnorm(pt,0.01,0.01)
  pleiosnp <- sample(truesnp,pt)
  df$Y11 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 50+50
  pt <- 100 #
  gamma <- rnorm(pt,0.01,0.01)
  pleiosnp <- c(sample(truesnp,pt/2),sample(setdiff(snp1,truesnp),pt/2))
  df$Y12 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

  # 0+100
  pt <- 100 #
  gamma <- rnorm(pt,0.01,0.01)
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
  

  snp <- df_iv%>%filter(p<1e-3)%>%pull(snp)
  X <- df$X
  Y <- df%>%pull(paste0("Y",i))
  N <- nrow(df)
  J <- length(snp)
  Z <- as.matrix(df[,snp],dim=c(N,J))
  
  # sisVIVE
  # sis <- cv.sisVIVE(Y,X,Z,K=5)
  # sisvive_est <- sis$beta-beta
  # lam <- sis$lambda
  # boot <- function(){
  #   dfn <- df %>% sample_n(n(),replace = TRUE)
  #   Z <- as.matrix(df[,snp],dim=c(N,J))
  #   sis <- cv.sisVIVE(Y,X,Z,lambdaSeq=lam,K=1)
  #   return(sis$beta-beta)
  # }
  # samplers <- replicate(100,boot())
  # sisvive_se <- sd(samplers)

  # CAUSE
  betaX <- c()
  betaY <- c()
  sebetaX <- c()
  sebetaY <- c()
  px <- c()
  py <- c()
  for(isnp in 1:J){
    regX <- lm(X ~ Z[,isnp])
    regY <- lm(Y ~ Z[,isnp])
    betaX[isnp] <- summary(regX)$coefficients[2,1]
    sebetaX[isnp] <- summary(regX)$coefficients[2,2]
    px[isnp] <- summary(regX)$coefficients[2,4]
    betaY[isnp] <- summary(regY)$coefficients[2,1]
    sebetaY[isnp] <- summary(regY)$coefficients[2,2]
    py[isnp] <- summary(regY)$coefficients[2,4]
  }
  df_summary <- tibble(snp=snp,beta_hat_1=betaX,seb1=sebetaX,p1=px,
                       beta_hat_2=betaY,seb2=sebetaY,p2=py)

  params <- est_cause_params(df_summary, snp)
  r2_thresh = 0.01
  pval_thresh = 1e-3
  # df_clump <- df_summary %>%
  #           rename(rsid = snp,
  #                   pval = p1) %>%
  #           ieugwasr::ld_clump(dat = .,
  #                     clump_r2 = r2_thresh,
  #                     clump_p = pval_thresh,
  #                     #plink_bin = genetics.binaRies::get_plink_binary(), 
  #                     pop = "EUR")
  # top_vars <- df_clump$rsid 
  top_vars <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise
  res <- cause(X=df_summary, variants = top_vars, param_ests = params,force=TRUE)
  rest <- res$elpd
  res_sum <- summary(res, ci_size=0.95)
  cause_e <- res_sum$tab[2,2]
  nums <- as.numeric(unlist(strsplit(gsub("[\\(\\),]", "", cause_e), " ")))
  cause_est <- as.numeric(nums[1])-beta
  cause_se <- (as.numeric(nums[3])-as.numeric(nums[2]))/(2*qnorm(0.975))
  cause_cr <- as.numeric(as.numeric(nums[3])>=beta & as.numeric(nums[2])<=beta)
  
  # JAMMR
  # eafs <- sapply(1:J,function(j) (2*sum(Z[,j]==2)+sum(Z[,j]==1))/(2*nrow(Z)))
  # G.cor <- cor(Z)
  es <- JAMMR(betaX,sebetaX,betaY,sebetaY,N1=dim(Z)[1],G.matrix = Z,n.grid=10,iter = 20000,jam.seed=123)#eafs=eafs,G.cor = G.cor,)
  jam_est <- es$causal-beta
  jam_se <- es$se
  jam_cr <- as.numeric(jam_est-qnorm(0.975)*jam_se<=0 & jam_est+qnorm(0.975)*jam_se>=0)
  
  
  # CIIV
  IV <- df_iv%>%filter(IV3==1)%>%pull(snp) #ld stepwise
  D <- as.numeric(df$X)
  Y <- as.numeric(df%>%pull(paste0("Y",i)))
  s <- IV
  N <- nrow(df)
  J <- length(s)
  Z <- as.matrix(df[,s],dim=c(N,J))
  
  
  ciiv <- CIIV(Y,D,Z)
  ci_est <- ciiv$Coefficients_CIM[1]-beta
  cigmm_est <- ciiv$Coefficients_CIM_GMM[1]-beta
  ci_se <- ciiv$sd_CIM[1]
  cigmm_se <- ciiv$sd_CIM_GMM[1]
  ci_ci <- unname(ciiv$ci_CIM)
  ci_cr <- as.numeric(as.numeric(ci_ci[2])>=beta & as.numeric(ci_ci[1])<=beta)
  cigmm_ci <- unname(ciiv$ci_CIM_GMM)
  cigmm_cr <- as.numeric(as.numeric(cigmm_ci[2])>=beta & as.numeric(cigmm_ci[1])<=beta)

  
  #bnmr
  re <- bnmr(df,snp,"X",paste0("Y",i),nsam=1000,psam=100,selectNum=20,repeats=1000,prior="lasso",n.iter=2000,n.chain=4)
  bnmr_est <- re$mean-beta
  bnmr_se <- re$se
  L <- re$lower
  U <- re$upper
  bnmr_cr <- as.numeric((L<=beta)&(U>=beta))

  re <- c(i,r,ivw,ivw_se,ivw_cr,egger,egger_se,egger_cr,median,median_se,median_cr,
          mode,mode_se,mode_cr,ci_est,ci_se,ci_cr,cigmm_est,cigmm_se,cigmm_cr,
          cause_est,cause_se,cause_cr,jam_est,jam_se,jam_cr,sisvive_est,
          TSLS,TSLS_se,TSLS_cr,liml,liml_se,liml_cr, bnmr_est,bnmr_se,coverage) #sisvive_se,
  names(re) <- NULL
  return(re)
}

k <- 13
reps <- 100
re <- foreach(i = seq(1,k),.combine = 'rbind',
              .packages = loaded_packages) %:%
              foreach(r = seq(1,reps),
                 .combine = 'rbind',.packages = loaded_packages) %dopar% {
                  FUN2(r, i)
                 }
stopImplicitCluster()

saveRDS(re,"sim_supp1.rds")

re <- as.matrix(re)
colnames(re) <- c("outcome","repeats","ivw_est","ivw_se","ivw_cr","egger_est","egger_se","egger_cr",
                  "median_est","median_se","median_cr","mode_est","mode_se","mode_cr",
                  "ci_est","ci_se","ci_cr","cigmm_est","cigmm_se","cigmm_cr",
                  "cause_est","cause_se","cause_cr","jam_est","jam_se","jam_cr",
                  "sisvive_est","TSLS_est","TSLS_se","TSLS_cr","liml_est","liml_se","liml_cr"
                  "bnmr_est","bnmr_se","bnmr_cr") #sisvive_se,

re %>% as_tibble() %>% write_csv("mr_re_type_supp2.csv")
stopImplicitCluster()



# FUN3 <- function(i){
#   rep <- 100
#   res <- foreach(r = seq(1,rep), .combine = 'rbind',.packages = loaded_packages) %dopar% {
#     FUN2(r, i)
#   }
#   return(res)
# }
# k <- 100
# re <- foreach(i = seq(1,k),.combine = 'rbind',
#               .packages = loaded_packages) %dopar% {
#   FUN3(i)
# }
# k <- 13
# items <- seq(1,13)
# cores <- detectCores(logical = FALSE)
# cl <- makeCluster(cores)
# resL <- parLapply(cl, items, FUN2)
# stopCluster(cl)
# dfre <- data.frame(t(matrix(unlist(resL),ncol=length(items))))

