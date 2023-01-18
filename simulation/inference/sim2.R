library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("../bn.R")
source("mr.R")

set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")

IV <- df_iv%>%filter(IV==1)%>%pull(snp)
IV2 <- df_iv%>%filter(IV2==1)%>%pull(snp)
IV3 <- df_iv%>%filter(IV3==1)%>%pull(snp)
IV4 <- paste0("PC",seq(1,40))
IV5 <- df_iv%>%filter(IV5==1)%>%pull(snp)


FUN <- function(i,IV,IVmethod){
  library(MendelianRandomization)
  library(AER)
  library(ivmodel)
  library(tidyverse)
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  exposure <- df$X
  outcome <- df%>%pull(paste0("Y",i))
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
  
  oggetto <- mr_input(bx = as.numeric(betaX),
                      bxse = as.numeric(sebetaX),
                      by = as.numeric(betaY),
                      byse = as.numeric(sebetaY),
                      correlation = cor(Z),
                      exposure = "X ", outcome = "Y",
                      snps = colnames(Z))
  
  risultato <- mr_allmethods(oggetto, method = "median")
  median <- risultato$Values[3,2]
  median_se <- risultato$Values[3,3]
  
  risultato <- mr_allmethods(oggetto, method = "egger")
  egger <- risultato$Values[7,2]
  egger_se <- risultato$Values[7,3]
  
  risultato <- mr_allmethods(oggetto, method = "ivw")
  ivw <- risultato$Values[4,2]
  ivw_se <- risultato$Values[4,3]
  
  risultato <- mr_mbe(oggetto,weighting = "weighted")
  mode <- risultato@Estimate
  mode_se <- risultato@StdError
  
  risultato <-  summary(ivreg(outcome~exposure|Z))
  TSLS <-  risultato$coefficients[2,1]
  TSLS_se <-  risultato$coefficients[2,2]
  
  m <- ivmodel(outcome,exposure,Z)
  risultato <- LIML(m)
  liml <- risultato$point.est
  liml_se <- risultato$std.err
  
  re <- mr(df,IV,"X",paste0("Y",i),prior="horseshoe",n.iter=5000,n.chain=4)
  R <- re$Rhat
  # while(R>1.1){
  #   re <- mr(df,IV,"X",paste("Y",i,sep=""),prior="horseshoe",n.iter=2000,n.chain=4)
  #   R <- re$Rhat
  # }
  bnmr <- re$mean
  bnmr_se <- re$se
  bnmr_sd <- re$sd
  L <- re$lower
  U <- re$upper
  res <- c(IVmethod,paste0("Y",i),median,median_se,egger,egger_se,ivw,ivw_se,mode,mode_se,TSLS,TSLS_se,liml,liml_se,bnmr,bnmr_se,bnmr_sd,L,U,R)
  return(res)
}

k <- 13
library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV"))
resL <- parLapply(cl, seq(1,k), FUN, IV,"bn")
stopCluster(cl)
# resL <- lapply(seq(1,k), FUN, IV,"bn")
dfre <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre) <- c("method","outcome","median","median_se","egger","egger_se","ivw","ivw_se","mode","mode_se","TSLS","TSLS_se","liml","liml_se","bnmr","bnmr_se","bnmr_sd","L","U","R")
# write_csv(dfre,"mr_re_50_bn_b.csv")
dfre

library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV2"))
resL <- parLapply(cl, seq(1,k), FUN, IV2,"gwas")
stopCluster(cl)
# resL <- lapply(seq(1,k), FUN, IV2,"gwas")
dfre1 <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre1) <- c("method","outcome","median","median_se","egger","egger_se","ivw","ivw_se","mode","mode_se","TSLS","TSLS_se","liml","liml_se","bnmr","bnmr_se","bnmr_sd","L","U","R")
#write_csv(dfre,"mr_re_50_gwas_b.csv")
dfre <- dfre %>% bind_rows(dfre1)
dfre

library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV3"))
resL <- parLapply(cl, seq(1,k), FUN, IV3,"LDprune")
stopCluster(cl)
# resL <- lapply(seq(1,k), FUN, IV3,"LDprune")
dfre1 <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre1) <- c("method","outcome","median","median_se","egger","egger_se","ivw","ivw_se","mode","mode_se","TSLS","TSLS_se","liml","liml_se","bnmr","bnmr_se","bnmr_sd","L","U","R")
#write_csv(dfre,"mr_re_50_ldprune_b.csv")
dfre <- dfre %>% bind_rows(dfre1)
dfre

library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV4"))
resL <- parLapply(cl, seq(1,k), FUN, IV4,"PCA")
stopCluster(cl)
# resL <- lapply(seq(1,k), FUN, IV4,"PCA")
dfre1 <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre1) <- c("method","outcome","median","median_se","egger","egger_se","ivw","ivw_se","mode","mode_se","TSLS","TSLS_se","liml","liml_se","bnmr","bnmr_se","bnmr_sd","L","U","R")
#write_csv(dfre,"mr_re_50_ldprune_b.csv")
dfre <- dfre %>% bind_rows(dfre1)
dfre

library(parallel)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV5"))
resL <- parLapply(cl, seq(1,k), FUN, IV4,"LASSO")
stopCluster(cl)
# resL <- lapply(seq(1,k), FUN, IV4,"LASSO")
dfre1 <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre1) <- c("method","outcome","median","median_se","egger","egger_se","ivw","ivw_se","mode","mode_se","TSLS","TSLS_se","liml","liml_se","bnmr","bnmr_se","bnmr_sd","L","U","R")
#write_csv(dfre,"mr_re_50_ldprune_b.csv")
dfre <- dfre %>% bind_rows(dfre1)
dfre

dfre %>% write_csv("mr_re_type8_7.csv")
