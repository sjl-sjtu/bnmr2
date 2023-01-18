library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("../bn.R")
source("mr.R")

set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")

snp <- df_iv%>%filter(p<1e-2)%>%pull(snp)

stepPrune <- function(df,dfgwas,pthres,cutoff){
  dfgwas <- dfgwas %>% filter(p<pthres)
  snps <- dfgwas %>% pull(snp)
  candidate <- c()
  while(length(snps)>0){
    dfgwas <- dfgwas%>%filter(snp %in% snps)
    j <- dfgwas[which.min(dfgwas$p),"snp"] %>% as.character()
    candidate <- c(candidate,j)
    snps <- setdiff(snps,j)
    calcLD <- function(i,j,df){
      return(cor(df[,i],df[,j]))  # Rogers and Huff (2008)
    }
    LDs <- sapply(snps,calcLD,j,df) 
    snps <- snps[which(LDs<cutoff)]
  }
  return(candidate)
}
top_vars <- stepPrune(df,df_iv,1e-3,0.1)
prune_s <- stepPrune(df,df_iv,1e-2,0.8)



FUN2 <- function(i){
  library(tidyverse)
  
  #CAUSE
  exposure <- df$X
  outcome <- df%>%pull(paste0("Y",i))
  s <- grep("^g",colnames(df),value = TRUE)
  N <- nrow(df)
  J <- length(s)
  X <- array(exposure,dim=N)
  Y <- array(outcome,dim=N)
  Z <- as.matrix(df[,s],dim=c(N,J))
  
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
  df_summary <- tibble(snp=s,beta_hat_1=betaX,seb1=sebetaX,p1=px,beta_hat_2=betaY,seb2=sebetaY,p2=py)
  
  library(cause)
  params <- est_cause_params(df_summary, s)
  # r2_thresh = 0.01
  # pval_thresh = 1e-3
  # df_clump <- df_summary %>%
  #           rename(rsid = snp,
  #                   pval = p1) %>%
  #           ieugwasr::ld_clump(dat = .,
  #                     clump_r2 = r2_thresh,
  #                     clump_p = pval_thresh,
  #                     #plink_bin = genetics.binaRies::get_plink_binary(), 
  #                     pop = "EUR")
  # top_vars <- df_clump$rsid
  
  res <- cause(X=df_summary, variants = top_vars, param_ests = params,force=TRUE)
  rest <- res$elpd
  res_sum <- summary(res, ci_size=0.95)
  cause_est <- res_sum$tab[2,2]
  
  # sisVIVE
  library(sisVIVE)
  J <- length(snp)
  Z <- as.matrix(df[,snp],dim=c(N,J))
  sis <- cv.sisVIVE(outcome,exposure,Z,K=10)
  sisvive_est <- sis$beta
  
  source("../TSHT.R")
  source("../JMCode.R")
  source("../CIIV_Functions.R")
  
  X <- as.matrix(array(X,dim=N))
  Y <- as.matrix(array(Y,dim=N))
  Z <- as.matrix(df[,prune_s],dim=c(N,J))
  # HT
  ht <- TSHT(Y,X,Z)
  print(ht)
  ht_est <- ht$betaHat
  ht_se <- sqrt(ht$betaVarHat)
  
  # CI
  ciiv <- CIIV(Y,X,Z)
  print(ciiv)
  ci_est <- ciiv$Coefficients_CIM[1]
  cigmm_est <- ciiv$Coefficients_CIM_GMM[1]
  ci_se <- ciiv$sd_CIM[1]
  cigmm_se <- ciiv$sd_CIM_GMM[1]
  
  #JAM-MR
  library(R2BGLiMS)
  J <- length(prune_s)
  betaX <- c()
  betaY <- c()
  sebetaX <- c()
  sebetaY <- c()
  for(isnp in 1:J){
    regX <- lm(X ~ Z[,isnp])
    regY <- lm(Y ~ Z[,isnp])
    betaX[isnp] <- summary(regX)$coefficients[2,1]
    sebetaX[isnp] <- summary(regX)$coefficients[2,2]
    betaY[isnp] <- summary(regY)$coefficients[2,1]
    sebetaY[isnp] <- summary(regY)$coefficients[2,2]
  }
  # eafs <- sapply(1:J,function(j) (2*sum(Z[,j]==2)+sum(Z[,j]==1))/(2*nrow(Z)))
  # G.cor <- cor(Z)
  es <- JAMMR(betaX,sebetaX,betaY,sebetaY,N1=dim(Z)[1],G.matrix = Z,n.grid=10,iter = 10000,jam.seed=123)#eafs=eafs,G.cor = G.cor,)
  jam_est <- es$causal
  jam_se <- es$se
  
  re <- c(paste0("Y",i),cause_est,sisvive_est,ht_est,ht_se,ci_est,ci_se,cigmm_est,
          cigmm_se,jam_est,jam_se)
  return(re)
}


k <- 13
#cores <- detectCores(logical = FALSE)
#cl <- makeCluster(cores)
#clusterExport(cl,varlist=c("df","mr","snp","top_vars","prune_s"))
resL <- lapply(seq(1,k), FUN2)
#stopCluster(cl)
# resL <- lapply(seq(1,k), FUN2)
dfre <- data.frame(t(matrix(unlist(resL),ncol=k)))
colnames(dfre) <- c("outcome","cause_est",
                    "sisvive_est","ht_est","ht_se","ci_est","ci_se","cigmm_est",
                    "cigmm_se","jam_est","jam_se")
write_csv(dfre,"mr_re_type8_8.csv")
