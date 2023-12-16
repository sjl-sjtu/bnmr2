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
library(parallel)
library(tidyverse)
# library(AER)
# library(ivmodel)
# library(sisVIVE)
# library(MendelianRandomization)
# library(cause)
# library(R2BGLiMS)
# source("../JMCode.R")
# source("../CIIV_Functions.R")
source("bnmr.R")

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

# set.seed(0)
df <- read_csv("flexibleIV_n.csv")
df_iv <- read_csv("snp_iv.csv")
dfgwas <- df_iv
truesnp <- df_iv %>% filter(label==1) %>%pull(snp)
snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)
# dfs <- read_csv("score.csv")
# IV <- dfs %>% arrange(desc(score)) %>% slice(1:20) %>% pull(snp)
# IV2 <- dfgwas %>% arrange(p) %>% slice(1:20) %>% pull(snp)

FUN2 <- function(r,i){
  # set.seed(0)
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

  
  #bnmr
  re <- bnmr(df,snp1,"X",paste0("Y",i),nsam=1000,psam=100,repeats=1000,prior="horseshoe",n.iter=2000,n.chain=4, init=0.5,selectNum=20)
  bnmr_est <- re$mean-beta
  bnmr_se <- re$se
  L <- re$lower
  U <- re$upper
  coverage <- as.numeric((L<=beta)&(U>=beta))

  # re <- mr(df,IV2,"X",paste0("Y",i),prior="horseshoe",n.iter=3000,n.chain=4, init=0.5)
  # bnmr_est2 <- re$mean-beta
  # bnmr_se2 <- re$se
  # L2 <- re$lower
  # U2 <- re$upper
  # coverage2 <- as.numeric((L2<=beta)&(U2>=beta))

  re <- c(i,r,bnmr_est,bnmr_se,coverage) #
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

saveRDS(re,"sim_supp_bnmr_horseshoe.rds")

re <- as.matrix(re)
colnames(re) <- c("outcome","repeats",
                  "bnmr_est","bnmr_se","bnmr_cr")

re %>% as_tibble() %>% write_csv("mr_re_type_supp_bnmr_horseshoe6.csv")
stopImplicitCluster()


# dfre <- data.frame(t(matrix(unlist(resL),ncol=length(items))))

