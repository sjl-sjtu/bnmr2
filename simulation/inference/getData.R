library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("../bn.R")

set.seed(23)
library(MendelianRandomization)

linkage <- function(r2,pi1,pi2){
  #X is linked with A,p1=P(X|A),p2=P(X|a)
  q <- 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2))
  return(q)
}

linkGeno <- function(i,G1,q){
  a <- G1[i]
  if(a==0){
    x <- sample(c(0,1,2),1,prob=c(q^2,2*q*(1-q),(1-q)^2))
  }else if(a==1){
    x <- sample(c(0,1,2),1,prob=c(q*(1-q),q^2+(1-q)^2,q*(1-q)))
  }else if(a==2){
    x <- sample(c(0,1,2),1,prob=c((1-q)^2,2*q*(1-q),q^2))
  }
  return(x)
}

simLD <- function(j,g,freq,m){
  library(purrr)
  G1 <- g[,j]
  pi1 <- freq[j]
  pi2 <- runif(m,0.05,0.95)
  r2 <- runif(m,0.2,0.8)
  qAX <- map2_dbl(pi2, r2, function(pi2,r2) 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2)))
  while (TRUE %in% (qAX>=1)) {
    pi2 <- runif(m,0.05,0.95)
    r2 <- runif(m,0.2,0.8)
    qAX <- map2_dbl(pi2, r2, function(pi2,r2) 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2)))
  }
  G2 <- sapply(qAX, function(q) sapply(seq(1,n),linkGeno,G1,q))
}

n <- 4000
p <- 50
m <- 199
pt <- 300

freq <- runif(p,0.05,0.95)
g <- sapply(freq,function(pi) sample(c(0,1,2),n,prob = c(pi^2,2*pi*(1-pi),(1-pi)^2),replace = TRUE))
g2 <- map(seq(1,p),simLD,g,freq,m) 
g2 <- do.call(cbind,g2)
df <- cbind(g,g2)
df <- df[,sample(ncol(df))]
colnames(df) <- paste0("g",seq(1,p*(m+1)))
df <- as_tibble(df)   
snp <- paste0("g",seq(1,p*(m+1)))

U <- rnorm(n,0,1)
df$U <- U

truesnpInd <- sample(seq(1,ncol(df)),pt)
mu_a <- 0.1
alpha <- mu_a+abs(rnorm(pt,0,0.05))
alpha_0 <- 10
sigma_x <- 1 #
delta_x <- 0.1
epsilon_x <- rnorm(n,0,sigma_x)
df$X <- as.numeric(alpha_0+as.matrix(df[,truesnpInd])%*%alpha+delta_x*U+epsilon_x)

truesnp <- snp[truesnpInd]
write_delim(tibble(num=truesnpInd,truesnp=truesnp),"truesnp.txt")
df %>% write_csv("fullsnp.csv")

gwas <- function(g,x,df){
    lm_sum <- summary(lm(get(x)~get(g),data=df))
    p <- lm_sum$coefficients[2,4]
    fstat <- lm_sum$fstatistic[1]
    return(tibble(snp=g,p=p,fstat=fstat))
}
library(purrr)
dfgwas <- map_dfr(snp,gwas,"X",df)
dfgwas %>% write_csv("dfgwas.csv")

snp1 <- dfgwas%>%filter(p<1e-3)%>%pull(snp)
print(length(snp1))
df1 <- df%>%dplyr::select(all_of(snp1),X,U) %>% write_csv("snp_filter.csv")

sum(truesnp%in%snp1)

ns <- 4000
ps <- 120
r <- 1000
pt <- 40
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: X, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("score.csv")

IV <- dfs %>% arrange(desc(score)) %>% slice(1:40) %>% pull(snp)

IV2 <- dfgwas %>% arrange(p) %>% slice(1:40) %>% pull(snp)

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
IV3 <- stepPrune(df,dfgwas,1e-3,0.03)
print(length(IV3))

# PCA
library(stats)
pc <- prcomp(df[,snp1],scale. = T)
K <- which(cumsum(pc$sdev^2/sum(pc$sdev^2))>0.4)[1] # K is number of principal components to include in analysis
plot(1:length(pc$sdev),cumsum(pc$sdev^2/sum(pc$sdev^2)))
print(K)
K <- 40
df <- df %>% bind_cols(pc$x[,1:K])
IV4 <- paste0("PC",seq(1,K))

# penalized
library(glmnet)
f <- cv.glmnet(x=as.matrix(df[,snp1]),y=df$X,type.measure='mse',nfolds=10,pmax=70)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IV5 <- setdiff(variables, '(Intercept)')
print(length(IV5))

df_iv <- dfs %>% left_join(dfgwas,by="snp") %>%
  mutate(IV = ifelse(snp%in%IV,1,0),IV2 = ifelse(snp%in%IV2,1,0),
         IV3 = ifelse(snp%in%IV3,1,0),IV4 = ifelse(snp%in%IV4,1,0),
         IV5 = ifelse(snp%in%IV5,1,0)) %>% as_tibble() %>% write_csv("snp_iv.csv")

set.seed(1896)
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

df %>% write_csv("flexibleIV.csv")

####
IV <- df_iv%>%dplyr::arrange(desc(score))%>%slice_head(n=50)%>%pull(snp)

n <- nrow(df)
# balanced
sigma_y <- 1 #
delta_y <- 0.1
epsilon_y <- rnorm(n,0,sigma_y)
beta0 <- 5
beta <- 2.0

snp <- grep("^g",colnames(df),value = TRUE)
snp2 <- dfgwas%>%filter(p<1e-1)%>%pull(snp)
pm <- 100
snpy <- sample(reduce(list(snp,truesnp,snp2),setdiff),pm)
gamma1 <- rnorm(pm,0.1,0.05)

pk <- 50
snpc <- sample(setdiff(snp1,IV),pk)
gamma2 <- rnorm(pk,0,0.05)

#
df$Y1 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y2 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

# 
rho <- 0.2 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y3 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y4 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.4 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y5 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y6 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.6 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y7 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y8 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.8 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y9 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y10 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

# directional
gamma2 <- rnorm(pk,0.05,0.05)
# 
rho <- 0.2 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.05,0.05)
pleiosnp <- sample(IV,pt)
df$Y11 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y12 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.4 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.05,0.05)
pleiosnp <- sample(IV,pt)
df$Y13 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y14 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.6 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.05,0.05)
pleiosnp <- sample(IV,pt)
df$Y15 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y16 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)


rho <- 0.8 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.05,0.05)
pleiosnp <- sample(IV,pt)
df$Y17 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y18 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)


df %>% write_csv("fixedIV_50_0.05.csv")

##########

#### 20 ivs
df_iv <- read_csv("snp_iv.csv")
IV <- df_iv%>%dplyr::arrange(desc(score))%>%slice_head(n=20)%>%pull(snp)

n <- nrow(df)
# balanced
sigma_y <- 1 #
delta_y <- 0.1
epsilon_y <- rnorm(n,0,sigma_y)
beta0 <- 5
beta <- 2.0

snp <- grep("^g",colnames(df),value = TRUE)
snp2 <- dfgwas%>%filter(p<1e-1)%>%pull(snp)
pm <- 100
snpy <- sample(reduce(list(snp,truesnp,snp2),setdiff),pm)
gamma1 <- rnorm(pm,0.1,0.05)

pk <- 50
snpc <- sample(setdiff(snp1,IV),pk)
gamma2 <- rnorm(pk,0,0.05)

#
df$Y1 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y2 <- as.numeric(beta0+beta*df$X+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

# 
rho <- 0.2 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y3 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y4 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.4 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y5 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y6 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.6 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y7 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y8 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.8 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0,0.05)
pleiosnp <- sample(IV,pt)
df$Y9 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y10 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

# directional
# 

# pk <- 50
# snpc <- sample(setdiff(snp1,IV),pk)
gamma2 <- rnorm(pk,0.1,0.05)


rho <- 0.2 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.1,0.05)
pleiosnp <- sample(IV,pt)
df$Y11 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y12 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.4 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.1,0.05)
pleiosnp <- sample(IV,pt)
df$Y13 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y14 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)

rho <- 0.6 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.1,0.05)
pleiosnp <- sample(IV,pt)
df$Y15 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y16 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)


rho <- 0.8 #
pt <- round(rho*length(IV)) #
gamma <- rnorm(pt,0.1,0.05)
pleiosnp <- sample(IV,pt)
df$Y17 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1)

df$Y18 <- as.numeric(beta0+beta*df$X+as.matrix(df[,pleiosnp])%*%gamma+delta_y*df$U+epsilon_y+as.matrix(df[,snpy])%*%gamma1+as.matrix(df[,snpc])%*%gamma2)


df %>% write_csv("fixedIV_20_0.1.csv")
