library(tidyverse)
library(data.table)
set.seed(0)
dfgeno <- fread("genedata/merged.raw")
n <- nrow(dfgeno)
colnames(dfgeno) <- gsub("_[ATCG]*","",colnames(dfgeno))
source("../bn.R")
#
# # print("1")
# #
# # gwas <- function(g,v,df){
# #   lm_sum <- summary(RcppEigen::fastLmPure(df$get(v),df$get(g)))
# #   p <- lm_sum$coefficients[2,4]
# #   fstat <- lm_sum$fstatistic[1]
# #   return(tibble(snp=g,p=p,fstat=fstat))
# # }
# # library(purrr)
snps <- colnames(dfgeno)[-(1:6)]
#
p <- 10000
snp <- sample(snps,p)

write_delim(tibble(snp=snp),"snp.txt",col_names=F)

# df <- (dfgeno%>%as_tibble)[,snp]
dfgeno <- dfgeno%>%as_tibble()
dfgeno <- bind_cols(dfgeno[,1:6],dfgeno[,snp])
pt <- 500
z_d <- sample(7:ncol(dfgeno),pt)
truesnp <- colnames(dfgeno)[z_d]
write_delim(tibble(truesnp=truesnp),"truesnp.txt")



n <- 5000
rownameID <- sample(seq(1,nrow(dfgeno)),n)
df <- dfgeno[rownameID,]
u <- rnorm(n,0,1)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)
x <- (alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x) %>% as.numeric()
df$x <- x
df %>% write_csv("data_5000.csv")
tibble(FID=df$FID,IID=df$IID,x=df$x)%>%write_delim("pheno_5000.txt")

n <- 10000
rownameID <- sample(seq(1,nrow(dfgeno)),n)
df <- dfgeno[rownameID,]
u <- rnorm(n,0,1)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)
x <- (alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x) %>% as.numeric()
df$x <- x
df %>% write_csv("data_10000.csv")
tibble(FID=df$FID,IID=df$IID,x=df$x)%>%write_delim("pheno_10000.txt")

n <- 50000
rownameID <- sample(seq(1,nrow(dfgeno)),n)
df <- dfgeno[rownameID,]
u <- rnorm(n,0,1)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)
x <- (alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x) %>% as.numeric()
df$x <- x
df %>% write_csv("data_50000.csv")
tibble(FID=df$FID,IID=df$IID,x=df$x)%>%write_delim("pheno_50000.txt")

n <- 2000
rownameID <- sample(seq(1,nrow(dfgeno)),n)
df <- dfgeno[rownameID,]
u <- rnorm(n,0,1)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)
x <- (alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x) %>% as.numeric()
df$x <- x
df %>% write_csv("data_2000.csv")
tibble(FID=df$FID,IID=df$IID,x=df$x)%>%write_delim("pheno_2000.txt")

# PLINK GWAS

truesnp <- read_delim("truesnp.txt",delim=" ") %>% pull(truesnp)


dfgwas <- read_delim("genedata/gwas_2000.x.glm.linear")
snp1 <- dfgwas%>%filter(P<1e-4)%>%pull(ID)
df <- read_csv("data_2000.csv")
df <- df[,c(snp1,"x")] %>% mutate_if(is.integer,as.numeric)

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)
dfs %>% write_csv("x_score_n2000_2.csv")


dfgwas <- read_delim("genedata/gwas_5000.x.glm.linear")
snp1 <- dfgwas%>%filter(P<1e-4)%>%pull(ID)
df <- read_csv("data_5000.csv")
df <- df[,c(snp1,"x")] %>% mutate_if(is.integer,as.numeric)

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)
dfs %>% write_csv("x_score_n5000_2.csv")



dfgwas <- read_delim("genedata/gwas_10000.x.glm.linear")
snp1 <- dfgwas%>%filter(P<1e-6)%>%pull(ID)
df <- read_csv("data_10000.csv")
df <- df[,c(snp1,"x")] %>% mutate_if(is.integer,as.numeric)

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)
dfs %>% write_csv("x_score_n10000_2.csv")



dfgwas <- read_delim("genedata/gwas_50000.x.glm.linear")
snp1 <- dfgwas%>%filter(P<1e-10)%>%pull(ID)
df <- read_csv("data_50000.csv")
df <- df[,c(snp1,"x")] %>% mutate_if(is.integer,as.numeric)
truesnp <- read_delim("truesnp.txt",delim=" ") %>% pull(truesnp)

ns <- 4000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x_score_n50000_2.csv")


################
library(tidyverse)
library(data.table)
set.seed(0)
dfgeno <- fread("merged2.raw")
n <- nrow(dfgeno)
colnames(dfgeno) <- gsub("_[ATCG]*","",colnames(dfgeno))
source("../bn.R")

gwas <- function(g,v,df){
  lm_sum <- summary(lm(get(v)~get(g),data=df))
  p <- lm_sum$coefficients[2,4]
  fstat <- lm_sum$fstatistic[1]
  return(tibble(snp=g,p=p,fstat=fstat))
}
library(purrr)
snps <- colnames(dfgeno)[-(1:6)]

p <- 10000
snp <- sample(snps,p)
df <- (dfgeno%>%as_tibble)[,snp] %>% mutate_if(is.integer,as.numeric)
u <- rnorm(n,0,1)
pt <- 500
z_d <- sample(1:ncol(df),pt)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)

x <- alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x
df$x <- x

truesnp <- colnames(df)[z_d]

dfgwas <- map_dfr(snp,gwas,"x",df)
snp1 <- dfgwas%>%filter(p<1e-4)%>%pull(snp)
print(length(snp1))

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

# plot(roc_object)
dfs %>% write_csv("x_score_p10000_3.csv")

p <- 50000
snp <- sample(snps,p)
df <- (dfgeno%>%as_tibble)[,snp] %>% mutate_if(is.integer,as.numeric)
u <- rnorm(n,0,1)
pt <- 500
z_d <- sample(1:ncol(df),pt)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)

x <- alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x
df$x <- x

truesnp <- colnames(df)[z_d]


dfgwas <- map_dfr(snp,gwas,"x",df)
snp1 <- dfgwas%>%filter(p<1e-4)%>%pull(snp)

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

# plot(roc_object)
dfs %>% write_csv("x_score_p50000_3.csv")


p <- 20000
snp <- sample(snps,p)
df <- (dfgeno%>%as_tibble)[,snp] %>% mutate_if(is.integer,as.numeric)
u <- rnorm(n,0,1)
pt <- 500
z_d <- sample(1:ncol(df),pt)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)

x <- alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x
df$x <- x

truesnp <- colnames(df)[z_d]


dfgwas <- map_dfr(snp,gwas,"x",df)
snp1 <- dfgwas%>%filter(p<1e-4)%>%pull(snp)

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

# plot(roc_object)
dfs %>% write_csv("x_score_p20000_3.csv")


p <- 100000
snp <- sample(snps,p)
df <- (dfgeno%>%as_tibble)[,snp] %>% mutate_if(is.integer,as.numeric)
u <- rnorm(n,0,1)
pt <- 500
z_d <- sample(1:ncol(df),pt)
mu_a <- 0.1
alpha <- rnorm(pt,mu_a,0.05)
delta_x <- 0.1
alpha_0 <- 2.0
sigma_x <- 0.25
epsilon_x <- rnorm(n,0,sigma_x)

x <- alpha_0+as.matrix(df[,z_d])%*%alpha+delta_x*u+epsilon_x
df$x <- x

truesnp <- colnames(df)[z_d]

dfgwas <- map_dfr(snp,gwas,"x",df)
snp1 <- dfgwas%>%filter(p<1e-4)%>%pull(snp)
print(length(snp1))

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp1,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
sprintf("finish: x, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

# plot(roc_object)
dfs %>% write_csv("x_score_p100000_3.csv")
