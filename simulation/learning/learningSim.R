source("bn.R")

library(tidyverse)
library("bnlearn")
library("plyr")
library("dplyr")
library("parallel")
library(pROC)

linkage <- function(r2,pi1,pi2){
  #X is linked with A,p1=P(X|A),p2=P(X|a)
  # model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
  #                             F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
  # q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rqu,pai1,pai2))$root)
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
  r2 <- runif(m,0.05,0.95)
  qAX <- map2_dbl(pi2, r2, function(pi2,r2) 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2)))
  qAX[which(qAX>1)] <- qAX[which(qAX>1)]-1
  qAX[which(qAX==1)] <- qAX[which(qAX==1)]/2
  # while (TRUE %in% (qAX>=1)) {
  #   pi2 <- runif(m,0.05,0.95)
  #   r2 <- runif(m,0.01,0.99)
  #   qAX <- map2_dbl(pi2, r2, function(pi2,r2) 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2)))
  # }
  G2 <- sapply(qAX, function(q) sapply(seq(1,n),linkGeno,G1,q))
}

n <- 10000
p <- 50
m <- 199

library(tidyverse)
freq <- runif(p,0.05,0.95)
g <- sapply(freq,function(pi) sample(c(0,1,2),n,prob = c(pi^2,2*pi*(1-pi),(1-pi)^2),replace = TRUE))
g2 <- map(seq(1,p),simLD,g,freq,m) 
g2 <- do.call(cbind,g2)
df <- cbind(g,g2)
df <- df[,sample(ncol(df))]
# df <- sample(c(0,1,2),n*p,replace=TRUE) %>% matrix(nrow=n,ncol=p)
colnames(df) <- paste0("g",seq(1,p*(m+1)))
df <- as_tibble(df)
df %>% write_csv("dfgeno.csv")

U <- rnorm(n,0,1)
df$U <- U

pt <- 300 #
truesnp1 <- sample(seq(1,ncol(df)),pt)
mu_a <- 0.1
alpha <- mu_a+abs(rnorm(pt,0,0.05))
alpha_0 <- 10
sigma_x <- 1 #
delta_x <- 1
epsilon_x <- rnorm(n,0,sigma_x)
df$X1 <- as.numeric(alpha_0+as.matrix(df[,truesnp1])%*%alpha+delta_x*U+epsilon_x)
hsqu_x1 <- var(as.matrix(df[,truesnp1])%*%alpha)/var(df$X1)
print(hsqu_x1)    

pt <- 300 #
truesnp2 <- sample(seq(1,ncol(df)),pt)
mu_a <- 0.1
alpha <- mu_a+abs(rnorm(pt,0,0.05))
alpha_0 <- 10
sigma_x <- 2 #
delta_x <- 1
epsilon_x <- rnorm(n,0,sigma_x)
df$X2 <- as.numeric(alpha_0+as.matrix(df[,truesnp2])%*%alpha+delta_x*U+epsilon_x)
hsqu_x2 <- var(as.matrix(df[,truesnp2])%*%alpha)/var(df$X2)
print(hsqu_x2)    

df

set.seed(100)

df %>% write_csv("simulated.csv")

snp <- paste0("g",seq(1,p*(m+1)))
write_delim(tibble(num=truesnp1,truesnp=snp[truesnp1]),"truesnp_x1.txt")
write_delim(tibble(num=truesnp2,truesnp=snp[truesnp2]),"truesnp_x2.txt")

truesnp1 <- snp[truesnp1]
truesnp2 <- snp[truesnp2]

gwas <- function(g,x,df){
    lm_sum <- summary(lm(get(x)~get(g),data=df))
    p <- lm_sum$coefficients[2,4]
    fstat <- lm_sum$fstatistic[1]
    return(tibble(snp=g,p=p,fstat=fstat))
}
library(purrr)
dfgwas1 <- map_dfr(snp,gwas,"X1",df)
dfgwas1 %>% write_csv("dfgwas_x1.csv")

dfgwas2 <- map_dfr(snp,gwas,"X2",df)
dfgwas2 %>% write_csv("dfgwas_x2.csv")


####
library(data.table)
df <- fread("simulated.csv")
snp <- grep("^g",colnames(df),value = TRUE)
df <- df%>%as_tibble()%>%mutate_if(is.integer,as.numeric)
ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1.csv")
timestart <- Sys.time()
dfre <- bn(df,snp,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1.csv")


dfgwas1 <- read_csv("dfgwas_x1.csv")
dfgwas2 <- read_csv("dfgwas_x2.csv")
snp1 <- dfgwas1%>%filter(p<1e-2)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-2)%>%pull(snp)
df1 <- df%>%dplyr::select(all_of(snp1),X1,X2,U) %>% write_csv("sim_filter_x1_1e-2.csv")
df2 <- df%>%dplyr::select(any_of(snp2),X1,X2,U) %>% write_csv("sim_filter_x2_1e-2.csv")
truesnp1 <- read_delim("truesnp_x1.txt")%>%pull(truesnp)
truesnp2 <- read_delim("truesnp_x2.txt")%>%pull(truesnp)
ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1e-2.csv")
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1e-2.csv")


#1e-3
snp1 <- dfgwas1%>%filter(p<1e-3)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-3)%>%pull(snp)

df1 <- df%>%dplyr::select(all_of(snp1),X1,X2,U) %>% write_csv("sim_filter_x1.csv")
df2 <- df%>%dplyr::select(any_of(snp2),X1,X2,U) %>% write_csv("sim_filter_x2.csv")



length(snp1)
length(snp2)
sum(truesnp1%in%snp1)
sum(truesnp1%in%snp1)/300
sum(truesnp2%in%snp2)
sum(truesnp2%in%snp2)/300

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x1",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x1, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x1_score.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx1 <- df1 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx1_iv.csv")


ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x2",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x2, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x2_score.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx2 <- df %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx2_iv.csv")

#1e-4 threshold
snp1 <- dfgwas1%>%filter(p<1e-4)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-4)%>%pull(snp)

df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U) %>% write_csv("sim_filter_x1_1e-4.csv")
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U) %>% write_csv("sim_filter_x2_1e-4.csv")

length(snp1)
length(snp2)
sum(truesnp1%in%snp1)
sum(truesnp1%in%snp1)/300
sum(truesnp2%in%snp2)
sum(truesnp2%in%snp2)/300

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x1",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x1, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x1_score_1e-4.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx1 <- df1 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx1_iv_1e-4.csv")


ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x2",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x2, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x2_score_1e-4.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx2 <- df2 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx2_iv_1e-4.csv")


#1e-6 threshold
snp1 <- dfgwas1%>%filter(p<1e-6)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-6)%>%pull(snp)

df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U) %>% write_csv("sim_filter_x1_1e-6.csv")
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U) %>% write_csv("sim_filter_x2_1e-6.csv")

length(snp1)
length(snp2)
sum(truesnp1%in%snp1)
sum(truesnp1%in%snp1)/300
sum(truesnp2%in%snp2)
sum(truesnp2%in%snp2)/300

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x1",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x1, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x1_score_1e-6.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx1 <- df1 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx1_iv_1e-6.csv")


ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x2",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x2, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x2_score_1e-6.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx2 <- df2 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx2_iv_1e-6.csv")


#1e-8
snp1 <- dfgwas1%>%filter(p<1e-8)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-8)%>%pull(snp)

df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U) %>% write_csv("sim_filter_x1_1e-8.csv")
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U) %>% write_csv("sim_filter_x2_1e-8.csv")

length(snp1)
length(snp2)
sum(truesnp1%in%snp1)
sum(truesnp1%in%snp1)/300
sum(truesnp2%in%snp2)
sum(truesnp2%in%snp2)/300

ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x1",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x1, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x1_score_1e-8.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx1 <- df1 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx1_iv_1e-8.csv")


ns <- 2000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
# re <- re%>%add_row(exposure="x2",ns=ns,ps=ps,r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: x2, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

dfs %>% write_csv("x2_score_1e-8.csv")

iv <- dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(snp)
dfivx2 <- df2 %>% dplyr::select(any_of(iv),X1,X2,U) %>% write_csv("dfx2_iv_1e-8.csv")

# LD step prune
df <- fread("simulated.csv")
snp <- grep("^g",colnames(df),value = TRUE)
df <- df%>%as_tibble()%>%mutate_if(is.integer,as.numeric)
# df1 <- read_csv("sim_filter_x1.csv")
# df2 <- read_csv("sim_filter_x2.csv")
df1 <- df%>%dplyr::select(all_of(snp),X1,X2,U)
df2 <- df%>%dplyr::select(all_of(snp),X1,X2,U)
truesnp1 <- read_delim("truesnp_x1.txt") %>% pull(truesnp)
truesnp2 <- read_delim("truesnp_x2.txt") %>% pull(truesnp)
snp1 <- grep("^g",colnames(df1), value = TRUE)
snp2 <- grep("^g",colnames(df2), value = TRUE)
# length(snp1)
# length(snp2)
dfgwas1 <- read_csv("dfgwas_x1.csv")
dfgwas2 <- read_csv("dfgwas_x2.csv")

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
      return(cor(df[,i],df[,j])^2)  # Rogers and Huff (2008)
    }
    LDs <- sapply(snps,calcLD,j,df) 
    snps <- snps[which(LDs<cutoff)]
    # snps <- setdiff(snps,names(LDs[which(LDs>=cutoff)]))
  }
  return(candidate)
}

timestart <- Sys.time()
prunex1_1 <- stepPrune(df1,dfgwas1,1e-8,0.001)
print(prunex1_1)
prunex1_2 <- stepPrune(df1,dfgwas1,1e-8,0.01)
print(prunex1_2)
prunex1_3 <- stepPrune(df1,dfgwas1,1e-8,0.1)
print(prunex1_3)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
len1 <- c(length(prunex1_1),length(prunex1_2),length(prunex1_3))
rec1 <- c(sum(prunex1_1%in%truesnp1)/length(prunex1_1),
          sum(prunex1_2%in%truesnp1)/length(prunex1_2),
          sum(prunex1_2%in%truesnp1)/length(prunex1_3))
print(sum(prunex1_1%in%truesnp1)/length(prunex1_1))
print(sum(prunex1_2%in%truesnp1)/length(prunex1_2))
print(sum(prunex1_3%in%truesnp1)/length(prunex1_3))

timestart <- Sys.time()
prunex2_1 <- stepPrune(df2,dfgwas2,1e-8,0.001)
print(prunex2_1)
prunex2_2 <- stepPrune(df2,dfgwas2,1e-8,0.01)
print(prunex2_2)
prunex2_3 <- stepPrune(df2,dfgwas2,1e-8,0.1)
print(prunex2_3)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")

len2 <- c(length(prunex2_1),length(prunex2_2),length(prunex2_3))
rec2 <- c(sum(prunex2_1%in%truesnp2)/length(prunex2_1),
          sum(prunex2_2%in%truesnp2)/length(prunex2_2),
          sum(prunex2_3%in%truesnp2)/length(prunex2_3))
print(sum(prunex2_1%in%truesnp2)/length(prunex2_1))
print(sum(prunex2_2%in%truesnp2)/length(prunex2_2))
print(sum(prunex2_3%in%truesnp2)/length(prunex2_3))

ldre1 <- tibble(pthres=1e-8,cut=c(0.001,0.01,0.1),len=len1,rec=rec1,t=t1/3)
ldre2 <- tibble(pthres=1e-8,cut=c(0.001,0.01,0.1),len=len2,rec=rec2,t=t2/3)


timestart <- Sys.time()
prunex1_1 <- stepPrune(df1,dfgwas1,1e-6,0.001)
print(prunex1_1)
prunex1_2 <- stepPrune(df1,dfgwas1,1e-6,0.01)
print(prunex1_2)
prunex1_3 <- stepPrune(df1,dfgwas1,1e-6,0.1)
print(prunex1_3)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
len1 <- c(length(prunex1_1),length(prunex1_2),length(prunex1_3))
rec1 <- c(sum(prunex1_1%in%truesnp1)/length(prunex1_1),
          sum(prunex1_2%in%truesnp1)/length(prunex1_2),
          sum(prunex1_2%in%truesnp1)/length(prunex1_3))
print(sum(prunex1_1%in%truesnp1)/length(prunex1_1))
print(sum(prunex1_2%in%truesnp1)/length(prunex1_2))
print(sum(prunex1_3%in%truesnp1)/length(prunex1_3))

timestart <- Sys.time()
prunex2_1 <- stepPrune(df2,dfgwas2,1e-6,0.001)
print(prunex2_1)
prunex2_2 <- stepPrune(df2,dfgwas2,1e-6,0.01)
print(prunex2_2)
prunex2_3 <- stepPrune(df2,dfgwas2,1e-6,0.1)
print(prunex2_3)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")

len2 <- c(length(prunex2_1),length(prunex2_2),length(prunex2_3))
rec2 <- c(sum(prunex2_1%in%truesnp2)/length(prunex2_1),
          sum(prunex2_2%in%truesnp2)/length(prunex2_2),
          sum(prunex2_3%in%truesnp2)/length(prunex2_3))
print(sum(prunex2_1%in%truesnp2)/length(prunex2_1))
print(sum(prunex2_2%in%truesnp2)/length(prunex2_2))
print(sum(prunex2_3%in%truesnp2)/length(prunex2_3))

ldre1 <- ldre1 %>% bind_rows(tibble(pthres=1e-6,cut=c(0.001,0.01,0.1),len=len1,rec=rec1,t=t1/3))
ldre2 <- ldre2 %>% bind_rows(tibble(pthres=1e-6,cut=c(0.001,0.01,0.1),len=len2,rec=rec2,t=t2/3))

timestart <- Sys.time()
prunex1_1 <- stepPrune(df1,dfgwas1,1e-4,0.001)
print(prunex1_1)
prunex1_2 <- stepPrune(df1,dfgwas1,1e-4,0.01)
print(prunex1_2)
prunex1_3 <- stepPrune(df1,dfgwas1,1e-4,0.1)
print(prunex1_3)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
len1 <- c(length(prunex1_1),length(prunex1_2),length(prunex1_3))
rec1 <- c(sum(prunex1_1%in%truesnp1)/length(prunex1_1),
          sum(prunex1_2%in%truesnp1)/length(prunex1_2),
          sum(prunex1_2%in%truesnp1)/length(prunex1_3))
print(sum(prunex1_1%in%truesnp1)/length(prunex1_1))
print(sum(prunex1_2%in%truesnp1)/length(prunex1_2))
print(sum(prunex1_3%in%truesnp1)/length(prunex1_3))

timestart <- Sys.time()
prunex2_1 <- stepPrune(df2,dfgwas2,1e-4,0.001)
print(prunex2_1)
prunex2_2 <- stepPrune(df2,dfgwas2,1e-4,0.01)
print(prunex2_2)
prunex2_3 <- stepPrune(df2,dfgwas2,1e-4,0.1)
print(prunex2_3)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")
len2 <- c(length(prunex2_1),length(prunex2_2),length(prunex2_3))
rec2 <- c(sum(prunex2_1%in%truesnp2)/length(prunex2_1),
          sum(prunex2_2%in%truesnp2)/length(prunex2_2),
          sum(prunex2_3%in%truesnp2)/length(prunex2_3))
print(sum(prunex2_1%in%truesnp2)/length(prunex2_1))
print(sum(prunex2_2%in%truesnp2)/length(prunex2_2))
print(sum(prunex2_3%in%truesnp2)/length(prunex2_3))

ldre1 <- ldre1 %>% bind_rows(tibble(pthres=1e-4,cut=c(0.001,0.01,0.1),len=len1,rec=rec1,t=t1/3)) #%>% write_csv("idprune_x1.csv")
ldre2 <- ldre2 %>% bind_rows(tibble(pthres=1e-4,cut=c(0.001,0.01,0.1),len=len2,rec=rec2,t=t2/3)) #%>% write_csv("idprune_x2.csv")

timestart <- Sys.time()
prunex1_1 <- stepPrune(df1,dfgwas1,1e-2,0.001)
print(prunex1_1)
prunex1_2 <- stepPrune(df1,dfgwas1,1e-2,0.01)
print(prunex1_2)
prunex1_3 <- stepPrune(df1,dfgwas1,1e-2,0.1)
print(prunex1_3)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
len1 <- c(length(prunex1_1),length(prunex1_2),length(prunex1_3))
rec1 <- c(sum(prunex1_1%in%truesnp1)/length(prunex1_1),
          sum(prunex1_2%in%truesnp1)/length(prunex1_2),
          sum(prunex1_2%in%truesnp1)/length(prunex1_3))
print(sum(prunex1_1%in%truesnp1)/length(prunex1_1))
print(sum(prunex1_2%in%truesnp1)/length(prunex1_2))
print(sum(prunex1_3%in%truesnp1)/length(prunex1_3))

timestart <- Sys.time()
prunex2_1 <- stepPrune(df2,dfgwas2,1e-2,0.001)
print(prunex2_1)
prunex2_2 <- stepPrune(df2,dfgwas2,1e-2,0.01)
print(prunex2_2)
prunex2_3 <- stepPrune(df2,dfgwas2,1e-2,0.1)
print(prunex2_3)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")
len2 <- c(length(prunex2_1),length(prunex2_2),length(prunex2_3))
rec2 <- c(sum(prunex2_1%in%truesnp2)/length(prunex2_1),
          sum(prunex2_2%in%truesnp2)/length(prunex2_2),
          sum(prunex2_3%in%truesnp2)/length(prunex2_3))
print(sum(prunex2_1%in%truesnp2)/length(prunex2_1))
print(sum(prunex2_2%in%truesnp2)/length(prunex2_2))
print(sum(prunex2_3%in%truesnp2)/length(prunex2_3))

ldre1 <- ldre1 %>% bind_rows(tibble(pthres=1e-2,cut=c(0.001,0.01,0.1),len=len1,rec=rec1,t=t1/3)) 
ldre2 <- ldre2 %>% bind_rows(tibble(pthres=1e-2,cut=c(0.001,0.01,0.1),len=len2,rec=rec2,t=t2/3)) 

ldre1 %>% write_csv("idprune_x1.csv")
ldre2 %>% write_csv("idprune_x2.csv")


#penalize

#1e-4
snp1 <- dfgwas1%>%filter(p<1e-4)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-4)%>%pull(snp)
library(glmnet)
f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_3 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_3 <- setdiff(variables, '(Intercept)')

len1 <- c(length(IVx1_1),length(IVx1_2),length(IVx1_3))
rec1 <- c(sum(IVx1_1%in%truesnp1)/length(IVx1_1),
          sum(IVx1_2%in%truesnp1)/length(IVx1_2),
          sum(IVx1_3%in%truesnp1)/length(IVx1_3))

len2 <- c(length(IVx2_1),length(IVx2_2),length(IVx2_3))
rec2 <- c(sum(IVx2_1%in%truesnp2)/length(IVx2_1),
          sum(IVx2_2%in%truesnp2)/length(IVx2_2),
          sum(IVx2_3%in%truesnp2)/length(IVx2_3))

ldre1 <- tibble(pthres=1e-4,alpha=c(1,0.5,0),len=len1,rec=rec1) %>% write_csv("penalized_x1.csv")
ldre2 <- tibble(pthres=1e-4,alpha=c(1,0.5,0),len=len2,rec=rec2) %>% write_csv("penalized_x2.csv")


#1e-6
snp1 <- dfgwas1%>%filter(p<1e-6)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-6)%>%pull(snp)
library(glmnet)
f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_3 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_3 <- setdiff(variables, '(Intercept)')

len1 <- c(length(IVx1_1),length(IVx1_2),length(IVx1_3))
rec1 <- c(sum(IVx1_1%in%truesnp1)/length(IVx1_1),
          sum(IVx1_2%in%truesnp1)/length(IVx1_2),
          sum(IVx1_3%in%truesnp1)/length(IVx1_3))

len2 <- c(length(IVx2_1),length(IVx2_2),length(IVx2_3))
rec2 <- c(sum(IVx2_1%in%truesnp2)/length(IVx2_1),
          sum(IVx2_2%in%truesnp2)/length(IVx2_2),
          sum(IVx2_3%in%truesnp2)/length(IVx2_3))

ldre1 <- tibble(pthres=1e-6,alpha=c(1,0.5,0),len=len1,rec=rec1) %>% write_csv("penalized_x1.csv",append = TRUE)
ldre2 <- tibble(pthres=1e-6,alpha=c(1,0.5,0),len=len2,rec=rec2) %>% write_csv("penalized_x2.csv",append = TRUE)

#1e-8
snp1 <- dfgwas1%>%filter(p<1e-8)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-8)%>%pull(snp)
library(glmnet)
f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_3 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_3 <- setdiff(variables, '(Intercept)')

len1 <- c(length(IVx1_1),length(IVx1_2),length(IVx1_3))
rec1 <- c(sum(IVx1_1%in%truesnp1)/length(IVx1_1),
          sum(IVx1_2%in%truesnp1)/length(IVx1_2),
          sum(IVx1_3%in%truesnp1)/length(IVx1_3))

len2 <- c(length(IVx2_1),length(IVx2_2),length(IVx2_3))
rec2 <- c(sum(IVx2_1%in%truesnp2)/length(IVx2_1),
          sum(IVx2_2%in%truesnp2)/length(IVx2_2),
          sum(IVx2_3%in%truesnp2)/length(IVx2_3))

ldre1 <- tibble(pthres=1e-8,alpha=c(1,0.5,0),len=len1,rec=rec1) %>% write_csv("penalized_x1.csv",append = TRUE)
ldre2 <- tibble(pthres=1e-8,alpha=c(1,0.5,0),len=len2,rec=rec2) %>% write_csv("penalized_x2.csv",append = TRUE)

#1e-2
snp1 <- dfgwas1%>%filter(p<1e-2)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-2)%>%pull(snp)
library(glmnet)
f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df1[,snp1]),y=df1$X1,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx1_3 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=1,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_1 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0.5,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_2 <- setdiff(variables, '(Intercept)')

f <- cv.glmnet(x=as.matrix(df2[,snp2]),y=df2$X2,alpha=0,type.measure='mse',nfolds=10)
c <- coef(f,s='lambda.1se',exact=TRUE)
inds <- which(c!=0)
variables <- row.names(c)[inds]
IVx2_3 <- setdiff(variables, '(Intercept)')

len1 <- c(length(IVx1_1),length(IVx1_2),length(IVx1_3))
rec1 <- c(sum(IVx1_1%in%truesnp1)/length(IVx1_1),
          sum(IVx1_2%in%truesnp1)/length(IVx1_2),
          sum(IVx1_3%in%truesnp1)/length(IVx1_3))

len2 <- c(length(IVx2_1),length(IVx2_2),length(IVx2_3))
rec2 <- c(sum(IVx2_1%in%truesnp2)/length(IVx2_1),
          sum(IVx2_2%in%truesnp2)/length(IVx2_2),
          sum(IVx2_3%in%truesnp2)/length(IVx2_3))

ldre1 <- tibble(pthres=1e-2,alpha=c(1,0.5,0),len=len1,rec=rec1) %>% write_csv("penalized_x1.csv",append = TRUE)
ldre2 <- tibble(pthres=1e-2,alpha=c(1,0.5,0),len=len2,rec=rec2) %>% write_csv("penalized_x2.csv",append = TRUE)


############plot correlation with F
library(data.table)
df <- fread("simulated.csv")
snp <- grep("^g",colnames(df),value = TRUE)
df <- df%>%as_tibble()%>%mutate_if(is.integer,as.numeric)

truesnp1 <- read_delim("truesnp_x1.txt")%>%pull(truesnp)
truesnp2 <- read_delim("truesnp_x2.txt")%>%pull(truesnp)
ns <- 5000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df,snp,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1_5000.csv")
timestart <- Sys.time()
dfre <- bn(df,snp,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1_5000.csv")


dfgwas1 <- read_csv("dfgwas_x1.csv")
dfgwas2 <- read_csv("dfgwas_x2.csv")

snp1 <- dfgwas1%>%filter(p<1e-2)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-2)%>%pull(snp)
df1 <- df%>%dplyr::select(all_of(snp1),X1,X2,U) 
df2 <- df%>%dplyr::select(any_of(snp2),X1,X2,U) 

ns <- 5000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1e-2_5000.csv")
timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1e-2_5000.csv")


snp1 <- dfgwas1%>%filter(p<1e-4)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-4)%>%pull(snp)
df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U)
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U)
ns <- 5000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1e-4_5000.csv")

timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1e-4_5000.csv")


snp1 <- dfgwas1%>%filter(p<1e-6)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-6)%>%pull(snp)
df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U)
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U)
ns <- 5000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1e-6_5000.csv")

timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1e-6_5000.csv")

snp1 <- dfgwas1%>%filter(p<1e-8)%>%pull(snp)
snp2 <- dfgwas2%>%filter(p<1e-8)%>%pull(snp)
df1 <- df1%>%dplyr::select(all_of(snp1),X1,X2,U)
df2 <- df2%>%dplyr::select(any_of(snp2),X1,X2,U)
ns <- 5000
ps <- 150
r <- 5000
pt <- 50
timestart <- Sys.time()
dfre <- bn(df1,snp1,"X1",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
dfs %>% write_csv("x1_score_1e-8_5000.csv")

timestart <- Sys.time()
dfre <- bn(df2,snp2,"X2",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp2,1,0))
dfs %>% write_csv("x2_score_1e-8_5000.csv")
