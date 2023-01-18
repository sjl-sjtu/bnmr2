bn <- function(df,snp,exposureName,bn_method="hc",repeats=1000,selectNum=50,nsam=1000,psam=100,sample_replace=TRUE){
  library("bnlearn")
  library("plyr")
  library("dplyr")
  library("parallel")
  learnBN <- function(iter,df,snp,exposureName,nsam,psam,bn_method,sample_replace){
    n <- nrow(df)
    if(sample_replace==TRUE){
      iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    }else{
      if(nsam>n){
        stop("subsample size is larger than the original sample size")
      }else{
        iSam <- sample(seq(1,n),size = nsam,replace=FALSE)
      }
    }
    if(psam>length(snp)){
      stop("subsample features is more than the number of features")
    }else{
      jSam <- sample(snp,size = psam, replace = FALSE)
    }
    dfSam <- df[iSam,c(jSam,exposureName)]
    rmFlag <- 0
    if(bn_method=="pc.stable"){
      model <- pc.stable(dfSam,undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="gs"){
      model <- gs(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="iamb"){
      model <- iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="fast.iamb"){
      model <- fast.iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="inter.iamb"){
      model <- inter.iamb(dfSam, undirected = TRUE)
      rmFlag = 1
    }else if(bn_method=="iamb.fdr"){
      cores1 <- detectCores(logical = FALSE)
      cl1 <- makeCluster(cores1)
      model <- iamb.fdr(dfSam, cl1, undirected = TRUE)
      rmFlag = 1
      stopCluster(cl1)
    }else if(bn_method=="hc"){
      model <- hc(dfSam)
    }else if(bn_method=="tabu"){
      model <- tabu(dfSam)
    }else if(bn_method=="mmhc"){
      model <- mmhc(dfSam)
    }else if(bn_method=="rsmax2"){
      model <- rsmax2(dfSam)
    }else if(bn_method=="h2pc"){
      model <- h2pc(dfSam)
    }else if(bn_method=="mmpc"){
      model <- mmpc(dfSam, undirected = TRUE)
    }else if(bn_method=="si.hiton.pc"){
      model <- si.hiton.pc(dfSam, undirected = TRUE)
    }else if(bn_method=="hpc"){
      model <- hpc(dfSam, undirected = TRUE)
    }else if(bn_method=="chow.liu"){
      model <- chow.liu(dfSam)
    }else if(bn_method=="aracne"){
      model <- aracne(dfSam)
    }else{
      stop("no this bn learning method")
    }
    dfarc <- data.frame(model$arcs)
    if(rmFlag==1){
      dfarc <- rmBidire(dfarc)
    }
    return(dfarc)
  }
  
  rmBidire <- function(df){
    df <- arrange(df,from,to)
    for(i in 1:nrow(df)){
      p <- which(df$from==df$from[i]|df$from==df$to[i])
      for(j in setdiff(p,i)){
        if(all(df[j,]%in%df[i,])){
          df <- df[-j,]
        }
      }
    }
    return(df)
  }
  
  getscore <- function(dfre,exposureName,snp,repeats){
    #exposureName is a str, snp is a vector of str.
    calc <- function(sn) {
      count1 <- dfre[which(dfre$from==sn&dfre$to==exposureName),"count"]
      count2 <- dfre[which(dfre$from==exposureName&dfre$to==sn),"count"]
      if(length(count1)==0){
        count1 <- 0
      }
      if(length(count2)==0){
        count2 <- 0
      }
      return((count1+count2)/repeats)
    }
    score <- sapply(snp, calc)
    dfscore <- data.frame(snp=snp,score=score)
    dfscore$snp <- as.character(dfscore$snp)
    return(dfscore)
  }
  
  BNbootstrap <- function(df,snp,exposureName,repeats,nsam,psam,bn_method,sample_replace){
    cores <- detectCores(logical = FALSE)
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {library("bnlearn")
      library("plyr")
      library("dplyr")
      library("parallel")
    })
    clusterExport(cl,deparse(substitute(learnBN)),envir=environment())
    arcsL <- parLapply(cl,seq(1,repeats),learnBN,df,snp,exposureName,nsam,psam,bn_method,sample_replace)
    stopCluster(cl)
    # arcsL <- replicate(repeats,learnBN(df,nsam,bn_method),simplify = FALSE)
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }
  
  df <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]
  dfre <- BNbootstrap(df,snp,exposureName,repeats,nsam,psam,bn_method,sample_replace)
  
  dfscore <- getscore(dfre,exposureName,snp,repeats)
  dfscore <- arrange(dfscore,desc(score))
  selectsnp <- dfscore[1:selectNum,"snp"]
  
  re <- list(IV=selectsnp,score=dfscore)
  return(re)
}

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
  while (TRUE %in% (qAX>=1)) {
    pi2 <- runif(m,0.05,0.95)
    r2 <- runif(m,0.2,0.8)
    qAX <- map2_dbl(pi2, r2, function(pi2,r2) 0.5+0.5*sqrt(r2)*sqrt(pi1*(1-pi1)/pi2*(1-pi2)))
  }
  G2 <- sapply(qAX, function(q) sapply(seq(1,n),linkGeno,G1,q))
}

n <- 5000
p <- 20
m <- 99
pt <- 100

getSimulate <- function(n,p,m,pt){
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
    truesnp <- sample(seq(1,ncol(df)),pt)
    mu_a <- 0.1
    alpha <- c(rnorm(pt/2,mu_a,0.01),rnorm(pt/2,-mu_a,0.01))
    alpha_0 <- 2.0
    sigma_x <- 0.25
    df$x <- alpha_0+as.matrix(df[,truesnp])%*%alpha+rnorm(n,0,sigma_x)
    
    return(list(df=df,truesnp=colnames(df)[truesnp]))
}

# nlist <- c(100,500,1000)
# plist <- c(100,200,500,1000)
# 
# re <- tibble(n=rep(nlist,each=length(plist)),p=rep(plist,times=length(nlist)),t=NA)

# for(n in c(100,500,1000,2000)){
#     for(p in c(100,200,500,1000)){
#         df <- getSimulate(n,p)
#         snp <- paste0("g",seq(1,p))
#         timestart <- Sys.time()
#         dfre <- bn_multi(df,snp,"x",bn_method="hc",selectNum=pt,repeats=100,nsam=n)
#         timeend <- Sys.time()
#         t <- difftime(timeend,timestart,units="min")
#         re[which(re$n==n&re$p==p),"t"] <- t
#         sprintf("finish: n=%d, p=%d, t=%.2f",n,p,t)
#     }
# }
# write_csv(re,"timeCompare.csv")

set.seed(43)
simd <- getSimulate(n,p,m,pt)
df <- simd$df
truesnp <- simd$truesnp
snp <- paste0("g",seq(1,p*(m+1)))


re <- tibble(ns=0,ps=0,t=0,recall=0,recall1=0,auc=0)

ns <- 2000
ps <- 30
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

ns <- 2000
ps <- 50
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)


ns <- 2000
ps <- 80
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)


ns <- 2000
ps <- 100
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

ns <- 2000
ps <- 120
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

ns <- 2000
ps <- 150
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)


ns <- 500
ps <- 100
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)


ns <- 1000
ps <- 100
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

ns <- 3000
ps <- 100
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

ns <- 5000
ps <- 100
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=1000,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re <- re%>%add_row(ns=ns,ps=ps,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,t,recal,recal1,AUC)

re %>% write_csv("compareTime_1.csv")

#repeats
re2 <- tibble(r=0,t=0,recall=0,recall1=0,auc=0)

ns <- 2000
ps <- 100

r <- 500
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

r <- 1000
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)


r <- 1500
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

r <- 2000
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)


r <- 3000
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)

r <- 5000
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re2 <- re2%>%add_row(r=r,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",ns,ps,r,t,recal,recal1,AUC)


re2 %>% write_csv("compareRepeats_1.csv")

# methods
ns <- 1000
ps <- 100
r <- 1000
re3 <- tibble(method=NA,t=0,recall=0,recall1=0,auc=0)

mt <- "hc"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "tabu"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "mmhc"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "rsmax2"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "gs"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "fast.iamb"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

mt <- "pc.stable"
timestart <- Sys.time()
dfre <- bn(df,snp,"x",bn_method=mt,selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t <- difftime(timeend,timestart,units="min")
dfs <- dfre$score %>% mutate(label=ifelse(snp%in%truesnp,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)
#re[which(re$n==n&re$p==p),"t"] <- t
re3 <- re3%>%add_row(method=mt,t=as.numeric(t),recall=recal,recall1=recal1,auc=as.numeric(AUC))
sprintf("finish: mt=%s, ns=%d, ps=%d, r=%d, t=%.2f, recall=%.2f, recall-1=%.2f, auc=%.2f",mt,ns,ps,r,t,recal,recal1,AUC)

re3 %>% write_csv("compareMethods_1.csv")

# write_csv(re,"timeCompare.csv")