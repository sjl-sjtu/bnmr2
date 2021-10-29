set.seed(111)
bn_multi <- function(df,snp,exposureName,bn_method="hr",cutoff=0.7,repeats=100,nsam=500){
  library("bnlearn")
  library("plyr")
  library("dplyr")
  learnBN <- function(df,nsam,bn_method){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    rmFlag <- 0
    if(bn_method=="pc.stable"){
      model <- pc.stable(dfSam)
      rmFlag = 1
    }else if(bn_method=="gs"){
      model <- gs(dfSam)
      rmFlag = 1
    }else if(bn_method=="iamb"){
      model <- iamb(dfSam)
      rmFlag = 1
    }else if(bn_method=="fast.iamb"){
      model <- fast.iamb(dfSam)
      rmFlag = 1
    }else if(bn_method=="inter.iamb"){
      model <- inter.iamb(dfSam)
      rmFlag = 1
    }else if(bn_method=="iamb.fdr"){
      model <- iamb.fdr(dfSam)
      rmFlag = 1
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
      model <- mmpc(dfSam)
    }else if(bn_method=="si.hiton.pc"){
      model <- si.hiton.pc(dfSam)
    }else if(bn_method=="hpc"){
      model <- hpc(dfSam)
    }else if(bn_method=="chow.liu"){
      model <- chow.liu(dfSam)
    }else if(bn_method=="aracne"){
      model <- aracne(dfSam)
    }else{
      return(message("no this bn learning method"))
    }
    dfarc <- data.frame(model$arcs)
    dfarc$from <- as.character(dfarc$from)
    dfarc$to <- as.character(dfarc$to)
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
  
  BNbootstrap <- function(df,repeats,nsam,bn_method){
    arcsL <- replicate(repeats,learnBN(df,nsam,bn_method),simplify = FALSE)
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }
  
  getscore <- function(dfre,exposureName,snp,repeats){
    #exposureName is a vector of str, snp is a vector of str.
    dfscoreList <- list()
    for(j in 1:length(exposureName)){
      expos <- exposureName[j]
      score <- rep(0,length(snp))
      for(i in 1:length(snp)){
        sn <- snp[i]
        count1 <- dfre[which(dfre$from==sn&dfre$to==expos),"count"]
        count2 <- dfre[which(dfre$from==expos&dfre$to==sn),"count"]
        if(length(count1)==0){
          count1 <- 0
        }
        if(length(count2)==0){
          count2 <- 0
        }
        score[i] <- (count1+count2)/repeats
      }
      dfscore <- data.frame(snp=snp,score=score)
      dfscore$snp <- as.character(dfscore$snp)
      dfscoreList[[j]] <- dfscore
    }
    return(dfscoreList)
  }
  
  df1 <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]
  
  dfre <- BNbootstrap(df1,repeats,nsam,bn_method)
  dfscore <- getscore(dfre,exposureName,snp,repeats)
  selectsnp <- list()
  for(i in 1:length(exposureName)){
    selectsnp1 <- dfscore[[i]][which(dfscore[[i]]$score>=cutoff),"snp"]
    selectsnp <- append(selectsnp,selectsnp1)
  }
  
  re <- list(IV=selectsnp,score=dfscore)
  return(re)
}

simulaMDdata <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                         thetaA,thetaB,thetaA2,thetaM,rsAX,rsBY,
                         n,FM1,SE1,FM2,SE2,sigma2E){
  library(rootSolve)
  library(bnlearn)
  GA <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiA)^2,2*paiA*(1-paiA),paiA^2))
  GB <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiB)^2,2*paiB*(1-paiB),paiB^2))
  phenotype <- function(i,G1,G2,alpha,theta1,theta2,FM,SE){
    a <- G1[i]
    b <- G2[i]
    if(a==0&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==0&b==1){
      d <- rnorm(1,FM*alpha*(1+theta2),SE)
    }else if(a==0&b==2){
      d <- rnorm(1,FM*alpha*(1+theta2)^2,SE)
    }else if(a==1&b==0){
      d <- rnorm(1,FM*alpha*(1+theta1),SE)
    }else if(a==1&b==1){
      d <- rnorm(1,FM*alpha*(1+theta1)*(1+theta2),SE)
    }else if(a==1&b==2){
      d <- rnorm(1,FM*alpha*(1+theta1)*(1+theta2)^2,SE)
    }else if(a==2&b==0){
      d <- rnorm(1,FM*alpha*(1+theta1)^2,SE)
    }else if(a==2&b==1){
      d <- rnorm(1,FM*alpha*(1+theta1)^2*(1+theta2),SE)
    }else if(a==2&b==2){
      d <- rnorm(1,FM*alpha*(1+theta1)^2*(1+theta2)^2,SE)
    }
    return(d)
  }
  
  F1 <- sapply(seq(1,n),phenotype,GA,GB,alpha,thetaA,thetaB,FM1,SE1)
  
  linkage <- function(rsqu,pai1,pai2){
    #X is linked with A,p1=P(X|A),p2=P(X|a)
    model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
                                F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
    q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rsqu,pai1,pai2))$root)
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
  
  qAX <- linkage(rsAX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(rsBY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,thetaA2,thetaM,FM2,SE2)
  
  pais <- runif(5,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(F1,F2,GA,GB,GM,GX,GY,geno)
  return(df)
}

simulaMDdata2 <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                          theta,theta2,rsAX,rsBY,
                          n,FM1,SE1,FM2,SE2,sigma2E){
  library(rootSolve)
  GA <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiA)^2,2*paiA*(1-paiA),paiA^2))
  GB <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiB)^2,2*paiB*(1-paiB),paiB^2))
  phenotype <- function(i,G1,G2,alpha,theta,FM,SE){
    a <- G1[i]
    b <- G2[i]
    if(a==0&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==0&b==1){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==0&b==2){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1&b==1){
      d <- rnorm(1,FM*alpha*(1+theta),SE)
    }else if(a==1&b==2){
      d <- rnorm(1,FM*alpha*(1+theta)^2,SE)
    }else if(a==2&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==2&b==1){
      d <- rnorm(1,FM*alpha*(1+theta)^2,SE)
    }else if(a==2&b==2){
      d <- rnorm(1,FM*alpha*(1+theta)^4,SE)
    }
    return(d)
  }
  
  F1 <- sapply(seq(1,n),phenotype,GA,GB,alpha,theta,FM1,SE1)
  
  linkage <- function(rsqu,pai1,pai2){
    #X is linked with A,p1=P(X|A),p2=P(X|a)
    model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
                                F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
    q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rsqu,pai1,pai2))$root)
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
  
  qAX <- linkage(rsAX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(rsBY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
  pais <- runif(5,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(F1,F2,GA,GB,GM,GX,GY,geno)
  return(df)
}

simulaMDdata3 <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                          theta,theta2,rsAX,rsBY,
                          n,FM1,SE1,FM2,SE2,sigma2E){
  library(rootSolve)
  GA <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiA)^2,2*paiA*(1-paiA),paiA^2))
  GB <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiB)^2,2*paiB*(1-paiB),paiB^2))
  phenotype <- function(i,G1,G2,alpha,theta,FM,SE){
    a <- G1[i]
    b <- G2[i]
    if(a==0&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==0&b==1){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==0&b==2){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1&b==1){
      d <- rnorm(1,FM*alpha*(1+theta),SE)
    }else if(a==1&b==2){
      d <- rnorm(1,FM*alpha*(1+theta),SE)
    }else if(a==2&b==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==2&b==1){
      d <- rnorm(1,FM*alpha*(1+theta),SE)
    }else if(a==2&b==2){
      d <- rnorm(1,FM*alpha*(1+theta),SE)
    }
    return(d)
  }
  
  F1 <- sapply(seq(1,n),phenotype,GA,GB,alpha,theta,FM1,SE1)
  
  linkage <- function(rsqu,pai1,pai2){
    #X is linked with A,p1=P(X|A),p2=P(X|a)
    model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
                                F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
    q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rsqu,pai1,pai2))$root)
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
  
  qAX <- linkage(rsAX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(rsBY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
  epsilon <- runif(n,0,sigma2E)
  pai <- runif(1,0.1,0.4)
  X <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
  D <- 1.2*F1+3*F2+epsilon
  
  pais <- runif(5,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(F1,F2,GA,GB,GM,GX,GY,geno)
  return(df)
}

snpnameA <- sapply(seq(1,15),paste,"GA",sep="")
snpnameB <- sapply(seq(1,15),paste,"GB",sep="")
snpnameM <- sapply(seq(1,15),paste,"GM",sep="")
snpnameX <- sapply(seq(1,15),paste,"GX",sep="")
snpnameY <- sapply(seq(1,15),paste,"GY",sep="")
phenoname1 <- sapply(seq(1,15),paste,"F1",sep="")
phenoname2 <- sapply(seq(1,15),paste,"F2",sep="")
xname1 <- sapply(seq(1,15),paste,"X1",sep="")
xname2 <- sapply(seq(1,15),paste,"X2",sep="")
xname3 <- sapply(seq(1,15),paste,"X3",sep="")
xname4 <- sapply(seq(1,15),paste,"X4",sep="")
xname5 <- sapply(seq(1,15),paste,"X5",sep="")

df <- data.frame(ID=seq(1,5000))
for(i in 1:5){
  paiA<-runif(1,0.12,0.3)
  paiB <-runif(1,0.12,0.3)
  paiM <-runif(1,0.12,0.3)
  paiX <-runif(1,0.12,0.3)
  paiY <-runif(1,0.12,0.3)
  alpha <- runif(1,0.15,0.3)
  alpha2 <- runif(1,0.15,0.3)
  thetaA <- runif(1,0.2,0.5)
  thetaB <- runif(1,0.2,0.5)
  thetaA2 <- runif(1,0.2,0.5)
  thetaM <- runif(1,0.2,0.5)
  rsAX <- runif(1,0.3,0.7)
  rsBY <- runif(1,0.3,0.7)
  n <- 5000
  FM1 <- sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                      thetaA,thetaB,thetaA2,thetaM,rsAX,rsBY,
                      n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i],xname1[i],xname2[i],
                   xname3[i],xname4[i],xname5[i])
  df <- cbind(df,dfs)
}
for(i in 6:10){
  paiA<-runif(1,0.12,0.3)
  paiB <-runif(1,0.12,0.3)
  paiM <-runif(1,0.12,0.3)
  paiX <-runif(1,0.12,0.3)
  paiY <-runif(1,0.12,0.3)
  alpha <- runif(1,0.15,0.3)
  alpha2 <- runif(1,0.15,0.3)
  theta <- runif(1,0.25,0.6)
  theta2 <- runif(1,0.25,0.6)
  rsAX <- runif(1,0.3,0.7)
  rsBY <- runif(1,0.3,0.7)
  n <- 5000
  FM1 <- sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata2(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                       theta,theta2,rsAX,rsBY,
                       n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i],xname1[i],xname2[i],
                   xname3[i],xname4[i],xname5[i])
  df <- cbind(df,dfs)
}
for(i in 11:15){
  paiA<-runif(1,0.12,0.3)
  paiB <-runif(1,0.12,0.3)
  paiM <-runif(1,0.12,0.3)
  paiX <-runif(1,0.12,0.3)
  paiY <-runif(1,0.12,0.3)
  alpha <- runif(1,0.15,0.3)
  alpha2 <- runif(1,0.15,0.3)
  theta <- runif(1,0.25,0.6)
  theta2 <- runif(1,0.25,0.6)
  rsAX <- runif(1,0.3,0.7)
  rsBY <- runif(1,0.3,0.7)
  n <- 5000
  FM1 <- sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata3(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                       theta,theta2,rsAX,rsBY,
                       n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i],xname1[i],xname2[i],
                   xname3[i],xname4[i],xname5[i])
  df <- cbind(df,dfs)
}
df$F1 <- df[,"1F1"]
df$F2 <- df[,"1F2"]
for(i in 2:15){
  c <- phenoname1[i]
  d <- phenoname2[i]
  df$F1 <- df$F1*df[,c]
  df$F2 <- df$F2*df[,d]
}
df$F1 <- abs(df$F1)^(1/15)
df$F2 <- abs(df$F2)^(1/15)
df1 <- df[,c("F1","F2",snpnameA,snpnameB,snpnameM,snpnameX,snpnameY,xname1,
             xname2,xname3,xname4,xname5)]

snp <- as.character(setdiff(colnames(df1),c("F1","F2")))
phenoName <- as.character(c("F1","F2"))

timestart <- Sys.time()
dfre1 <- bn_multi(df1,snp,phenoName,bn_method="hc",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")

timestart <- Sys.time()
dfre2 <- bn_multi(df1,snp,phenoName,bn_method="iamb",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")

timestart <- Sys.time()
dfre3 <- bn_multi(df1,snp,phenoName,bn_method="mmhc",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t3 <- difftime(timeend,timestart,units="min")
print("time")
print(c(t1,t2,t3))

dfs1_1 <- dfre1$score[[1]]
dfs1_2 <- dfre2$score[[1]]
dfs1_3 <- dfre3$score[[1]]
dfs2_1 <- dfre1$score[[2]]
dfs2_2 <- dfre2$score[[2]]
dfs2_3 <- dfre3$score[[2]]

truesnp1 <- c(snpnameA,snpnameB)
truesnp2 <- c(snpnameA,snpnameM)
lse1_1 <- c()
lsp1_1 <- c()
lse1_2 <- c()
lsp1_2 <- c()
lse1_3 <- c()
lsp1_3 <- c()
lse2_1 <- c()
lsp2_1 <- c() 
lse2_2 <- c()
lsp2_2 <- c()
lse2_3 <- c()
lsp2_3 <- c()
cutoff <- seq(0,1,0.05)
for(i in 1:length(cutoff)){
  c <- cutoff[i]
  f1snp_1 <- dfs1_1[which(dfs1_1$score>=c),"snp"]
  f1snp_2 <- dfs1_2[which(dfs1_2$score>=c),"snp"]
  f1snp_3 <- dfs1_3[which(dfs1_3$score>=c),"snp"]
  f2snp_1 <- dfs2_1[which(dfs2_1$score>=c),"snp"]
  f2snp_2 <- dfs2_2[which(dfs2_2$score>=c),"snp"]
  f2snp_3 <- dfs2_3[which(dfs2_3$score>=c),"snp"]
  tp1_1 <- sum(f1snp_1 %in% truesnp1)
  fp1_1 <- length(f1snp_1)-tp1_1
  fn1_1 <- length(setdiff(truesnp1,f1snp_1))
  tn1_1 <- length(setdiff(snp,f1snp_1))-fn1_1
  se1_1 <- tp1_1/(tp1_1+fn1_1) #敏感度
  sp1_1 <- tn1_1/(tn1_1+fp1_1) #特异度
  lse1_1[i] <- se1_1
  lsp1_1[i] <- sp1_1
  tp1_2 <- sum(f1snp_2 %in% truesnp1)
  fp1_2 <- length(f1snp_2)-tp1_2
  fn1_2 <- length(setdiff(truesnp1,f1snp_2))
  tn1_2 <- length(setdiff(snp,f1snp_2))-fn1_2
  se1_2 <- tp1_2/(tp1_2+fn1_2) #敏感度
  sp1_2 <- tn1_2/(tn1_2+fp1_2) #特异度
  lse1_2[i] <- se1_2
  lsp1_2[i] <- sp1_2
  tp1_3 <- sum(f1snp_3 %in% truesnp1)
  fp1_3 <- length(f1snp_3)-tp1_3
  fn1_3 <- length(setdiff(truesnp1,f1snp_3))
  tn1_3 <- length(setdiff(snp,f1snp_3))-fn1_3
  se1_3 <- tp1_3/(tp1_3+fn1_3) #敏感度
  sp1_3 <- tn1_3/(tn1_3+fp1_3) #特异度
  lse1_3[i] <- se1_3
  lsp1_3[i] <- sp1_3
  tp2_1 <- sum(f2snp_1 %in% truesnp2)
  fp2_1 <- length(f2snp_1)-tp2_1
  fn2_1 <- length(setdiff(truesnp2,f2snp_1))
  tn2_1 <- length(setdiff(snp,f2snp_1))-fn2_1
  se2_1 <- tp2_1/(tp2_1+fn2_1) #敏感度
  sp2_1 <- tn2_1/(tn2_1+fp2_1) #特异度
  lse2_1[i] <- se2_1
  lsp2_1[i] <- sp2_1
  tp2_2 <- sum(f2snp_2 %in% truesnp2)
  fp2_2 <- length(f2snp_2)-tp2_2
  fn2_2 <- length(setdiff(truesnp2,f2snp_2))
  tn2_2 <- length(setdiff(snp,f2snp_2))-fn2_2
  se2_2 <- tp2_2/(tp2_2+fn2_2) #敏感度
  sp2_2 <- tn2_2/(tn2_2+fp2_2) #特异度
  lse2_2[i] <- se2_2
  lsp2_2[i] <- sp2_2
  tp2_3 <- sum(f2snp_3 %in% truesnp2)
  fp2_3 <- length(f2snp_3)-tp2_3
  fn2_3 <- length(setdiff(truesnp2,f2snp_3))
  tn2_3 <- length(setdiff(snp,f2snp_3))-fn2_3
  se2_3 <- tp2_3/(tp2_3+fn2_3) #敏感度
  sp2_3 <- tn2_3/(tn2_3+fp2_3) #特异度
  lse2_3[i] <- se2_3
  lsp2_3[i] <- sp2_3
}
df <- data.frame(cutoff=cutoff,se1_1=lse1_1,se1_2=lse1_2,se1_3=lse1_3,se2_1=lse2_1,
      se2_2=lse2_2,se2_3=lse2_3,sp1_1=lsp1_1,sp1_2=lsp1_2,sp1_3=lsp1_3,sp2_1=lsp2_1,
      sp2_2=lsp2_2,sp2_3=lsp2_3)
write.csv(df,"re1.csv")

