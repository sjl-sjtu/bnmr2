set.seed(120)
bn <- function(df,snp,exposureName,bn_method="hr",cutoff=0.7,repeats=100,nsam=500){
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
    #exposureName is a str, snp is a vector of str.
    score <- rep(0,length(snp))
    for(i in 1:length(snp)){
      sn <- snp[i]
      count1 <- dfre[which(dfre$from==sn&dfre$to==exposureName),"count"]
      count2 <- dfre[which(dfre$from==exposureName&dfre$to==sn),"count"]
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
    return(dfscore)
  }
  
  df1 <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]
  
  dfre <- BNbootstrap(df1,repeats,nsam,bn_method)
  dfscore <- getscore(dfre,exposureName,snp,repeats)
  selectsnp <- dfscore[which(dfscore$score>=cutoff),"snp"]
  
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
  
  phenotype1 <- function(i,G1,alpha,theta1,FM,SE){
    a <- G1[i]
    if(a==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1){
      d <- rnorm(1,FM*alpha*(1+theta1),SE)
    }else if(a==2){
      d <- rnorm(1,FM*alpha*(1+theta1)^2,SE)
    }
  }
  
  #F2 <- sapply(seq(1,n),phenotype1,GM,alpha2,thetaM,FM2,SE2)
  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,thetaA2,thetaM,FM2,SE2)
  
  epsilon <- runif(n,0,sigma2E)
  pai <- runif(1,0.1,0.4)
  X <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
  D <- 1.2*F1+3*F2+20*GA+30*X+epsilon
  
  pais <- runif(100,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(D,F1,F2,GA,GB,GM,GX,GY,geno)
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
  
  phenotype1 <- function(i,G1,alpha,theta,FM,SE){
    a <- G1[i]
    if(a==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1){
      d <- rnorm(1,FM*alpha*(1+theta1),SE)
    }else if(a==2){
      d <- rnorm(1,FM*alpha*(1+theta1)^2,SE)
    }
  }
  
  #F2 <- sapply(seq(1,n),phenotype1,GM,alpha2,thetaM,FM2,SE2)
  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
  epsilon <- runif(n,0,sigma2E)
  pai <- runif(1,0.1,0.4)
  X <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
  D <- 3*F2+epsilon
  
  pais <- runif(100,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(D,F1,F2,GA,GB,GM,GX,GY,geno)
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
  
  phenotype1 <- function(i,G1,alpha,theta1,FM,SE){
    a <- G1[i]
    if(a==0){
      d <- rnorm(1,FM*alpha,SE)
    }else if(a==1){
      d <- rnorm(1,FM*alpha*(1+theta1),SE)
    }else if(a==2){
      d <- rnorm(1,FM*alpha*(1+theta1)^2,SE)
    }
  }
  
  #F2 <- sapply(seq(1,n),phenotype1,GM,alpha2,thetaM,FM2,SE2)
  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
  epsilon <- runif(n,0,sigma2E)
  pai <- runif(1,0.1,0.4)
  X <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
  D <- 3*F2+epsilon
  
  pais <- runif(100,0.1,0.4)
  generaG <- function(pai){
    return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
  }
  genoList <- lapply(pais,generaG)
  geno <- matrix(unlist(genoList),nrow=n)
  df <- data.frame(D,F1,F2,GA,GB,GM,GX,GY,geno)
  return(df)
}

paiA <- 0.12
paiB <- 0.08
paiM <- 0.25
paiX <- 0.3
paiY <- 0.21
alpha <- 0.2
alpha2 <- 0.25
thetaA <- 0.4
thetaB <- 0.5
thetaA2 <- 0.25
thetaM <- 0.35
rsAX <- 0.5
rsBY <- 0.32
n <- 5000
FM1 <- 700
SE1 <- sqrt(FM1)
FM2 <- 50
SE2 <- sqrt(FM2)/3
sigma2E <- 20

#df <- simulaMDdata(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,thetaA,thetaB,thetaA2,thetaM,rsAX,rsBY,n,FM1,SE1,FM2,SE2,sigma2E)

theta <- 0.4
theta2 <- 0.25
#df <- simulaMDdata2(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,theta,theta2,rsAX,rsBY,n,FM1,SE1,FM2,SE2,sigma2E)

df <- simulaMDdata3(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,theta,theta2,rsAX,rsBY,n,FM1,SE1,FM2,SE2,sigma2E)

df1 <- df[,-which(colnames(df)=="D")]

phenoName <- "F1"
snp <- colnames(df)[-c(1,2,3)]
snphc <- bn(df,snp,phenoName,"hc",0.9,100,500)
snppc <- bn(df,snp,phenoName,"iamb",0.9,100,500)
snpmmhc <- bn(df,snp,phenoName,"mmhc",0.9,100,500)
print(phenoName)
print("hc")
snphc
print("iamb")
snppc
print("mmhc")
snpmmhc

phenoName <- "F2"
snp <- colnames(df)[-c(1,2,3)]
snphc <- bn(df,snp,phenoName,"hc",0.9,100,500)
snppc <- bn(df,snp,phenoName,"iamb",0.9,100,500)
snpmmhc <- bn(df,snp,phenoName,"mmhc",0.9,100,500)
print(phenoName)
print("hc")
snphc
print("iamb")
snppc
print("mmhc")
snpmmhc