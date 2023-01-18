library(tidyverse)
set.seed(0)

simulaMDdata <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                         thetaA,thetaB,thetaA2,thetaM,r2AX,r2BY,
                         n,FM1,SE1,FM2,SE2,sigma2E){
  # library(rootSolve)
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
  
  qAX <- linkage(r2AX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(r2BY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,thetaA2,thetaM,FM2,SE2)

#   pais <- runif(5,0.1,0.4)
#   generaG <- function(pai){
#     return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
#   }
#   genoList <- lapply(pais,generaG)
#   geno <- matrix(unlist(genoList),nrow=n)
  
  df <- tibble(F1,F2,GA,GB,GM,GX,GY)
  return(df)
}

simulaMDdata2 <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                          theta,theta2,r2AX,r2BY,
                          n,FM1,SE1,FM2,SE2,sigma2E){
  # library(rootSolve)
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
  
  linkage <- function(r,pi1,pi2){
    #X is linked with A,p1=P(X|A),p2=P(X|a)
    # model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
    #                             F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
    # q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rqu,pai1,pai2))$root)
    q <- 0.5+0.5*r*sqrt(pi1*(1-pi1)/pi2*(1-pi2))
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
  
  qAX <- linkage(r2AX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(r2BY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
#   pais <- runif(5,0.1,0.4)
#   generaG <- function(pai){
#     return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
#   }
#   genoList <- lapply(pais,generaG)
#   geno <- matrix(unlist(genoList),nrow=n)
  df <- tibble(F1,F2,GA,GB,GM,GX,GY)#,geno)
  return(df)
}

simulaMDdata3 <- function(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                          theta,theta2,r2AX,r2BY,
                          n,FM1,SE1,FM2,SE2,sigma2E){
  # library(rootSolve)
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
  
  linkage <- function(r,pi1,pi2){
    #X is linked with A,p1=P(X|A),p2=P(X|a)
    # model <- function(x,parms)c(F1=(x[1]-x[2])^2*(parms[2]*(1-parms[2]))/(parms[3]*(1-parms[3]))-parms[1],
    #                             F2=x[1]+x[2]-1)  #assumption:P(X|A)=q,P(X|a)=1-q,paiA<paiX<1-paiA
    # q <- max(multiroot(f=model,start=c(0.5,0.5),parms=c(rqu,pai1,pai2))$root)
    q <- 0.5+0.5*r*sqrt(pi1*(1-pi1)/pi2*(1-pi2))
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
  
  qAX <- linkage(r2AX,paiA,paiX)
  GX <- sapply(seq(1,n),linkGeno,GA,qAX)
  qBY <- 1-linkage(r2BY,paiB,paiY)
  GY <- sapply(seq(1,n),linkGeno,GB,qBY)
  
  GM <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-paiM)^2,2*paiM*(1-paiM),paiM^2))

  F2 <- sapply(seq(1,n),phenotype,GA,GM,alpha2,theta2,FM2,SE2)
  
  # epsilon <- runif(n,0,sigma2E)
  # pai <- runif(1,0.1,0.4)
  # X <- sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2))
  # D <- 1.2*F1+3*F2+epsilon
  
#   pais <- runif(5,0.1,0.4)
#   generaG <- function(pai){
#     return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
#   }
#   genoList <- lapply(pais,generaG)
#   geno <- matrix(unlist(genoList),nrow=n)
  df <- tibble(F1,F2,GA,GB,GM,GX,GY)
  return(df)
}

# snpnameA <- sapply(seq(1,15),paste,"GA",sep="")
# snpnameB <- sapply(seq(1,15),paste,"GB",sep="")
# snpnameM <- sapply(seq(1,15),paste,"GM",sep="")
# snpnameX <- sapply(seq(1,15),paste,"GX",sep="")
# snpnameY <- sapply(seq(1,15),paste,"GY",sep="")
# phenoname1 <- sapply(seq(1,15),paste,"F1",sep="")
# phenoname2 <- sapply(seq(1,15),paste,"F2",sep="")
# xname1 <- sapply(seq(1,15),paste,"X1",sep="")
# xname2 <- sapply(seq(1,15),paste,"X2",sep="")
# xname3 <- sapply(seq(1,15),paste,"X3",sep="")
# xname4 <- sapply(seq(1,15),paste,"X4",sep="")
# xname5 <- sapply(seq(1,15),paste,"X5",sep="")

r <- 20
snpnameA <- paste0("GA",seq(1,3*r))
snpnameB <- paste0("GB",seq(1,3*r))
snpnameM <- paste0("GM",seq(1,3*r))
snpnameX <- paste0("GX",seq(1,3*r))
snpnameY <- paste0("GY",seq(1,3*r))
phenoname1 <- paste("F1",seq(1,3*r),sep="_")
phenoname2 <- paste("F2",seq(1,3*r),sep="_")


df <- data.frame(ID=seq(1,5000))
for(i in 1:r){
  paiA <-0.15#runif(1,0.12,0.3)
  paiB <-0.2#runif(1,0.12,0.3)
  paiM <-0.16#runif(1,0.12,0.3)
  paiX <-0.25#runif(1,0.12,0.3)
  paiY <-0.3#runif(1,0.12,0.3)
  alpha <- 0.2#runif(1,0.15,0.3)
  alpha2 <- 0.15#runif(1,0.15,0.3)
  thetaA <- 0.3#runif(1,0.2,0.5)
  thetaB <- 0.5#runif(1,0.2,0.5)
  thetaA2 <- 0.4#runif(1,0.2,0.5)
  thetaM <- 0.4#runif(1,0.2,0.5)
  r2AX <- 0.49#runif(1,0.3,0.7)  #r2_AX=0.49
  r2BY <- 0.36#runif(1,0.3,0.7)  #r2_BY=0.36
  n <- 5000
  FM1 <- 500#sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- 150#sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                      thetaA,thetaB,thetaA2,thetaM,r2AX,r2BY,
                      n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i])#,xname1[i],xname2[i],xname3[i],xname4[i],xname5[i])
  df <- bind_cols(df,dfs)
}
for(i in (r+1):(2*r)){
  paiA <-0.15#runif(1,0.12,0.3)
  paiB <-0.2#runif(1,0.12,0.3)
  paiM <-0.16#runif(1,0.12,0.3)
  paiX <-0.25#runif(1,0.12,0.3)
  paiY <-0.3#runif(1,0.12,0.3)
  alpha <- 0.2#runif(1,0.15,0.3)
  alpha2 <- 0.15#runif(1,0.15,0.3)
  theta <- 0.3#runif(1,0.25,0.6)
  theta2 <- 0.5#runif(1,0.25,0.6)
  r2AX <- 0.49#runif(1,0.3,0.7)
  r2BY <- 0.36#runif(1,0.3,0.7)
  n <- 5000
  FM1 <- 500#sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- 150#sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata2(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                       theta,theta2,r2AX,r2BY,
                       n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i])#,xname1[i],xname2[i],xname3[i],xname4[i],xname5[i])
  df <- bind_cols(df,dfs)
}
for(i in (2*r+1):(3*r)){
  paiA <-0.15#runif(1,0.12,0.3)
  paiB <-0.2#runif(1,0.12,0.3)
  paiM <-0.16#runif(1,0.12,0.3)
  paiX <-0.25#runif(1,0.12,0.3)
  paiY <-0.3#runif(1,0.12,0.3)
  alpha <- 0.2#runif(1,0.15,0.3)
  alpha2 <- 0.15#runif(1,0.15,0.3)
  theta <- 0.3#runif(1,0.25,0.6)
  theta2 <- 0.5#runif(1,0.25,0.6)
  r2AX <- 0.49#runif(1,0.3,0.7)
  r2BY <- 0.36#runif(1,0.3,0.7)
  n <- 5000
  FM1 <- 500#sample(2000:9000,size=1)
  SE1 <- sqrt(FM1)
  FM2 <- 150#sample(10:200,size=1)
  SE2 <- sqrt(FM2)/3
  sigma2E <- 20
  dfs <- simulaMDdata3(paiA,paiB,paiM,paiX,paiY,alpha,alpha2,
                       theta,theta2,r2AX,r2BY,
                       n,FM1,SE1,FM2,SE2,sigma2E)
  colnames(dfs)<-c(phenoname1[i],phenoname2[i],snpnameA[i],snpnameB[i],
                   snpnameM[i],snpnameX[i],snpnameY[i])#,xname1[i],xname2[i],xname3[i],xname4[i],xname5[i])
  df <- bind_cols(df,dfs)
}
df$F1 <- rowMeans(df[,phenoname1])
df$F2 <- rowMeans(df[,phenoname2])

n_ge <- 900
ename <- paste0("E",seq(1,n_ge))
pais <- runif(n_ge,0.1,0.4)
generaG <- function(pai){
  return(sample(c(0,1,2),n,replace=TRUE,prob = c((1-pai)^2,2*pai*(1-pai),pai^2)))
}
geno <- sapply(pais,generaG)
colnames(geno) <- ename
geno <- as_tibble(geno)

df %>% bind_cols(geno) %>% select(F1,F2,any_of(c(snpnameA,snpnameB,snpnameM,snpnameX,snpnameY,ename))) %>% write_csv("simData.csv")

# for(i in 2:15){
#   c <- phenoname1[i]
#   d <- phenoname2[i]
#   df$F1 <- df$F1+df[,c]
#   df$F2 <- df$F2+df[,d]
# }
# df$F1 <- df$F1/15
# df$F2 <- df$F2/15
# for(i in 2:15){
#   c <- phenoname1[i]
#   d <- phenoname2[i]
#   df$F1 <- df$F1*df[,c]
#   df$F2 <- df$F2*df[,d]
# }
# df$F1 <- abs(df$F1)^(1/15)
# df$F2 <- abs(df$F2)^(1/15)

# write_csv(df,"simData.csv")