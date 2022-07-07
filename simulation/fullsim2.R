
bn_multi <- function(df,snp,exposureName,bn_method="hr",cutoff=0.7,repeats=100,nsam=500){
  library("bnlearn")
  library("plyr")
  library("dplyr")
  library("parallel")
  learnBN <- function(iter,df,nsam,bn_method){
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
    cores <- detectCores(logical = FALSE)
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {library("bnlearn")
                      library("plyr")
                      library("dplyr")
                      })
    clusterExport(cl,deparse(substitute(learnBN)),envir=environment())
    arcsL <- parLapply(cl,seq(repeats),learnBN,df,nsam,bn_method)
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

library(tidyverse)
df <- read_csv("simData.csv")

library(bnlearn)
library(plyr)
library(dplyr)
library(parallel)

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

timestart <- Sys.time()
dfre4 <- bn_multi(df1,snp,phenoName,bn_method="gs",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t4 <- difftime(timeend,timestart,units="min")

timestart <- Sys.time()
dfre5 <- bn_multi(df1,snp,phenoName,bn_method="rsmax2",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t5 <- difftime(timeend,timestart,units="min")

timestart <- Sys.time()
dfre6 <- bn_multi(df1,snp,phenoName,bn_method="pc.stable",cutoff=0.7,repeats=100,nsam=500)
timeend <- Sys.time()
t6 <- difftime(timeend,timestart,units="min")

dftime <- data.frame(algorithm=c("hc","iamb","mmhc","gs","rsmax2","pc.stable"),time=c(t1,t2,t3,t4,t5,t6))
write.csv(dftime,"time.csv")

#F1
dfs1_1 <- dfre1$score[[1]]
dfs1_2 <- dfre2$score[[1]]
dfs1_3 <- dfre3$score[[1]]
dfs1_4 <- dfre4$score[[1]]
dfs1_5 <- dfre5$score[[1]]
dfs1_6 <- dfre6$score[[1]]
#F2
dfs2_1 <- dfre1$score[[2]]
dfs2_2 <- dfre2$score[[2]]
dfs2_3 <- dfre3$score[[2]]
dfs2_4 <- dfre4$score[[2]]
dfs2_5 <- dfre5$score[[2]]
dfs2_6 <- dfre6$score[[2]]

write.csv(dfs1_1,"dfs1_1.csv",row.names=FALSE)
write.csv(dfs1_2,"dfs1_2.csv",row.names=FALSE)
write.csv(dfs1_3,"dfs1_3.csv",row.names=FALSE)
write.csv(dfs1_4,"dfs1_4.csv",row.names=FALSE)
write.csv(dfs1_5,"dfs1_5.csv",row.names=FALSE)
write.csv(dfs1_6,"dfs1_6.csv",row.names=FALSE)

write.csv(dfs2_1,"dfs2_1.csv",row.names=FALSE)
write.csv(dfs2_2,"dfs2_2.csv",row.names=FALSE)
write.csv(dfs2_3,"dfs2_3.csv",row.names=FALSE)
write.csv(dfs2_4,"dfs2_4.csv",row.names=FALSE)
write.csv(dfs2_5,"dfs2_5.csv",row.names=FALSE)
write.csv(dfs2_6,"dfs2_6.csv",row.names=FALSE)

truesnp1 <- c(snpnameA,snpnameB)
truesnp2 <- c(snpnameA,snpnameM)
lse1_1 <- c()
lsp1_1 <- c()
lse1_2 <- c()
lsp1_2 <- c()
lse1_3 <- c()
lsp1_3 <- c()
lse1_4 <- c()
lsp1_4 <- c()
lse1_5 <- c()
lsp1_5 <- c()
lse1_6 <- c()
lsp1_6 <- c()


lse2_1 <- c()
lsp2_1 <- c() 
lse2_2 <- c()
lsp2_2 <- c()
lse2_3 <- c()
lsp2_3 <- c()
lse2_4 <- c()
lsp2_4 <- c() 
lse2_5 <- c()
lsp2_5 <- c()
lse2_6 <- c()
lsp2_6 <- c()

cutoff <- seq(0,1,0.05)
for(i in 1:length(cutoff)){
  c <- cutoff[i]
  f1snp_1 <- dfs1_1[which(dfs1_1$score>=c),"snp"]
  f1snp_2 <- dfs1_2[which(dfs1_2$score>=c),"snp"]
  f1snp_3 <- dfs1_3[which(dfs1_3$score>=c),"snp"]
  f1snp_4 <- dfs1_4[which(dfs1_4$score>=c),"snp"]
  f1snp_5 <- dfs1_5[which(dfs1_5$score>=c),"snp"]
  f1snp_6 <- dfs1_6[which(dfs1_6$score>=c),"snp"]

  f2snp_1 <- dfs2_1[which(dfs2_1$score>=c),"snp"]
  f2snp_2 <- dfs2_2[which(dfs2_2$score>=c),"snp"]
  f2snp_3 <- dfs2_3[which(dfs2_3$score>=c),"snp"]
  f2snp_4 <- dfs2_4[which(dfs2_4$score>=c),"snp"]
  f2snp_5 <- dfs2_5[which(dfs2_5$score>=c),"snp"]
  f2snp_6 <- dfs2_6[which(dfs2_6$score>=c),"snp"]

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
  
  tp1_4 <- sum(f1snp_4 %in% truesnp1)
  fp1_4 <- length(f1snp_4)-tp1_4
  fn1_4 <- length(setdiff(truesnp1,f1snp_4))
  tn1_4 <- length(setdiff(snp,f1snp_4))-fn1_4
  se1_4 <- tp1_4/(tp1_4+fn1_4) #敏感度
  sp1_4 <- tn1_4/(tn1_4+fp1_4) #特异度
  lse1_4[i] <- se1_4
  lsp1_4[i] <- sp1_4

  tp1_5 <- sum(f1snp_5 %in% truesnp1)
  fp1_5 <- length(f1snp_5)-tp1_5
  fn1_5 <- length(setdiff(truesnp1,f1snp_5))
  tn1_5 <- length(setdiff(snp,f1snp_5))-fn1_5
  se1_5 <- tp1_5/(tp1_5+fn1_5) #敏感度
  sp1_5 <- tn1_4/(tn1_4+fp1_4) #特异度
  lse1_5[i] <- se1_5
  lsp1_5[i] <- sp1_5

  tp1_6 <- sum(f1snp_6 %in% truesnp1)
  fp1_6 <- length(f1snp_6)-tp1_6
  fn1_6 <- length(setdiff(truesnp1,f1snp_6))
  tn1_6 <- length(setdiff(snp,f1snp_6))-fn1_6
  se1_6 <- tp1_6/(tp1_6+fn1_6) #敏感度
  sp1_6 <- tn1_6/(tn1_6+fp1_6) #特异度
  lse1_6[i] <- se1_6
  lsp1_6[i] <- sp1_6


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
  
  tp2_4 <- sum(f2snp_4 %in% truesnp2)
  fp2_4 <- length(f2snp_4)-tp2_4
  fn2_4 <- length(setdiff(truesnp2,f2snp_4))
  tn2_4 <- length(setdiff(snp,f2snp_4))-fn2_4
  se2_4 <- tp2_4/(tp2_4+fn2_4) #敏感度
  sp2_4 <- tn2_4/(tn2_4+fp2_4) #特异度
  lse2_4[i] <- se2_4
  lsp2_4[i] <- sp2_4

  tp2_5 <- sum(f2snp_5 %in% truesnp2)
  fp2_5 <- length(f2snp_5)-tp2_5
  fn2_5 <- length(setdiff(truesnp2,f2snp_5))
  tn2_5 <- length(setdiff(snp,f2snp_5))-fn2_5
  se2_5 <- tp2_5/(tp2_5+fn2_5) #敏感度
  sp2_5 <- tn2_5/(tn2_5+fp2_5) #特异度
  lse2_5[i] <- se2_5
  lsp2_5[i] <- sp2_5

  tp2_6 <- sum(f2snp_6 %in% truesnp2)
  fp2_6 <- length(f2snp_6)-tp2_6
  fn2_6 <- length(setdiff(truesnp2,f2snp_6))
  tn2_6 <- length(setdiff(snp,f2snp_6))-fn2_6
  se2_6 <- tp2_6/(tp2_6+fn2_6) #敏感度
  sp2_6 <- tn2_6/(tn2_6+fp2_6) #特异度
  lse2_6[i] <- se2_6
  lsp2_6[i] <- sp2_6

}
dfre <- data.frame(cutoff=cutoff,
                   se1_1=lse1_1,se1_2=lse1_2,se1_3=lse1_3,
                   se1_4=lse1_4,se1_5=lse1_5,se1_6=lse1_6,
                   se2_1=lse2_1,se2_2=lse2_2,se2_3=lse2_3,
                   se2_4=lse2_4,se2_5=lse2_5,se2_6=lse2_6,
                   sp1_1=lsp1_1,sp1_2=lsp1_2,sp1_3=lsp1_3,
                   sp1_4=lsp1_4,sp1_5=lsp1_5,sp1_6=lsp1_6,
                   sp2_1=lsp2_1,sp2_2=lsp2_2,sp2_3=lsp2_3,
                   sp2_4=lsp2_4,sp2_5=lsp2_5,sp2_6=lsp2_6)
write.csv(dfre,"re1.csv",row.names=FALSE)


