#' Title Getting suitable genetic IVs through random graph forest, which is based on Bayesian network structure learning
#'
#' @param df a data frame which contains data of SNPs and specified exposure.
#' @param snp a vector of string belonging to column names of df, which is the name of SNPs included in BN structure learning.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param bn_method method for BN structure learning. Possible values are the function name of structure learning algorithm implemented in bnlearn. Default is "hc".
#' @param repeats an integer standing for the number of subsamples or bootstraps. Default is 1000.
#' @param selectNum the number of instrument to select. Default is 50.
#' @param nsam the size of individuals in each subsample of random graph forest. Default is 1000.
#' @param psam the size of variants in each subsample of random graph forest. Default is 100.
#' @param sample_replace is a boolean value to determine the sampling methods for individuals. TRUE with replacement and FALSE without replacement. Default is TRUE.
#'
#' @return a list containing:
#'
#'   selectsnp: a vector of string containing the colnames of df corresponding to
#'
#'   selected SNPs.
#'
#'   dfscore: a data frame containing the score calculated for each SNP.
#'   
#' @export
#'
#' @examples
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

