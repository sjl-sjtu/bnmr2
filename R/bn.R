#' @title Getting suitable genetic IVs through Bayesian network learning
#'
#' @description is used to get the suitable SNPs as instrumental variables (IVs) of specified exposure by Bayesian network (BN) structure learning.
#'
#' @param df a data frame which contains data of SNPs and specified exposure.
#' @param snp a vector of string belonging to colnames of df, which is the name of SNPs included in BN structure learning.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param bn_method method for BN structure learning. Possible values are "hc", "pc" or "mix". Default is "hc".
#' @param cutoff a numeric between 0 to 1. Those SNPs with score larger than "cutoff" will be chosen as IVs. Default is 0.7.
#' @param repeats an integer standing for the times of bootstraping. Default is 100.
#' @param nsam an integer standing for the sample size for bootstraping sampling. Default is 500.
#'
#' @return a list containing:
#'   selectsnp: a vector of string containing the colnames of df corresponding to
#'   selected SNPs.
#'   dfscore: a data frame containing the score calculated for each SNP.
#' @export
#'
#' @examples
#'
#'
bn <- function(df,snp,exposureName,bn_method="hr",cutoff=0.7,repeats=100,nsam=500){
  library("bnlearn")

  learnBNpc <- function(df,nsam){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    model <- pc.stable(dfSam)
    dfarc <- data.frame(model$arcs)
    return(dfarc)
  }

  learnBNhc <- function(df,nsam){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    model <- hc(dfSam)
    dfarc <- data.frame(model$arcs)
    return(dfarc)
  }

  rmBidire <- function(df){
    delL <- c()
    for(i in 1:nrow(df)){
      for(j in i+1:nrow(df)){
        if(all(df[j,]%in%df[i,])){
          delL <- append(delL,j)
        }
      }
    }
    df <- df[-delL,]
    return(df)
  }

  BNbootstrapPC <- function(df,repeats,nsam){
    arcsL <- replicate(repeats,learnBNpc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL <- rmBidire(arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }

  BNbootstrapHC <- function(df,repeats,nsam){
    arcsL <- replicate(repeats,learnBNhc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }

  BNbootstrapBi <- function(df,repeats,nsam){
    arcsL <- replicate(repeats/2,learnBNpc(df,nsam),simplify = FALSE)
    arcsL2 <- replicate(repeats/2,learnBNhc(df,nsam),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL <- rmBidire(arcsL)
    arcsL2 <- do.call(rbind.fill,arcsL2)
    arcsL <- rbind(arcsL,arcsL2)
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
    return(dfscore)
  }

  df1 <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]


  if(bn_method=="pc"){
    dfre <- BNbootstrapPC(df1,repeats,nsam)
  }else if(bn_method=="hc"){
    dfre <- BNbootstrapHC(df1,repeats,nsam)
  }else if(bn_method=="mix"){
    dfre <- BNbootstrapBi(df1,repeats,nsam)
  }else{
    return(message("no this bn learning method"))
  }

  dfscore <- getscore(dfre,exposureName,snp,repeats)
  selectsnp <- dfscore[which(dfscore$score>=cutoff),"snp"]

  re <- list(IV=selectsnp,score=dfscore)
  return(re)
}
