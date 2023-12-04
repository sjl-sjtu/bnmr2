#' Getting suitable genetic IVs through random graph forest, which is based on Bayesian network structure learning
#'
#' @param df a data frame which contains data of SNPs and specified exposure. The values of snps in the data frame should be either numeric or factors (not integers) for BN learning.
#' @param snp a vector of string belonging to column names of df, which is the name of SNPs included in BN structure learning.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param bn_method method for BN structure learning. Possible values are the function name of structure learning algorithm implemented in bnlearn. Default is "hc".
#' @param repeats an integer standing for the number of subsamples or bootstraps. Default is 1000.
#' @param selectNum the number of instrument to select. Default is NA.
#' @param alpha a number between 0 and 1 to specify the threshold for IV selection. We will use a threshold for variant selection as alpha*psam/length(snp). If selectNum is specified, the parameter will not be used. Default is 0.9.
#' @param nsam the size of individuals in each subsample of random graph forest. Default is 1000.
#' @param psam the size of variants in each subsample of random graph forest. Default is 100.
#' @param sample_replace is a boolean value to determine the sampling methods for individuals. TRUE with replacement and FALSE without replacement. Default is TRUE.
#'
#' @return a list containing:
#'   \item{selectsnp}{a vector of string containing the colnames of df corresponding to selected SNPs.}
#'   \item{dfscore}{a data frame containing the score calculated for each SNP.}
#'
#' @export
#'
#' @examples
#' n <- 2000
#' p <- 200
#' snps <- replicate(p,sample(1:3,n,replace = TRUE))
#' snps <- apply(snps,2,as.numeric)
#' snpname <- paste0("g",1:p)
#' df <- as.data.frame(snps)
#' colnames(df) <- snpname
#' truesnp <- paste0("g",sample(1:p,50))
#' df$x <- as.matrix(df[,truesnp])%*%rnorm(50,0.05,0.05)+rnorm(n,0,1)
#' df$y <- 0.5*df$x+rnorm(n,0,1)
#'
#' model <- bn(df,snpname,"x")
#'
bn <- function(df,snp,exposureName,bn_method="hc",repeats=1000,selectNum=NA,alpha=0.9,nsam=1000,psam=100,sample_replace=TRUE){
  library("bnlearn")
  library("plyr")
  library("dplyr")
  library("parallel")

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- 2L
  } else {
    # use all cores in devtools::test()
    cores <- parallel::detectCores(logical = FALSE)
  }

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
      model <- iamb.fdr(dfSam, undirected = TRUE)
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

  if(!is.na(selectNum)){
    selectsnp <- dfscore[1:selectNum,"snp"]
  }else{
    selectsnp <- dfscore%>%filter(score>=alpha*psam/length(snp))%>%pull(snp)
  }

  re <- list(IV=selectsnp,score=dfscore)
  return(re)
}

