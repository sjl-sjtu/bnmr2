library(rstan)
library(MendelianRandomization)
library(AER)
library(ivmodel)
library(parallel)

source("inference/mr.R")

library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
df <- read_csv("inference/fixedIV_2_n12.csv")
df_iv <- read_csv("inference/snp_iv.csv")

snp <- grep("^g",colnames(df))
exposureName <- "X"
outcomeName <- "Y7"
IV <- df_iv %>% arrange(desc(score)) %>% slice_head(n=20) %>% pull(snp)



FUN <- function(iter){
    library(MendelianRandomization)
    library(AER)
    library(ivmodel)
    library(tidyverse)
    library(rstan)
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    
    timestart <- Sys.time()
    re <- mr(df,IV,"X","Y3",prior="horseshoe",n.iter=iter,n.chain=4)
    timeend <- Sys.time()
    t <- difftime(timeend,timestart,units="min")
    b <- re$mean
    s <- re$sd
    r <- re$Rhat
    L <- re$lower
    U <- re$upper
    res <- c(iter,t,b,s,L,U,r)
}

iterlist <- c(1000,2000,5000,10000,20000)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","IV"))
resL <- parLapply(cl, iterlist, FUN)
stopCluster(cl)

dfre <- data.frame(t(matrix(unlist(resL),ncol=length(iterlist))))
colnames(dfre) <- c("iter","time","est","se","L","U","R")
write_csv(dfre,"compareIter2.csv")





