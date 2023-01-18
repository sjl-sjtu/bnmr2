library(rstan)
library(MendelianRandomization)
library(AER)
library(ivmodel)
library(parallel)

source("inference/mr.R")

library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
df <- read_csv("inference/flexibleIV_n.csv")
df_iv <- read_csv("inference/snp_iv.csv")

snp <- grep("^g",colnames(df))
exposureName <- "X"
outcomeName <- "Y5"



FUN <- function(num){
    library(MendelianRandomization)
    library(AER)
    library(ivmodel)
    library(tidyverse)
    library(rstan)
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    IV <- df_iv %>% arrange(desc(score)) %>% slice_head(n=num) %>% pull(snp)
    timestart <- Sys.time()
    re <- mr(df,IV,"X","Y5",prior="horseshoe",n.iter=5000,n.chain=4)
    timeend <- Sys.time()
    t <- difftime(timeend,timestart,units="min")
    b <- re$mean
    s <- re$sd
    r <- re$Rhat
    L <- re$lower
    U <- re$upper
    res <- c(num,t,b,s,L,U,r)
}

numlist <- c(10,20,40,60,100,150)
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
clusterExport(cl,varlist=c("df","mr","df_iv"))
resL <- parLapply(cl, numlist, FUN)
stopCluster(cl)

dfre <- data.frame(t(matrix(unlist(resL),ncol=length(numlist))))
colnames(dfre) <- c("IVnum","time","est","se","L","U","R")
write_csv(dfre,"compareIV.csv")





