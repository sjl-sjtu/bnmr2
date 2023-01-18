

library(tidyverse)
setwd("./SunJianle/BNMR/simulation2211/")
df <- read_csv("simData.csv")
source("bn.R")
set.seed(0)

library(bnlearn)
library(plyr)
library(dplyr)
library(parallel)

r <- 20
snpnameA <- paste0("GA",seq(1,3*r))
snpnameB <- paste0("GB",seq(1,3*r))
snpnameM <- paste0("GM",seq(1,3*r))
snpnameX <- paste0("GX",seq(1,3*r))
snpnameY <- paste0("GY",seq(1,3*r))
n_ge <- 900
ename <- paste0("E",seq(1,n_ge))


snp <- as.character(setdiff(colnames(df),c("F1","F2")))
phenoName <- "F1"#as.character(c("F1","F2"))
truesnp1 <- c(snpnameA,snpnameB)


ns <- 1000
ps <- 100
r <- 1000
pt <- 50

timestart <- Sys.time()
dfre1 <- bn(df,snp,phenoName,bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
print(t1)

dfs <- dfre1$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)

timestart <- Sys.time()
dfre2 <- bn(df,snp,phenoName,bn_method="iamb",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")
print(t2)

timestart <- Sys.time()
dfre3 <- bn(df,snp,phenoName,bn_method="mmhc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t3 <- difftime(timeend,timestart,units="min")
print(t3)

timestart <- Sys.time()
dfre4 <- bn(df,snp,phenoName,bn_method="gs",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t4 <- difftime(timeend,timestart,units="min")
print(t4)

timestart <- Sys.time()
dfre5 <- bn(df,snp,phenoName,bn_method="rsmax2",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t5 <- difftime(timeend,timestart,units="min")
print(t5)

timestart <- Sys.time()
dfre6 <- bn(df,snp,phenoName,bn_method="pc.stable",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t6 <- difftime(timeend,timestart,units="min")
print(t6)

timestart <- Sys.time()
dfre7 <- bn(df,snp,phenoName,bn_method="tabu",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t7 <- difftime(timeend,timestart,units="min")
print(t7)

x <- list(dfre1$score,dfre2$score,dfre3$score,dfre4$score,dfre5$score,dfre6$score,dfre7$score)
x <- lapply(x,as_tibble)
dfre <- reduce(x,left_join,by="snp")
colnames(dfre)[-1] <- paste0("score",1:7)
dfre <- dfre %>% mutate(label=ifelse(snp%in%truesnp1,1,0))


dftime <- tibble(algorithm=c("hc","iamb","mmhc","gs","rsmax2","pc.stable","tabu"),time=c(t1,t2,t3,t4,t5,t6,t7))
dftime$auc <- map_dbl(1:7, function(i) auc(roc(dfre$label,dfre%>%pull(paste0("score",i)))))
dftime$pr_50 <- map_dbl(1:7, function(i) (dfre%>%dplyr::arrange(desc(paste0("score",i)))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt)
dftime$pr_25 <- map_dbl(1:7, function(i) (dfre%>%dplyr::arrange(desc(paste0("score",i)))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2))
dftime

write_csv(dftime,"epistas_sum.csv")
write_csv(dfre,"epistas_det.csv")


################
phenoName <- "F2"#as.character(c("F1","F2"))
truesnp1 <- c(snpnameA,snpnameM)


ns <- 1000
ps <- 100
r <- 1000
pt <- 50

timestart <- Sys.time()
dfre1 <- bn(df,snp,phenoName,bn_method="hc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t1 <- difftime(timeend,timestart,units="min")
print(t1)

dfs <- dfre1$score %>% mutate(label=ifelse(snp%in%truesnp1,1,0))
roc_object <- roc(dfs$label,dfs$score)
recal <- (dfs%>%arrange(desc(score))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt
recal1 <- (dfs%>%arrange(desc(score))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2)
AUC <- auc(roc_object)

timestart <- Sys.time()
dfre2 <- bn(df,snp,phenoName,bn_method="iamb",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t2 <- difftime(timeend,timestart,units="min")
print(t2)

timestart <- Sys.time()
dfre3 <- bn(df,snp,phenoName,bn_method="mmhc",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t3 <- difftime(timeend,timestart,units="min")
print(t3)

timestart <- Sys.time()
dfre4 <- bn(df,snp,phenoName,bn_method="gs",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t4 <- difftime(timeend,timestart,units="min")
print(t4)

timestart <- Sys.time()
dfre5 <- bn(df,snp,phenoName,bn_method="rsmax2",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t5 <- difftime(timeend,timestart,units="min")
print(t5)

timestart <- Sys.time()
dfre6 <- bn(df,snp,phenoName,bn_method="pc.stable",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t6 <- difftime(timeend,timestart,units="min")
print(t6)

timestart <- Sys.time()
dfre7 <- bn(df,snp,phenoName,bn_method="tabu",selectNum=pt,repeats=r,nsam=ns,psam=ps)
timeend <- Sys.time()
t7 <- difftime(timeend,timestart,units="min")
print(t7)

x <- list(dfre1$score,dfre2$score,dfre3$score,dfre4$score,dfre5$score,dfre6$score,dfre7$score)
x <- lapply(x,as_tibble)
dfre <- reduce(x,left_join,by="snp")
colnames(dfre)[-1] <- paste0("score",1:7)
dfre <- dfre %>% mutate(label=ifelse(snp%in%truesnp1,1,0))



dftime <- tibble(algorithm=c("hc","iamb","mmhc","gs","rsmax2","pc.stable","tabu"),time=c(t1,t2,t3,t4,t5,t6,t7))
dftime$auc <- map_dbl(1:7, function(i) auc(roc(dfre$label,dfre%>%pull(paste0("score",i)))))
dftime$pr_50 <- map_dbl(1:7, function(i) (dfre%>%dplyr::arrange(desc(paste0("score",i)))%>%slice_head(n=pt)%>%pull(label)%>%sum())/pt)
dftime$pr_25 <- map_dbl(1:7, function(i) (dfre%>%dplyr::arrange(desc(paste0("score",i)))%>%slice_head(n=round(pt/2))%>%pull(label)%>%sum())/round(pt/2))
dftime

write_csv(dftime,"epistas_sum_f2.csv")
write_csv(dfre,"epistas_det_f2.csv")


###############

associaT<-function(j,pheno){
    Z <- df%>%pull(j)
    X <- df%>%pull(pheno)
    f <- lm(X~Z)
    return(summary(f)$coefficients[2,4])
}
d1 <- data.frame(snp=snp,p=sapply(snp,associaT,"X1"))
write.csv(d1,"assoF1.csv",row.names=FALSE)
d2 <- data.frame(snp=snp,p=sapply(snp,associaT,"X2"))
write.csv(d2,"assoF1.csv",row.names=FALSE)




# truesnp1 <- c(snpnameA,snpnameB)
# truesnp2 <- c(snpnameA,snpnameM)
# lse1_1 <- c()
# lsp1_1 <- c()
# lse1_2 <- c()
# lsp1_2 <- c()
# lse1_3 <- c()
# lsp1_3 <- c()
# lse1_4 <- c()
# lsp1_4 <- c()
# lse1_5 <- c()
# lsp1_5 <- c()
# lse1_6 <- c()
# lsp1_6 <- c()


# lse2_1 <- c()
# lsp2_1 <- c() 
# lse2_2 <- c()
# lsp2_2 <- c()
# lse2_3 <- c()
# lsp2_3 <- c()
# lse2_4 <- c()
# lsp2_4 <- c() 
# lse2_5 <- c()
# lsp2_5 <- c()
# lse2_6 <- c()
# lsp2_6 <- c()

# cutoff <- seq(0,1,0.05)
# for(i in 1:length(cutoff)){
#   c <- cutoff[i]
#   f1snp_1 <- dfs1_1[which(dfs1_1$score>=c),"snp"]
#   f1snp_2 <- dfs1_2[which(dfs1_2$score>=c),"snp"]
#   f1snp_3 <- dfs1_3[which(dfs1_3$score>=c),"snp"]
#   f1snp_4 <- dfs1_4[which(dfs1_4$score>=c),"snp"]
#   f1snp_5 <- dfs1_5[which(dfs1_5$score>=c),"snp"]
#   f1snp_6 <- dfs1_6[which(dfs1_6$score>=c),"snp"]

#   f2snp_1 <- dfs2_1[which(dfs2_1$score>=c),"snp"]
#   f2snp_2 <- dfs2_2[which(dfs2_2$score>=c),"snp"]
#   f2snp_3 <- dfs2_3[which(dfs2_3$score>=c),"snp"]
#   f2snp_4 <- dfs2_4[which(dfs2_4$score>=c),"snp"]
#   f2snp_5 <- dfs2_5[which(dfs2_5$score>=c),"snp"]
#   f2snp_6 <- dfs2_6[which(dfs2_6$score>=c),"snp"]

#   tp1_1 <- sum(f1snp_1 %in% truesnp1)
#   fp1_1 <- length(f1snp_1)-tp1_1
#   fn1_1 <- length(setdiff(truesnp1,f1snp_1))
#   tn1_1 <- length(setdiff(snp,f1snp_1))-fn1_1
#   se1_1 <- tp1_1/(tp1_1+fn1_1) #敏感度
#   sp1_1 <- tn1_1/(tn1_1+fp1_1) #特异度
#   lse1_1[i] <- se1_1
#   lsp1_1[i] <- sp1_1

#   tp1_2 <- sum(f1snp_2 %in% truesnp1)
#   fp1_2 <- length(f1snp_2)-tp1_2
#   fn1_2 <- length(setdiff(truesnp1,f1snp_2))
#   tn1_2 <- length(setdiff(snp,f1snp_2))-fn1_2
#   se1_2 <- tp1_2/(tp1_2+fn1_2) #敏感度
#   sp1_2 <- tn1_2/(tn1_2+fp1_2) #特异度
#   lse1_2[i] <- se1_2
#   lsp1_2[i] <- sp1_2

#   tp1_3 <- sum(f1snp_3 %in% truesnp1)
#   fp1_3 <- length(f1snp_3)-tp1_3
#   fn1_3 <- length(setdiff(truesnp1,f1snp_3))
#   tn1_3 <- length(setdiff(snp,f1snp_3))-fn1_3
#   se1_3 <- tp1_3/(tp1_3+fn1_3) #敏感度
#   sp1_3 <- tn1_3/(tn1_3+fp1_3) #特异度
#   lse1_3[i] <- se1_3
#   lsp1_3[i] <- sp1_3
  
#   tp1_4 <- sum(f1snp_4 %in% truesnp1)
#   fp1_4 <- length(f1snp_4)-tp1_4
#   fn1_4 <- length(setdiff(truesnp1,f1snp_4))
#   tn1_4 <- length(setdiff(snp,f1snp_4))-fn1_4
#   se1_4 <- tp1_4/(tp1_4+fn1_4) #敏感度
#   sp1_4 <- tn1_4/(tn1_4+fp1_4) #特异度
#   lse1_4[i] <- se1_4
#   lsp1_4[i] <- sp1_4

#   tp1_5 <- sum(f1snp_5 %in% truesnp1)
#   fp1_5 <- length(f1snp_5)-tp1_5
#   fn1_5 <- length(setdiff(truesnp1,f1snp_5))
#   tn1_5 <- length(setdiff(snp,f1snp_5))-fn1_5
#   se1_5 <- tp1_5/(tp1_5+fn1_5) #敏感度
#   sp1_5 <- tn1_4/(tn1_4+fp1_4) #特异度
#   lse1_5[i] <- se1_5
#   lsp1_5[i] <- sp1_5

#   tp1_6 <- sum(f1snp_6 %in% truesnp1)
#   fp1_6 <- length(f1snp_6)-tp1_6
#   fn1_6 <- length(setdiff(truesnp1,f1snp_6))
#   tn1_6 <- length(setdiff(snp,f1snp_6))-fn1_6
#   se1_6 <- tp1_6/(tp1_6+fn1_6) #敏感度
#   sp1_6 <- tn1_6/(tn1_6+fp1_6) #特异度
#   lse1_6[i] <- se1_6
#   lsp1_6[i] <- sp1_6


#   tp2_1 <- sum(f2snp_1 %in% truesnp2)
#   fp2_1 <- length(f2snp_1)-tp2_1
#   fn2_1 <- length(setdiff(truesnp2,f2snp_1))
#   tn2_1 <- length(setdiff(snp,f2snp_1))-fn2_1
#   se2_1 <- tp2_1/(tp2_1+fn2_1) #敏感度
#   sp2_1 <- tn2_1/(tn2_1+fp2_1) #特异度
#   lse2_1[i] <- se2_1
#   lsp2_1[i] <- sp2_1

#   tp2_2 <- sum(f2snp_2 %in% truesnp2)
#   fp2_2 <- length(f2snp_2)-tp2_2
#   fn2_2 <- length(setdiff(truesnp2,f2snp_2))
#   tn2_2 <- length(setdiff(snp,f2snp_2))-fn2_2
#   se2_2 <- tp2_2/(tp2_2+fn2_2) #敏感度
#   sp2_2 <- tn2_2/(tn2_2+fp2_2) #特异度
#   lse2_2[i] <- se2_2
#   lsp2_2[i] <- sp2_2

#   tp2_3 <- sum(f2snp_3 %in% truesnp2)
#   fp2_3 <- length(f2snp_3)-tp2_3
#   fn2_3 <- length(setdiff(truesnp2,f2snp_3))
#   tn2_3 <- length(setdiff(snp,f2snp_3))-fn2_3
#   se2_3 <- tp2_3/(tp2_3+fn2_3) #敏感度
#   sp2_3 <- tn2_3/(tn2_3+fp2_3) #特异度
#   lse2_3[i] <- se2_3
#   lsp2_3[i] <- sp2_3
  
#   tp2_4 <- sum(f2snp_4 %in% truesnp2)
#   fp2_4 <- length(f2snp_4)-tp2_4
#   fn2_4 <- length(setdiff(truesnp2,f2snp_4))
#   tn2_4 <- length(setdiff(snp,f2snp_4))-fn2_4
#   se2_4 <- tp2_4/(tp2_4+fn2_4) #敏感度
#   sp2_4 <- tn2_4/(tn2_4+fp2_4) #特异度
#   lse2_4[i] <- se2_4
#   lsp2_4[i] <- sp2_4

#   tp2_5 <- sum(f2snp_5 %in% truesnp2)
#   fp2_5 <- length(f2snp_5)-tp2_5
#   fn2_5 <- length(setdiff(truesnp2,f2snp_5))
#   tn2_5 <- length(setdiff(snp,f2snp_5))-fn2_5
#   se2_5 <- tp2_5/(tp2_5+fn2_5) #敏感度
#   sp2_5 <- tn2_5/(tn2_5+fp2_5) #特异度
#   lse2_5[i] <- se2_5
#   lsp2_5[i] <- sp2_5

#   tp2_6 <- sum(f2snp_6 %in% truesnp2)
#   fp2_6 <- length(f2snp_6)-tp2_6
#   fn2_6 <- length(setdiff(truesnp2,f2snp_6))
#   tn2_6 <- length(setdiff(snp,f2snp_6))-fn2_6
#   se2_6 <- tp2_6/(tp2_6+fn2_6) #敏感度
#   sp2_6 <- tn2_6/(tn2_6+fp2_6) #特异度
#   lse2_6[i] <- se2_6
#   lsp2_6[i] <- sp2_6

# }
# dfre <- data.frame(cutoff=cutoff,
#                    se1_1=lse1_1,se1_2=lse1_2,se1_3=lse1_3,
#                    se1_4=lse1_4,se1_5=lse1_5,se1_6=lse1_6,
#                    se2_1=lse2_1,se2_2=lse2_2,se2_3=lse2_3,
#                    se2_4=lse2_4,se2_5=lse2_5,se2_6=lse2_6,
#                    sp1_1=lsp1_1,sp1_2=lsp1_2,sp1_3=lsp1_3,
#                    sp1_4=lsp1_4,sp1_5=lsp1_5,sp1_6=lsp1_6,
#                    sp2_1=lsp2_1,sp2_2=lsp2_2,sp2_3=lsp2_3,
#                    sp2_4=lsp2_4,sp2_5=lsp2_5,sp2_6=lsp2_6)
# write.csv(dfre,"re1.csv",row.names=FALSE)

