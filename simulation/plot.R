library(tidyverse)
library(stringr)
library(latex2exp)
library(ggthemes)
library(purrr)
library(ggsci)
library(RColorBrewer)

df <- read_csv("compareTime_1.csv")
df1 <- df%>% filter(ns==2000) %>% 
  arrange(ps) %>% dplyr::select(-ns) %>% 
  #mutate(recall=1-recall,recall1=1-recall1) %>%
  pivot_longer(cols = -c(ps,t), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(ps),as.factor) %>%
  mutate_at(vars(indices),factor,levels=c("recall","recall1","auc"))
ggplot(df1)+
  #geom_col(aes(x=ps,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  geom_col(aes(x=ps,y=t/10),width=0.6,color="black",fill=brewer.pal(9,"BuGn")[2],position = position_dodge(preserve = "total"))+
  geom_point(aes(x=ps,y=performance,colour=indices,group=indices))+
  geom_line(aes(x=ps,y=performance,colour=indices,group=indices))+
  xlab(TeX("$p_s$"))+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_aaas(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*10,name = "running time (min)"))+
  scale_fill_nejm()
  #geom_point(aes(x=ps,y=t/10,group=1))+
  
ggsave("g1.png", width=5.5, height=5, dpi = 300)

df2 <- df%>% filter(ps==100) %>% 
  arrange(ps) %>% dplyr::select(-ps) %>% 
  #mutate(recall=1-recall,recall1=1-recall1) %>%
  pivot_longer(cols = -c(ns,t), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(ns),as.factor) %>%
  mutate_at(vars(indices),factor,levels=c("recall","recall1","auc"))
ggplot(df2)+
  #geom_col(aes(x=ns,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  geom_col(aes(x=ns,y=t/10),width=0.6,fill=brewer.pal(9,"BuGn")[2],color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=ns,y=performance,colour=indices,group=indices))+
  geom_line(aes(x=ns,y=performance,colour=indices,group=indices))+
  xlab(TeX("$n_s$"))+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_aaas(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*10,name = "running time (min)"))
  # geom_point(aes(x=ns,y=t/10,group=1))+
  # geom_line(aes(x=ns,y=t/10,group=1),linewidth=1,color="black")
ggsave("g2.png", width=5.5, height=5, dpi = 300)

df3 <- read_csv("compareRepeats_1.csv") %>% dplyr::slice(-1) %>% arrange(r) %>% 
  #mutate(recall=1-recall,recall1=1-recall1) %>%
  pivot_longer(cols = -c(r,t), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(r),as.factor) %>%
  mutate_at(vars(indices),factor,levels=c("recall","recall1","auc"))
ggplot(df3)+
  geom_col(aes(x=r,y=t/10),width=0.6,fill=brewer.pal(9,"BuGn")[2],color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=r,y=performance,colour=indices,group=indices))+
  geom_line(aes(x=r,y=performance,colour=indices,group=indices))+
  # geom_col(aes(x=r,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  xlab(TeX("$r$"))+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_aaas(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*10,name = "running time (min)"))
  # geom_point(aes(x=r,y=t/10,group=1))+
  # geom_line(aes(x=r,y=t/10,group=1),linewidth=1,color="black")
ggsave("g3.png", width=5.5, height=5,dpi = 300)

df4 <- read_csv("compareMethods.csv") %>% dplyr::slice(-1) %>%
  #mutate(recall=1-recall,recall1=1-recall1) %>%
  pivot_longer(cols = -c(method,t), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(indices),factor,levels=c("recall","recall1","auc")) %>%
  mutate_at(vars(method),factor,levels=c("hc","tabu","pc.stable","gs","fast.iamb","mmhc","rsmax2"))
ggplot(df4)+
  geom_col(aes(x=method,y=t/12),width=0.6,fill=brewer.pal(9,"BuGn")[2],color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=method,y=performance,colour=indices,group=indices))+
  geom_line(aes(x=method,y=performance,colour=indices,group=indices))+
  # geom_col(aes(x=method,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  xlab("learning method")+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_aaas(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*12,name = "running time (min)"))
  # geom_point(aes(x=method,y=t/12,group=1))+
  # geom_line(aes(x=method,y=t/12,group=1),linewidth=1,color="black")
ggsave("g4.png", width=5.5, height=5, dpi = 300)


################
df <- read_csv("idprune_x1.csv") %>% mutate_at(vars(pthres,cut),as.factor)
df1 <- read_csv("x1_score_1e-4.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x1_score_1e-6.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x1_score_1e-8.csv") %>% dplyr::arrange(desc(score))
truesnp <- read_delim("truesnp_x1.txt") %>% pull(truesnp)
a <- c(100,50,25)
dfgwas <- read_csv("dfgwas_x1.csv") %>% arrange(p)
library(purrr)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s1 <- map_dbl(a,function(x) sum((df1%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s2 <- map_dbl(a,function(x) sum((df2%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s3 <- map_dbl(a,function(x) sum((df3%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)

ggplot(df)+
  geom_col(aes(x=cut,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  #geom_line(aes(x=cut,y=rec,colour=pthres,group=pthres))+
  xlab(TeX("$\\rho$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_aaas()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_1$"))#+
  # scale_y_continuous(sec.axis = sec_axis(~.*1000,name = "selected number"))+
  # geom_point(aes(x=cut,y=len/1000,group=1))+
  # geom_line(aes(x=cut,y=len/1000,group=1),linewidth=1,color="black")
ggsave("g5.png", width=5.5, height=5, dpi = 300)

dfbn <- tibble(num=a,th1=bn_s1,th2=bn_s2,th3=bn_s3) %>% 
  pivot_longer(cols = -c(num), names_to = "pthres", values_to = "accuracy") %>%
  mutate_at(vars(num),as.factor) %>%
  mutate_at(vars(pthres),dplyr::recode,th1=1e-4,th2=1e-6,th3=1e-8) %>%
  mutate_at(vars(pthres),as.factor)
ggplot(dfbn)+
  geom_col(aes(x=num,y=accuracy,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_stata()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1.2))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_1$"))
ggsave("g6.png", width=5.5, height=5, dpi = 300)

dfg <- tibble(num=a,accuracy=gwas_s) %>%
  mutate_at(vars(num),as.factor)
ggplot(dfg)+
  geom_col(aes(x=num,y=accuracy),fill=8,position = position_dodge(preserve = 'single'),color="black")+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  scale_fill_stata()+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(TeX("$X_1$"))+
  coord_cartesian(ylim = c(0,1))
ggsave("g7.png", width=5.5, height=5, dpi = 300)


df <- read_csv("idprune_x2.csv") %>% mutate_at(vars(pthres,cut),as.factor)
df1 <- read_csv("x2_score_1e-4.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x2_score_1e-6.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x2_score_1e-8.csv") %>% dplyr::arrange(desc(score))
truesnp <- read_delim("truesnp_x2.txt") %>% pull(truesnp)
a <- c(100,50,25)
dfgwas <- read_csv("dfgwas_x2.csv") %>% arrange(p)
library(purrr)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s1 <- map_dbl(a,function(x) sum((df1%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s2 <- map_dbl(a,function(x) sum((df2%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s3 <- map_dbl(a,function(x) sum((df3%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)

ggplot(df)+
  geom_col(aes(x=cut,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  xlab(TeX("$\\rho$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_stata()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_2$"))#+
# scale_y_continuous(sec.axis = sec_axis(~.*1000,name = "selected number"))+
# geom_point(aes(x=cut,y=len/1000,group=1))+
# geom_line(aes(x=cut,y=len/1000,group=1),linewidth=1,color="black")
ggsave("g8.png", width=5.5, height=5, dpi = 300)

dfbn <- tibble(num=a,th1=bn_s1,th2=bn_s2,th3=bn_s3) %>% 
  pivot_longer(cols = -c(num), names_to = "pthres", values_to = "accuracy") %>%
  mutate_at(vars(num),as.factor) %>%
  mutate_at(vars(pthres),dplyr::recode,th1=1e-4,th2=1e-6,th3=1e-8) %>%
  mutate_at(vars(pthres),as.factor)
ggplot(dfbn)+
  geom_col(aes(x=num,y=accuracy,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_stata()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1.2))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_2$"))
ggsave("g9.png", width=5.5, height=5, dpi = 300)

dfg <- tibble(num=a,accuracy=gwas_s) %>%
  mutate_at(vars(num),as.factor)
ggplot(dfg)+
  geom_col(aes(x=num,y=accuracy),fill=8,position = position_dodge(preserve = 'single'),color="black")+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  scale_fill_stata()+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_2$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(ylim = c(0,1))
ggsave("g10.png", width=5.5, height=5, dpi = 300)

df <- read_csv("penalized_x1.csv") %>% mutate_at(vars(pthres,alpha),as.factor)
ggplot(df)+
  geom_col(aes(x=alpha,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  xlab(TeX("$\\alpha$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_stata()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_1$"))
ggsave("g11.png", width=5.5, height=5, dpi = 300)

df <- read_csv("penalized_x2.csv") %>% mutate_at(vars(pthres,alpha),as.factor)
ggplot(df)+
  geom_col(aes(x=alpha,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),color="black")+
  xlab(TeX("$\\alpha$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995))+
  scale_fill_stata()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_2$"))
ggsave("g12.png", width=5.5, height=5, dpi = 300)


##########################
df1 <- read_csv("penalized_x1.csv") %>% mutate_at(vars(pthres,alpha),as.factor)
df2 <- read_csv("penalized_x2.csv") %>% mutate_at(vars(pthres,alpha),as.factor)
df <- bind_rows(df1,df2) %>% mutate(exposure=rep(c("X1","X2"),each=nrow(df1)))
ggplot(df)+facet_grid(cols=vars(exposure))+
  geom_col(aes(x=alpha,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),
           color="black",alpha = 0.3)+
  xlab(TeX("$\\alpha$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  scale_fill_aaas()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle("Penalized Regression")
ggsave("g13.png", width=7, height=5.6, dpi = 300)

df1 <- read_csv("idprune_x1.csv") %>% mutate_at(vars(pthres,cut),as.factor)
df2 <- read_csv("idprune_x2.csv") %>% mutate_at(vars(pthres,cut),as.factor)
df <- bind_rows(df1,df2) %>% mutate(exposure=rep(c("X1","X2"),each=nrow(df1)))
ggplot(df)+facet_grid(cols=vars(exposure))+
  geom_col(aes(x=cut,y=rec,fill=pthres),position = position_dodge(preserve = 'single'),
           color="black",alpha=0.3)+
  xlab(TeX("$r^2$"))+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  scale_fill_aaas()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle("LD Pruning")
ggsave("g14.png", width=7, height=5.6, dpi = 300)

df1 <- read_csv("x1_score_1_5000.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x1_score_1e-2_5000.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x1_score_1e-4_5000.csv") %>% dplyr::arrange(desc(score))
df4 <- read_csv("x1_score_1e-6_5000.csv") %>% dplyr::arrange(desc(score))
df5 <- read_csv("x1_score_1e-8_5000.csv") %>% dplyr::arrange(desc(score))
truesnp <- read_delim("truesnp_x1.txt") %>% pull(truesnp)
a <- c(60,40,20)
dfgwas <- read_csv("dfgwas_x1.csv") %>% arrange(p)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s1 <- map_dbl(a,function(x) sum((df1%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s2 <- map_dbl(a,function(x) sum((df2%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s3 <- map_dbl(a,function(x) sum((df3%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s4 <- map_dbl(a,function(x) sum((df4%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s5 <- map_dbl(a,function(x) sum((df5%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)

# alpha <- c(0.95,0.99,0.999)
# ps <- 150; pt <- 10000
# bn_s1 <- map_dbl(alpha,function(x) sum((df1%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df1%>%filter(score>=x*ps/pt)))
# bn_s2 <- map_dbl(alpha,function(x) sum((df2%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df2%>%filter(score>=x*ps/pt)))
# bn_s3 <- map_dbl(alpha,function(x) sum((df3%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df3%>%filter(score>=x*ps/pt)))
# bn_s4 <- map_dbl(alpha,function(x) sum((df4%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df4%>%filter(score>=x*ps/pt)))
# bn_s5 <- map_dbl(alpha,function(x) sum((df5%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df5%>%filter(score>=x*ps/pt)))
dfbn1 <- tibble(num=a,th1=bn_s1,th2=bn_s2,th3=bn_s3,th4=bn_s4,th5=bn_s5) %>% 
  pivot_longer(cols = -c(num), names_to = "pthres", values_to = "accuracy") %>%
  mutate_at(vars(num),as.factor) %>%
  mutate_at(vars(pthres),dplyr::recode,th1=1,th2=1e-2,th3=1e-4,th4=1e-6,th5=1e-8) %>%
  mutate_at(vars(pthres),as.factor)


df1 <- read_csv("x2_score_1_5000.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x2_score_1e-2_5000.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x2_score_1e-4_5000.csv") %>% dplyr::arrange(desc(score))
df4 <- read_csv("x2_score_1e-6_5000.csv") %>% dplyr::arrange(desc(score))
df5 <- read_csv("x2_score_1e-8_5000.csv") %>% dplyr::arrange(desc(score))
truesnp <- read_delim("truesnp_x2.txt") %>% pull(truesnp)
a <- c(60,40,20)
dfgwas <- read_csv("dfgwas_x2.csv") %>% arrange(p)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s1 <- map_dbl(a,function(x) sum((df1%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s2 <- map_dbl(a,function(x) sum((df2%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s3 <- map_dbl(a,function(x) sum((df3%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s4 <- map_dbl(a,function(x) sum((df4%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
bn_s5 <- map_dbl(a,function(x) sum((df5%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)

# alpha <- c(0.95,0.99,0.999)
# ps <- 150; pt <- 10000
# bn_s1 <- map_dbl(alpha,function(x) sum((df1%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df1%>%filter(score>=x*ps/pt)))
# bn_s2 <- map_dbl(alpha,function(x) sum((df2%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df2%>%filter(score>=x*ps/pt)))
# bn_s3 <- map_dbl(alpha,function(x) sum((df3%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df3%>%filter(score>=x*ps/pt)))
# bn_s4 <- map_dbl(alpha,function(x) sum((df4%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df4%>%filter(score>=x*ps/pt)))
# bn_s5 <- map_dbl(alpha,function(x) sum((df5%>%filter(score>=x*ps/pt)%>%pull(snp))%in%truesnp)/nrow(df5%>%filter(score>=x*ps/pt)))
dfbn2 <- tibble(num=a,th1=bn_s1,th2=bn_s2,th3=bn_s3,th4=bn_s4,th5=bn_s5) %>% 
  pivot_longer(cols = -c(num), names_to = "pthres", values_to = "accuracy") %>%
  mutate_at(vars(num),as.factor) %>%
  mutate_at(vars(pthres),dplyr::recode,th1=1,th2=1e-2,th3=1e-4,th4=1e-6,th5=1e-8) %>%
  mutate_at(vars(pthres),as.factor)

df <- bind_rows(dfbn1,dfbn2) %>% mutate(exposure=rep(c("X1","X2"),each=nrow(dfbn1)))
dfbn <- df%>%filter(pthres%in%c(1e-2,1e-4,1e-6,1e-8))
ggplot(dfbn)+facet_grid(cols=vars(exposure))+
  geom_col(aes(x=num,y=accuracy,fill=pthres),position = position_dodge(preserve = 'single'),
           color="black",alpha=0.3)+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  scale_fill_aaas()+
  labs(fill="GWAS threshold")+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle("BN")
ggsave("g15.png", width=7, height=5.6, dpi = 300)

dfgwas <- read_csv("dfgwas_x1.csv") %>% arrange(p)
truesnp <- read_delim("truesnp_x1.txt") %>% pull(truesnp)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
dfg1 <- tibble(num=a,accuracy=gwas_s) %>%
  mutate_at(vars(num),as.factor)
dfgwas <- read_csv("dfgwas_x2.csv") %>% arrange(p)
truesnp <- read_delim("truesnp_x2.txt") %>% pull(truesnp)
gwas_s <- map_dbl(a,function(x) sum((dfgwas%>%dplyr::slice(1:x)%>%pull(snp))%in%truesnp)/x)
dfg2 <- tibble(num=a,accuracy=gwas_s) %>%
  mutate_at(vars(num),as.factor)
df <- bind_rows(dfg1,dfg2) %>% mutate(exposure=rep(c("X1","X2"),each=nrow(dfg1))) %>% mutate(color="x")
ggplot(df)+facet_grid(cols=vars(exposure))+
  geom_col(aes(x=num,y=accuracy),position = position_dodge(preserve = 'single'),
           color="black",fill=brewer.pal(9,"BuPu")[5],alpha=0.4)+
  xlab("selected number")+
  ylab("TDR")+
  theme_few()+
  scale_fill_npg()+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  ggtitle(TeX("$X_2$"))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("GWAS filtering")+
  guides(fill="none")
ggsave("g16.png", width=7, height=5.6, dpi = 300)



###############
dfgwas <- read_csv("dfgwas_x1.csv") %>% arrange(p)
df1 <- read_csv("x1_score_1_5000_1.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x1_score_1e-2_5000_1.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x1_score_1e-4_5000_1.csv") %>% dplyr::arrange(desc(score))
#df0 <- read_csv("x1_score_1e-2.csv") %>% dplyr::arrange(desc(score))
dfss3 <- df3 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
#qqplot(dfss$score,dfss$fstat,xlab="adjacency score",ylab ="F statistics" )
dfss2 <- df2 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
#qqplot(dfss$score,dfss$fstat,xlab="adjacency score",ylab ="F statistics" )
dfss1 <- df1 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
#qqplot(dfss$score,dfss$fstat,xlab="adjacency score",ylab ="F statistics" )
#dfss0 <- df0 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
dfss_x1 <- map_dfr(list(dfss1,dfss2,dfss3),bind_rows,.id = "threshold") %>% mutate(exposure="x1")

dfgwas <- read_csv("dfgwas_x2.csv") %>% arrange(p)
#df0 <- read_csv("x2_score_1e-2.csv") %>% dplyr::arrange(desc(score))
df1 <- read_csv("x2_score_1_5000_1.csv") %>% dplyr::arrange(desc(score))
df2 <- read_csv("x2_score_1e-2_5000_1.csv") %>% dplyr::arrange(desc(score))
df3 <- read_csv("x2_score_1e-4_5000_1.csv") %>% dplyr::arrange(desc(score))
dfss3 <- df3 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
dfss2 <- df2 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
dfss1 <- df1 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
#dfss0 <- df0 %>% left_join(dfgwas,by="snp")# %>% slice_max(order_by = fstat,n=100)
dfss_x2 <- map_dfr(list(dfss1,dfss2,dfss3),bind_rows,.id = "threshold") %>% mutate(exposure="x2")

dfss <- dfss_x1 %>% bind_rows(dfss_x2) %>% mutate_at(vars(threshold),factor,labels=c("1","1e-2","1e-4"))

library(ggpmisc)
ggplot(dfss,aes(fstat,score))+
  facet_grid(threshold~exposure,scales = "free_y")+
  geom_point()+
  geom_smooth(method = "lm", se=TRUE)+
  stat_poly_eq(aes(label = paste(after_stat(rr.label),after_stat(p.value.label), sep = "*\", \"*")), 
               parse = TRUE)+
  theme_few()+
  theme(axis.text.y = element_text(angle=0, hjust=1, vjust=.5),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  xlab("F statistics")+
  ylab("adjacency score")
ggsave("com.jpg", width=7.5, height=7.6, dpi = 300)

########################
df <- read_csv("epistas_sum.csv")
df[which(df$algorithm=="iamb"),"algorithm"] <- "fast.iamb"
df1 <- df%>% 
  pivot_longer(cols = -c(algorithm,time), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(algorithm),factor,levels=c("hc","tabu","pc.stable","gs","fast.iamb","mmhc","rsmax2")) %>%
  mutate_at(vars(indices),factor,levels=c("pr_50","pr_25","auc"))
ggplot(df1)+
  geom_col(aes(x=algorithm,y=time/5),width=0.6,fill=8,color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=algorithm,y=performance,color=indices,group=indices))+
  geom_line(aes(x=algorithm,y=performance,color=indices,group=indices))+
  # geom_col(aes(x=algorithm,y=performance,fill=indices),position = position_dodge(preserve = 'single'),color="black")+
  xlab("learning method")+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_stata(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*5,name = "running time (min)"))+
  ggtitle(TeX("$X_1$"))
  # geom_point(aes(x=algorithm,y=time/5,group=1))+
  # geom_line(aes(x=algorithm,y=time/5,group=1),linewidth=1,color="black")
ggsave("epis1.png", width=5.5, height=5, dpi = 300)

df <- read_csv("epistas_sum_f2.csv")
df[which(df$algorithm=="iamb"),"algorithm"] <- "fast.iamb"
df2 <- df%>% 
  pivot_longer(cols = -c(algorithm,time), names_to = "indices", values_to = "performance") %>% 
  mutate_at(vars(algorithm),factor,levels=c("hc","tabu","pc.stable","gs","fast.iamb","mmhc","rsmax2")) %>%
  mutate_at(vars(indices),factor,levels=c("pr_50","pr_25","auc"))
ggplot(df)+
  geom_col(aes(x=algorithm,y=time/5),width=0.6,fill=8,color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=algorithm,y=performance,color=indices,group=indices))+
  geom_line(aes(x=algorithm,y=performance,color=indices,group=indices))+
  # geom_col(aes(x=algorithm,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  xlab("learning method")+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank())+
  scale_color_stata(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*5,name = "running time (min)"))+
  ggtitle(TeX("$X_2$"))
  # geom_point(aes(x=algorithm,y=time/5,group=1))+
  # geom_line(aes(x=algorithm,y=time/5,group=1),linewidth=1,color="black")
ggsave("epis2.png", width=5.5, height=5, dpi = 300)

df <- bind_rows(df1,df2) %>% mutate(exposure=rep(c("X1","X2"),each=nrow(df1)))
ggplot(df)+facet_grid(cols=vars(exposure))+
  #geom_col(aes(x=algorithm,y=time/5),width=0.6,fill=8,color="black",position = position_dodge(preserve = "total"))+
  geom_point(aes(x=algorithm,y=performance,color=indices,group=indices))+
  geom_line(aes(x=algorithm,y=performance,color=indices,group=indices))+
  # geom_col(aes(x=algorithm,y=performance,fill=indices),position = position_dodge(preserve = 'single'))+
  xlab("learning method")+
  ylab("performance")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank(),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+
  scale_color_stata(labels = c("TDR in top\n50 variants", "TDR in top\n25 variants","AUC"))+
  coord_cartesian(ylim = c(0,1))+
  #scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*5,name = "running time (min)"))+
  ggtitle("Epistas")
  
ggsave("epis.png", width=6, height=5, dpi = 300)

################
df <- read_delim("../application2211_binary/pheno.txt")

df <- df%>%mutate_at(vars(SCZ),factor,labels=c("non-SCZ","SCZ")) %>% mutate_at(vars(PTSD),factor,labels=c("non-PTSD","PTSD"))
ggplot(df)+
  facet_grid(cols=vars(SCZ))+
  geom_histogram(aes(lymphocyte,..density..,fill=SCZ),bins=30)+
  geom_freqpoly(aes(lymphocyte,..density..,color=SCZ),bins=30,size=0.7)+
  theme_few()+
  scale_color_gdocs()
ggplot(df)+
  facet_grid(cols=vars(PTSD))+
  geom_histogram(aes(lymphocyte,..density..,fill=PTSD),bins=30)+
  geom_freqpoly(aes(lymphocyte,..density..,color=PTSD),bins=30,size=0.7)+
  theme_few()+
  scale_color_gdocs()


df1 <- df%>%filter(SCZ==1|(SCZ==0&PTSD==0))%>%
  mutate(outcome=if_else(SCZ==1,1,0))%>%
  mutate_at(vars(outcome),factor,labels=c("Control","Case")) %>%
  mutate(label="SCZ") %>%
  select(lymphocyte,outcome,label)
df2 <- df%>%filter(PTSD==1|(SCZ==0&PTSD==0))%>%
  mutate(outcome=if_else(PTSD==1,1,0))%>%
  mutate_at(vars(outcome),factor,labels=c("Control","Case")) %>%
  mutate(label="PTSD") %>%
  select(lymphocyte,outcome,label)


dfc <- bind_rows(df1,df2)
ggplot(dfc)+
  facet_grid(col=vars(label))+
  geom_boxplot(aes(outcome,lymphocyte,fill=outcome),outlier.colour = NULL,outlier.size = 1)+
  geom_signif(aes(outcome,lymphocyte),comparisons = list(c("Control","Case")),map_signif_level = T,test = t.test)+
  theme_few()+
  xlab("group")+
  scale_fill_stata()+
  theme(legend.title = element_blank(),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))
ggsave("app.png", width=6.5, height=5, dpi = 300)

library(ggsignif)
ggplot(df)+
  geom_boxplot(aes(SCZ,lymphocyte,fill=SCZ),outlier.colour = NULL,outlier.size = 1)+
  geom_signif(aes(SCZ,lymphocyte),comparisons = list(c("non-SCZ","SCZ")),map_signif_level = T,test = t.test)+
  theme_few()+
  xlab("group")+
  scale_fill_stata()
ggsave("app1.png", width=5.5, height=5, dpi = 300)
ggplot(df)+
  geom_boxplot(aes(PTSD,lymphocyte,fill=PTSD),outlier.colour = NULL,outlier.size = 1)+
  geom_signif(aes(PTSD,lymphocyte),comparisons = list(c("non-PTSD","PTSD")),map_signif_level = T,test = t.test)+
  theme_few()+
  xlab("group")+
  scale_fill_stata()
ggsave("app2.png", width=5.5, height=5, dpi = 300)



#############3
df <- read_csv("compareIter.csv") %>% mutate(bia=est-2) %>% mutate_at(vars(iter),as.factor)
ggplot(df)+
  geom_point(aes(iter,bia))+
  geom_line(aes(iter,bia,group=1))+
  geom_errorbar(aes(iter,bia,ymin=L-2,ymax=U-2), width=.2)+
  xlab("Iteration")+
  ylab("Bias")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  scale_fill_aaas()+
  coord_cartesian(ylim = c(-0.25,0.25))
ggsave("sens1.png", width=3.5, height=3, dpi = 300)

ggplot(df)+geom_point(aes(iter,R))+
  geom_line(aes(iter,R,group=1))

df <- read_csv("compareMCMC_5000.csv")
ggplot(df)+
  geom_point(aes(prior,bias))+
  geom_line(aes(prior,bias,group=1))+
  geom_errorbar(aes(prior,bias,ymin=low-2,ymax=up-2), width=.2)+
  xlab("Prior")+
  ylab("Bias")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+
  scale_fill_aaas()+
  coord_cartesian(ylim = c(-0.25,0.25))
ggsave("sens2.png", width=4, height=3, dpi = 300)

ggplot(df)+geom_point(aes(iter,R))+
  geom_line(aes(iter,R,group=1))

df <- read_csv("compareIV.csv") %>% mutate(bias=est-2) %>% mutate_at(vars(IVnum),as.factor)
ggplot(df)+
  geom_point(aes(IVnum,bias))+
  geom_line(aes(IVnum,bias,group=1))+
  geom_errorbar(aes(IVnum,bias,ymin=L-2,ymax=U-2), width=.2)+
  xlab("IV number")+
  ylab("Bias")+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.99,0.995), legend.justification=c(0.99,0.995),
        strip.text.x = element_text(face = "bold.italic"),
        strip.background = element_rect(fill="grey"))+
  scale_color_aaas()+
  coord_cartesian(ylim = c(-0.25,0.25))
ggsave("sens3.png", width=3.5, height=3, dpi = 300)

ggplot(df)+geom_point(aes(IVnum,R))+
  geom_line(aes(IVnum,R,group=1))

library(pROC)
df1 <- read_csv("x_score_p10000_3.csv")
df2 <- read_csv("x_score_p20000_3.csv")
df3 <- read_csv("x_score_p50000_3.csv")
df4 <- read_csv("x_score_p100000_3.csv")
roc_object1 <- roc(df1$label,df1$score)
roc_object2 <- roc(df2$label,df2$score)
roc_object3 <- roc(df3$label,df3$score)
roc_object4 <- roc(df4$label,df4$score)
# auc_df <- data.frame(auc = c(auc(roc_object1), auc(roc_object2), auc(roc_object3), auc(roc_object4)),
#                      roc_label = c("ROC1", "ROC2", "ROC3", "ROC4"))
# ggroc(list(roc_object1,roc_object2,roc_object3,roc_object4),legacy.axes = T)+theme_few()+
#   geom_segment(aes(x=0,xend=1,y=0,yend=1),color="grey",linetype=4)+
#   annotate("text",x=0.75,y=0.25,label=paste("AUC = ",round(roc_object1$auc,3)))+
#   annotate("text",x=0.75,y=0.2,label=paste("AUC = ",round(roc_object2$auc,3)))

roc_data <- rbind(
  data.frame(spe = roc_object1$specificities, sen = roc_object1$sensitivities, dataset = "p = 10000", stringsAsFactors = FALSE),
  data.frame(spe = roc_object2$specificities, sen = roc_object2$sensitivities, dataset = "p = 20000", stringsAsFactors = FALSE),
  data.frame(spe = roc_object3$specificities, sen = roc_object3$sensitivities, dataset = "p = 50000", stringsAsFactors = FALSE),
  data.frame(spe = roc_object4$specificities, sen = roc_object4$sensitivities, dataset = "p = 100000", stringsAsFactors = FALSE)
)

# Plot ROC curves with annotations
library(ggsci)
ggplot(roc_data, aes(x = 1-spe, y = sen, color = dataset)) +
  geom_line() +
  scale_color_aaas()+
  theme_classic() +
  # scale_color_manual(values = c("red", "blue", "green", "orange")) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  annotate("text", x = 0.7, y = 0.3, label = paste0("p = 10000, AUC = ", sprintf("%0.3f",auc(roc_object1)))) +
  annotate("text", x = 0.7, y = 0.25, label = paste0("p = 20000, AUC = ", sprintf("%0.3f",auc(roc_object2)))) +
  annotate("text", x = 0.7, y = 0.2, label = paste0("p = 50000, AUC = ", sprintf("%0.3f",auc(roc_object3)))) +
  annotate("text", x = 0.7, y = 0.15, label = paste0("p = 100000, AUC = ", sprintf("%0.3f",auc(roc_object4))))+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),color="grey",linetype=4)
ggsave("g_p.png", width=5.5, height=5, dpi = 300)
