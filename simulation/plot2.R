library(tidyverse)
df1 <- read_csv("C:\\Users\\Administrator\\Desktop\\sim.csv")
df <- df1[,-1]
df <- df %>% pivot_longer(cols = -c("pleio","id","item"),names_to = "method")
result <- df %>%
  group_by(id, method) %>%
  summarise(
    value_mse = sum(value[item %in% c("bias", "se")]^2)
  ) %>%
  mutate(item = "mse")


library(ggsci)
library(ggthemes)

df <- df[-which(df$method=="TSHT"),]
df[which(df$item=="se"),"item"] <- "variance"
df[which(df$item=="bias"),"item"] <- "bias2"


# df1 <- df %>% filter(index %in% seq(1,7)) %>% mutate(id = case_when(
#   index==1~0,
#   index==2~1,
#   index==3~2,
#   index==4~3,
#   index==5~4,
#   index==6~5,
#   index==7~6
# )) %>%
#   mutate_at(vars(id),as.factor) %>%
#   mutate_at(vars(item),as.factor)
g1 <- df %>% 
  mutate(value = value^2) %>%
  mutate_at(vars(id),factor,
            labels=c("0+0","50+0","25+25","0+50","100+0","50+50","0+100")) %>%
  group_by(id, method, pleio) %>% 
  mutate(mse = cumsum(value)) %>% 
  ggplot(aes(id, mse, fill = method)) + 
  facet_wrap(~ pleio, nrow = 2, 
             labeller = labeller(pleio = c(balanced = "Balanced Pleiotropy", directional = "Directional Pleiotropy"))) +
  geom_col(data = . %>% filter(item=="bias2"), position = position_dodge(width = 0.9), alpha = 1,color="black") +
  geom_col(data = . %>% filter(item=="variance"), position = position_dodge(width = 0.9), alpha = 0.4,color="black") +
  geom_tile(aes(y=NA_integer_, alpha = item)) + 
  scale_alpha_manual(values = c(1,0.4)) +
  theme_classic() +
  scale_fill_stata() +
  theme(axis.text = element_text(size = 12),  # 调整刻度标签的文字大小
        axis.title = element_text(size = 14,face="bold"),
        strip.text = element_text(size = 14,face = "bold.italic"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 11)) +
  xlab("Scenarios (#{True pleiotropic variants} + #{Correlated pleiotropic variants})") +
  ylab("MSE")

setwd("C:\\Users\\Administrator\\Desktop")
ggsave("1.tiff",g1,width=13,height=10)



df1 <- read_csv("C:\\Users\\Administrator\\Desktop\\sim2.csv")
df <- df1 %>% pivot_longer(cols = -c("pleio","id","item"),names_to = "method")
df[which(df$item=="se"),"item"] <- "variance"
df[which(df$item=="bias"),"item"] <- "bias2"
g2 <- df %>% 
  mutate(value = value^2) %>%
  mutate_at(vars(id),factor,
            labels=c("0+0","50+0","25+25","0+50","100+0","50+50","0+100")) %>%
  mutate_at(vars(method),factor,levels=c("BNMR","BMR")) %>%
  group_by(id, method, pleio) %>% 
  mutate(mse = cumsum(value)) %>% 
  ggplot(aes(id, mse, fill = method)) + 
  facet_wrap(~ pleio, nrow = 2, 
             labeller = labeller(pleio = c(balanced = "Balanced Pleiotropy", directional = "Directional Pleiotropy"))) +
  geom_col(data = . %>% filter(item=="bias2"), position = position_dodge(width = 0.9), alpha = 1,color="black") +
  geom_col(data = . %>% filter(item=="variance"), position = position_dodge(width = 0.9), alpha = 0.4,color="black") +
  geom_tile(aes(y=NA_integer_, alpha = item)) + 
  scale_alpha_manual(values = c(1,0.4)) +
  theme_classic() +
  scale_fill_stata() +
  theme(axis.text = element_text(size = 12),  # 调整刻度标签的文字大小
        axis.title = element_text(size = 14,face="bold"),
        strip.text = element_text(size = 14,face = "bold.italic"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 11)) +
  xlab("Scenarios (#{True pleiotropic variants} + #{Correlated pleiotropic variants})") +
  ylab("MSE")

setwd("C:\\Users\\Administrator\\Desktop")
ggsave("2.tiff",g2,width=8,height=10)



library(tidyverse)
library(readxl)
library(forestplot)
df <- read_excel("C:\\Users\\Administrator\\Desktop\\res.xlsx")
#df <- df[9:11,]
df <- df%>%mutate(name=paste(exposure,outcome,sep = "→"))

df <- df%>%mutate(CI=sprintf("%.2f (%.2f,%.2f)",mean,low,up))

# png(filename="fore.png",width = 8.00,height = 4, units = "in", res = 300)
p1 <- forestplot::forestplot(labeltext=as.matrix(df[,c("name","CI")]),
           mean=df$mean,
           lower=df$low,
           upper=df$up,
           zero=0,
           boxsize=0.1,
           lineheight = unit(7,'mm'),
           colgap=unit(7,'mm'),
           lwd.zero=1.5,
           lwd.ci=1.5,
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = 'red'),
           xlab="beta",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 1),
           lty.ci = "solid",
           title = "Strict pre-filtering P + BNMR",
           line.margin = 1,
           graph.pos=3,
           xticks=seq(-1,11))
# dev.off()

df <- read_excel("C:\\Users\\Administrator\\Desktop\\res1.xlsx")
#df <- df[9:11,]
df <- df%>%mutate(name=paste(exposure,outcome,sep = "→"))

df <- df%>%mutate(CI=sprintf("%.2f (%.2f,%.2f)",mean,low,up))
p2 <- forestplot::forestplot(labeltext=as.matrix(df[,c("CI")]),
                 mean=df$mean,
                 lower=df$low,
                 upper=df$up,
                 zero=0,
                 boxsize=0.1,
                 lineheight = unit(7,'mm'),
                 colgap=unit(7,'mm'),
                 lwd.zero=1.5,
                 lwd.ci=1.5,
                 col=fpColors(box='#458B00',
                              summary='#8B008B',
                              lines = 'black',
                              zero = 'red'),
                 xlab="beta",
                 lwd.xaxis =1,
                 txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                                  xlab  = gpar(cex = 0.8),
                                  cex = 1),
                 lty.ci = "solid",
                 title = "Loose pre-filtering P with LD clumping + BNMR",
                 line.margin = 1,
                 graph.pos=2,
                 xticks=seq(-1,11))

library(gridExtra)
# grob1 <- ggplotGrob(p1)
# grob2 <- ggplotGrob(p2)
library(ggplotify)
library(patchwork)

pr1 <- grid.grabExpr(print(p1))
pr2 <- grid.grabExpr(print(p2))
combined_plot <- gridExtra::grid.arrange(pr1, pr2,nrow=1)
ggsave("C:\\Users\\Administrator\\Desktop\\p1.tiff",combined_plot,width=12,height=6,dpi=300)

# p1 <- grid2grob(print(p1))
# p2 <- grid2grob(print(p2))
# combined_plot <- grid.arrange(p1, p2)
# arrangeGrob(p1, p2)

library(patchwork)
p1+p2

# my_plot <- ggplot(df, aes(x = mean, y = name)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = low, xmax = up), height = 0.2) +
#   geom_text(aes(label = paste0(mean, " (", low, ", ", up, ")")), hjust = -0.2)



df <- read_excel("C:\\Users\\Administrator\\Desktop\\lymm.xlsx")
df1 <- df%>%filter(method=="BNMR")
df1 <- df1%>%mutate(name=paste(exposure,outcome,sep = "→"))
df1 <- df1%>%mutate(CI=sprintf("%.2f (%.2f,%.2f)",mean,low,up))
p1 <- forestplot::forestplot(labeltext=as.matrix(df1[,c("name","CI")]),
                             mean=df1$mean,
                             lower=df1$low,
                             upper=df1$up,
                             zero=0,
                             boxsize=0.1,
                             lineheight = unit(7,'mm'),
                             colgap=unit(7,'mm'),
                             lwd.zero=1.5,
                             lwd.ci=1.5,
                             col=fpColors(box='#458B00',
                                          summary='#8B008B',
                                          lines = 'black',
                                          zero = 'red'),
                             xlab="beta",
                             lwd.xaxis =1,
                             txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                                              xlab  = gpar(cex = 0.8),
                                              cex = 1),
                             lty.ci = "solid",
                             title = "BNMR",
                             line.margin = 1,
                             graph.pos=3)

df2 <- df%>%filter(method=="MR-Egger")
df2 <- df2%>%mutate(name=paste(exposure,outcome,sep = "→"))
df2 <- df2%>%mutate(CI=sprintf("%.2f (%.2f,%.2f)",mean,low,up))
p2 <- forestplot::forestplot(labeltext=as.matrix(df2[,c("name","CI")]),
                             mean=df2$mean,
                             lower=df2$low,
                             upper=df2$up,
                             zero=0,
                             boxsize=0.1,
                             lineheight = unit(7,'mm'),
                             colgap=unit(7,'mm'),
                             lwd.zero=1.5,
                             lwd.ci=1.5,
                             col=fpColors(box='#458B00',
                                          summary='#8B008B',
                                          lines = 'black',
                                          zero = 'red'),
                             xlab="beta",
                             lwd.xaxis =1,
                             txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                                              xlab  = gpar(cex = 0.8),
                                              cex = 1),
                             lty.ci = "solid",
                             title = "MR-Egger",
                             line.margin = 1,
                             graph.pos=3)
pr1 <- grid.grabExpr(print(p1))
pr2 <- grid.grabExpr(print(p2))
combined_plot <- gridExtra::grid.arrange(pr1, pr2,nrow=2)
ggsave("C:\\Users\\Administrator\\Desktop\\p2_1.tiff",combined_plot,width=7,height=6,dpi=300)


df1 <- read_excel("C:\\Users\\Administrator\\Desktop\\res.xlsx") %>% 
  select(exposure,outcome,mean,low,up) %>%
  mutate(method="Strict pre-filtering P + BNMR")
df2 <- read_excel("C:\\Users\\Administrator\\Desktop\\res1.xlsx") %>% 
  select(exposure,outcome,mean,low,up) %>%
  mutate(method="Loose pre-filtering P with LD clumping + BNMR")
df <- rbind(df1,df2)
df <- df%>%mutate(name=paste(exposure,outcome,sep = "→"))
df <- df%>%mutate(CI=sprintf("%.2f (%.2f,%.2f)",mean,low,up))
tbex <- tibble(name=c(NA,unique(df$name)),CI=c(NA,df%>%group_by(name)%>%summarize(CI=str_c(CI,collapse = "\n"))%>%pull(CI)))
forestplot::forestplot(labeltext=as.matrix(tbex[,c("name","CI")]),
                       align = c("l","l","l","l"),
                       legend =  c( "Strict pre-filtering P + BNMR", 
                                    "Loose pre-filtering P with LD clumping + BNMR"),
                       # legend_args = fpLegend(pos = list(x=.4, y=0.02),
                       #                        gp=gpar(col="#CCCCCC", fill="#F9F9F9")),
                       mean=cbind(c(NA,
                                    df%>%filter(method=="Strict pre-filtering P + BNMR")%>%pull(mean)
                                    ),
                                  c(NA,
                                    df%>%filter(method=="Loose pre-filtering P with LD clumping + BNMR")%>%pull(mean))
                                  ),
                       lower=cbind(c(NA,
                                     df%>%filter(method=="Strict pre-filtering P + BNMR")%>%pull(low)
                                     ),
                                   c(NA,
                                     df%>%filter(method=="Loose pre-filtering P with LD clumping + BNMR")%>%pull(low))
                                   ),
                       upper=cbind(c(NA,
                                     df%>%filter(method=="Strict pre-filtering P + BNMR")%>%pull(up)
                                     ),
                                   c(NA,
                                     df%>%filter(method=="Loose pre-filtering P with LD clumping + BNMR")%>%pull(up))
                                   ),
                       zero=0,
                       boxsize=0.1,
                       # lineheight = unit(2,'cm'),
                       colgap=unit(7,'mm'),
                       lwd.zero=1.5,
                       lwd.ci=1.5,
                       # col=fpColors(box='#458B00',
                       #              summary='#8B008B',
                       #              lines = 'black',
                       #              zero = 'red'),
                       col=fpColors(box=c('gray', "black"), 
                                    lines = c('gray', "black"),
                                    zero = 'red'),
                       xlab="beta",
                       lwd.xaxis =1,
                       txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                                        xlab  = gpar(cex = 0.8),
                                        cex = 1),
                       lty.ci = "solid",
                       title = "MR-Egger",
                       line.margin = 1,
                       graph.pos=3)

# df |> mutate_at(vars(method),as.factor) |>
#   group_by(method) |> 
ls <- unique(df$name)
df <- df %>% group_by(name) %>% 
  group_modify(~{.x%>%mutate(C=.x%>%summarize(CI=str_c(CI,collapse = "\n"))%>%pull(CI))}) %>% 
  ungroup() %>%
  mutate_at(vars(name),factor,levels=ls) %>%
  arrange(name)
p <- df |> group_by(method) |> forestplot(clip = c(-1, 11),
             labeltext = c(name,C),
             mean=mean,
             lower=low,
             upper=up,
             zero=0,
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             xlab = "beta",
             lwd_Ci = 1.5,
             txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                              xlab  = gpar(cex = 1),
                              label = gpar(cex = 0.9),
                              legend = gpar(cex = 1),
                              cex = 1),
             xticks=seq(-1,11),
             col=fpColors(zero = 'black')) |>
  fp_add_lines("steelblue") |>
  # fp_add_header("Variable") |>
  fp_set_style(box = c("darkblue", "darkred") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) |>
  fp_set_zebra_style("#F5F9F9")
pr <- grid.grabExpr(print(p))
ggsave("C:\\Users\\Administrator\\Desktop\\p.tiff",pr,width=12,height=6,dpi=300)
