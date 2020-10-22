library(ggplot2)
library(easyGgplot2)
#library(extrafont)
library(latex2exp)
#font_import()
scaleFUN <- function(x) sprintf("%.1f", x)


data1 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsGreedyResult_2.txt",header = TRUE)
kk = dim(data1)[1];

df1 <- data.frame(
  k = c(1:kk),
  data1,
  type = 'Greedy'
)

data2 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsRandResult_2.txt",header = TRUE)
df2 <- data.frame(
  k = c(1:kk),
  data2,
  type = 'Random'
)

data3 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsMaxDResult_2.txt",header = TRUE)
df3 <- data.frame(
  k = c(1:kk),
  data3,
  type = 'Max-Degree'
)

df1.new = df1[seq(1, nrow(df1), 20), ]
df2.new = df2[seq(1, nrow(df2), 20), ]
df3.new = df3[seq(1, nrow(df3), 20), ]

df = rbind(df1.new, df2.new, df3.new)

plot1 <- ggplot(df,
       aes(x = k , y = Infections, colour=type, shape=type)) + theme_bw() +
      ylab("Approx of Expected Infections")  + scale_y_continuous(labels=scaleFUN) + xlab("k") + geom_line() + geom_point(size=3)+
      theme(legend.position = c(0.65,0.85),legend.background = element_rect(fill = "transparent",colour = NA),
            legend.title = element_blank(), legend.key.size = unit(4,'cm'), legend.key.height = unit(0.7,'cm'),
            legend.key.width = unit(2,'cm'), legend.text = element_text(size = 16,face="bold"),
            axis.text = element_text(size = 16,face="bold", color = 'black'), 
            axis.title.x = element_text(size = 20,face="bold"),
            axis.title.y = element_text(size = 20,face="bold")) +
      theme(panel.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA)) + 
      theme(panel.border = element_rect(size=1, colour = "black"))


data1 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsGreedyResultReal_2.txt",header = TRUE)
kk = dim(data1)[1];

df1 <- data.frame(
  k = c(1:kk),
  data1,
  type = 'Greedy'
)

data2 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsRandResultReal_2.txt",header = TRUE)
df2 <- data.frame(
  k = c(1:kk),
  data2,
  type = 'Random'
)

data3 <- read.table("D:\\LirenShan\\西北\\SIRcontrol\\SIRControl-main\\code\\results\\dsMaxDResultReal_2.txt",header = TRUE)
df3 <- data.frame(
  k = c(1:kk),
  data3,
  type = 'Max-Degree'
)

df1.new = df1[seq(1, nrow(df1), 20), ]
df2.new = df2[seq(1, nrow(df2), 20), ]
df3.new = df3[seq(1, nrow(df3), 20), ]

df = rbind(df1.new, df2.new, df3.new)

plot2 <- ggplot(df,
                aes(x = k , y = Infections, colour=type, shape=type)) + theme_bw() +
  ylab("Expected Infections")  + scale_y_continuous(labels=scaleFUN) + xlab("k") + geom_line() + geom_point(size=3)+
  theme(legend.position = c(0.65,0.85),legend.background = element_rect(fill = "transparent",colour = NA),
        legend.title = element_blank(), legend.key.size = unit(4,'cm'), legend.key.height = unit(0.7,'cm'),
        legend.key.width = unit(2,'cm'), legend.text = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16,face="bold", color = 'black'), 
        axis.title.x = element_text(size = 20,face="bold"),
        axis.title.y = element_text(size = 20,face="bold")) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) + 
  theme(panel.border = element_rect(size=1, colour = "black"))

# 
ggplot2.multiplot(plot1,plot2, cols=2)
# 
library(Cairo)
cairo_ps("D-SIR_ER.eps", width = 11.50, height = 4.50)
ggplot2.multiplot(plot1,plot2, cols=2)
dev.off()



