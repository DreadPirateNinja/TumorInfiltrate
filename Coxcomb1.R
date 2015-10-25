df1 <- read.csv("~/R-tables/JL-M128_Counts.csv")

#Subsets matrix values
df1[,1] <- NULL
Other <- df1$CD45.Cells - rowSums(df1[,2:ncol(df1)])
df2 <- cbind(df1, Other)
colnames(df2) <- c("CD45.Cells", "B cells", "CD4 T-cells", "CD8 T-cells", "TAMs", "Monocytes", "Other")

#Normalizes table to CD45 count
df3 <- sweep(df2, 1, df2$CD45.Cells, "/")
df3$CD45.Cells <- NULL
df3 <- df3 * 100

#Define Treatment Groups
Tx <- rep(c("KO", "CTRL"),times=c(3,4))
df4 <- cbind(Tx, df3)

#Reshape table
library(reshape2)
df5 <- melt(df4, id.vars="Tx")
colnames(df5)<- c("Tx", "Subset", "Count")

#Summarize table
library(plyr)
df6 <- ddply(df5, .(Tx, Subset), summarise,mean=mean(Count),sd=sd(Count))
  
#######################################################################################################################
#Makes Coxcomb plot

#Splits data frame by treatment group, appends fold change from CTRL
df6_CTRL <- df6[ which(df6$Tx=='CTRL'),]
df6_KO <- df6[ which(df6$Tx=='KO'),]
df6_CTRL_Fold <- df6_CTRL$mean/df6_CTRL$mean
df6_KO_Fold <- df6_KO$mean/df6_CTRL$mean

cbind(df6_CTRL, df6_CTRL_Fold)
cbind(df6_KO, df6_KO_Fold)


#Overlays graphs
library (ggplot2)
v <- df6_CTRL$mean
w <- df6_KO$mean
pos1 <- 0.5 * (cumsum(v) + cumsum(c(0, v[-length(v)])))
p <- ggplot(df6_CTRL, aes(x = pos1)) +
  ylim(0,2.3) + 
  geom_vline(xintercept = pos1, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_bar(aes(y = df6_KO_Fold, fill = Subset), stat = "identity",  width = v, alpha = 0.7, color = "black") +
  scale_fill_brewer(palette="Spectral") +
  geom_bar(aes(y = df6_CTRL_Fold, fill = Subset), stat = "identity",  width = v, alpha = 0.3, color = "black", linetype = "dashed") +
  coord_polar(theta = "x") + scale_x_continuous(labels = df6_CTRL$Subset, breaks = pos1) +
  guides(fill = guide_legend(title="CD45+ Subsets", override.aes = list(colour = NULL))) +
  geom_text(aes(label = "0.5", y = 0.5, x = 0), color = "black") +
  geom_text(aes(label = "1.0", y = 1, x = 0), color = "black") +
  geom_text(aes(label = "1.5", y = 1.5, x = 0), color = "black") +
  geom_text(aes(label = "2.0", y = 2, x = 0), color = "black") +
  geom_hline(yintercept = 1, size = .75) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16),
    legend.key = element_rect(colour = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_blank())

################################################################################################################################
#Makes spie chart (p) of df6_CTRL
v <- df6_CTRL$mean
pos1 <- 0.5 * (cumsum(v) + cumsum(c(0, v[-length(v)])))
p <- ggplot(df6_CTRL, aes(x = pos1)) + 
  geom_bar(aes(y = df6_CTRL_Fold), stat = "identity",  width = v, color = "black") +
  coord_polar(theta = "x") + scale_x_continuous(labels = df6_CTRL$Subset, breaks = pos1)

#Makes spie chart(q) of df6_KO
w <- df6_KO$mean
pos2 <- 0.5 * (cumsum(w) + cumsum(c(0, w[-length(w)])))
q <- ggplot(df6_KO, aes(x = pos1)) + 
  geom_bar(aes(y = df6_KO_Fold), stat = "identity",  width = w, color = "black") +
  coord_polar(theta = "x") + scale_x_continuous(labels = df6_CTRL$Subset, breaks = pos1)

################################################################################################################################

#Makes stacked bar plot

ggplot(df6, aes(x = Tx)) +
  geom_bar(aes(y = mean, fill = Subset), stat = "identity") +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), stat = "identity", width = 0.5) +
  scale_fill_brewer(palette="Spectral") +
  ylab("Percent of CD45+ Immune cells") + xlab(NULL) + 
  guides(fill = guide_legend(title="CD45+ Subsets", override.aes = list(colour = NULL))) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size=16),
    axis.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=16),
    legend.title = element_text(size=18))

#########################################################################################
#Makes radar plot
library(ggplot2)
ggplot(df6, aes(x = Subset, y = mean)) + 
  geom_bar(color = "black", stat = "identity",  width = 1)
