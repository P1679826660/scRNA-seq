

#install.packages("ggplot2")


library(ggplot2)           
logFCfilter=0.5             
adj.P.Val.Filter=0.05                #有要求的话，改阈值
inputFile="all.txt"        
setwd("C:\\Users\\X\\Desktop")       #这里修改地址


rt=read.table(inputFile, header=T, sep="\t", check.names=F)

Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
rt=cbind(rt, Sig=Sig)

summary(rt$logFC)  
summary(-log10(rt$adj.P.Val))  


p = ggplot(rt, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(col = Sig)) +
  scale_color_manual(values = c("#1976D2", "black", "red")) +
  xlim(-6, 6) +
  labs(title = " ") +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), col = "#8E24AA", linewidth = 1, linetype = 2) +  # Replace cex with linewidth
  geom_hline(yintercept = -log10(adj.P.Val.Filter), col = "#8E24AA", linewidth = 1, linetype = 2) +  # Replace cex with linewidth
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

p=p+theme_bw()


pdf(file="volcano.pdf", width=6, height=5.1)
print(p)
dev.off()

