library(phangorn)
library(tidyverse)


#read amino acid alignments
dmrt1L.alignment <- read.phyDat("dmrt1L.bivalve.aligned.trim04.faa", format = "fasta", type = "AA")
dmrt2.alignment <- read.phyDat("dmrt2.bivalve.aligned.trim04.faa", format = "fasta", type = "AA")
dmrt3.alignment <- read.phyDat("dmrt3.bivalve.aligned.trim04.faa", format = "fasta", type = "AA")
dmrt45.alignment <- read.phyDat("dmrt45.bivalve.aligned.trim04.faa", format = "fasta", type = "AA")

#compute pairwise amino acid distances
dmrt1L.dist <- dist.ml(dmrt1L.alignment, model = "JTT")
dmrt2.dist <- dist.ml(dmrt2.alignment, model = "JTT")
dmrt3.dist <- dist.ml(dmrt3.alignment, model = "JTT")
dmrt45.dist <- dist.ml(dmrt45.alignment, model = "JTT")

#get the list of computed pairwise amino acid distance
dmrt1L.dist.list <- dmrt1L.dist[1:length(dmrt1L.dist)]
dmrt2.dist.list <- dmrt2.dist[1:length(dmrt2.dist)]
dmrt3.dist.list <- dmrt3.dist[1:length(dmrt3.dist)]
dmrt45.dist.list <- dmrt45.dist[1:length(dmrt45.dist)]

#prepare the dataframe for the boxplot and statistical tests
dist.al.distr <- data.frame(c(rep("DmrtL1",length(dmrt1L.dist.list)),
                              rep("Other Dmrts",length(dmrt2.dist.list)),
                              rep("Other Dmrts",length(dmrt3.dist.list)),
                              rep("Other Dmrts",length(dmrt45.dist.list))),
                            c(dmrt1L.dist.list,
                              dmrt2.dist.list,
                              dmrt3.dist.list,
                              dmrt45.dist.list))

colnames(dist.al.distr) <- c("V1","V2")

#test for normality
dist.al.distr.lm <- lm(V2 ~ V1, data = dist.al.distr)
qqnorm(resid(dist.al.distr.lm))
shapiro.test(resid(dist.al.distr.lm))

#KW test
kruskal.test(dist.al.distr$V2, dist.al.distr$V1)

#Wilcoxon test
wilcox.test(dist.al.distr$V2~dist.al.distr$V1)

#generate the boxplot
dist.al.plot <- ggplot(dist.al.distr, aes(as.factor(V1), V2, fill = as.factor(V1), color = as.factor(V1))) +
  geom_boxplot(lwd = 1.3) +
  ylim(0, 3.33) +
  geom_segment(aes(x = "DmrtL1", xend = "Other Dmrts", y = 3.26, yend = 3.26), lwd = 1, col = "#7f878f") +
  annotate("text", label = "****", x = 1.5, y = 3.33, col = "#7f878f", size = 6) +
  scale_fill_manual(values = alpha(c("#3174bc","#bbc5d0"), 0.5)) +
  scale_color_manual(values = c("#3174bc","#bbc5d0")) +
  ylab("Amino acid distance") +
  xlab("") +
  theme_light() +
  theme(axis.title.y = element_text(size = 13, vjust = 2.5, color = "#4d4d4d"),
        axis.text.x = element_text(size = 13, vjust = -0.5),
        legend.position = "none")

ggsave("dist.al.plot.pdf", dist.al.plot, width = 3.5, height = 6, units = "in", dpi = 400)
