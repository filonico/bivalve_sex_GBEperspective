library(dplyr)
library(viridis)
library(phangorn)
library(ggpubr)

##################################################################
### GENERATE THE PLOT WITH PRESENCE/ABSENCE DATA ON DMRT GENES ###
##################################################################

# import and format data
compilation.raw <- read.table("dmrt.compilation.tsv", header = TRUE, check.names = FALSE)
colnames(compilation.raw)[5] <- "Dmrt4/5"

compilation.tibble <- dplyr::as_tibble(compilation.raw) %>%
  tidyr::pivot_longer(!Species, names_to = "Gene", values_to = "Count")

compilation.tibble$Gene["Dmrt45"]

# Species_order <- factor(compilation.tibble$Species, levels = c("Mmar","Mner","Pstr","Scon","Dpol","Dros","Amar","Csin","Mmer","Rphi","Lfor","Bpla","Mphi","Pvir","Mcor","Mgal","Medu","Sbro","Pyes","Pmax","Apur","Airi","Airc","Pmar","Pfum","Pfuc","Sglo","Cvir","Cari","Cgig","Chon"))
#species ored 4evolution
Species_order <- factor(compilation.tibble$Species, levels = c("Mmar","Pstr","Scon","Dros","Dpol","Amar","Csin","Rphi","Mmer","Sbro","Pyes","Pmax","Apur","Airi","Airc","Pmar","Sglo","Cvir","Cari","Cgig","Pvir","Mcor","Mgal","Medu"))
gene_order <- factor(compilation.tibble$Gene, levels = c("Dmrt2","Dmrt3","Dmrt4/5","Dmrt1L"))
                     
plot.compilation <- ggplot(compilation.tibble, aes(x = Species_order,
                                                   y = fct_rev(gene_order),
                                                   fill = factor(Count))) +
  geom_tile(color = "grey20", size = 0.4) +
  geom_text(aes(label = Count), size = 2) +
  scale_fill_manual(values = c("white", "#D6DBE2", "#9CA8BA", "#58677E")) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  labs(fill = "Count") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 0),
        axis.text.y = element_text(face = 'bold.italic',
                                   colour = c("#0092ff", "#bbc5d0","#bbc5d0","#bbc5d0")),
        panel.border = element_rect(size = 1),
        # axis.text.y = element_text(size = 5),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_equal(expand = 0)

plot.compilation


#boxplot
Dmrt1L.alignment <- read.phyDat("Dmrt1L.aligned.trim06.faa", format = "fasta", type = "AA")
Dmrt1L.dist <- dist.ml(Dmrt1L.alignment, model = "JTT")
Dmrt1L.dist.list <- Dmrt1L.dist[1:length(Dmrt1L.dist)]

Dmrt2.alignment <- read.phyDat("Dmrt2.aligned.trim06.faa", format = "fasta", type = "AA")
Dmrt2.dist <- dist.ml(Dmrt2.alignment, model = "JTT")
Dmrt2.dist.list <- Dmrt2.dist[1:length(Dmrt2.dist)]

Dmrt3.alignment <- read.phyDat("Dmrt3.aligned.trim06.faa", format = "fasta", type = "AA")
Dmrt3.dist <- dist.ml(Dmrt3.alignment, model = "JTT")
Dmrt3.dist.list <- Dmrt3.dist[1:length(Dmrt3.dist)]

Dmrt45.alignment <- read.phyDat("Dmrt45.aligned.trim06.faa", format = "fasta", type = "AA")
Dmrt45.dist <- dist.ml(Dmrt45.alignment, model = "JTT")
Dmrt45.dist.list <- Dmrt45.dist[1:length(Dmrt45.dist)]

dist.al.distr <- data.frame(c(rep("Dmrt2",length(Dmrt2.dist.list)),
                              rep("Dmrt3",length(Dmrt3.dist.list)),
                              rep("Dmrt4/5",length(Dmrt45.dist.list)),
                              rep("Dmrt1L",length(Dmrt1L.dist.list))),
                            c(Dmrt2.dist.list,
                              Dmrt3.dist.list,
                              Dmrt45.dist.list,
                              Dmrt1L.dist.list))

colnames(dist.al.distr) <- c("Gene","Dist")

dist.al.distr$Gene <- factor(dist.al.distr$Gene, levels = c("Dmrt2", "Dmrt3", "Dmrt4/5", "Dmrt1L"))

#statistics related to boxplots
dist.al.distr.lm <- lm(Dist ~ Gene, data = dist.al.distr)
qqnorm(resid(dist.al.distr.lm))
shapiro.test(resid(dist.al.distr.lm)) #not normally-distributed

kruskal.test(dist.al.distr$Dist, dist.al.distr$Gene) #there are significative differences

pairwise.wilcox.test(dist.al.distr$Dist, dist.al.distr$Gene, p.adjust.method = "bonferroni")


boxplot <- ggplot(dist.al.distr, aes(Gene, Dist,
                                     fill = as.factor(Gene), color = as.factor(Gene))) +
  geom_boxplot(lwd = 1) +
  # ylim(0, 3.5) +
  # geom_segment(aes(x = "Dmrt1L", xend = "Dmrt2", y = 2.8, yend = 2.8), lwd = 0.7, col = "#7f878f") +
  # geom_segment(aes(x = "Dmrt1L", xend = "Dmrt3", y = 3, yend = 3), lwd = 0.7, col = "#7f878f") +
  # geom_segment(aes(x = "Dmrt1L", xend = "Dmrt4/5", y = 3.2, yend = 3.2), lwd = 0.7, col = "#7f878f") +
  # geom_segment(aes(x = "Dmrt3", xend = "Dmrt4/5", y = 2.6, yend = 2.6), lwd = 0.7, col = "#7f878f") +
  # geom_segment(aes(x = "Dmrt2", xend = "Dmrt4/5", y = 2.4, yend = 2.4), lwd = 0.7, col = "#7f878f") +
  # annotate("text", label = "****", x = 2.5, y = 2.85, col = "#7f878f", size = 4) +
  # annotate("text", label = "****", x = 3, y = 3.05, col = "#7f878f", size = 4) +
  # annotate("text", label = "****", x = 3.5, y = 3.25, col = "#7f878f", size = 4) +
  # annotate("text", label = "****", x = 2.5, y = 2.65, col = "#7f878f", size = 4) +
  # annotate("text", label = "***", x = 2, y = 2.45, col = "#7f878f", size = 4) +
  scale_fill_manual(values = alpha(c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff"), 0.5)) +
  scale_color_manual(values = c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff")) +
  ylab("Amino acid distance") +
  xlab("") +
  theme_light() +
  theme(axis.title.y = element_text(size = 13, vjust = 2.5, color = "#4d4d4d"),
        axis.text.x = element_text(size = 13, vjust = -0.5, face = 'bold.italic',
                                   colour = c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff")),
        legend.position = "none")

boxplot


panel <- ggarrange(plot.compilation, boxplot,
                   ncol = 1, nrow = 2, heights = c(0.4,1), align = "v",
                   labels = c("AB","C"))

panel

ggsave("plot.panel.4evolution.pdf", plot = panel, device = "pdf",
       dpi = 300, height = 6, width = 5.5, units = ("in"))

