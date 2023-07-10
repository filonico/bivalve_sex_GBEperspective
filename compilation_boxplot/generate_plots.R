library(dplyr)
library(ggplot2)

##################################################################
### GENERATE THE PLOT WITH PRESENCE/ABSENCE DATA ON DMRT GENES ###
##################################################################

# import and format data
compilation.raw <- read.table("dmrt.compilation.tsv", header = TRUE, check.names = FALSE)
colnames(compilation.raw)[5] <- "Dmrt4/5"

compilation.tibble <- as_tibble(compilation.raw) %>%
  tidyr::pivot_longer(!Species, names_to = "Gene", values_to = "Count")

# define the species and gene order that should be displayed in the plot
species_order <- factor(compilation.tibble$Species, levels = c("Mmar","Pstr","Scon","Dros","Dpol","Amar","Csin","Rphi","Mmer","Sbro","Pyes","Pmax","Apur","Airi","Airc","Pmar","Sglo","Cvir","Cari","Cgig","Pvir","Mcor","Mgal","Medu"))
gene_order <- factor(compilation.tibble$Gene, levels = c("Dmrt1L","Dmrt4/5","Dmrt3","Dmrt2"))
                       
plot.compilation <- ggplot(compilation.tibble,
                           aes(x = species_order, y = gene_order, fill = factor(Count))) +
  
  # plot tile and gene counts and scale colors
  geom_tile(color = "grey20", linewidth = 0.4) +
  geom_text(aes(label = Count), size = 2) +
  scale_fill_manual(values = c("white", "#D6DBE2")) +
  
  theme_bw() +
  scale_x_discrete(position = "top") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 0),
        axis.text.y = element_text(face = 'bold.italic', colour = c("#0092ff", "#bbc5d0","#bbc5d0","#bbc5d0")),
        panel.border = element_rect(linewidth = 1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_equal(expand = 0)

plot.compilation

##############################################################################
### GENERATE THE BOXPLOT ABOUT PAIRWISE AMINO ACID DISTANCES OF DMRT GENES ###
##############################################################################

# define a function to calculate pairwise amino acid distances of each Dmrt gene
calculate_distance <- function(gene_name) {
  file_name <- paste(gene_name, ".aligned.trim06.faa", sep = "")
  alignment <- phangorn::read.phyDat(file_name, format = "fasta", type = "AA")
  dist <- phangorn::dist.ml(alignment, model = "JTT")
  dist_list <- dist[1:length(dist)]
  
  data.frame(Gene = rep(gene_name, length(dist_list)), Dist = dist_list)
}

# define the dataframe that will be used as input to ggplot
dist.all.distr <- data.frame()
genes <- c("Dmrt2", "Dmrt3", "Dmrt45", "Dmrt1L")

# calculate pairwise amino acid distances of each Dmrt gene and put the results in the dataframe
for (gene in genes) {
  dist_data <- calculate_distance(gene)
  dist.all.distr <- rbind(dist.all.distr, dist_data)
}

# transform Gene column into factors
dist.all.distr$Gene <- factor(dist.all.distr$Gene, levels = genes)

### STATISTICS ABOUT DIFFERENCES IN AMINO ACID DISTANCES BETWEEN DMRT GROUPS ###

# linear model fo pairwise amino acid distance
dist.all.distr.lm <- lm(Dist ~ Gene, data = dist.all.distr)

# Q-Q plot 
qqnorm(resid(dist.all.distr.lm))

# shapiro test for normality of data
shapiro.test(resid(dist.all.distr.lm)) #not normally-distributed: W = 0.88544, p-value < 2.2e-16

# non-parametric KW test for significant differences in data
kruskal.test(dist.all.distr$Dist, dist.all.distr$Gene) #there are significant differences: KW chi-squared = 234.9, df = 3, p-value < 2.2e-16

# pairwise non-parametric Wilcoxon test
pairwise.wilcox.test(dist.all.distr$Dist, dist.all.distr$Gene, p.adjust.method = "bonferroni")

#         Dmrt2    Dmrt3    Dmrt45 
# Dmrt3   6.3e-06  -        -      
# Dmrt45  0.0058   1.8e-08  -      
# Dmrt1L  < 2e-16  < 2e-16  < 2e-16

# generate boxplot with p-values of wilcoxon test for Dmrt1L comparisons
boxplot <- ggplot(dist.all.distr,
                  aes(Gene, Dist, fill = as.factor(Gene), color = as.factor(Gene))) +
  
  # plot boxplots and scale their colors
  geom_boxplot(lwd = 1) +
  scale_fill_manual(values = alpha(c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff"), 0.5)) +
  scale_color_manual(values = c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff")) +
  
  # plot comparisons
  ylim(0, 3.75) +
  geom_segment(aes(x = "Dmrt1L", xend = "Dmrt2", y = 3.3, yend = 3.3), lwd = 0.7, col = "#7f878f") +
  geom_segment(aes(x = "Dmrt1L", xend = "Dmrt3", y = 3.5, yend = 3.5), lwd = 0.7, col = "#7f878f") +
  geom_segment(aes(x = "Dmrt1L", xend = "Dmrt45", y = 3.7, yend = 3.7), lwd = 0.7, col = "#7f878f") +
  annotate("text", label = "****", x = 2.5, y = 3.35, col = "#7f878f", size = 4) +
  annotate("text", label = "****", x = 3, y = 3.55, col = "#7f878f", size = 4) +
  annotate("text", label = "****", x = 3.5, y = 3.75, col = "#7f878f", size = 4) +

  ylab("Amino acid distance") +
  xlab("") +
  theme_light() +
  theme(axis.title.y = element_text(size = 13, vjust = 2.5, color = "#4d4d4d"),
        axis.text.x = element_text(size = 13, vjust = -0.5, face = 'bold.italic', colour = c("#bbc5d0","#bbc5d0","#bbc5d0","#0092ff")),
        legend.position = "none")

boxplot

# generate panel of compilation + boxplot
panel <- ggpubr::ggarrange(plot.compilation, boxplot,
                           ncol = 1, nrow = 2, heights = c(0.4,1), align = "v",
                           labels = c("B","C"))

panel

ggsave("plot.panel.pdf", plot = panel, device = "pdf",
	        dpi = 300, height = 6, width = 5.5, units = ("in"))
