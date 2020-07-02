#### beta div for the parasite data  ####

library(ggplot2)

uwunif <- read.delim('../../data/qiime2_190820-nochloro/beta_diversity/unweighted_unifrac_distance_matrix.txt',
                     header = T, row.names = 1, check.names = F, sep = '\t')
wtunif <- read.delim('../../data/qiime2_190820-nochloro/beta_diversity/weighted_unifrac_distance_matrix.txt',
                     header = T, row.names = 1, check.names = F, sep = '\t')
map.beta <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_190820.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

beta.samps <- intersect(rownames(uwunif), rownames(map.beta))
beta.samps <- sort(beta.samps)
uwunif <- uwunif[beta.samps, beta.samps]
wtunif <- wtunif[beta.samps, beta.samps]
map.beta <- map.beta[beta.samps, ]

unw.pca <- as.data.frame(cmdscale(uwunif, k=5))
pcnames = c()
for(i in 1:ncol(unw.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(unw.pca) = pcnames
unw.pca.meta <- merge(unw.pca, map.beta, by=0)
# The PCA plot
ggplot(unw.pca.meta, aes(x=PC1, y=PC2, color=factor(CaptiveWild))) +
  geom_point(size=4, alpha=0.6, aes(shape=SpeciesCommonName)) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC1 (49%)') + ylab('PC2 (5.9%)') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title=""), shape=guide_legend(title="Primate Species")) +
  scale_color_manual(values=c("red","blue"))
  # stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/unweighted_unifrac_captivity.png', dpi=300, height=4, width=5.5)





