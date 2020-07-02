#### beta div for the DOUC parasite data  ####

library(ggplot2)

uwunif.h <- read.delim('../../data/qiime2_howlers_190820//beta_diversity/unweighted_unifrac_distance_matrix.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')
wtunif.h <- read.delim('../../data/qiime2_howlers_190820//beta_diversity/weighted_unifrac_distance_matrix.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')
map.beta <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_190820.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

beta.samps.h <- intersect(rownames(uwunif.h), rownames(map.beta))
beta.samps.h <- sort(beta.samps.h)
uwunif.h <- uwunif[beta.samps.h, beta.samps.h]
wtunif.h <- wtunif[beta.samps.h, beta.samps.h]
map.beta.h <- map.beta[beta.samps.h, ]

unw.h.pca <- as.data.frame(cmdscale(uwunif.h, k=5))
pcnames = c()
for(i in 1:ncol(unw.h.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(unw.h.pca) = pcnames
unw.h.pca.meta <- merge(unw.h.pca, map.beta.h, by=0)
# The PCA plot
unw.h.pca.meta$all_parasites_semiquant <- factor(unw.h.pca.meta$all_parasites_semiquant,
                                                 levels = c("None","Low","Moderate"),
                                                 ordered = T)
ggplot(unw.h.pca.meta, aes(x=PC1, y=PC2, color=factor(all_parasites_semiquant))) +
  geom_point(size=4, alpha=0.6) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC1 (20%)') + ylab('PC2 (4.9%)') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title="All Howler\nSemi-Quantitative\nParasite Load")) +
  scale_color_manual(values=c("blue","purple","red"))
  # stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/howler-all_parasitesemiquant_unwunif_beta.png', dpi=300, height=4, width=5.5)

ggplot(unw.h.pca.meta[unw.h.pca.meta$CaptiveWild=="Wild",], aes(x=PC1, y=PC2, color=factor(all_parasites_semiquant))) +
  geom_point(size=4, alpha=0.6) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC1') + ylab('PC2') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title="Wild Howler\nSemi-Quantitative\nParasite Load")) +
  scale_color_manual(values=c("blue","purple","red")) +
  stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/howler-wild_parasitesemiquant_unwunif_beta.png', dpi=300, height=4, width=5.5)

