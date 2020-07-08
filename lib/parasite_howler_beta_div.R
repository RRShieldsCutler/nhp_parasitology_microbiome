#### beta div for the DOUC parasite data  ####

library(ggplot2)

# WILD ONLY

uwunif.h <- read.delim('../../data/qiime_200602_wild-howlers/beta_diversity/unweighted_unifrac_distance_matrix.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')
wtunif.h <- read.delim('../../data/qiime_200602_wild-howlers/beta_diversity/weighted_unifrac_distance_matrix.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')
map.beta <- read.delim('../../data/howler_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

beta.samps.h <- intersect(rownames(uwunif.h), rownames(map.beta))
beta.samps.h <- sort(beta.samps.h)
uwunif.h <- uwunif[beta.samps.h, beta.samps.h]
wtunif.h <- wtunif[beta.samps.h, beta.samps.h]
map.beta.h <- map.beta[beta.samps.h, ]

unw.h.pca <- as.data.frame(cmdscale(uwunif.h, k=5))
unw.h.stat <- cmdscale(uwunif.h, k=5, eig = T, add = T)
round(unw.h.stat$eig*100/sum(unw.h.stat$eig),1)

pcnames = c()
for(i in 1:ncol(unw.h.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(unw.h.pca) = pcnames
unw.h.pca.meta <- merge(unw.h.pca, map.beta.h, by=0)
# The PCA plot
unw.h.pca.meta$all_parasites_semiquant <- factor(unw.h.pca.meta$all_parasites_semiquant,
                                                 levels = c("none","low","high"),
                                                 ordered = T)
ggplot(unw.h.pca.meta, aes(x=PC2, y=PC3, color=factor(all_parasites_semiquant))) +
  geom_point(size=4, alpha=0.6) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC2 (6.7%)') + ylab('PC3 (5.7%)') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title="Wild Howler\nSemi-Quantitative\nParasite Load")) +
  scale_color_manual(values=c("blue","purple","red"))
  # stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/howler-wild_parasitesemiquant_unwunif_beta.png', dpi=300, height=4, width=5.5)



wt.h.pca <- as.data.frame(cmdscale(wtunif.h, k=5))
wt.h.stat <- cmdscale(wtunif.h, k=5, eig = T, add = T)
round(wt.h.stat$eig*100/sum(wt.h.stat$eig),1)

pcnames = c()
for(i in 1:ncol(wt.h.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(wt.h.pca) = pcnames
wt.h.pca.meta <- merge(wt.h.pca, map.beta.h, by=0)
# The PCA plot
wt.h.pca.meta$all_parasites_semiquant <- factor(wt.h.pca.meta$all_parasites_semiquant,
                                                 levels = c("none","low","high"),
                                                 ordered = T)

ggplot(wt.h.pca.meta, aes(x=PC1, y=PC2, color=factor(all_parasites_semiquant))) +
  geom_point(size=4, alpha=0.6) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC1 (18.9%)') + ylab('PC2 (11.8%)') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title="Wild Howler\nSemi-Quantitative\nParasite Load")) +
  scale_color_manual(values=c("blue","purple","red"))
  # stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/howler-wild_parasitesemiquant_wtunif_beta.png', dpi=300, height=4, width=5.5)

