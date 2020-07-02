#### beta div for the DOUC parasite data  ####

library(ggplot2)

uwunif.d <- read.delim('../../data/qiime2_doucs_190820/beta_diversity/unweighted_unifrac_distance_matrix.txt',
                     header = T, row.names = 1, check.names = F, sep = '\t')
wtunif.d <- read.delim('../../data/qiime2_doucs_190820/beta_diversity/weighted_unifrac_distance_matrix.txt',
                     header = T, row.names = 1, check.names = F, sep = '\t')
map.beta <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_190820.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

beta.samps.d <- intersect(rownames(uwunif.d), rownames(map.beta))
beta.samps.d <- sort(beta.samps.d)
uwunif.d <- uwunif[beta.samps.d, beta.samps.d]
wtunif.d <- wtunif[beta.samps.d, beta.samps.d]
map.beta.d <- map.beta[beta.samps.d, ]

unw.d.pca <- as.data.frame(cmdscale(uwunif.d, k=5))
pcnames = c()
for(i in 1:ncol(unw.d.pca)){
  pcnames[i] <- paste0("PC",i)
}
colnames(unw.d.pca) = pcnames
unw.d.pca.meta <- merge(unw.d.pca, map.beta.d, by=0)
# The PCA plot
unw.d.pca.meta$all_parasites_semiquant <- factor(unw.d.pca.meta$all_parasites_semiquant,
                                                 levels = c("None","Low","Moderate"),
                                                 ordered = T)
ggplot(unw.d.pca.meta, aes(x=PC1, y=PC2, color=factor(all_parasites_semiquant))) +
  geom_point(size=4, alpha=0.6) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  xlab('PC1 (10%)') + ylab('PC2 (6.1%)') + 
  # xlim(-32, 32) +
  guides(color=guide_legend(title="Douc\nSemi-Quantitative\nParasite Load")) +
  scale_color_manual(values=c("blue","purple","red")) +
  stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/douc_parasitesemiquant_unwunif_beta.png', dpi=300, height=4, width=5.5)


# Now just the binary
unw.d.pca.meta$all_parasites_presence <- as.character(unw.d.pca.meta$all_parasites_presence)

ggplot(unw.d.pca.meta, aes(x=PC1, y=PC2, group=all_parasites_presence)) +
  geom_point(size=4, alpha=0.8, aes(shape=all_parasites_presence)) +
  theme_classic() + theme(aspect.ratio = 1, axis.text = element_text(color='black', size=10)) +
  scale_shape_manual(values = c('0'=19,'1'=1)) +
  xlab('PC1 (10%)') + ylab('PC2 (6.1%)') + 
  guides(shape=guide_legend(title="Douc\nParasite Presence"))
  # stat_ellipse(level = 0.95, linetype=2)

ggsave(filename = '../../results/douc_parasites-presence_unwunif_beta.png', dpi=300, height=4, width=5.5)

vegan::adonis2(dist(uwunif.d)~map.beta.d$all_parasites_presence, permutations = 99)
# Not significant


