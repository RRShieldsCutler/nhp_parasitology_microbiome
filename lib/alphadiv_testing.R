library(tidyverse)
library(ggplot2)
library(polycor)
library(beeswarm)

# Correlations

L6_otus <- read.delim('../../data/qiime2_old/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
map <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_190612.txt', header=1, row.names = 1, sep = '\t', check.names = F)

#CLR
otu.t = t(L6_otus); eps = 0.5
otu.t = otu.t*(1 - rowSums(otu.t==0)*eps/rowSums(otu.t))
otu.t[otu.t==0]=eps
otu.t = sweep(otu.t,1,rowSums(otu.t),'/');
ls = log(otu.t)
otu.t = t(ls - rowMeans(ls))
L6_otus = otu.t[,!is.nan(colSums(otu.t))]

sample.ids <- intersect(colnames(L6_otus), rownames(map))
sample.ids <- sort(sample.ids)
# otus <- otus[, sample.ids]
L6_otus <- data.frame(t(L6_otus))[sample.ids, ]
map <- map[sample.ids,]

alphadiv <- read.delim('../../data/qiime2_190705/alpha_diversity/alphadiv.txt', header=1, row.names = 1, check.names = F, sep = '\t')
alphadiv <- alphadiv[sample.ids, ]

alphamap <- merge(map, alphadiv, by=0)
# alphamap$Trichostrongylus_semiquant <- factor(alphamap$Trichostrongylus_semiquant, levels = c("None","Low","Moderate"), ordered = T)
# alphamap$all_parasites_semiquant <- factor(alphamap$all_parasites_semiquant, levels = c("None","Low","Moderate"), ordered = T)

boxplot(alphamap$shannon ~ alphamap$SpeciesCommonName)

table(alphamap$Trichostrongylus_semiquant, alphamap$Trichuris_semiquant)
table(alphamap$Trichuris_semiquant, alphamap$SpeciesCommonName)

cor.test(alphamap$observed_otus, alphamap$Trichostrongylus)

boxplot(alphamap$faith_pd ~ alphamap$all_parasites_semiquant)

boxplot(alphamap$observed_otus ~ alphamap$all_parasites_semiquant)
boxplot(alphamap$pielou_e ~ alphamap$all_parasites_semiquant)
boxplot(alphamap$shannon ~ alphamap$all_parasites_semiquant)

# Significantly higher diversity in samples with higher parasite loads
pairwise.wilcox.test(alphamap$faith_pd, alphamap$all_parasites_semiquant, paired = F, p.adjust.method = 'fdr')

alphamap.howler <- alphamap[alphamap$SpeciesCommonName == "mantled howling monkey",]
pairwise.wilcox.test(alphamap.howler$faith_pd, alphamap.howler$Strongyloides_semiquant, paired = F, p.adjust.method = 'fdr')

boxplot(alphamap.howler$observed_otus ~ alphamap.howler$all_parasites_presence)
boxplot(alphamap.howler$faith_pd ~ alphamap.howler$Strongyloides_semiquant)
boxplot(alphamap.howler$faith_pd ~ alphamap.howler$Trichuris_semiquant)


alphamap.douc <- alphamap[alphamap$SpeciesCommonName == "red-shanked douc",]
pairwise.wilcox.test(alphamap.douc$faith_pd, alphamap.douc$Trichuris_semiquant, paired = F, p.adjust.method = 'fdr')

boxplot(alphamap.douc$observed_otus ~ alphamap.douc$all_parasites_presence)
boxplot(alphamap.douc$faith_pd ~ alphamap.douc$Strongyloides_semiquant)
boxplot(alphamap.douc$faith_pd ~ alphamap.douc$Trichuris_semiquant)

pairwise.wilcox.test(alphamap.douc$faith_pd, alphamap.douc$distinct_parasite_taxa, paired = F, p.adjust.method = 'fdr')


