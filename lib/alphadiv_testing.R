library(tidyverse)
library(ggplot2)
library(polycor)
library(beeswarm)

# Correlations

L6_otus <- read.delim('../../data/qiime2_old/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
map <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_180720.txt', header=1, row.names = 1, sep = '\t', check.names = F)

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
otus <- otus[, sample.ids]
L6_otus <- data.frame(t(L6_otus))[sample.ids, ]
map <- map[sample.ids,]

alphadiv <- read.delim('../../data/qiime2_180309/alpha_diversity/alphadiv.txt', header=1, row.names = 1, check.names = F, sep = '\t')
alphadiv <- alphadiv[sample.ids, ]

alphamap <- merge(map, alphadiv, by=0)

boxplot(alphamap$faith_pd ~ alphamap$SpeciesCommonName)

table(alphamap$Trichostrongylus_semiquant, alphamap$Trichuris_semiquant)
table(alphamap$Trichuris_semiquant, alphamap$SpeciesCommonName)

pairwise.wilcox.test(alphamap$faith_pd, alphamap$Trichuris_semiquant, paired = F, p.adjust.method = 'fdr')

cor.test(alphamap$observed_otus, alphamap$Trichostrongylus)

alphamap.howler <- alphamap[alphamap$SpeciesCommonName == "mantled howling monkey",]

boxplot(alphamap.howler$faith_pd ~ alphamap.howler$Trichostrongylus_semiquant)
boxplot(alphamap.howler$faith_pd ~ alphamap.howler$Strongyloides_semiquant)
boxplot(alphamap.howler$faith_pd ~ alphamap.howler$Trichuris_semiquant)
boxplot(alphamap.howler$faith_pd ~ alphamap.howler$all_parasites_semiquant)

