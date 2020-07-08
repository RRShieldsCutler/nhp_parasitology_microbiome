## OTU table processing
# QC by shi7
# shi7.py -i raw_fastq/ -o shi7_stitch -t 48 --adaptor Nextera --trim_qual 34 --filter_qual 37 -s .,1 --min_overlap 233 --max_overlap 281 --allow_outies False
# Alignment with BURST12
# burst12 -r /project/flatiron2/sop/gg97.edx -a /project/flatiron2/sop/gg97.acx -b /project/flatiron2/sop/gg97.tax -q shi7_stitch/combined_seqs.fna -o douc_howler_stitch_gg_95.b6 -n -m CAPITALIST -bs -i 0.95 -fr -sa -t 48


library(ggplot2)
library(vegan)
library(reshape2)
library(gplots)
library(viridis)

otus <- read.delim('../../data/douc_howler_gg95id_otu.txt', header=1, row.names = 1, check.names = F, sep = '\t')
colnames(otus) <- unlist(lapply(colnames(otus), function(xx) strsplit(xx, '.', fixed = T)[[1]][1]))

sum(otus)  # 5.3 million reads, so we'll drop anything with less than 5 counts

otus <- otus[(rowSums(otus) >= 5), ]  # Remove rare

otus = otus[rowMeans(otus > 0) >= 0.05, ]  # Remove rare OTUs in less than 5% of samples
depths = colSums(otus)
hist(depths, breaks=40)
summary(depths)
### Remove samples with low depth?
sort(depths)[1:50]
# All look ok!
# Rarefaction level could be 21990
dim(otus)

otu.counts = colSums(otus > 0)
hist(otu.counts, breaks=30)

top_otu <- rowSums(otus)
tail(sort(top_otu))

# write the filtered table
sink('../../data/filtered_nhp_gg_otu.txt'); cat('#OTU ID\t')
write.table(otus, file = '../../data/filtered_nhp_gg_otu.txt', sep = '\t', row.names = T, quote = F, append = T)
sink(NULL)

# CLR abundance OTU table
otus.c <- t(otus); eps <- 0.5
otus.c <- otus.c * (1 - rowSums(otus.c==0) * eps / rowSums(otus.c))
otus.c[otus.c == 0] <- eps
otus.c <- sweep(otus.c, 1, rowSums(otus.c), '/')
ls <- log(otus.c)
otus.c <- t(ls - rowMeans(ls))
otus.c <- otus.c[, !is.nan(colSums(otus.c))]
otus.clr <- as.data.frame(otus.c)

# Write the CLR table
sink('../../data/filtered_CLR_nhp_gg_otu.txt'); cat('#OTU ID\t')
write.table(otus, file = '../../data/filtered_CLR_nhp_gg_otu.txt', sep = '\t', row.names = T, quote = F, append = T)
sink(NULL)

L6_otus <- read.delim('../../data/qiime2_190705/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
L3_otus <- read.delim('../../data/qiime2_190705/otu_table_L3.txt', header = 1, sep = '\t', row.names = 1, skip=1)
map <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_190612.txt', header=1, row.names = 1, sep = '\t', check.names = F)

sample.ids <- intersect(colnames(otus), rownames(map))
sample.ids <- sort(sample.ids)
otus <- otus[, sample.ids]
L6_otus <- data.frame(t(L6_otus))[sample.ids, ]
map <- map[sample.ids,]

trich.cor.result <- psych::corr.test(L6_otus, y = map[, c(24,25,26)], adjust = 'fdr')
sum(trich.cor.result$p < 0.05)
colnames(trich.cor.result$p)
sort(trich.cor.result$p[,1])[1:20]
sort(trich.cor.result$p[,1])[1:20]
sort(trich.cor.result$p[,1])[1:20]



##### Heatmap for L3 #######

L3_otus <- L3_otus[,sample.ids]

otus.c <- t(L3_otus); eps <- 0.5
otus.c <- otus.c * (1 - rowSums(otus.c==0) * eps / rowSums(otus.c))
otus.c[otus.c == 0] <- eps
otus.c <- sweep(otus.c, 1, rowSums(otus.c), '/')
ls <- log(otus.c)
otus.c <- t(ls - rowMeans(ls))
otus.c <- otus.c[, !is.nan(colSums(otus.c))]
L3_otus.clr <- as.data.frame(otus.c)

map2 <- map[order(map$all_parasites_ct), ]
ctorder <- as.character(rownames(map2))

L3_otus.clr <- L3_otus.clr[, ctorder]

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")

heatmap.2(as.matrix(L3_otus.clr), Colv = FALSE, symm = F, symbreaks = T, scale="none", 
          symkey = F, dendrogram="row",trace="none", margin=c(5,20),
          hclust=hclustfunc,distfun=distfunc, col=viridis(n = 300),
          breaks=seq(-3.44, 9.48, length.out = 301))

ctsumm <- summary(map2$all_parasites_ct)
ctlength <- ctsumm[[6]][1] - ctsumm[[1]][1] + 1
rowpal <- colorRampPalette(c("yellow","red"))(n=ctlength)
names(rowpal) <- as.character(seq(0,143,by=1))
rowpal.in <- rowpal[as.character(map2$all_parasites_ct)]

lowday <- min(map$relative_day)
maxday <- max(map$relative_day)
ctlength <- maxday - lowday + 1
rowpal <- colorRampPalette(c("yellow","red"))(n=ctlength)
names(rowpal) <- as.character(seq(lowday, maxday, by=1))
rowpal.in <- rowpal[as.character(map$relative_day)]



heatmap.2(as.matrix(L3_otus.clr), Colv = FALSE, symm = F, symbreaks = T, scale="none", 
          symkey = F, dendrogram="row",trace="none", margin=c(5,20),
          hclust=hclustfunc,distfun=distfunc, col=viridis(n = 300),
          breaks=seq(-3.44, 9.48, length.out = 301), ColSideColors = rowpal.in)
nrow(as.matrix(L3_otus.clr))
length(rowpal.in)


###### Summer 2020 Updates #########
# Let's use boruta to check for taxa that might track with infection.
# Howlers
otus.L6.h <- read.delim('../../data/qiime_200602_wild-howlers/otu_table_L6.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.h.rf <- read.delim('../../data/howler_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.h <- intersect(colnames(otus.L6.h), rownames(map.h.rf))
samps.h <- sort(samps.h)
otus.L6.h <- otus.L6.h[,samps.h]
map.h.rf <- map.h.rf[samps.h, ]

otus.L6.h.c <- t(otus.L6.h); eps <- 0.5
otus.L6.h.c <- otus.L6.h.c * (1 - rowSums(otus.L6.h.c==0) * eps / rowSums(otus.L6.h.c))
otus.L6.h.c[otus.L6.h.c == 0] <- eps
otus.L6.h.c <- sweep(otus.L6.h.c, 1, rowSums(otus.L6.h.c), '/')
ls <- log(otus.L6.h.c)
otus.L6.h.c <- t(ls - rowMeans(ls))
otus.L6.h.c <- otus.L6.h.c[, !is.nan(colSums(otus.L6.h.c))]
otus.L6.h.clr <- as.data.frame(otus.L6.h.c)

#### Boruta feature selection by RF for L6####
library(Boruta)

L6.h.clr.t <- t(otus.L6.h.clr)[samps.h,]
map.h.rf <- map.h.rf[samps.h, ]
boruta.rfimp <- Boruta(x=L6.h.clr.t, y=factor(map.h.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
length(sig.taxa)


# Try it with L4
otus.L4.h <- read.delim('../../data/qiime_200602_wild-howlers/otu_table_L4.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.h.rf <- read.delim('../../data/howler_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.h <- intersect(colnames(otus.L4.h), rownames(map.h.rf))
samps.h <- sort(samps.h)
otus.L4.h <- otus.L4.h[,samps.h]
map.h.rf <- map.h.rf[samps.h, ]

otus.L4.h.c <- t(otus.L4.h); eps <- 0.5
otus.L4.h.c <- otus.L4.h.c * (1 - rowSums(otus.L4.h.c==0) * eps / rowSums(otus.L4.h.c))
otus.L4.h.c[otus.L4.h.c == 0] <- eps
otus.L4.h.c <- sweep(otus.L4.h.c, 1, rowSums(otus.L4.h.c), '/')
ls <- log(otus.L4.h.c)
otus.L4.h.c <- t(ls - rowMeans(ls))
otus.L4.h.c <- otus.L4.h.c[, !is.nan(colSums(otus.L4.h.c))]
otus.L4.h.clr <- as.data.frame(otus.L4.h.c)

#### Boruta feature selection by RF for L4 ####

L4.h.clr.t <- t(otus.L4.h.clr)[samps.h,]
map.h.rf <- map.h.rf[samps.h, ]
boruta.rfimp <- Boruta(x=L4.h.clr.t, y=factor(map.h.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
## None


## Try with L5

otus.L5.h <- read.delim('../../data/qiime_200602_wild-howlers/otu_table_L5.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.h.rf <- read.delim('../../data/howler_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.h <- intersect(colnames(otus.L5.h), rownames(map.h.rf))
samps.h <- sort(samps.h)
otus.L5.h <- otus.L5.h[,samps.h]
map.h.rf <- map.h.rf[samps.h, ]

otus.L5.h.c <- t(otus.L5.h); eps <- 0.5
otus.L5.h.c <- otus.L5.h.c * (1 - rowSums(otus.L5.h.c==0) * eps / rowSums(otus.L5.h.c))
otus.L5.h.c[otus.L5.h.c == 0] <- eps
otus.L5.h.c <- sweep(otus.L5.h.c, 1, rowSums(otus.L5.h.c), '/')
ls <- log(otus.L5.h.c)
otus.L5.h.c <- t(ls - rowMeans(ls))
otus.L5.h.c <- otus.L5.h.c[, !is.nan(colSums(otus.L5.h.c))]
otus.L5.h.clr <- as.data.frame(otus.L5.h.c)

#### Boruta feature selection by RF for L5 ####

L5.h.clr.t <- t(otus.L5.h.clr)[samps.h,]
map.h.rf <- map.h.rf[samps.h, ]
boruta.rfimp <- Boruta(x=L5.h.clr.t, y=factor(map.h.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
## None


## Doucs
otus.L6.d <- read.delim('../../data/qiime_200602_doucs/otu_table_L6.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.d.rf <- read.delim('../../data/douc_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.d <- intersect(colnames(otus.L6.d), rownames(map.d.rf))
samps.d <- sort(samps.d)
otus.L6.d <- otus.L6.d[,samps.d]
map.d.rf <- map.d.rf[samps.d, ]

otus.L6.d.c <- t(otus.L6.d); eps <- 0.5
otus.L6.d.c <- otus.L6.d.c * (1 - rowSums(otus.L6.d.c==0) * eps / rowSums(otus.L6.d.c))
otus.L6.d.c[otus.L6.d.c == 0] <- eps
otus.L6.d.c <- sweep(otus.L6.d.c, 1, rowSums(otus.L6.d.c), '/')
ls <- log(otus.L6.d.c)
otus.L6.d.c <- t(ls - rowMeans(ls))
otus.L6.d.c <- otus.L6.d.c[, !is.nan(colSums(otus.L6.d.c))]
otus.L6.d.clr <- as.data.frame(otus.L6.d.c)

#### Boruta feature selection by RF for L6####
library(Boruta)

L6.d.clr.t <- t(otus.L6.d.clr)[samps.d,]
map.d.rf <- map.d.rf[samps.d, ]
boruta.rfimp <- Boruta(x=L6.d.clr.t, y=factor(map.d.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
# NONE


# Try it with L4
otus.L4.d <- read.delim('../../data/qiime_200602_doucs/otu_table_L4.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.d.rf <- read.delim('../../data/douc_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.d <- intersect(colnames(otus.L4.d), rownames(map.d.rf))
samps.d <- sort(samps.d)
otus.L4.d <- otus.L4.d[,samps.d]
map.d.rf <- map.d.rf[samps.d, ]

otus.L4.d.c <- t(otus.L4.d); eps <- 0.5
otus.L4.d.c <- otus.L4.d.c * (1 - rowSums(otus.L4.d.c==0) * eps / rowSums(otus.L4.d.c))
otus.L4.d.c[otus.L4.d.c == 0] <- eps
otus.L4.d.c <- sweep(otus.L4.d.c, 1, rowSums(otus.L4.d.c), '/')
ls <- log(otus.L4.d.c)
otus.L4.d.c <- t(ls - rowMeans(ls))
otus.L4.d.c <- otus.L4.d.c[, !is.nan(colSums(otus.L4.d.c))]
otus.L4.d.clr <- as.data.frame(otus.L4.d.c)

#### Boruta feature selection by RF for L4 ####

L4.d.clr.t <- t(otus.L4.d.clr)[samps.d,]
map.d.rf <- map.d.rf[samps.d, ]
boruta.rfimp <- Boruta(x=L4.d.clr.t, y=factor(map.d.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))

length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative'])
tent.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
#no tent

sig.taxa # "k__Bacteria;p__Tenericutes;c__RF3;o__ML615J;28"



## Try with L5

otus.L5.d <- read.delim('../../data/qiime_200602_doucs/otu_table_L5.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.d.rf <- read.delim('../../data/douc_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.d <- intersect(colnames(otus.L5.d), rownames(map.d.rf))
samps.d <- sort(samps.d)
otus.L5.d <- otus.L5.d[,samps.d]
map.d.rf <- map.d.rf[samps.d, ]

otus.L5.d.c <- t(otus.L5.d); eps <- 0.5
otus.L5.d.c <- otus.L5.d.c * (1 - rowSums(otus.L5.d.c==0) * eps / rowSums(otus.L5.d.c))
otus.L5.d.c[otus.L5.d.c == 0] <- eps
otus.L5.d.c <- sweep(otus.L5.d.c, 1, rowSums(otus.L5.d.c), '/')
ls <- log(otus.L5.d.c)
otus.L5.d.c <- t(ls - rowMeans(ls))
otus.L5.d.c <- otus.L5.d.c[, !is.nan(colSums(otus.L5.d.c))]
otus.L5.d.clr <- as.data.frame(otus.L5.d.c)

#### Boruta feature selection by RF for L5 ####

L5.d.clr.t <- t(otus.L5.d.clr)[samps.d,]
map.d.rf <- map.d.rf[samps.d, ]
boruta.rfimp <- Boruta(x=L5.d.clr.t, y=factor(map.d.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed'])
# sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
tent.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                    FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))



tent.taxa
# "k__Bacteria;p__Tenericutes;c__RF3;o__ML615J;28;f__"

L5.d.clr.t.df <- data.frame(L5.d.clr.t)

boxplot(L5.d.clr.t.df$k__Bacteria.p__Tenericutes.c__RF3.o__ML615J.28.f__ ~ map.d.rf$all_parasites_presence)
wilcox.test(L5.d.clr.t.df$k__Bacteria.p__Tenericutes.c__RF3.o__ML615J.28.f__ ~ map.d.rf$all_parasites_presence, paired=F)
# p = 0.006
# Should go to L3 because of annotations


# L3
otus.L3.d <- read.delim('../../data/qiime_200602_doucs/otu_table_L3.txt', header=1,
                        row.names = 1, check.names = F, skip=1, sep = '\t')
map.d.rf <- read.delim('../../data/douc_animal_parasite_map_200602.txt',
                       header = T, row.names = 1, check.names = F, sep = '\t')

samps.d <- intersect(colnames(otus.L3.d), rownames(map.d.rf))
samps.d <- sort(samps.d)
otus.L3.d <- otus.L3.d[,samps.d]
map.d.rf <- map.d.rf[samps.d, ]

otus.L3.d.c <- t(otus.L3.d); eps <- 0.5
otus.L3.d.c <- otus.L3.d.c * (1 - rowSums(otus.L3.d.c==0) * eps / rowSums(otus.L3.d.c))
otus.L3.d.c[otus.L3.d.c == 0] <- eps
otus.L3.d.c <- sweep(otus.L3.d.c, 1, rowSums(otus.L3.d.c), '/')
ls <- log(otus.L3.d.c)
otus.L3.d.c <- t(ls - rowMeans(ls))
otus.L3.d.c <- otus.L3.d.c[, !is.nan(colSums(otus.L3.d.c))]
otus.L3.d.clr <- as.data.frame(otus.L3.d.c)

#### Boruta feature selection by RF for L3 ####

L3.d.clr.t <- t(otus.L3.d.clr)[samps.d,]
map.d.rf <- map.d.rf[samps.d, ]
boruta.rfimp <- Boruta(x=L3.d.clr.t, y=factor(map.d.rf$all_parasites_presence))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed'])
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
tent.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Tentative']),
                    FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
sig.taxa
# "k__Bacteria;p__Tenericutes;c__RF3"

L3.d.clr.t.df <- data.frame(L3.d.clr.t)
boxplot(L3.d.clr.t.df$k__Bacteria.p__Tenericutes.c__RF3 ~ map.d.rf$all_parasites_presence)
wilcox.test(L3.d.clr.t.df$k__Bacteria.p__Tenericutes.c__RF3 ~ map.d.rf$all_parasites_presence, paired=F)
p = 0.000984

L3.metamerge <- merge(L3.d.clr.t.df, map.d.rf, by=0)

ggplot(L3.metamerge, aes(x=all_parasites_presence,
                         y=k__Bacteria.p__Tenericutes.c__RF3,
                         group=all_parasites_presence)) +
  xlab("") + ylab("Tenericutes clade, class RF3\n(CLR abundance)") +
  geom_violin(aes(fill=all_parasites_presence)) +
  geom_jitter(height = 0, width = 0.15, size=0.6, alpha=0.7) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  guides(fill=F) +
  annotate("text", x=1.5, y=2.1, label="italic(p) == 0.001", size=3.5, parse=T)
ggsave('../../results/douc_tenericutes_violin.png', height = 4, width = 4, dpi=300)  
