library(tidyverse)
library(ggplot2)
library(polycor)
library(beeswarm)

# Correlations

L6_otus <- read.delim('../../data/qiime2_190705/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
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

trich.cor.result <- psych::corr.test(L6_otus, y = map[, c(24,25,26)], adjust = 'fdr')
sum(trich.cor.result$p < 0.05)
colnames(trich.cor.result$p)
sort(trich.cor.result$p[,1])[1:20]
sort(trich.cor.result$p[,1])[1:20]
sort(trich.cor.result$p[,1])[1:20]

# Or, for heatmap

# set x and y to be your dataframes of interest, rows as subjects and columns as variables of interest
y <- L6_otus
sum(is.na(y))
dim(y)
x <- map[, c(24,25,26)]
sum(is.na(x))
dim(x)

# create the dataframe that will hold your correlations
cors <- NULL

for (i in colnames(x)) {
  for (j in colnames(y)) {
    a <- as.numeric(x[,i])
    b <- as.numeric(y[,j])
    tmp <- c(i, j, cor(a, b, method = "pearson"), cor.test(a,b, method = "pearson", exact = F)$p.value)   # you can tweak this line to suppress warnings if you like
    cors <- rbind(cors,tmp)
  }
}

cors_stash_b_genus <- cors
# save(cors_stash_b_genus, file = '../data/metabolomics/nmr_cors_genus_stool-post.rdata')
cors<- data.frame(row.names = NULL, cors, stringsAsFactors = FALSE)

colnames(cors) <- c("parasite", "organism", "correlation", "pvalue")  
cors$correlation <- as.numeric(cors$correlation)
cors$pvalue <- as.numeric(cors$pvalue)
cors$qvalue <- p.adjust(cors$pvalue, method = "fdr")
cors$Significance<-cut(cors$qvalue, breaks=c(-Inf, 0.05, 0.25, Inf), label=c("**","*", ""))


# prune for complete cases with significant correlations
dat.c<-cors[complete.cases(cors),]
xkeep<-unique(dat.c$parasite[dat.c$qvalue<=0.25])  # adjust this q-value for your liking, 1 will show everything, <0.25, only signif.
xassign<-xkeep
dat.w<-dat.c[dat.c$parasite %in% xassign,] 
ykeep <- unique(dat.w$organism[dat.w$qvalue<=0.25]) # adjust this one too as above
dat.m<-dat.w[dat.w$organism %in% ykeep,]


#reorder variables for plotting if more than one observation
order_y <- dat.m %>% select(parasite, organism, correlation) %>% na.omit()
order_y <- order_y %>% spread(parasite, correlation)
order_y <- remove_rownames(order_y)
order_y <- column_to_rownames(order_y, "organism")
yorder <- hclust((dist(1-cor(order_y))/2))$order


# reorder the other variables for plotting
order_x <- dat.m %>% select(organism, parasite, correlation) %>% na.omit()
order_x <- order_x %>% spread(organism, correlation)
order_x <- remove_rownames(order_x)
order_x<- column_to_rownames(order_x, "parasite")
xorder <- hclust((dist(1-cor(order_x))/2))$order

dat.m$parasite <- as.factor(dat.m$parasite)
dat.m$parasite_ordered <- factor(dat.m$parasite, levels(dat.m$parasite)[yorder])

# dat.m$organism <- as.factor(dat.m$organism)
# dat.m$organism <- factor(dat.m$organism, levels(dat.m$organism)[xorder])

dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism), FUN = function(xx) gsub(x = xx, pattern = '.p__', replacement = ';p__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.c__', replacement = ';c__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.o__', replacement = ';o__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.f__', replacement = ';f__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.g__', replacement = ';', fixed = T)))

dat.m$organism_short <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) tail(strsplit(x=xx, split = as.character(';'), fixed = T)[[1]], n=1)))

dat.m$organism_short <- as.factor(dat.m$organism_short)
dat.m$organism_short_ordered <- factor(dat.m$organism_short, levels(dat.m$organism_short)[xorder])

summary(dat.m$correlation)
summary(dat.m$qvalue)
table(dat.m$Significance)
#plot the correlations
myheatmap <- ggplot(data = dat.m, aes(x=parasite_ordered, y=organism_short_ordered, fill=correlation)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low="yellow", mid="black", high="blue", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Pearson") +
  geom_text(data = dat.m, aes(label=Significance), color="white", size=4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color='black'), 
        axis.text.y = element_text(size = 8, color='black'),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.ticks = element_blank(),
        legend.key.size = unit(0.1, "inch"),
        plot.title = element_text(size=9),
        panel.background = element_rect(fill='white')) +
  labs(x = "Parasite", y = "Genus or higher (collapsed to L6)", title='Genus-level correlations\nwith three parasites\n(microscopy egg counts)') +
  # theme(strip.text.y = element_text(angle = 0, face = "italic"), 
  #       strip.background = element_rect(color="white", fill = "white")) +
  coord_fixed(ratio=0.4)

myheatmap
ggsave(myheatmap, filename = '../../results/microbiome_parasite_correlations.png', dpi = 300, height = 10, width = 4.5)


# Just do within the Doucs or howlers
douc.ids <- rownames(map[map$SpeciesCommonName == 'red-shanked douc', ])
L6.douc <- L6_otus[douc.ids,]
map.douc <- map[douc.ids, ]

howler.ids <- rownames(map[map$SpeciesCommonName == 'mantled howling monkey', ])
L6.howler <- L6_otus[howler.ids,]
map.howler <- map[howler.ids, ]


L3_otus <- read.delim('../../data/qiime2_180309/otu_table_L3.txt', header = 1, sep = '\t', row.names = 1, skip=1)
L3.douc <- L3_otus[,douc.ids]
L3.howler <- L3_otus[,howler.ids]


## For the howler:
# set x and y to be your dataframes of interest, rows as subjects and columns as variables of interest
y <- data.frame(t(L6.howler))
sum(is.na(y))
dim(y)
x <- map.howler[, c(24,26)]
sum(is.na(x))
dim(x)

# create the dataframe that will hold your correlations
cors <- NULL

for (i in colnames(x)) {
  for (j in colnames(y)) {
    a <- as.numeric(x[,i])
    b <- as.numeric(y[,j])
    tmp <- c(i, j, cor(a, b, method = "pearson"), cor.test(a,b, method = "pearson", exact = F)$p.value)   # you can tweak this line to suppress warnings if you like
    cors <- rbind(cors,tmp)
  }
}

cors_stash_b_genus <- cors
# save(cors_stash_b_genus, file = '../data/metabolomics/nmr_cors_genus_stool-post.rdata')
cors<- data.frame(row.names = NULL, cors, stringsAsFactors = FALSE)

colnames(cors) <- c("parasite", "organism", "correlation", "pvalue")  
cors$correlation <- as.numeric(cors$correlation)
cors$pvalue <- as.numeric(cors$pvalue)
cors$qvalue <- p.adjust(cors$pvalue, method = "fdr")
cors$Significance<-cut(cors$qvalue, breaks=c(-Inf, 0.05, 0.25, Inf), label=c("**","*", ""))


# prune for complete cases with significant correlations
dat.c<-cors[complete.cases(cors),]
xkeep<-unique(dat.c$parasite[dat.c$qvalue<=0.5])  # adjust this q-value for your liking, 1 will show everything, <0.25, only signif.
xassign<-xkeep
dat.w<-dat.c[dat.c$parasite %in% xassign,] 
ykeep <- unique(dat.w$organism[dat.w$qvalue<=0.5]) # adjust this one too as above
dat.m<-dat.w[dat.w$organism %in% ykeep,]


#reorder variables for plotting if more than one observation
order_y <- dat.m %>% select(parasite, organism, correlation) %>% na.omit()
order_y <- order_y %>% spread(parasite, correlation)
order_y <- remove_rownames(order_y)
order_y <- column_to_rownames(order_y, "organism")
yorder <- hclust((dist(1-cor(order_y))/2))$order


# reorder the other variables for plotting
order_x <- dat.m %>% select(organism, parasite, correlation) %>% na.omit()
order_x <- order_x %>% spread(organism, correlation)
order_x <- remove_rownames(order_x)
order_x<- column_to_rownames(order_x, "parasite")
xorder <- hclust((dist(1-cor(order_x))/2))$order

dat.m$parasite <- as.factor(dat.m$parasite)
dat.m$parasite_ordered <- factor(dat.m$parasite, levels(dat.m$parasite)[yorder])

# dat.m$organism <- as.factor(dat.m$organism)
# dat.m$organism <- factor(dat.m$organism, levels(dat.m$organism)[xorder])

dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism), FUN = function(xx) gsub(x = xx, pattern = '.p__', replacement = ';p__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.c__', replacement = ';c__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.o__', replacement = ';o__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.f__', replacement = ';f__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.g__', replacement = ';', fixed = T)))

dat.m$organism_short <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) tail(strsplit(x=xx, split = as.character(';'), fixed = T)[[1]], n=1)))

dat.m$organism_short <- as.factor(dat.m$organism_short)
dat.m$organism_short_ordered <- factor(dat.m$organism_short, levels(dat.m$organism_short)[xorder])

summary(dat.m$correlation)
summary(dat.m$qvalue)
table(dat.m$Significance)
#plot the correlations
myheatmap <- ggplot(data = dat.m, aes(x=parasite_ordered, y=organism_short_ordered, fill=correlation)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low="yellow", mid="black", high="blue", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Pearson") +
  geom_text(data = dat.m, aes(label=Significance), color="white", size=4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color='black'), 
        axis.text.y = element_text(size = 8, color='black'),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.ticks = element_blank(),
        legend.key.size = unit(0.1, "inch"),
        plot.title = element_text(size=9),
        panel.background = element_rect(fill='white')) +
  labs(x = "Parasite", y = "Genus or higher (collapsed to L6)", title='HOWLER\nGenus-level correlations\nwith three parasites\n(microscopy egg counts)') +
  # theme(strip.text.y = element_text(angle = 0, face = "italic"), 
  #       strip.background = element_rect(color="white", fill = "white")) +
  coord_fixed(ratio=.75)

myheatmap
ggsave(myheatmap, filename = '../../results/howler_microbiome_parasite_correlations_L3.png', dpi = 300, height = 4.5, width = 3.2)

table(x$Trichostrongylus > 0)


## Adding some metadata updates after Tom's comment about using counts as a proxy for intensity/load

map$Strongyloides_semiquant <- cut(map$Strongyloides, right = F,
                                   breaks=c(-Inf, 1, 10, 100, 1000), label=c("None","Low", "Moderate","High"))
map$Trichuris_semiquant <- cut(map$Trichuris, right = F,
                               breaks=c(-Inf, 1, 10, 100, 1000), label=c("None","Low", "Moderate","High"))
map$Trichostrongylus_semiquant <- cut(map$Trichostrongylus, right = F,
                                      breaks=c(-Inf, 1, 10, 100, 1000), label=c("None","Low", "Moderate","Moderate"))

map$Strongyloides_semiquant <- factor(map$Strongyloides_semiquant, levels = c('None','Low','Moderate','High'), ordered = T)
map$Trichuris_semiquant <- factor(map$Trichuris_semiquant, levels = c('None','Low','Moderate','High'), ordered = T)
map$Trichostrongylus_semiquant <- factor(map$Trichostrongylus_semiquant, levels = c('None','Low','Moderate','High'), ordered = T)

map$Trichuris_semiquant <- droplevels(map$Trichuris_semiquant)
map$Trichostrongylus_semiquant <- droplevels(map$Trichostrongylus_semiquant)

map.douc <- map[douc.ids, ]
map.howler <- map[howler.ids, ]

L6.douc <- data.frame(t(L6.douc))
L6.howler <- data.frame(t(L6.howler))

y <- data.frame(t(L3.douc))
sum(is.na(y))
dim(y)
x <- map.douc[, c(25,28)]
sum(is.na(x))
dim(x)

# create the dataframe that will hold your correlations
cors <- NULL

for (i in colnames(x)) {
  for (j in colnames(y)) {
    a <- as.numeric(x[,i])
    b <- as.numeric(y[,j])
    tmp <- c(i, j, cor(a, b, method = "pearson"), cor.test(a,b, method = "pearson", exact = F)$p.value)   # you can tweak this line to suppress warnings if you like
    cors <- rbind(cors,tmp)
  }
}

cors_stash_b_genus <- cors
# save(cors_stash_b_genus, file = '../data/metabolomics/nmr_cors_genus_stool-post.rdata')
cors<- data.frame(row.names = NULL, cors, stringsAsFactors = FALSE)

colnames(cors) <- c("parasite", "organism", "correlation", "pvalue")  
cors$correlation <- as.numeric(cors$correlation)
cors$pvalue <- as.numeric(cors$pvalue)
cors$qvalue <- p.adjust(cors$pvalue, method = "fdr")
cors$Significance<-cut(cors$qvalue, breaks=c(-Inf, 0.05, 0.25, Inf), label=c("**","*", ""))


# prune for complete cases with significant correlations
dat.c<-cors[complete.cases(cors),]
xkeep<-unique(dat.c$parasite[dat.c$qvalue<=1])  # adjust this q-value for your liking, 1 will show everything, <0.25, only signif.
xassign<-xkeep
dat.w<-dat.c[dat.c$parasite %in% xassign,] 
ykeep <- unique(dat.w$organism[dat.w$qvalue<=0.98]) # adjust this one too as above
dat.m<-dat.w[dat.w$organism %in% ykeep,]


#reorder variables for plotting if more than one observation
order_y <- dat.m %>% select(parasite, organism, correlation) %>% na.omit()
order_y <- order_y %>% spread(parasite, correlation)
order_y <- remove_rownames(order_y)
order_y <- column_to_rownames(order_y, "organism")
yorder <- hclust((dist(1-cor(order_y))/2))$order


# reorder the other variables for plotting
order_x <- dat.m %>% select(organism, parasite, correlation) %>% na.omit()
order_x <- order_x %>% spread(organism, correlation)
order_x <- remove_rownames(order_x)
order_x<- column_to_rownames(order_x, "parasite")
xorder <- hclust((dist(1-cor(order_x))/2))$order

dat.m$parasite <- as.factor(dat.m$parasite)
dat.m$parasite_ordered <- factor(dat.m$parasite, levels(dat.m$parasite)[yorder])

# dat.m$organism <- as.factor(dat.m$organism)
# dat.m$organism <- factor(dat.m$organism, levels(dat.m$organism)[xorder])

dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism), FUN = function(xx) gsub(x = xx, pattern = '.p__', replacement = ';p__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.c__', replacement = ';c__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.o__', replacement = ';o__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.f__', replacement = ';f__', fixed = T)))
dat.m$organism_2 <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) gsub(x = xx, pattern = '.g__', replacement = ';', fixed = T)))

dat.m$organism_short <- unlist(lapply(X=as.character(dat.m$organism_2), FUN = function(xx) tail(strsplit(x=xx, split = as.character(';'), fixed = T)[[1]], n=1)))

dat.m$organism_short <- as.factor(dat.m$organism_short)
dat.m$organism_short_ordered <- factor(dat.m$organism_short, levels(dat.m$organism_short)[xorder])

summary(dat.m$correlation)
summary(dat.m$qvalue)
table(dat.m$Significance)
#plot the correlations
myheatmap <- ggplot(data = dat.m, aes(x=parasite_ordered, y=organism_short_ordered, fill=correlation)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low="yellow", mid="black", high="blue", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Pearson") +
  geom_text(data = dat.m, aes(label=Significance), color="white", size=4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1, color='black'), 
        axis.text.y = element_text(size = 8, color='black'),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.ticks = element_blank(),
        legend.key.size = unit(0.1, "inch"),
        plot.title = element_text(size=9),
        panel.background = element_rect(fill='white')) +
  labs(x = "Parasite", y = "Genus or higher (collapsed to L6)", title='DOUC\nGenus-level correlations\nwith three parasites\n(microscopy egg counts)') +
  # theme(strip.text.y = element_text(angle = 0, face = "italic"), 
  #       strip.background = element_rect(color="white", fill = "white")) +
  coord_fixed(ratio=.75)

myheatmap
ggsave(myheatmap, filename = '../../results/douc_microbiome_parasite_correlations_L3.png', dpi = 300, height = 4.5, width = 3.2)


## NOTHING INTERESTING







####### not using this part currently
# Grp.Pvals = rep(1,nrow(L6.douc)) # Initialize p values for groupwise sig. tests\
# Grp.Corrs = rep(0,nrow(L6.douc)) # Initialize correlation chamber (groupwise)\
# Wld.Pvals = rep(1,nrow(L6.douc)) # Initialize p values for wild vs captive tests\
# for (m.ix in 1:nrow(L6.douc)) {  # Loop through all the rows (picrust pathways)
#   try({
#     ps = polyserial(L6.douc[m.ix,],map.douc$Trichuris_semiquant,ML=T,std.err = T)
#     if (is.na(ps$rho)) next
#     Grp.Corrs[m.ix] = ps$rho
#     Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) 
#     },silent=T)
#   Wld.Pvals[m.ix] = wilcox.test(as.numeric(L6.douc[m.ix,]) ~ (map.douc$Trichuris>10), exact=F)$p.value
#   }
# 
# # Adjust for multiple tests
# Grp.Pvals = p.adjust(Grp.Pvals, method = 'fdr')
# Wld.Pvals = p.adjust(Wld.Pvals, method = 'fdr')
# 
# # make and sort a data frame with these columns
# df = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals, row.names = rownames(L6.douc))
# select = abs(df$Grp.Corrs) > 0.25 | df$Grp.Pvals < 0.25 | df$Wld.Pvals < 0.25
# df = df[select,]
# df = df[order(df$Grp.Corrs),]


#### Approach with pairwise t ####
pairwise.t.pvals <- rep(1,nrow(L6.douc))
for (pw.ix in 1:nrow(L6.douc)) {
  pairwise.t.pvals[pw.ix] <- min(pairwise.t.test(as.numeric(L6.douc[pw.ix,]), g=map.douc$Trichuris_semiquant, p.adjust.method = 'fdr')$p.value, na.rm=T)
}
sum(pairwise.t.pvals < 0.05)
d.res.trichuris <- data.frame(rownames(L6.douc),pairwise.t.pvals, row.names = rownames(L6.douc))
d.L6hits <- rownames(d.res.trichuris[d.res.trichuris$pairwise.t.pvals<0.05,])
# clostridium.pw.douc <- pairwise.t.test(as.numeric(L6.douc[51,]), g=map.douc$Trichuris_semiquant, p.adjust.method = 'fdr')$p.value

beeswarm::beeswarm(as.numeric(L6.douc[d.L6hits[1],]) ~ map.douc$Trichuris_semiquant)
d.L6hits[1]  # Is interesting

# Save this one
pdf(file='../../results/Douc_trichuris_taxa_v2_boxplot.pdf')
for (dt in c(1,2)) {
  boxplot(as.numeric(L6.douc[d.L6hits[dt],]) ~ map.douc$Trichuris_semiquant,
           xlab='Trichuris semi-quantitative index', ylab='CLR abundance',
           main=d.L6hits[dt], cex.main=0.65, cex.axis=0.7)
}
dev.off()

# Now with howler
pairwise.t.pvals <- rep(1,nrow(L6.howler))
for (pw.ix in 1:nrow(L6.howler)) {
  pairwise.t.pvals[pw.ix] <- min(pairwise.t.test(as.numeric(L6.howler[pw.ix,]), g=map.howler$Trichostrongylus_semiquant, p.adjust.method = 'fdr')$p.value, na.rm=T)
}
sum(pairwise.t.pvals < 0.05)
howler.res.trichostrong <- data.frame(rownames(L6.howler),pairwise.t.pvals, row.names = rownames(L6.howler))
h.L6hits <- rownames(howler.res.trichostrong[howler.res.trichostrong$pairwise.t.pvals<0.05,])

beeswarm::beeswarm(as.numeric(L6.howler[h.L6hits[5],]) ~ map.howler$Trichostrongylus_semiquant)
# 1,2,6,8,12 are potentially interesting
# Save these ones
pdf(file='../../results/Howler_trichostrongylus_taxa_v2_boxplots.pdf')
for (ht in 1:length(h.L6hits)) {
  boxplot(as.numeric(L6.howler[h.L6hits[ht],]) ~ map.howler$Trichostrongylus_semiquant,
           xlab='Trichostrongylus semi-quantitative index', ylab='CLR abundance',
           main=h.L6hits[ht], cex.main=0.65, cex.axis=0.7)
}
dev.off()


# Now do it for all of them!
map$all_parasites_ct <- rowSums(map[, 22:32])
map$all_parasites_semiquant <- cut(map$all_parasites_ct, right = F,
                                   breaks=c(-Inf, 1, 10, 100, 1000), label=c("None","Low", "Moderate","Moderate"))
# Changed the high to moderate because there is only 1 with 100+ counts
map$all_parasites_semiquant <- factor(map$all_parasites_semiquant,
                                      levels = c('None','Low','Moderate'), ordered = T)

map <- map[sample.ids, ]
L6_otus.pwt <- data.frame(t(L6_otus[sample.ids, ]))

pairwise.t.pvals <- rep(1,nrow(L6_otus.pwt))
for (pw.ix in 1:nrow(L6_otus.pwt)) {
  pairwise.t.pvals[pw.ix] <- min(pairwise.t.test(as.numeric(L6_otus.pwt[pw.ix,]), g=map$all_parasites_semiquant,
                                                 p.adjust.method = 'fdr')$p.value, na.rm=T)
}
sum(pairwise.t.pvals < (0.002))
all.res <- data.frame(rownames(L6_otus.pwt),pairwise.t.pvals, row.names = rownames(L6_otus.pwt))
all.L6hits <- rownames(all.res[all.res$pairwise.t.pvals<0.001,])

pdf(file='../../results/all_parasites_vs_taxa_boxplots.pdf')
for (ht in 1:length(all.L6hits)) {
  boxplot(as.numeric(L6_otus.pwt[all.L6hits[ht],]) ~ map$all_parasites_semiquant,
          xlab='All parasites semi-quantitative index', ylab='CLR abundance',
          main=all.L6hits[ht], cex.main=0.65, cex.axis=0.7)
}
dev.off()

# Save the full map file
# map.w <- tibble::rownames_to_column(map, var = '#SampleID')
# write.table(map.w, '../../data/nhp_microbe-parasite_map_w-groups_180720.txt', quote=F,
#             sep = '\t', row.names = F)
