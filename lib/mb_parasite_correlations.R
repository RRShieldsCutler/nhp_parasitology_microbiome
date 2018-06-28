# Correlations

L6_otus <- read.delim('../../data/qiime2_work/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
map <- read.delim('../../data/nhp_microbe_parasite_mapfile_v1.txt', header=1, row.names = 1, sep = '\t', check.names = F)

sample.ids <- intersect(colnames(L6_otus), rownames(map))
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

# Or, for heatmap

library(tidyverse)

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
