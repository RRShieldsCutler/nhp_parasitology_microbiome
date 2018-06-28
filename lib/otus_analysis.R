## OTU table processing
# QC by shi7
# shi7.py -i raw_fastq/ -o shi7_stitch -t 48 --adaptor Nextera --trim_qual 34 --filter_qual 37 -s .,1 --min_overlap 233 --max_overlap 281 --allow_outies False
# Alignment with BURST12
# burst12 -r /project/flatiron2/sop/gg97.edx -a /project/flatiron2/sop/gg97.acx -b /project/flatiron2/sop/gg97.tax -q shi7_stitch/combined_seqs.fna -o douc_howler_stitch_gg_95.b6 -n -m CAPITALIST -bs -i 0.95 -fr -sa -t 48


library(ggplot2)
library(vegan)
library(reshape2)

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

