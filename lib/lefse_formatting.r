######### Formatting data for lefse analysis #########

L6_otus <- read.delim('../../data/qiime2_180309/otu_table_L6.txt', header = 1, sep = '\t', row.names = 1, skip=1)
map <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_180720.txt', header=1, row.names = 1, sep = '\t', check.names = F)

sum(L6_otus)

sample.ids <- intersect(colnames(L6_otus), rownames(map))
sample.ids <- sort(sample.ids)
L6_otus <- L6_otus[, sample.ids]
map <- map[sample.ids,]

#CLR
otu.t = t(L6_otus); eps = 0.5
otu.t = otu.t*(1 - rowSums(otu.t==0)*eps/rowSums(otu.t))
otu.t[otu.t==0]=eps
otu.t = sweep(otu.t,1,rowSums(otu.t),'/');
ls = log(otu.t)
otu.t = t(ls - rowMeans(ls))
L6_otus.c = otu.t[,!is.nan(colSums(otu.t))]

# RA

L6.ra <- sweep(L6_otus, 2, colSums(L6_otus), '/')

map.sub <- map[, c("SpeciesCommonName", "all_parasites_semiquant")]

L6.lefse <- data.frame(t(L6.ra))
L6.lefse <- merge(map.sub, L6.lefse, by=0)
L6.lefse <- L6.lefse[, c(3,2,1,4:ncol(L6.lefse))]
colnames(L6.lefse)[1:3] <- c("parasite_load","primate_species","subject_id")
L6.lefse <- tibble::rownames_to_column(data.frame(t(L6.lefse)), var = 'rows')
write.table(L6.lefse, file = '../../data/lefse/l6_ra_subclass.txt', col.names = F, row.names = F, quote = F, sep = '\t')


