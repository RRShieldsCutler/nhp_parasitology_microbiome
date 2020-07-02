## Aggregate the data in the mapping file by Animal ID

animal.meta <- read.delim('../../data/nhp_per-animal_metadata.txt',
                     header = T, sep = '\t', check.names = F)
worm.cts <- read.delim('../../data/nhp_animal_parasite_counts.txt',
                       header = T, sep = '\t', check.names = F)
dim(worm.cts)
worm.mean.counts <- aggregate(worm.cts[,2:12], list(worm.cts$AnimalID), mean)
hist(worm.mean.counts$all_parasites_ct)
colnames(worm.mean.counts)[1] <- "AnimalID"
worm.mean.counts$all_parasites_presence <- worm.mean.counts$all_parasites_ct > 0
table(worm.mean.counts$all_parasites_presence)

nrow(worm.mean.counts[worm.mean.counts$all_parasites_ct == 0,])
nrow(worm.mean.counts[worm.mean.counts$all_parasites_ct > 0 & worm.mean.counts$all_parasites_ct <= 10,])
nrow(worm.mean.counts[worm.mean.counts$all_parasites_ct > 10,])

worm.mean.counts$all_parasites_presence[worm.mean.counts$all_parasites_presence==TRUE] <- "positive"
worm.mean.counts$all_parasites_presence[worm.mean.counts$all_parasites_presence==FALSE] <- "negative"

worm.mean.counts$all_parasites_semiquant <- cut(worm.mean.counts$all_parasites_ct,
                                                breaks = c(-Inf,0,10,Inf), right=T,
                                                labels = c("none","low","high"))
table(worm.mean.counts$all_parasites_semiquant)

animal.meta <- animal.meta[!duplicated(animal.meta$AnimalID),]

nhp.parasite.meta <- merge(worm.mean.counts, animal.meta, by='AnimalID')

colnames(nhp.parasite.meta)[1] <- "sampleid"  # Because QIIME2 will want it this way

write.table(nhp.parasite.meta, file = '../../data/nhp_animal_parasite_map_200602.txt', quote = F, sep = '\t', row.names = F)
# make the howler and douc ones

howler.parasite.meta <- nhp.parasite.meta[nhp.parasite.meta$SpeciesCommonName == "mantled-howling-monkey",]
douc.parasite.meta <- nhp.parasite.meta[nhp.parasite.meta$SpeciesCommonName == "red-shanked-douc",]
write.table(howler.parasite.meta, file = '../../data/howler_animal_parasite_map_200602.txt', quote = F, sep = '\t', row.names = F)
write.table(douc.parasite.meta, file = '../../data/douc_animal_parasite_map_200602.txt', quote = F, sep = '\t', row.names = F)
wild.howler.parasite.meta <- howler.parasite.meta[howler.parasite.meta$CaptiveWild == 'Wild',]
write.table(wild.howler.parasite.meta, file = '../../data/wild-howler_animal_parasite_map_200602.txt', quote = F, sep = '\t', row.names = F)
dim(wild.howler.parasite.meta)

