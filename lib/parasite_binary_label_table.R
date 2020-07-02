# Making a simple table

map.table <- read.delim('../../data/nhp_microbe-parasite_map_w-groups_180720.txt', header=1, row.names = 1, sep = '\t', check.names = F)

# What's there
#colSums(map.table[,c(22:32)])

map.table$trich.b <- map.table$Trichuris_semiquant != "None"
map.table$tricho.b <- map.table$Trichostrongylus_semiquant != "None"
map.table$strony.b <- map.table$Strongyloides_semiquant != "None"
map.table$Oes.b <- map.table$Oesophagostumum > 0
map.table$Ascaris.b <- map.table$Ascaris > 0 
map.table$fasc.b <- map.table$Fasciola > 0
map.table$digen.b <- map.table$`Digenean trematode eggs` > 0

map.table$trich.b[map.table$trich.b==TRUE] <- "Trichuris"
map.table$tricho.b[map.table$tricho.b==TRUE] <- "Trichostrongylus"
map.table$strony.b[map.table$strony.b==TRUE] <- "Strongyloides"
map.table$Oes.b[map.table$Oes.b==TRUE] <- "Oesophagostumum"
map.table$Ascaris.b[map.table$Ascaris.b==TRUE] <- "Ascaris"
map.table$fasc.b[map.table$fasc.b==TRUE] <- "Fasciola"
map.table$digen.b[map.table$digen.b==TRUE] <- "Digenean trematode"

map.table$trich.b[map.table$trich.b==FALSE] <- NA
map.table$tricho.b[map.table$tricho.b==FALSE] <- NA
map.table$strony.b[map.table$strony.b==FALSE] <- NA
map.table$Oes.b[map.table$Oes.b==FALSE] <- NA
map.table$Ascaris.b[map.table$Ascaris.b==FALSE] <- NA
map.table$fasc.b[map.table$fasc.b==FALSE] <- NA
map.table$digen.b[map.table$digen.b==FALSE] <- NA

map.table2 <- unite(map.table, all_parasites_labels, 39:45, sep = "|", remove = FALSE)
# map.table2$all_parasites_labels <- gsub(x=map.table2$all_parasites_labels, pattern = "|NA|", replacement = "", fixed = T)
map.table2$all_parasites_labels <- gsub(x=map.table2$all_parasites_labels, pattern = "|NA", replacement = "", fixed = T)
map.table2$all_parasites_labels <- gsub(x=map.table2$all_parasites_labels, pattern = "NA|", replacement = "", fixed = T)
map.table2$all_parasites_labels <- gsub(x=map.table2$all_parasites_labels, pattern = "NA", replacement = "", fixed = T)
map.table2$all_parasites_labels[map.table2$all_parasites_labels==""] <- "None"

table(map.table2$all_parasites_labels, map.table2$SpeciesCommonName)
