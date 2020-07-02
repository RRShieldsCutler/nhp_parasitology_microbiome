# README for updated parasite OTU processing
## June 2020 and forward
## Robin, Shiv, Jonathan

All of this is synced to GitHub repo

(qiime2-2020.2) [13:10] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table group --i-table otutable-noMitoChloro.qza --p-axis 'sample' --m-metadata-file nhp_microbe-parasite_map_w-groups_200602.txt --m-metadata-column 'AnimalID' --p-mode 'sum' --o-grouped-table animal-sum_otutable_200602.qza
Saved FeatureTable[Frequency] to: animal-sum_otutable_200602.qza

(qiime2-2020.2) [13:22] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table filter-samples --i-table animal-sum_otutable_200602.qza --m-metadata-file wild-howler_animal_parasite_map_200602.txt --o-filtered-table wild-howler_animal-sum_otutable_200602.qza
Saved FeatureTable[Frequency] to: wild-howler_animal-sum_otutable_200602.qza
(qiime2-2020.2) [13:26] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table filter-samples --i-table animal-sum_otutable_200602.qza --m-metadata-file howler_animal_parasite_map_200602.txt --o-filtered-table howler_animal-sum_otutable_200602.qza
Saved FeatureTable[Frequency] to: howler_animal-sum_otutable_200602.qza
(qiime2-2020.2) [13:26] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table filter-samples --i-table animal-sum_otutable_200602.qza --m-metadata-file douc_animal_parasite_map_200602.txt --o-filtered-table douc_animal-sum_otutable_200602.qza
Saved FeatureTable[Frequency] to: douc_animal-sum_otutable_200602.qza

(qiime2-2020.2) [13:28] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table summarize --i-table animal-sum_otutable_200602.qza --m-sample-metadata-file nhp_animal_parasite_map_200602.txt --o-visualization animal_sum_tablesummary.qzv
Saved Visualization to: animal_sum_tablesummary.qzv

(qiime2-2020.2) [13:30] ~/Dropbox/Research/Knights_Lab/knights_box/nhp_projects/parasitology_microbiome/data $ qiime feature-table heatmap --i-table animal-sum_otutable_200602.qza --m-sample-metadata-file nhp_animal_parasite_map_200602.txt --o-visualization animal_sum_heatmap_v1_species.qzv --m-sample-metadata-column 'SpeciesCommonName' --p-color-scheme 'viridis'


# Qiime2 automater for all samples
python ~/bin/qiime2_automater/qiime2_2018.11_automater.py animal-sum_otutable_200602.qza nhp_animal_parasite_map_200602.txt ~/bin/gg97_taxonomy.qza 14825 ~/bin/gg97_q2tree.qza NONE

# Qiime2 automater for howlers
python ~/bin/qiime2_automater/qiime2_2018.11_automater.py ../howler_animal-sum_otutable_200602.qza ../howler_animal_parasite_map_200602.txt ~/bin/gg97_taxonomy.qza 14825 ~/bin/gg97_q2tree.qza NONE

# For doucs
python ~/bin/qiime2_automater/qiime2_2018.11_automater.py ../douc_animal-sum_otutable_200602.qza ../douc_animal_parasite_map_200602.txt ~/bin/gg97_taxonomy.qza 14825 ~/bin/gg97_q2tree.qza NONE

# For wild howlers
python ~/bin/qiime2_automater/qiime2_2018.11_automater.py ../wild-howler_animal-sum_otutable_200602.qza ../wild-howler_animal_parasite_map_200602.txt ~/bin/gg97_taxonomy.qza 14825 ~/bin/gg97_q2tree.qza NONE


