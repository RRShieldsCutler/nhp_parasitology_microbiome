# Attempt at doing DE NOVO OTU CLUSTERING using BURST
# 2/14/18

# Things to NOTE:
#   1. project_030 contains only IMP samples (removed FC samples, shared run)
#   2. project_032 was a full IMP run (contains some Longitudinal M1 repeats from 030)
#   3. project_032_redo contains only 3 samples (repeats that failed from 032)
#   4. project_037 contains only 10 samples from IMP (shared run)
#   5. project_040 was a full IMP run

# run Gabe's manicure fasta script which will run shi7 on a light mode:
#       -no quality control
#       -stitches R1 and R2
#       -trims adaptors/primers
#       -filters for desired amplicon length
#       -this produces a filtered fasta file at varying levels of D
#       -and a manicured.fna file (keeps full length UMGC sample names in tact, will be useful for identifying runs later)
cd /project/flatiron2/pj/imp/all_seq_runs/projects/project_040/
bash ../manicure.fasta.sh ./raw/

cd /project/flatiron2/pj/imp/all_seq_runs/projects/project_037/
bash ../manicure.fasta.sh ./raw/

# for all other projects, run on MSI with PBS scripts (see ~/imp_denovo folder)

# union all sequences from all runs that appeared > 100 times, to form our rep_set
cd /project/flatiron2/pj/imp/all_seq_runs/pick_otus/denovo
source activate myR
R
    dirs <- c("/project/flatiron2/pj/imp/all_seq_runs/projects/project_030/D100_filt.fa",
        "/project/flatiron2/pj/imp/all_seq_runs/projects/project_032/D100_filt.fa",
        "/project/flatiron2/pj/imp/all_seq_runs/projects/project_032redo/D100_filt.fa",
        "/project/flatiron2/pj/imp/all_seq_runs/projects/project_037/D100_filt.fa",
        "/project/flatiron2/pj/imp/all_seq_runs/projects/project_040/D100_filt.fa")
    all_reps <- NULL
    for(i in 1:length(dirs)){
        reps <- read.table(dirs[i], comment=">", as.is=T, quote="")
        print(length(reps[,1]))
        all_reps <- union(all_reps, reps[,1])
    }
    write.table(data.frame(paste0(">", 1:length(all_reps)), all_reps), file="rep_set.fa", col.names=F, row.names=F, quote=F, sep="\n")

# concatenate all manicured fasta files into one prior to running BURST
#13752 rep_set.fa
cd /project/flatiron2/pj/imp/all_seq_runs/projects
cat project_030/manicured.fna project_032/manicured.fna project_032redo/manicured.fna project_037/manicured.fna project_040/manicured.fna > manicured.allruns.fna 
wc -l manicured.allruns.fna
# 69354240 manicured.allruns.fna

# create the rep_set required BURST database files
cd /project/flatiron2/pj/imp/all_seq_runs/pick_otus/denovo
# BURST now has automatic creation of all databases (no need for intermediate steps)
burst12 -d -o rep_set.edx -a rep_set.acx -r rep_set.fa -f

# Align representative sequences to PROK database at 0% so that EVERY sequence gets taxonomy assigned
burst12 -q rep_set.fa -r /project/flatiron2/sop/PROK_170704.edx -a /project/flatiron2/sop/PROK_170704.acx -b /project/flatiron2/sop/PROK_170704.tax -n -m CAPITALIST -bs -i 0.0 -f -o ./burst_rep_set.b6
# grab the 1st and last columns of the b6 file
cut -d$'\t' -f 1,13 burst_rep_set.b6 > rep_set.tax

# run BURST at 99% of all manicured sequences against newly created repset
    burst12 -r rep_set.edx -a rep_set.acx -b rep_set.tax -q /project/flatiron2/pj/imp/all_seq_runs/projects/manicured.allruns.fna -o /project/flatiron2/pj/imp/all_seq_runs/pick_otus/denovo/embalmer_output_99.b6 -n -m CAPITALIST -bs -i 0.99 -f -sa
    embalmulate embalmer_output_99.b6 otutable_99.txt taxatable_99.txt GGtrim
    #Parsed 30764141 reads [706 samples, 555 taxa, 6871 refs]. Collating...
    # 30764141/(69354240/2)) == 88.7% (vs 88.37% previously)

# run BURST at 98%
burst12 -r rep_set.edx -a rep_set.acx -b rep_set.tax -q /project/flatiron2/pj/imp/all_seq_runs/projects/manicured.allruns.fna -o /project/flatiron2/pj/imp/all_seq_runs/pick_otus/denovo/embalmer_output_98.b6 -n -m CAPITALIST -bs -i 0.98 -f -sa
embalmulate embalmer_output_98.b6 otutable_98.txt taxatable_98.txt GGtrim
    # Parsed 32437414 reads [706 samples, 604 taxa, 6871 refs]. Collating...
    # 32437414/(69354240/2)) = 93.54%
    # this means, that the remaining 6.5% did not hit our rep_set database at 98%:
    #       --> appeared few times < 100 anyway (else, they would have been a member of the rep_set)
    #       --> or they weren't 98% similar enough to any of the reps
    # ref-based otu picking = 89.7%
    # Transfer otu and taxa files to local machine
        # handed off rep_set.fa to Gabe to build a rep_set.tre!

# HAND-PICK samples from runs to include for further processing
    # run code in map.sample.runs.r to decide which samples to drop and keep between runs, and to rename all samples to their real sample IDs
    # note that this step is somewhat manual
    # transfer newly renamed taxatable_98.txt to teraminx in /clr folder to do the below

# run CLR
    source activate myR3.3
    cd /project/flatiron2/pj/imp/all_seq_runs/pick_otus/denovo/clr
    Rscript ./impute.clr.r -t taxatable_98.txt -f ./clr -d 1000

    # transfer RData to local ~/Dropbox/UMN/KnightsLab/IMP/ANALYSES/ALL_SEQ_RUNS/denovo/clr
    # run transform locally
    Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/embalmer-post-processing/transform.clr.r -i clr

# use CLR for ordination and stats only
# use sqrt rel abundance (of regular taxa tables) for actually showing taxa (beeswarm, taxa summary)

# Run QIIME generation as normal WITH phylogenetic features (Gabe generated one with our rep set)
    # this step is useful for alpha diversity generation, OTU filtering etc, and having relative abundance taxa files
    bash /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/embalmer_to_qiime_output.sh otutable_98.txt taxatable_98.txt ./qiime_files rep_set.tre
        # filter and rarefy depth: 2143
        
# Summarize filtered and formatted taxa table
    # removed first line, and replace UNKNOWN with k__Other
      Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/embalmer-post-processing/gabe.burst.summarize.r -i ./qiime_files/taxa0_s2_f.txt -o ./qiime_files/taxatable
    
# use taxa table to generate alpha diversity (non-phylo since using taxa, no tree available with taxa alone)
    cd /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/ALL_RUNS/denovo/qiime_files
    biom summarize-table -i taxa0_s2.biom -o taxa.stats.txt
    single_rarefaction.py -i taxa0_s2_f.biom -o taxa0_s2_f_rare.biom -d 2143 # use same depth as for OTU
    alpha_diversity.py -i taxa0_s2_f_rare.biom -o alpha.taxa.txt -m chao1,observed_otus,shannon,simpson
