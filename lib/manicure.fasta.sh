###########################################
## DE NOVO DRAFT PROTOCOL - PREP
# Reqs.: ninja_filter (ninja_filter_linux)
# From Gabe 2/13/18
###########################################

INPUTDIR=$1
TOO_SHORT=250 # note, pass these into line 20 if possible!
TOO_LONG=255

# 1. SHI7 without quality control. (optionally, learn beforehand for adaptor determination)
shi7.py -i $INPUTDIR -o shi7 --adaptor Nextera --allow_outies F --filter_qual 1 --trim_qual 1

# 2. Filter for and remove primers. 
grep -B1 --no-group-separator '.*GCCGCGGTAA.*ATTAGA.ACCC.*' shi7/combined_seqs.fna | sed 's/.*GCCGCGGTAA//' | sed 's/ATTAGA.ACCC.*//' > combined_seqs_np.fna

# 3. Filter for desired amplicon length(s)
  #(optional: check lengths):
  awk "!(NR % 2) {print length}" combined_seqs_np.fna | sort | uniq -c | sort -rh | head -50
awk "{if (!(NR % 2) && length > $TOO_SHORT && length < $TOO_LONG) {print x; print \$0}; x=\$0}" combined_seqs_np.fna > manicured.fna

#3.1 To see the per-sample alignment rate by querying the original fasta.
grep '^>' combined_seqs_np.fna | cut -f1 -d_ | sort | uniq -c | sort -hr
grep '^>' manicured.fna | cut -f1 -d_ | sort | uniq -c | sort -hr


# 4. Create some "seed sets" to determine proper ref set cutoff (refine if desired)
for i in 002 003 005 010 025 050 100; do ninja_filter_linux manicured.fna D$i D $i && rm D$i.db; done
