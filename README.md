# Antimicrobial Peptide
This repository holds the files and code used to predict and annotate antimicrobial peptides in Chiroptera and identify gene gains/losses while also trying to characterize species-specific gene clusters (Castellanos et al. 2023).

The project attempts to get a better understanding of gene-family evolution of antimicrobial peptides, specifically the defensins (α and β) and cathelicidin families across the order Chiroptera.

## Pipeline

## Ortholog identification in bat proteomes

Search was done using [orthofisher](https://github.com/JLSteenwyk/orthofisher#quick-start), which initiaties a HMMER search and further incorporates a best-hit BUSCO pipeline. DEFAs with e-values < 0.001 and 80% best hits were kept, whereas 75% was used for DEFBs.

Refer to the orthofisher documentation to familiarize with the contents of `fasta_args.txt` used in the following loop.

```
#for each line of the fasta_args.txt file.. do
cat $DIR/fasta_args.txt | while read LINE; do

#grab the second column of the file
SPECIES=$(echo ${LINE} | awk '{print $2}')
#create a directory with the name of the species
mkdir ${SPECIES}_DEFA
#move into that directory
(cd ${SPECIES}_DEFA

#copy files to the current directory
cp proteomes/${SPECIES}_nodups.fasta .
cp defa_outgroup.hmm .

echo ${LINE} > ${SPECIES}_dir.txt
sed -i 's/ /\t/g' ${SPECIES}_dir.txt

# Activate conda environment
conda activate hmmer
orthofisher -m $DIR/hmms.txt -f ${SPECIES}_dir.txt -b 0.80

# Activate conda seqkit
conda activate seqkit
seqkit fx2tab orthofisher_output/all_sequences/defa_outgroup.hmm.orthofisher -C c | awk '($4 >= 2)' | seqkit tab2fx > ${SPECIES}_defa_outgroup.ortho.fa);
done
```
You can either write one loop for every AMP family or write a loop that iterates through multiple AMP families.


