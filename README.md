# Antimicrobial Peptide
This repository holds the files and code used to predict and annotate antimicrobial peptides in Chiroptera and identify gene gains/losses while also trying to characterize species-specific gene clusters (Castellanos et al. 2023).

The project attempts to get a better understanding of gene-family evolution of antimicrobial peptides, specifically the defensins (α and β) and cathelicidin families across the order Chiroptera.

# Table of Contents

- [Ortholog identification in bat proteomes](#Ortholog-identification-in-bat-proteomes)
- [AMP probability estimation](#AMP-probability-estimation)
- [Gene search in genomes](#Gene-search-in-genomes)
  - [The AMPlify Pipeline](#The-AMPlify-Pipeline)
  - [Databases](#Databases)
    - [Positive training set](#Positive-training-set)
    - [Negative training set](#Negative-training-set)
  - [Aligning AMPs precursors to the genomes](##Aligning-AMPs-precursors-to-the-genomes)
    - [GMAP](###Installing-and-running-GMAP)
    - [Aligning AMPs proteins to the genomes](##Aligning-AMPs-proteins-to-the-genomes)
  - [MAKER2](##MAKER2)


# Ortholog identification in bat proteomes

Search was done using [orthofisher](https://github.com/JLSteenwyk/orthofisher#quick-start), which initiaties a HMMER search and further incorporates a best-hit BUSCO pipeline. DEFAs with e-values < 0.001 and 80% best hits were kept, whereas 75% was used for DEFBs.

Refer to the orthofisher documentation to familiarize with the contents of `fasta_args.txt` used in the following loop.

```bash
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

# Activate hmmer
conda activate hmmer
orthofisher -m $DIR/hmms.txt -f ${SPECIES}_dir.txt -b 0.80

# Activate seqkit
conda activate seqkit
seqkit fx2tab orthofisher_output/all_sequences/defa_outgroup.hmm.orthofisher -C c | awk '($4 >= 2)' | seqkit tab2fx > ${SPECIES}_defa_outgroup.ortho.fa);
done
```
You can either write one loop for every AMP family or write a loop that iterates through multiple AMP families.

# AMP probability estimation

After quering bat proteomes and keeping sequences based on criteria mentioned in Castellanos et al. (2023), resulting datasets were subjected to an AMP probability test using the [ampir](https://github.com/Legana/AMP_pub) package.

First a training dataset was built following the [How to train your model](https://cran.r-project.org/web/packages/ampir/vignettes/train_model.html) vignette. After the training was completed, a subsampling of the dataset (50 random seqs) is obtained and the model is tested:


```r
# Create train and test sets
trainIndex <-createDataPartition(y=bats_features$Label, p=.7, list = FALSE)
bats_featuresTrain <-bats_features[trainIndex,]
bats_featuresTest <-bats_features[-trainIndex,]

#train model
my_bat_svm_model <- train(Label~.,
                          data = bats_featuresTrain[,-1],
                          method="svmRadial",
                          trControl = trctrl_prob,
                          preProcess = c("center", "scale"))

#test model in a subsample of the dataset
test_bat_AMPs <- predict_amps(bat_test_set, min_len = 5, model = my_bat_svm_model)
```

# GENE SEARCH IN GENOMES
## The AMPlify Pipeline

This workflow is taken from [Li et al. 2022](https://bmcgenomics-biomedcentral-com.lib-e2.lib.ttu.edu/articles/10.1186/s12864-022-08310-4#Sec9), and has been modified to adjust to our needs.

## Databases

### Positive training set

Two databases were created for the **positive training set**: [APD3](https://aps.unmc.edu/) containing 3172 AMP aminoacid sequences from a variety of Phyla. The latest version (Jan. 2020) of entire database was downloaded [here](https://aps.unmc.edu/downloads). It also includes DEFA, DEFB, and CTHL fasta files that have been manually curated from NCBI and ENSEMBL of my outgroup. Finally, the manually curated database of DEFA,DEFB and CTHl of bats proteomes. **Total 3694**.

To achieve this, the fasta file was first edited to avoid conflicts due to headers.

```bash
#After I merged the bats, vertebrate and APD3 in a single file, I used

seqkit fx2tab apd_bats_laura.fasta -l -i -H -C BJOUXZ | awk '{ if ($3 <= 200 && $4 < 1 ) {print $1 "\t" $2} }' | seqkit tab2fx | seqkit rmdup -s -o clean_apd_bats_laura.fasta -D duplicated_apd_bats_laura.txt
```

### Negative training set

For the **negative traning set** I downloaded the whole data file for the Reviewed SwissProt from [here](https://www.uniprot.org/help/downloads). I applied a custom script to filter based on the same keywords as in the [AMPlify paper](https://bmcgenomics-biomedcentral-com.lib-e2.lib.ttu.edu/articles/10.1186/s12864-022-08310-4#Sec9)

```bash
grep '^AC\|^KW' uniprot_sprot.dat | awk '{printf "%s%s", /^KW/?OFS:(NR>1)?"\n":"", $0} END{print ""}' | \
sed "s/;/\t/g" | grep -iv "antimicrobial\|antibiotic\|antibacterial\|antiviral\|antifungal\|antimalarial\|antiparasitic\|anti-protist\|anticancer\|defense\|defensin\|cathelicidin\|histatin\|bacteriocin\|microbicidal\|fungicide" | \
awk -F'\t' '$3 != ""'| awk '{print $2}' > list_uniprot_non_amps.list
```
For the **potential AMPs dataset** I used a variation of the code above

```bash=
grep '^AC\|^KW' uniprot_sprot.dat | awk '{printf "%s%s", /^KW/?OFS:(NR>1)?"\n":"", $0} END{print ""}' | \
sed "s/;/\t/g" | grep -i "antimicrobial\|antibiotic\|antibacterial\|antiviral\|antifungal\|antimalarial\|antiparasitic\|anti-protist\|anticancer\|defense\|defensin\|cathelicidin\|histatin\|bacteriocin\|microbicidal\|fungicide" | \
awk '{print $2}' > list_uniprot_potential_amps.list
```
Once I got that file then I downloaded it using a modified script from one I found [online]().

```bash
cat list_uniprot_non_amps.list | \
xargs -n 1 -P 32 -I % curl -s 'https://rest.uniprot.org/uniprotkb/search?query=%&format=fasta' \
>> uniprot_non_amps.fasta
```
I didn´t edit the headers as for the [positive training set](https://hackmd.io/iah38oSGSu-_1qGy9XMVZw?both#-Positive-training-set), however, I deleted sequences that contained non-standard aminoacids (B,J,O,U,X and Z), sequences larger than 200 AAs, and lastly deleted duplicates with a single command:

```bash
seqkit fx2tab uniprot_potential_amps.fasta -l -i -H -C BJOUXZ | \
awk '{ if ($3 <= 200 && $4 < 1 ) {print $1 "\t" $2} }' | seqkit tab2fx | \
seqkit rmdup -s -o clean_uniprot_potential_amps.fasta -D duplicated_potential.txt
```
This resulted in **132,180 aminoacid seqs**

## Aligning AMPs precursors to the genomes

### GMAP

GMAP is available [here](https://github.com/juliangehring/GMAP-GSNAP). And can be run like:

```bash
BAT=<NAME>
GMAP=/home/frcastel/software/gmap-2021-12-17/bin/
GENOMES=/lustre/scratch/frcastel/defensins/dataset/amplify/genomes/
CDNA=/lustre/scratch/frcastel/defensins/dataset/amplify/evidence/

#Build gmap database
perl $GMAP/gmap_build -d $BAT -k 15 $GENOMES/$BAT -s none

#Align cDNA from mammals to gmap database
$GMAP/gmap -d $BAT -A --max-intronlength-ends=200000 -f gff3_gene -O -n 20 --nofails -B 3 -k 15 -t 10 $CDNA/mammalia_cDNA_AMPs.fasta

#Get scaffolds which contain hits from gmap alignment
grep -A 1000000 "gff-version" "$BAT"_gmap_gff3.*.out > maker/gff/"$BAT"_gmap.gff3
```

## MAKER2

