[![DOI](https://zenodo.org/badge/660285082.svg)](https://zenodo.org/badge/latestdoi/660285082)

# Antimicrobial Peptides Annotation and Prediction in Chiroptera
This repository holds the code used to predict and annotate antimicrobial peptides in 20 species of the Order Chiroptera, and to identify gene gains/losses. The main resulting files of this pipeline are included in the `results` folder.

This project attempts to get a better understanding of gene-family evolution of antimicrobial peptides, specifically the defensins (α and β) and cathelicidin families across the order Chiroptera.

# Table of Contents

- [Software requirements](#Software-requirements)
- [Ortholog identification in bat proteomes](#Ortholog-identification-in-bat-proteomes)
- [AMP probability estimation](#AMP-probability-estimation)
- [Gene search in genomes](#Gene-search-in-genomes)
  - [The AMPlify Pipeline](#The-AMPlify-Pipeline)
  - [Databases](#Databases)
    - [Positive training set](#Positive-training-set)
    - [Negative training set](#Negative-training-set)
    - [Potential AMPS](#Potential-AMPs-dataset)
  - [Aligning AMPs precursors to the genomes with GMAP](#Aligning-AMPs-precursors-to-the-genomes-with-GMAP)
  - [MAKER2](#MAKER2)
- [Functional Annotation](#functional-annotation)
  - [Extract defensins and cathelicidins from Interpro annotation](#Extract-defensins-and-cathelicidins-from-Interpro-annotation)
  - [Bat defensins/cathelicidins concatenation](#Bat-defensins/cathelicidins-concatenation)
  - [Extract any other AMPs from Interpro annotation](#Extract-any-other-AMPs-from-Interpro-annotation)
- [Final AMPs Prediction](#Final-AMPs-Prediction)
  - [Defensins and cathelicidins prediction](#Defensins-and-cathelicidins-prediction)
  - [Any other AMP prediction](#Any-other-AMP-prediction)
- [Functional annotations stats](#Functional-annotations-stats)
  - [All other AMPs](#All-other-AMPs)
- [Extract genes from gffs](#Extract-genes-from-gffs)
  - [Overview stats of genes](#Overview-stats-of-genes)
  - [Extract genes](#Extract-genes)
- [Annotate proteins with defensin subfamilies](#Annotate-proteins-with-defensin-subfamilies)

- [Citation](#Citation)

# Software requirements

You should have the following software installed in your system:

- [MAKER2](http://www.yandell-lab.org/software/maker.html)
- [GMAP](http://research-pub.gene.com/gmap/)
- [orthofisher](https://jlsteenwyk.com/orthofisher/install/index.html)
- [Interproscan](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [seqtk](https://github.com/lh3/seqtk)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)

The installation of MAKER2 and its dependencies is explained in the [MAKER2 documentation](http://www.yandell-lab.org/software/maker.html). The installation of GMAP is detailed [here](https://github.com/juliangehring/GMAP-GSNAP).
The rest of the software can be installed using [conda](https://docs.conda.io/en/latest/) environments.

# Ortholog identification in bat proteomes

Search was done using [orthofisher](https://jlsteenwyk.com/orthofisher/), which initiaties a HMMER search and further incorporates a best-hit BUSCO pipeline. DEFAs with e-values < 0.001 and 80% best hits were kept, whereas 75% was used for DEFBs.

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
cp proteomes/${SPECIES}.fasta .
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

First a training dataset was built following the [How to train your model](https://cran.r-project.org/web/packages/ampir/vignettes/train_model.html) vignette. After the training was completed, a subsampling of the dataset (70 random sequences) is obtained and the model is tested:

```r
# Create train and test sets
trainIndex <- createDataPartition(y = bats_features$Label, p = 0.7, list = FALSE)
bats_featuresTrain <- bats_features[trainIndex,]
bats_featuresTest <- bats_features[-trainIndex,]

#train model
my_bat_svm_model <- train(Label~.,
                          data = bats_featuresTrain[,-1],
                          method ="svmRadial",
                          trControl = trctrl_prob,
                          preProcess = c("center", "scale"))

#test model in a subsample of the dataset
test_bat_AMPs <- predict_amps(bat_test_set, min_len = 5, model = my_bat_svm_model)
```

# GENE SEARCH IN GENOMES
## The AMPlify Pipeline

This workflow is taken from [Li et al. 2022](https://bmcgenomics-biomedcentral-com.lib-e2.lib.ttu.edu/articles/10.1186/s12864-022-08310-4#Sec9), and has been modified and adjusted to our needs.

## Databases

### Positive training set

The **training set** included:
- [APD3 database](https://aps.unmc.edu/) containing 3,172 AMP aminoacid sequences from a variety of Phyla. The latest version (Jan. 2020) of entire database was downloaded [here](https://aps.unmc.edu/downloads).
- DEFA, DEFB, and CTHL fasta files that were manually curated from NCBI and ENSEMBL of the outgroup chosen for this work.
- Manually curated database of DEFA,DEFB and CTHl of bats proteomes.

After I concatenated the file that contained the bats, vertebrate and APD3 sequences (3,694 proteins), I cleaned it using [seqkit](https://bioinf.shenwei.me/seqkit/) as follows:

```bash
seqkit fx2tab apd_bats_laura.fasta -l -i -H -C BJOUXZ | awk '{ if ($3 <= 200 && $4 < 1 ) {print $1 "\t" $2} }' | seqkit tab2fx | seqkit rmdup -s -o clean_apd_bats_laura.fasta -D duplicated_apd_bats_laura.txt
```

### Negative training set

For the **negative traning set** I downloaded the Reviewed SwissProt database from [here](https://www.uniprot.org/help/downloads). I applied a custom script to filter the sequences using the same keywords as in the [AMPlify paper](https://bmcgenomics-biomedcentral-com.lib-e2.lib.ttu.edu/articles/10.1186/s12864-022-08310-4#Sec9)

```bash
grep '^AC\|^KW' uniprot_sprot.dat | awk '{printf "%s%s", /^KW/?OFS:(NR>1)?"\n":"", $0} END{print ""}' | \
sed "s/;/\t/g" | grep -iv "antimicrobial\|antibiotic\|antibacterial\|antiviral\|antifungal\|antimalarial\|antiparasitic\|anti-protist\|anticancer\|defense\|defensin\|cathelicidin\|histatin\|bacteriocin\|microbicidal\|fungicide" | \
awk -F'\t' '$3 != ""'| awk '{print $2}' > list_uniprot_non_amps.list
```

### Potential AMPs dataset

For the **potential AMPs dataset** I used a variation of the code above:

```bash
grep '^AC\|^KW' uniprot_sprot.dat | awk '{printf "%s%s", /^KW/?OFS:(NR>1)?"\n":"", $0} END{print ""}' | \
sed "s/;/\t/g" | grep -i "antimicrobial\|antibiotic\|antibacterial\|antiviral\|antifungal\|antimalarial\|antiparasitic\|anti-protist\|anticancer\|defense\|defensin\|cathelicidin\|histatin\|bacteriocin\|microbicidal\|fungicide" | \
awk '{print $2}' > list_uniprot_potential_amps.list
```
Once I got the names of the sequences that are potential AMPS, I downloaded them using this script modified from one I found [online]().

```bash
cat list_uniprot_non_amps.list | \
xargs -n 1 -P 32 -I % curl -s 'https://rest.uniprot.org/uniprotkb/search?query=%&format=fasta' \
>> uniprot_non_amps.fasta
```
I didn´t edit the headers as in the [positive training set](https://hackmd.io/iah38oSGSu-_1qGy9XMVZw?both#-Positive-training-set), however, I deleted sequences that contained non-standard aminoacids (B,J,O,U,X and Z), sequences larger than 200 AAs, and lastly deleted duplicates:

```bash
seqkit fx2tab uniprot_potential_amps.fasta -l -i -H -C BJOUXZ | \
awk '{ if ($3 <= 200 && $4 < 1 ) {print $1 "\t" $2} }' | seqkit tab2fx | \
seqkit rmdup -s -o clean_uniprot_potential_amps.fasta -D duplicated_potential.txt
```
The result is a database of **132,180 aminoacid sequences**.

## Aligning AMPs precursors to the genomes with GMAP

GMAP was ran in this pipeline by submitting an individual script for each species. You can replace the `<NAME>` of each species using a `sed` command and iterating the names of the species from a text file.

```bash
BAT=<NAME> #Insert name of the species here
GMAP=/your/path/to/gmap/bin/
GENOMES=/your/path/to/genomes/
CDNA=/your/path/to/cdna/

#Build gmap database
perl $GMAP/gmap_build -d $BAT -k 15 $GENOMES/$BAT -s none

#Align cDNA from mammals to gmap database
$GMAP/gmap -d $BAT -A --max-intronlength-ends=200000 -f gff3_gene -O -n 20 --nofails -B 3 -k 15 -t 10 $CDNA/mammalia_cDNA_AMPs.fasta

#Get scaffolds which contain hits from gmap alignment
grep -A 1000000 "gff-version" "$BAT"_gmap_gff3.*.out > maker/gff/"$BAT"_gmap.gff3
```

## MAKER2

Running MAKER2 is straightfoward after the installation is done. The major changes performed in this work were in the MAKER Behavior Options section in the `maker_opts.ctl` file, as follows:

```bash
#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=1000 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=10 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

# Functional annotation

I used [interproscan](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html) in a conda environment along with opendjk version 11 for this step. The `bat_list_amps.txt` in this section contains the abbreviations for the species I am working with.

```bash
conda activate interpro

for i in $(cat bat_list_amps.txt)
do
interproscan.sh -i $i/"$i".maker.output/curated_"$i"_proteins.fasta -f tsv -dp --cpu 10 -goterms -o $i/"$i".maker.output/curated_"$i"_proteins.tsv
done
```

## Extract defensins and cathelicidins from Interpro annotation
Once Intepro was run, I extracted α-β defensins, and cathelicidins from the `curated_"$i"_proteins.tsv` that I got in the previous step:

```bash
conda activate seqkit

for i in $(cat bat_list_amps.txt)
do
grep -i "Beta_defensin\|Beta-defensin" "$i"/"$i".maker.output/curated_"$i"_proteins.tsv | cut -f1 | uniq > "$i"/"$i".maker.output/"$i"_beta_defensins_headers.txt
seqtk subseq "$i"/"$i".maker.output/curated_"$i"_proteins.fasta "$i"/"$i".maker.output/"$i"_beta_defensins_headers.txt > "$i"/"$i".maker.output/"$i"_beta_defensins.fasta
#Get alpha defensins files
grep -i "Alpha_defensin\|Alpha-defensin" "$i"/"$i".maker.output/curated_"$i"_proteins.tsv | cut -f1 | uniq > "$i"/"$i".maker.output/"$i"_alpha_defensins_headers.txt
seqtk subseq "$i"/"$i".maker.output/curated_"$i"_proteins.fasta "$i"/"$i".maker.output/"$i"_alpha_defensins_headers.txt > "$i"/"$i".maker.output/"$i"_alpha_defensins.fasta
#Get cathelicidin files
grep -i "cath" "$i"/"$i".maker.output/curated_"$i"_proteins.tsv | cut -f1 | uniq > "$i"/"$i".maker.output/"$i"_cathelicidins_headers.txt
seqtk subseq "$i"/"$i".maker.output/curated_"$i"_proteins.fasta "$i"/"$i".maker.output/"$i"_cathelicidins_headers.txt > "$i"/"$i".maker.output/"$i"_cathelicidins.fasta
done
```
## Bat defensins/cathelicidins concatenation

I ran the following command to incorporate the type of AMP and species name in the headers to make it easier for some analyses. However, I removed them later when dealing with the MAKER gff files.

```bash
for i in $(cat bat_list_amps.txt)
do
# This is just for the defb, do the same for the other families
sed "s/^>/>defb_"$i"_/g" "$i"/"$i".maker.output/"$i"_beta_defensins.fasta >> bats.predicted.maker.beta.defensins.fasta
done
```

## Extract any other AMPs from Interpro annotation

I filtered out the defensins and cathelicidins using the headers I extracted in the previous step, from the file that contains all the proteins annotated with MAKER in each species. This will allow me to use ampir to predict which of these are in fact AMPs.

```bash
for bat in $(cat bat_list_amps.txt)
do
DIR=/path/to/MAKER/results/"$bat"/"$bat".maker.output/
mkdir -p "$DIR"/all_other_amps/
cat "$DIR"/*headers.txt > "$DIR"/all_other_amps/"$bat"_def_cath_headers.txt
grep -Fvwf "$DIR"/all_other_amps/"$bat"_def_cath_headers.txt "$DIR"/curated_"$bat"_proteins.tsv | awk '{print $1}' | sort | uniq > "$DIR"/all_other_amps/"$bat"_all_other_amps_headers.txt
seqtk subseq "$DIR"/"$bat".all.maker.proteins.fasta "$DIR"/all_other_amps/"$bat"_all_other_amps_headers.txt > "$DIR"/all_other_amps/"$bat"_all_other_amps.fasta
sed -i "s/^>/"$bat"_/g" "$DIR"/all_other_amps/"$bat"_all_other_amps.fasta
done
```
Then I concatenated all these files into a single `bats_all_other_amps.fasta` file for a final prediction round.

# Final AMPs Prediction

## Defensins and cathelicidins prediction

To make the defensins/cathelicidins prediction of the proteins obtained in MAKER, I subsetted the potential amps dataset created [here](###Potential-AMPs-dataset), by only selecting the defensins and cathelicidins sequences first:

```bash
grep ">" uniprot_potential_amps.fasta | grep -i "cathel\|defens" | sed 's/>//g' > defens_cath_headers.txt && seqtk subseq uniprot_potential_amps.fasta defens_cath_headers.txt > uniprot_potential_defens_cath.fasta
```
And then I kept sequences greater than 6 AAs long:

```bash
seqkit fx2tab uniprot_potential_defens_cath.fasta -l | awk '{ if ($3 >= 6 ) {print $1 "\t" $2} }' | seqkit tab2fx > uniprot_potential_defens_cath_over6aas.fasta
```
Now, to obtain the lengths of the sequences in the datasets I used:

```bash
seqkit fx2tab uniprot_potential_defens_cath_over6aas.fasta -l | awk '{print $1"\t"$3}' > uniprot_potential_defens_cath_over6aas_length.txt
seqkit fx2tab uniprot_non_amps.fasta -l | awk '{print $1"\t"$3}' > uniprot_non_amps_length.txt
```
Finally, I approximated the peptide length distribution of the **negative dataset** to the **positive dataset** using the custom R script shown below.

```R
# Import data
positive <- read.table("uniprot_potential_defens_cath_over6aas_length.txt",
                       header = T, sep = "\t") %>% na.omit()
negative <- read.table("uniprot_non_amps_length.txt",
                       header = T, sep = "\t") %>% na.omit()

# Loop to approximate the negative dataset to the positive dataset
iterations <- seq(1, 210, 10)
final_distribution <- c()

for (i in 1 : length(iterations)){

tmp <- negative %>% filter(length >= iterations[i] & length < iterations[i+1])
tmp <- tmp[sample(nrow(tmp),
                 length(which(positive$length>= iterations[i] & positive$length < iterations[i+1]))), ]

message("sequence length ", iterations[i])
print(length(tmp$length))
final_distribution <- bind_rows(final_distribution,tmp)
}

write.table(final_distribution, "distribution_non_amps.txt",
            quote = FALSE, row.names = F, col.names = F)
```
The resulting distributions are shown below:

![](https://i.imgur.com/9VET2Yp.jpg)

For alpha-beta defensins and cathelicidins prediction, I repeated the same steps as explained [previously](#AMP-probability-estimation). However, I varied the dataset used as the training dataset, using the same that I built for the AMPlify pipeline, but keeping only sequences that contain more than **10 aminoacids** and I fixed the length distribution as well.

I made a final multisequence alignment and got rid of sequences that did not have the 4-6 cystein motif. So the dataset went from 337 predicted proteins to 333 for the **final defensins+cathelicidins dataset**.

This file is found in the `fasta_files` folder under the name `predicted_bat_defensins_cathelicidins.fasta`.

## Any other AMP prediction

To estimate the probability of all the other AMPs besides defensins and cathelicidins, I used the built-in `predict_mature` function in ampir using the `bats_all_other_amps.fasta` file I generated above.

```r
bats_all_other_amps <- read_faa(file="bats_all_other_amps.fasta")
bats_all_other_amps <- remove_nonstandard_aa(bats_all_other_amps)

# Predict with built-in model
bat_other_amps_mature <- predict_amps(bats_all_other_amps, model = "mature")
hist(bat_other_amps_pred$prob_AMP)

## Verify
bat_other_amps_mature %>%
  ggplot() +
  geom_density(aes(x=prob_AMP)) +
  geom_point(aes(y=0,x=prob_AMP)) + xlim(0,1)

## Filter
filtered_all_other_amps_mature <- bat_other_amps_mature %>% filter(prob_AMP >= 0.8) %>%
  select (seq_name, seq_aa, prob_AMP)

## Plot
filtered_all_other_amps_mature %>%
  ggplot() +
  geom_density(aes(x=prob_AMP))

## Save file
df_to_faa(filtered_all_other_amps_mature, "filtered.ampir.bats.all.other.amps.fasta")
```
This file is found in the `fasta_files` folder under the name `predicted_bat_any_other_AMP.fasta`.

# Functional annotations stats

## All other AMPs

First I obtained the headers of the proteins predicted by ampir:

```bash
# Get tsv files of all the bats
for bat in $(cat bat_list_amps.txt)
do
grep ">${bat}" filtered.ampir.bats.all.other.amps.fasta | sed "s/>${bat}_//g" | awk '{print $1}' >> filtered.ampir.bats.all.other.amps.headers.txt
awk -v species="$bat" '{print species "\t" $0}' "$bat"/"$bat".maker.output/curated_"$bat"_proteins.tsv >> bats_curated_proteins.tsv
done

#Get the annotations for the predicted proteins
grep -Fwf filtered.ampir.bats.all.other.amps.headers.txt bats_curated_proteins.tsv > all_other_amps_functional_annotation.tsv
```

I also got a list of AMP names from the APD3 database to compare with my results.

```bash
grep ">0" apd_and_vertebrates_AMPs.fasta | sed "s/|/\t/g ; s/(/\t/g; s/,/\t/g" | awk '{print $2}' | sort | uniq > list_of_APD3_amps_names.txt
```

Finally, I got rid of proteins that were not AMPs and counted predicted AMP families per species with the code below and saved it in the `interpro_results.tsv` file.

```bash
for bat in $(cat bat_list_amps.txt)
do
column -t -s $'\t' "$bat"/"$bat".maker.output/curated_"$bat"_proteins.tsv | grep -v MobiDBLite | grep -v PRINTS | awk '{print $1 "\t" $6 "\t" $7}' | sed "s/,//g"  | awk '{ if ($3 >= 0 ) {print $1 "\t" $2} }' | awk 'NR==1 {print $0}; NR>1 {if(cat[$1])cat[$1]=cat[$1]", "$2; else cat[$1]=$2;}; END{j=1; for (i in cat) print i, cat[i]}' | sed "s/^/"$bat"\t/g"  | sed "s/,/\t/g" | column -t -s $'\t' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' >> interpro.results.tsv
done
```
# Extract genes from gffs

Extracting $\alpha$, $\beta$-defensin, and cathelicidin genes (nucleotides) from the gff files that MAKER creates was necessary for the TE analysis.

## Overview stats of genes

I used the `predicted_bat_defensins_cathelicidins.fasta` file from the ampir prediction and the MAKER gffs to have an overview of the lengths of the genes. The bash script to estimate gene length and number of exons is shown below:

```bash
#Small edits were done for the beta and cathelicidin genes
for i in $(cat bat_list_amps.txt)
do
while read -r line
do
  gene_name=$(echo "$line" | awk '{print $2}')
  count=$(grep "$gene_name" "$i".gff | awk '$3 == "exon"' | wc -l)
  echo -e "$line\t$count"
done < <(grep -Fwf "$i"_alpha_defensins_headers.txt "$i".gff | sed "s/=/\t/g" | awk -v species="$i" '$3 == "gene" {print species"\t"$11"\t"$5-$4}') >> bats_alpha_defensins_genes_stats.txt
done
```
I had to do some manual editing in the resulting file because of a bug where some exon counts were >3 due to grep was not matching exact patterns even though I was using the `-Fw` flags.

In the $\beta$-defensins I manually edited one header because of the same issue found in the $\alpha$. Also the `maker-AnoCau_scaffold_2984-exonerate_protein2genome-gene-0.1` was deleted because it was not properly annotated by MAKER and had two concatenated defensin sequences in one.

Finally, for the cathelicidins I deleted one sequence because it lacked the cystein residues and another small gene of tsaurophihla.

## Extract genes

The headers and the `./Gff3ToBed.sh` script I created for this purpose will only retrieve the genes that are >200 bp.

```bash
#!/bin/bash

# Usage: ./Gff3ToBed.sh headers.txt input.gff3 output.bed

# Check if input file exists
if [ ! -f $1 ] && [ -f $2 ]; then
    echo "Input files do not exist."
    exit 1
fi

# Check if output file exists
if [ -f $3 ]; then
    echo "Output file already exists. Overwriting.."
    rm $3
fi

# Look for the header in the gff and continue with the code
# only if the gene is bigger than 200 bp.
# Then print some columns of the gff and include the length of
# the gene in the last column of the bed file.
grep -Fwf "$1" "$2" | awk 'BEGIN {FS="\t"; OFS="\t"}
     $2!~/^#/ && $3=="gene" && $5-$4>=200 {
         split($9,a,";");
         split(a[1],b,"=");
         print $1,$2,$3,$4-1,$5,$6,$7,b[2],$5-$4}' > "$3.tmp"

# Rename temporary file to output file
mv "$3.tmp" "$3"
echo "Conversion done!"
```

Then I used bedtools to retrieve the genes from the genomes. The script requires the paths to several files created in this tutorial which will vary based on your working directory. Furthermore, the `bedtools getfasta` command requires the path to the **softmasked genome** of each species.

```bash
conda activate bedtools

#This is  for the cathelicidins but I just had to change the names of the subfamilies accordingly
for i in $(cat bat_list_amps.txt)
do
./Gff3ToBed.sh "$i"_cathelicidins_headers.txt "$i".gff "$i"_cathelicidins_genes.bed
bedtools getfasta -fi /path/to/softmasked/genome/"$i" -bed "$i"_cathelicidins_genes.bed -fo "$i"_cathelicidins_genes.fasta
done
```

# Annotate proteins with defensin subfamilies

The Interpro output doesn't exactly states which defensin subfamilies I may have retrieved, thus, I blasted the dataset I used to predict AMPs in ampir to get functional names of the proteins.

First I have to copy the files to the directory and get rid of any prefix I gave to the `filtered.ampir.bats.maker.predicted.amps.fasta` file.

```bash
# Copy files I need for this pipeline
cp /your/directory/to/filtered.ampir.bats.maker.predicted.amps.fasta .
cp /your/directory/to/uniprot_potential_defens_cath_over6aas.fasta .

# Delete prefixes
for i in $(cat bat_list_amps.txt)
do
sed -i "s/^>defa_"$i"_\|^>defb_"$i"_\|^>cath_"$i"_/>/g" filtered.ampir.bats.maker.predicted.amps.fasta
done
```
Then I blasted the potential amps file from UNIPROT that I have been using.

```bash
#Create database before blasting
makeblastdb -in /lustre/scratch/frcastel/defensins/dataset/amplify/databases/UNIPROT/uniprot_potential_defens_cath_over6aas.fasta -dbtype prot -out uniprotdb

#BLASTP
blastp -query filtered.ampir.bats.maker.predicted.amps.fasta -db uniprotdb -evalue 1e-4 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp
```
Then some looping to change the headers for the genes that were hits from BLAST.

```bash
for i in $(cat /lustre/scratch/frcastel/defensins/dataset/maker/bat_list_amps.txt)
do
echo "BLASTING $i"
#Get gff for the species
cp /lustre/scratch/frcastel/defensins/dataset/maker/"$i"/"$i".maker.output/"$i".gff .

#Change headers
maker_map_ids --prefix "$i" --justify 8 "$i".gff > "$i".contig.map

map_fasta_ids "$i".contig.map filtered.ampir.bats.maker.predicted.amps.fasta
#Get new header names for the species
grep "$i" filtered.ampir.bats.maker.predicted.amps.fasta | sed "s/>//g" > "$i".headers.txt

#Map blast results to the header
map_data_ids "$i".contig.map output.blastp

#Get hits from blast fot this species
grep "$i" output.blastp | sed "s/|/\t/g" | sed "s/_/\t/g"  | cut -f1,4,6 | sed "s/\t/_/g" > "$i".hits.txt

#Concatenate all the results
cat "$i".hits.txt >> all.hits.txt

#Remove some files
rm "$i".gff "$i".contig.map "$i".headers.txt
done
```

The result of this would be a change in the headers of the `filtered.ampir.bats.maker.predicted.amps.fasta` file. I moved the hits results to `/lustre/scratch/frcastel/defensins/dataset/cafe/blast_hits_bats/`

After I annotated the headers with the gene names, I had to get the genes for the outgroup and do some cleaning for < 200 AAs and excluding non traditional aminoacids.


```bash
conda activate seqkit
#Here I am filtering for the outgroups and get sequences using headers
for FILE in $(ls /lustre/scratch/frcastel/defensins/dataset/outgroup_proteins/no_homo_mus/)
do
grep -i "taurus\|btaurus\|equus\|ecaballus\|scrofa\|sscrofa\|familiaris\|cfamiliaris" /lustre/scratch/frcastel/defensins/dataset/outgroup_proteins/no_homo_mus/$FILE >> outgroup_amps_headers.txt
grep -i "taurus\|btaurus\|equus\|ecaballus\|scrofa\|sscrofa\|familiaris\|cfamiliaris" uniprot_potential_defens_cath_over6aas.fasta >> outgroup_amps_headers.txt
sed -i "s/^>//g" outgroup_amps_headers.txt
seqtk subseq /lustre/scratch/frcastel/defensins/dataset/outgroup_proteins/no_homo_mus/$FILE outgroup_amps_headers.txt >> outgroup_amps.fasta
done

seqtk subseq uniprot_potential_defens_cath_over6aas.fasta outgroup_amps_headers.txt >> outgroup_amps.fasta

#Do some cleaning of the file with the mentioned criteria
seqkit fx2tab outgroup_amps.fasta -l -i -H -C BJOUXZ | awk '{ if ($3 <= 200 && $4 < 1 ) {print $1 "\t" $2} }' | seqkit tab2fx | seqkit rmdup -s -o clean_outgroup_amps.fasta -D duplicated_outgroup_amps.txt

#Merge bat amps with outgroups
cat filtered.ampir.bats.maker.predicted.amps.fasta clean_outgroup_amps.fasta > outgroup_bats_amps.fasta

```
After the automatic cleaning I removed some "*" that were by the end of some ecaballus sequences only. I am reducing the complexity of the headers with some code:

```bash
sed -E 's/>([a-z]+)_\S+\|(DEF[^|]+).*/>\1_\2/' outgroup_bats_amps.fasta  > outgroup_bats_amps_edit_headers.fasta
sed -i "s/_Paneth_cellspecific_alphadefensin/DEFA/g ; s/_precursor//g ; s/_alphadefensin_/DEFA_/g; s/_neutrophil_defensin/DEFA/g; s/_PE=.*//g; s/OX=[*]//g" outgroup_bats_amps_edit_headers.fasta
```
I also had to look for the AMPs names in some sequences that were not writen in the headers, mostly sequences of cow, dog and cow corresponding to cathelicidins.

Also, since the changes of the `maker_map_id` get rid of the annotation of the somewhat informative headers, I am merging the results of the blast hits with the headers to make it more informative for downstream analyses.

```bash
# Change the headers in outgroup_bats_amps.fasta
hits_file="all.hits.txt"
fasta_file="outgroup_bats_amps_edit_headers.fasta"

# Loop through each line in the hits file
while read line; do
    # Extract the sequence ID and desired header
    seq_id=$(echo "$line" | cut -d"_" -f1)
    header=$(echo "$line" | cut -d"_" -f2-)

    # Use sed to replace the header in the fasta file
    sed -i "s/^>$seq_id/>${seq_id}_${header}/" $fasta_file
done < $hits_file
```
Final editing to the new headers, just for simplicity

```bash
sed -i "s/_Alpha-defensin//g ; s/_Defensin-.*//g ; s/_Beta-defensin//g ; s/_Cathel.*//g; s/_Defensin//g; s/_Putative//g; s/DB/DEFB/g; s/CAMP/CTHL/g; s/D103A/DEFB103A/g; s/DEF5/DEFA5/g" outgroup_bats_amps_edit_headers.fasta
```
I also did a manual search for the proteins that I manually separated from pdiscolor and acaudifer because since I added a "_dup" to the maker header, then they were not retrieved or recognized when I ran the loop to change the header's names. BLAST showed that both belonged to DEFB109 like proteins so I manually changed the header name to their species and DEFB, adding a "_dup" at the end of it.
I also had to manually look for some sequences that did not have the AMP family in them:
|Accession Number|NCBI ID |
|----------------|--------|
| XP_024838922.1| cathelicidin-1-like   |
|XP_024839052.1|cathelicidin-4|
|XP_020937599.1|protegrin-1|
|NP_777250.1|cathelicidin-1 precursor|
|NP_001123448.1|antibacterial peptide PMAP-23 precursor|
|NP_001123437.1|antibacterial peptide PMAP-36 precursor|
|NP_999615.1|antibacterial protein PR-39 precursor|
|NP_777250.1|cathelicidin-1 precursor|
|NP_777251.1|cathelicidin-2 precursor|
|NP_776426.1|cathelicidin-3 precursor|
|NP_776935.1|cathelicidin-5 precursor|
|NP_777256.1|cathelicidin-7 precursor|
|NP_777257.1|cathelicidin-6 precursor|
|NP_001003359.1|cathelicidin antimicrobial peptide precursor|
|NP_001075338.1|myeloid cathelicidin 2 precursor|
|NP_001075399.1|myeloid cathelicidin 3 precursor|
|NP_001116621.1|protegrin-1 precursor|
|NP_001123438.1|protegrin-2 precursor|
|NP_001116622.1|protegrin-3 precursor|

# Citation (Under revision)

If you use this code or any dataset click in "*Cite this repository*" at the top of this page to get the citation in APA and BibTeX formats.