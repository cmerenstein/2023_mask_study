#!/bin/bash

## load the qiime conda environment (qiime2-2022.2)
conda activate qiime2-2022.2
PATH=$PATH:/home/cmeren/rotation_2020/utilities/ 

## set the directory in which you're working.
cd /media/lorax/users/cmeren/mask/
wd=/media/lorax/users/cmeren/mask/
mkdir qiime_data

## Import the data. Each sample has it's own fastq file, one for forward one for reverse.
qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-path scripts/FIXED_qiime_manifest.txt --output-path qiime_data/MASK.qza --input-format PairedEndFastqManifestPhred33V2

## get the summary vizualizations for the data
qiime demux summarize --i-data qiime_data/MASK.qza --o-visualization qiime_data/dmux.qzv

## export the summary so I can look at it locally
qiime tools export --input-path qiime_data/dmux.qzv --output-path qiime_data/summary_export

## run DADA2 denoising to call ASV
## Needs it to use the version of R that the qiime conda environment wants
PATH=/home/cmeren/miniconda3/envs/qiime2-2022.2/lib/R/bin/:$PATH

## the folder qiime_data/summary_export has the QC information needed to pick the trimming cutoffs.
## The html file "quality-plot.html" shows the quality score cutoffs for both the forward and reverse reads
qiime dada2 denoise-paired --i-demultiplexed-seqs qiime_data/MASK.qza --p-trunc-len-f 197 --p-trunc-len-r 250 --p-trim-left-r 13 --p-trim-left-f 13 --p-max-ee-f 2 --p-max-ee-r 2 --p-n-threads 8 --o-table qiime_data/dada2_table.qza --o-representative-sequences qiime_data/dada2_features.qza --o-denoising-stats qiime_data/dada2_denoising.qza --verbose

## export the results from DADA2
mkdir qiime_data/DADA2/

qiime tools export --input-path qiime_data/dada2_denoising.qza --output-path qiime_data/DADA2/
qiime tools export --input-path qiime_data/dada2_table.qza --output-path qiime_data/DADA2/
qiime tools export --input-path qiime_data/dada2_features.qza --output-path qiime_data/DADA2/

## convert the feature table to plain text
biom convert -i qiime_data/DADA2/feature-table.biom -o qiime_data/DADA2/feature_table.txt --to-tsv

## assign taxonomy
mkdir qiime_data/ASV
mkdir qiime_data/ASV/NB

## Naive Bayes Classifier. I downloaded this from the QIIME2 website
qiime feature-classifier classify-sklearn --i-reads qiime_data/dada2_features.qza --i-classifier qiime_data/silva-138-99-nb-classifier.qza --o-classification qiime_data/ASV/NB/silva_138_99.qza --verbose

qiime tools export --input-path qiime_data/ASV/NB/silva_138_99.qza --output-path qiime_data/ASV/NB/

## there might be some way via QIIME do do this, but this is just a python script I wrote to get
## counts matricies at each taxonomy level.
counts_by_taxon.py -t "$wd"/qiime_data/ASV/NB/taxonomy.tsv -l 7 -f "$wd"/qiime_data/DADA2/feature_table.txt -o "$wd"/qiime_data/ASV/NB/

## generate phylogenetic tree
mkdir qiime_data/tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences qiime_data/dada2_features.qza  --o-alignment qiime_data/tree/aligned-rep-seqs.qza --o-masked-alignment qiime_data/tree/masked-aligned-rep-seqs.qza --o-tree qiime_data/tree/unrooted-tree.qza --o-rooted-tree qiime_data/tree/rooted-tree.qza

## Export tree to do UniFrac in R
qiime tools export --input-path qiime_data/tree/rooted-tree.qza --output-path qiime_data


## Use this phylogenetic tree to get the UniFrac and other metrics
#qiime diversity core-metrics-phylogenetic  --i-phylogeny qiime_data/tree/rooted-tree.qza --i-table qiime_data/dada2_table.qza --p-sampling-depth 1000  --output-dir qiime_data/diversity

## For one analysis we are considering using 97% OTUs, because ASVs maybe are splitting up taxa
## that shouldn't be split up (at least, splitting them doesn't seem clinically important).
mkdir qiime_data/OTU
qiime vsearch cluster-features-de-novo --i-table qiime_data/dada2_table.qza --i-sequences qiime_data/dada2_features.qza  --p-perc-identity 0.97 --o-clustered-table qiime_data/OTU/table-dn-97.qza --o-clustered-sequences qiime_data/OTU/rep-seqs-dn-97.qza

qiime tools export --input-path qiime_data/OTU/table-dn-97.qza --output-path qiime_data/OTU/
biom convert -i qiime_data/OTU/feature-table.biom -o qiime_data/OTU/feature_table.txt --to-tsv

## classify
qiime feature-classifier classify-sklearn --i-reads qiime_data/OTU/rep-seqs-dn-97.qza --i-classifier /home/cmeren/program_files/qiime/classifiers/silva-138-99-nb-classifier.qza --o-classification qiime_data/OTU/OTU_97_silva_138_99.qza --verbose

qiime tools export --input-path qiime_data/OTU/OTU_97_silva_138_99.qza --output-path qiime_data/OTU/





