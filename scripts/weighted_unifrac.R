library(ape)
library(phyloseq)
library(vegan)
library(tidyverse)

## ------- Get phylogenetic tree from QIIME ---------------
tree = read_tree("../data/tree.nwk")

meta = read.csv("../data/all_meta.csv")
rownames(meta) = meta$sample

## load and filter ASVs
asv_df = read.table("../data/DADA2/feature_table_clean.txt", sep = "\t", header = T)
asv = asv_df[,2:ncol(asv_df)]
rownames(asv) = asv_df$OTU.ID

## 20 count cutoff for individual ASVs, 1000 reads minimum for samples
asv[asv < 20] = 0
asv = asv[,colSums(asv) > 1000]

## phyloseq requires its own otu_table object
asv_table = otu_table(asv, taxa_are_rows = T)
phylo = phyloseq(asv_table , tree)

## rarefied weighted UniFrac
set.seed(02130)
rare = rrarefy(t(asv), 1000)
rare_phylo = phyloseq( otu_table(rare, taxa_are_rows = F), tree)
rare_weighted = UniFrac(rare_phylo, weighted = T)

saveRDS(rare_weighted, "../data/pcoa/rare_weighted_unifrac.rds") 

## ------------------- PERMANOVA by sample type -------------
unifrac_mat = as.matrix(rare_weighted)
adonis2( unifrac_mat ~ type, data = meta[rownames(unifrac_mat),], permutations = 9999)

# -------------------- PCOA plot of weighted UniFrac --------------
set.seed(19002)
ordinates = pcoa(rare_weighted)
vectors = ordinates$vectors

pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    type = meta[rownames(vectors), "type"],
                    days = as.character(meta[rownames(vectors), "n_days_this_mask"]))
pcoa_df$days[is.na(pcoa_df$days)] = "ctrl"


## ---- remove positive controls and improve label readability -----
pcoa_df = pcoa_df[ !(rownames(pcoa_df) %in%  c("POS.MASK", "mockdna.2", "mockdna.3")),]

## Change the names of type and make them look nice
pcoa_df$type = ifelse( pcoa_df$days != "ctrl", pcoa_df$type,
                ifelse( grepl("SKIN", pcoa_df$sample), "Unused skin swab",
                  ifelse( grepl("throat", pcoa_df$sample), "Unused throat swab",
                    ifelse( grepl("MASKBLANK", pcoa_df$sample), "Unused mask", "Extraction blank"))))

pcoa_df$type = sapply(pcoa_df$type, function(type){
                    paste(toupper( substr(type, 1,1)), tolower( substr(type, 2, nchar(type))), sep = "" )})

pcoa_df$type = factor(pcoa_df$type, levels = c("Mask", "Throat", "Skin",
                               "Unused throat swab", "Unused skin swab", "Unused mask", "Extraction blank") )

pdf("../figures/sample_type_rarefied_pcoa.pdf")
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = type,  fill = days) ) +
    theme_classic() +
    geom_point(shape = 21, size = 4, stroke = 2) +
    scale_fill_manual(values = c("white", "grey", "white")) +
    scale_color_manual(values = c("skyblue", "salmon", "palegreen4", "plum2",
                                "lightpink3", "lightpink1", "darkorchid4"))

dev.off()


### ------------------ Distance between sample types -------------------------
unifrac = as.matrix(rare_weighted)
unifrac_df = data.frame(sample.x = rownames(unifrac), unifrac) %>% 
                    gather("sample.y", "distance", -sample.x)

unifrac_df$comparison = paste( meta[unifrac_df$sample.x, "type"],
                               meta[unifrac_df$sample.y, "type"], sep = "-")

## only care about mask-skin, mask-throat, mask-blank
filtered = filter(unifrac_df, grepl("MASK-", comparison))  

TukeyHSD(aov( distance ~ comparison, data = filtered)) 







