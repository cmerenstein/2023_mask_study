library(ape)
library(phyloseq)
library(vegan)
library(tidyverse)

tree = read_tree("../data/tree.nwk")

meta = read.csv("../data/all_meta.csv")
rownames(meta) = meta$sample

## load and filter ASVs
asv_df = read.table("../data/DADA2/feature_table_clean.txt", sep = "\t", header = T)
asv = asv_df[,2:ncol(asv_df)]
rownames(asv) = asv_df$OTU.ID 

## Filter low abundance
asv[asv < 20] = 0
asv = asv[ rowSums(asv) > 100 ,colSums(asv) > 1000]

## convert to percent
percent = t(asv)
percent = percent / rowSums(percent)

## ----- find significant associations with controls -----
controls = meta[rownames(percent), "type"] == "CTRL"
control_pval_list = apply(percent, 2, function(col){
                    test = wilcox.test(col ~ controls)
                    diff = mean( col[controls] ) - mean( col[!(controls)] )
                    return( data.frame( W = test$statistic , p = test$p.value, diff ) )
                })
control_pval = do.call("rbind", control_pval_list)
contaminants = control_pval[control_pval$p < 0.01 & control_pval$diff > 0,]

asv_no_contam = asv[ !(rownames(asv) %in% rownames(contaminants)), ]
asv_no_contam = asv_no_contam[ , colSums(asv_no_contam) > 1000 ] 

## --------

## phyloseq requires its own otu_table object
asv_table = otu_table(asv_no_contam, taxa_are_rows = T)
new_tree <- ape::multi2di(tree) ## need to fix the tree after pruning contaminats
phylo = phyloseq(asv_table , new_tree) 

set.seed(19143)
unweighted = UniFrac(phylo, weighted = F)

## rarefied
set.seed(02130)
rare = rrarefy(t(asv_no_contam), 1000)
rare_phylo = phyloseq( otu_table(rare, taxa_are_rows = F), new_tree)
rare_weighted = UniFrac(rare_phylo, weighted = T)

saveRDS(unweighted, "../data/pcoa/no_contam_unweighted_unifrac.rds")
saveRDS(rare_weighted, "../data/pcoa/no_contam_rare_weighted_unifrac.rds") 


## -------------------
### unweighted PCOA of mask samples colored by days worn
set.seed(19002)
unweighted = as.matrix(unweighted)
unifrac_mask = unweighted[grepl("MASK.+MASK", rownames(unweighted)),
                               grepl("MASK.+MASK", colnames(unweighted))]
ordinates = pcoa(unifrac_mask)
vectors = ordinates$vectors

## make df for plotting pcoa
pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    visit = meta[rownames(vectors), "visit"],
                    subject = meta[rownames(vectors), "subject"])
pcoa_df$visit = factor(pcoa_df$visit, levels = c("V1", "V2", "V3"))
pcoa_df = arrange(pcoa_df, visit)

## plot by number of days wearing this mask
pdf("figures/unweighted/no_contam_mask_pcoa_by_days.pdf")
df_day = left_join(pcoa_df, meta, by = "sample") %>%
    mutate( days_of_mask = ifelse(visit.x == "V3" & group == "reuse group", "7 days", "1 day"))

ggplot(df_day, aes(x = Axis.1, y = Axis.2, color = days_of_mask)) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        geom_point(size = 2) +
        stat_ellipse()
dev.off()

adonis2( unifrac_mask[df_day$sample, df_day$sample] ~ df_day$days_of_mask)

## --------------------------------------------------------------------
### weighted PCOA of mask samples colored by days worn
set.seed(19002)
rare_weighted = as.matrix(rare_weighted)
unifrac_mask = rare_weighted[grepl("MASK.+MASK", rownames(rare_weighted)),
                               grepl("MASK.+MASK", colnames(rare_weighted))]
ordinates = pcoa(unifrac_mask)
vectors = ordinates$vectors

## make df for plotting pcoa
pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    visit = meta[rownames(vectors), "visit"],
                    subject = meta[rownames(vectors), "subject"])
pcoa_df$visit = factor(pcoa_df$visit, levels = c("V1", "V2", "V3"))
pcoa_df = arrange(pcoa_df, visit)

## plot by number of days wearing this mask
pdf("figures/no_contam_mask_pcoa_by_days.pdf")
df_day = left_join(pcoa_df, meta, by = "sample") %>%
    mutate( days_of_mask = ifelse(visit.x == "V3" & group == "reuse group", "7 days", "1 day"))

ggplot(df_day, aes(x = Axis.1, y = Axis.2, color = days_of_mask)) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        geom_point(size = 2) +
        stat_ellipse()
dev.off()

adonis2( unifrac_mask[df_day$sample, df_day$sample] ~ df_day$days_of_mask)

