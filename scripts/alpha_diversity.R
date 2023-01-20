library(ape)
library(vegan)
library(tidyverse)

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
asv_no_contam = asv[, colSums(asv) > 1000]

## ------------ calculate diversity metrics --------------

richness = colSums( asv_no_contam> 0 )
simpson = diversity(t(asv_no_contam), "simpson")
shannon = diversity(t(asv_no_contam), "shannon")
invsimpson = diversity(t(asv_no_contam), "invsimpson")

diversity_all = gather( data.frame(sample = names(simpson), richness, simpson, shannon, invsimpson),
                            "metric", "diversity", -sample)

## ------------- plot plot plot -----------------------------
## ---- MASK -----

masks = meta[meta$type == "MASK", ]
masks = merge(diversity_all, masks, by = "sample")

pdf("../figures/mask_alpha_diversity.pdf")
ggplot(masks, aes( x = as.factor(n_days_this_mask), y = diversity) ) + 
    geom_boxplot() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    facet_wrap(~metric, scale = "free") + 
    xlab("Days of this Mask") +
    ylab("Diversity")

dev.off()

## ---- THROAT ------------

throat = meta[meta$type == "THROAT", ]
throat = merge(diversity_all, throat, by = "sample")

## ~~ Compare each subject to their own baseline ~~
V1_throat = meta %>% filter(type == "THROAT" & visit == "V1") %>%
                inner_join( diversity_all, by = "sample")
V3_throat = meta %>% filter(type == "THROAT" & visit == "V3") %>%
                inner_join( diversity_all, by = "sample")

first_last = merge(V1_throat, V3_throat, by = "subject") %>% 
                filter(metric.x == metric.y)

first_last$change = (first_last$diversity.y - first_last$diversity.x)

for (metric in unique(first_last$metric.x)){
    metric_change = first_last[ first_last$metric.x == metric, ]
    print(metric)
    print(wilcox.test( metric_change$change ~ metric_change$group.y ))
}

pdf("../figures/throat_alpha_diversity.pdf")
ggplot( filter(first_last, metric.x %in% c("richness", "simpson")),
  aes( x = as.factor(n_days_this_mask.y), y = change) ) + 
    geom_boxplot() + 
    geom_point() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    facet_wrap(~metric.x, scale = "free") + 
    xlab("Days of this Mask") +
    ylab("Change in Diversity") +
    ggtitle("Throat samples")

dev.off()

## ---- SKIN --------

skin = meta[meta$type == "SKIN", ]
skin = merge(diversity_all, skin, by = "sample")

V1_skin = meta %>% filter(type == "SKIN" & visit == "V1") %>%
                inner_join( diversity_all, by = "sample")
V3_skin = meta %>% filter(type == "SKIN" & visit == "V3") %>%
                inner_join( diversity_all, by = "sample")

skin_first_last = merge(V1_skin, V3_skin, by = "subject") %>% 
                    filter(metric.x == metric.y)

skin_first_last$change = (skin_first_last$diversity.y - skin_first_last$diversity.x)

for (metric in unique(skin_first_last$metric.x)){
    metric_change = skin_first_last[ skin_first_last$metric.x == metric, ]
    print(metric)
    print(wilcox.test( metric_change$change ~ metric_change$group.y ))
}

pdf("../figures/skin_alpha_diversity.pdf")

ggplot(filter( skin_first_last, metric.x %in% c("richness", "simpson")), 
  aes( x = as.factor(n_days_this_mask.y), y = change) ) + 
    geom_boxplot() + 
    geom_point() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    facet_wrap(~metric.x, scale = "free") + 
    xlab("Days of this Mask") +
    ylab("Change in Diversity") +
    ggtitle("Skin samples")

dev.off()














