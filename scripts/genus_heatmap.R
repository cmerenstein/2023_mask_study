library(pheatmap)
library(tidyverse)

## read and filter
genus = read.csv("../data/qiime_data/NB/counts_taxa_level_5.csv")

totals = rowSums(genus[,2:ncol(genus)])

filtered = genus[totals >= 100,]
rownames(filtered) = filtered$sample

## make metadata from names
meta= read.csv("../data/all_meta.csv")
rownames(meta) = meta$sample
meta = meta[rownames(filtered),]

## get matrix in percent format
filtered_mat = filtered[,2:ncol(filtered)]
percent = filtered_mat / rowSums(filtered_mat)

## plot common taxa heatmap
threshold = 0.05
common = percent[,colSums(percent > threshold) > 1] ## over threshold percent in at least one sample

ann_colors = list( type = c(CTRL = "grey80", MASK = "skyblue", SKIN = "goldenrod", THROAT = "chocolate"), 
                    n_days_this_mask = c(na = "grey80", "1" = "palegreen", "7" = "forestgreen"))

## change NA to "na" for plotting
meta[ is.na(meta$n_days_this_mask), "n_days_this_mask"] = "na"

pdf("../figures/genus_heatmap.pdf", width = 10, height = 6)
pheatmap(t(common), annotation_col = meta[,c("type", "n_days_this_mask")] , 
                    labels_col = character(nrow(common)),
                    annotation_colors = ann_colors, treeheight_row = 0)
dev.off()

##  ---------- plot heatmap by sample type ----------------

controls = meta[ meta$type == "CTRL", "sample"]
throat = meta[ meta$type == "THROAT", "sample"]
mask = meta[ meta$type == "MASK", "sample"]
skin = meta[ meta$type == "SKIN", "sample"]

## we still want clustering within each sample type
get_order <- function(samples) { return(hclust(dist(common[ rownames(common) %in% samples,]))$order)}

ord_controls = controls[get_order(controls)]
ord_throat = throat[get_order(throat)]
ord_mask = mask[get_order(mask)]
ord_skin = skin[get_order(skin)]

## reorder common
reordered = c(ord_controls, ord_mask, ord_skin, ord_throat)

pdf("../figures/genus_heatmap_by_type.pdf", width = 10, height = 6)
pheatmap(t(common[reordered,]), annotation_col = meta[reordered,c("type", "n_days_this_mask")], 
            labels_col = character(nrow(common)), cluster_cols = F,
            annotation_colors = ann_colors, treeheight_row = 0)
dev.off()















