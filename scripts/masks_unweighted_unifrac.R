library(vegan)
library(tidyverse)
library(ape)

meta = read.csv("../data/all_meta.csv")
rownames(meta) = meta$sample

unweighted = as.matrix(readRDS( "../data/pcoa/no_contam_unweighted_unifrac.rds"))

## samples without contamination is a smaller list, because some got removed
meta = meta[rownames(unweighted),]

##
V1_mask = meta %>% filter(type == "MASK" & visit == "V1")
V3_mask = meta %>% filter(type == "MASK" & visit == "V3")

first_last = merge(V1_mask, V3_mask, by = "subject")
first_last$distance = apply(first_last, 1, function(row){ 
                             first = row["sample.x"]
                             last = row["sample.y"]
                             return( unweighted[first, last] ) })

pdf("../figures/unweighted_distance_to_first_mask.pdf", height = 5, width = 5)
ggplot(first_last, aes(x = group.x, y = distance)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Weighted UniFrac distance, Day 7 vs Day 0") + 
    xlab("Mask Study Group")
dev.off()

t.test(first_last$distance~group.x, first_last)


## -------------- PCOA plot by number of days worn --------
set.seed(19002)
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
pdf("../figures/unweighted_mask_pcoa_by_days.pdf")
df_day = left_join(pcoa_df, meta, by = "sample") %>%
    mutate( days_of_mask = ifelse(visit.x == "V3" & group == "reuse group", "7 days", "1 day"))

ggplot(df_day, aes(x = Axis.1, y = Axis.2, color = days_of_mask)) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        geom_point(size = 4, stroke = 2, shape = 21) +
        stat_ellipse() +
        scale_color_manual(values = c("skyblue", "navy"))
dev.off()


adonis2( unifrac_mask[df_day$sample, df_day$sample] ~ df_day$days_of_mask)




