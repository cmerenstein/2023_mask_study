library(vegan)
library(tidyverse)
library(ape)


meta = read.csv("../data/all_meta.csv")
rownames(meta) = meta$sample

unweighted = as.matrix(readRDS( "../data/pcoa/no_contam_unweighted_unifrac.rds"))

## filter
meta = meta[rownames(unweighted),]

##
V1_skin = meta %>% filter(type == "SKIN" & visit == "V1")
V2_skin = meta %>% filter(type == "SKIN" & visit == "V2")
V3_skin = meta %>% filter(type == "SKIN" & visit == "V3")

first_last = merge(V1_skin, V3_skin, by = "subject")
first_last$distance = apply(first_last, 1, function(row){ 
                             first = row["sample.x"]
                             last = row["sample.y"]
                             return( unweighted[first, last] ) })

pdf("../figures/unweighted_distance_to_first_skin.pdf", height = 5, width = 5)
ggplot(first_last, aes(x = group.x, y = distance)) +
    geom_boxplot() +
    theme_classic() + 
    theme(text = element_text(size = 20))
dev.off()

t.test(first_last$distance~group.x, first_last)


### ----------------- PCOA ----------------------------
set.seed(19002)
unifrac_skin = unweighted[grepl("MASK.+SKIN", rownames(unweighted)),
                               grepl("MASK.+SKIN", colnames(unweighted))]
ordinates = pcoa(unifrac_skin)
vectors = ordinates$vectors

## plot pcoa
pcoa_df = data.frame(sample = rownames(vectors), Axis.1 = vectors[,1],
                    Axis.2 = vectors[,2], Axis.3 = vectors[,3], Axis.4 = vectors[,4],
                    visit = meta[rownames(vectors), "visit"],
                    subject = meta[rownames(vectors), "subject"])
pcoa_df$visit = factor(pcoa_df$visit, levels = c("V1", "V2", "V3"))
pcoa_df = arrange(pcoa_df, visit)


## by number of days 
pdf("../figures/unweighted_skin_pcoa_by_days.pdf", height = 4.5, width = 5)
df_day = left_join(pcoa_df, meta, by = "sample") %>%
    mutate( days_of_mask = ifelse(visit.x == "V3" & group == "reuse group", "7 days", "1 day"))

ggplot(df_day, aes(x = Axis.1, y = Axis.2, color = days_of_mask)) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        geom_point(size = 4, shape = 21, stroke = 2) +
        scale_color_manual(values = c("lightgreen", "forestgreen")) +
        stat_ellipse() +
        theme(legend.position = "none")
dev.off()

adonis2( unifrac_skin[df_day$sample, df_day$sample] ~ df_day$days_of_mask)




