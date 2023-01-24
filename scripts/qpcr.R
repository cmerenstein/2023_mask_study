library(tidyverse)

## load and clean
qpcr = read.csv("../data/qPCR_ct_values.csv")

meta = read.csv("../data/all_meta.csv")

## remove failed sample and make numeric
qpcr_all = merge(qpcr, meta, by.x = "SampleID", by.y = "sample", all.x = T, all.y = T)
qpcr_all$Ct_ave = as.numeric(qpcr_all$Ct_ave)
qpcr_all = qpcr_all[ !(is.na(qpcr_all$Ct_ave)),]

## set group for controls
qpcr_all[qpcr_all$SampleID == "POS.MASK", "group"] = "Pos.control"
qpcr_all[qpcr_all$SampleID == "mockdna.3", "group"] = "Pos.control"
qpcr_all[ is.na(qpcr_all$group), "group"] = "Neg.control"

## set control days 0
qpcr_all[ grepl("MASKBLANK", qpcr_all$SampleID), "n_days_this_mask"] = 0
qpcr_masks = qpcr_all[ !(is.na(qpcr_all$n_days_this_mask)) ,]

## -------------- mask CT by visit number -----------------------------
pdf("../figures/ct_by_visit.pdf")
ggplot(qpcr_masks, aes(x = visit, y = Ct_ave, fill = group)) + 
    theme_classic() +  
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    theme(text = element_text(size = 20)) #+ 
    #scale_fill_manual(values = c("grey50", "grey90"))

ggplot(qpcr_masks, aes(x = visit, y = Ct_ave, color = group)) + 
    theme_classic() + 
    geom_point() + 
    geom_line(aes(group = subject))
dev.off()

TukeyHSD(aov(Ct_ave~visit, data = qpcr_masks))

## -------------- mask CT by number of days worn ----------------------

pdf("../figures/ct_by_days.pdf")
ggplot(qpcr_masks, aes(x = n_days_this_mask, y = Ct_ave)) + 
    theme_classic() +  
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    geom_point() + 
    theme(text = element_text(size = 20)) +
    xlab("Number of days wearing this mask")#+ 
    #scale_fill_manual(values = c("grey50", "grey90"))
dev.off()

group_by(qpcr_masks, as.character(n_days_this_mask)) %>% summarize( mean(Ct_ave))

wilcox.test( Ct_ave ~ as.character(n_days_this_mask), data = qpcr_masks[qpcr_masks$n_days_this_mask != 0,])
wilcox.test( Ct_ave ~ as.character(n_days_this_mask), data = qpcr_masks[qpcr_masks$n_days_this_mask != 7,])
wilcox.test( Ct_ave ~ as.character(n_days_this_mask), data = qpcr_masks[qpcr_masks$n_days_this_mask != 1,])

## correlation within subject
V1 = qpcr_masks[qpcr_masks$visit == "V1", c("subject", "Ct_ave", "group")]
V3 = qpcr_masks[qpcr_masks$visit == "V3", c("subject", "Ct_ave")]
first_last = merge(V1, V3, by = "subject") %>% filter(!(is.na(subject)))

## --------------- Correlation between V1 and V3 masks ---------------

pdf("../figures/ct_correlation_v1_v3.pdf")
ggplot(first_last, aes(x = Ct_ave.x, y = Ct_ave.y)) + 
    theme_classic() + 
    geom_point() + 
    stat_smooth(method = "lm") +
    ylab("Visit 3 CT") + 
    xlab("Visit 1 CT") + 
    theme(text = element_text(size = 20)) +
    coord_fixed()+
    ylim(c( min( c(first_last$Ct_ave.x, first_last$Ct_ave.y)), 
            max( c(first_last$Ct_ave.x, first_last$Ct_ave)))) + 
    xlim(c( min( c(first_last$Ct_ave.x, first_last$Ct_ave.y)), 
            max( c(first_last$Ct_ave.x, first_last$Ct_ave)))) 
dev.off()

summary(lm(Ct_ave.x~Ct_ave.y, data = first_last))






