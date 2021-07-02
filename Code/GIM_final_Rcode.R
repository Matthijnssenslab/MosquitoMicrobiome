############### *** load library ################
necessary_packages <- c('ggpubr', 'PMCMR', 'corrplot','Hmisc', 'ape','ade4','phyloseq','ggplot2','plyr',
                        'made4','mclust','mixtools','data.table','reshape2','biomformat','DESeq2','gridExtra',
                        'grid','plyr','vegan','dplyr','magrittr','scales','reshape2','randomForest', 'knitr', 'ggthemes',
                        'gplots', 'Heatplus', 'vegan', 'RColorBrewer', 'ComplexHeatmap', 'circlize', 'viridis', 'ggsci',
                        'akmedoids', 'wesanderson', 'tidyr','VennDiagram','viridis')
packages <- lapply(necessary_packages, library, character.only = TRUE)
library(viridis)
shap=c(15,16,17,7,8,10,11,3,4,13,18,12)

############### *** Load Eukaryotic Virome Phyloseq Project #######################
setwd("/Users/shichenyan/Desktop/mosq_Guad_infection/")
GP.M1 <- readRDS("/Users/shichenyan/Desktop/mosq_Guad_infection/RDS/GIM_virome_raw.RDS")

# otu <- as.data.frame(otu_table(GP.M1_AA))
# tax <- as.data.frame(tax_table(GP.M1_AA))
# df <- merge(tax, otu, by="row.names")
# write.table(df, file = "/Users/shichenyan/Desktop/mosq_Guad_infection/AA_tax_abundance.txt", quote = FALSE, sep = "\t", row.names = F)
#################### *** 1 EuV Aedes aegypti *** ##############################
#################### * color ####################
col = plasma(9,alpha = 1, begin=0.9, end = 0, direction = 1)
#################### * 1-1) subset Aedes aegypti Eukaryotic Virome #################### 
GP.M1_AA <- subset_samples(GP.M1, mosq_species  == "Ae_aegypti")
GP.M1_AA
GP.M1_AA <- subset_samples(GP.M1_AA, group.plaque.qPCR.1000. != "")
GP.M1_AA
GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Guadeloupe Culex tymo-like virus")

GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Chuvirus Mos8Chu0")
GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Kaiowa virus")
GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Cumbaru virus")
GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Guato virus")
GP.M1_AA <- subset_taxa(GP.M1_AA, species != "Trichoplusia ni TED virus")

GP.M1_AA
tax <- as.data.frame(tax_table(GP.M1_AA))
tax_table(GP.M1_AA)

sum(sample_sums(GP.M1_AA))
sample_variables(GP.M1_AA)
sample_data(GP.M1_AA)$Treatment <- factor(sample_data(GP.M1_AA)$Treatment, levels = c("AA-sucrose", "AA-blood", "AA-ZIKV"))
sample_data(GP.M1_AA)$Info1 <- factor(sample_data(GP.M1_AA)$Info1, levels = c("AA-7dpi-sucrose", "AA-7dpi-blood", "AA-7dpi-ZIKV",   
                                                                              "AA-21dpi-sucrose","AA-21dpi-blood","AA-21dpi-ZIKV"))
sample_data(GP.M1_AA)$group.plaque.qPCR.1000. <- factor(sample_data(GP.M1_AA)$group.plaque.qPCR.1000., 
                                                        levels = c("AA-7dpi-ZIKV(-/+)","AA-7dpi-ZIKV(-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                   "AA-21dpi-ZIKV(+/+)","AA-21dpi-ZIKV(-/+)","AA-21dpi-ZIKV(-/-)","AA-21dpi-blood","AA-21dpi-sucrose"))
sample_data(GP.M1_AA)$group.plaque.qPCR.0. <- factor(sample_data(GP.M1_AA)$group.plaque.qPCR.0., 
                                                     levels = c("AA-7dpi-ZIKV(-/+)","AA-7dpi-ZIKV(-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                "AA-21dpi-ZIKV(+/+)","AA-21dpi-ZIKV(-/+)","AA-21dpi-ZIKV(-/-)","AA-21dpi-blood","AA-21dpi-sucrose"))
sample_data(GP.M1_AA)$group.plaque.qPCR.1000.body.1000. <- factor(sample_data(GP.M1_AA)$group.plaque.qPCR.1000.body.1000.,
                                                                  levels = c("AA-7dpi-ZIKV(-/+/+)","AA-7dpi-ZIKV(-/-/+)","AA-7dpi-ZIKV(-/-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                             "AA-21dpi-ZIKV(+/+/+)","AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(-/-/+)","AA-21dpi-ZIKV(-/-/-)",
                                                                             "AA-21dpi-blood","AA-21dpi-sucrose"))
sample_data(GP.M1_AA)$group.plaque.qPCR.0.body.0. <- factor(sample_data(GP.M1_AA)$group.plaque.qPCR.0.body.0, 
                                                            levels = c("AA-7dpi-sucrose", "AA-7dpi-blood","AA-7dpi-ZIKV(-/-/-)", "AA-7dpi-ZIKV(-/-/+)", "AA-7dpi-ZIKV(-/+/+)",
                                                                       "AA-21dpi-sucrose", "AA-21dpi-blood", "AA-21dpi-ZIKV(-/-/-)", "AA-21dpi-ZIKV(-/-/+)","AA-21dpi-ZIKV(-/+/-)",
                                                                       "AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(+/+/+)"))
sample_data(GP.M1_AA)$plaque_assay <- factor(sample_data(GP.M1_AA)$plaque_assay,
                                             levels = c("AA-7dpi-sucrose", "AA-7dpi-blood", "AA-7dpi-ZIKV(-)",
                                                        "AA-21dpi-sucrose", "AA-21dpi-blood", "AA-21dpi-ZIKV(-)","AA-21dpi-ZIKV(+)"))

sample_data(GP.M1_AA)$dpi <- factor(sample_data(GP.M1_AA)$dpi, levels = c("7dpi", "21dpi"))
levels(sample_data(GP.M1_AA)$Treatment)
levels(sample_data(GP.M1_AA)$Info1)
levels(sample_data(GP.M1_AA)$group.plaque.qPCR.1000.)
levels(sample_data(GP.M1_AA)$group.plaque.qPCR.0.)
levels(sample_data(GP.M1_AA)$group.plaque.qPCR.1000.body.1000.)
levels(sample_data(GP.M1_AA)$group.plaque.qPCR.0.body.0.)
levels(sample_data(GP.M1_AA)$plaque_assay)

#################### * 1-2) merge on species level  ########
GP.M1_species <- GP.M1_AA %>%
  tax_glom(taxrank = "species")
GP.M1_species

# sample
GP.M1_samselect = prune_samples(sample_sums(GP.M1_species) >0, GP.M1_species)
GP.M1_samselect                                

# species
GP.M1_species_select_0 = prune_taxa(taxa_sums(GP.M1_samselect) > 0, GP.M1_samselect)
GP.M1_species_select_0
summary(taxa_sums(GP.M1_species_select_0))
summary(sample_sums(GP.M1_species_select_0))

GP.M1_species_select = prune_taxa(taxa_sums(GP.M1_samselect) > 500, GP.M1_samselect)
GP.M1_species_select
summary(taxa_sums(GP.M1_species_select))
summary(sample_sums(GP.M1_species_select))
#################### * 1-3) HEATMAP on viral species level ############
#### get table - reads for each viral species ####
OTU_species <- merge(otu_table(GP.M1_species_select), tax_table(GP.M1_species_select), by=0)
rownames(OTU_species) <- OTU_species$species
rownames(OTU_species)
colnames(OTU_species)

#### order viral species by familiy #### 
OTU_species$family <- as.character(OTU_species$family)
OTU_species$family
OTU_species$family[is.na(OTU_species$family)] <- "unclassified"
OTU_species_m<- melt(OTU_species)
OTU_species_m_agg <- aggregate(value~family+species, FUN=sum, data=OTU_species_m)
OTU_species_m_agg1 <- aggregate(value~family, FUN=sum, data=OTU_species_m)
family_ord <- OTU_species_m_agg1[order(-OTU_species_m_agg1$value),]$family
family_ord <- c(family_ord[-3],"unclassified")
family_ord
OTU_species_m_agg$index <- 1
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Totiviridae")] <- 2
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Flaviviridae")] <- 3
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Orthomyxoviridae")] <- 4
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Xinmoviridae")] <- 5
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Metaviridae")] <- 6
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Virgaviridae")] <- 7
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Rhabdoviridae")] <- 8
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Reoviridae")] <- 9
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "unclassified")] <- 10
OTU_species_m_agg
OTU_species_m_agg <- OTU_species_m_agg[order(OTU_species_m_agg$index, -OTU_species_m_agg$value),]
OTU_species_m_agg
species_ord <- as.character(OTU_species_m_agg$species)
species_ord
which(species_ord == "Zika virus")
species_ord_1 <- c("Zika virus", species_ord[-5])
species_ord_1

#### transform the table - reads for each viral species ##### 
rm <- c("Row.names", "superkingdom", "phylum", "class", "order", "family", "genus", "species")
OTU_species <- OTU_species[,!colnames(OTU_species) %in% rm]

#log
OTU_species_trans <- log10(OTU_species)
is.na(OTU_species_trans) <- sapply(OTU_species_trans, is.infinite)
OTU_species_trans[is.na(OTU_species_trans)]<-0
df.original <- OTU_species
df <- OTU_species_trans

#### use meta data for top annotation #####
data.frame(sample_data(OTU_species))
meta = read.csv("GIM_metadata_new.csv", header=TRUE, row.names=1, sep=",")
meta_s = meta[!meta$mosq_species == "",]

meta_AA <- meta[colnames(df),]

#### order sample by ZIKV reads within each group & column annotation  #### 
meta_AAz <- merge(meta_AA, as.data.frame(t(df)),by=0)
meta_AAz
meta_AAz$Treatment <- as.character(meta_AAz$Treatment)
meta_AAz$Treatment <- factor(meta_AAz$Treatment, levels = c("AA-ZIKV", "AA-blood", "AA-sucrose"))
meta_AAz$Info1 <- as.character(meta_AAz$Info1)
meta_AAz$Info1 <- factor(meta_AAz$Info1, levels = c("AA-7dpi-ZIKV", "AA-7dpi-blood", "AA-7dpi-sucrose",
                                                    "AA-21dpi-ZIKV", "AA-21dpi-blood", "AA-21dpi-sucrose"))
meta_AAz$group.plaque.qPCR.1000. <- as.character(meta_AAz$group.plaque.qPCR.1000.)
meta_AAz$group.plaque.qPCR.1000. <- factor(meta_AAz$group.plaque.qPCR.1000.,levels = c("AA-7dpi-ZIKV(-/+)","AA-7dpi-ZIKV(-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                                       "AA-21dpi-ZIKV(+/+)","AA-21dpi-ZIKV(-/+)","AA-21dpi-ZIKV(-/-)","AA-21dpi-blood","AA-21dpi-sucrose"))
meta_AAz$group.plaque.qPCR.0. <- as.character(meta_AAz$group.plaque.qPCR.0.)
meta_AAz$group.plaque.qPCR.0. <- factor(meta_AAz$group.plaque.qPCR.0.,levels = c("AA-7dpi-ZIKV(-/+)","AA-7dpi-ZIKV(-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                                 "AA-21dpi-ZIKV(+/+)","AA-21dpi-ZIKV(-/+)","AA-21dpi-ZIKV(-/-)","AA-21dpi-blood","AA-21dpi-sucrose"))
meta_AAz$group.plaque.qPCR.1000.body.1000. <- as.character(meta_AAz$group.plaque.qPCR.1000.body.1000.)
meta_AAz$group.plaque.qPCR.1000.body.1000. <- factor(meta_AAz$group.plaque.qPCR.1000.body.1000.,
                                                     levels = c("AA-7dpi-ZIKV(-/+/+)","AA-7dpi-ZIKV(-/-/+)","AA-7dpi-ZIKV(-/-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                "AA-21dpi-ZIKV(+/+/+)","AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(-/-/+)","AA-21dpi-ZIKV(-/-/-)",
                                                                "AA-21dpi-blood","AA-21dpi-sucrose"))
meta_AAz$group.plaque.qPCR.0.body.0. <- as.character(meta_AAz$group.plaque.qPCR.0.body.0.)
meta_AAz$group.plaque.qPCR.0.body.0. <- factor(meta_AAz$group.plaque.qPCR.0.body.0.,
                                               levels = c("AA-7dpi-ZIKV(-/+/+)","AA-7dpi-ZIKV(-/-/+)","AA-7dpi-ZIKV(-/-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                          "AA-21dpi-ZIKV(+/+/+)","AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(-/+/-)","AA-21dpi-ZIKV(-/-/+)",
                                                          "AA-21dpi-ZIKV(-/-/-)","AA-21dpi-blood","AA-21dpi-sucrose"))
meta_AAz$plaque_assay <- as.character(meta_AAz$plaque_assay)
meta_AAz$plaque_assay <- factor(meta_AAz$plaque_assay, levels = c("AA-7dpi-ZIKV(-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                  "AA-21dpi-ZIKV(+)","AA-21dpi-ZIKV(-)","AA-21dpi-blood","AA-21dpi-sucrose"))

levels(meta_AAz$Treatment)
levels(meta_AAz$Info1)
levels(meta_AAz$group.plaque.qPCR.1000.)
levels(meta_AAz$group.plaque.qPCR.0.)
levels(meta_AAz$plaque_assay)

meta_AAz$qPCR_head._ZIKV_WNV
meta_AAz$qPCR_head._ZIKV_WNV <- log10(meta_AAz$qPCR_head._ZIKV_WNV)
meta_AAz$qPCR_head._ZIKV_WNV

meta_AAz$qPCR_body._ZIKV_WNV
meta_AAz$qPCR_body._ZIKV_WNV <- log10(meta_AAz$qPCR_body._ZIKV_WNV)
meta_AAz$qPCR_body._ZIKV_WNV

which(colnames(meta_AAz)=="qPCR_head._ZIKV_WNV")
which(colnames(meta_AAz)=="qPCR_body._ZIKV_WNV")

is.na(meta_AAz[,16:17])<-sapply(meta_AAz[,16:17], is.infinite)
meta_AAz[,16:17][is.na(meta_AAz[,16:17])]<-0
meta_AAz[,16:17]

####
meta_AAz$plaque_assay
meta_AAz$group_n1 <- meta_AAz$plaque_assay
meta_AAz$plaque_assay <- factor(meta_AAz$plaque_assay, levels = c("AA-7dpi-ZIKV(-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                  "AA-21dpi-ZIKV(+)","AA-21dpi-ZIKV(-)","AA-21dpi-blood","AA-21dpi-sucrose"), ordered = TRUE)
meta_AAz$plaque_assay

meta_AAz_order1 <- meta_AAz[meta_AAz$dpi== "7dpi",][with(meta_AAz[meta_AAz$dpi== "7dpi",], 
                                                         order(plaque_assay, -qPCR_head._ZIKV_WNV, -qPCR_body._ZIKV_WNV, -`Zika virus`, -`Guadeloupe mosquito virus`,-`Phasi Charoen-like phasivirus`,
                                                               -`Aedes aegypti toti-like virus`,-`Aedes anphevirus`,-`Guadeloupe mosquito quaranja-like virus 1`)), ]
meta_AAz_order2 <- meta_AAz[meta_AAz$dpi== "21dpi",][with(meta_AAz[meta_AAz$dpi== "21dpi",], order(plaque_assay, -qPCR_head._ZIKV_WNV,  -qPCR_body._ZIKV_WNV, -`Zika virus`, -`Guadeloupe mosquito virus`,-`Phasi Charoen-like phasivirus`,
                                                                                                   -`Aedes aegypti toti-like virus`,-`Aedes anphevirus`,-`Guadeloupe mosquito quaranja-like virus 1`)), ]
meta_AAz_order <- rbind(meta_AAz_order1,meta_AAz_order2)

meta_AAz$Row.names
meta_AAz_order$Row.names
df <- df[,meta_AAz_order$Row.names]
df
col_ha = HeatmapAnnotation(qPCR = anno_lines(meta_AAz_order[,16:17], gp = gpar(col = c("black","darkgrey")), 
                                             height = unit(1.5, "cm"),ylim = c(-0.5,8), axis_param = list(at = c(0,2,4,6,8)),
                                             add_points = TRUE, pt_gp = gpar(col = c("black","darkgrey")), pch = c(16, 17)))

meta_AAz_order$group_n1 

#### column split ####
df_1 <- df
df_1$ZIKV <- "others"
df_1$ZIKV[which(rownames(df_1)=="Zika virus")] <- "aZIKV"

#### draw figure #####
heatmap_col_tp2 <- colorRampPalette(brewer.pal(9, "YlGnBu"), space = "rgb")(100)

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Heatmap_AA_EuVspecies_reads_sampleSum+0_speciesSum+500.new.pdf", width = 14,height = 6,useDingbats = FALSE)
Heatmap(as.matrix(df),
        name = "log10 reads number", #title of legend
        column_title = "Samples", row_title = "viral species",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 6),
        column_title_gp = gpar(fontsize = 12),
        cluster_rows = F,
        cluster_columns= F,
        row_order = species_ord_1,
        # column_order = meta_AAz_order$Row.names,
        column_labels = meta_AAz_order$mosq,
        top_annotation = col_ha,
        col = heatmap_col_tp2,
        column_split = meta_AAz_order$group_n1,
        row_split = df_1$ZIKV,
        # column_title_gp = gpar(fontsize = 8),
        border = TRUE,
        # gap = unit(10, "mm"),
        rect_gp = gpar(col = "white", lwd = 0.5)) 
dev.off()



#################### * 1-4) alpha diversity on species level ####################
plot_richness <- function (physeq, x = "samples", color = NULL, shape = NULL, 
                           title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL, 
                           sortby = NULL) 
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i, 
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures, 
      ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, 
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby, 
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE, 
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color, 
                            shape = shape)
  p = ggplot(mdf, richness_map) + geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
    # geom_point(size=1)
    geom_jitter(position=position_jitter(0.25), size = 3.5, alpha=0.7)
  # if (any(!is.na(mdf[, "se"]))) {
  #   p = p + geom_errorbar(aes(ymax = value + se, ymin = value - 
  #                               se), width = 0.1)
  # }
  # p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, 
                                           # hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

sample_data(GP.M1_species_select_0)$Treatment
GP.M1_species_select_0_7dpi = subset_samples(GP.M1_species_select_0, dpi == "7dpi")
GP.M1_species_select_0_21dpi = subset_samples(GP.M1_species_select_0, dpi == "21dpi")

#### diet 7dpi ##### 
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Alpha-diversity_AA_EuVspecies_diet_7dpiv2.pdf", width = 5,height = 7,useDingbats = FALSE)
plot_richness(GP.M1_species_select_0_7dpi, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") +
  scale_color_manual(values = c(col[9],col[5],col[2])) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of eukaryotic virus in Ae. aegypti at 7dpi")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GP.M1_species_select_0_7dpi, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GP.M1_species_select_0_7dpi)$Treatment, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_species_select_0_7dpi)$Treatment, p.adjust.method = "BH")

#### diet 21dpi ##### 
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Alpha-diversity_AA_EuVspecies_diet_21dpiv2.pdf", width = 5,height = 7,useDingbats = FALSE)
plot_richness(GP.M1_species_select_0_21dpi, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") +
  scale_color_manual(values = c(col[9],col[5],col[2])) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of eukaryotic virus in Ae. aegypti at 21dpi")+
  ylab("Alpha Diversity Value")+
  stat_compare_means()
dev.off()

erich <- estimate_richness(GP.M1_species_select_0_21dpi, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GP.M1_species_select_0_21dpi)$Treatment, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_species_select_0_21dpi)$Treatment, p.adjust.method = "BH")
Kruskal−Wallis
PlantGrowth %>% kruskal_test(weight ~ group)
#### diet all ##### 
sample_data(GP.M1_species_select_0)$Info1

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Alpha-diversity_AA_EuVspecies_diet_allv2_reorder.pdf", width = 6,height = 7,useDingbats = FALSE)
pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Alpha-diversity_AA_EuVspecies_diet_allv2_reorder_rmEVE.pdf", width = 6,height = 7,useDingbats = FALSE)
plot_richness(GP.M1_species_select_0, x="Info1",measures=c("Shannon"),col="Treatment", shape = "dpi") + 
  scale_color_manual(values = rev(c(col[9],col[5],col[2]))) +
  scale_shape_manual(values = c(16, 17))+
  # facet_wrap(~dpi)+
  ylim(-0.001, 2.5)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 12, colour = "black", angle = 90, hjust = 1),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of eukaryotic virus in Ae. aegypti")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GP.M1_species_select_0, measures = c("Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GP.M1_species_select_0)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_species_select_0)$Info1, p.adjust.method = "BH")

#### infectious status all #####

sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000.
sample_data(GP.M1_species_select_0)$hd <- as.character(sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000.)

sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-7dpi-ZIKV(-/+/+)"] <- "7dpe-Head(+)/Body(+)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-7dpi-ZIKV(-/-/+)"] <- "7dpe-Head(-)/Body(+)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-7dpi-ZIKV(-/-/-)"] <- "7dpe-Head(-)/Body(-)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-21dpi-ZIKV(+/+/+)"] <- "21dpe-Head(+)/Body(+)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-21dpi-ZIKV(-/+/+)"] <- "21dpe-Head(+)/Body(+)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-21dpi-ZIKV(-/-/+)"] <- "21dpe-Head(-)/Body(+)"
sample_data(GP.M1_species_select_0)$hd[sample_data(GP.M1_species_select_0)$group.plaque.qPCR.1000.body.1000. == "AA-21dpi-ZIKV(-/-/-)"] <- "21dpe-Head(-)/Body(-)"

sample_data(GP.M1_species_select_0)$hd <- factor(sample_data(GP.M1_species_select_0)$hd,
                                             levels = c("AA-7dpi-sucrose","AA-7dpi-blood","7dpe-Head(-)/Body(-)","7dpe-Head(-)/Body(+)", "7dpe-Head(+)/Body(+)", 
                                                        "AA-21dpi-sucrose","AA-21dpi-blood","21dpe-Head(-)/Body(-)","21dpe-Head(-)/Body(+)", "21dpe-Head(+)/Body(+)"))

sample_data(GP.M1_species_select_0)$hd

pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Alpha-diversity_AA_EuVspecies_ZIKVhd_all_reorder_rmEVE.pdf", width = 7,height = 7,useDingbats = FALSE)
# plot_richness(GP.M1_species_select_0_7dpi, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") + 
plot_richness(GP.M1_species_select_0, x="hd", measures=c("Shannon"), col="hd", shape = "dpi") + 
  scale_color_manual(values = rev(c(col[9],col[7],col[5],col[3],col[1],col[9],col[7],col[5],col[3],col[1]))) +
  scale_shape_manual(values = c(16, 17))+
  ylim(-0.0001,2.5)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust=0.5),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of eukaryotic virus in Ae. aegypti")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GP.M1_species_select_0, measures = c("Shannon"))
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_species_select_0)$hd, p.adjust.method = "BH")
summary(erich$Shannon)
#################### * 1-5) PCoA on species level #####################
#### 7dpi ##### 
GP.M1_species_select_0_7dpi = subset_samples(GP.M1_species_select_0, dpi == "7dpi")
sample_data(GP.M1_species_select_0_7dpi)$group.plaque.qPCR.1000.body.1000.
sample_data(GP.M1_species_select_0_7dpi)$Treatment
GP.M1_species_select_0_rmZIKV <- subset_taxa(GP.M1_species_select_0_7dpi, !taxa_names(GP.M1_species_select_0_7dpi) == "GIM628_NODE_1_length_10744_cov_806.867254")
GP.M1_species_select_0_rmZIKV
sample_variables(GP.M1_species_select_0_rmZIKV)

GP.ord.GP.M1_species_select_0_rmZIKV <- ordinate(GP.M1_species_select_0_rmZIKV, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of EuV on species AA excludeZIKV 7dpi label.pdf", width = 7,height = 6, useDingbats = FALSE)
pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/PCoA of EuV on species AA excludeZIKV 7dpi label rmEVE.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_species_select_0_rmZIKV, GP.ord.GP.M1_species_select_0_rmZIKV, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.1000.body.1000.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_shape_manual(values=rev(shap[2:6])) +
  scale_color_manual(values = c(col[2],col[5],col[9]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of eukaryotic virus in Ae. aegypti at 7dpe excluding ZIKV") +
  stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_species_select_0_rmZIKV, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_species_select_0_rmZIKV))
# Adonis test
adonis(GP.M1_species_select_bray ~ group.plaque.qPCR.1000.body.1000., data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$group.plaque.qPCR.1000.body.1000.)
permutest(beta)


#### 21dpi ##### 
GP.M1_species_select_0_21dpi = subset_samples(GP.M1_species_select_0, dpi == "21dpi")
sample_data(GP.M1_species_select_0_21dpi)$group.plaque.qPCR.1000.
GP.M1_species_select_0_rmZIKV <- subset_taxa(GP.M1_species_select_0_21dpi, !taxa_names(GP.M1_species_select_0_21dpi) == "GIM628_NODE_1_length_10744_cov_806.867254")
GP.M1_species_select_0_rmZIKV
sample_variables(GP.M1_species_select_0_rmZIKV)

GP.ord.GP.M1_species_select_0_rmZIKV <- ordinate(GP.M1_species_select_0_rmZIKV, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of EuV on species AA excludeZIKV 21dpi label.pdf",width = 7,height = 6, useDingbats = FALSE)
pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/PCoA of EuV on species AA excludeZIKV 21dpi label rmEVE.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_species_select_0_rmZIKV, GP.ord.GP.M1_species_select_0_rmZIKV, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.1000.body.1000.")+
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_shape_manual(values=c(shap[1:4],shap[6],shap[5])) +
  # scale_shape_manual(values=c(shap[2:4],shap[6],shap[8])) + 
  scale_color_manual(values = c(col[2],col[5],col[9]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of eukaryotic virus in Ae. aegypti at 21dpe excluding ZIKV")
# stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_species_select_0_rmZIKV, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_species_select_0_rmZIKV))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$group.plaque.qPCR.1000.body.1000.)
permutest(beta)


#################### *** 2 EuV Cx. quinquefasciatus *** ##############################
#################### * 2-1) subset Cx. quinquefasciatus Eukaryotic Virome ####################
GP.M1_CQ <- subset_samples(GP.M1, mosq_species  == "Cx_quinquefasciatus")
GP.M1_CQ
sum(sample_sums(GP.M1_CQ))
sample_variables(GP.M1_CQ)

otu <- as.data.frame(otu_table(GP.M1_CQ))
tax <- as.data.frame(tax_table(GP.M1_CQ))
df <- merge(tax, otu, by="row.names")
write.table(df, file = "/Users/shichenyan/Desktop/mosq_Guad_infection/CQ_tax_abundance.txt", quote = FALSE, sep = "\t", row.names = F)

sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0. <- factor(sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0., 
                                                            levels = c("7dpi-WNV(+/)", "7dpi-WNV(-/+)", "7dpi-WNV(-/-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                                       "14dpi-WNV(+/)","14dpi-WNV(-/+)","14dpi-WNV(-/-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))
sample_data(GP.M1_CQ)$Info1 <- factor(sample_data(GP.M1_CQ)$Info1, levels = c("CQ-7dpi-WNV", "CQ-7dpi-blood", "CQ-7dpi-sucrose", "CQ-14dpi-WNV", "CQ-14dpi-blood", "CQ-14dpi-sucrose"))

levels(sample_data(GP.M1_CQ)$Info1)
levels(sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0.)
#################### * 2-2) merge on species level ####################
GP.M1_species <- GP.M1_CQ %>%
  tax_glom(taxrank = "species")
GP.M1_species

summary(sample_sums(GP.M1_species))
which(sample_sums(GP.M1_species) == 0)
GP.M1_samselect = prune_samples(sample_sums(GP.M1_species) > 0, GP.M1_species)
GP.M1_samselect                                

GP.M1_species_select_0 = prune_taxa(taxa_sums(GP.M1_samselect) > 0, GP.M1_samselect)
GP.M1_species_select_0
summary(taxa_sums(GP.M1_species_select_0))
GP.M1_species_select_0 <- subset_taxa(GP.M1_species_select_0, species != "Aedes aegypti totivirus")
GP.M1_species_select_0

GP.M1_species_select = prune_taxa(taxa_sums(GP.M1_samselect) > 500, GP.M1_samselect)
GP.M1_species_select
summary(taxa_sums(GP.M1_species_select))

GP.M1_species_select <- subset_taxa(GP.M1_species_select, species != "Aedes aegypti totivirus")
GP.M1_species_select
#################### * 2-3) HEATMAP on viral species level ############
#### get table - reads for each viral species ####
OTU_species <- merge(otu_table(GP.M1_species_select), tax_table(GP.M1_species_select), by=0)
rownames(OTU_species) <- OTU_species$species
rownames(OTU_species)
colnames(OTU_species)

#### order viral species by familiy #### 
OTU_species$family <- as.character(OTU_species$family)
OTU_species$family
OTU_species$family[is.na(OTU_species$family)] <- "unclassified"
OTU_species_m<- melt(OTU_species)
OTU_species_m_agg <- aggregate(value~family+species, FUN=sum, data=OTU_species_m)
OTU_species_m_agg1 <- aggregate(value~family, FUN=sum, data=OTU_species_m)
family_ord <- OTU_species_m_agg1[order(-OTU_species_m_agg1$value),]$family
family_ord
family_ord <- c(family_ord[-1],"unclassified")
family_ord
OTU_species_m_agg$index <- 1
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Phenuiviridae")] <- 2
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Metaviridae")] <- 3
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Circoviridae")] <- 4
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Picobirnaviridae")] <- 5
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Iridoviridae")] <- 6
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "Bacilladnaviridae")] <- 7
OTU_species_m_agg$index[which(OTU_species_m_agg$family == "unclassified")] <- 8

OTU_species_m_agg
OTU_species_m_agg <- OTU_species_m_agg[order(OTU_species_m_agg$index, -OTU_species_m_agg$value),]
OTU_species_m_agg

species_ord <- as.character(OTU_species_m_agg$species)
species_ord
which(species_ord == "West Nile virus")
species_ord_1 <- c("West Nile virus", species_ord[-1])
species_ord_1

#### transform the table - reads for each viral species ##### 
rm <- c("Row.names", "superkingdom", "phylum", "class", "order", "family", "genus", "species")
OTU_species <- OTU_species[,!colnames(OTU_species) %in% rm]

#log
OTU_species_trans <- log10(OTU_species)
is.na(OTU_species_trans) <- sapply(OTU_species_trans, is.infinite)
OTU_species_trans[is.na(OTU_species_trans)]<-0
df.original <- OTU_species
df <- OTU_species_trans

#### use meta data for top annotation #####
data.frame(sample_data(OTU_species))
meta = read.csv("GIM_metadata_new.csv", header=TRUE, row.names=1, sep=",")
meta_s = meta[!meta$mosq_species == "",]

meta_AA <- meta_s[colnames(df),]

#### order sample by WNV reads within each group & column annotation  #### 
meta_AAz <- merge(meta_AA, as.data.frame(t(df)),by=0)

sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0. <- factor(sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0., 
                                                            levels = c("7dpi-WNV(+/)", "7dpi-WNV(-/+)", "7dpi-WNV(-/-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                                       "14dpi-WNV(+/)","14dpi-WNV(-/+)","14dpi-WNV(-/-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))
sample_data(GP.M1_CQ)$Info1 <- factor(sample_data(GP.M1_CQ)$Info1, levels = c("CQ-7dpi-WNV", "CQ-7dpi-blood", "CQ-7dpi-sucrose", "CQ-14dpi-WNV", "CQ-14dpi-blood", "CQ-14dpi-sucrose"))

meta_AAz$Info1 <- as.character(meta_AAz$Info1)
meta_AAz$Info1 <- factor(meta_AAz$Info1,levels = c("CQ-7dpi-WNV", "CQ-7dpi-blood", "CQ-7dpi-sucrose", "CQ-14dpi-WNV", "CQ-14dpi-blood", "CQ-14dpi-sucrose"))

meta_AAz$group.plaque.qPCR.0.body.0 <- as.character(meta_AAz$group.plaque.qPCR.0.body.0)
meta_AAz$group.plaque.qPCR.0.body.0 <- factor(meta_AAz$group.plaque.qPCR.0.body.0,
                                              levels = c("7dpi-WNV(+/)", "7dpi-WNV(-/+)", "7dpi-WNV(-/-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                         "14dpi-WNV(+/)","14dpi-WNV(-/+)","14dpi-WNV(-/-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))

meta_AAz$plaque_assay <- as.character(meta_AAz$plaque_assay)
meta_AAz$plaque_assay <- factor(meta_AAz$plaque_assay,
                                levels = c("7dpi-WNV(+)", "7dpi-WNV(-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                           "14dpi-WNV(+)","14dpi-WNV(-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))

levels(meta_AAz$Info1)
levels(meta_AAz$group.plaque.qPCR.0.body.0)
levels(meta_AAz$plaque_assay)

meta_AAz$qPCR_body._ZIKV_WNV
meta_AAz$qPCR_body._ZIKV_WNV <- log10(meta_AAz$qPCR_body._ZIKV_WNV)
meta_AAz$qPCR_body._ZIKV_WNV

# which(colnames(meta_AAz)=="qPCR_head._ZIKV_WNV")
which(colnames(meta_AAz)=="qPCR_body._ZIKV_WNV")

is.na(meta_AAz[,17])<-sapply(meta_AAz[,17], is.infinite)
meta_AAz[,17][is.na(meta_AAz[,17])]<-0
meta_AAz[,17]

####
meta_AAz$group.plaque.qPCR.0.body.0
meta_AAz$group_n1 <- meta_AAz$group.plaque.qPCR.0.body.0
meta_AAz$group.plaque.qPCR.0.body.0 <- factor(meta_AAz$group.plaque.qPCR.0.body.0, 
                                              levels = c("7dpi-WNV(+/)", "7dpi-WNV(-/+)", "7dpi-WNV(-/-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                         "14dpi-WNV(+/)","14dpi-WNV(-/+)","14dpi-WNV(-/-)","CQ-14dpi-blood","CQ-14dpi-sucrose"), ordered = TRUE)
meta_AAz$group.plaque.qPCR.0.body.0

meta_AAz_order1 <- meta_AAz[meta_AAz$dpi== "7dpi",][with(meta_AAz[meta_AAz$dpi== "7dpi",], 
                                                         order(group.plaque.qPCR.0.body.0, -qPCR_body._ZIKV_WNV, -`West Nile virus`, -`Guadeloupe Culex tymo-like virus`,
                                                               -`Wenzhou sobemo-like virus 3`, -`Phasi Charoen-like phasivirus`)), ]
meta_AAz_order2 <- meta_AAz[meta_AAz$dpi== "14dpi",][with(meta_AAz[meta_AAz$dpi== "14dpi",], 
                                                          order(group.plaque.qPCR.0.body.0, -qPCR_body._ZIKV_WNV, -`West Nile virus`, -`Guadeloupe Culex tymo-like virus`,
                                                                -`Wenzhou sobemo-like virus 3`, -`Phasi Charoen-like phasivirus`)), ]
meta_AAz_order <- rbind(meta_AAz_order1,meta_AAz_order2)

meta_AAz$Row.names
meta_AAz_order$Row.names
df <- df[,meta_AAz_order$Row.names]
df
col_ha = HeatmapAnnotation(qPCR = anno_lines(meta_AAz_order[,17], gp = gpar(col = c("darkgrey")), 
                                             height = unit(1.5, "cm"),ylim = c(-0.5,6), axis_param = list(at = c(0,2,4,6)),
                                             add_points = TRUE, pt_gp = gpar(col = c("darkgrey")), pch = c(17)))

meta_AAz_order$group_n1 
meta_AAz_order$plaque_assay
# meta_AAz_order$group_n1 <- factor(meta_AAz_order$group_n1,levels = c("AA-7dpi-ZIKV(-/+)","AA-7dpi-ZIKV(-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
#                                                                                            "AA-21dpi-ZIKV(+/+)","AA-21dpi-ZIKV(-/+)","AA-21dpi-ZIKV(-/-)","AA-21dpi-blood","AA-21dpi-sucrose"),ordered = F)


#### column split ####
df_1 <- df
df_1$ZIKV <- "others"
df_1$ZIKV[which(rownames(df_1)=="West Nile virus")] <- "aZIKV"

#### draw figure #####
# col = c(brewer.pal(4, "Blues")[3],brewer.pal(4, "Reds")[2:4],brewer.pal(12, "Paired")[3:4], brewer.pal(n=4,name="Purples")[3],brewer.pal(12, "Paired")[10])
heatmap_col_tp2 <- colorRampPalette(brewer.pal(9, "YlGnBu"), space = "rgb")(100)

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Heatmap_CQ_EuVspecies_reads_sampleSum+0_speciesSum+500.pdf", width = 10,height = 3.5,useDingbats = FALSE)
Heatmap(as.matrix(df),
        name = "log10 reads number", #title of legend
        column_title = "Samples", row_title = "viral species",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 6),
        column_title_gp = gpar(fontsize = 12),
        cluster_rows = F,
        cluster_columns= F,
        row_order = species_ord_1,
        # column_order = meta_AAz_order$Row.names,
        column_labels = meta_AAz_order$mosq,
        top_annotation = col_ha,
        col = heatmap_col_tp2,
        column_split = meta_AAz_order$plaque_assay,
        row_split = df_1$ZIKV,
        # column_title_gp = gpar(fontsize = 8),
        border = TRUE,
        # gap = unit(10, "mm"),
        rect_gp = gpar(col = "white", lwd = 0.5)) 
dev.off()



#################### *** 3 DESeq2 differental virus *** #######################
#################### * 3-1）Aedes ####################
df_otu <- as.data.frame(GP.M1_species_select_0@otu_table)
df_tax <- as.data.frame(GP.M1_species_select_0@tax_table)
df_meta <- data.frame(sample_data(GP.M1_species_select_0))
df <- merge(df_otu, df_tax, by = "row.names", all=T)
head(df)
#### treatment ####
AAdiet <- readRDS("./RDS/DESeq2_AA_treatment.RDS")
AAdiet

#cutoff padj<0.05 basemean >50
sp <- c("Guadeloupe mosquito virus", "Beihai barnacle virus 12", "Aedes aegypti totivirus", "Aedes aegypti toti-like virus")

df_sp <- df[df$species %in% sp,]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "species")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0

df_merge$species <- factor(df_merge$species, levels = sp); levels(df_merge$species)

##  plot
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/AA_treatment_DESeq.pdf", width = 9, height = 6,useDingbats = FALSE)
ggplot(df_merge, aes(x = Treatment, y = log10reads, fill = Treatment)) + 
  facet_wrap(~species,nrow =1, scales="fixed") + 
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(col[9], col[5], col[2]))+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, hjust = 0.5),
        legend.position = "bottom",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 raw reads")
dev.off()

#### infectious status ####
AAinfec <- readRDS("./RDS/DESeq2_AA_infectious.RDS")
AAinfec

# "Aedes aegypti toti-like virus" 7headZIKVposvsneg
# "Beihai barnacle virus 12"      21plaqueposvsneg
# " Aedes aegypti toti-like virus" 21headZIKVposvsneg

cog <- headZIKVposvsneg
sp2 <- c("Aedes aegypti toti-like virus")
, "Beihai barnacle virus 12",
sp <- c("Guadeloupe mosquito virus", "Beihai barnacle virus 12", "Aedes aegypti totivirus", "Aedes aegypti toti-like virus")

df_sp <- df[df$species == "Aedes aegypti toti-like virus",]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "species")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0

df_merge <- df_merge[df_merge$head_ZIKV.1000 %in% c("AA-7dpi-ZIKV(-)","AA-7dpi-ZIKV(+)", "AA-21dpi-ZIKV(-)", "AA-21dpi-ZIKV(+)"),]
df_merge$head_ZIKV.1000 <- factor(df_merge$head_ZIKV.1000, levels = c("AA-7dpi-ZIKV(+)","AA-7dpi-ZIKV(-)", "AA-21dpi-ZIKV(+)", "AA-21dpi-ZIKV(-)"))

##  plot
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/AA_infectious_DESeq.pdf", width = 4.5, height = 6,useDingbats = FALSE)
ggplot(df_merge, aes(x = head_ZIKV.1000, y = log10reads, fill = head_ZIKV.1000)) + 
  facet_wrap(~species,nrow =1, scales="fixed") +
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(col[9], col[3], col[8], col[2]))+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, hjust = 0.5),
        legend.position = "bottom",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 raw reads")
dev.off()

#################### * 3-2）Culex ####################
df_otu <- as.data.frame(GP.M1_species_select_0@otu_table)
df_tax <- as.data.frame(GP.M1_species_select_0@tax_table)
df_meta <- data.frame(sample_data(GP.M1_species_select_0))
df <- merge(df_otu, df_tax, by = "row.names", all=T)
head(df)

#### treatment ####
CQdiet <- readRDS("./RDS/DESeq2_CQ_treatment.RDS")
CQdiet

#cutoff padj<0.05 basemean >50
sp <- c("Guadeloupe Culex tymo-like virus")

df_sp <- df[df$species %in% sp,]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "species")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0

df_merge$Treatment <- factor(df_merge$Treatment, levels = c("CQ-WNV", "CQ-blood", "CQ-sucrose"))
# df_merge$species <- factor(df_merge$species, levels = sp); levels(df_merge$species)

##  plot
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/CQ_treatment_DESeq.pdf", width = 3, height = 6,useDingbats = FALSE)
ggplot(df_merge, aes(x = Treatment, y = log10reads, fill = Treatment)) + 
  facet_wrap(~species,nrow =1, scales="fixed") + 
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(col[9], col[5], col[2]))+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, hjust = 0.5),
        legend.position = "bottom",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 raw reads")
dev.off()

#### infectious status ####
CQinfec <- readRDS("./RDS/DESeq2_CQ_infectious.RDS")
CQinfec

df_sp <- df[df$species == "Guadeloupe Culex tymo-like virus",]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "species")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0

df_merge <- df_merge[df_merge$plaque_assay %in% c("7dpi-WNV(+)","7dpi-WNV(-)", "14dpi-WNV(+)", "14dpi-WNV(-)"),]
df_merge$plaque_assay <- factor(df_merge$plaque_assay, levels = c("7dpi-WNV(+)","7dpi-WNV(-)", "14dpi-WNV(+)", "14dpi-WNV(-)"))

##  plot
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/CQ_infectious_DESeq.pdf", width = 4.5, height = 6,useDingbats = FALSE)
ggplot(df_merge, aes(x = plaque_assay, y = log10reads, fill = plaque_assay)) + 
  facet_wrap(~species,nrow =1, scales="fixed") +
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(col[9], col[3], col[8], col[2]))+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, hjust = 0.5),
        legend.position = "bottom",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 raw reads")
dev.off()

#################### *** 4 qRT-PCR analysis *** ####################
setwd("/Users/shichenyan/Desktop/mosq_Guad_infection/")
copies = read.table("./qRT-PCR/qRTPCR_copies.csv", header=TRUE, dec=".", sep=",")
colnames(copies)

## log transform
copies_log <-copies
copies_log <- cbind(copies[,1:7], log10(copies[,8:15]),copies[,16:23])
is.na(copies_log) <- sapply(copies_log, is.infinite)
copies_log[is.na(copies_log)] <- 0
colnames(copies_log)

########### all core viruses #########
copies_log$plaque_assay <- gsub("CQ-14dpi-","CQ-",copies_log$plaque_assay,fixed = T)
copies_log$plaque_assay <- gsub("CQ-7dpi-","CQ-",copies_log$plaque_assay,fixed = T)
copies_log$plaque_assay <- gsub("14dpi-","CQ-",copies_log$plaque_assay,fixed = T)
copies_log$plaque_assay <- gsub("7dpi-","CQ-",copies_log$plaque_assay,fixed = T)
levels(copies_log$plaque_assay)
copies_log$plaque_assay <- factor(copies_log$plaque_assay, levels = c("AA-sucrose", "AA-blood", "AA-ZIKV(-)", "AA-ZIKV(+)", 
                                                                      "CQ-sucrose", "CQ-blood", "CQ-WNV(-)", "CQ-WNV(+)"))
levels(copies_log$plaque_assay)

copies_log$PCLPV_copies[copies_log$mosq_species == "Cx_quinquefasciatus"] <- NA
copies_log$GMV_copies[copies_log$mosq_species == "Cx_quinquefasciatus"] <- NA
copies_log$ATV_copies[copies_log$mosq_species == "Cx_quinquefasciatus"] <- NA
copies_log$AANV_copies[copies_log$mosq_species == "Cx_quinquefasciatus"] <- NA
copies_log$GCTLV_copies[copies_log$mosq_species == "Ae_aegypti"] <- NA
copies_log$WSLV3_copies[copies_log$mosq_species == "Ae_aegypti"] <- NA

copies_log_m <- melt(copies_log)
copies_log_m_sub <- copies_log_m[copies_log_m$variable %in% c("PCLPV_copies", "GMV_copies", "ATV_copies","AANV_copies","GCTLV_copies", "WSLV3_copies"),]
copies_log_m_sub <- na.omit(copies_log_m_sub)

copies_log_m_sub$variable
copies_log_m_sub$variable <- factor(copies_log_m_sub$variable, levels=c("GMV_copies", "PCLPV_copies","ATV_copies","AANV_copies","GCTLV_copies", "WSLV3_copies"))

copies_log_m_sub$virus[copies_log_m_sub$variable == "GMV_copies"] <- "Guadeloupe mosquito virus"
copies_log_m_sub$virus[copies_log_m_sub$variable == "PCLPV_copies"] <- "Phasi Charoen−like phasivirus"
copies_log_m_sub$virus[copies_log_m_sub$variable == "ATV_copies"] <- "Aedes aegypti totivirus"
copies_log_m_sub$virus[copies_log_m_sub$variable == "AANV_copies"] <- "Aedes anphevirus"
copies_log_m_sub$virus[copies_log_m_sub$variable == "GCTLV_copies"] <- "Guadeloupe Culex tymo−like virus"
copies_log_m_sub$virus[copies_log_m_sub$variable == "WSLV3_copies"] <- "Wenzhou sobemo−like virus 3"
copies_log_m_sub$virus <- factor(copies_log_m_sub$virus, levels=c("Guadeloupe mosquito virus", "Phasi Charoen−like phasivirus","Aedes aegypti totivirus",
                                                                  "Aedes anphevirus","Guadeloupe Culex tymo−like virus", "Wenzhou sobemo−like virus 3"))

copies_log_m_sub$dpi <- gsub("dpi", "dpe", copies_log_m_sub$dpi)
copies_log_m_sub$dpi<- factor(copies_log_m_sub$dpi, levels = c("7dpe", "21dpe", "14dpe"))


copies_log_m_sub_rmATV21 <- copies_log_m_sub[!(copies_log_m_sub$dpi=="21dpe"&copies_log_m_sub$variable=="ATV_copies"), ] 
copies_log_m_sub_rmGCTLV7 <- copies_log_m_sub_rmATV21[!(copies_log_m_sub_rmATV21$dpi=="7dpe"&copies_log_m_sub_rmATV21$variable=="GCTLV_copies"), ] 
copies_log_m_sub_rmWSLV37 <- copies_log_m_sub_rmGCTLV7[!(copies_log_m_sub_rmGCTLV7$dpi=="7dpe"&copies_log_m_sub_rmGCTLV7$variable=="WSLV3_copies"), ] 

col = plasma(8,alpha = 1, begin=0.9, end = 0, direction = 1)

pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/qRTPCR_allcoreviruses_plaque_dpe_reorder_rmdpe.pdf", width = 12,height =9,useDingbats = FALSE)
# pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/qRTPCR_allcoreviruses_plaque.pdf", width = 13,height =9,useDingbats = FALSE)
ggplot(copies_log_m_sub_rmWSLV37, aes(x=plaque_assay, y=value, col=plaque_assay)) +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA) + geom_jitter(shape=19, position=position_jitter(0.25), size = 2.2, alpha=0.7)+
  theme_bw()+
  facet_wrap(~ virus+part, scales="free", ncol = 6)+
  # geom_line(color = "darkgrey", size = 0.7 ) +
  # geom_point(aes(shape=variable, color=variable),alpha=0.8, size = 2)+
  scale_color_manual(values = col) +
  # scale_color_manual(values = c(col[9], col[7], col[4],col[1])) +
  # scale_y_continuous(trans="log10")+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=10), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12))+
  ggtitle("Viral genome copies determined by qRT-PCR")+
  ylab("log10 copies")+
  scale_y_continuous(limits = c(-0.01,9), breaks=seq(0, 9, by = 1))
dev.off()

copies_log_m_sub_rmATV7 <- copies_log_m_sub[!(copies_log_m_sub$dpi=="7dpe"&copies_log_m_sub$variable=="ATV_copies"), ] 
copies_log_m_sub_rmGCTLV21 <- copies_log_m_sub_rmATV7[!(copies_log_m_sub_rmATV7$dpi=="14dpe"&copies_log_m_sub_rmATV7$variable=="GCTLV_copies"), ] 
copies_log_m_sub_rmWSLV321 <- copies_log_m_sub_rmGCTLV21[!(copies_log_m_sub_rmGCTLV21$dpi=="14dpe"&copies_log_m_sub_rmGCTLV21$variable=="WSLV3_copies"), ] 

copies_log_m_sub_ATV <- copies_log_m_sub[copies_log_m_sub$variable=="ATV_copies", ]
copies_log_m_sub_ATV$dpi <- "7dpe+21dpe"
copies_log_m_sub_GCTLV <- copies_log_m_sub[copies_log_m_sub$variable=="GCTLV_copies", ]
copies_log_m_sub_GCTLV$dpi <- "7dpe+14dpe"
copies_log_m_sub_WSLV3 <- copies_log_m_sub[copies_log_m_sub$variable=="WSLV3_copies", ]
copies_log_m_sub_WSLV3$dpi <- "7dpe+14dpe"

copies_log_m_sub_new <- rbind(copies_log_m_sub_rmWSLV321,copies_log_m_sub_ATV)
copies_log_m_sub_new <- rbind(copies_log_m_sub_new,copies_log_m_sub_GCTLV)
copies_log_m_sub_new <- rbind(copies_log_m_sub_new,copies_log_m_sub_WSLV3)

pdf("/Users/shichenyan/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/qRTPCR_allcoreviruses_plaque_dpe_reorder_all_new.pdf", width = 9,height =16,useDingbats = FALSE)
ggplot(copies_log_m_sub_new, aes(x=plaque_assay, y=value, col=plaque_assay)) +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA) + geom_jitter(shape=19, position=position_jitter(0.25), size = 2.2, alpha=0.7)+
  theme_bw()+
  facet_wrap(~ virus+dpi+part, scales="free", ncol = 4)+
  # geom_line(color = "darkgrey", size = 0.7 ) +
  # geom_point(aes(shape=variable, color=variable),alpha=0.8, size = 2)+
  scale_color_manual(values = col) +
  # scale_color_manual(values = c(col[9], col[7], col[4],col[1])) +
  # scale_y_continuous(trans="log10")+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=10), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12))+
  ggtitle("Viral genome copies determined by qRT-PCR")+
  ylab("log10 copies")+
  scale_y_continuous(limits = c(-0.01,9), breaks=seq(0, 9, by = 1))
dev.off()

########### pairwise.wilcox.test #######
head(copies_log)
  
pairwise.wilcox.test(copies_log$GMV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GMV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GMV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GMV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")

pairwise.wilcox.test(copies_log$ATV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$ATV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$ATV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$ATV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")

pairwise.wilcox.test(copies_log$AANV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$AANV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$AANV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$AANV_copies[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Ae_aegypti" & copies_log$dpi=="21dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")

pairwise.wilcox.test(copies_log$GCTLV_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GCTLV_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GCTLV_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$GCTLV_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")

pairwise.wilcox.test(copies_log$WSLV3_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$WSLV3_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="7dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$WSLV3_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="body"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="body"], 
                     p.adjust.method = "BH")
pairwise.wilcox.test(copies_log$WSLV3_copies[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="head"], 
                     copies_log$plaque_assay[copies_log$mosq_species=="Cx_quinquefasciatus" & copies_log$dpi=="14dpi" & copies_log$part=="head"], 
                     p.adjust.method = "BH")
#################### *** 5 Phageome Phyloseq Project *** #######################
setwd("/Users/shichenyan/Desktop/mosq_Guad_infection/")
GuadPhage <- readRDS("./RDS/GIM_phagome_raw.RDS")
#################### * 5-1) length versus abundance distribution #################
df <- as.data.frame(GuadPhage@otu_table)
head(df)
rownames(df) <- gsub("GIM319_GIM319_", "GIM319_", rownames(df), fixed = T)
df$name <- rownames(df)
library(tidyr)
df <- df %>% separate(name, c("a","b","c","d","length","f","g"), sep ="_"); head(df)
df <- df[,!colnames(df) %in% c("a","b","c","d","f","g")]; head(df)
df_m <- melt(df)
df$length

df_meta <- data.frame(sample_data(GuadPhage))
head(df_meta)
df_meta <- df_meta[,6:7]; head(df_meta)
df_meta$variable <- rownames(df_meta); head(df_meta)
df_merge <- merge(df_m, df_meta, by="variable"); head(df_merge)
df_merge$length <-  as.numeric(df_merge$length)
summary(df_merge$length)
df_merge$logvalue <- log10(df_merge$value)
df_merge$logvalue[is.infinite(df_merge$logvalue)] <- 0
df_merge_rm <- df_merge[df_merge$logvalue>0,]

col = viridis(8,alpha = 0.6, begin=0.8, end = 0.1, direction = 1)

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Phage_length_vs_abundance_distribution.pdf", width = 8,height =7,useDingbats = FALSE)
ggplot(df_merge_rm, aes(x=log10(length), y=logvalue, shape=mosq_species, color=mosq_species)) +
  geom_point(size=1.8)+
  geom_smooth()+
  theme_bw()+
  scale_color_manual(values = c(col[7],col[3]))+
  ylab("log10 Abundance")+
  xlab("log10 Length of phage contigs")
dev.off()  

#################### * 5-2) HEATMAP of phageome ############
#### load phyloseq project #### 
Phage_host = readRDS("./RDS/GIM_phage_host.RDS")
df <- as.data.frame(otu_table(Phage_host))
df_tax <- as.data.frame(tax_table(Phage_host))
df_meta <- data.frame(sample_data(Phage_host))
# complete <- read.table("~/OneDrive/1. Project + Guadeloupe mosq/Infection/checkv/completeness.tsv", row.names = 1, header=TRUE, sep="\t"); head(complete)
# complete <- complete[ ,4:5]
# rownames(complete) <- gsub("GIM319_GIM319_", "GIM319_", rownames(complete))
# df_meta_complete <- merge(df_meta, complete , by="row.names" ,all =T)
# rownames(df_meta_complete) <- df_meta_complete$Row.names
# df_meta_complete <-  df_meta_complete[,-1]
# write.csv(df_meta_complete, "./PhageHostPrediction/Phage_host_completeness.csv", quote = F)

Phage_host_complete40 <- subset_samples(Phage_host, aai_completeness >=40)
sum(sample_sums(Phage_host_complete40))
sum(sample_sums(Phage_host))


cat <- c("AA-ZIKV", "AA-blood", "AA-sucrose", "CQ-WNV", "CQ-blood", "CQ-sucrose")
#### mat1 log10 reads per treatment #####
Phage_host_treatment <- Phage_host %>%
  tax_glom(taxrank = "Treatment")

Phage_host_treatment_merge <- merge(as.data.frame(Phage_host_treatment@otu_table), as.data.frame(Phage_host_treatment@tax_table), by="row.names")
Phage_host_treatment_merge
rownames(Phage_host_treatment_merge) <- Phage_host_treatment_merge$Treatment
Phage_host_treatment_merge <- Phage_host_treatment_merge[, !colnames(Phage_host_treatment_merge) %in% c("Row.names", "mosq_species", "Treatment","dpi", "Info1")] 
colnames(Phage_host_treatment_merge)
mat1 <- log10(Phage_host_treatment_merge)
is.na(mat1) <- sapply(mat1, is.infinite)
mat1[is.na(mat1)]<-0 
mat1[1:5,1:5]

rownames(mat1)
mat1 <- mat1[cat,]
rownames(mat1)

#### mat2 presence% per treatment #####
otu_table(Phage_host) <- data.frame(sample_data(Phage_host))
Phage_host_pre <- Phage_host %>%
  tax_glom(taxrank = "Treatment")
Phage_host_pre_merge <- merge(as.data.frame(Phage_host_pre@otu_table), as.data.frame(Phage_host_pre@tax_table), by="row.names")
Phage_host_pre_merge
rownames(Phage_host_pre_merge) <- Phage_host_pre_merge$Treatment
Phage_host_pre_merge <- Phage_host_pre_merge[, !colnames(Phage_host_pre_merge) %in% c("Row.names", "mosq_species", "Treatment","dpi", "Info1")] 
colnames(Phage_host_pre_merge)
mat2 <- Phage_host_pre_merge
mat2 <- mat2[cat,]
rownames(mat2)

table(df_tax$Treatment)
mat2[1,] <- mat2[1,]/74*100
mat2[2,] <- mat2[2,]/10*100
mat2[3,] <- mat2[3,]/10*100
mat2[4,] <- mat2[4,]/40*100
mat2[5,] <- mat2[5,]/10*100
mat2[6,] <- mat2[6,]/10*100


##abundance log transform
dflog <- log10(df)
is.na(dflog) <- sapply(dflog, is.infinite)
dflog[is.na(dflog)]<-0

####  order rows-sample and columns-contigs #### 
##column - phage host
host_or <- data.frame(sort(table(df_meta$host_genus), decreasing=T)); host_or
df_meta$host_genus <- factor(df_meta$host_genus, levels = host_or$Var1); levels(df_meta$host_genus)

##column - viral annotation
family_or <- data.frame(sort(table(df_meta$family), decreasing=T)); family_or
df_meta$family <- factor(df_meta$family, levels = family_or$Var1); levels(df_meta$family)

##column - sample order-length
df_meta$log10length <- log10(df_meta$length); df_meta$log10length
contigs_ordered <- df_meta[with(df_meta, order(host_genus, -log10length)),]
df_meta$length1k <- df_meta$length/1000

########  draw the map ########
heatmap_col_tp1 <- colorRamp2(c(0:5), brewer.pal(6, "YlOrBr"), space = "RGB") #heatmap colour
heatmap_col_tp2 <- colorRamp2(c(0,20,40,60,80,100), brewer.pal(6, "YlGnBu"), space = "RGB")
heatmap_col_tp3 <- colorRamp2(c(0,20,40,60,80,100), brewer.pal(6, "Greens"), space = "RGB")

#### host block ####  
colh = c("#818687FF", "#2D5C21FF", "#2CA0A9FF", "#A1E4BDFF", "#323B34FF", "#C8D656FF")
colhost = HeatmapAnnotation(host = anno_block(gp = gpar(fill = colh)),show_annotation_name = F, annotation_label = F)  ## HeatmapAnnotation default annotaion is col

#### viral annotation ####  
colf = c("#818687FF", "#66C2A5", "#FC8D62", "#8DA0CB", "#FFD92F", "#E5C494")
colfamily = HeatmapAnnotation(condition = df_meta$family,
                              col = list(condition = c("no annotation"=colf[1], "Siphoviridae"=colf[2], "Myoviridae"=colf[3],
                                                       "Podoviridae"=colf[4], "Microviridae"=colf[5], "Herelleviridae"=colf[6]),
                                                       show_annotation_name = T, annotation_label = F, 
                                                       border = T))

#### contigs length ####  
collen = HeatmapAnnotation(length = anno_barplot(df_meta$length1k, ylim = c(0,50),
                                                 height = unit(2, "cm")))
df_meta$aai_completeness[is.na(df_meta$aai_completeness)] <- 0
complete = HeatmapAnnotation(complete = df_meta$aai_completeness, col = list(complete=heatmap_col_tp3))

#### legend ####  
lgd_host = Legend(labels = levels(df_meta$host_genus), title = "Predicted_host", legend_gp = gpar(fill = colh))
lgd_virus = Legend(labels = levels(df_meta$family), title = "Phage_annotation", legend_gp = gpar(fill = colf))
lgd_reads = Legend(title = "Log10 reads number", col = heatmap_col_tp1, at = c(0:6), labels = c("0", "1", "2", "3", "4", "5", "6"))
lgd_pro = Legend(title = "Present proportion", col = heatmap_col_tp2, at = c(0,20,40,60,80,100), labels = c("0", "20", "40", "60", "80", "100"))
pd = packLegend(lgd_host, lgd_virus, 
                # lgd_reads, lgd_pro, 
                direction = "vertical")

#### plot ####  
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/phage_host_wheatmap_complete.pdf", width = 16,height = 6,useDingbats = FALSE)
ht1 <- Heatmap(as.matrix(mat1),
               name = "log10 reads number", #title of legend
               row_title = "Samples",
               col = heatmap_col_tp1,
               column_split = df_meta$host_genus, 
               # row_split = df_tax$Treatment,
               row_names_gp = gpar(fontsize = 12),
               column_names_gp = gpar(fontsize = 0),
               column_order = contigs_ordered$name,
               row_order = cat,
               # column_title_gp = gpar(fontsize = 12),
               cluster_rows = F,
               cluster_columns= F,
               rect_gp = gpar(col = "white", lwd = 0.5))
ht2 <- Heatmap(as.matrix(mat2),
               name = "Present proportion", #title of legend
               row_title = "Samples",
               col = heatmap_col_tp2,
               column_split = df_meta$host_genus, 
               # row_split = df_tax$Treatment,
               row_names_gp = gpar(fontsize = 12),
               column_names_gp = gpar(fontsize = 0),
               column_order = contigs_ordered$name,
               row_order = cat,
               # column_title_gp = gpar(fontsize = 12),
               cluster_rows = F,
               cluster_columns= F,
               rect_gp = gpar(col = "white", lwd = 0.5))
 
ht_list = colhost %v% colfamily  %v% complete %v% collen %v% ht1 %v% ht2
# draw(ht_list, annotation_legend_list = pd)
ht_draw = draw(ht_list, annotation_legend_list = pd)
decorate_annotation("length", {
  grid.lines(c(0, 280), c(10, 10), gp = gpar(lty = 1, col = "lightgrey"), default.units = "native")
  grid.lines(c(0, 280), c(20, 20), gp = gpar(lty = 1, col = "lightgrey"), default.units = "native")
  grid.yaxis(at = c(0, 10, 20, 30, 40, 50))
})
dev.off()


#################### *** 6 Phageome Aedes aegypti *** #######################
GuadPhage_samselect = prune_samples(sample_sums(GuadPhage) > 0 , GuadPhage)
GuadPhage_samselect

factor(sample_data(GuadPhage_samselect)$Treatment)

sample_data(GuadPhage_samselect)$Treatment <- factor(sample_data(GuadPhage_samselect)$Treatment, 
                                                     levels = c("AA-sucrose","AA-blood", "AA-ZIKV", "CQ-sucrose", "CQ-blood", "CQ-WNV"))

phage <- scan(file = "~/OneDrive/1. Project + Guadeloupe mosq/phage_contigs_complete+20.txt", what="", sep="\n")
my_subset <- subset(otu_table(GuadPhage_samselect), rownames(otu_table(GuadPhage_samselect)) %in% phage)
new_GuadPhage_samselect <- merge_phyloseq(my_subset, tax_table(GuadPhage_samselect), sample_data(GuadPhage_samselect))

pdf("~/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Alpha_diversity_phage_all_diet_complete+20reorder.pdf", width = 7,height = 6,useDingbats = FALSE)
plot_richness(new_GuadPhage_samselect, x="Treatment",measures=c("Shannon"),col="Treatment", shape = "mosq_species") + 
  scale_color_manual(values = rev(c(col[9],col[5],col[2],col[9],col[5],col[2]))) +
  theme_bw()+
  facet_wrap(~mosq_species, scales = "free")+
  scale_shape_manual(values = c(16, 17))+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.ticks.x = element_blank()) +
  ggtitle("Alpha diversity of phages")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(new_GuadPhage_samselect, measures = c("Shannon"))
pairwise.wilcox.test(erich$Shannon,sample_data(new_GuadPhage_samselect)$Treatment, p.adjust.method = "BH")

#################### * color ####################
col = plasma(9,alpha = 1, begin=0.9, end = 0, direction = 1)
#### Phyloseq Project #### 
GuadPhage_Ae_aegypti <- subset_samples(GuadPhage, mosq_species  == "Ae_aegypti" ) 
GuadPhage_Ae_aegypti
GP.M1_AA = GuadPhage_Ae_aegypti
GP.M1_AA
sum(sample_sums(GP.M1_AA))
sample_data(GP.M1_AA)$Treatment <- factor(sample_data(GP.M1_AA)$Treatment, levels = c("AA-ZIKV","AA-blood","AA-sucrose"))

# a <- as.data.frame(otu_table(GuadPhage_Ae_aegypti))
# c <- as.data.frame(otu_table(GuadPhage_Cx_quinquefasciatus))
# head(a)
# head(c)
# 
# al <- merge(a, c, by="row.names")
# write.table(al, file = "phage.txt", quote = FALSE, sep = "\t", row.names = T)
# colnames(a)

#### subset #### 
sample_sums(GP.M1_AA)
sample_sums(GP.M1_AA)[which(sample_sums(GP.M1_AA) == 0)]
GP.M1_AA_samselect = prune_samples(sample_sums(GP.M1_AA) > 0 , GP.M1_AA)
GP.M1_AA_samselect

summary(taxa_sums(GP.M1_AA_samselect))
GP.M1_AA_select_0 = prune_taxa(taxa_sums(GP.M1_AA_samselect) > 0, GP.M1_AA_samselect)
GP.M1_AA_select_0

#################### * 6-1) alpha diversity ####################
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Alpha_diversity_phage_AA_dietv2.pdf", width = 6,height = 7,useDingbats = FALSE)
plot_richness(GP.M1_AA_samselect, x="Info1",measures=c("Shannon"),col="Treatment", shape = "mosq_species") + 
  # scale_color_manual(values = c(col[9],col[5],col[2],col[9],col[5],col[2])) + 
  theme_bw()+
  facet_wrap(~mosq_species, scales = "free")+
  scale_shape_manual(values = c(16, 17))+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.ticks.x = element_blank()) +
  ggtitle("Alpha diversity of phages in Ae. aegypti")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GP.M1_AA_samselect, measures = c("Shannon"))
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_AA_samselect)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_AA_samselect)$Treatment, p.adjust.method = "BH")

#################### * 6-2) PCoA #####################
#### all ####
GP.ord.GP.M1_AA_samselect <- ordinate(GP.M1_AA_samselect, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level AA raw reads.pdf",
    width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_AA_samselect, GP.ord.GP.M1_AA_samselect, type = "samples", col = "Treatment")+#,label= "name") +
  theme_bw() +
  geom_point(size = 4, alpha=0.7)  + 
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  # scale_shape_manual(values=c(shap[1],shap[3],shap[5])) +
  scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in all Ae. aegypti samples")
  # stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()
set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_CQ_select_0, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

##### 7dpi #####
GP.M1_AA_select_0_7dpi = subset_samples(GP.M1_AA_samselect, dpi == "7dpi")
GP.M1_AA_select_7dpi = prune_taxa(taxa_sums(GP.M1_AA_select_0_7dpi) > 0, GP.M1_AA_select_0_7dpi)
levels(sample_data(GP.M1_AA_select_7dpi)$group.plaque.qPCR.1000.body.1000.)

GP.ord.GP.M1_AA_select_7dpi <- ordinate(GP.M1_AA_select_7dpi, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level AA raw reads ovals 7dpi v2.pdf",
    width = 7.2,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_AA_select_7dpi, GP.ord.GP.M1_AA_select_7dpi, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.1000.body.1000.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4, alpha=0.7)  + 
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  # scale_shape_manual(values=c(shap[1:4],shap[6])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in Ae. aegypti at 7dpe")
  # stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment), color="black")
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_AA_select_0_7dpi, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0_7dpi))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

##### 21dpi #####
GP.M1_AA_select_0_21dpi = subset_samples(GP.M1_AA_samselect, dpi == "21dpi")
GP.M1_AA_select_21dpi = prune_taxa(taxa_sums(GP.M1_AA_select_0_21dpi) > 0, GP.M1_AA_select_0_21dpi)
levels(sample_data(GP.M1_AA_select_0_21dpi)$group.plaque.qPCR.1000.body.1000.)

GP.ord.GP.M1_AA_select_21dpi <- ordinate(GP.M1_AA_select_21dpi, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level AA raw reads 21dpe.pdf",
    width = 7.2,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_AA_select_21dpi, GP.ord.GP.M1_AA_select_21dpi, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.0.body.0.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4, alpha=0.7)  +  
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  # scale_shape_manual(values=c(shap[1:4],shap[6])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in Ae. aegypti at 21dpe")
  # stat_ellipse(type = "norm", linetype = 2,level=0.8,aes(group = Treatment), color="black")
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_CQ_select_0_14dpi, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0_14dpi))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

#################### *** 7 Phageome Cx. quinquefasciatus *** #######################
#################### * color ####################
col = plasma(9,alpha = 1, begin=0.9, end = 0, direction = 1)
#### Phyloseq Project #### 
GuadPhage_Cx_quinquefasciatus <- subset_samples(GuadPhage, mosq_species  == "Cx_quinquefasciatus" ) 
GuadPhage_Cx_quinquefasciatus
GP.M1_CQ = GuadPhage_Cx_quinquefasciatus
GP.M1_CQ
sum(sample_sums(GP.M1_CQ))

sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0. <- factor(sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0., 
                                                            levels = c("7dpi-WNV(+/)", "7dpi-WNV(-/+)", "7dpi-WNV(-/-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                                       "14dpi-WNV(+/)","14dpi-WNV(-/+)","14dpi-WNV(-/-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))
sample_data(GP.M1_CQ)$Info1 <- factor(sample_data(GP.M1_CQ)$Info1, levels = c("CQ-7dpi-WNV", "CQ-7dpi-blood", "CQ-7dpi-sucrose", "CQ-14dpi-WNV", "CQ-14dpi-blood", "CQ-14dpi-sucrose"))
sample_data(GP.M1_CQ)$Treatment <- factor(sample_data(GP.M1_CQ)$Treatment, levels = c("CQ-WNV","CQ-blood","CQ-sucrose"))
levels(sample_data(GP.M1_CQ)$Info1)
levels(sample_data(GP.M1_CQ)$group.plaque.qPCR.0.body.0.)

#### subset #### 
sample_sums(GP.M1_CQ)
sample_sums(GP.M1_CQ)[which(sample_sums(GP.M1_CQ) == 0)]
GP.M1_CQ_samselect = prune_samples(sample_sums(GP.M1_CQ) > 0 , GP.M1_CQ)
GP.M1_CQ_samselect

summary(taxa_sums(GP.M1_CQ_samselect))
GP.M1_CQ_select_0 = prune_taxa(taxa_sums(GP.M1_CQ_samselect) > 0, GP.M1_CQ_samselect)
GP.M1_CQ_select_0

#################### * 7-1) alpha diversity ####################
pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/Alpha_diversity_phage_CQ_dietv2.pdf", width = 6,height = 7,useDingbats = FALSE)
plot_richness(GP.M1_CQ_samselect, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") + 
  scale_color_manual(values = c(col[9],col[5],col[2])) + 
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.ticks.x = element_blank()) +
  ggtitle("Alpha diversity of phages in Cx. quinquefasciatus")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GP.M1_CQ_samselect, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GP.M1_CQ_samselect)$Treatment, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GP.M1_CQ_samselect)$Treatment, p.adjust.method = "BH")

#################### * 7-2) PCoA #####################
#### all ####
GP.ord.GP.M1_CQ_samselect <- ordinate(GP.M1_CQ_samselect, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level CQ raw reads ovals v2.pdf",
    width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_CQ_samselect, GP.ord.GP.M1_CQ_samselect, type = "samples", col = "Treatment")+#,label= "name") +
  theme_bw() +
  geom_point(size = 4, alpha=0.7)  + 
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  # scale_shape_manual(values=c(shap[1],shap[3],shap[5])) +
  scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in all Cx. quinquefasciatus samples")+
  stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()
set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_CQ_select_0, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

##### 7dpi #####
GP.M1_CQ_select_0_7dpi = subset_samples(GP.M1_CQ_samselect, dpi == "7dpi")
GP.M1_CQ_select_7dpi = prune_taxa(taxa_sums(GP.M1_CQ_select_0_7dpi) > 0, GP.M1_CQ_select_0_7dpi)
levels(sample_data(GP.M1_CQ_select_7dpi)$group.plaque.qPCR.0.body.0.)

GP.ord.GP.M1_CQ_select_7dpi <- ordinate(GP.M1_CQ_select_7dpi, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level CQ raw reads ovals 7dpi v2.pdf",
    width = 7.2,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_CQ_select_7dpi, GP.ord.GP.M1_CQ_select_7dpi, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.0.body.0.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4, alpha=0.7)  + 
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  scale_shape_manual(values=c(shap[1:4],shap[6])) +
  scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in Cx. quinquefasciatus at 7dpi")+
  stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment), color="black")
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_CQ_select_0_7dpi, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0_7dpi))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

##### 14dpi #####
GP.M1_CQ_select_0_14dpi = subset_samples(GP.M1_CQ_samselect, dpi == "14dpi")
GP.M1_CQ_select_14dpi = prune_taxa(taxa_sums(GP.M1_CQ_select_0_14dpi) > 0, GP.M1_CQ_select_0_14dpi)
levels(sample_data(GP.M1_CQ_select_14dpi)$group.plaque.qPCR.0.body.0.)

GP.ord.GP.M1_CQ_select_14dpi <- ordinate(GP.M1_CQ_select_14dpi, "PCoA", "bray")

pdf("/Users/shichenyan/Desktop/mosq_Guad_infection/Rfigure/PCoA of Phage contigs level CQ raw reads ovals 14dpi v2.pdf",
    width = 7.2,height = 6, useDingbats = FALSE)
plot_ordination(GP.M1_CQ_select_14dpi, GP.ord.GP.M1_CQ_select_14dpi, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.0.body.0.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4, alpha=0.7)  +  
  # geom_text(aes(label= name),size = 2, vjust = 2.5, hjust = 0) +
  scale_shape_manual(values=c(shap[1:4],shap[6])) +
  scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of phages in Cx. quinquefasciatus at 14dpi")+
  stat_ellipse(type = "norm", linetype = 2,level=0.8,aes(group = Treatment), color="black")
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GP.M1_species_select_bray <- phyloseq::distance(GP.M1_CQ_select_0_14dpi, method = "bray")
# make a data frame from the sample_data
GP.M1_species_select_sampledf <- data.frame(sample_data(GP.M1_CQ_select_0_14dpi))
# Adonis test
adonis(GP.M1_species_select_bray ~ Treatment, data = GP.M1_species_select_sampledf)
# # Call:
beta <- betadisper(GP.M1_species_select_bray, GP.M1_species_select_sampledf$Treatment)
permutest(beta)

#################### *** 8 Load 16s rRNA Phyloseq Project *** #######################
setwd("/Users/shichenyan/Desktop/mosq_Guad_infection/16s_rRNA/")
GIMBac_decont <- readRDS("GIM_Bacteriome_decont_SLV132.rds")
GIMBac_decont
col_6_Ae <-c("#8D634A", "#CCC591", "#84A0A5", "#972D15", "#C27D38", "#1D769C") 
col_6_Cq <-c("#7A7358",  "#85A789", "#C6CDF7","#29211F", "#376138","#7294D4") 

####  Library Sizes ####
df <- as.data.frame(sample_data(GIMBac_decont)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(GIMBac_decont)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
pdf("./Rfigure/total counts for each sample new.pdf", width = 8,height = 6, useDingbats = FALSE)
ggplot(data=df, aes(x=Index, y=log10(LibrarySize), color=mosq_species)) + geom_point()
dev.off()


#################### *** 9 Bacteriome Aedes aegypti  *** ####################
GIMBac.Aa <- subset_samples(GIMBac_decont, mosq_species == "Ae_aegypti")
GIMBac.Aa <- subset_taxa(GIMBac.Aa, taxa_sums(GIMBac.Aa) > 0 )
GIMBac.Aa
sort(rowSums(otu_table(GIMBac.Aa)))

col_6_Ae <-c("#8D634A", "#CCC591", "#84A0A5", "#972D15", "#C27D38", "#1D769C") 
col_6_Cq <-c("#7A7358",  "#85A789", "#C6CDF7","#29211F", "#376138","#7294D4") 
#################### * 9-1) Library Sizes #################### 
df <- as.data.frame(sample_data(GIMBac.Aa)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(GIMBac.Aa)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
head(df)

df$Info1 <- factor(df$Info1, levels = c("AA-7dpi-ZIKV", "AA-7dpi-blood", "AA-7dpi-sucrose", 
                                        "AA-21dpi-ZIKV", "AA-21dpi-blood", "AA-21dpi-sucrose"))
levels(df$Info1)

pdf("./Rfigure/total counts for Aa new.pdf", width = 9,height = 6, useDingbats = FALSE)
ggplot(data=df, aes(x=Index, y=log10(LibrarySize), color=Info1)) + geom_point() +
  facet_wrap(~ Info1)+
  scale_color_manual(values = col_6_Ae)+
  scale_y_continuous(breaks=seq(2,5,0.5))
dev.off()
#################### * 9-2) Genera proportion per sample #################### 
#### collape on Genus #### 
GIMBac.Aa.Genus <- GIMBac.Aa %>% #raw counts
  tax_glom(taxrank = "Genus")

GIMBac.Aa.Genus.per <- GIMBac.Aa.Genus  %>% 
  transform_sample_counts(function(x) {100*(x/sum(x))} ) %>%  
  psmelt()

head(GIMBac.Aa.Genus.per)
GIMBac.Aa.Genus

levels(factor(GIMBac.Aa.Genus.per$group.plaque.qPCR.1000.body.1000.))
GIMBac.Aa.Genus.per$group.plaque.qPCR.1000.body.1000. <- as.character(GIMBac.Aa.Genus.per$group.plaque.qPCR.1000.body.1000.)
GIMBac.Aa.Genus.per$group.plaque.qPCR.1000.body.1000.[which(is.na(GIMBac.Aa.Genus.per$group.plaque.qPCR.1000.body.1000.))] <- "NC"
levels(GIMBac.Aa.Genus.per$plaque_assay)

group.plaque.qPCR.1000.body.1000.order <- c("AA-7dpi-ZIKV(-/+/+)","AA-7dpi-ZIKV(-/-/+)","AA-7dpi-ZIKV(-/-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                            "AA-21dpi-ZIKV(+/+/+)","AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(-/-/+)","AA-21dpi-ZIKV(-/-/-)",
                                            "AA-21dpi-blood","AA-21dpi-sucrose","NC")

#### 把一些低丰度的anno改成"others" #### 
AaBac_agg <- aggregate(Abundance ~ mosq + plaque_assay + Genus, FUN=sum, data=GIMBac.Aa.Genus.per)
AaBac_wide <- reshape2::dcast(AaBac_agg, mosq + plaque_assay ~ Genus, value.var = 'Abundance')
AaBac_wide <- AaBac_wide[order(-AaBac_wide$g_Asaia),]
sample_order <- AaBac_wide$mosq
# head(as.data.frame(sort(colSums(AaBac_wide[,3:ncol(AaBac_wide)]) , decreasing=T)), n =10)

AaBac_name_10 <- colnames(AaBac_wide[,3:ncol(AaBac_wide)][,colSums(AaBac_wide[,3:ncol(AaBac_wide)])<10])
AaBac_name_10

# GIMBac.Aa.Genus.per$Genus <- as.character(GIMBac.Aa.Genus.per$Genus)
nc <- which(colnames(GIMBac.Aa.Genus.per) == "Genus")
nc
GIMBac.Aa.Genus.per[GIMBac.Aa.Genus.per$Genus %in% AaBac_name_10, nc] <- "Others"

GIMBac.Aa.Genus.per$Genus <- as.factor(GIMBac.Aa.Genus.per$Genus)
levels(GIMBac.Aa.Genus.per$Genus)

###--- repeat ---###
AaBac_agg <- aggregate(Abundance ~ mosq + plaque_assay + Genus, FUN=sum, data=GIMBac.Aa.Genus.per)
AaBac_wide <- reshape2::dcast(AaBac_agg, mosq + plaque_assay ~ Genus, value.var = 'Abundance')
as.data.frame(sort(colSums(AaBac_wide[,3:ncol(AaBac_wide)]) , decreasing=T))
order <- rev(rownames(as.data.frame(sort(colSums(AaBac_wide[,3:ncol(AaBac_wide)]) , decreasing=T))))
which(order == "Others" | order == "uc_f_Enterobacteriaceae")
order.1 <- c("Others", "uc_f_Enterobacteriaceae", order[-which(order == "Others" | order == "uc_f_Enterobacteriaceae")])
order.1

head(AaBac_agg)
AaBac_agg$Genus <- factor(AaBac_agg$Genus, levels = order.1)
levels(AaBac_agg$plaque_assay)
AaBac_agg$mosq <- factor(AaBac_agg$mosq, levels = sample_order)

col <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
         "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
         "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861","grey")

pdf("./Rfigure/Ae_Genera_proportion_newv2.pdf", width = 16,height = 8, useDingbats = FALSE)
ggplot(AaBac_agg, aes(x = mosq, y = Abundance, fill = Genus)) + 
  facet_grid(~ plaque_assay,scales="free",space = "free") +  #### Lay out panels in a grid
  geom_bar(stat = "identity")+ #### stat = "identity" : the height of the bar will represent the value in a column of the data frame
  scale_fill_manual(values = rev(col)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5,hjust = 1),
        axis.title.y = element_text(size = 12,vjust = 0.5),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        strip.text.x = element_text(size =8, colour = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()


#################### * 9-3) rarefy  ####################
sort(rowSums(otu_table(GIMBac.Aa)))

set.seed(100)
GIMBac.Aa_rare = rarefy_even_depth(GIMBac.Aa, sample.size = 11428, replace=FALSE, trimOTUs = F)

# rm <- c("GIM-549","GIM-319","GIM-318","GIM-320","GIM-317","GIM-30","GIM-257")
GIMBac.Aa_rm = prune_samples(sample_sums(GIMBac.Aa)<11428, GIMBac.Aa)
GIMBac.Aa_rm

GIMBac.Aa_rareall <- merge_phyloseq(GIMBac.Aa_rare, GIMBac.Aa_rm)
sort(rowSums(otu_table(GIMBac.Aa_rareall)))

#### collape on Genus #### 
GIMBac.Aa_rare.Genus <- GIMBac.Aa_rareall %>% #raw counts
  tax_glom(taxrank = "Genus")

#################### * 9-4) alpha diversity on genus  ####################
plot_richness <- function (physeq, x = "samples", color = NULL, shape = NULL,
                           title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL,
                           sortby = NULL)
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs),
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i,
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures,
      ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[,
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby,
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE,
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color,
                            shape = shape)
  p = ggplot(mdf, richness_map) + geom_boxplot(na.rm = TRUE,outlier.shape = NA)+geom_jitter(shape=19, position=position_jitter(0.25), size = 2.2, alpha=0.7)
  # if (any(!is.na(mdf[, "se"]))) {
  #   p = p + geom_errorbar(aes(ymax = value + se, ymin = value -
  #                               se), width = 0.1)
  # }
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5,
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

GIMBac.Aa.Genus
GIMBac.Aa_rare.Genus
#################### * 9-4-1) alpha diversity on genus rarefy counts ####################
GIMBac.Aa_rare.Genus
sample_data(GIMBac.Aa_rare.Genus)$Info1 <- factor(sample_data(GIMBac.Aa_rare.Genus)$Info1, levels = c( "AA-7dpi-sucrose", "AA-7dpi-blood", "AA-7dpi-ZIKV", 
                                                                                                      "AA-21dpi-sucrose","AA-21dpi-blood", "AA-21dpi-ZIKV"))
sample_data(GIMBac.Aa_rare.Genus)$Treatment <- factor(sample_data(GIMBac.Aa_rare.Genus)$Treatment, levels = c("AA-sucrose", "AA-blood", "AA-ZIKV"))
sample_data(GIMBac.Aa_rare.Genus)$dpi <- factor(sample_data(GIMBac.Aa_rare.Genus)$dpi, levels = c("7dpi", "21dpi"))
GIMBac.Aa_rare.Genus.7d = subset_samples(GIMBac.Aa_rare.Genus, dpi == "7dpi")
GIMBac.Aa_rare.Genus.21d = subset_samples(GIMBac.Aa_rare.Genus, dpi == "21dpi")

#### all ##### 
pdf("~/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Bacteriome_Aedes_alpha_diversity_rarefy_reorder.pdf", width = 6,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Aa_rare.Genus, x="Info1",measures=c("Shannon"),col="Treatment", shape = "dpi") +
  scale_color_manual(values = rev(col_6_Ae[c(4:6)])) +
  theme_bw() +
  ylim(-0.001,1.9)+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust=0.5),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Ae. aegypti on genus level")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Aa_rare.Genus, measures = c("Shannon"))
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa_rare.Genus)$Info1, p.adjust.method = "BH")


pdf("./Rfigure/Aedes_alpha_diversity_rarefy.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Aa_rare.Genus, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Ae. aegypti on genus level")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Aa_rare.Genus, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Aa_rare.Genus)$Treatment, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa_rare.Genus)$Treatment, p.adjust.method = "BH")


#### 7dpi ##### 
pdf("./Rfigure/Aedes_alpha_diversity_rarefy_7dpi.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Aa_rare.Genus.7d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Ae. aegypti at 7dpi")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Aa_rare.Genus.7d, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Aa_rare.Genus.7d)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa_rare.Genus.7d)$Info1, p.adjust.method = "BH")

#### 21dpi #####
pdf("./Rfigure/Aedes_alpha_diversity_rarefy_21dpi.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Aa_rare.Genus.21d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Ae. aegypti at 21dpi")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Aa_rare.Genus.21d, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Aa_rare.Genus.21d)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa_rare.Genus.21d)$Info1, p.adjust.method = "BH")




# # #################### alpha diversity on genus raw counts ####################
# # GIMBac.Aa.Genus
# # sample_data(GIMBac.Aa.Genus)$Info1 <- factor(sample_data(GIMBac.Aa.Genus)$Info1, levels = c("AA-7dpi-ZIKV", "AA-7dpi-blood", "AA-7dpi-sucrose", 
# #                                                                                                       "AA-21dpi-ZIKV", "AA-21dpi-blood", "AA-21dpi-sucrose"))
# # GIMBac.Aa.Genus.7d = subset_samples(GIMBac.Aa.Genus, dpi == "7dpi")
# # GIMBac.Aa.Genus.21d = subset_samples(GIMBac.Aa.Genus, dpi == "21dpi")
# # 
# # #### - 7dpi ##### 
# # pdf("./Rfigure/Aedes_alpha_diversity_7dpi_rarefy1000.pdf", width = 7,height = 8, useDingbats = FALSE)
# # plot_richness(GIMBac.Aa.Genus.7d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
# #   # scale_color_manual(values = col_6_Ae[c(1:2)]) +
# #   theme(axis.title.x = element_blank(),
# #         # axis.text.x  = element_text(size = 12, colour = "black")
# #         strip.text.x = element_text(size = 12, colour = "black"),
# #         plot.title = element_text(hjust = 0.5, size = 12)) +
# #   ggtitle("Bacteriome in Ae. aegypti at 7dpi on genus level")+
# #   ylab("Alpha Diversity Value")
# # dev.off()
# # 
# # 
# # erich <- estimate_richness(GIMBac.Aa.Genus.7d, measures = c("Chao1", "Shannon"))
# # pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Aa.Genus.7d)$Info1, p.adjust.method = "BH")
# # pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa.Genus.7d)$Info1, p.adjust.method = "BH")
# # 
# # #### - 21dpi ##### 
# # pdf("./Rfigure/Aedes_alpha_diversity_21dpi_rarefy1000.pdf", width = 8,height = 8, useDingbats = FALSE)
# # plot_richness(GIMBac.Aa.Genus.21d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
# #   scale_color_manual(values = col_6_Ae[c(4:6)]) +
# #   theme(axis.title.x = element_blank(),
# #         # axis.text.x  = element_text(size = 12, colour = "black")
# #         strip.text.x = element_text(size = 12, colour = "black"),
# #         plot.title = element_text(hjust = 0.5, size = 12)) +
# #   ggtitle("Bacteriome in Ae. aegypti at 21dpi on genus level")+
# #   ylab("Alpha Diversity Value")
# # dev.off()
# # 
# # 
# # erich <- estimate_richness(GIMBac.Aa.Genus.21d, measures = c("Chao1", "Shannon"))
# # pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Aa.Genus.21d)$Info1, p.adjust.method = "BH")
# # pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Aa.Genus.21d)$Info1, p.adjust.method = "BH")
# # 

#################### * 9-5) Beta diversity on genus rarefy counts ####################
GIMBac.Aa_rare.GenusS <- subset_samples(GIMBac.Aa_rare.Genus, group.plaque.qPCR.1000.body.1000. != "")
sample_data(GIMBac.Aa_rare.GenusS)$group.plaque.qPCR.1000.body.1000. <- factor(sample_data(GIMBac.Aa_rare.GenusS)$group.plaque.qPCR.1000.body.1000.,
                                                                  levels = c("AA-7dpi-ZIKV(-/+/+)","AA-7dpi-ZIKV(-/-/+)","AA-7dpi-ZIKV(-/-/-)","AA-7dpi-blood","AA-7dpi-sucrose",
                                                                             "AA-21dpi-ZIKV(+/+/+)","AA-21dpi-ZIKV(-/+/+)","AA-21dpi-ZIKV(-/-/+)","AA-21dpi-ZIKV(-/-/-)",
                                                                             "AA-21dpi-blood","AA-21dpi-sucrose"))

GIMBac.Aa_rare.Genus.7dS = subset_samples(GIMBac.Aa_rare.GenusS, dpi == "7dpi")
GIMBac.Aa_rare.Genus.21dS = subset_samples(GIMBac.Aa_rare.GenusS, dpi == "21dpi")

#### 7dpi ##### 
levels(sample_data(GIMBac.Aa_rare.Genus.7dS)$group.plaque.qPCR.1000.body.1000.)
levels(sample_data(GIMBac.Aa_rare.Genus.7dS)$Treatment)

GP.ord.GIMBac.Aa_rare.Genus.7d <- ordinate(GIMBac.Aa_rare.Genus.7dS, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome  AArarefy 7dpi.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Aa_rare.Genus.7dS, GP.ord.GIMBac.Aa_rare.Genus.7d, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.1000.body.1000.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  scale_shape_manual(values=rev(shap[2:6])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Ae. aegypti at 7dpi") +
  stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Aa_rare.Genus.7dS_bray <- phyloseq::distance(GIMBac.Aa_rare.Genus.7dS, method = "bray")
# make a data frame from the sample_data
GIMBac.Aa_rare.Genus.7dS_sampledf <- data.frame(sample_data(GIMBac.Aa_rare.Genus.7dS))
# Adonis test
adonis(GIMBac.Aa_rare.Genus.7dS_bray ~ Treatment, data = GIMBac.Aa_rare.Genus.7dS_sampledf)
# #                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2    2.8781 1.43907  10.953 0.35966  0.001 ***
# # Residuals 39    5.1242 0.13139         0.64034           
# # Total     41    8.0024                 1.00000   
# Call:
beta <- betadisper(GIMBac.Aa_rare.Genus.7dS_bray, GIMBac.Aa_rare.Genus.7dS_sampledf$Treatment)
permutest(beta)


#### 21dpi ##### 
levels(sample_data(GIMBac.Aa_rare.Genus.21dS)$group.plaque.qPCR.1000.body.1000.)
levels(sample_data(GIMBac.Aa_rare.Genus.21dS)$Treatment)

GP.ord.GIMBac.Aa_rare.Genus.21d <- ordinate(GIMBac.Aa_rare.Genus.21dS, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome  AArarefy 21dpi.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Aa_rare.Genus.21dS, GP.ord.GIMBac.Aa_rare.Genus.21d, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.1000.body.1000.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  scale_shape_manual(values= c(shap[1:4],shap[6],shap[5])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Ae. aegypti at 21dpi")
  # stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Aa_rare.Genus.21dS_bray <- phyloseq::distance(GIMBac.Aa_rare.Genus.21dS, method = "bray")
# make a data frame from the sample_data
GIMBac.Aa_rare.Genus.21dS_sampledf <- data.frame(sample_data(GIMBac.Aa_rare.Genus.21dS))
# Adonis test
adonis(GIMBac.Aa_rare.Genus.21dS_bray ~ Treatment, data = GIMBac.Aa_rare.Genus.21dS_sampledf)
# #                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2   0.69415 0.34707  9.1003 0.27915  0.001 ***
# # Residuals 47   1.79252 0.03814         0.72085           
# # Total     49   2.48667                 1.00000   
# Call:
beta <- betadisper(GIMBac.Aa_rare.Genus.21dS_bray, GIMBac.Aa_rare.Genus.21dS_sampledf$Treatment)
permutest(beta)


#### all ####
GP.ord.GIMBac.Aa_rare.Genus <- ordinate(GIMBac.Aa_rare.Genus, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome AArarefy.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Aa_rare.Genus, GP.ord.GIMBac.Aa_rare.Genus, type = "samples", color = "Treatment")+
  # shape = "group.plaque.qPCR.0.body.0.")+,label= "name") +
  theme_bw()+
  geom_point(size = 4,  pch = 16, alpha=0.7)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  # scale_shape_manual(values=shap) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Ae. aegypti level")
# stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Aa_rare.Genus_bray <- phyloseq::distance(GIMBac.Aa_rare.Genus, method = "bray")
# make a data frame from the sample_data
GIMBac.Aa_rare.Genus_sampledf <- data.frame(sample_data(GIMBac.Aa_rare.Genus))
# Adonis test
adonis(GIMBac.Aa_rare.Genus_bray ~ Treatment, data = GIMBac.Aa_rare.Genus_sampledf)
# #         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2    0.5371 0.26853  1.4751 0.0985  0.193
# # Residuals 27    4.9153 0.18205         0.9015       
# # Total     29    5.4524                 1.0000  
# Call:
beta <- betadisper(GIMBac.Aa_rare.Genus_bray, GIMBac.Aa_rare.Genus_sampledf$Treatment)
permutest(beta)


#################### *** 10 Bacteriome Culex quinquefasciatus *** ####################
GIMBac.Cq <- subset_samples(GIMBac_decont, mosq_species == "Cx_quinquefasciatus")
GIMBac.Cq <- subset_taxa(GIMBac.Cq, taxa_sums(GIMBac.Cq) > 0 )
GIMBac.Cq
sample_sums(GIMBac.Cq)
sort(sample_sums(GIMBac.Cq))

sample_data(GIMBac.Cq)$group.plaque.qPCR.0.body.0. <- factor(sample_data(GIMBac.Cq)$group.plaque.qPCR.0.body.0., 
                                                             levels = c("CQ-7dpi-sucrose","CQ-7dpi-blood","7dpi-WNV(-/-)", "7dpi-WNV(-/+)", "7dpi-WNV(+/)", 
                                                                        "CQ-14dpi-sucrose","CQ-14dpi-blood","14dpi-WNV(-/-)","14dpi-WNV(-/+)","14dpi-WNV(+/)"))
sample_data(GIMBac.Cq)$Info1 <- factor(sample_data(GIMBac.Cq)$Info1, levels = c("CQ-7dpi-sucrose", "CQ-7dpi-blood", "CQ-7dpi-WNV", "CQ-14dpi-sucrose", "CQ-14dpi-blood", "CQ-14dpi-WNV"))
sample_data(GIMBac.Cq)$Treatment <- factor(sample_data(GIMBac.Cq)$Treatment, levels = c("CQ-sucrose", "CQ-blood", "CQ-WNV"))
#################### * 10-1) Library Sizes #################### 
df <- as.data.frame(sample_data(GIMBac.Cq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(GIMBac.Cq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
head(df)

df$Info1 <- factor(df$Info1, levels = c("CQ-7dpi-WNV", "CQ-7dpi-blood", "CQ-7dpi-sucrose", 
                                        "CQ-14dpi-WNV", "CQ-14dpi-blood", "CQ-14dpi-sucrose"))

pdf("./Rfigure/total counts for Cq new.pdf", width = 9,height = 6, useDingbats = FALSE)
ggplot(data=df, aes(x=Index, y=log10(LibrarySize), color=Info1)) + geom_point() +
  facet_wrap(~ Info1) +
  scale_color_manual(values = col_6_Cq)+
  scale_y_continuous(breaks=seq(2,5,0.5))
dev.off()

#################### * 10-2) Genera proportion per sample ####################
#### collape on Genus #### 
GIMBac.Cq.Genus <- GIMBac.Cq %>% #raw counts
  tax_glom(taxrank = "Genus")

GIMBac.Cq.Genus.per <- GIMBac.Cq.Genus  %>% 
  transform_sample_counts(function(x) {100*(x/sum(x))} ) %>%  
  psmelt()

colnames(GIMBac.Cq.Genus.per)
levels(GIMBac.Cq.Genus.per$plaque_assay)
GIMBac.Cq.Genus.per$plaque_assay <- factor(GIMBac.Cq.Genus.per$plaque_assay, 
                                           levels = c("7dpi-WNV(+)", "7dpi-WNV(-)","CQ-7dpi-blood","CQ-7dpi-sucrose",
                                                      "14dpi-WNV(+)","14dpi-WNV(-)","CQ-14dpi-blood","CQ-14dpi-sucrose"))
levels(GIMBac.Cq.Genus.per$plaque_assay)
#### 把一些低丰度的anno改成“others”  ####
CqBac_agg <- aggregate(Abundance ~ mosq + plaque_assay + Genus, FUN=sum, data=GIMBac.Cq.Genus.per)
CqBac_wide <- dcast(CqBac_agg, mosq + plaque_assay ~ Genus, value.var = 'Abundance')
CqBac_wide_ordered <- CqBac_wide[with(CqBac_wide, order(plaque_assay,-g_Wolbachia,-g_Serratia)), ]
sample_order <- CqBac_wide_ordered$mosq

CqBac_name_10 <- colnames(CqBac_wide[,3:ncol(CqBac_wide)][,colSums(CqBac_wide[,3:ncol(CqBac_wide)])<10])
CqBac_name_10

# GIMBac.Cq.Genus.per$Genus <- as.character(GIMBac.Cq.Genus.per$Genus)
nc <- which(colnames(GIMBac.Cq.Genus.per) == "Genus")
nc
GIMBac.Cq.Genus.per[GIMBac.Cq.Genus.per$Genus %in% CqBac_name_10, nc] <- "Others"

GIMBac.Cq.Genus.per$Genus <- as.factor(GIMBac.Cq.Genus.per$Genus)
levels(GIMBac.Cq.Genus.per$Genus)

###--- repeat ---### 
CqBac_agg <- aggregate(Abundance ~ mosq + plaque_assay + Genus, FUN=sum, data=GIMBac.Cq.Genus.per)
CqBac_wide <- dcast(CqBac_agg, mosq + plaque_assay ~ Genus, value.var = 'Abundance')
as.data.frame(sort(colSums(CqBac_wide[,3:ncol(CqBac_wide)]) , decreasing=T))
order <- rev(rownames(as.data.frame(sort(colSums(CqBac_wide[,3:ncol(CqBac_wide)]) , decreasing=T))))
which(order == "Others")
order.1 <- c("Others", order[-which(order == "Others")])
order.1

CqBac_agg$mosq <- factor(CqBac_agg$mosq, levels = sample_order)
CqBac_agg$Genus <- factor(CqBac_agg$Genus, levels = order.1)
levels(CqBac_agg$plaque_assay)

#### plot Genera proportion bar plot ####
col <- c("#5F7FC7", "#AD6F3B", "#CBD588", "#673770", "orange", "#CD9BCD", "#652926", "#C84248", "grey")
pdf("./Rfigure/Cq_Genera_proportion newv2.pdf", width = 13,height = 8, useDingbats = FALSE)
ggplot(CqBac_agg, aes(x = mosq, y = Abundance, fill = Genus)) + 
  facet_grid(~ plaque_assay,scales="free",space = "free") +  #### Lay out panels in a grid
  geom_bar(stat = "identity") + #### stat = "identity" : the height of the bar will represent the value in a column of the data frame
  scale_fill_manual(values = rev(col)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5,hjust = 1),
        axis.title.y = element_text(size = 12,vjust = 0.5),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        strip.text.x = element_text(size =8, colour = "black"),
        plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))
dev.off()

#### boxplot ####
CqBac_agg$dpi <- CqBac_agg$plaque_assay
CqBac_agg$dpi <- gsub("CQ-", "", CqBac_agg$dpi)

CqBac_agg_v <- CqBac_agg %>% separate(dpi, c("dpi", "food"), sep ="-"); head(CqBac_agg_v)
# Serratia Asaia Wolbachia
CqBac_agg_v$Genus <- gsub("g_", "", CqBac_agg_v$Genus); head(CqBac_agg_v)
CqBac_agg_v <- CqBac_agg_v[CqBac_agg_v$Genus %in% c("Asaia", "Wolbachia"), ] #"Serratia",
CqBac_agg_v$dpi <- factor(CqBac_agg_v$dpi, levels= c("7dpi", "14dpi"))


wilcox.test(CqBac_agg_v[CqBac_agg_v$Genus=="Serratia"&CqBac_agg_v$dpi=="7dpi",]$Abundance,
                     CqBac_agg_v[CqBac_agg_v$Genus=="Serratia"&CqBac_agg_v$dpi=="14dpi",]$Abundance, 
                     p.adjust.method = "BH")

wilcox.test(CqBac_agg_v[CqBac_agg_v$Genus=="Asaia"&CqBac_agg_v$dpi=="7dpi",]$Abundance,
            CqBac_agg_v[CqBac_agg_v$Genus=="Asaia"&CqBac_agg_v$dpi=="14dpi",]$Abundance, 
            p.adjust.method = "BH")

wilcox.test(CqBac_agg_v[CqBac_agg_v$Genus=="Wolbachia"&CqBac_agg_v$dpi=="7dpi",]$Abundance,
            CqBac_agg_v[CqBac_agg_v$Genus=="Wolbachia"&CqBac_agg_v$dpi=="14dpi",]$Abundance, 
            p.adjust.method = "BH")

CqBac_agg_v$Food_sources[CqBac_agg_v$food == "WNV(+)"] <- "WNV(plaque+)"
CqBac_agg_v$Food_sources[CqBac_agg_v$food == "WNV("] <- "WNV(plaque-)"
CqBac_agg_v$Food_sources[CqBac_agg_v$food == "sucrose"] <- "sucrose"
CqBac_agg_v$Food_sources[CqBac_agg_v$food == "blood"] <- "blood"
CqBac_agg_v$Food_sources <- factor(CqBac_agg_v$Food_sources, levels=c("sucrose", "blood", "WNV(plaque-)", "WNV(plaque+)"))

pdf("~/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Bacteriome_Cq_Genera_proportion_sepcific_v2.pdf", width = 12,height = 6, useDingbats = FALSE)
ggplot(CqBac_agg_v, aes(x=Genus, y=Abundance,fill=dpi)) + 
  geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
  # geom_jitter(size = 0.5, position=position_jitter(0.1))+
  facet_wrap(~Food_sources, nrow = 1)+
  # geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1)+
  # scale_shape_manual(values = c(9,19))+
  theme_bw()+
  scale_fill_manual(values = col_6_Ae[c(3,1)])+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12,vjust = 0.5),
        # plot.margin = unit(c(1, 1, 1, 1), "cm"),
        strip.text.x = element_text(size =12, colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  ylab("Proportion")+
  stat_compare_means()
dev.off()

#################### * 10-3) rarefy  ####################
sort(rowSums(otu_table(GIMBac.Cq)))

set.seed(100)
GIMBac.Cq_rare = rarefy_even_depth(GIMBac.Cq, sample.size = 24168, replace=FALSE, trimOTUs = F)

# rm <- c("GIM-549","GIM-319","GIM-318","GIM-320","GIM-317","GIM-30","GIM-257")
GIMBac.Cq_rm = prune_samples(sample_sums(GIMBac.Cq)<24168, GIMBac.Cq)
GIMBac.Cq_rm

GIMBac.Cq_rareall <- merge_phyloseq(GIMBac.Cq_rare, GIMBac.Cq_rm)
sort(rowSums(otu_table(GIMBac.Cq_rareall)))

#### collape on Genus #### 
GIMBac.Cq_rare.Genus <- GIMBac.Cq_rareall %>% #raw counts
  tax_glom(taxrank = "Genus")

#################### * 10-4) alpha diversity on genus  ####################
plot_richness <- function (physeq, x = "samples", color = NULL, shape = NULL,
                           title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL,
                           sortby = NULL)
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs),
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i,
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures,
      ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[,
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby,
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE,
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color,
                            shape = shape)
  p = ggplot(mdf, richness_map) + geom_boxplot(na.rm = TRUE,outlier.shape = NA)+geom_jitter(shape=19, position=position_jitter(0.25), size = 2.2, alpha=0.7)
  # if (any(!is.na(mdf[, "se"]))) {
  #   p = p + geom_errorbar(aes(ymax = value + se, ymin = value -
  #                               se), width = 0.1)
  # }
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5,
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

GIMBac.Cq.Genus
GIMBac.Cq_rare.Genus
#################### * 10-4-1) alpha diversity on genus rarefy counts ####################
GIMBac.Cq_rare.Genus
sample_data(GIMBac.Cq_rare.Genus)$Info1
sample_data(GIMBac.Cq_rare.Genus)$Treatment
sample_data(GIMBac.Cq_rare.Genus)$group.plaque.qPCR.0.body.0.
GIMBac.Cq_rare.Genus.7d = subset_samples(GIMBac.Cq_rare.Genus, dpi == "7dpi")
GIMBac.Cq_rare.Genus.14d = subset_samples(GIMBac.Cq_rare.Genus, dpi == "14dpi")
sample_data(GIMBac.Cq_rare.Genus)$dpi <- factor(sample_data(GIMBac.Cq_rare.Genus)$dpi, levels = c("7dpi", "14dpi"))
#### all ##### 
pdf("~/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Bacteriome_Culex_alpha_diversity_rarefy_reorder.pdf", width = 6,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Cq_rare.Genus, x="Info1",measures=c("Shannon"),col="Treatment", shape = "dpi") +
  scale_color_manual(values = rev(col_6_Ae[c(4:6)])) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust=0.5),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Cx. quinquefasciatus on genus level")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Cq_rare.Genus, measures = c("Shannon"))
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Cq_rare.Genus)$Info1, p.adjust.method = "BH")

pdf("./Rfigure/Culex_alpha_diversity_rarefy.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Cq_rare.Genus, x="Treatment",measures=c("Chao1", "Shannon"),col="Treatment") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Cx. quinquefasciatus on genus level")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Cq_rare.Genus, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Cq_rare.Genus)$Treatment, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Cq_rare.Genus)$Treatment, p.adjust.method = "BH")


#### 7dpi ##### 
pdf("./Rfigure/Culex_alpha_diversity_rarefy_7dpi.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Cq_rare.Genus.7d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Cx. quinquefasciatus at 7dpi")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Cq_rare.Genus.7d, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Cq_rare.Genus.7d)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Cq_rare.Genus.7d)$Info1, p.adjust.method = "BH")

#### 14dpi #####
pdf("./Rfigure/Culex_alpha_diversity_rarefy_14dpi.pdf", width = 5,height = 7, useDingbats = FALSE)
plot_richness(GIMBac.Cq_rare.Genus.14d, x="Info1",measures=c("Chao1", "Shannon"),col="Info1") +
  scale_color_manual(values = col_6_Ae[c(4:6)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Alpha diversity of bacteriome in Cx. quinquefasciatus at 14dpi")+
  ylab("Alpha Diversity Value")
dev.off()

erich <- estimate_richness(GIMBac.Cq_rare.Genus.14d, measures = c("Chao1", "Shannon"))
pairwise.wilcox.test(erich$Chao1,sample_data(GIMBac.Cq_rare.Genus.14d)$Info1, p.adjust.method = "BH")
pairwise.wilcox.test(erich$Shannon,sample_data(GIMBac.Cq_rare.Genus.14d)$Info1, p.adjust.method = "BH")




#################### * 10-5) Beta diversity on genus rarefy counts ####################
GIMBac.Cq_rare.Genus
GIMBac.Cq_rare.Genus.7d
GIMBac.Cq_rare.Genus.14d
#### all ##### 
GP.ord.GIMBac.Cq_rare.Genus <- ordinate(GIMBac.Cq_rare.Genus, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome CQrarefy.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Cq_rare.Genus, GP.ord.GIMBac.Cq_rare.Genus, type = "samples", color = "Treatment")+
                # shape = "group.plaque.qPCR.0.body.0.")+,label= "name") +
  theme_bw()+
  geom_point(size = 4,  pch = 16)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  # scale_shape_manual(values=shap) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Cx. quinquefasciatus on genus level")
# stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Cq_rare.Genus.7d_bray <- phyloseq::distance(GIMBac.Cq_rare.Genus.7d, method = "bray")
# make a data frame from the sample_data
GIMBac.Cq_rare.Genus.7d_sampledf <- data.frame(sample_data(GIMBac.Cq_rare.Genus.7d))
# Adonis test
adonis(GIMBac.Cq_rare.Genus.7d_bray ~ Treatment, data = GIMBac.Cq_rare.Genus.7d_sampledf)
# #         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2    0.5371 0.26853  1.4751 0.0985  0.193
# # Residuals 27    4.9153 0.18205         0.9015       
# # Total     29    5.4524                 1.0000  
# Call:
beta <- betadisper(GIMBac.Cq_rare.Genus.7d_bray, GIMBac.Cq_rare.Genus.7d_sampledf$Treatment)
permutest(beta)


#### 7dpi ##### 
GP.ord.GIMBac.Cq_rare.Genus.7d <- ordinate(GIMBac.Cq_rare.Genus.7d, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome CQrarefy 7dpi.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Cq_rare.Genus.7d, GP.ord.GIMBac.Cq_rare.Genus.7d, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.0.body.0.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  scale_shape_manual(values=c(shap[1:4],shap[6])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Cx. quinquefasciatus at 7dpi")
  # stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Cq_rare.Genus.7d_bray <- phyloseq::distance(GIMBac.Cq_rare.Genus.7d, method = "bray")
# make a data frame from the sample_data
GIMBac.Cq_rare.Genus.7d_sampledf <- data.frame(sample_data(GIMBac.Cq_rare.Genus.7d))
# Adonis test
adonis(GIMBac.Cq_rare.Genus.7d_bray ~ Treatment, data = GIMBac.Cq_rare.Genus.7d_sampledf)
# #         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2    0.5371 0.26853  1.4751 0.0985  0.193
# # Residuals 27    4.9153 0.18205         0.9015       
# # Total     29    5.4524                 1.0000  
# Call:
beta <- betadisper(GIMBac.Cq_rare.Genus.7d_bray, GIMBac.Cq_rare.Genus.7d_sampledf$Treatment)
permutest(beta)


#### 14dpi ##### 
GP.ord.GIMBac.Cq_rare.Genus.14d <- ordinate(GIMBac.Cq_rare.Genus.14d, "PCoA", "bray")

pdf("./Rfigure/PCoA of bacteriome CQrarefy 14dpi.pdf", width = 7,height = 6, useDingbats = FALSE)
plot_ordination(GIMBac.Cq_rare.Genus.14d, GP.ord.GIMBac.Cq_rare.Genus.14d, type = "samples", color = "Treatment", shape = "group.plaque.qPCR.0.body.0.")+#,label= "name") +
  theme_bw()+
  geom_point(size = 4)  + 
  # geom_text(aes(label= mosq),size = 2, vjust = 0, hjust = -0.5) +
  scale_color_manual(values = col_6_Ae[c(4:6)])+
  scale_shape_manual(values=c(shap[1:4],shap[6])) +
  # scale_color_manual(values = c(col[9],col[5],col[2]))+
  theme(legend.title=element_blank(),
        # legend.text= element_text(),
        legend.position = "right",
        plot.title = element_text(size = 12,hjust = 0.5))+
  ggtitle("PCoA of bacteriome in Cx. quinquefasciatus at 14dpi")
# stat_ellipse(type = "norm", linetype = 2,level=0.6,aes(group = Treatment))
dev.off()

set.seed(1)
# Calculate bray curtis distance matrix
GIMBac.Cq_rare.Genus.14d_bray <- phyloseq::distance(GIMBac.Cq_rare.Genus.14d, method = "bray")
# make a data frame from the sample_data
GIMBac.Cq_rare.Genus.14d_sampledf <- data.frame(sample_data(GIMBac.Cq_rare.Genus.14d))
# Adonis test
adonis(GIMBac.Cq_rare.Genus.14d_bray ~ Treatment, data = GIMBac.Cq_rare.Genus.14d_sampledf)
# #         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Treatment  2    0.4523 0.22613  1.1465 0.07828  0.323
# # Residuals 27    5.3254 0.19724         0.92172       
# # Total     29    5.7777                 1.00000 
# Call:
beta <- betadisper(GIMBac.Cq_rare.Genus.14d_bray, GIMBac.Cq_rare.Genus.14d_sampledf$Treatment)
permutest(beta)


#################### *** 11 DESeq2 differential bacteria *** #######################
#################### * 11-1）Aedes * ####################
df_otu <- as.data.frame(t(GIMBac.Aa_rare.Genus@otu_table)); head(df_otu)
df_tax <- as.data.frame(GIMBac.Aa_rare.Genus@tax_table); head(df_tax)
df_meta <- data.frame(sample_data(GIMBac.Aa_rare.Genus))
df <- merge(df_otu, df_tax, by = "row.names", all=T)
head(df)
#################### * 11-1-1) Aedes-Diet * ####################
AAdiet <- readRDS("../RDS/bacteria_DESeq2_AA_treatment.RDS")
AAdiet
gp <- as.data.frame(sort(table(c(AAdiet[[1]]$Genus, AAdiet[[2]]$Genus, AAdiet[[3]]$Genus)), decreasing = T))
gp <- gp[-c(8,10,12),]
order <- c("Chryseobacterium", "Novosphingobium", "Acidovorax", "Cupriavidus", "Roseococcus",
           "Acinetobacter", "Klebsiella", "Brevundimonas", "Sphingobacterium") #"Pseudomonas",

df_sp <- df[df$Genus %in% gp$Var1,]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "Genus")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T); head(df_merge)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0
df_merge$Genus <- gsub("g_", "", df_merge$Genus)
df_merge$Genus <- factor(df_merge$Genus, levels = order); levels(df_merge$Genus)

df_merge_aa <- df_merge
df_merge_aa <- df_merge_aa[,colnames(df_merge_aa) %in% c("Treatment", "Genus", "variable", "log10reads")]

####  plot  ####
pdf("./Rfigure/Bacteria_AA_treatment_DESeq.pdf", width = 12, height = 10,useDingbats = FALSE)
ggplot(df_merge, aes(x = Treatment, y = log10reads, fill = Treatment)) + 
  facet_wrap(~Genus,nrow =2, scales="fixed") + 
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = rev(col_6_Ae[c(4:6)]))+
  theme(axis.title.y = element_text(size = 14, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        strip.text.x =element_text(size = 14, hjust = 0.5)) +
  ylab("log10 rarefied counts")
dev.off()

#################### * 11-1-2) Aedes-Infectious * ####################
AAinfec <- readRDS("../RDS/Bacteria_DESeq2_AA_infectious.RDS")
AAinfec

df_sp <- df[df$Genus == "g_Acinetobacter",]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "Genus")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T); head(df_merge)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0
df_merge$Genus <- gsub("g_", "", df_merge$Genus)
# df_merge$Genus <- factor(df_merge$Genus, levels = order); levels(df_merge$Genus)

df_merge <- df_merge[df_merge$head_ZIKV.1000 %in% c("AA-7dpi-ZIKV(-)","AA-7dpi-ZIKV(+)", "AA-21dpi-ZIKV(-)", "AA-21dpi-ZIKV(+)"),]
df_merge$head_ZIKV.1000 <- factor(df_merge$head_ZIKV.1000, levels = c("AA-7dpi-ZIKV(+)","AA-7dpi-ZIKV(-)", "AA-21dpi-ZIKV(+)", "AA-21dpi-ZIKV(-)"))

####  plot ####  
# pdf("/Rfigure/AA_infectious_DESeq.pdf", width = 4.5, height = 6,useDingbats = FALSE)
ggplot(df_merge, aes(x = plaque_assay, y = log10reads, fill = plaque_assay)) + 
  facet_wrap(~Genus+dpi,nrow =1, scales="fixed") +
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(col[9], col[3], col[8], col[2]))+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, hjust = 0.5),
        legend.position = "bottom",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 raw reads")
dev.off()

#################### * 11-2）Culex * ####################
df_otu <- as.data.frame(t(GIMBac.Cq_rare.Genus@otu_table)); head(df_otu)
df_tax <- as.data.frame(GIMBac.Cq_rare.Genus@tax_table); head(df_tax)
df_meta <- data.frame(sample_data(GIMBac.Cq_rare.Genus))
df <- merge(df_otu, df_tax, by = "row.names", all=T)
head(df)
#################### * 11-2-1) Culex-Diet * ####
CQdiet <- readRDS("../RDS/Bacteria_DESeq2_CQ_treatment.RDS")
CQdiet

gp <- as.data.frame(sort(table(c(CQdiet[[1]]$Genus, CQdiet[[2]]$Genus, CQdiet[[3]]$Genus)), decreasing = T))
gp
order <- c("Serratia", "Pseudomonas", "Bacillus")

df_sp <- df[df$Genus %in% c("g_Serratia", "g_Pseudomonas", "g_Bacillus"),]
df_sp <- df_sp[, colnames(df_sp) %in% c(colnames(df_otu), "Genus")]
df_spm <- melt(df_sp); head(df_spm)
df_meta$variable <- rownames(df_meta)
df_merge <- merge(df_spm, df_meta, by = "variable", all =T); head(df_merge)

df_merge$log10reads <- log10(df_merge$value)
is.na(df_merge$log10reads) <- sapply(df_merge$log10reads, is.infinite)
df_merge$log10reads[is.na(df_merge$log10reads)] <- 0
df_merge$Genus <- gsub("g_", "", df_merge$Genus)
df_merge$Genus <- factor(df_merge$Genus, levels = order); levels(df_merge$Genus)

df_merge_cq <- df_merge
df_merge_cq <- df_merge_cq[,colnames(df_merge_cq) %in% c("Treatment", "Genus", "variable", "log10reads")]


####  plot  ####
pdf("./Rfigure/Bacteria_CQ_treatment_DESeq.pdf", width = 3.5, height = 10,useDingbats = FALSE)
ggplot(df_merge, aes(x = Treatment, y = log10reads, fill = Treatment)) + 
  facet_wrap(~Genus,nrow =2, scales="fixed") + 
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = col_6_Ae[c(4:6)])+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        strip.text.x =element_text(size = 10, hjust = 0.5)) +
  ylab("log10 rarefied counts")+
  stat_compare_means()
dev.off()

#################### * 11-1）Aedes & Culex * ####################
df_merge_all <- rbind(df_merge_aa, df_merge_cq)
df_merge_all$Genus
df_merge_all$Treatment
df_merge_all$Treatment <- factor(df_merge_all$Treatment , levels = c("AA-sucrose", "AA-blood", "AA-ZIKV", "CQ-sucrose", "CQ-blood", "CQ-WNV"))

####  plot  ####
pdf("~/OneDrive/1. Project + Guadeloupe mosq/Infection/Rfigure/Bacteria_AA+CQ_treatment_DESeq_v2.pdf", width = 8.5, height = 12,useDingbats = FALSE)
ggplot(df_merge_all, aes(x = Treatment, y = log10reads, fill = Treatment)) + 
  facet_wrap(~Genus,ncol = 3, scales="free_x") + 
  theme_bw() +
  geom_boxplot(na.rm = TRUE,outlier.shape = NA, alpha = 0.7)+
  geom_jitter(shape=19, position=position_jitter(0.2), size = 2,alpha = 0.5)+
  scale_fill_manual(values = c(rev(col_6_Ae[c(4:6)]), rev(col_6_Ae[c(4:6)])))+
  theme(axis.title.y = element_text(size = 14, vjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right",
        strip.text.x =element_text(size = 14, hjust = 0.5)) +
  ylab("log10 rarefied counts")
  # stat_compare_means()
dev.off()

#################### *** 12 Phage correlate with bacteria *** #######################
####  phage  ####
Phage_host = readRDS("../RDS/GIM_phage_host.RDS")
Phage_host@sam_data
Phage_host_subset <- subset_samples(Phage_host, host_genus != "no result")
Phage_host_subset@sam_data
pOTU <- as.data.frame(t(Phage_host_subset@otu_table)); head(pOTU)
pTAX <- data.frame(sample_data(Phage_host_subset))
pTAX <- pTAX[,15:16]
pTAX$phageID <- rownames(pTAX)
pOTU$phageID <- rownames(pOTU)
phage_merge <- merge(pOTU, pTAX, by = "phageID"); head(phage_merge)

####  bacteria  ####
Bac = readRDS("GIM_Bacteriome_decont_SLV132.rds")
BacGenus <- Bac %>% #raw counts
  tax_glom(taxrank = "Genus")

BacOTU <- as.data.frame(t(BacGenus@otu_table)); BacOTU[1:5,1:5]
colnames(BacOTU) <- paste("Bac_", colnames(BacOTU), sep = ""); BacOTU[1:5,1:5]
BacTAX <- as.data.frame(BacGenus@tax_table); BacTAX[1:5,1:6]
BacOTU$bacID <- rownames(BacOTU); head(BacOTU)
BacTAX$bacID <- rownames(BacTAX); head(BacTAX)
Bac_merge <- merge(BacOTU, BacTAX[,5:7], by = "bacID"); head(Bac_merge)
Bac_merge$Genus <- gsub("g_", "", Bac_merge$Genus)

####  merge phage + bacteria  ####
head(phage_merge)
head(Bac_merge)
colnames(phage_merge) <- gsub("host_genus", "Genus", colnames(phage_merge))

table(phage_merge$Genus)
levels(factor(phage_merge$Genus))
levels(factor(Bac_merge$Genus))

# c(phage_merge$Genus, "uc_f_Enterobacteriaceae", 
Bac_merge_slec <- Bac_merge[Bac_merge$Genus %in% c(phage_merge$Genus, "uc_f_Enterobacteriaceae"),]
phage_merge$Genus_og <- phage_merge$Genus
phage_merge$Genus <- gsub("Enterobacter", "uc_f_Enterobacteriaceae", phage_merge$Genus)
phage_merge$Genus <- gsub("Phytobacter", "uc_f_Enterobacteriaceae", phage_merge$Genus)

levels(factor(phage_merge$Genus))
levels(factor(Bac_merge_slec$Genus))

phage_Bac <- merge(phage_merge, Bac_merge_slec, all = T, by ="Genus") 
colnames(phage_Bac)
phage_Bac <- phage_Bac[, !colnames(phage_Bac) %in% c("Family", "host_species", "bacID", "Bac_NC-p07", "Bac_NC-p06", "Bac_GIM-NC1", "Bac_GIM-NC2")]
colnames(phage_Bac)
phage_Bac_melt <- melt(phage_Bac)
head(phage_Bac_melt)
phage_Bac_melt <- phage_Bac_melt %>% separate(variable, c("cat","Sample"), sep ="GIM"); head(phage_Bac_melt)
levels(factor(phage_Bac_melt$cat))
levels(factor(phage_Bac_melt$Sample))

phage_Bac_melt$Sample <- gsub("-", "", phage_Bac_melt$Sample)
phage_Bac_melt$Sample <- paste("GIM", phage_Bac_melt$Sample, sep = "")
phage_Bac_melt$cat <- ifelse(phage_Bac_melt$cat == "", "phage", "bacteria"); head(phage_Bac_melt)

phage_Bac_melt_wide <- dcast(phage_Bac_melt, phageID + Sample + Genus + Genus_og ~ cat, value.var = 'value')
head(phage_Bac_melt_wide)

####  add meta  ####
head(phage_Bac_melt)
phage_Bac_melt_wide$Genus <- as.factor(phage_Bac_melt_wide$Genus)

phage_Bac_melt_wide$log10bacteria <- log10(phage_Bac_melt_wide$bacteria)
is.na(phage_Bac_melt_wide$log10bacteria) <- sapply(phage_Bac_melt_wide$log10bacteria, is.infinite)
phage_Bac_melt_wide$log10bacteria[is.na(phage_Bac_melt_wide$log10bacteria)]<-0

phage_Bac_melt_wide$log10phage <- log10(phage_Bac_melt_wide$phage)
is.na(phage_Bac_melt_wide$log10phage) <- sapply(phage_Bac_melt_wide$log10phage, is.infinite)
phage_Bac_melt_wide$log10phage[is.na(phage_Bac_melt_wide$log10phage)]<-0
phage_Bac_melt_wide

samplemeta <- data.frame(BacGenus@sam_data); colnames(samplemeta)
samplemeta <- samplemeta[ ,colnames(samplemeta) %in% c("mosq", "Index", "mosq_species", "Treatment", "dpi", "plaque_assay", "Info1")]
samplemeta
samplemeta$Sample <- gsub("-", "", rownames(samplemeta)); samplemeta

phage_Bac_melt_wide_meta <- merge(phage_Bac_melt_wide, samplemeta, by= "Sample")

####  plot  ####
levels(phage_Bac_melt_wide_meta$Genus)
phage_Bac_melt_wide_metaSle <- phage_Bac_melt_wide_meta[phage_Bac_melt_wide_meta$Genus %in% c("Wolbachia", "Serratia"),]

phage_Bac_melt_wide_metaSle$Genus <- factor(phage_Bac_melt_wide_metaSle$Genus, levels = c("Wolbachia", "Serratia"))
phage_Bac_melt_wide_metaSle$mosq_species <- factor(phage_Bac_melt_wide_metaSle$mosq_species, levels = c("Ae_aegypti", "Cx_quinquefasciatus"))

pdf("./Rfigure/Phage_vs_WolbachiaSerratia_abundance_distribution.pdf", width = 8,height =5,useDingbats = FALSE)
ggscatter(phage_Bac_melt_wide_metaSle, x = "log10bacteria", y = "log10phage",
          color = "mosq_species", palette = c("#7294D4", "#C27D38"), size = 2.2, alpha = 0.8,
          conf.int = T, cor.coef = TRUE, cor.method = "pearson")+
  facet_wrap(~Genus+dpi+mosq_species)+
  geom_smooth(color="black")+
  theme_bw()+
  theme(axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_text(size = 12, hjust=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        strip.text.x =element_text(size = 14, hjust = 0.5)) +
  xlab("log10 host bacteria reads")+
  ylab("1og10 phage reads")
dev.off()

Wolbachia <- phage_Bac_melt_wide_metaSle[phage_Bac_melt_wide_metaSle$Genus=="Wolbachia",]
Serratia <- phage_Bac_melt_wide_metaSle[phage_Bac_melt_wide_metaSle$Genus=="Serratia",]
cor.test(Wolbachia$bacteria, Wolbachia$phage, method = "pearson")
cor.test(Serratia$bacteria, Serratia$phage, method = "pearson")
cor.test

Wolbachia_o <- Wolbachia[order(-Wolbachia$log10bacteria),]
Wolbachia_m <- melt(Wolbachia[,-5:-6]); head(Wolbachia_m)
Wolbachia_m$variable
Wolbachia_m$Sample <- factor(Wolbachia_m$Sample, levels = unique(Wolbachia_o$Sample))

phage_Bac_melt_wide_metaSle_Culex <- phage_Bac_melt_wide_metaSle[phage_Bac_melt_wide_metaSle$mosq_species == "Cx_quinquefasciatus",]
phage_Bac_melt_wide_metaSle_Culex$dpi <- factor(phage_Bac_melt_wide_metaSle_Culex$dpi, levels = c("7dpi", "14dpi"))

phage_Bac_melt_wide_metaSle_Culex <- phage_Bac_melt_wide_metaSle_Culex[,-5:-6]
phage_Bac_melt_wide_metaSle_Culex_m <- melt(phage_Bac_melt_wide_metaSle_Culex); head(phage_Bac_melt_wide_metaSle_Culex_m)
# pdf("./Rfigure/Cq_Genera_proportion_sepcific.pdf", width = 8,height = 6, useDingbats = FALSE)
ggplot(phage_Bac_melt_wide_metaSle_Culex_m, aes(x=dpi, y=value,col=variable)) + 
  # geom_point(shape = 1)+
  facet_grid(~Genus)+
  geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
  # facet_wrap(~ Genus)+
  # geom_jitter(size = 2, position=position_jitter(0.2))
  scale_shape_manual(values = c(9,19))+
  theme_bw()+
  # scale_color_manual(values = col_6_Ae[c(1:2,4:6)])+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12,vjust = 0.5),
        # plot.margin = unit(c(1, 1, 1, 1), "cm"),
        strip.text.x = element_text(size =12, colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  stat_compare_means()
dev.off()
# pdf("./Rfigure/Phage_vs_Bacteria_abundance_distribution.pdf", width = 30,height =20,useDingbats = FALSE)
# ggplot(phage_Bac_melt_wide_meta, aes(x=log10bacteria, y=log10phage)) + # shape=mosq_species, color=mosq_species
#   facet_wrap(~Genus+phageID)+
#   geom_point(size=1.8)+
#   # geom_smooth(method=lm)
#   geom_smooth()+
#   theme_bw()+
#   scale_color_manual(values = c(col[7],col[3]))+
#   ylab("log10 phage")+
#   xlab("log10 bacteria")
# dev.off()  