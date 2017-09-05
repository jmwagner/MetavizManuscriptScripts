library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

bangladesh <- msd16s[, which(pData(msd16s)$Country == "Bangladesh")]

bangladesh_filtered <- filterData(bangladesh, present = 5)

normed_bangladesh <-  cumNorm(bangladesh_filtered, p = 0.75)

aggregation_level <- "species"
aggregated_normed_bangladesh <- aggregateByTaxonomy(normed_bangladesh, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_bangladesh) <- normFactors(normed_bangladesh)

bangladesh_sample_data <-  pData(aggregated_normed_bangladesh)
mod <-  model.matrix(~1+Dysentery, data = bangladesh_sample_data)
results_bangladesh <-  fitFeatureModel(aggregated_normed_bangladesh, mod)
logFC_bangladesh <- MRcoefs(results_bangladesh, number = nrow(aggregated_normed_bangladesh))

speciesOfInterest <- c("Escherichia coli", "Escherichia sp. oral clone 3RH-30", "Granulicatella adiacens", "Streptococcus equinus", "Streptococcus mitis", "Streptococcus parasanguinis", "Streptococcus salivarius", 
                       "Haemophilus parainfluenzae", "Acinetobacter sp. SF6", "Collinsella sp. CB20", "Bacteroides fragilis", "Faecalibacterium prausnitzii", "Faecalibacterium sp. DJF_VR20", "Ruminococcus gnavus")

logFC_bangladesh[speciesOfInterest,]