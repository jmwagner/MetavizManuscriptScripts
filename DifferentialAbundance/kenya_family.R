library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

kenya <- msd16s[, which(pData(msd16s)$Country == "Kenya")]

kenya_filtered <- filterData(kenya, present = 5)

normed_kenya <-  cumNorm(kenya_filtered, p = 0.75)

aggregation_level <- "family"
aggregated_normed_kenya <- aggregateByTaxonomy(normed_kenya, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_kenya) <- normFactors(normed_kenya)

kenya_sample_data <-  pData(aggregated_normed_kenya)
mod <-  model.matrix(~1+Dysentery, data = kenya_sample_data)
results_kenya <-  fitFeatureModel(aggregated_normed_kenya, mod)
logFC_kenya <- MRcoefs(results_kenya, number = nrow(aggregated_normed_kenya))

familyOfInterest <- c("Veillonellaceae", "Campylobacteraceae", "Enterobacteriaceae", "Pasteurellaceae", "Prevotellaceae")

logFC_kenya[familyOfInterest,]