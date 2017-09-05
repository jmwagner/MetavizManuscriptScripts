library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

gambia <- msd16s[, which(pData(msd16s)$Country == "Gambia")]

gambia_filtered <- filterData(gambia, present = 5)

normed_gambia <-  cumNorm(gambia_filtered, p = 0.75)

aggregation_level <- "family"
aggregated_normed_gambia <- aggregateByTaxonomy(normed_gambia, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_gambia) <- normFactors(normed_gambia)

gambia_sample_data <-  pData(aggregated_normed_gambia)
mod <-  model.matrix(~1+Dysentery, data = gambia_sample_data)
results_gambia <-  fitFeatureModel(aggregated_normed_gambia, mod)
logFC_gambia <- MRcoefs(results_gambia, number = nrow(aggregated_normed_gambia))

familyOfInterest <- c("Actinomycetaceae", "Micrococcaceae", "Carnobacteriaceae", "Streptococcaceae", "Campylobacteraceae", "Enterobacteriaceae", "Pasteurellaceae", "Moraxellaceae", "Porphyromonadaceae", "Prevotellaceae", "Eubacteriaceae")

logFC_gambia[familyOfInterest,]