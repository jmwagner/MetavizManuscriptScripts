library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

mali <- msd16s[, which(pData(msd16s)$Country == "Mali")]

mali_filtered <- filterData(mali, present = 5)

normed_mali <-  cumNorm(mali_filtered, p = 0.75)

aggregation_level <- "order"
aggregated_normed_mali <- aggregateByTaxonomy(normed_mali, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_mali) <- normFactors(normed_mali)

mali_sample_data <-  pData(aggregated_normed_mali)
mod <-  model.matrix(~1+Dysentery, data = mali_sample_data)
results_mali <-  fitFeatureModel(aggregated_normed_mali, mod)
logFC_mali <- MRcoefs(results_mali, number = nrow(aggregated_normed_mali))

orderOfInterest <- c("Actinomycetales", "Neisseriales", "Fusobacteriales", "Enterobacteriales", "Pasteurellales", "Pseudomonadales")

logFC_mali[orderOfInterest,]