library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

bangladesh <- msd16s[, which(pData(msd16s)$Country == "Bangladesh")]

bangladesh_filtered <- filterData(bangladesh, present = 5)

normed_bangladesh <-  cumNorm(bangladesh_filtered, p = 0.75)

aggregation_level <- "family"
aggregated_normed_bangladesh <- aggregateByTaxonomy(normed_bangladesh, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_bangladesh) <- normFactors(normed_bangladesh)

bangladesh_sample_data <-  pData(aggregated_normed_bangladesh)
mod <-  model.matrix(~1+Dysentery, data = bangladesh_sample_data)
results_bangladesh <-  fitFeatureModel(aggregated_normed_bangladesh, mod)
logFC_bangladesh <- MRcoefs(results_bangladesh, number = nrow(aggregated_normed_bangladesh))

familyOfInterest <- c("Micrococcaceae", "Enterobacteriaceae", "Carnobacteriaceae", "Streptococcaceae", "Pasteurellaceae", "Moraxellaceae", "Coriobacteriaceae", "Bacteroidaceae", "Porphyromonadaceae", "Clostridiaceae", "Eubacteriaceae", "Lachnospiraceae", "Ruminococcaceae"
)

logFC_bangladesh[familyOfInterest,]