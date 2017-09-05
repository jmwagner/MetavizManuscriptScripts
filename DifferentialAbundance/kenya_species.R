library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

kenya <- msd16s[, which(pData(msd16s)$Country == "Kenya")]

kenya_filtered <- filterData(kenya, present = 5)

normed_kenya <-  cumNorm(kenya_filtered, p = 0.75)

aggregation_level <- "species"
aggregated_normed_kenya <- aggregateByTaxonomy(normed_kenya, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_kenya) <- normFactors(normed_kenya)

kenya_sample_data <-  pData(aggregated_normed_kenya)
mod <-  model.matrix(~1+Dysentery, data = kenya_sample_data)
results_kenya <-  fitFeatureModel(aggregated_normed_kenya, mod)
logFC_kenya <- MRcoefs(results_kenya, number = nrow(aggregated_normed_kenya))

speciesOfInterest <- c("Veillonella parvula", "Veillonella sp. HF9", "Veillonella sp. oral clone VeillC8", "Veillonella sp. oral clone VeillD5", 
                       "Enterobacter cancerogenus", "Enterobacter cloacae", "Escherichia coli", "Escherichia sp. oral clone 3RH-30", "Klebsiella pneumoniae", 
                       "Haemophilus haemolyticus", "Haemophilus parainfluenzae", "Prevotella copri", "Prevotella histicola", "Prevotella sp. BI-42", "Prevotella sp. DJF_B112", 
                       "Prevotella sp. DJF_B116", "Prevotella sp. DJF_RP53")

logFC_kenya[speciesOfInterest,]