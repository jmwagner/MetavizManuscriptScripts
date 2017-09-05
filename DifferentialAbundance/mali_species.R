library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

mali <- msd16s[, which(pData(msd16s)$Country == "Mali")]

mali_filtered <- filterData(mali, present = 5)

normed_mali <-  cumNorm(mali_filtered, p = 0.75)

aggregation_level <- "species"
aggregated_normed_mali <- aggregateByTaxonomy(normed_mali, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_mali) <- normFactors(normed_mali)

mali_sample_data <-  pData(aggregated_normed_mali)
mod <-  model.matrix(~1+Dysentery, data = mali_sample_data)
results_mali <-  fitFeatureModel(aggregated_normed_mali, mod)
logFC_mali <- MRcoefs(results_mali, number = nrow(aggregated_normed_mali))

speciesOfInterest <- c("Rothia mucilaginosa", "Citrobacter freundii", "Erwinia chrysanthemi", "Enterobacter cancerogenus", "Enterobacter cloacae", "Escherichia albertii", 
                       "Escherichia coli", "Escherichia sp. oral clone 3RH-30", "Klebsiella pneumoniae", "Shigella boydii", "Shigella sonnei", "Haemophilus parainfluenzae", 
                       "Haemophilus sp. oral clone BP2-46", "Acinetobacter sp. SF6", "Bifidobacterium longum", "Bacteroides fragilis", "Prevotella copri", "Prevotella histicola", 
                       "Prevotella sp. BI-42", "Prevotella sp. DJF_B112", "Prevotella sp. DJF_RP53")

logFC_mali[speciesOfInterest,]