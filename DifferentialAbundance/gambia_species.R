library(metagenomeSeq)
library(msd16s)

feature_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "OTU")
fData(msd16s) <- fData(msd16s)[feature_order]

gambia <- msd16s[, which(pData(msd16s)$Country == "Gambia")]

gambia_filtered <- filterData(gambia, present = 5)

normed_gambia <-  cumNorm(gambia_filtered, p = 0.75)

aggregation_level <- "species"
aggregated_normed_gambia <- aggregateByTaxonomy(normed_gambia, lvl=aggregation_level, norm=FALSE)

normFactors(aggregated_normed_gambia) <- normFactors(normed_gambia)

gambia_sample_data <-  pData(aggregated_normed_gambia)
mod <-  model.matrix(~1+Dysentery, data = gambia_sample_data)
results_gambia <-  fitFeatureModel(aggregated_normed_gambia, mod)
logFC_gambia <- MRcoefs(results_gambia, number = nrow(aggregated_normed_gambia))

speciesOfInterest <- c("Rothia mucilaginosa", "Granulicatella adiacens", "Granulicatella elegans", "Granulicatella sp. oral clone ASCG05", 
                     "Streptococcus mitis", "Streptococcus oralis", "Streptococcus parasanguinis", "Streptococcus sanguinis", "Streptococcus sp. C101", 
                     "Streptococcus sp. oral clone ASCC01", "Streptococcus sp. oral clone ASCE09", "Citrobacter freundii", "Erwinia chrysanthemi", 
                     "Escherichia coli", "Klebsiella pneumoniae", "Haemophilus haemolyticus", "Haemophilus parainfluenzae", "Haemophilus sp. oral clone BP2-46", 
                     "Prevotella copri", "Prevotella histicola", "Prevotella sp. BI-42", "Prevotella sp. DJF_B112", "Prevotella sp. DJF_B116", "Prevotella sp. DJF_LS16", 
                     "Prevotella sp. DJF_RP53", "Prevotella sp. oral clone BP1-28")

logFC_gambia[speciesOfInterest,]