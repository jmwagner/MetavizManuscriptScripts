# Load required libraries
library(devtools)
install_github("epiviz/metavizr")
library(metavizr)
library(metagenomeSeq)
library(etec16s)

# Start Metaviz App
app <- startMetaviz(verbose=TRUE)

# Load etec16s dataset and select only samples in first 11 days 
data(etec16s)
etec16s <- etec16s[,-which(pData(etec16s)$Day>9)]
collinsella_indices <- which(fData(etec16s)[,"Species"] == "Collinsella aerofaciens ")
fData(etec16s)$Species[collinsella_indices] <- "Collinsella aerofaciens"

# Add a Kingdom level to hiearchy
featureData(etec16s)$Kingdom <- "Bacteria"
feature_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU.ID")
fData(etec16s) <- fData(etec16s)[,feature_order]
fData(etec16s)[which(fData(etec16s)[,"Phylum"] == "NULL"),]$Kingdom <- "Not_Annotated_Kingdom"

#Fit a time series analysis to the etec16s data, use B=100 to set permutations
timeSeriesFits <- fitMultipleTimeSeries(obj=etec16s,
                                            formula = abundance~id + time*class + AntiGiven,
                                            class="AnyDayDiarrhea",
                                            id="SubjectID",
                                            time="Day",
                                            lvl="Species",
                                            featureOrder = feature_order,
                                            C=0.3,
                                            B=100,
                                            seed=12345)

# Use if loading time series from RData
#load("timeSeriesFits.RData")

# Select features that have significant intervals of at least 2 days 
ts_intervals_results <- sapply(timeSeriesFits, function(x){class(x[["timeIntervals"]])})
ts_intervals_indices <- which(ts_intervals_results == "matrix")
splines_to_plot <- sapply(names(ts_intervals_indices), function(x){
  (unname(timeSeriesFits[[x]]$timeIntervals[,"Interval end"]) - unname(timeSeriesFits[[x]]$timeIntervals[,"Interval start"])) >= 2})
splines_to_plot <- names(which(splines_to_plot == TRUE))

# Create an MRexperiment obj with features of interest and hierarchy aggregated to species level
feature_order_agg <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
splinesMRexp <- ts2MRexperiment(timeSeriesFits, 
                                featureData = featureData(aggregateByTaxonomy(etec16s, lvl="Species", featureOrder = feature_order_agg)),
                                sampleNames = timeSeriesFits[[2]]$fit$timePoints)

splinesMRsubselect <- splinesMRexp[splines_to_plot,]
fData(splinesMRsubselect) <- fData(splinesMRsubselect)[order(fData(splinesMRsubselect)[,"Species"]),]
assayData(splinesMRsubselect)$counts <- assayData(splinesMRsubselect)$counts[rownames(fData(splinesMRsubselect)),]

# Add FacetZoom to Metaviz workspace
ic_plot <- app$plot(splinesMRsubselect, datasource_name = "etec16_base", 
                    control=metavizControl(norm = FALSE, aggregateAtDepth = 7), 
                    feature_order = feature_order_agg)

splineObj <- app$data_mgr$add_measurements(splinesMRsubselect, 
                                           datasource_name = "splines", 
                                           control = metavizControl(norm=FALSE, aggregateAtDepth = 7, log=FALSE), feature_order = feature_order_agg)
splineMeasurements <- splineObj$get_measurements()
splineChart <- app$chart_mgr$visualize("LinePlot", splineMeasurements)

selectedSpecies <- rownames(fData(splinesMRsubselect))
etec16s_agg <- aggregateByTaxonomy(etec16s, lvl="Species", featureOrder = feature_order_agg)

# Remove samples with zero counts for all features that were selected
sub_etec16s_obj <- etec16s_agg[which(fData(etec16s_agg)$Species %in% selectedSpecies),]

fData(sub_etec16s_obj) <- fData(sub_etec16s_obj)[order(fData(sub_etec16s_obj)[,"Species"]),]
assayData(sub_etec16s_obj)$counts <- assayData(sub_etec16s_obj)$counts[rownames(fData(sub_etec16s_obj)),]

# update settings on splineChart 
settings <- list(yMin = -4, yMax = 2, colLabel="label")
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
app$chart_mgr$set_chart_settings(splineChart, settings=settings, colors = colors)

# Add the counts for selected samples and features
splineObj_stacked_plot <- app$data_mgr$add_measurements(sub_etec16s_obj, datasource_name = "etec16s_sub_stacked_plot", 
                                                        control=metavizControl(norm = FALSE, aggregateAtDepth = 7), 
                                                        feature_order = feature_order_agg)

# Add a stacked plot to the Metaviz workspace for one case sample
splineMeasurements_stacked_plot_case <- splineObj_stacked_plot$get_measurements()[which(pData(sub_etec16s_obj)$SubjectID == "E01JH0016")]
splineChart_stacked_plot_case <- app$chart_mgr$visualize("StackedLinePlot", splineMeasurements_stacked_plot_case)
app$chart_mgr$set_chart_settings(splineChart_stacked_plot_case, settings = list(title = "Case"))

# Add a stacked plot to the Metaviz workspace for one control sample
splineMeasurements_stacked_plot_control <- splineObj_stacked_plot$get_measurements()[which(pData(sub_etec16s_obj)$SubjectID == "E01JH0029")]
splineChart_stacked_plot_control <- app$chart_mgr$visualize("StackedLinePlot", splineMeasurements_stacked_plot_control)
app$chart_mgr$set_chart_settings(splineChart_stacked_plot_control, settings = list(title = "Control"))
