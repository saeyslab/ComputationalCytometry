library("flowCore") # For basic cytometry operations
library("PeacoQC") # For quality control
library("FlowSOM") # For FlowSOM clustering and functions
library("ggplot2") # For nice plots
library("ggpubr") # For nice plots

library("pheatmap") # For pretty heatmap plots
library("tidyr") # For nice data handling

source("Extra_functions.R") # To load in additional functions


# Set directories ----

dir_raw <- "Data/AdditionalExercise/" #where the raw data is located
dir_prepr <- "Data/AdditionalExercise/Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "Data/AdditionalExercise//Preprocessed/QC/" #where the data QC results will be stored
dir_RDS <- "RDS/AdditionalExercise/" #where the R objects will be stored
dir_results <- "Results/AdditionalExercise/" #where the results will be stored

for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  if(!dir.exists(path)) dir.create(path)
}

# Preprocessing ----
# NOTE: Do not apply removeMargins or compensation

files <- list.files(dir_raw, pattern = ".fcs", full.names = TRUE)
manual <- readRDS("Data/AdditionalExercise/manual.RDS")

ff <- read.FCS(files[1])

channels_of_interest <- colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                                       11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34,
                                       41, 26, 30, 28, 29, 25, 35)]
markers_of_interest <- GetMarkers(ff, channels_of_interest)

transformList <- transformList(channels_of_interest,
                               arcsinhTransform(a = 0, b = 1/5))

for (file in files){

  print(file)

  ff <- read.FCS(file)
  ff_t <- transform(ff, transformList)

  ff_QC <- PeacoQC::PeacoQC(ff = ff_t,
                            channels = channels_of_interest,
                            output_directory = dir_QC,
                            plot = TRUE,
                            save_fcs = TRUE,
                            remove_zeros = TRUE)
  p1 <- plotScatter(ff_after = ff_QC$FinalFF, 
                    ff_before = ff_t,
                    channels = GetChannels(ff, c("Time", "CCR7"), exact = FALSE),
                    adjust = 1,
                    title = "Remove low quality events")
  

  m <- manual[[basename(file)]]$matrix[exprs(ff_QC$FinalFF)[,"Original_ID"],]
  ff_subset <- ff_QC$FinalFF[m[,"CD45+CD66-"],]

  p2 <- plotScatter(ff_after = ff_subset, 
                    ff_before = ff_QC$FinalFF,
                    channels = GetChannels(ff, c("CD45", "CD66"), exact = FALSE),
                    adjust = 1,
                    title = "Select CD45 cells")
  

  png(paste0(dir_QC, sub(".fcs", ".png", basename(file))), height = 300, width = 600)
  print(ggarrange(p1,p2, nrow = 1))
  dev.off()

  write.FCS(ff_subset, file.path(dir_prepr, basename(file)))
}

# FlowSOM ----

files <- list.files(path = dir_prepr,
                    pattern = ".fcs",
                    full.names = TRUE)
files

set.seed(2022)
agg <- AggregateFlowFrames(fileNames = files,
                           cTotal = 500000)

agg_manual <- c()
i <- 1
for (file in files){
  agg_manual <- c(agg_manual,
                  as.character(manual[[basename(file)]][["manual"]][(exprs(agg)[,"Original_ID"])[exprs(agg)[,"File"]==i]]))
  i <- i + 1
}
agg_manual <- factor(agg_manual,
                     levels = c("Unknown",
                                "Naive CD4+ T cells", "Memory CD4+ T cells", 
                                "Naive CD8+ T cells", "Memory CD8+ T cells", 
                                "TCRgd T cells", "pDCs", "Macrophages", 
                                "B cells", "NK cells"))

fsom <- FlowSOM(input = agg,
                colsToUse = channels_of_interest,
                seed = 2022,
                nClus = 30,
                xdim = 15, ydim = 15)
saveRDS(fsom, paste0(dir_RDS, "fsom_additionalexercise.rds"))

FlowSOMmary(fsom,
            plotFile = file.path(dir_results, "FlowSOMmary.pdf"))

# Downstream analysis ----
fsom <- readRDS(paste0(dir_RDS, "fsom_additionalexercise.rds"))


PlotPies(fsom, fsom$data[,"File"])

fsom_1Unstim <- NewData(fsom,
                        "Data/AdditionalExercise/Preprocessed/Gates_PTLG021_Unstim_Control_1.fcs")
fsom_1Stim <- NewData(fsom,
                      "Data/AdditionalExercise/Preprocessed/Gates_PTLG021_IFNa_LPS_Control_1.fcs")

p1 <- PlotStars(fsom_1Unstim, list_insteadof_ggarrange = TRUE)
p2 <- PlotStars(fsom_1Stim, list_insteadof_ggarrange = TRUE)


diff  <-  100 * (GetPercentages(fsom_1Stim, level = "clusters") -
  GetPercentages(fsom_1Unstim, level = "clusters"))
diff[is.na(diff)] <- 0

p3 <- PlotVariable(fsom, diff, "Difference",
             colorPalette = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(9,"RdBu"))),
             equalNodeSize = TRUE)
p4 <- PlotPies(fsom = fsom, 
               cellTypes = agg_manual)

gridExtra::grid.arrange(nrow = 2,
                        p1$tree + ggtitle("Unstimulated"),
                        p2$tree + ggtitle("Stimulated"),
                        p3 + ggtitle("Difference"),
                        p4 + ggtitle("Manual labeling"))

