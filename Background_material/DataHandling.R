##### Data ComputationalCytometry ##############################################
dir.create("Background_material/Subsetted/")
dir.create("Background_material/Preprocessed//")
dir.create("Background_material/Batched/")
dir.create("Background_material/Final/")


#### List raw files ####
filesKO <- list.files(path = "Background_material/Raw", 
                      pattern = ".fcs", 
                      full.names = TRUE)[1:3]
filesWT <- list.files(path = "Background_material/Raw", 
                      pattern = ".fcs", 
                      full.names = TRUE)[4:7]

#### Subset data and give more informative names ####
subset_i <- matrix(NA,
                   nrow = 300000, ncol = 10,
                   dimnames = list(NULL, 
                                   c(basename(filesWT),
                                     basename(filesKO),
                                     basename(filesKO[1]),
                                     rep(basename(filesWT[1]), 2))))

day <- rep(rep(c(1,2), each = 2),2)

tube <- 1
set.seed(1)
for (file in filesWT){
  ff <- flowCore::read.FCS(file)
  i <- sort(sample(1:nrow(ff), 300000))
  flowCore::write.FCS(ff[i,], paste0("Background_material/Subsetted/Tube0", tube, "_WT_", "Day", day[tube], ".fcs"))
  subset_i[,tube] <- i
  tube <- tube + 1
}

for (file in filesKO){
  ff <- flowCore::read.FCS(file)
  i <- sort(sample(1:nrow(ff), 300000))
  flowCore::write.FCS(ff[i,], paste0("Background_material/Subsetted/Tube0", tube, "_KO_", "Day", day[tube], ".fcs"))
  subset_i[,tube] <- i
  tube <- tube + 1
}

# Sample other file out of first file (to have 2 KO samples in batch 2)
ff <- flowCore::read.FCS(filesKO[1])
i <- sort(sample(1:nrow(ff), 300000))  
flowCore::write.FCS(ff[i,], paste0("Background_material/Subsetted/Tube0", tube, "_KO_", "Day", day[tube], ".fcs"))
subset_i[,tube] <- i
tube <- tube + 1

# Sample 2 other file out of first WT file (to have control in each batch)
ff <- flowCore::read.FCS(filesWT[1])
for (ctrl in 1:2){
  i <- sort(sample(1:nrow(ff), 300000))  
  flowCore::write.FCS(ff[i,], paste0("Background_material/Subsetted/Control_Day", ctrl, ".fcs"))
  subset_i[,tube] <- i
  tube <- tube + 1
}


saveRDS(subset_i, "Background_material/Subsetted/subsetIndices.rds")

#### Get manual labels ####
gating <- FlowSOM::GetFlowJoLabels(files = basename(c(filesWT, filesKO)),
                          wspFile = "Background_material/Raw/General_panel.wsp")
files <- list.files(path = "Background_material/Subsetted", 
                    pattern = ".fcs")[c(3:10, 1:2)]

manual_labeling <- list()
for (col in 1:ncol(subset_i)){
  manual_labeling[[files[col]]] <- gating[[colnames(subset_i)[col]]][["manual"]][subset_i[,col]]
}

saveRDS(manual_labeling, "Background_material/Subsetted/ManualLabeling.rds")

#### Preprocess files ####
files <- list.files(path = "Background_material/Subsetted",
                    pattern = ".fcs",
                    full.names = TRUE)
p5 <- NULL
for (file in files){
  ff <- flowCore::read.FCS(file)
  comp <- ff@description$SPILL
  ff_c <- flowCore::compensate(ff, comp)
  cols_to_transform <- colnames(comp)
  t <- flowCore::arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=1/150)
  t_list <- flowCore::transformList(cols_to_transform, t)
  ff_t <- flowCore::transform(ff_c, t_list)
  SSCA <- ff_t@exprs[,"SSC-A"]
  if (is.null(p5)){ # Calculate percentiles on first file and apply to all
    p5 <- quantile(SSCA, 0.05)
    p95 <- quantile(SSCA, 0.95)
  }
  SSCA <- 10 * (SSCA - p5) / (p95 - p5)
  ff_t@exprs[,"SSC-A"] <- SSCA
  flowCore::write.FCS(ff_t, sub("Subsetted", "Preprocessed", file))
}

#### Introduce batch effects ####
files <- list.files(path = "Background_material/Preprocessed",
                    pattern = ".fcs",
                    full.names = TRUE)
files_batch <- list.files(path = "Background_material/Preprocessed",
                          pattern = "2.fcs",
                          full.names = TRUE)

channels_of_interest <- c("PE-Cy7-A", "PE-Cy5-A", "PE-A", "Alexa Fluor 700-A", 
                          "APC-A", "SSC-A", "PerCP-Cy5-5-A", "Pacific Blue-A", 
                          "BV605-A","BV711-A", "BV786-A")
n <- length(channels_of_interest)
transformation <- matrix(c(0.5, -0.5, 0.3, 1, 0.9, 1.1),
                         nrow = 3,
                         ncol = 2,
                         dimnames = list(c("PE-A", "PerCP-Cy5-5-A", "APC-A"), # Other transformation values for each channel
                                         c("shift", "scale"))) # 2 transformations
for (file in files_batch){
  ff <- flowCore::read.FCS(file)
  for (channel in channels_of_interest){
    if (channel == "PE-A"){ # Only affect positive population for CD3
      gate <- flowDensity::deGate(ff, channel)
      pos <- ff@exprs[,channel] > gate
      ff@exprs[pos,channel] <- ff@exprs[pos,channel] + transformation[channel, "shift"]

    } else if (channel == "PerCP-Cy5-5-A"){ #Only affect MHCII expression of DCs
      DC <- manual_labeling[[basename(file)]] == "DCs"
      ff@exprs[DC,channel] <- (ff@exprs[DC,channel] + transformation[channel, "shift"])
      ff@exprs[DC,channel] <- scales::rescale(ff@exprs[DC,channel],
                                                 to=c(min(ff@exprs[DC,channel]),
                                                      max(ff@exprs[DC,channel]*transformation[channel, "scale"])))
    } else if (channel == "APC-A"){# Affect whole population for CD161
      ff@exprs[,channel] <- ff@exprs[,channel] + transformation[channel, "shift"]
      ff@exprs[,channel] <- scales::rescale(ff@exprs[,channel],
                                            to=c(min(ff@exprs[,channel]),
                                                 max(ff@exprs[,channel])*transformation[channel, "scale"]))
    }
  }
  flowCore::write.FCS(ff, sub("Preprocessed", "Batched", file))
}

#### Backtransform and save files ####
files_batched <- list.files(path = "Background_material/Batched",
                            pattern = "2.fcs",
                            full.names = TRUE)
files_or <- list.files(path = "Background_material/Subsetted",
                       pattern = "1.fcs",
                       full.names = TRUE)

for (file in files_batched){
  ff <- flowCore::read.FCS(file)
  
  SSCA <- ff@exprs[,"SSC-A"]
  SSCA <- SSCA * (p95 - p5) / 10 + p5
  ff@exprs[,"SSC-A"] <- SSCA
  
  ff@exprs[,cols_to_transform] <- apply(X = ff@exprs[,cols_to_transform], 
                                        MARGIN = 2,
                                        FUN = function(x) sinh(x)*150)
  ff_or <- flowCore::read.FCS(sub("Batched", "Subsetted", file))
  ff_or@exprs <- ff@exprs
  
  ff_or <- flowCore::decompensate(ff_or, ff_or@description$SPILL)
  flowCore::write.FCS(ff_or, sub("Batched", "Final", file))
}

for (file in files_or){
  ff <- flowCore::read.FCS(file)
  flowCore::write.FCS(ff, sub("Subsetted", "Final", file))
}


##### Data ComputationalCytometry Beginner #####################################


##### Data AdditionalExercise ##################################################
files <- list.files(path = "Data/AdditionalExercise",
                    pattern = ".fcs",
                    full.names = TRUE)
for (file in files){
  ff <- read.FCS(file)
  write.FCS(ff, file)
}

  