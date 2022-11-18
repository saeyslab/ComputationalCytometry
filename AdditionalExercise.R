library("flowCore") # For basic cytometry operations
library("PeacoQC") # For quality control
library("FlowSOM") # For FlowSOM clustering and functions
library("ggplot2") # For nice plots
library("ggpubr") # For nice plots

library("pheatmap") # For pretty heatmap plots
library("tidyr") # For nice data handling

source("Extra_functions.R") # To load in additional functions


# Set directories ----


# Preprocessing ----
manual <- readRDS("Data/AdditionalExercise/manual.RDS")


channels_of_interest <- colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33, 17,
                                       11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34,
                                       41, 26, 30, 28, 29, 25, 35)]


for (file in files){
  
  m <- manual[[basename(file)]]$matrix[exprs(ff_QC$FinalFF)[,"Original_ID"],]
  ff_subset <- ff_QC$FinalFF[m[,"CD45+CD66-"],] #ff_QC is the result of PeacoQC

  write.FCS(ff_subset, file.path(dir_prepr, basename(file)))
}

# FlowSOM ----


# Downstream analysis ----
