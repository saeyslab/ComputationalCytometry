files <- list.files(path = "Data/Preprocessed", pattern = ".fcs", full.names = TRUE)
files

channels <- c(`SSC-A` = "SSC-A", `MHCII#PerCp-Cy5-5` = "PerCP-Cy5-5-A", `CD49b#eFluor450` = "Pacific Blue-A", 
              `CD11b#BV605` = "BV605-A", `CD64#BV711` = "BV711-A", `FcERI#BV786` = "BV786-A", 
              `CD161#APC` = "APC-A", `Ly-6G#AF700` = "Alexa Fluor 700-A", `CD3#PE` = "PE-A", 
              `CD19#PE-Cy5` = "PE-Cy5-A", `CD11c#PE-Cy7` = "PE-Cy7-A")
markers <- c("SSC_A", "MHCII", "CD49b", "CD11b", "CD64", "FcERI", "CD161", 
             "Ly_6G", "CD3", "CD19", "CD11c")
day <- sub(".*(Day.).*", "\\1", files)

agg <- AggregateFlowFrames(fileNames = files, 
                           channels = c(channels, "Original_ID"), 
                           cTotal = 1000000)

#### Plot all data ####
df <- data.frame(agg@exprs)
colnames(df) <- c(markers, "Original_ID", "File", "File_scattered", "Original_ID2")
df$Day <- factor(day[df$File])

plotlist <- list()
for (marker in markers){
  df_sub <- df[,c(marker, "File", "Day")]
  p <- ggplot(df_sub, aes_string(marker)) +
    geom_density(aes(color = Day, group = File)) +
    theme_minimal()
    
  plotlist[[length(plotlist)+1]] <- p
}

pdf("Background_material/BatchEffect_AllData.pdf", height = 12, width = 6)
ggarrange(plotlist = plotlist, ncol = 2, nrow = 6)
dev.off()

#### Plot DCs ####
manual <- readRDS("Data/Raw/ManualLabeling.rds")

labels <- rep(NA, nrow(df))
for (i in unique(df$File)){
  man <- as.character(manual[[basename(files[i])]])
  file_i <- df$File == i
  labels[file_i] <- man[df$Original_ID[file_i]]
}
df$CellType <- labels

plotlist <- list()
for (marker in markers){
  df_sub <- df[df$CellType == "DCs", c(marker, "File", "Day")]
  p <- ggplot(df_sub, aes_string(marker)) +
    geom_density(aes(color = Day, group = File)) +
    theme_minimal()
  
  plotlist[[length(plotlist)+1]] <- p
}

pdf("Background_material/BatchEffect_DCs.pdf", height = 12, width = 6)
ggarrange(plotlist = plotlist, ncol = 2, nrow = 6)
dev.off()



#### Plot B cells ####
plotlist <- list()
for (marker in markers){
  df_sub <- df[df$CellType == "B cells", c(marker, "File", "Day")]
  p <- ggplot(df_sub, aes_string(marker)) +
    geom_density(aes(color = Day, group = File)) +
    theme_minimal()
  
  plotlist[[length(plotlist)+1]] <- p
}

pdf("Background_material/BatchEffect_Bcells.pdf", height = 12, width = 6)
ggarrange(plotlist = plotlist, ncol = 2, nrow = 6)
dev.off()
