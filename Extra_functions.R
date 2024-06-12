plotScatter <- function (ff_after, 
                         ff_before = NULL, 
                         channels, 
                         adjust = 1, 
                         title = NULL)
{
  set.seed(1)
  if (!is.null(ff_before)){
    if (!"Original_ID" %in% colnames(ff_before)){
      original_ff_new <- PeacoQC:::AppendCellID(ff_before, 1:nrow(ff_before))
      to_plot_ff <- original_ff_new
    } else{
      to_plot_ff <- ff_before
    }
    idx <- sample(to_plot_ff@exprs[,"Original_ID"],20000 )
    ff_sub <- to_plot_ff[to_plot_ff@exprs[,"Original_ID"] %in% idx,]
    removed_meas <- setdiff(to_plot_ff@exprs[,"Original_ID"],
                            ff_after@exprs[,"Original_ID"])
    removed_meas <- which(ff_sub@exprs[,"Original_ID"] %in% removed_meas)
    
  } else{
    to_plot_ff <- ff_after
    idx <- sample(1:nrow(to_plot_ff), 20000)
    ff_sub <- to_plot_ff[idx,]
    
  }
  
  channels <- FlowSOM::GetMarkers(to_plot_ff, channels, exact = TRUE)
  
  df <- ff_sub@exprs[,names(channels)] %>% as.data.frame()
  colnames(df) <- gsub("-", "", colnames(df))
  
  colors <-  c("#0000BF", "#0000FF", "#0080FF", "#00FFFF", "#55FFAA", "#AAFF55",
               "#FFFF00", "#FF8000", "#FF0000", "#BF0000")
  
  
  p <- ggplot(data = df, aes(x = !!sym(colnames(df)[1]), y = !!sym(colnames(df)[2]))) +
    ggpointdensity::geom_pointdensity(adjust=adjust) +
    ggplot2::scale_color_gradientn(colours =
                                     colors) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab(channels[1]) +
    ylab(channels[2])
  
  if(!is.null(ff_before)){
    p <- p + geom_point(data = df[removed_meas,], color = "grey60")
  }
  
  if(!is.null(title))
    p <- p + ggtitle(title)
  return(p)
  
}

estimateGate <- function (ff, marker_values = list(`FSC-A` = c(min = 60000, max = 125000, 
                                                                range_min = 50000, range_max = 50000), 
                                                    `SSC-A` = c(max = 80000,
                                                                range_max = 50000)), 
                           plot = FALSE, plot_title = "Gate", tinypeak.removal = 1/25) 
{
  selection <- rep(TRUE, nrow(ff))
  for (m in names(marker_values)) {
    cutoffs <- c(flowDensity::deGate(ff, m, all.cuts = TRUE, 
                                     tinypeak.removal = tinypeak.removal), 
                 flowDensity::deGate(ff, m, use.upper = TRUE, upper = TRUE))
    for (type in c("min", "max")) {
      if (!is.na(marker_values[[m]][type])) {
        default <- marker_values[[m]][type]
        cutoff <- cutoffs[which.min(abs(cutoffs - default))]
        if (!is.na(marker_values[[m]][paste0("range_", 
                                             type)]) & abs(cutoff - default) > marker_values[[m]][paste0("range_", 
                                                                                                         type)]) {
          cutoff <- default
        }
        if (type == "min") 
          selection <- selection & ff@exprs[, m] > cutoff
        if (type == "max") 
          selection <- selection & ff@exprs[, m] < cutoff
        marker_values[[m]][type] <- cutoff
      }
    }
    if (!"plot_min" %in% names(marker_values[[m]])) 
      marker_values[[m]]["plot_min"] <- NA
    if (!"plot_max" %in% names(marker_values[[m]])) 
      marker_values[[m]]["plot_max"] <- NA
  }
  if (plot) {
    p <- list()
    for (i in seq(1, length(marker_values), 2)) {
      m1 <- names(marker_values)[i]
      if (i + 1 <= length(marker_values)) {
        m2 <- names(marker_values)[i + 1]
      }
      else {
        m1 <- colnames(ff)[1]
        m2 <- names(marker_values)[i]
      }
      if (!"Original_ID" %in% colnames(ff)) {
        ff <- flowCore::fr_append_cols(ff, matrix(seq_len(nrow(ff)), 
                                                  ncol = 1, dimnames = list(NULL, "Original_ID")))
      }
      p[[length(p) + 1]] <- filter_plot(ff, ff[selection, 
      ], "", m1, m2) + ggplot2::geom_vline(xintercept = marker_values[[m1]]["min"], 
                                           color = c("cyan")) + ggplot2::geom_vline(xintercept = marker_values[[m1]]["max"], 
                                                                                    color = c("red")) + ggplot2::geom_hline(yintercept = marker_values[[m2]]["min"], 
                                                                                                                            color = c("cyan")) + ggplot2::geom_hline(yintercept = marker_values[[m2]]["max"], 
                                                                                                                                                                     color = c("red")) + ggplot2::theme_minimal() + 
        ggplot2::ggtitle(plot_title)
      if (!is.na(marker_values[[m1]]["plot_min"]) | !is.na(marker_values[[m1]]["plot_max"])) 
        p[[length(p)]] <- p[[length(p)]] + ggplot2::xlim(marker_values[[m1]]["plot_min"], 
                                                         marker_values[[m1]]["plot_max"])
      if (!is.na(marker_values[[m2]]["plot_min"]) | !is.na(marker_values[[m2]]["plot_max"])) 
        p[[length(p)]] <- p[[length(p)]] + ggplot2::ylim(marker_values[[m2]]["plot_min"], 
                                                         marker_values[[m2]]["plot_max"])
    }
    return(list(selection = selection, plot = p))
  }
  else {
    return(selection)
  }
}

filter_plot <- function (ff_pre, ff_post, title, channel_x, channel_y, n = 10000) 
{
  df <- data.frame(x = flowCore::exprs(ff_pre)[, channel_x], 
                   y = flowCore::exprs(ff_pre)[, channel_y])
  i <- sample(nrow(df), min(n, nrow(df)))
  if (!"Original_ID" %in% colnames(flowCore::exprs(ff_pre))) {
    ff_pre <- flowCore::fr_append_cols(ff_pre, matrix(seq_len(nrow(df)), 
                                                      ncol = 1, dimnames = list(NULL, c("Original_ID"))))
  }
  p <- ggplot2::ggplot(df[i, ], ggplot2::aes(x = .data$x, y = .data$y)) + 
    ggplot2::geom_point(size = 0.5, color = ifelse(flowCore::exprs(ff_pre)[i, 
                                                                           "Original_ID"] %in% flowCore::exprs(ff_post)[, "Original_ID"], 
                                                   "blue", "red")) + ggplot2::xlab(paste(flowCore::getChannelMarker(ff_pre, 
                                                                                                                    channel_x), collapse = " ")) + ggplot2::ylab(paste(flowCore::getChannelMarker(ff_pre, 
                                                                                                                                                                                                  channel_y), collapse = " ")) + ggplot2::theme_minimal() + 
    ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(paste0(round(flowCore::nrow(ff_post)/flowCore::nrow(ff_pre) * 
                                                                               100, 2), "% ", title))
  return(p)
}

plotDens <- function(input,
                     channels,
                     level,
                     plotFile){
  plotlist <- list()
  for (channel in channels){
    plotlist[[length(plotlist)+1]] <- ggplot(data = input, aes(x = !!sym(channel),
                                                               color = !!sym(level))) +
      stat_density(aes(group = .data$File),
                   geom = "line", position = "identity",
                   alpha = 0.2) +
      stat_density(aes(group = !!sym(level)),
                   geom = "line", position = "identity",
                   alpha = 1) +
      xlab(GetMarkers(object = ff,
                      channels = channel)) +
      theme_minimal()
  }
  pdf(plotFile, width = 20, height = 12)
  print(ggarrange(plotlist = plotlist, common.legend = TRUE, nrow = 5, ncol = 5))
  dev.off()
}

plotUMAP <- function(input,
                     level,
                     plotFile){
  p <- ggplot(input, aes(x = UMAP1, y = UMAP2, color = !!sym(level))) +
    geom_point() +
    theme_minimal()
  pdf(plotFile, width = 20, height = 12)
  print(p)
  dev.off()
}

plotGate <- function(ff, 
                     gate_coord){
  gate_df <- data.frame(x = gate_coord[, 1],
                        y = gate_coord[, 2],
                        xend = c(gate_coord[-1,1], gate_coord[1,1]),
                        yend = c(gate_coord[-1,2], gate_coord[1,2]))
  
  channel_x <- colnames(gate_coord)[1]
  channel_y <- colnames(gate_coord)[2]
  
  p <- ggplot(as.data.frame(exprs(ff))[sample(nrow(exprs(ff)), 10000),],
              aes_string(x = paste0("`", channel_x, "`"), y = paste0("`", channel_y, "`"))) +
    geom_point() +
    geom_segment(data = gate_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "red", linewidth = 1.5) +
    geom_point(data = gate_df,
               aes(x = x, y = y),
               color = "red", size = 2.5) +
    theme_minimal()
  
  return(p)
}