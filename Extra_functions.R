plotDens <- function (ff, original_ff = NULL, markers, adjust = 1)
{

  if (!is.null(original_ff)){
    if (!"Original_ID" %in% colnames(original_ff)){
      original_ff_new <- PeacoQC:::AppendCellID(original_ff, 1:nrow(original_ff))
      to_plot_ff <- original_ff_new
    } else{
      to_plot_ff <- original_ff
    }
    idx <- sample(to_plot_ff@exprs[,"Original_ID"],20000 )
    ff_sub <- to_plot_ff[to_plot_ff@exprs[,"Original_ID"] %in% idx,]
    removed_meas <- setdiff(to_plot_ff@exprs[,"Original_ID"],
                            ff@exprs[,"Original_ID"])
    removed_meas <- which(ff_sub@exprs[,"Original_ID"] %in% removed_meas)

  } else{
    to_plot_ff <- ff
    idx <- sample(1:nrow(ff), 20000)
    ff_sub <- to_plot_ff[idx,]

  }



  channels <- GetChannels(to_plot_ff, markers, exact = FALSE)

  df <- ff_sub@exprs[,channels] %>% as.data.frame()
  colnames(df) <- markers
  colnames(df) <- gsub("-", "", colnames(df))

  colors <-  c("#0000BF", "#0000FF", "#0080FF", "#00FFFF", "#55FFAA", "#AAFF55",
               "#FFFF00", "#FF8000", "#FF0000", "#BF0000")


  p <- ggplot(data = df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
    ggpointdensity::geom_pointdensity(adjust=adjust) +
    ggplot2::scale_color_gradientn(colours =
                                     colors) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab(markers[1]) +
    ylab(markers[2])

  if(!is.null(original_ff)){
    p <- p + geom_point(data = df[removed_meas,], color = "black")
  }

  return(p)

}



estimateTransformation <- function(files, channels, cTotal = 5000*length(files)){
  if(length(names(channels))>0){
    channels <- unname(channels)
  }

  ff_agg <- AggregateFlowFrames(files, cTotal = cTotal)
  transformList <- estimateLogicle(ff_agg, channels)


  return(transformList)
}
