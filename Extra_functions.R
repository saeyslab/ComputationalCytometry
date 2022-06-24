plotDens <- function (ff, markers, adjust = 1)
{

  idx <- sample(1:nrow(ff),10000 )
  ff_sub <- ff[idx,]

  channels <- GetChannels(ff, markers, exact = FALSE)

  df <- ff_sub@exprs[,channels] %>% as.data.frame()
  colnames(df) <- markers

  colors <-  c("#0000BF", "#0000FF", "#0080FF", "#00FFFF", "#55FFAA", "#AAFF55",
               "#FFFF00", "#FF8000", "#FF0000", "#BF0000")


  p <- ggplot(data = df, aes_string(x = markers[1], y = markers[2])) +
    ggpointdensity::geom_pointdensity(adjust=adjust) +
    ggplot2::scale_color_gradientn(colours =
                                     colors) +
    theme_minimal() +
    theme(legend.position = "none")


  p

}



estimateTransformation <- function(files, channels, cTotal = 5000*length(files)){
  if(length(names(channels))>0){
    channels <- unname(channels)
  }

  ff_agg <- AggregateFlowFrames(files, cTotal = cTotal)
  transformList <- estimateLogicle(ff_agg, channels)


  return(transformList)
}
