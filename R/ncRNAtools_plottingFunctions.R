plotPairsProbabilityMatrix <- function(basePairProbsMatrix, probabilityThreshold=0.1, filename=NULL,
                                       colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep="")) {
  Var1 <- NULL
  Var2 <- NULL
  value <- NULL
  meltedMatrix <- melt(basePairProbsMatrix)
  plot <- ggplot(data = meltedMatrix, aes(x=Var1, y=Var2, fill=value)) +
    geom_abline(intercept=nrow(basePairProbsMatrix)+1,
                slope=-1,
                colour="darkgrey") +
    geom_hline(yintercept=(seq(nrow(basePairProbsMatrix), 0, by=-10)+1),
               color="lightgrey",
               size=0.2) +
    geom_vline(xintercept=seq(10, nrow(basePairProbsMatrix), by=10),
               color="lightgrey",
               size=0.2) +
    geom_tile() +
    scale_fill_gradientn(colours=c("#FFFFFF00", colorPalette),
                         values=c(0, seq(probabilityThreshold, 1, length.out=length(colorPalette))),
                         breaks=c(0, 0.5, 1),
                         labels=c(0, 0.5, 1),
                         limits=c(0, 1),
                         name="Probability") +
    scale_x_discrete(position="top",
                     breaks=levels(meltedMatrix$Var1)[seq(10, nrow(basePairProbsMatrix), by=10)]) +
    scale_y_discrete(breaks=levels(meltedMatrix$Var1)[seq(10, nrow(basePairProbsMatrix), by=10)],
                     limits=rev(levels(meltedMatrix$Var2))) +
    coord_equal() +
    theme(panel.border=element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill = 'white', colour="white"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(color="black",size=12),
          axis.text.x = element_text(angle = 90))
  if (!is.null(filename)) {
    ggsave(filename, plot)
  }
  return(plot)
}

plotCompositePairsMatrix <- function(basePairProbsMatrix, pairedBases, probabilityThreshold=0.1, filename=NULL,
                                     colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep="")) {
  Var1 <- NULL
  Var2 <- NULL
  value <- NULL
  compositeMatrix <- makeCompositeMatrix(basePairProbsMatrix, pairedBases)
  compositeMatrix[upper.tri(compositeMatrix) & (compositeMatrix == 1)] <- NA
  meltedMatrix <- melt(compositeMatrix)
  plot <- ggplot(data = meltedMatrix, aes(x=Var1, y=Var2, fill=value)) +
    geom_abline(intercept=nrow(basePairProbsMatrix)+1,
                slope=-1,
                colour="darkgrey") +
    geom_hline(yintercept=(seq(nrow(basePairProbsMatrix), 0, by=-10)+1),
               color="lightgrey",
               size=0.2) +
    geom_vline(xintercept=seq(10, nrow(basePairProbsMatrix), by=10),
               color="lightgrey",
               size=0.2) +
    geom_tile() +
    scale_fill_gradientn(colours=c("#FFFFFF00", colorPalette),
                         values=c(0, seq(probabilityThreshold, 1, length.out=length(colorPalette))),
                         breaks=c(0, 0.5, 1),
                         labels=c(0, 0.5, 1),
                         limits=c(0, 1),
                         name="Probability",
                         na.value="black") +
    scale_x_discrete(position="top",
                     breaks=levels(meltedMatrix$Var1)[seq(10, nrow(basePairProbsMatrix), by=10)]) +
    scale_y_discrete(breaks=levels(meltedMatrix$Var1)[seq(10, nrow(basePairProbsMatrix), by=10)],
                     limits=rev(levels(meltedMatrix$Var2))) +
    coord_equal() +
    theme(panel.border=element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill = 'white', colour="white"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(color="black",size=12),
          axis.text.x = element_text(angle = 90))
  if (!is.null(filename)) {
    ggsave(filename, plot)
  }
  return(plot)
}
