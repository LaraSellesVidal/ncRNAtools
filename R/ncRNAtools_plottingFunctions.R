plotPairsProbabilityMatrix <- function(basePairProbsMatrix, probabilityThreshold=0.1,
                                       colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep="")) {
  meltedData <- data.frame(nucleotide1=factor(rownames(basePairProbsMatrix)[row(basePairProbsMatrix)],
                                              levels=rownames(basePairProbsMatrix)),
                           nucleotide2=factor(colnames(basePairProbsMatrix)[col(basePairProbsMatrix)],
                                              levels=colnames(basePairProbsMatrix)),
                           probability=as.vector(basePairProbsMatrix), 
                           row.names=NULL, stringsAsFactors=TRUE)
  plot <- ggplot(data = meltedData, aes(x=nucleotide1, y=nucleotide2, fill=probability)) +
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
                         values=c(0, seq(probabilityThreshold, 1, 
                                         length.out=length(colorPalette))),
                         breaks=c(0, 0.5, 1),
                         labels=c(0, 0.5, 1),
                         limits=c(0, 1),
                         name="Probability") +
    scale_x_discrete(position="top",
                     breaks=levels(meltedData$nucleotide1)[seq(10, nrow(basePairProbsMatrix), by=10)]) +
    scale_y_discrete(breaks=levels(meltedData$nucleotide1)[seq(10, nrow(basePairProbsMatrix), by=10)],
                     limits=rev(levels(meltedData$nucleotide2))) +
    coord_equal() +
    theme(panel.border=element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill = 'white', colour="white"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(color="black",size=12),
          axis.text.x = element_text(angle = 90))
  return(plot)
}

plotCompositePairsMatrix <- function(basePairProbsMatrix, pairedBases, probabilityThreshold=0.1,
                                     colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep="")) {
  compositeMatrix <- makeCompositeMatrix(basePairProbsMatrix, pairedBases)
  compositeMatrix[upper.tri(compositeMatrix) & (compositeMatrix == 1)] <- NA
  meltedData <- data.frame(nucleotide1=factor(rownames(compositeMatrix)[row(compositeMatrix)],
                                              levels=rownames(compositeMatrix)),
                           nucleotide2=factor(colnames(compositeMatrix)[col(compositeMatrix)],
                                              levels=colnames(compositeMatrix)),
                           probability=as.vector(compositeMatrix), 
                           row.names=NULL, stringsAsFactors=TRUE)
  plot <- ggplot(data = meltedData, aes(x=nucleotide1, y=nucleotide2, fill=probability)) +
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
                         values=c(0, seq(probabilityThreshold, 1, 
                                         length.out=length(colorPalette))),
                         breaks=c(0, 0.5, 1),
                         labels=c(0, 0.5, 1),
                         limits=c(0, 1),
                         name="Probability",
                         na.value="black") +
    scale_x_discrete(position="top",
                     breaks=levels(meltedData$nucleotide1)[seq(10, nrow(basePairProbsMatrix), by=10)]) +
    scale_y_discrete(breaks=levels(meltedData$nucleotide1)[seq(10, nrow(basePairProbsMatrix), by=10)],
                     limits=rev(levels(meltedData$nucleotide2))) +
    coord_equal() +
    theme(panel.border=element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill = 'white', colour="white"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text=element_text(color="black",size=12),
          axis.text.x = element_text(angle = 90))
  return(plot)
}