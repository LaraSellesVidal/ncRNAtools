\name{plotPairsProbabilityMatrix}
\alias{plotPairsProbabilityMatrix}
\title{Plots a matrix of base pair probabilities
}
\description{
Generates a heatmap-like plot to visualize the probabilities that different
bases of an RNA molecule form a pair.
}
\usage{
plotPairsProbabilityMatrix(basePairProbsMatrix, probabilityThreshold=0.1,
                           colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep=""))
}
\arguments{
  \item{basePairProbsMatrix}{A symmetric square matrix containing the probabilities
  of pairs between different bases. Should be in the same format as output by
  the generatePairsProbabilityMatrix function.
  }
  \item{probabilityThreshold}{Threshold for representing the probability that
  two given bases form a pair. Pairs with a probability lower than the threshold
  will not be considered, and their corresponding cell in the plot will be left
  blank.
  }
  \item{colorPalette}{Color palette to be used for displaying the probabilities
  above the specified threshold.
  }
}
\value{
A ggplot object with a representation of the base pair probability matrix.
}
\examples{
# Load an example table of base pair probabilities, calculated with centroidFold:

basePairProbabilitiesTable <- read.csv(system.file("extdata", 
"exampleBasePairProbabilitiesTable.csv", package="ncRNAtools"))

# Generate a matrix of base pair probabilities:

probabilitiesMatrix <- generatePairsProbabilityMatrix(basePairProbabilitiesTable)

# Plot the probability matrix

plotPairsProbabilityMatrix(probabilitiesMatrix)
}
