\name{plotCompositePairsMatrix}
\alias{plotCompositePairsMatrix}
\title{Plots a composite matrix of base pair probabilities and paired bases
}
\description{
Generates a heatmap-like plot to visualize the probabilities that different
bases of an RNA molecule form a pair together with bases that are paired
in a given secondary structure for the same RNA. The top-right half of the plot
represents base pair probabilities, while in the bottom-left half paired bases
are represented with black squares in the corresponding cells. Usually, a
correspondance between paired bases and high-probability base pairs is observed.
}
\usage{
plotCompositePairsMatrix(basePairProbsMatrix, pairedBases, probabilityThreshold=0.1,
                         colorPalette=paste(rainbow(7, rev=TRUE), "FF", sep=""))
}
\arguments{
  \item{basePairProbsMatrix}{A symmetric square matrix containing the probabilities
  of pairs between different bases. Should be in the same format as output by
  the generatePairsProbabilityMatrix function.
  }
  \item{pairedBases}{A dataframe indicating paired bases, in the same format as
  output by the findPairedBases function.
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
A ggplot object with a representation of the base pair probability matrix and
bases paired in a given secondary structure.
}
\examples{
# Create a list with an RNA sequence, its secondary structure and a table of 
# base pair probabilities, calculated with centroidFold:

secondaryStructure <- list(sequence="GGGGAUGUAGCUCAUAUGGUAGAGCGCUCGCUUUGCAUGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA", 
secondaryStructure="(((((((..((((........)))).(((((.......))))).....(((((.......)))))))))))).", 
basePairProbabilities=read.csv(system.file("extdata", "exampleBasePairProbabilitiesTable.csv", package="ncRNAtools")))

# Generate a matrix of base pair probabilities:

probabilitiesMatrix <- generatePairsProbabilityMatrix(secondaryStructure$basePairProbabilities)

# Generate a dataframe with paired bases:

pairedBases <- findPairedBases(secondaryStructure$secondaryStructure, secondaryStructure$sequence)

# Plot base pair probabilities and paired bases:

plotCompositePairsMatrix(probabilitiesMatrix, pairedBases)
}
