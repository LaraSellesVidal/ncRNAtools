\name{generatePairsProbabilityMatrix}
\alias{generatePairsProbabilityMatrix}
\title{Generates a matrix of base pair probabilities
}
\description{
Generates a matrix of base pair probabilities from a table of base pair
probabilities in the format returned by predictSecondaryStructure when using
centroidFold or centroidHomFold as the prediction method.
}
\usage{
generatePairsProbabilityMatrix(basePairProbsTable)
}
\arguments{
  \item{basePairProbsTable}{A dataframe where each line corresponds to a nucleotide
  of the query RNA sequence. The first column indicates the position number, the
  second column indicates the corresponding nucleotide type and additional columns
  indicating the probability of forming a base pair with other nucleotides. The
  potentially pairing nucleotides and their corresponding probabilities should be
  provided as strings, with a colon separating both fields.
  }
}
\value{
A symmetric square matrix with a number of rows and columns equal to the number
of nucleotides in the corresponding RNA sequence, determined by the number of
rows of basePairProbsTable. The names of rows and columns are the nucleotide
type and position. The value of each cell of the matrix is the probability that
the nucleotides at positions given by the row and columns of the cell form a
base pair.
}
\examples{
# Predict the secondary structure of an RNA sequence with centroidFold:

structurePrediction <- predictSecondaryStructure("GGGGAUGUAGCUCAUAUGGUAGAGCGCUCG
CUUUGCAUGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA", "centroidFold")

# Generate a base pair probability matrix from the returned base pair probability
# table:

basePairProbabilityMatrix <- generatePairsProbabilityMatrix(structurePrediction$basePairProbabilities)
}
