\name{pairsToSecondaryStructure}
\alias{pairsToSecondaryStructure}
\title{Generates a string with the secondary structure of an RNA sequence
}
\description{
Generates a string with the secondary structure of an RNA sequence from a
table of paired bases.
}
\usage{
pairsToSecondaryStructure(pairedBases, sequence)
}
\arguments{
  \item{pairedBases}{A dataframe where each row contains the information of two
  bases that form a pair. The dataframe should contain columns named "Position1"
  and "Position2" indicating respectively the positions of the  5' and 3' bases
  involved in the base pair.
  }
  \item{sequence}{string with the RNA sequence corresponding to the provided
  table of paired bases. Should contain only standard RNA symbols (i.e., "A",
  "U", "G" and "C"), and no spaces or newlines.
  }
}
\value{
A string representing the secondary structure of the provided RNA in the Dot-Bracket
format.
}
\examples{
# Predict the secondary structure of an RNA sequence with centroidFold:

structurePrediction <- predictSecondaryStructure("GGGGAUGUAGCUCAUAUGGUAGAGCGCUCG
CUUUGCAUGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA", "centroidFold")

# Generate a dataframe of paired bases from the returned secondary structure
# string:

pairedBasesTable <- findPairedBases(sequence=structurePrediction$sequence,
secondaryStructureString=structurePrediction$secondaryStructure)

# Generate a secondary structure string from the table of paired bases:

secondaryStructureString <- pairsToSecondaryStructure(pairedBasesTable, "GGGGAUG
UAGCUCAUAUGGUAGAGCGCUCGCUUUGCAUGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA")

# Verify that the resulting secondary structure string is equal to the original
# prediction:

secondaryStructureString == structurePrediction$secondaryStructure
}
