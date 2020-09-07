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
# Read a Dot Bracket file with the secondary structure of an RNA sequence:

exampleDotBracketFile <- system.file("extdata", "exampleDotBracket.dot", 
                                     package="ncRNAtools")

exampleDotBracket <- readDotBracket(exampleDotBracketFile)

# Generate a dataframe of paired bases from the returned secondary structure
# string:

pairedBasesTable <- findPairedBases(secondaryStructureString=exampleDotBracket$secondaryStructure,
sequence=exampleDotBracket$sequence)

# Generate a secondary structure string from the table of paired bases:

secondaryStructureString <- pairsToSecondaryStructure(pairedBasesTable,
exampleDotBracket$sequence)

# Verify that the resulting secondary structure string is equal to the original
# prediction:

secondaryStructureString == exampleDotBracket$secondaryStructure
}
