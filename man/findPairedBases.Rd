\name{findPairedBases}
\alias{findPairedBases}
\title{Determines paired bases from secondary structure string
}
\description{
Determines bases that form pairs in a given RNA sequence from the secondary
structure string of the RNA sequence.
}
\usage{
findPairedBases(secondaryStructureString, sequence)
}
\arguments{
  \item{secondaryStructureString}{A string representing the secondary structure
  of the RNA sequence in the Dot-Bracket format.
  }
  \item{sequence}{string with the RNA sequence corresponding to the provided
  secondary structure string. Should contain only standard RNA symbols (i.e.,
  "A", "U", "G" and "C"), and no spaces or newlines.
  }
}
\value{
A dataframe where each row denotes a base pair. The dataframe comprises the
following 4 columns:
  \item{Position1}{position in the sequence of the nucleotide involved in the
  base pair closest to the 5' end of the RNA.
  }
  \item{Position2}{position in the sequence of the nucleotide involved in the
  base pair closest to the 3' end of the RNA.
  }
  \item{Nucleotide1}{nucleotide type of the base in Position1.
  }
  \item{Nucleotide2}{nucleotide type of the base in Position2.
  }
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
}
