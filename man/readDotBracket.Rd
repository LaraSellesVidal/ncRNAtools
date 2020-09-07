\name{readDotBracket}
\alias{readDotBracket}
\title{Reads a file with the sequence and secondary structure of an RNA in the
Dot-Bracket format
}
\description{
Reads a file with the secondary structure of an RNA in the CT format
(Connectivity Table).
}
\usage{
readDotBracket(filename)
}
\arguments{
  \item{filename}{A string indicating the path to the Dot-Bracket file to be
  read. A description of the Dot-Bracket format can be found at
  https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html.
  Both basic and extended Dot-Bracket notations are accepted. The Dot-Bracket
  file should always contain at least two lines: the first one with the RNA
  sequence, and the second one with the secondary structure representation. The
  free energy of the secondary structure can be optionally present at the end of
  the secondary structure line, in kcal/mol.
  An additional optional line can precede these two, containing ">" as its first
  character followed by the name or description of the sequence.
  }
}
\value{
A list with the following 4 elements:
  \item{sequenceName }{name of the RNA sequence. Will be set to NA if no string
  starting with ">" is present at the beginning of the Dot-Bracket file
  }
  \item{sequence }{RNA sequence
  }
  \item{secondaryStructure }{string representing the secondary structure of the
  RNA in the Dot-Bracket format
  }
  \item{freeEnergy }{free energy of the secondary structure of the RNA, usually
  in kcal/mol. If not present in the Dot-Bracket file, it will be set to NA
  }
}
\references{
https://software.broadinstitute.org/software/igv/RNAsecStructure

http://projects.binf.ku.dk/pgardner/bralibase/RNAformats.html
}
\examples{
exampleDotBracketFile <- system.file("extdata", "exampleDotBracket.dot", 
                                     package="ncRNAtools")

exampleDotBracket <- readDotBracket(exampleDotBracketFile)

# Since the file does not contain an initial line with the name of the sequence,
# sequenceName is set to NA
}
