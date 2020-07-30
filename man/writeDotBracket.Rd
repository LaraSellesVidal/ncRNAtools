\name{writeDotBracket}
\alias{writeDotBracket}
\title{Writes a file with the sequence and secondary structure of an RNA in the
Dot-Bracket format
}
\description{
Writes a file with the sequence and secondary structure of an RNA in the
Dot-Bracket format.
}
\usage{
writeDotBracket(filename, sequence, secondaryStructure, sequenceName="Sequence")
}
\arguments{
  \item{filename}{A string indicating the path to the Dot-Bracket file to be
  written. A description of the Dot-Bracket format can be found at
  https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html.
  The resulting file will contain three lines. The first line starts with ">",
  after which the name of the sequence follows. The second line contains the RNA
  sequence. The third line contains the secondary structure representation in
  Dot-Bracket notation.
  }
  \item{sequence}{A string with the full-length sequence of the RNA for which
  a Dot-Bracket file should be written.
  }
  \item{secondaryStructure}{A string representing the secondary structure of the
  RNA sequence in the Dot-Bracket format.
  }
  \item{sequenceName}{A string with the name of the RNA sequence for which
  a Dot-Bracket file should be written.
  }
}
\value{
Called for its effect of writing a Dot-Bracket file. Invisibly returns the status
code returned by close when closing the file connection. See documentation of close
for details.
}
\references{
https://software.broadinstitute.org/software/igv/RNAsecStructure

http://projects.binf.ku.dk/pgardner/bralibase/RNAformats.html
}
\examples{
# Write a Dot-Bracket file by providing an RNA sequence and its secondary structure in
# the Dot-Bracket format:

tempDir <- tempdir()
DotBracketFile <- paste(tempDir, "testDotBracketFile.dot", sep="")
writeDotBracket(DotBracketFile, "AGCGGGUUCCUGGUUCCCCAAGGUUGA",
secondaryStructure="...(((..((.))..))).((..))..", sequenceName="Test sequence")
}
