\name{writeCT}
\alias{writeCT}
\title{Writes a file with the secondary structure of an RNA in the CT format
}
\description{
Writes a file with the secondary structure of an RNA in the CT format
(Connectivity Table).
}
\usage{
writeCT(filename, sequence, secondaryStructure=NULL, sequenceName="Sequence", pairedBases=NULL)
}
\arguments{
  \item{filename}{A string indicating the path to the CT file to be written. A
  description of the CT format can be found at
  http://rna.urmc.rochester.edu/Text/File_Formats.html#CT.
  }
  \item{sequence}{A string with the full-length sequence of the RNA for which
  a CT file should be written.
  }
  \item{secondaryStructure}{A string representing the secondary structure of the
  RNA sequence in the Dot-Bracket format. A correct value must be provided for
  either secondaryStructure or pairedBases. If a table of paired bases is provided
  through pairedBases, the string provided through secondaryStructure will be
  ignored. If no table of paired bases is provided, one will be internally
  calculated from the string representation of the secondary structure.
  }
  \item{sequenceName}{A string with the name of the sequence to be written in
  the first line of the CT file.
  }
  \item{pairedBases}{A table of paired bases in the same format as output by the
  findPairedBases function.
  }
}
\value{
Called for its effect of writing a CT file. Invisibly returns the status code
returned by close when closing the file connection. See documentation of close
for details.
}
\references{
http://rna.urmc.rochester.edu/Text/File_Formats.html#CT.
}
\examples{
# Write a CT file by providing an RNA sequence and its secondary structure in
# the Dot-Bracket format:

tempDir <- tempdir()
CTFile <- paste(tempDir, "testCTFile.ct", sep="")
writeCT(CTFile, "AGCGGGUUCCUGGUUCCCCAAGGUUGA", 
secondaryStructure="...(((..((.))..))).((..))..",
sequenceName="Test sequence")
}
