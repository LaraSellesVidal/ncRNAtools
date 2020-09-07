\name{readCT}
\alias{readCT}
\title{Reads a file with the secondary structure of an RNA in the CT format
}
\description{
Reads a file with the secondary structure of an RNA in the CT format
(Connectivity Table).
}
\usage{
readCT(filename, sequence=NULL)
}
\arguments{
  \item{filename}{A string indicating the path to the CT file to be read. A
  description of the CT format can be found at
  http://rna.urmc.rochester.edu/Text/File_Formats.html#CT.
  The CT file should contain information for all nucleotides of the RNA
  sequence, including those that do not form base pairs, unless the full-length
  RNA sequence is provided through the sequence argument. In such case, the CT
  file may contain only information for nucleotides involved in base pairs.
  }
  \item{sequence}{A string with the full-length sequence of the RNA whose
  secondary structure is represented in the CT file. Such argument is optional,
  but if it is not provided, a complete CT file with information for all
  nucleotides of the RNA must be provided.
  }
}
\value{
A list with the following 4 elements:
  \item{sequenceName }{name of the RNA sequence
  }
  \item{sequence }{RNA sequence
  }
  \item{sequenceLength }{number of nucleotides of the RNA sequence
  }
  \item{pairsTable }{a dataframe indicating the paired bases, in the same format
  as that returned by the findPairedBases function
  }
}
\references{
http://rna.urmc.rochester.edu/Text/File_Formats.html#CT.
}
\examples{
exampleCTFile <- system.file("extdata", "exampleCT.ct", package="ncRNAtools")

# If the CT file does not contain information for unpaired nucleotides, the
# original sequence must also be supplied in order to read it (tmRNA of E. coli 
# encoded by the ssrA gene):

tmRNASequence <- "GGGGCUGAUUCUGGAUUCGACGGGAUUUGCGAAACCCAAGGUGCAUGCCGAGGGGCGGUUGG
CCUCGUAAAAAGCCGCAAAAAAUAGUCGCAAACGACGAAAACUACGCUUUAGCAGCUUAAUAACCUGCUUAGAGCCCUCU
CUCCCUAGCCUCCGCUCUUAGGACGGGGAUCAAGAGAGGUCAAACCCAAAAGAGAUCGCGUGGAAGCCCUGCCUGGGGUU
GAAGCGUUAAAACUUAAUCAGGCUAGUUUGUUAGUGGCGUGUCCGUCCGCAGCUGGCAAGCGAAUGUAAAGACUGACUAA
GCAUGUAGUACCGAGGAUGUAGGAAUUUCGGACGCGGGUUCAACUCCCGCCAGCUCCACCA"

tmRNASecondaryStructure <- readCT(exampleCTFile, tmRNASequence)
}
