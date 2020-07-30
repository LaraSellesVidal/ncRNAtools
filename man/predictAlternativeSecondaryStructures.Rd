\name{predictAlternativeSecondaryStructures}
\alias{predictAlternativeSecondaryStructures}
\title{Predicts alternative secondary structures of a given RNA sequence
}
\description{
Attempts to identify potential alternative secondary structures of the provided
RNA sequence using the RintW method, based on the decomposition of the base-pairing
probability matrix over the Hamming distance to a reference secondary structure.
}
\usage{
predictAlternativeSecondaryStructures(sequence, gammaWeight=4, inferenceEngine="BL")
}
\arguments{
  \item{sequence}{string with an RNA sequence whose secondary structure should
  be predicted. Should contain only standard RNA symbols (i.e., "A", "U", "G"
  and "C"), and no spaces or newlines.
  }
  \item{gammaWeight}{weight factor for predicted base pairs. It directly affects
  the number of predicted base pairs. A higher value leads to a higher number of
  base pairs predicted. It should be a positive number. In the default behavior,
  a value of 4 is used.
  }
  \item{inferenceEngine}{engine used to identify the optimal canonical secondary
  structure. Possible values are "BL", "Turner" and "CONTRAfold". In the first
  two cases, a McCaskill partition function is applied, using respectively the
  Boltzmann likelihood model or Turner's energy model. In the third case, the
  CONTRAfold engine, based on conditional log-linear models, is applied. In the
  default behavior, a McCaskill partition function with a Boltzmann likelihood
  model is used.
  }
}
\value{
A list of two-element lists, where each element of the upper level list represents
a potential secondary structure. The first top-level element always represents
the canonical secondary structure. If no alternative secondary structures are
found, simply a list of two elements is returned, comprising the query sequence
and the canonical secondary structure.

When alternative secondary structure elements are found, each top-level element
comprises the following two elements:
  \item{sequence }{Query RNA sequence
  }
  \item{secondaryStructure }{Predicted secondary structure
  }
}
\references{
Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Computational approaches
for RNA energy parameter estimation. RNA. 2010;16(12):2304-2318.
doi:10.1261/rna.1950510

Do CB, Woods DA, Batzoglou S. CONTRAfold: RNA secondary structure prediction
without physics-based models. Bioinformatics. 2006;22(14):e90-e98.
doi:10.1093/bioinformatics/btl246

Hagio T, Sakuraba S, Iwakiri J, Mori R, Asai K. Capturing alternative secondary
structures of RNA by decomposition of base-pairing probabilities.
BMC Bioinformatics. 2018;19(Suppl 1):38. Published 2018 Feb 19.
doi:10.1186/s12859-018-2018-4

Hamada M, Ono Y, Kiryu H, et al. Rtools: a web server for various secondary
structural analyses on single RNA sequences. Nucleic Acids Res. 2016;44(W1):W302-W307.
doi:10.1093/nar/gkw337

Mathews DH, Sabina J, Zuker M, Turner DH. Expanded sequence dependence of
thermodynamic parameters improves prediction of RNA secondary structure.
J Mol Biol. 1999;288(5):911-940. doi:10.1006/jmbi.1999.2700

http://rtools.cbrc.jp/
}
\examples{
# Predict alternative secondary structures of an RNA sequence:

alternativeStructures <- predictAlternativeSecondaryStructures("GGGGAUGUAGCUCAUA
UGGUAGAGCGCUCGCUUUGCAUGCGAGAGGC")

# Count the number of potential alternative structures identified:

length(alternativeStructures)
}
