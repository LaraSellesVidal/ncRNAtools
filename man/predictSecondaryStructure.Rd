\name{predictSecondaryStructure}
\alias{predictSecondaryStructure}
\title{Predicts the secondary structure of a given RNA sequence
}
\description{
Predicts the secondary structure of the provided RNA sequence, using the chosen
prediction method. Secondary structure predictions are carried out with the
rtools RNA Bioinformatics Web server.
}
\usage{
predictSecondaryStructure(sequence, method, gammaWeight=NULL, inferenceEngine=NULL,
alignmentEngine=NULL, eValueRfamSearch=NULL, numHomSeqsRfamSearch=NULL)
}
\arguments{
  \item{sequence}{string with an RNA sequence whose secondary structure should
  be predicted. Should contain only standard RNA symbols (i.e., "A", "U", "G"
  and "C").
  }
  \item{method}{method that should be used for the prediction of RNA secondary
  structure. Possible values are "centroidFold", "centroidHomFold" and "IPknot".
  Only IPknot is able to predict pseudoknots. For a detailed description of each
  method, see respectively \emph{Hamada et al., 2008}; \emph{Hamada et al., 2009},
  and \emph{Sato et al., 2011}.
  }
  \item{gammaWeight}{weight factor for predicted base pairs. It directly affects
  the number of predicted base pairs. A higher value leads to a higher number of
  base pairs predicted. It should be a positive number. In the default behavior,
  when no specific value is provided, the default value for each secondary
  structure prediction method in the rtools webserver is used (4 for centroidFold
  and IPknot, and 8 forcentroidHomFold).
  }
  \item{inferenceEngine}{engine used to identify optimal secondary structures.
  Possible values are "BL", "Turner" and "CONTRAfold". In the first two cases, a
  McCaskill partition function is applied, using respectively the Boltzmann
  likelihood model or Turner's energy model. In the third case, the CONTRAfold
  engine, based on conditional log-linear models, is applied. Additionally, if
  IPknot is chosen as the method for secondary structure prediction, "NUPACK" is
  also a possible value if the sequence has 100 nucleotides or less. In this case,
  the NUPACK scoring model is used. In the default behavior, when no specific
  value is provided, the default inference engine in the rtools webserver is,
  which is a McCaskill partition function with a Boltzmann likelihood model.
  }
  \item{alignmentEngine}{engine used to perform pairwise alignments of the query
  sequence and Rfam homologous sequences during the application of centroidHomFold.
  Possible values are "CONTRAlign" and "ProbCons". For details on each alignment
  engine, see \emph{Do et al., 2006} and \emph{Do et al., 2005} respectively.
  In the default behavior, when no value is specified, CONTRAlign is used.
  }
  \item{eValueRfamSearch}{e-value used to select homologous sequences from the
  Rfam database during the application of centroidHomFold. Should be a number
  equal to or greater than 0. In the default behavior, when no value is specified,
  a value of 0.01 is used.
  }
  \item{numHomSeqsRfamSearch}{maximum number of homologous sequences to be considered
  during the application of centroidHomFold. Should be a positive integer.In the
  default behavior, when no value is specified, a value of 30 is used.
  }
}
\value{
If either centroidFold or centroidHomFold are used for predicting secondary
structure, a list of three elements comprising the query RNA sequence, the predicted
secondary structure and a table of base pair probabilities. The secondary structure
is represented as a string in the Dot-Bracket format. The three elements of the list
are:
  \item{sequence }{Query RNA sequence
  }
  \item{secondaryStructure }{Predicted secondary structure
  }
  \item{basePairProbabilities }{A dataframe where each line corresponds to a nucleotide
  of the query RNA sequence. The first column indicates the position number, the
  second column indicates the corresponding nucleotide type and additional columns
  indicating the probability of forming a base pair with other nucleotides. The
  potentially pairing nucleotides and their corresponding probabilities are provided
  as strings, where a colon separates both fields
  }

If IPknot is used for predicting secondary structure, no table of base pair
probabilities is returned. Therefore, the output is a list of only two elements
(sequence and secondaryStructure). Additionally, the secondary structure is
provided in the extended Dot-Bracket format if required to represent pseudoknots
unambiguously.
}
\references{
Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Computational approaches
for RNA energy parameter estimation. RNA. 2010;16(12):2304-2318.
doi:10.1261/rna.1950510

Do C.B., Gross S.S., Batzoglou S. CONTRAlign: Discriminative Training for
Protein Sequence Alignment. In: Apostolico A., Guerra C., Istrail S., Pevzner P.A.,
Waterman M. (eds) Research in Computational Molecular Biology. 2006.
Lecture Notes in Computer Science, vol 3909. Springer, Berlin, Heidelberg
doi:10.1007/11732990_15

Do CB, Mahabhashyam MS, Brudno M, Batzoglou S. ProbCons: Probabilistic
consistency-based multiple sequence alignment. Genome Res. 2005;15(2):330-340.
doi:10.1101/gr.2821705

Do CB, Woods DA, Batzoglou S. CONTRAfold: RNA secondary structure prediction
without physics-based models. Bioinformatics. 2006;22(14):e90-e98.
doi:10.1093/bioinformatics/btl246

Hamada M, Kiryu H, Sato K, Mituyama T, Asai K. Prediction of RNA secondary
structure using generalized centroid estimators. Bioinformatics. 2009;25(4):465-473.
doi:10.1093/bioinformatics/btn601

Hamada M, Ono Y, Kiryu H, et al. Rtools: a web server for various secondary
structural analyses on single RNA sequences. Nucleic Acids Res. 2016;44(W1):W302-W307.
doi:10.1093/nar/gkw337

Hamada M, Yamada K, Sato K, Frith MC, Asai K. CentroidHomfold-LAST: accurate
prediction of RNA secondary structure using automatically collected homologous
sequences. Nucleic Acids Res. 2011;39(Web Server issue):W100-W106.
doi:10.1093/nar/gkr290

Mathews DH, Sabina J, Zuker M, Turner DH. Expanded sequence dependence of
thermodynamic parameters improves prediction of RNA secondary structure.
J Mol Biol. 1999;288(5):911-940. doi:10.1006/jmbi.1999.2700

Sato K, Kato Y, Hamada M, Akutsu T, Asai K. IPknot: fast and accurate prediction
of RNA secondary structures with pseudoknots using integer programming.
Bioinformatics. 2011;27(13):i85-i93. doi:10.1093/bioinformatics/btr215

http://rtools.cbrc.jp/
}
\examples{
# Predict the secondary structure of an RNA sequence with IPknot:

structurePrediction <- predictSecondaryStructure("UGCGAGAGGCACAGGGUUCGAUUCCCUGCA
UCUCCA", "IPknot", inferenceEngine = "NUPACK")

# Extract the string representing the secondary structure prediction:

structurePrediction$secondaryStructure
}
