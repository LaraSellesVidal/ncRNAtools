\name{flattenDotBracket}
\alias{flattenDotBracket}
\title{Converts an extended Dot-Bracket secondary structure to basic Dot-Bracket
notation
}
\description{
Generates a string with the secondary structure of an RNA sequence in the basic
Dot-Bracket notation from a string in the extended Dot-Bracket notation.
}
\usage{
flattenDotBracket(extendedDotBracketString)
}
\arguments{
  \item{extendedDotBracketString}{A string with a secondary structure representation
  of an RNA molecule in the extended Dot-Bracket notation. The extended Dot-Bracket
  notation uses dots (".") for unpaired bases, and multiple pairs of symbols ("("-")",
  "["-"]", "{"-"}", "<"-">", "A"-"a", "B"-"b", "C"-"c" and "D"-"d") to indicate
  paired bases. This allows the format to unambiguously represent highly nested
  structures, including pseudoknots.
  }
}
\value{
A string representing the secondary structure of the provided RNA in the basic
Dot-Bracket format. In such format, dots (".") indicate unpaired bases, and
paired bases are represented with pairs of round brackets ("("-")").
}
\references{
https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html
}
\examples{
# A secondary structure string with characters other than dots and round 
# brackets, representing unambiguously pseudo-knots:

extendedDotBracketString <- "((((....((((.((((([[[[...)))))]]]]..)))).....))))."

# Convert it to the basic Dot-Bracket format:

basicDotBracketString <- flattenDotBracket(extendedDotBracketString)
}
