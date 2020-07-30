\name{rnaCentralGenomicCoordinatesSearch}
\alias{rnaCentralGenomicCoordinatesSearch}
\title{Retrieves annotated non-coding RNA in a set of genomic ranges
}
\description{
Retrieves RNAcentral entries corresponding to non-coding RNA present within
a specified set of genomic coordinates
}
\usage{
rnaCentralGenomicCoordinatesSearch(genomicRanges, species)
}
\arguments{
  \item{genomicRanges}{GRanges object specifying the genomic coordinates for
  which known non-coding RNA should be retrieved. Each sequence name must be of
  the format chr? or simply ?, where ? can be any number or the "X" or "Y"
  characters. These denote, respectively, standard somatic chromosomes and
  sexual chromosomes. Additionally, sequence name can also be "MT", which refers
  to the mitochondrial DNA. Multiple ranges can be provided in the same GRanges
  object, but they must all belong to the same organism.
  }
  \item{species}{string specifying the species to which the ranges supplied with
  the genomicRanges argument refer. The scientific name should be provided. The
  list of species for which search of RNAcentral by genomic coordinates is
  supported can be found at https://rnacentral.org/help/genomic-mapping.
  }
}
\value{
A nested list whose length matches the number of genomic ranges provided through
genomicRanges. Each top-level element is a list containing all hits found for
a given genomic range. In turn, each of the elements of such list representing
a hit is a list comprising the following elements that describe the hit:
  \item{rnaCentralID }{RNAcentral ID of the non-coding RNA that was identified
  as a hit
  }
  \item{species }{species whose genome was searched
  }
  \item{description }{brief description of the RNA sequence of the hit
  }
  \item{RNAType }{category of the RNA sequence of the hit
  }
  \item{genomicCoordinates }{GRanges object containing a single range, within
  which the hit was found
  }
}
\references{
https://rnacentral.org/help/genomic-mapping

https://rnacentral.org/api
}
\examples{
# Generate a GRanges object with 2 genomic ranges specifying coordinates of
# the human genome:

genomicCoordinates <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(c("chr3", "chr4")),
ranges=IRanges::IRanges(rep(39745816, 2), rep(39847679, 2)))

# Retrieve known annotated non-coding RNA present in the specified genomic
# ranges:

knownNonCodingRNA <- rnaCentralGenomicCoordinatesSearch(genomicCoordinates, "Homo sapiens")
}
