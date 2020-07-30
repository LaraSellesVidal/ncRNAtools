\name{rnaCentralRetrieveEntry}
\alias{rnaCentralRetrieveEntry}
\title{Retrieves information for a specific entry of the RNAcentral database
}
\description{
Retrieves information for a specific entry of the RNAcentral database. The
retrieved information includes the corresponding sequence, the sequence length,
a brief description of the sequence, the species where the sequence is found,
the NCBI taxid of the species and the RNA type.
}
\usage{
rnaCentralRetrieveEntry(rnaCentralID)
}
\arguments{
  \item{rnaCentralID}{string indicating the RNAcentral ID of the entry that
  should be retrieved. It must start with "URS".
  }
}
\value{
A list containing the following elements that correspond to the RNAcentral entry
associated with the provided RNAcentral ID:
  \item{rnaCentralID }{Query RNAcentral ID
  }
  \item{sequence }{RNA sequence of the retrieved entry
  }
  \item{sequenceLength }{length of the RNA sequence of the retrieved entry
  }
  \item{description }{short description of the retrieved entry
  }
  \item{species }{species where the sequence of the retrieved entry is found
  }
  \item{ncbiTaxID }{NCBI taxonomy identifier of the species corresponding to the
  retrieved entry
  }
  \item{RNATypes }{categories of RNA associated to the RNA sequence of the
  retrieved entry
  }
}
\references{
https://rnacentral.org/about-us

https://rnacentral.org/api
}
\examples{
# Retrieve information of the RNAcentral entry associated with ID URS00007C2D83_224308:

rnaCentralEntry <- rnaCentralRetrieveEntry("URS00007C2D83_224308")

# Extract the corresponding RNA sequence

rnaCentralEntry$sequence
}
