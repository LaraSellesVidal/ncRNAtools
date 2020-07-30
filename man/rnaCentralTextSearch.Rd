\name{rnaCentralTextSearch}
\alias{rnaCentralTextSearch}
\title{Searches the RNAcentral database with a text query
}
\description{
Searches the RNAcentral database of non-coding RNA with a text query.
}
\usage{
rnaCentralTextSearch(query)
}
\arguments{
  \item{query}{string to be used as the text query to search the RNAcentral
  database. The most basic usage involves providing a word or group of words to
  be matched against any of the fields of entries of the RNAcentral database.
  More refined searches can be done by using the field_name:"field value" syntax.
  Several search terms can be combined with logical operators (and, or, not).
  Asterisks ("*") can be added to signal any number of any characters. If not
  added to the query, only exact matches are produced.
  A detailed description of the syntax can be found at https://rnacentral.org/help/text-search.
  The entire query should be provided as a single string, and double quotes
  enclosing field values must be escaped with a backslash (``\``).
  It should also be noted that logical operators should be capitalized.
  }
}
\value{
A character vector where each element is the RNAcentral ID of an entry that
matched the query. A maximum of 15 IDs will be returned.
}
\references{
https://rnacentral.org/about-us

https://rnacentral.org/help/text-search
}
\examples{
# Find RNAcentral entries that correspond to the HOTAIR long noncoding RNA:

rnaCentralTextSearch("HOTAIR")

# Find RNAcentral entries that correspond to FMN riboswitches in Bacillus subtilis:

rnaCentralTextSearch("FMN AND species:\"Bacillus subtilis\"")
}
