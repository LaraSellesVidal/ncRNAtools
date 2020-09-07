rnaCentralTextSearch <- function(query) {
  if (!is.character(query)) {
    stop("Input query is not a string")
  }
  result <- GET(rnaCentralEbiApiURL,
                query=list(query=query,
                           format="idlist"))
  matchList <- splitString(content(result, as="text"), split="\n")
  return(matchList[grep("^URS", matchList)])
}

rnaCentralRetrieveEntry <- function(rnaCentralID) {
  checkMultipleQuery(rnaCentralID)
  if (!grep("^URS", rnaCentralID)) {
    stop("Invalid RNAcentral ID provided.")
  }
  result <- GET(paste(rnaCentralApiURL, "/", rnaCentralID, sep=""),
                accept_json())
  resultContent <- content(result)
  parsedResult <- list(rnaCentralID=resultContent$rnacentral_id,
                       sequence=resultContent$sequence,
                       sequenceLength=resultContent$length,
                       description=resultContent$description,
                       species=resultContent$species,
                       ncbiTaxID=resultContent$taxid,
                       RNATypes=unlist(resultContent$ncrna_types))
  return(parsedResult)
}

rnaCentralGenomicCoordinatesSearch <- function(genomicRanges, species) {
  if (!is(genomicRanges, "GRanges")) {
    stop("Please provide a valid GRanges object with genomic coordinates.")
  }
  if (length(strsplit(species, split=" ")[[1]]) != 2) {
    stop("Please enter a valid species name. See documentation for details.")
  }
  speciesString <- gsub(" ", "_", tolower(species))
  chromosomes <- gsub("chr", "", as.character(seqnames(genomicRanges)))
  if (any(grepl("[^[:alnum:]]", chromosomes) & chromosomes != "MT")) {
    stop(strwrap("Please enter valid chromosome names as seqnames for the GRanges
                 object. See documentation for details.", initial="", prefix="\n"))
  }
  startPoints <- start(genomicRanges)
  endPoints <- end(genomicRanges)
  annotatedRNA <- vector(mode="list", length=length(genomicRanges))
  for (i in seq_len(length(genomicRanges))) {
    result <- GET(paste(rnaCentralRangeSearchURL, "/", speciesString, "/", chromosomes[i],
                        ":", startPoints[i], "-", endPoints[i], sep=""),
                  accept_json())
    parsedResult <- content(result)
    parsedResult <- parsedResult[grepl("^URS", unlist(lapply(parsedResult, `[`, "external_name")))]
    annotatedRNA[[i]] <- lapply(parsedResult, function(hit) list(rnaCentralID=splitString(hit$ID, split="@")[1],
                                                                 species=species,
                                                                 description=hit$description,
                                                                 RNAType=hit$biotype,
                                                                 genomicCoordinates=GRanges(seqnames=Rle(chromosomes[i]),
                                                                                            ranges=IRanges(hit$start, hit$end),
                                                                                            strand=splitString(hit$ID, split=":")[2])))
  }
  if (length(genomicRanges == 1)) {
    annotatedRNA <- unlist(annotatedRNA, recursive=FALSE)
  }
  return(annotatedRNA)
}