library(RUnit)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)

## Test rnaCentralTextSearch

checkTrue(grepl("^URS", rnaCentralTextSearch("HOTAIR")[1]))

## Test rnaCentralRetrieveEntry

checkEquals("URS000075C808_9606", rnaCentralRetrieveEntry("URS000075C808_9606")$rnaCentralID)

## Test rnaCentralGenomicCoordinates

checkTrue(grepl("^URS", rnaCentralGenomicCoordinatesSearch(GenomicRanges::GRanges(seqnames=S4Vectors::Rle("chr3"), ranges=IRanges::IRanges(39788104, 39795326)), "Homo sapiens")[[1]][[1]]$rnaCentralID))
