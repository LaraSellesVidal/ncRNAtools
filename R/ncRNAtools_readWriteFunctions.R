readCT <- function(filename, sequence=NULL) {
  inputFile <- file(filename, "r")
  firstLine <- readLines(inputFile, n=1)
  connectivityTable <- read.table(inputFile)
  close(inputFile)
  firstLineFields <- splitString(firstLine, split="\\s", perl=TRUE)
  sequenceLength <- as.integer(firstLineFields[1])
  sequenceName <- firstLineFields[length(firstLineFields)]
  reconstructedSequence <- paste(connectivityTable[, 2], collapse="")
  if (nchar(reconstructedSequence) != sequenceLength) {
    if (is.null(sequence) || !(is.character(sequence) & nchar(sequence) == sequenceLength)) {
      stop(strwrap(prefix="\n", initial="", "CT file does not contain information
                   for all nucleotides. and a valid full-length sequence has not
                   been provided. Please provide either a complete CT file, or a
                   CT file with information only for paired bases and the
                   corresponding full-length RNA sequence."))
    }
    reconstructedSequence <- sequence
  }
  pairsTable <- connectivityTable[(connectivityTable[, 5] != 0) & (connectivityTable[, 5] > connectivityTable[, 1]), ]
  pairsTable <- pairsTable[, c(1, 5, 2)]
  complementaryNucleotides <- splitString(sequence, split="")[pairsTable[, 2]]
  pairsTable <- cbind(pairsTable, complementaryNucleotides)
  colnames(pairsTable) <- c("Position1", "Position2", "Nucleotide1", "Nucleotide2")
  resultsList <- list(sequenceName=sequenceName, sequence=sequence, 
                      sequenceLength=sequenceLength, pairsTable=pairsTable)
  return(resultsList)
}

readDotBracket <- function(filename) {
  inputFile <- file(filename, "r")
  firstLine <- readLines(inputFile, n=1)
  if (startsWith(firstLine, ">")) {
    sequenceName <- substring(firstLine, 2)
    sequence <- readLines(inputFile, n=1)
    structureLine <- readLines(inputFile, n=1)
    secondaryStructure <- splitString(structureLine, " ")[1]
    freeEnergy <- as.numeric(gsub("[^-+[:digit:]\\.]", "",
                                  splitString(structureLine, split=" ")[2]))
  }
  else {
    sequenceName <- NA
    sequence <- firstLine
    structureLine <- readLines(inputFile, n=1)
    secondaryStructure <- splitString(structureLine, " ")[1]
    freeEnergy <- as.numeric(gsub("[^-+[:digit:]\\.]", "",
                                  splitString(structureLine, " ")[2]))
  }
  close(inputFile)
  resultsList <- list(sequenceName=sequenceName, sequence=sequence,
                      secondaryStructure=secondaryStructure, freeEnergy=freeEnergy)
  return(resultsList)
}

writeCT <- function(filename, sequence, secondaryStructure=NULL, sequenceName="Sequence", pairedBases=NULL) {
  if (is.null(secondaryStructure) & is.null(pairedBases)) {
    stop("Please provide a secondary structure string or dataframe of paired bases.")
  } else if (!is.null(secondaryStructure) & is.null(pairedBases)) {
    pairsTable <- findPairedBases(secondaryStructure, sequence)
  } else {
    pairsTable <- pairedBases
  }
  outputFile <- file(filename, "w")
  writeLines(paste(nchar(sequence), sequenceName, sep="\t"), outputFile)
  for (i in seq_len(nchar(sequence))) {
    if (i %in% pairsTable[, 1]) {
      rowToWrite <- which(pairsTable[, 1] == i)
      writeLines(paste(i, substr(sequence, i, i), i-1, i+1, 
                       pairsTable[rowToWrite, 2], i, sep="\t"), 
                 outputFile)
    } else if (i %in% pairsTable[, 2]) {
      rowToWrite <- which(pairsTable[, 2] == i)
      writeLines(paste(i, substr(sequence, i, i), i-1, i+1, 
                       pairsTable[rowToWrite, 1], i, sep="\t"), 
                 outputFile)
    } else {
      writeLines(paste(i, substr(sequence, i, i), i-1, i+1, 0, i, sep="\t"), 
                 outputFile)
    }
  }
  close(outputFile)
}

writeDotBracket <- function(filename, sequence, secondaryStructure, 
                            sequenceName="Sequence") {
  outputFile <- file(filename, "w")
  writeLines(paste(">", sequenceName, sep=""), outputFile)
  writeLines(sequence, outputFile)
  writeLines(secondaryStructure, outputFile)
  close(outputFile)
}
