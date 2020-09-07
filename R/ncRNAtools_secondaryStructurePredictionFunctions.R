predictSecondaryStructure <- function(sequence, method, gammaWeight=NULL, inferenceEngine=NULL,
                                      alignmentEngine=NULL, eValueRfamSearch=NULL, numHomSeqsRfamSearch=NULL) {
  sequence <- removeNewLines(sequence)
  checkRNAString(sequence)
  sendQueryResult <- sendSecondaryStructureQuery(sequence, method, gammaWeight=gammaWeight, inferenceEngine=inferenceEngine,
                                                 alignmentEngine=alignmentEngine, eValueRfamSearch=eValueRfamSearch,
                                                 numHomSeqsRfamSearch=numHomSeqsRfamSearch)
  predictionFinished <- FALSE
  while (!predictionFinished) {
    Sys.sleep(1)
    queryRunning <- checkSecondaryStructureQuery(sendQueryResult[1])
    if (!queryRunning) {
      predictionFinished <- TRUE
    }
  }
  predictionResult <- retrieveSecondaryStructureResults(sendQueryResult[1], sendQueryResult[2])
  return(predictionResult)
}

predictAlternativeSecondaryStructures <- function(sequence, gammaWeight=4, inferenceEngine="BL") {
  sequence <- removeNewLines(sequence)
  checkRNAString(sequence)
  sendQueryResult <- sendAlternativeSecondaryStructureQuery(sequence, gammaWeight=gammaWeight, inferenceEngine=inferenceEngine)
  predictionFinished <- FALSE
  while (!predictionFinished) {
    Sys.sleep(3)
    queryRunning <- checkSecondaryStructureQuery(sendQueryResult)
    if (!queryRunning) {
      predictionFinished <- TRUE
    }
  }
  predictionResult <- retrieveAlternativeSecondaryStructureResults(sendQueryResult)
  return(predictionResult)
}

generatePairsProbabilityMatrix <- function(basePairProbsTable) {
  sequenceLength <- nrow(basePairProbsTable)
  sequenceNucleotides <- paste0(basePairProbsTable[, 2], seq_len(sequenceLength))
  basePairProbsMatrix <- matrix(0, nrow=sequenceLength, ncol=sequenceLength, 
                                dimnames=list(sequenceNucleotides, sequenceNucleotides))
  for (row in seq_len(nrow(basePairProbsTable))) {
    basePairProbabilitiesList <- lapply(basePairProbsTable[row, c(-1, -2)], 
                                        FUN=extractBasePairProbability)
    basePairProbabilitiesList <- basePairProbabilitiesList[lengths(basePairProbabilitiesList) != 0]
    pairedBasesPositions <- unlist(lapply(basePairProbabilitiesList, `[[`, 1))
    basePairProbabilities <- unlist(lapply(basePairProbabilitiesList, `[[`, 2))
    basePairProbsMatrix[row, pairedBasesPositions] <- basePairProbabilities
  }
  basePairProbsMatrix[lower.tri(basePairProbsMatrix)] <- t(basePairProbsMatrix)[lower.tri(basePairProbsMatrix)]
  return(basePairProbsMatrix)
}

findPairedBases <- function(secondaryStructureString, sequence) {
  sequence <- removeNewLines(sequence)
  checkBasicDotBracketString(secondaryStructureString)
  secondaryStructureCharacters <- splitString(secondaryStructureString, split="")
  sequenceCharacters <- splitString(sequence, split="")
  pairedBases <- data.frame(Position1=integer(), Position2=integer(),
                            Nucleotide1=character(), Nucleotide2=character())
  basePairsCount <- 0
  openBracketsPositions <- integer()
  openBracketsNucleotides <- character()
  for (i in seq(1, length(secondaryStructureCharacters))) {
    if (secondaryStructureCharacters[i] == "(") {
      openBracketsPositions <- c(openBracketsPositions, i)
      openBracketsNucleotides <- c(openBracketsNucleotides, sequenceCharacters[i])
    } else if (secondaryStructureCharacters[i] == ")") {
      basePairsCount <- basePairsCount + 1
      pairedBases[basePairsCount, "Position1"] <- openBracketsPositions[length(openBracketsPositions)]
      pairedBases[basePairsCount, "Position2"] <- i
      pairedBases[basePairsCount, "Nucleotide1"] <- openBracketsNucleotides[length(openBracketsNucleotides)]
      pairedBases[basePairsCount, "Nucleotide2"] <- sequenceCharacters[i]
      openBracketsPositions <- openBracketsPositions[-length(openBracketsPositions)]
      openBracketsNucleotides <- openBracketsNucleotides[-length(openBracketsNucleotides)]
    }
  }
  pairedBases <- pairedBases[order(pairedBases$Position1),]
  rownames(pairedBases) <- NULL
  return(pairedBases[order(pairedBases$Position1),])
}

pairsToSecondaryStructure <- function(pairedBases, sequence) {
  sequence <- removeNewLines(sequence)
  secondaryStructureString <- rep(".", nchar(sequence))
  secondaryStructureString[pairedBases[, "Position1"]] <- "("
  secondaryStructureString[pairedBases[, "Position2"]] <- ")"
  return(paste(secondaryStructureString, collapse=""))
}
