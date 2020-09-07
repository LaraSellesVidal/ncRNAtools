sendSecondaryStructureQuery <- function(sequence, method, gammaWeight, inferenceEngine,
                                        alignmentEngine, eValueRfamSearch, numHomSeqsRfamSearch) {
  registerURL <- paste(rtoolsBaseURL, "register.cgi", sep="")
  methodCode <- match(method, c("centroidFold", "centroidHomFold", "IPknot"))
  if (is.na(methodCode)) {
    stop("Invalid method for secondary structure prediction")
  }
  formData <- list(query=paste(">", "Sequence", "\n", sequence, sep="", collapse=""),
                   methodSelection="checked",
                   model_str_1="CONTRAfold",
                   gamma_1="4",
                   model_str_2="BL",
                   model_aln_2="CONTRAlign",
                   gamma_2="8",
                   eval_2="0.01",
                   homnum_2="30",
                   model_str_3="BL",
                   gamma_3="4")
  if (!is.null(gammaWeight)) {
    checkGammaWeight(gammaWeight)
    formData[c("gamma_1", "gamma_2", "gamma_3")] <- as.character(gammaWeight)
  }
  if (!is.null(inferenceEngine)) {
    checkInferenceEngine(inferenceEngine, method, sequence)
    formData[c("model_str_1", "model_str_2", "model_str_3")] <- inferenceEngine
  }
  if (!is.null(alignmentEngine)) {
    checkAlignmentEngine(alignmentEngine, method)
    formData["model_aln_2"] <- alignmentEngine
  }
  if (!is.null(eValueRfamSearch)) {
    checkEValueRfamSearch(eValueRfamSearch, method)
    formData["eval_2"] <- as.character(eValueRfamSearch)
  }
  if (!is.null(numHomSeqsRfamSearch)) {
    checkNumHomSeqsRfamSearch(numHomSeqsRfamSearch, method)
    formData["homnum_2"] <- as.character(numHomSeqsRfamSearch)
  }
  names(formData)[2] <- paste("exec_sel_", methodCode, sep="")
  response <- POST(registerURL,
                   body=formData,
                   encode="form",
                   add_headers("Host"="rtools.cbrc.jp",
                               "Content-Type"="application/x-www-form-urlencoded",
                               "Upgrade-Insecure-Requests"="1",
                               "Accept"="*/*"))
  redirectURL <- response$all_headers[[1]]$headers$location
  requestID <- splitString(redirectURL, split="req_id=")[2]
  return(c(requestID, methodCode))
}

sendAlternativeSecondaryStructureQuery <- function(sequence, gammaWeight, inferenceEngine) {
  registerURL <- paste(rtoolsBaseURL, "register.cgi", sep="")
  formData <- list(query=paste(">", "Sequence", "\n", sequence, sep="", collapse=""),
                   exec_sel_8="checked",
                   model_str_8=inferenceEngine,
                   gamma_8=gammaWeight)
  response <- POST(registerURL,
                   body=formData,
                   encode="form",
                   add_headers("Host"="rtools.cbrc.jp",
                               "Content-Type"="application/x-www-form-urlencoded",
                               "Upgrade-Insecure-Requests"="1",
                               "Accept"="*/*"))
  redirectURL <- response$all_headers[[1]]$headers$location
  requestID <- splitString(redirectURL, split="req_id=")[2]
  return(requestID)
}

checkSecondaryStructureQuery <- function(requestID) {
  secondaryStructureURL <- paste(rtoolsBaseURL, "cgi-bin/result.cgi?req_id=", 
                                 requestID, sep="")
  queryRunning <- !any(grepl("adding", names(GET(secondaryStructureURL)$headers)))
  if (queryRunning) {
    message("Secondary structure prediction is running, please wait.")
    return(queryRunning)
  }
  else if (!queryRunning) {
    message("Secondary structure prediction completed.")
    return(queryRunning)
  }
  else {
    stop("Malformed query or server unavailable. Please try again.")
  }
}

retrieveSecondaryStructureResults <-function(requestID, methodCode) {
  secondaryStructureURL <- paste(rtoolsBaseURL, "work/", requestID, "/", 
                                 methodCode, "/structure.txt", sep="")
  secondaryStructureResponseContent <- content(GET(secondaryStructureURL))
  secondaryStructure <- list(sequence=splitString(secondaryStructureResponseContent, split="\n")[2],
                             secondaryStructure=splitString(splitString(secondaryStructureResponseContent, split="\n")[3], split=" ")[1])
  if(methodCode %in% c(1, 2)) {
    basePairProbsURL <- paste(rtoolsBaseURL, "work/", requestID, "/", 
                              methodCode, "/base-pairing-prob.txt", sep="")
    basePairProbsResponseContent <- content(GET(basePairProbsURL))
    maxFieldsNumber <- max(count.fields(textConnection(basePairProbsResponseContent), sep=" "))
    basePairProbsTable <- read.table(textConnection(basePairProbsResponseContent), 
                                     header=FALSE, 
                                     col.names=c("Position", "Nucleotide", paste0("Pairing", seq(1, maxFieldsNumber-2))), 
                                     fill = TRUE, sep=" ")
    return(list(sequence=secondaryStructure[["sequence"]], 
                secondaryStructure=secondaryStructure[["secondaryStructure"]],
                basePairProbabilities=basePairProbsTable[, -ncol(basePairProbsTable)]))
  }
  else {
    return(secondaryStructure)
  }
}

retrieveAlternativeSecondaryStructureResults <- function(requestID) {
  numberAltStructures <- sum(gregexpr("range of Hamming distance",
                                      xml_text(content(GET(paste(rtoolsBaseURL, "cgi-bin/result.cgi?req_id=", requestID, sep="")))),
                                      fixed=TRUE)[[1]] > 0)
  if (numberAltStructures == 1) {
    message("No alternative structures were found. 
            Returning canonical structure.")
    canonicalStructureURL <- paste(rtoolsBaseURL, "work/", requestID, "/", "8", 
                                   "/rintw.range.", 1, ".ss.txt", sep="")
    canonicalStructureResponseContent <- content(GET(canonicalStructureURL))
    canonicalStructure <- list(sequence=splitString(canonicalStructureResponseContent, split="\n")[1],
                               secondaryStructure=splitString(canonicalStructureResponseContent, split="\n")[2])
    return(canonicalStructure)
  }
  alternativeStructures <- vector(mode="list", length=numberAltStructures)
  for (i in seq_len(numberAltStructures)){
    altStructureURL <- paste(rtoolsBaseURL, "work/", requestID, "/", "8", 
                             "/rintw.range.", i, ".ss.txt", sep="")
    altStructureResponseContent <- content(GET(altStructureURL))
    alternativeStructures[[i]] <- list(sequence=splitString(altStructureResponseContent, split="\n")[1],
                                       secondaryStructure=splitString(altStructureResponseContent, split="\n")[2])
  }
  names(alternativeStructures) <- c("alternativeStructure1-Canonical", 
                                    paste("alternativeStructure", seq(2, numberAltStructures), sep=""))
  return(alternativeStructures)
}

extractBasePairProbability <- function(basePairProbsTableField) {
  if (basePairProbsTableField != "") {
    basePairProb <- as.numeric(splitString(basePairProbsTableField, split=":"))
    names(basePairProb) <- c("targetPosition", "probability")
    return(basePairProb)
  }
}

makeCompositeMatrix <-function(basePairProbsMatrix, pairedBases) {
  compositeMatrix <- basePairProbsMatrix
  compositeMatrix[upper.tri(compositeMatrix)] <- 0
  for (row in seq_len(nrow(pairedBases))) {
    compositeMatrix[pairedBases[row, "Position1"], pairedBases[row, "Position2"]] <- 1
  }
  return(compositeMatrix)
}

flattenDotBracket <- function(extendedDotBracketString) {
  extendedDotBracketVector <- splitString(extendedDotBracketString, "")
  simpleDotBracketVector <- character(length=length(extendedDotBracketVector))
  for (i in seq_len(length(extendedDotBracketVector))) {
    if (extendedDotBracketVector[i] %in% c("(", "[", "<", "{", "A", "B", "C", "D")) {
      simpleDotBracketVector[i] <- "("
    } else if (extendedDotBracketVector[i] %in% c(")", "]", ">", "}", "a", "b", "c", "d")) {
      simpleDotBracketVector[i] <- ")"
    } else if (extendedDotBracketVector[i] == ".") {
      simpleDotBracketVector[i] <- "."
    } else {
      stop("Invalid characters present in extended Dot-Bracket string")
    }
  }
  return(paste(simpleDotBracketVector, collapse=""))
}

splitString <- function(string, split, ...) {
  return(unlist(strsplit(string, split=split, ...)))
}

removeNewLines <- function(string) {
  return(gsub("[\r\n]", "", string))
}