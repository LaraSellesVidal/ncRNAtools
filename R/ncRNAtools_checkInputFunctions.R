checkMultipleQuery <- function(inputString) {
  if (length(inputString) > 1) {
    stop("Multiple queries at a time not supported")
  }
}

checkInferenceEngine <- function(inferenceEngine, method, sequence) {
  if (method == "centroidFold" & !(inferenceEngine %in% c("CONTRAfold", "BL", "Turner"))) {
    stop(strwrap("Please select a secondary structure inference engine 
                 compatible with CentroidFold. See documentation for details", 
                 initial="", prefix="\n"))
  }
  else if (method == "centroidHomFold" & !(inferenceEngine %in% c("CONTRAfold", "BL", "Turner"))) {
    stop(strwrap("Please select a secondary structure inference engine 
                 compatible with CentroidHomfold. See documentation for
                 details", initial="", prefix="\n"))
  }
  else if (method == "IPknot" & !(inferenceEngine %in% c("CONTRAfold", "BL", "Turner", "NUPACK"))) {
    stop(strwrap("Please select a secondary structure inference engine 
                 compatible with IPknot. See documentation for details", 
                 initial="", prefix="\n"))
  }
  else if (method == "IPknot" & inferenceEngine == "NUPACK" & nchar(sequence)>100) {
    stop("NUPACK model is only available for sequences of 100 or less nucleotides.")
  }
}

checkAlignmentEngine <- function(alignmentEngine, method) {
  if (method != "centroidHomFold") {
    warning(strwrap("Alignment engine ignored, since CentroidHomfold was not
                    used for secondary structure prediction.", 
                    initial="", prefix="\n"))
  }
  else if (!alignmentEngine %in% c("CONTRAlign", "ProbCons")) {
    stop(strwrap("Please select a valid engine for pairwise alignment during 
                 CentroidHomfold prediction. See documentation for details.", 
                 initial="", prefix="\n"))
  }
}

checkGammaWeight <- function(gammaWeight) {
  if (!is.numeric(gammaWeight) || gammaWeight <= 0) {
    stop("Gamma weight must be a positive number.")
  }
}

checkEValueRfamSearch <- function(eValueRfamSearch, method) {
  if (method != "centroidHomFold") {
    warning(strwrap("E-value for searching the Rfam database ignored, since 
                    CentroidHomfold was not used for secondary structure prediction.",
                    initial="", prefix="\n"))
  }
  else if (!is.numeric(eValueRfamSearch) || eValueRfamSearch <= 0) {
    stop("E-value must be a positive number.")
  }
}

checkNumHomSeqsRfamSearch <- function(numHomSeqsRfamSearch, method) {
  if (method != "centroidHomFold") {
    warning(strwrap("Number of homolog sequences for searching the Rfam database
                    ignored, since CentroidHomfold was not used for secondary 
                    structure prediction.", initial="", prefix="\n"))
  }
  else if (!is.numeric(numHomSeqsRfamSearch) || numHomSeqsRfamSearch <= 0) {
    stop("The number of homolog sequences must be a positive number.")
  }
}

checkRNAString <- function(inputString) {
  if(!grepl('^[AUGCaugc]+$', gsub("[\r\n]", "", inputString))) {
    stop("Please provide a sequence string containing only standard RNA symbols")
  }
}

checkBasicDotBracketString <- function(inputString) {
  if(!grepl('^[\\(\\)\\.]+$', gsub("[\r\n]", "", inputString))) {
    stop("Please provide a secondary structure string in the basic Dot-Bracket 
         notation.")
  }
}
