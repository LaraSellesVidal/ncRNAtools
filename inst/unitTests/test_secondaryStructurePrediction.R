library(RUnit)

testPrediction <- predictSecondaryStructure("UGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA", "IPknot")
testAlternativePrediction <- predictAlternativeSecondaryStructures("A")
testBasePairProbabilitiesTable <- read.csv(system.file("extdata", "exampleBasePairProbabilitiesTable.csv", package="ncRNAtools"))
testPairedBases <- findPairedBases(testPrediction$secondaryStructure, testPrediction$sequence)

## Test predictSecondaryStructure

checkEquals(testPrediction$sequence, "UGCGAGAGGCACAGGGUUCGAUUCCCUGCAUCUCCA")

## Test predictAlternativeSecondaryStructures

checkEquals(testAlternativePrediction$sequence, "A")

## Test generatePairsProbabilityMatrix

checkTrue(is.matrix(generatePairsProbabilityMatrix(testBasePairProbabilitiesTable)))

## Test findPairedBases

checkTrue(is.data.frame(testPairedBases))

## Test pairsToSecondaryStructure

checkEquals(pairsToSecondaryStructure(testPairedBases, testPrediction$sequence), testPrediction$secondaryStructure)
