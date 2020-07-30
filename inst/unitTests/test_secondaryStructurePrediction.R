library(RUnit)

testPrediction <- predictSecondaryStructure("AAAAAAAAAAAAAGGGGGGGGGGGUUUUUUUUUUUUU", "centroidFold")
testAlternativePrediction <- predictAlternativeSecondaryStructures("AAAAAAAAAAAAAGGGGGGGGGGGUUUUUUUUUUUUUCCCCCCCCCCC")
testPairedBases <- findPairedBases(testPrediction$secondaryStructure, testPrediction$sequence)

## Test predictSecondaryStructure

checkEquals(testPrediction$sequence, "AAAAAAAAAAAAAGGGGGGGGGGGUUUUUUUUUUUUU")

## Test predictAlternativeSecondaryStructures

checkEquals(testAlternativePrediction$alternativeStructure2$sequence, "AAAAAAAAAAAAAGGGGGGGGGGGUUUUUUUUUUUUUCCCCCCCCCCC")

## Test generatePairsProbabilityMatrix

checkTrue(is.matrix(generatePairsProbabilityMatrix(testPrediction$basePairProbabilities)))

## Test findPairedBases

checkTrue(is.data.frame(testPairedBases))

## Test pairsToSecondaryStructure

checkEquals(pairsToSecondaryStructure(testPairedBases, testPrediction$sequence), testPrediction$secondaryStructure)

## Test flattenDotBracket

checkEquals(flattenDotBracket("...((((..[[.))))]]"), "...((((..((.))))))")
