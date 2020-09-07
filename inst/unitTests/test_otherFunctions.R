library(RUnit)

## Test flattenDotBracket

checkEquals(flattenDotBracket("...((((..[[.))))]]"), "...((((..((.))))))")

## Test readCT

testCTFile <- system.file("extdata", "exampleCT.ct", package="ncRNAtools")
tmRNASequence <- "GGGGCUGAUUCUGGAUUCGACGGGAUUUGCGAAACCCAAGGUGCAUGCCGAGGGGCGGUUGGCCUCGUAAAAAGCCGCAAAAAAUAGUCGCAAACGACGAAAACUACGCUUUAGCAGCUUAAUAACCUGCUUAGAGCCCUCUCUCCCUAGCCUCCGCUCUUAGGACGGGGAUCAAGAGAGGUCAAACCCAAAAGAGAUCGCGUGGAAGCCCUGCCUGGGGUUGAAGCGUUAAAACUUAAUCAGGCUAGUUUGUUAGUGGCGUGUCCGUCCGCAGCUGGCAAGCGAAUGUAAAGACUGACUAAGCAUGUAGUACCGAGGAUGUAGGAAUUUCGGACGCGGGUUCAACUCCCGCCAGCUCCACCA"

checkTrue(nrow(readCT(testCTFile, tmRNASequence)$pairsTable) == 15)

## Test readDotBracket

testDotBracketFile <- system.file("extdata", "exampleDotBracket.dot", package="ncRNAtools")

checkTrue(readDotBracket(testDotBracketFile)$freeEnergy == -41.2)
