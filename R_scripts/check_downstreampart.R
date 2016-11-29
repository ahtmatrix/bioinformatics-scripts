library(Biostrings)


#----------------------------TESTS IF THE START CODON AND STOP COODN IS IN THE RIGHT PLACE
test = readDNAStringSet('viral30upstream_STOP_200downstream.fasta')
stop.codon.count <- 0
start.codon.count <- 0
start.errors <- NULL
stop.errors <- NULL
for (i in 1:length(test)) {
  
  
  atg.codon.subseq <- subseq(test[i], start = 31, end = 33)
  
  stop.codon.subseq <-subseq(test[i],start = test@ranges@width[i]-2,end = test@ranges@width[i])
  
  if (atg.codon.subseq == "ATG") {
    start.codon.count = start.codon.count + 1
  } else {
    start.errors <- rbind(start.errors, toString(atg.codon.subseq))
  }
  
  if (stop.codon.subseq == "TAA" || stop.codon.subseq == "TAG" || stop.codon.subseq == "TGA") {
    stop.codon.count = stop.codon.count + 1
  } else {
    stop.errors <- rbind(stop.errors, toString(stop.codon.subseq))
  }
}
