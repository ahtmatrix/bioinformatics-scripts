library(Biostrings)

#make a pwm of size 13
human = readDNAStringSet('annotated_extracted_TIS_rna_30upstream_.CDS.TIS.fasta')
human.kozak = DNAStringSet(human)
human.pwm    = PWM(human.kozak, type = 'log2probratio')

query = readDNAStringSet('viral_30upstream_.CDS.fasta')


# B) Find out the dinucleotides profile of +1/2/3 frames
# +1 means dint at positions 1&2 of a codon
# +2 means dint at positions 2&3 of a codon
# +3 means dint spans the 3rd position of a codon, and the 1st nt of the next codon
# Compare these 3 profiles to see if there's skewness in dinucleotide occurrences.




#this generates DINUCLEOTIDE FREQUENCIES FOR EVERYTHING AFTER START UNTIL 99 AFTER START
#THEN TAKES DINUCLEOTIDE FREQUENCIES OF THE WHOLE SEQUENCE AND MAKES MEANS
dinucleotide.freqs = NULL
for (x in 1:length(query)) {
  for (i in range(99)) {
    per.seq.diNT.freq <- dinucleotideFrequency(query)
    
    dinucleotide.freqs = rbind(dinucleotide.freqs, per.seq.diNT.freq)
    
    
  }
}



plus.1.diNT.freq <-
  dinucleotideFrequency(subseq(query, start = 34, end = 99), as.prob = TRUE)

totals.plus.1.diNT.freq <- colSums(plus.1.diNT.freq)


all.seq.diNT.freq <- dinucleotideFrequency(query)

totals.all.seq.diNT.freq <- colSums(all.seq.diNT.freq)

normalized.means.diNT.freq <-
  log10(totals.plus.1.diNT.freq / totals.all.seq.diNT.freq)


df <- data.frame(totals.plus.1.diNT.freq, totals.all.seq.diNT.freq)

write.csv(df, "probabilities.diNT.csv")

barplot(
  normalized.means.diNT.freq,
  xlab = "diNT possibilities",
  ylab = "log10(means)",
  main = "diNT vs log10(means)"
)



#THIS GENERATES DINUCLEOTIDE FREQUENCIES PER CODON AFTER START
#start with 34 which is +1 after start codon

#step = 3 to preseve window frame

#change the start position to get different frames




#----------------------------TESTS IF THE START CODON AND STOP COODN IS IN THE RIGHT PLACE
test = readDNAStringSet('rna30upstream_and_CDS.fasta')
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



#checks if the query length may be too short
count <- 0
short.count <- 0
for (i in 1:length(query)) {
  if (query@ranges@width[i] < 99) {
    short.count = short.count + 1
  } else if (query@ranges@width[i] > 99) {
    count = count + 1
  }
}



all.data.frame = NULL
query = readDNAStringSet('rna30upstream_and_CDS.fasta')

smallest.seq <- min(query@ranges@width)

for (i in 1:3) {
  #1 2 or 3 frame changer
  frame = i
  #----------------------
  #calculate diNT freq of after start to 99
  frame.start.pos = 33
  
  #allows easy change of frame
  subseq.start = frame.start.pos + frame
  
  #make a subsequence starting after start codon and ending at 99
  after.start.seq = subseq(query, start = subseq.start, end = query@ranges@width - 3)
  
  #calculates diNT freq of every +<frame> position until end = 99
  after.start.diNT.freq <-dinucleotideFrequency(after.start.seq, step = 3)
  
  #sum all frequencies
  totals.after.start.diNT.freq <- colSums(after.start.diNT.freq)
  
  sum(totals.after.start.diNT.freq)
  
  
  #num.codons.after.start[1] = floor(after.start.seq@ranges@width / 3)+1
  
  #-----------------------
  #calculate diNT freq of full length sequence WITHOUT START AND WITHOUT STOP
  
  #filter out the start and stop codon from full length sequences
  
  no.start.stop.full.seq <-subseq(query, start = 34, end = query@ranges@width - 3)
  
  no.start.stop.full.seq.diNT.freq <-dinucleotideFrequency(no.start.stop.full.seq)
  
  totals.no.start.stop.full.seq.diNT.freq <-
    colSums(no.start.stop.full.seq.diNT.freq)
  sum(totals.no.start.stop.full.seq.diNT.freq)
  #diNT.freq.means <- totals.after.start.diNT.freq/totals.no.start.stop.full.diNT.freq
  
  #diNT.freq.means <-scale(totals.after.start.diNT.freq)/scale(totals.no.start.stop.full.seq.diNT.freq)
  
  #num.codons.full <- floor(no.start.stop.full.seq@ranges@width / 3)
  
  norm.fact <-
    (
      sum(totals.after.start.diNT.freq) / sum(totals.no.start.stop.full.seq.diNT.freq)
    )
  
  final.means <-
    log10(
      totals.after.start.diNT.freq / (totals.no.start.stop.full.seq.diNT.freq *
                                        norm.fact)
    )
  
  barplot(final.means, main = paste(c("frame = +"), frame))
  
  
}

data <-
  data.frame(
    time = seq(0, 23),
    noob = rnorm(24),
    plus = runif(24),
    extra = rpois(24, lambda = 1)
  )
