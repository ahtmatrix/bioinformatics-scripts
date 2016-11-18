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


query = readDNAStringSet('viral_30upstream_.CDS.fasta')

for(i in 1:3){
#1 2 or 3 frame changer
frame = 2
#----------------------
#calculate diNT freq of after start to 99
frame.start.pos = 33

#allows easy change of frame
subseq.start = frame.start.pos+frame

#make a subsequence starting after start codon and ending at 99
after.start.diNT.freq = subseq(query, start = subseq.start, end = 99)

#calculates diNT freq of every +<frame> position until end = 99
diNT.freq <- dinucleotideFrequency(after.start.diNT.freq, step = 3)

#sum all frequencies
totals.after.start.diNT.freq <- colSums(diNT.freq)

scale(totals.after.start.diNT.freq)



#-----------------------
#calculate diNT freq of full length sequence WITHOUT START AND WITHOUT STOP

#filter out the start and stop codon from full length sequences

no.start.stop.full <- xscat(subseq(query, start = frame, end = 30), subseq(query, start = 34, end = query@ranges@width-3))

no.start.stop.full.diNT.freq <- dinucleotideFrequency(no.start.stop.full, step = 3)

totals.no.start.stop.full.diNT.freq <-colSums(no.start.stop.full.diNT.freq)

#need to normalize denominator

scale(totals.no.start.stop.full.diNT.freq)




diNT.freq.means <- scale(totals.after.start.diNT.freq)/scale(totals.no.start.stop.full.diNT.freq)



plot(log10(diNT.freq.means), main = paste(c("frame = +"), frame))

