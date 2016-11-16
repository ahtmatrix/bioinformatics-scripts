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



dinucleotide.freqs = NULL
for (x in 1:length(query)) {
  for (i in range(99)) {
    
      per.seq.diNT.freq <- dinucleotideFrequency(query)

      dinucleotide.freqs = rbind(dinucleotide.freqs, codon.diNT.freq)

    
  }
}

subseq(query, start = 34)

plut.1.diNT.freq <- dinucleotideFrequency( subseq(query, start = 34, end = 99))

totals.plus.1.diNT.freq <- colSums(diNT_freq)


all.seq.diNT.freq <- dinucleotideFrequency(query)

totals.all.seq.diNT.freq <- colSums(all_seq_diNT)

normalized.means.diNT.freq <- log10(totals.plus.1.diNT.freq/totals.all.seq.diNT.freq)

barplot(normalized.means.diNT.freq, xlab = "diNT possibilities",ylab = "normalized means",  main = "diNT vs normalized means")





#start with 34 which is +1 after start codon

#step = 3 to preseve window frame

#change the start position to get different frames


diNT.dataframe = NULL
#i should be 3 every time
for (i in 0:99) {
  if (i %% 3 == 0 ){
    
    subseq.start = 36 + i
    subseq.end = 38 + i
    
    codon.window = subseq(query, start = subseq.start, end = subseq.end)
    
    
    
    diNT.freq <- dinucleotideFrequency(codon.window, step = 3)
    
    codon.diNT.totals <- colSums(diNT.freq)
    
    diNT.dataframe = cbind(diNT.dataframe, codon.diNT.totals)
    
    }
}

barplot(diNT.dataframe)

write.table(diNT.dataframe, "diNT_+3frame.csv")

