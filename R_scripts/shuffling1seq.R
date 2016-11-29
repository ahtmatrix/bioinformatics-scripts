# Our next step is to score translation start sites by translation initiation mechanism.
#
# Attached is a sample code that does similar thing but not exactly.
#Here are the things you should do in your codes:
# 1. Collect all coding sequences with a bit of 5'UTRs, say 30 bps upstream from ATG.
# 2. Use viral genes that use scanning mechanism to build a PWM by R's Biostrings package
#   (google for documentation, there is only one site)
# 3. Score a sequence by scanning from 5' to 3', including 5'utr & coding region (ORF).
#  It should produce a score profile where y-axis is the score, x-axis is position.
#
# In theory, viral genes using scanning mechanism should
# show a peak at the annotated translation start site.
#  Whereas, viral genes using other mechanisms should not show a peak. Fingers crossed!

library(Biostrings)
library(matrixStats)
library(data.table)

# virus = readDNAStringSet('annotated_extracted_TIS_viral_30upstream_.CDS.TIS.fasta')
# virus.kozak = DNAStringSet(virus)
# virus.cons   = consensusString(virus.kozak)
# virus.pwm    = PWM(virus.kozak, type = 'log2probratio')

human.filename = "rna30upstream_and_CDS.fasta"
query.filename = "finalviral.fasta"


#make a pwm of size 13
human = readDNAStringSet(human.filename)
human.kozak = DNAStringSet(subseq(human, start = 22, end= 34))
human.pwm    = PWM(human.kozak, type = 'log2probratio')

query = readDNAStringSet(query.filename)

#print(PWMscoreStartingAt(human.pwm, query[[1]], starting.at = 99))

#smallest.sequence <- min(query@ranges@width)-13

#loop through ever single element
pwm.scores = NULL

size.select = 140

  #atg.codon.subseq <- subseq(query[i], start = 31, end = 33)
  
i  = 1
if(query@ranges@width[i] > size.select){
    #if (grepl("Scanning", query@ranges@NAMES[i]) == TRUE) {
    
    score <- PWMscoreStartingAt(human.pwm, query[[i]], starting.at = 1:130)
    
    #longest length the score can do is 168 b/c that is the shortest sequence
    pwm.scores = cbind(pwm.scores, score)
}


shift <- function(d, k) rbind( tail(d,k), head(d,-k), deparse.level = 0 )

osciliation.generator.df <- data.frame()
osciliation.generator.df <- rbind(osciliation.generator.df, pwm.scores)


for(i in 1:2528){
  
  random.chance = sample(1:129, 1, replace = T)
  shifted <- shift(osciliation.generator.df[1], random.chance)
  osciliation.generator.df <- cbind(osciliation.generator.df, shifted)
  
}

osctest.gen.means <- rowMeans(osciliation.generator.df)

plot(
  osctest.gen.means,
  type = "l",
  xlab = "position",
  ylab = "PWM Score",
  ylim = c(0.1, 0.8),
  main =  paste(human.filename, "-->PWM-->", query.filename)
)

