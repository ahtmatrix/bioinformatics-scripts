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

# virus = readDNAStringSet('annotated_extracted_TIS_viral_30upstream_.CDS.TIS.fasta')
# virus.kozak = DNAStringSet(virus)
# virus.cons   = consensusString(virus.kozak)
# virus.pwm    = PWM(virus.kozak, type = 'log2probratio')

#make a pwm of size 13
human = readDNAStringSet('annotated_extracted_TIS_rna_30upstream_.CDS.TIS.fasta')
human.kozak = DNAStringSet(human)
human.pwm    = PWM(human.kozak, type = 'log2probratio')

query = readDNAStringSet('viral_30upstream_.CDS.fasta')

#print(PWMscoreStartingAt(human.pwm, query[[1]], starting.at = 99))


#loop through ever single element
pwm.scores = NULL

print(length(query))

for (i in 1:length(query)) {
  #if (grepl("Scanning", query@ranges@NAMES[i]) == TRUE) {
  score <-
    PWMscoreStartingAt(human.pwm, query[[i]], starting.at = 1:99)
  
  pwm.scores = cbind(pwm.scores, score)
  
  #}
}

pwm.score.means <- rowMeans(pwm.scores)

plot(
  pwm.score.means,
  type = "l",
  xlab = "position",
  ylab = "PWM Score",
  main = ""
)


#pick particular virus
#small group
#lower sample size
#12 sequences
#

#random sequences should expect steady state

#stddev near true ATG should be 0

#starndard devation per row




# ## Use human PWM to score viral seqs
# ## size of kozak sequence
# ## score only the beginning 336 bps (same as internal atg file)
# ksize = 13
# max.len = 336 - ksize + 1
# seqfiles = c('')
# outfiles = c('viral_orf_u36.tsv', 'viral_internal_atg.tsv')
# #seqfiles = c('viral_test.fasta','test_internal.fasta')
# #outfiles = c('viral_test.tsv','test_internal.tsv')
#
#
# for (idx in 1:length(seqfiles)) {
#   fname = seqfiles[idx]
#   print(fname)
#
#   virus = readDNAStringSet(fname)
#   all.scores = c()
#   all.scores.rownames = c()
#
#   for (i in 1:length(virus)) {
#     s.name = names(virus[i])
#     print(s.name)
#     s = virus[[i]]
#
#     if ((length(s) - ksize + 1) >= max.len) {
#       all.scores.rownames = c(all.scores.rownames, names(virus[i]))
#
#       s.scores = c()
#
#       for (st in 1:max.len) {
#         score = PWMscoreStartingAt(hs.pwm, s[st:(st + ksize - 1)])
#         s.scores = c(s.scores, score)
#       }
#
#       all.scores = rbind(all.scores, s.scores)
#     }
#   }
#   pos.names = c(seq(-27,-1, 1), seq(1, 297, 1))
#   colnames(all.scores) = pos.names
#   rownames(all.scores) = all.scores.rownames
#   write.table(all.scores,
#               file = outfiles[idx],
#               col.names = T,
#               sep = '\t')
# }
#
#
# ## Plot graphs
# real = read.table(outfiles[1], head = T, sep = '\t')
# internal = read.table(outfiles[2], head = T, sep = '\t')
# plot(
#   colMeans(real),
#   xlab = 'Position from ATG',
#   ylim = c(0, 1),
#   ylab = 'PWM scores',
#   type = 'l',
#   axes = F
# )
# pos.names = c(seq(-27,-1, 10), seq(1, 297, 10))
# axis(side = 1, at = pos.names)
# axis(side = 2, at = seq(0, 1, 0.1))
# #lines(colMeans(internal),col='blue')