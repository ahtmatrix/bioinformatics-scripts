#randomly pick score out of human PWM

# 
# randomly pick a S_NT_position and sum up each 1 million times
# 
# 
# 
# graph score vs freq histogram
# 
# should be exponential / poisson distrubtion
# 





human.filename = "rna30upstream_and_CDS.fasta"

#make a pwm of size 13
human = readDNAStringSet(human.filename)
human.kozak = DNAStringSet(subseq(human, start = 22, end = 34))

for (x in 1:2) {
  human.pwm = NULL
  raw.data = NULL
  unlisted.data = NULL
  if (x == 1) {
    PWM.type = "log2probratio"
  } else if (x == 2) {
    PWM.type = "prob"
  }
  
  human.pwm    = PWM(human.kozak, type = PWM.type)
  
  num.random.samples = 50000
  
  raw.data <- list()
  #loop unrolling
  for (i in 1:num.random.samples) {
    list.scores <- c()
    list.scores <- c(list.scores, sample(human.pwm[1:4], 1))
    list.scores <- c(list.scores, sample(human.pwm[5:8], 1))
    list.scores <- c(list.scores, sample(human.pwm[9:12], 1))
    list.scores <- c(list.scores, sample(human.pwm[13:16], 1))
    list.scores <- c(list.scores, sample(human.pwm[17:20], 1))
    list.scores <- c(list.scores, sample(human.pwm[21:24], 1))
    list.scores <- c(list.scores, sample(human.pwm[25:28], 1))
    list.scores <- c(list.scores, sample(human.pwm[29:32], 1))
    list.scores <- c(list.scores, sample(human.pwm[33:36], 1))
    list.scores <- c(list.scores, sample(human.pwm[37:40], 1))
    list.scores <- c(list.scores, sample(human.pwm[41:44], 1))
    list.scores <- c(list.scores, sample(human.pwm[45:48], 1))
    list.scores <- c(list.scores, sample(human.pwm[49:52], 1))
    
    raw.data <- c(raw.data, sum(list.scores))
  }
  
  unlisted.data <- unlist(raw.data)
  
  hist(
    unlisted.data,
    breaks = num.random.samples,
    xlab = "scores",
    main = paste(num.random.samples, "random samples| PWM type = ", PWM.type)
  )
}