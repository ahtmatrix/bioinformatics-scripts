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
human.kozak = DNAStringSet(subseq(human, start = 22, end= 34))
human.pwm    = PWM(human.kozak, type = 'log2probratio')



list.scores <- c()

x = seq(from = 1, to = 52, by = 4)
y = seq(from = 4, to = 52, by = 4)

x
y

list.scores <- c(list.scores, sample(human.pwm[x:y],1))

list.scores <- c(list.scores, sample(human.pwm[1:4],1))
list.scores <- c(list.scores, sample(human.pwm[5:8],1))
list.scores <- c(list.scores, sample(human.pwm[9:12],1))
list.scores <- c(list.scores, sample(human.pwm[13:16],1))
list.scores <- c(list.scores, sample(human.pwm[17:20],1))
list.scores <- c(list.scores, sample(human.pwm[21:24],1))
list.scores <- c(list.scores, sample(human.pwm[25:28],1))
list.scores <- c(list.scores, sample(human.pwm[29:32],1))
list.scores <- c(list.scores, sample(human.pwm[33:36],1))
list.scores <- c(list.scores, sample(human.pwm[37:40],1))
list.scores <- c(list.scores, sample(human.pwm[41:44],1))
list.scores <- c(list.scores, sample(human.pwm[45:48],1))
list.scores <- c(list.scores, sample(human.pwm[49:52],1))


raw.data <- data.frame()
for (i in 1:1000000000) {
  raw.data <- rbind(raw.data, sum(list.scores))
}



list.scores
  
length(human.pwm)

dank <- human.pwm[1:4]
dank <- human.pwm[5:8]
dank


sample(dank, 1)
