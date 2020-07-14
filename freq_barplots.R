table = read.table(file = 'freq.tsv', sep = '\t', header = TRUE)

codons = table$Codon
table = table[-c(1,2,3)]
ncol = length(colnames(table))
table = table[(ncol/2+1): ncol]
ncol = ncol/2

p_reps = vector(mode = "numeric", length = length(codons))
for(i in 1:length(p_reps)) {
  p_reps[i] = 0
  for(j in c(3,ncol)) {
    p_reps[i] = p_reps[i] + table[i,j]}
  p_reps[i] = p_reps[i]/(ncol-2)}

png("freq_barplots.png")
par( mfrow=c(1, 3))
for(i in c(1,2)) {
  barplot(unlist(table[i]), names.arg = codons, main = colnames[i], 
          xlab = "Codon", ylab = "Frequency", ylim = c(0,0.02))
}
barplot(p_reps, names.arg = codons, main = "p_reps.Frequency", 
        xlab = "Codon", ylab = "Frequency", ylim = c(0,0.02))
dev.off()
