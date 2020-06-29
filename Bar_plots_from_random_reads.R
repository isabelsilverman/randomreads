table = read.table(file = 'freq.tsv', sep = '\t', header = TRUE)

p_reps = list()
for(i in 1:length(table$p_rep1.Frequency)) {p_reps[i] = (table$p_rep1.Frequency[i] + table$p_rep2.Frequency[i])/2}
p_reps = as.numeric(p_reps)

par( mfrow=c(1, 3))
barplot(table$plasmid.Frequency, names.arg = table$Codon, main = "plasmid", 
        xlab = "Codon", ylab = "Frequency", ylim = c(0,0.02))
barplot(table$p0.Frequency, names.arg = table$Codon, main = "p0",
        xlab = "Codon", ylab = "Frequency", ylim = c(0,0.02))
barplot(p_reps, names.arg = table$Codon, main = "p reps 1 and 2",
        xlab = "Codon", ylab = "Frequency", ylim = c(0,0.02))

