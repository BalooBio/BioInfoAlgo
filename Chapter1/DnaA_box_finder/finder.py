from Bio import SeqIO

record = SeqIO.read(input('File Name: '), 'fasta')
seq = record.seq

seq2 = input('pattern')

positions = []
k = len(seq2)
s = len(seq)
for p in range(s-k+1):
    if seq[p:p+k] == seq2:
        positions.append(p)
print(positions)
