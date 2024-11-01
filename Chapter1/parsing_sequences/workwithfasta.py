from Bio import SeqIO

with open('smallestgenome.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(record.id)
