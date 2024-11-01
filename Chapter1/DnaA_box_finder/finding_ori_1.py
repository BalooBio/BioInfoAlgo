from supporting_functions import minimum_skew, dnabox_search, final_frequent_words
from Bio import SeqIO

def main():
    file_name = input('File Name: ').strip()

    # Get genome sequence from FASTA file using seqRecord object
    record = SeqIO.read(file_name, 'fasta')
    seq = record.seq

    # Find location(s) with minimum skew
    minima = minimum_skew(seq)
    print(minima)
    #Search for DnaA boxes in 50bp region centered on each minima
    for i in minima:
        frame = seq[i-250:i+249]
        print(final_frequent_words(frame, 9, 1))


if __name__ == '__main__':
    main()
