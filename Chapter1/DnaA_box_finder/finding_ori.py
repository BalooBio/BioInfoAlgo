from supporting_functions import minimum_skew, dnaabox_search, final_frequent_words, ori_test
from Bio import SeqIO

def main():
    file_name = input('File Name: ').strip()

    # Get genome sequence from FASTA file using seqRecord object
    record = SeqIO.read(file_name, 'fasta')
    seq = record.seq

    # Find index location(s) in genome with minimum skew
    minima = minimum_skew(seq)

    #Create list, ori_candidates, to hold all possible index positions for oric, ranging from GC-minima - 500bp to GC-minima + 499 bp"""

    start = minima[0]-500
    ori_candidates = []
    for i in range(1000):
        ori_candidates.append(start+i)

    """ Identify the general index position of ori by passing ori_candidates to the ori_test(seq, candidates, n) funciton, which passes the 500bp
    sequence centered at each ori_candidate to the final frequent words funciton and returns those sequences if any of the most common 9-mers
    match the DnaA box consensus sequence  5' TTATC[CA]A[CA]A 3'."""

    ori_index = ori_test(seq, ori_candidates, 100)
    ori_seq = seq[ori_index[0]-250:ori_index[0]+249]
    dnaabox = (dnaabox_search(ori_seq))

    print(f'\nProposed Origin:\n\nPosition = {ori_index[0]-250}:{ori_index[0]+249}\n\nSequence = {ori_seq}')


    number = 1
    print('\nProposed DnaA Boxes:\n')
    for key, value in dnaabox.items():
        print(f'{number}: Position = {ori_index[0] - 250 + key}, Sequence = {value}\n')
        number += 1


if __name__ == '__main__':
    main()
