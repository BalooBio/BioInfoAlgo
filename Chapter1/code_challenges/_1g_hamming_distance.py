"""
Implementation of code challenge 1G on page 29 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def main():
    seq1 = input('Sequence 1: ').upper()
    seq2 = input('Sequence 2: ').upper()
    print(hamming_distance(seq1, seq2))

def hamming_distance(x, y):
    """ Returns the hamming distance between two DNA sequences, x and y, where
    hamming distance is defined as the number of mismatched bases

    :param x: DNA sqeunce, y: DNA sequence
    :type x: str, y: str
    :return: Hamming distance between x and y
    :rtype: int

    """

    h = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            h += 1
    return h

if __name__ == '__main__':
    main()
