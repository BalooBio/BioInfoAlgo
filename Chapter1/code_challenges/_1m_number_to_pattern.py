11111111112222222222333333333344444444445555555555666666666677777777778888888888
"""
Implementation of pseudocode 1M on page 44 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
import math

def main():
    index = int(input('Index: '))
    k = int(input('kmer length: '))
    print(number_to_pattern(index, k))

def number_to_pattern(i, k):
    """ Returns the DNA sequence, pattern, stored at index, i, in the
    lexicographically ordered list, sequences, containing all possible k-mers of
    length, k. This is done using a recursive algorithm based on the following
    observation: When repeatedly dividing i by 4, the remainder of each division
    represents the index of the ordered list ('A', 'C', 'G', 'T') contaning the
    nucleotide stored at pattern[k-d], where d is the number of times i has been
    divided by 4. A visual example from the book is adapted below:
    i    DividedBy   Quotient   Remainder   Nucleotide
    619      4          154         3          'T'
    154      4           38         2          'G'
    38       4            9         2          'G'
    9        4            2         1          'C'
    2        4            0         0          'A'
    The pattern stored at patterns[619], containing k-mers of length k is
    'ACGGT'

    :param i: index of the list, k: length of k-mers contained in the list
    :type i: int, k: int
    :return: The DNA sequence, pattern, stored at patterns[i]
    :rtype: str

    """
    if k == 1:
        return number_to_nucleotide(i)
    quo = math.floor(i/4)
    rem = i % 4
    nucleotide = number_to_nucleotide(rem)
    prefix_seq = number_to_pattern(quo, k-1)
    return str(prefix_seq) + nucleotide

def number_to_nucleotide(r):
    """ Returns the index of the list nucleotides('A', 'C', 'G', 'T') containing the
    nucleotide, r

    :param r: A DNA nucleotide
    :type r: str
    :return: The index of list, nucleotides, containing r
    :rtype: int

    """
    if r == 0:
        return 'A'
    if r == 1:
        return 'C'
    if r == 2:
        return 'G'
    if r == 3:
        return 'T'


if __name__ == '__main__':
    main()
