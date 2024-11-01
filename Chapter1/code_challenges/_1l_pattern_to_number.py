"""
Implementation of pseudocode 1L on page 43 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def main():
     print(pattern_to_number(input('k-mer: ').strip().upper()))


def pattern_to_number(pattern):
    """ Returns the index of a DNA sequence, pattern, in the lexicographically
    ordered list containing all possible k-mers of length, len(pattern). This
    is done using a recursive algorithm that takes advantacge of the observation
    that the number of k mers occuring before pattern in theh list is equal to
    4 times the number of k-1 mers occuring before pattern[0:k-2] plus the
    number of 1-mers occuring before pattern[k-1], in their respective
    lexicographically orderd lists of k-1 mers and 1 mers. This observation is
    simplified in the equation below.

    pattern_to_number(pattern) = 4 * pattern_to_number(pattern[0:k-2]) +
    pattern_to_number[k-1]

    :param pattern: A DNA sequence with no spaces and only uppercase characters. Example: 'ATGC'
    :type pattern: str
    :return: The index of pattern in the lexicographically ordered list containing all possible
    k-mers of length, len(pattern)
    :rtype: int

    """
    if pattern == 0:
        return 0
    k = len(pattern)
    symbol = pattern[k-1]
    prefix = pattern[0:k-1]
    if k == 1:
        prefix = 0
    return(4 * pattern_to_number(prefix) + ACGT_to_number(symbol))

def ACGT_to_number(n):
    """ Returns the index a nucleatide, n, in the lexicographically
    ordered list containing the four DNA nucletides, A, C, G, and T.

    :param n: A sinlge DNA nucleotide
    :type pattern: str of length 1
    :return: The index of n in the list ('A', 'C', 'G', 'T')
    :rtype: int

    """
    if n == 'A':
        return 0
    if n == 'C':
        return 1
    if n == 'G':
        return 2
    if n == 'T':
        return 3


if __name__ == '__main__':
    main()
