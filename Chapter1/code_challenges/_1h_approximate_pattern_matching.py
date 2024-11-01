"""
Implementation of code challenge 1H on page 29 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""

from _1g_hamming_distance import hamming_distance

def main():
    seq = input('Sequence: ').upper()
    pattern = input('Pattern: ').upper()
    d = int(input('Allowable Mismatches: '))
    print(approximate_pattern_matching(seq, pattern, d))

def approximate_pattern_matching(seq, pattern, d):
    """ Returns starting positions of pattern in larger string, sequence, allowing for
    d number of mismatches.
    :param seq: DNA sequence , pattern: k-mer sequence, d: Allowable hamming distance
    :type seq: str, pattern: str , d:int
    :return all starting positions of pattern in seq
    :rtype list

    """

    positions = []
    k = len(pattern)
    s = len(seq)
    for p in range(s-k+1):
        if hammingDistance(pattern, seq[p:p+k]) < d+1:
            positions.append(p)
    return positions

if __name__ == '__main__':
    main()
