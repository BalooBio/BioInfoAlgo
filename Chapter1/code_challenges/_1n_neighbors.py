"""
Implementation of code challenge 1N on page 51 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from _1g_hamming_distance import hamming_distance
def main():
    pattern = input('Pattern: ').strip().upper()
    d = int(input('Hamming Distance: '))
    print(neighbors(pattern, d))

def neighbors(pattern, d):
    """ Returns the neighborhood of a k-mer, pattern, based on the specified
    hamming distance, d, where neighborhood is defined as a set of all k-mers,
    including pattern, whose hamming distance from pattern is <=d. This is done
    using a recursive algorithm built on the observation that the neightbors of
    pattern[1:] (referred to as suffix_pattern) will always have a hamming distance
    from pattern <=d. Thus, they can be used to fill out the neighbors of
    pattern by adding pattern[0] to the front of suffix_pattern's whose hamming
    distance from pattern > d and adding 'A', 'C', 'G', and 'T' (making 4 new kmers)
    to the front of suffix_pattern's whose hamming distance from pattern is = d. An
    example from the book is illustrated below:

    For neighbors(CAA, 1), first compute neighbors (AA, 1)

        For neighbors(AA, 1), first compute neighbors(A, 1)

            neighbors(A, 1) = ['A', 'C', 'G', 'T']

            Now, to compute neighbors(AA, 1), add prefix_pattern or ['A', 'C', 'G', 'T'] to front of each
            component of neighbors(A, 1) based on the rules described above.
            - Hamming distance of 'A' from 'AA' < d, so prefix with ['A', 'C', 'G', 'T']
                    ['AA', 'CA', 'GA', 'TA']
            - Hamming distance of 'C', 'G', and 'T' from 'AA' == d, so prefix with prefix_pattern
                    ['AC', 'AG', 'AT']

        neighbors(AA, 1) = ['AA', 'CA', 'GA', 'TA', 'AC', 'AG', 'AT']

        Now, to compute neighbors(CAA, 1), add prefix_pattern or ['A', 'C', 'G', 'T'] to front of each
            component of neighbors(AA, 1) based on the rules described above.
            - Hamming distance of 'AA' from 'CAA' < d, so prefix with ['A', 'C', 'G', 'T']
                    ['AAA', 'CAA', 'GAA', 'TAA']
            - Hamming distance of 'CA', 'GA', 'TA', 'AC', 'AG', and 'AT' fron 'CAA' == d, so prefix with
              prefix_pattern
                    ['CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT']

    Thus neightbors(CAA, 1) = ['AAA', 'CAA', 'GAA', 'TAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT']

    :param pattern: kmer sequence of DNA, d: hamming distance to use for generating
    the neighbors
    :type pattern: str, d: integer
    :return: All neighbor kmers for pattern, including pattern
    :rtype: set

    """
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A', 'G', 'C', 'T']
    nucleotides = ['A', 'G', 'C', 'T']
    neighborhood = []
    prefix = pattern[0]
    suffix_pattern = pattern[1:]
    suffix_neighbors = neighbors(suffix_pattern, d)
    for suffix in suffix_neighbors:
        if hamming_distance(suffix, suffix_pattern) < d:
            for n in nucleotides:
                neighborhood.append(n+suffix)
        else:
            neighborhood.append(prefix+suffix)
    return set(neighborhood)

def specific_neighbors(pattern, d):
    """ Modification of neighbors that generates all k-mers of Hamming distance
    exactly d from pattern.
    """

    if d == 0:
        return pattern
    if len(pattern) == 1:
        if d == 1:
            return nucleotide_not(pattern)
    nucleotides = ['A', 'G', 'C', 'T']
    neighborhood = []
    prefix = pattern[0]
    suffix_pattern = pattern[1:]
    suffix_neighbors = neighbors(suffix_pattern, d)
    for suffix in suffix_neighbors:
        if hamming_distance(suffix, suffix_pattern) == d:
            neighborhood.append(prefix+suffix)
        if hamming_distance(suffix, suffix_pattern) == d-1:
            for n in nucleotide_not(prefix):
                neighborhood.append(n+suffix)
    return set(neighborhood)

def nucleotide_not(n):
    """
    Returns list containing nucleotides that are not n. For example,
    nucleotide_not('A') will return ['C', 'G', 'T']

    """
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(n)
    return nucleotides

if __name__ == '__main__':
    main()
