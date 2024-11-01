"""
Implementation of code challenge 1A on page 8 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def pattern_count(sequence, pattern):
    """ Returns number of instances of a string, pattern, in a larger string, sequence.

    :param sequence: larger string to search through, pattern: smaller string to search for
    :type sequence: str , pattern: str
    :return: Number of times pattern is found in sequence
    :rtype: int

    """
    k = len(pattern)
    count = 0
    for p in range(len(sequence)-k+1):
        if sequence[p:p+k] == pattern:
            count += 1
    return count

"""
Implementation of code challenge 1B on page 8 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def frequent_words(s, k):
    """ Returns set of k-mers with the most frequent occurance in a DNA sequence

    :param s: DNA sequence , k: Size of k-mers
    :type s: str, k: int
    :return: set of most frequent strings of len k in larger string s
    :rtype set

    """
    frequentPatterns = []
    count = []
    for i in range(len(s)-k+1):
        pattern = s[i:i+k]
        count.append(pattern_count(s, pattern))

    maxCount = max(count)

    for i in range(len(count)):
        if count[i] == maxCount:
            frequentPatterns.append(s[i:i+k])

    return(set(frequentPatterns))

"""
Implementation of code challenge 1C on page 12 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def reverse_complement(seq):
    """ Returns the reverse complement of a DNA sequence, seq

    :param seq: DNA sequence
    :type seq: str
    :return: A DNA sequence that is the reverse complement of seq
    :rtype str

    """
    reverse = []
    position = len(seq)-1
    for p in range(position+1):
        reverse.append(seq[position-p])
    for p in range(len(reverse)):
        if reverse[p] == 'A':
            reverse[p] = 'T'
        elif reverse[p] == 'T':
            reverse[p] = 'A'
        elif reverse[p] == 'G':
            reverse[p] = 'C'
        elif reverse[p] == 'C':
            reverse[p] = 'G'

    return ''.join(reverse)

"""
Implementation of code challenge 1D on page 13 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def pattern_matching(seq, pattern):
    """ Returns the starting positions of all instances of a pattern in a DNA sequence, seq

    :param seq: DNA sequence , pattern: Smaller DNA sequence to search for
    :type seq: str, pattern: str
    :return: the starting positions of each instance of pattern in the larger
    sequence, seq
    :rtype list
    """

    positions = []
    k = len(pattern)
    s = len(seq)
    for p in range(s-k+1):
        if seq[p:p+k] == pattern:
            positions.append(p)
    return positions

"""
Implementation of code challenge 1E on page 15 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def clump_finder(genome, L, k, t):
    """ Returns a dictionary, clumps, for which each key is a pattern of length k that forms a
    'clump' in the genome, which is defined as t instances of a pattern existing within a length, L,
    of the genome. Each key has an associated list as the value, containing t intergers, each representing
    the starting position each of the t patterns forming the 'clump'

    :params genome: DNA sequence, L: length of span of geomone that can define a clump, k: k-mer length,
    t: Number of times a k-mer must occur in a spane of the genome to be defined as a clump
    :type genome: string, L: int, k: int, t: int
    :return: k-mers that form clumps in the genome and their startign positions in the genome
    :rtype: dict

    """

    gL= len(genome)
    kmers = {}
    for p in range(gL-k+1):
        if genome[p:p+3] in kmers:
            kmers[genome[p:p+3]].append(p)
        else:
            kmers[genome[p:p+3]] = [p]
    clumps= {}
    for key in kmers:
        listLen = len(kmers[key])
        if listLen >= t:
            for i in range(listLen-t+1):
                if kmers[key][i+t-1] <= L-k+1:
                    clumps[key] = (kmers[key][i:i+t])
    return clumps

"""
Implementation of code challenge 1F on page 27 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def minimum_skew(seq):
    """ Returns the position(s) in a genome, seq, where skew is at its minimum.

    :param seq: string representing a bacterial genome
    :return: list, minskew, containging indicies of genome, seq, in whick skew is minimized
    :rtype list

    """
    skew_dict = {}
    skew = 0
    for i in range(len(seq)):
        if seq[i] == 'G':
            skew += 1
            skew_dict[i] = skew
        elif seq[i] == 'C':
            skew -= 1
            skew_dict[i] = skew
        else:
            skew_dict[i] = skew
    minimum = (min(skew_dict.values()))
    minskew = []
    for k, v in skew_dict.items():
        if v == minimum:
            minskew.append(k)
    return minskew

"""
Implementation of code challenge 1G on page 29 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
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


"""
Implementation of code challenge 1H on page 29 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def approximate_pattern_matching(seq, pattern, d):
    """ Finds all approximate occurances of a k-mer, pattern, in a larger DNA
    sequence, seq, where an approximate occurance is defined as a k-mer whose
    hamming distance from pattern is <= d

    :param seq: DNA sequence , pattern: k-mer DNA sequence , d: Hamming distance
    to use to identify approximate patterns.
    :type seq: str, pattern: str , d: Hamming distance
    :return: Starting position in seq for all approximate occurances of pattern
    :rtype: list

    """
    positions = []
    k = len(pattern)
    s = len(seq)
    for p in range(s-k+1):
        if hammingDistance(pattern, seq[p:p+k]) < d+1:
            positions.append(p)
    return positions


"""
Implementation of code challenge 1I on page 30 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from BioInfoAlgo.Chapter1.supporting_functions.neighbors import neighbors
from _1l_pattern_to_number import pattern_to_number
from _1m_number_to_pattern import number_to_pattern

def frequent_words_with_mismatches(text, k, d):
    """ Returns set of k-mers with the most frequent occurance in a DNA
    sequence, s, allowing for up to d number of mismatches

    :param text: DNA sequence, k: k-mer lenght , d: Hamming Distance
    :type text: string containing only uppercase representatives of DNA
                nucleotides, k: int, d: int
    :return: A set of most frequent k-mers in a DNA sequence, s
    :rtype words: set
    """

    """
    Create list, neighborhoods, containing the neighborhood for each k-mer in
    text, where a neighborhood includes the k-mer and all possible k-mers
    containing up to d number of mismatches. For example, the neighborhood for
    'AA' if d = 1 would be ['AA, 'CA', 'GA', 'TA', 'AC', 'AG', 'AT']
    """
    neighborhoods = []
    for i in range(len(text)-k+1):
        neighborhoods.extend(neighbors(text[i:i+k], d))

    """
    Create list, index, mirroring the list, neighborhoods, but replacing each
    k-mer at a given indext with its representative number using the
    pattern_to_number funciton

    Create list, count, mirroring the list, neighborhoods, but replacing each
    k-mer at a given index with the number 1
    """
    index = []
    count = []
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        index.append(pattern_to_number(pattern))
        count.append(1)

    """
    Creat list, sorted_index, to be a sorted version of the list, index, so that
    all occurances of a given pattern are adjacent.

    Loop through sorted_index, updating count at the next index position to be
    equal to count[i] + 1 if count[i] == count[i+1]
    """
    sorted_index = sorted(index)
    for i in range(len(neighborhoods)-1):
        if sorted_index[i] == sorted_index[i+1]:
            count[i+1] = count[i] + 1

    """
    Determine the largest number in the list, count, and set it equal to new
    variable, max_count. Loop through count, adding the pattern represented by the
    number at the same index in sorted_index to a new list, freq_words, if
    count[i] == max
    """
    freq_words = []
    max_count = max(count)
    print(max_count)
    for i in range(len(neighborhoods)):
        if count[i] == max_count:
            pattern = number_to_pattern(sorted_index[i], k)
            freq_words.append(pattern)
    return freq_words


"""
Implementation of code challenge 1J on page 31 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def final_frequent_words(text, k, d):
    """ Returns set of k-mers with the most frequent occurance in a DNA
    sequence, s, allowing for up to d number of mismatches and considering
    reverse complements

    :param text: DNA sequence, k: k-mer lenght , d: Hamming Distance
    :type text: string containing only uppercase representatives of DNA
                nucleotides, k: int, d: int
    :return: A set of most frequent k-mers in a DNA sequence, s
    :rtype words: set
    """

    """
    Create list, neighborhoods, containing the neighborhood for each k-mer, and
    its reverse complement, in text, where a neighborhood includes the k-mer
    and all possible k-mers containing up to d number of mismatches. For
    example, the neighborhood for 'AA' if d == 1 would be ['AA, 'CA', 'GA', 'TA',
    'AC', 'AG', 'AT']
    """
    neighborhoods = []
    for i in range(len(text)-k+1):
        neighborhoods.extend(neighbors(text[i:i+k], d))
        neighborhoods.extend(neighbors(reverse_complement(text[i:i+k]), d))

    """
    Create list, index, mirroring the list, neighborhoods, but replacing each
    k-mer at a given indext with its representative number using the
    pattern_to_number funciton

    Create list, count, mirroring the list, neighborhoods, but replacing each
    k-mer at a given index with the number 1
    """
    index = []
    count = []
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        index.append(pattern_to_number(pattern))
        count.append(1)

    """
    Creat list, sorted_index, to be a sorted version of the list, index, so that
    all occurances of a given pattern are adjacent.

    Loop through sorted_index, updating count at the next index position to be
    equal to count[i] + 1 if count[i] == count[i+1]
    """
    sorted_index = sorted(index)
    for i in range(len(neighborhoods)-1):
        if sorted_index[i] == sorted_index[i+1]:
            count[i+1] = count[i] + 1

    """
    Determine the largest number in the list, count, and set it equal to new
    variable, max_count. Loop through count, adding the pattern represented by the
    number at the same index in sorted_index to a new list, freq_words, if
    count[i] == max
    """
    freq_words = []
    max_count = max(count)
    for i in range(len(neighborhoods)):
        if count[i] == max_count:
            pattern = number_to_pattern(sorted_index[i], k)
            freq_words.append(pattern)
    return freq_words


"""
Implementation of code challenge 1k on page 41 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""

"""
Implementation of code challenge 1L on page 43 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def pattern_to_number(pattern):
    """ Returns the index of a DNA sequence, pattern, in the lexicographically
    ordered list containing all possible k-mers of length, len(pattern). This
    is done using a recursive algorithm that takes advantacge of the observation
    that the number of k mers occuring before pattern in the list is equal to
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
"""
Implementation of code challenge 1M on page 44 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
import math

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

"""
Implementation of code challenge 1N on page 51 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from _1g_hamming_distance import hamming_distance
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
    """ Modification of neighbors(pattern, d) that generates all k-mers with a Hamming distance
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

# Fupporting functions
def approximate_pattern_count(text, pattern, d):
    """ Returns...

    :param text: , pattern: , d:
    :type s: str, k: int
    :rtype
    :return:

    """

    count = 0
    for i in range(len(text)-len(pattern)+1):
        newPattern = text[i:i+len(pattern)]
        if hamming_distance(pattern, newPattern) < d+1:
            count += 1
    return count
