"""
Implementation of pseudocode 1e on page 46 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from _1k_computing_frequencies import computing_frequencies
from _1m_number_to_pattern import number_to_pattern

def main():
    genome = input('Genome: ').strip().upper()
    L = int(input('L: '))
    k = int(input('k: '))
    t = int(input('t: '))
    print(clump_finder(genome, L, k, t))

def clump_finder(genome, L, k, t):
    """ Returns a set, frequent_patterns, containing all patterns of length k that
    form a "clump" in the the genome, where "clump" is defined as t number of
    repeats in an L length region of the genome.

    :params genome: A DNA sequence with no spaces and only uppercase characters,
            L: Length of clump window,
            k: k-mer length,
            t: Number of repeats of a pattern required to form a clump
    :type: genome: str, L: int, k: int, t: int
    :return: frequent_patterns
    :rtype: set

    """
    # Create list, clump, of length 4^k-1 and set all values to 0
    clump = []
    for i in range(4**k):
        clump.append(0)

    """
    Running L-lenth window down genome one base at a time, compute frequency_array
    for each L-len window, and for each frequnecy_array produced, if an
    index position in frequency_array holds a value greater-than or equal to t,
    update that same index position in clump from 0 to 1.
    """
    gL= len(genome)
    for i in range(gL - L + 1):
        text = genome[i:i+L]
        frequency_array = computing_frequencies(text, k)
        for j in range(4**k):
            if frequency_array[j] >= t:
                clump[j] = 1

    """
    For each index in list, clump, if the integer held in that position == 1, compute
    number_to_pattern for that index (number_to_pattern(i, k)), and add the
    retuned pattern to the set, frequent patterns
    """
    frequent_patterns = []
    for i in range(4**k):
        if clump[i] == 1:
            frequent_patterns.append(number_to_pattern(i, k))
    return set(frequent_patterns)

if __name__ == '__main__':
    main()
