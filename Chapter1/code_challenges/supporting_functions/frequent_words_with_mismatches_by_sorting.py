"""
Implementation of pseudocode page 52 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from BioInfoAlgo.Chapter1.supporting_functions.neighbors import neighbors
from _1l_pattern_to_number import pattern_to_number
from _1m_number_to_pattern import number_to_pattern

def main():
    seq = input('Sequence: ').strip().upper()
    k = int(input('K-mer Length: '))
    d = int(input('Hamming Distance: '))
    print(frequent_words(seq, k, d))

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


if __name__ == "__main__":
    main()
