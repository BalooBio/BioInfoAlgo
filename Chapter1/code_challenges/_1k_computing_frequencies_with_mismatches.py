"""
Implementation of a variant of the pseudocode 1k on page 41 of Bioinformatics
Algorithms: An Active Learnign Approach, 3rd Edition, by Phillip Compeau and
Paval Pevzner. The variant pseudocode is on page 49.
"""
from _1l_pattern_to_number import pattern_to_number

def main():
    seq = input('Sequence: ').strip().upper()
    k = int(input('k-mer Length: '))
    print(computing_frequencies(seq, k))

def computing_frequencies(text, k):
    """Returns list of size 4**k-1 called frequency_array. Each index of
    frequency_array stores the the number of times the k-mer, produced by
    computing number_to_pattern(index), appears in text.

    :param text: A DNA sequence with no spaces and only uppercase characters, k:
    Length of the patterns you're searching for
    :type text: string, k: int
    :return: List called frquency_array containing number of occurances of all
    possible k-mers of length k, in text
    :rtype: list

    """

    #Create list, frequency_array, of size 4**k-1 and set all index positions to = 0
    frequency_array = []
    for i in range(4**k):
        frequency_array.append(0)

    """
    For all patterns of k-length in text, compute the patterns index position in
    frequency_array using pattern_to_number(text, k), and update the value at that
    position in frequency_array += 1
    """
    for i in range(len(text)-k):
        j = pattern_to_number(text[i:i+k])
        frequency_array[j] += 1

    return(frequency_array)

if __name__ == '__main__':
    main()
