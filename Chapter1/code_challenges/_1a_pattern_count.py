"""
Implementation of code challenge 1A on page 8 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""

def main():
    sequence = input('Sequence: ').upper()
    pattern = input('Pattern: ').upper()
    print(pattern_count(sequence, pattern))



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

if __name__ == '__main__':
    main()

