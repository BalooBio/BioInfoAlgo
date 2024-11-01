"""
Implementation of code challenge 1B on page 8 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
from _1a_pattern_count import pattern_count
from BioInfoAlgo.Chapter1.parser import parser

def main():
    seq = parser(input('filename').strip())
    newseq = seq[:3000]
    kmer = int(input('kmer: '))
    print(frequent_words(newseq, kmer))

def frequent_words(s, k):
    """ Returns set of strings of len k with the most frequent occurance in a larger string s

    :param s: , k:
    :type s: str, k: int
    :rtype set
    :return: set of most frequent strings of len k in larger string s
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

if __name__ == '__main__':
    main()
