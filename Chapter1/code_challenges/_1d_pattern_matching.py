"""
Implementation of code challenge 1D on page 13 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""
def main():
    seq = input('Sequence: ')
    pattern = input('Pattern: ')
    print(pattern_matching(seq, pattern))

def pattern_matching(seq, pattern):
    """ Returns ...

    :param seq: , pattern:
    :type seq: str, pattern:
    :rtype list
    :return: ....
    """

    positions = []
    k = len(pattern)
    s = len(seq)
    for p in range(s-k+1):
        if seq[p:p+k] == pattern:
            positions.append(p)
    return positions

if __name__ == '__main__':
    main()
