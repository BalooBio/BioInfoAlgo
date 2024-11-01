"""
Implementation of code challenge 1C on page 12 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""

def main():
    seq = input('Sequence: ').upper()
    print(reverse_complement(seq))

def reverse_complement(seq):
    """ Returns the reverse complement of a DNA sequence, seq.

    :param seq: A DNA sequence containing only uppercase characters ('A', 'C', 'G', 'T') and
    no empty spaces
    :type seq: str
    :return: The reverse complement of seq
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

if __name__ == '__main__':
    main()

