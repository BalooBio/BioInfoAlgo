from _1g_hamming_distance import hamming_distance

def main():
    seq = input('Sequence:').upper()
    kmer = input('Pattern: ').upper()
    d = int(input('Hamming Distance: ').upper())
    print(approximate_pattern_count(seq, kmer, d))


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


if __name__ == '__main__':
    main()
