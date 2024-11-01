import math
from BioInfoAlgo.Chapter1.parser import parser
from _1j_frequent_words_with_mismatches_and_reverse_complements import frequent_words_platinum
import matplotlib.pyplot as plt
    plt.plot(skewlist)
    plt.show()

def main():
    seq = parser(input('filename').strip())
    #seq = input('Genome: ').upper()
    frame = int(input('Size of skew frame: '))
    framed_seq = frames(seq, frame)
    minskew = minimum_skew(framed_seq)
    for key in minskew:
        print(frequent_words_platinum(framed_seq[key], 9, 1))


def frames(seq, frame):
    """ Returns a list, frames, containing the genome, seq, split into segments of length, frame.

    :param seq: String representing a bacterial genome , frame: integer specifying the size of the interval in which seq should be split into
    :type seq: str, frame: int
    :rtype: list
    :return: list containing bacterial genome, seq, split into sements of length, frame.

    """


def minimum_skew(seq, frame):
    """ Returns the indicies of a list, frames, which hold the intervals of the genome, seq, for which skew is minimized

    :param seq: string representing a bacterial genome , frame: integer specifying the size of the interval in which seq should be split into
    :type s: str, k: int
    :rtype set
    :return: set of most frequent strings of len k in larger string s

    """

    # Generate list, frames, containing the genome, seq, in segments of length, frame
    # For example
    #   if seq = AAATTTCCCGGG and frame = 3
    #   then frames = ('AAA', 'TTT', 'CCC', 'GGG')
    s = len(seq)
    frameNumber = math.ceil(s/frame)
    frames = []
    for i in range(frameNumber):
        newList = seq[i*frame:i*frame+frame]
        frames.append(newList)
    # Generate dictionary, frameSkew, pairing the indicies of list, frames, with the skew
    # value of the genome segment at that index
    frameSkew = {}
    for i in range(len(frames)):
        score = 0
        for p in frames[i]:
            if p == 'G':
                score += 1
            elif p == 'C':
                score -= 1
            else:
                pass
        frameSkew[frames[i]] = score
    minskew = {}
    minimum = (min(frameSkew.values()))
    for k, v in frameSkew.items():
        if v == minimum:
            minskew[frames.index(k)] = v

    return minskew


if __name__ == '__main__':
    main()
