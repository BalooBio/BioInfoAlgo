import matplotlib.pyplot as plt
from BioInfoAlgo.Chapter1.parser import parser


def main():
    #minimum_skew(parser(input('filename').strip()))
    minimum_skew(input('Genome: ').upper())

def minimum_skew(seq):
    """ Returns the position(s)in genome, seq, where skew is at its minimum which hold the intervals of the genome, seq, for which skew is minimized

    :param seq: string representing a bacterial genome
    :rtype list
    :return: list, minskew, containging indicies of genome, seq, in whick skew is minimized

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
    plt.plot(minskew)
    plt.show()
    return minskew


if __name__ == '__main__':
    main()

