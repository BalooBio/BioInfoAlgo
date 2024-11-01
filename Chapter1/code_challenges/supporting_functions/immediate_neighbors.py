"""
Implementation of pseudocode page 50 of Bioinformatics Algorithms: An
Active Learnign Approach, 3rd Edition, by Phillip Compeau and Paval Pevzner
"""

def main():
    pattern = input('Pattern: ').strip().upper()
    print(immediate_neighbors(pattern))

def immediate_neighbors(pattern):
    """
    Returns set containing pattern plus all DNA squences with a hamming distance
    of 1 compared to pattern. 
    """
    # Create list containing only pattern
    neighborhood = [pattern]
    # Convert pattern to list
    plist = list(pattern)
    for i in range(len(plist)):
        n_not = nucleotide_not(plist[i])
        for j in n_not:
            new_plist = plist
            new_plist[i] = j
            new_pattern = ''.join(new_plist)
            neighborhood.append(new_pattern)
    return set(neighborhood)

def nucleotide_not(n):
    """
    Returns list containing nucleotides that are not n. For example,
    nucleotide_not('A') will return ['C', 'G', 'T']

    """
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(n)
    return nucleotides

if __name__ == '__main__':
    main()
