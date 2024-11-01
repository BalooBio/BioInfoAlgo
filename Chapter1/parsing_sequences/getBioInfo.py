def main():
    print(get_DNAseq('Sequence'))

def int(prompt):
    while True:
        try:
            return int(input(prompt))
        except ValueError():
            print(f'{prompt} must be an integer')

def get_DNAseq(prompt):
    while True:
        try:
            seq = input(prompt).strip().upper()
            nucleotide_check(seq)
            return seq
        except ValueError:
            print("Sequence must be a string consisting only of A's, G's, C's, and T's")


def nucleotide_check(seq):
    nucleotides = ['A', 'G', 'C', 'T']
    for i in seq:
        if i not in nucleotides:
            return -1

    pass

if __name__ == '__main__':
    main()
