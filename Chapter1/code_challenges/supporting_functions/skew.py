def main():
    print(skew(input('Sequence: ').upper()))

def skew(seq):
    s = 0
    for p in seq:
        if p == 'G':
            s += 1
        elif p == 'C':
            s -= 1
        else:
            pass

    return s


if __name__ == '__main__':
    main()
