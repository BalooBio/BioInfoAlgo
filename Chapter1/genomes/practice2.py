#Ways of getting DNA sequence from text file

with open('prac.txt') as file:
    seq = file.readlines() # get all lines
    seq = [x.strip() for x in seq] # rmv new line chars
    print(len(''.join(seq)))
    print(''.join(seq))


with open('prac.txt') as file:
    seq = file.read()
    seq = [x.strip() for x in seq]
    print(''.join(seq))
