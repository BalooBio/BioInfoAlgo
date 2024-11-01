filename = input('Filename: ')

with open(filename) as f:
    all_lines = f.readlines() # get all lines
    all_lines = [x.strip() for x in all_lines] # rmv new line chars
    print(len(''.join(all_lines)))

