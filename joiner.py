import os, sys

input = sys.argv[1]

with open(input, 'r') as infile:
	lines = infile.readlines()

#print(lines)

lines = list(map(lambda x: '_'.join(x.split()), lines[:-1]))

#print(lines)

with open(input, 'w') as outfile:
	for element in lines:
		outfile.write(element + '\n')

print('bingo')
