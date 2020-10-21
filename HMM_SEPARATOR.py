import os, sys

try:
    big_hmmfile = sys.argv[1]
except:
    print("Give me a big HMM file to separate out!")
    sys.exit(420)

f = open(big_hmmfile)

names = []
files = []

for line in f.readlines():
    if line.startswith('HMMER3'):
        files.append([])
    if line.startswith('NAME'):
        try:
            names.append(line.rstrip('\n').split('\t')[1])
        except:
            try:
                names.append(line.rstrip('\n').split('  ')[1])
            except:
                print(line.rstrip('\n').split('\t'))
    files[-1].append(line)

for index, file in enumerate(files):
    name = names[index]
    #if 'rhodopsin' not in name.lower():
    #    continue
    filename = name + '.hmm'
    with open(filename, 'w') as handle:
        list(map(lambda x: handle.write(x), file))



#os.system('Rscript ~/scripts/dale.R')
