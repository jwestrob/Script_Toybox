import os, time, sys

for filename in os.listdir('.'):
    if filename.split('.')[-1] == "fa" or filename.split('.')[-1] == "fasta":
        print(filename)
        os.system('grep ">" ' + filename + ' | wc -l')
print("boogie")
