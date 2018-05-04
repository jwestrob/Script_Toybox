import os, sys, csv
import argparse

parser=argparse.ArgumentParser(description='Do some things to all the files in your directory, NERD')

parser.add_argument('-gzip', metavar='gunzip all your files!', nargs='?', default=None, help='Do you wanna do it or not?')
parser.add_argument('-q', metavar='-q option for fastq_quality_filter', help='GIMME AN INT!')
parser.add_argument('-p', metavar='-p option for fastq_quality_filter', help='GIMME AN INT!')

args=parser.parse_args()

q_opt = int(args.q)
p_opt = int(args.p)
if args.gzip != None:
	gzip = True
else:
	gzip = False


def gunzipper():
	for filename in os.listdir('.'):
		if filename.split('.')[-1] == 'gz':
			os.system('gunzip ' + filename)
	return


def fastqq(q_opt, p_opt):
	for filename in os.listdir('.'):
		if filename.split('.')[-1] == 'fq':
			os.system('fastq_quality_filter -Q 33 -q ' + q_opt + ' -p ' + p_opt + ' -z -i ' + filename + ' -o ' + filename.split('.')[0] + '_TREATED.fastq')

	return

def main():
	if gzip:
		gunzipper()
	
	fastqq(q_opt, p_opt)
	print("boogie")
if __name__ == "__main__":
	main()

