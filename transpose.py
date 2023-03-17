#!/bin/bash python3


import sys
import argparse
import numpy 
import re
#import shlex # splitting while respecting quotes
#import re # for splitting by 2 delimiters

####### Set up #############
parser = argparse.ArgumentParser(description='For transposing data table or matrix')
parser.add_argument('-i', '--input', type=str
                    ,help='File or directory to transpose')
parser.add_argument('-d', '--delimiter', type=str
                    ,help='Input delimiter')
parser.add_argument('-e', '--exportdelimiter', type=str
                    ,help='Export delimiter [ Default: input delimiter ]')
parser.add_argument('-l', '--large', type=bool, default=False, nargs='?'
                    ,help='Is this a very large array? [Default: False]')
parser.add_argument('-o', '--output', type=str, default="transposed.txt"
                    ,help='Transposed file')

args=parser.parse_args()

dat = args.input
d = args.delimiter
e = args.exportdelimiter
output = args.output
large = args.large

if d is None:
	print("\n\nError: Please specify input delimiter\n\n")
	sys.exit()
	
if e is None:
	e = d

####### Transpose #############

f = open(dat, 'r')
contents = [re.split(d,i.strip()) for i in f.readlines()]
f.close()


# Check if array is as you'd expect
if len(contents[0])==1:
	print('WARNING: First row has only a single value. Check that delimiter is correct')

if not large:
	contents_trans = numpy.transpose(contents)
	##### Print ##############
	toPrint = ''

	for i in contents_trans:
		tempLine = ''
		for v in i:
			tempLine += str(v) + "\t"
		tempLine = tempLine.strip()
		tempLine += "\n"
		toPrint += tempLine
		
	toPrint = toPrint.strip()

	writefile = open(output, "w")
	writefile.write(toPrint)
	writefile.close()
else:
	totLen = len(contents[0])
	writefile = open(output, 'w')
	##### Print ##############
	for i in range(totLen):
		tempLine = ''
		for v in [x[i] for x in contents]:
			tempLine += str(v) + "\t"
		tempLine = tempLine.strip()
		if i != (totLen-1):
			tempLine += "\n"
		writefile.write(tempLine)
	writefile.close()
		


