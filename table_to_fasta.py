
import sys
import argparse
import numpy 
import re
#import shlex # splitting while respecting quotes
#import re # for splitting by 2 delimiters

####### Set up #############
parser = argparse.ArgumentParser(description='For converting a table with two columns (ESVID, sequence) into a fasta file')
parser.add_argument('-i', '--input', type=str
                    ,help='Filepath to table')
parser.add_argument('-d', '--delimiter', type=str
                    ,help='Input delimiter')

parser.add_argument('-o', '--output', type=str, default="seq.fasta"
                    ,help='Fasta file')

args=parser.parse_args()

inputFP = args.input
outputFP = args.output
delim = args.delimiter

####### Load in file #######

f = open(inputFP, 'r')
newfasta = open(outputFP, 'w')
linenum = 1
for l in f.readlines():
    if linenum>1:
        temp = l.split(delim)
        newfasta.write('>' + str(temp[0]) + '\n' + str(temp[1]))
    linenum+=1
    # newfasta += '>' + str(temp[0]) + '\n' + str(temp[1])
    # newfasta = newfasta.strip()
f.close()
newfasta.close()





