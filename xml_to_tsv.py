#!bin/bash python3


#from os.path import isfile, join
import sys
import argparse
from os import listdir, mkdir
from os.path import isfile, isdir, exists
import shlex # splitting while respecting quotes
import re # for splitting by 2 delimiters

####### Set up #############
parser = argparse.ArgumentParser(description='Input for converting cml to tsv from EBI/ENA')
parser.add_argument('-i', '--input', type=str
                    ,help='File or directory with XML files (files must end in .xml, directories should not have / at end)')
parser.add_argument('-o', '--output', type=str, default="tsv_outputs"
                    ,help='Directory to put tsv files')

args=parser.parse_args()



if isdir(args.input):
	allFiles_temp = listdir(args.input)
	allFiles = [args.input+'/'+f for f in allFiles_temp if f.endswith('.xml')]
elif isfile(args.input) and args.input.endswith('.xml') :
	allFiles = [args.input]
else:
	print('Not file or directory')
	sys.exit() 


# FOR TESTING
#allFiles_temp = listdir('00_xml_meta')
#allFiles = ['00_xml_meta'+'/'+f for f in allFiles_temp if f.endswith('.xml')]

################################
#currF='00_xml_meta/SAMN07360285.xml'
allInfoDict = {}
allInfoDict_table = {}
for currF in allFiles:

	f = open(currF, 'r')
	contents = [i.strip() for i in f.readlines()]
	f.close()

	## Get sample accession
	sampacc = re.sub('.xml','', re.sub('^.*/','',currF)) ### NEED TO TEST STILL
#	sampacc = re.sub('.xml','', re.sub('00_xml_meta/','',currF))
	

	# Get sample info

	sampleinfo = [keep.split('=') for keep in shlex.split([i for i in contents if i.startswith('<SAMPLE ')][0]) if '<' not in keep]

	# get identifier info
	#sampleIdentifiers = [j[0:2] for j in [list(filter(None, i)) for i in [re.split('>|<', i) for i in [contents[i] for i in list(range(contents.index("<IDENTIFIERS>")+1, contents.index("</IDENTIFIERS>")))]]]]
	#sampleIds = [[re.sub(' namespace=.+', '', i) for i in j] for j in sampleIdentifiers]

	sampleIdentifiers = [[re.sub(' namespace=.+', '', i) for i in k] for k in [j[0:2] for j in [list(filter(None, i)) for i in [re.split('>|<', i) for i in [contents[i] for i in list(range(contents.index("<IDENTIFIERS>")+1, contents.index("</IDENTIFIERS>")))]]]]]

	tags = [re.sub('<[^>]+>', '', i) for i in contents if '<TAG>' in i]
	values = [re.sub('<[^>]+>', '', j) for j in contents if '<VALUE>' in j]
	
	# get DB and IDs
	DB = [re.sub('<[^>]+>', '', i) for i in contents if '<DB>' in i]
	IDs = [re.sub('<[^>]+>', '', j) for j in contents if '<ID>' in j]

	tagvalues = [[tags[i],values[i]] for i in range(len(tags))]
	dbIDs = [[DB[i],IDs[i]] for i in range(len(DB))]
	
	### Combine info into dictionary
	singleDic = {sampacc : {i[0]:i[1] for i in sampleinfo + sampleIdentifiers + tagvalues+dbIDs}}
	allInfoDict.update(singleDic)
	# Make full list of ALL headers collected
	allInfoDict_table.update({k:[] for k in singleDic[sampacc].keys() if k not in allInfoDict_table.keys()})
	

##### Re-iterate to make format that guarantees all items present
allInfoDict_table.update({"SampleAccession":[]})
for samp in allInfoDict.keys():
	for val in allInfoDict_table.keys():
		if val == "SampleAccession":
			allInfoDict_table["SampleAccession"].append(samp)
		elif val in allInfoDict[samp].keys():
			allInfoDict_table[val].append(allInfoDict[samp][val])
		else:
			allInfoDict_table[val].append('')
#### Make into printable format

#print(allInfoDict_table)

######## FIGURE OUT HOW TO WRITE #############
toWrite = ''

# header
for header in allInfoDict_table:
	toWrite = toWrite + header + '\t'
toWrite.strip()
#toWrite = toWrite + '\n'
# body
for r in range(len(allInfoDict_table['SampleAccession'])):
	toWrite = toWrite + '\n'
	for header in allInfoDict_table:
		toWrite = toWrite + allInfoDict_table[header][r] + '\t'
	toWrite.strip()

# Print
out = open(args.output, 'w')
out.write(toWrite)
out.close()

print('DONE COLLAPSING XML TO TSV')




