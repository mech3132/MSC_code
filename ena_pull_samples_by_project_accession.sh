#!bin/bash

#### Script to get list of samples from ENAEBI #####

################################################################################################
#### Help ####
################################################################################################

Help()
{
   # Display Help
   echo
   echo
   echo "This script takes a list of ENA/EBI project accessions and extracts all samples found in these projects. First input should be a .txt file where each line is an accession number. Second input should be directory to save this into"
#   echo "-m 
#   echo "options:"
#   echo "g     Print the GPL license notification."
#   echo "h     Print this Help."
#   echo "v     Verbose mode."
#   echo "V     Print software version and exit."
   echo
}


while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

################################################################################################
#### Get study accessions ####
################################################################################################

########## Get all studies in list
# Check if output directory input exists
if [ ${#2} == 0 ]
then
output="SAMPLES"
else
output="$2"
fi

# Make output folder
mkdir -p $output
mkdir -p $output/study_acc

while read s; do
echo "Prcessing study $s"

if [ -e $output/study_acc/$s.tsv ]
then
echo "Already downloaded $s"
else
wget -q "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${s}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,experiment_alias,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" -O $output/study_acc/$s.tsv
fi

done < $1

######### Remove studies that are not in list

for f in $output/study_acc/*
do
studysearch=$(echo $f | sed 's/\.tsv//g' | sed 's/^.*\///g')

if grep -q $studysearch $1
then
echo "keep $studysearch"
else
echo "remove $studysearch"
rm $output/study_acc/$studysearch.tsv
fi

done 



################################################################################################
#### Get metadata xml files ####
################################################################################################

###### Get metadata xml files #######
mkdir -p $output/xml_meta

# sample accession is always column 2
echo 'GETTING ALL SAMPLE ACCESSIONS'
:>$output/all_sample_accessions.txt
for filename in $output/study_acc/*.tsv; do
	echo ${filename}
	awk '{ print $2 }' $filename | grep -v 'sample_accession' >> $output/all_sample_accessions.txt
done

# Look through all sample accessions and download xml file
nSamps=$(cat $output/all_sample_accessions.txt | wc -l) 
currSamp=1
while read sampacc
do
	
let "tempPerc=($currSamp*100/$nSamps)"
let "endBar=($tempPerc*20/100)"
echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"

if [ ! -e "$output/xml_meta/${sampacc}.xml" ]
then
#	echo "$sampacc exists"
#else
curl -s -C - -X POST "https://www.ebi.ac.uk/ena/browser/api/xml?accessions=${sampacc}&download=true" -H "accept: application/xml" -o $output/xml_meta/${sampacc}.xml
fi

let "currSamp++"

done < $output/all_sample_accessions.txt



##### Remove any sample accessions that are NOT in list of studies ####
for f in $output/xml_meta/*
do
sampsearch=$(echo $f | sed 's/\.xml//g' | sed 's/^.*\///g')

if grep -q $sampsearch $output/all_sample_accessions.txt
then
echo -ne "\r Checking things to keep"
else
echo "remove $sampsearch"
rm $output/xml_meta/$sampsearch.xml
fi

done 

echo ''
echo 'DONE'


################################################################################################
################################################################################################

