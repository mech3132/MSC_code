#!bin/bash

###### Script to download fasta files from EBI/ENA given a list of accession IDs ##########

################################################################################################
#### Help ####
################################################################################################

Help()
{
   # Display Help
   echo
   echo
   echo "This script takes a list of ENA/EBI sample accessions and downloads all the associated fastq files. First input should be a .txt file where each line is an accession number. Second input should be directory to save the fastq files into"
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
#### Setup ####
################################################################################################

BAR='####################'


# Troubleshotting #
#input='00_study_information/all_sample_run_accessions_short.txt'
#input=$1

# Check if output directory input exists
if [ ${#2} == 0 ]
then
output="fastq"
else
output="$2"
fi

# Make output folder
mkdir -p $output

nSamps=$(cat $1 | wc -l) 

currSamp=1



################################################################################################
#### Download files ####
################################################################################################


########################
while read sampacc
do
let "tempPerc=($currSamp*100/$nSamps)"
let "endBar=($tempPerc*20/100)"
echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"
# Get folders
fold1=${sampacc:0:6}
lastDig=${sampacc: 9}

if (( ${#lastDig}==1 ))
then
	nzeros=00
else
	nzeros=0
fi
# Testing
#echo "$fold1/$nzeros$lastDig/$sampacc"



if [ -e $output/$sampacc ] # check if sampacc folder exists
then
	echo "$sampacc exists"
else
	# Try both; for some reason they might be different?
	echo "$sampacc downloading..."
	wget -q -P $output/$sampacc/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$fold1/$nzeros$lastDig/$sampacc/*
	wget -q -P $output/$sampacc/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$fold1/$sampacc/*
# Store un-downloaded files here
: > $output/MANUAL_DOWNLOAD.txt
if [ ! -e $output/$sampacc ]
then
echo "Saving to manual download"
echo $sampacc >> $output/MANUAL_DOWNLOAD.txt
fi # manual downloade

fi # check if sampacc folder exists


let "currSamp++"
done < $1
####################


echo "DONE DOWNLOADING"

