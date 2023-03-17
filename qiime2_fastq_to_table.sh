#!bin/bash

######## This script takes a list of fastq file pathways and "merges" them into one ASV table ######
 
. activate qiime2-2021.11



################################################################################################
#### Help ####
################################################################################################

Help()
{
   # Display Help
   echo
   echo
   echo "This script takes a list of fastq pathways (one file== one sample) and:
   	(a) imports them as qza files
   	(b) denoises them with deblur
   	(c) merges them into one file
   	
   	Options:
   	[input .txt file] [output]
   	"
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

if [ ${#2} == 0 ];
then
	output='qiime2_process'
else
	output="$2"
fi

input=$1


################################################################################################
#### Setup ####
################################################################################################


mkdir -p $output
mkdir -p $output/all_manifests
mkdir -p $output/demux
mkdir -p $output/deblur
mkdir -p $output/tables
mkdir -p $output/repsets

: > $output/changed_accession_names.txt # for sample names that get changed because of _ - rule


# Progress bar
BAR='####################'

nSamps=$(cat $input | wc -l) 

currSamp=1
################################################################################################
#### Start ####
################################################################################################


while read samp_path;
do

# Progress bar
let "tempPerc=($currSamp*100/$nSamps)"
let "endBar=($tempPerc*20/100)"
echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"


# Get rid of before, suffix, and _1 if F/R present
sampacc=$(echo $samp_path | sed 's/^.*\///g' | sed 's/\.fastq\.gz//g' | sed 's/_1$//g')

### Run QIIME stuff

#### Check if imported; import if not

### For deblur, can't have underscores
sampacc2=$(echo $sampacc | sed 's/_/-/g')
# Save filenames that are not the same so they can be mapped later
echo -e "$sampacc\t$sampacc2" >> $output/changed_accession_names.txt


if [ -e $output/demux/$sampacc2.qza ]
then
echo "Already imported $sampacc"

else

echo -e "sample-id\tabsolute-filepath" > $output/all_manifests/${sampacc2}_manifest.txt
echo -e "$sampacc2\t$samp_path" >> $output/all_manifests/${sampacc2}_manifest.txt


# Import R1 as single demultiplexed sample
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path $output/all_manifests/${sampacc2}_manifest.txt \
--output-path $output/demux/$sampacc2.qza \
--input-format SingleEndFastqManifestPhred33V2 

fi # deblur-demux exists?



if [ -e $output/tables/$sampacc2-table.qza ] 
then
	echo "deblur already run for $sampacc"

else

# Re-name underscored sample names-- "deblur cannot operate on sample IDs that contain underscores"

# Run deblur -- default filter metrics, but included here for clarity
qiime deblur denoise-16S \
--i-demultiplexed-seqs $output/demux/$sampacc2.qza \
--p-trim-length 150 --p-min-reads 10 \
--p-min-size 2 \
--p-no-hashed-feature-ids \
--output-dir $output/deblur/$sampacc2 \
--quiet


mv $output/deblur/$sampacc2/table.qza $output/tables/$sampacc2-table.qza
mv $output/deblur/$sampacc2/representative_sequences.qza $output/repsets/$sampacc2-repset.qza

fi


#if [ -e allRawData/$sampacc/${sampacc}_2.fastq.gz ]
#then
#rm -r allRawData/$sampacc/${sampacc}_2.fastq.gz ######### REHIGHLIGHT when doing full batch to save HD space
#echo "Removed extra reverse fastq file $sampacc"
#fi


let "currSamp++"
done < $input




