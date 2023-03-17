#!/bin/bash python3

import argparse
import numpy as np # for expanding array
import re # for splitting via variable
from random import choices
from random import sample
import matplotlib.pyplot as plt
import os # to make directories
import math # for rounding
from progressbar import ProgressBar # Progress tracker
from time import sleep

########## Set up ##########
parser = argparse.ArgumentParser(description = "For conducting coverage-based rarefaction. Finds minimum coverage amongst samples and rarefies all samples to that depth")
parser.add_argument('-i', '--input', type=str
					, help='Text file of species x sample table. Tab or space delimited')
parser.add_argument('-d', '--delimiter', type=str
					, help='Delimiter for species table [ Default: tab ]', default="\t")
parser.add_argument('-m', '--min', type=int
					, help='Minimum rarefaction depth to try [ Default: 100 ]', default=10)
parser.add_argument('-M', '--max', type=int
					, help='Maximum rarefaction depth to try [ Default: 10000 ]', default=10000)
parser.add_argument('-b', '--by', type=int
					, help='Size between rarefaction steps [ Default: 500 ]', default=10)
parser.add_argument('-r', '--reps', type=int
					, help='Number of replications per rarefaction depth [ Default: 10 ]', default=10)
parser.add_argument('-c', '--coverage_min', type=float
					, help='Minimum coverage for a sample to not be thrown out [ Default: 0.85 ]', default=0.85)
parser.add_argument('-o', '--output', type=str
					, help='Directory to save final rarefied species table', default="rarefied_table.txt")				

args = parser.parse_args()
inputFP = args.input
d = args.delimiter
m = args.min
M = args.max
b = args.by
r = args.reps
c = args.coverage_min
outputFP = args.output

##### Functions ####


		
class speciestable:
	table=[]
	sampleList=[]
	sampleDepths=[]
	
	def __init__(self, table):
		self.table=table
		self.sampleDepths=[sum([float(i) for i in x[1:]]) for x in self.table[1:]]
		
	class sample:
		name=''
		speciesList=''
		def __init__(self, name, speciesList):
			self.name=name
			self.speciesList=speciesList
			# Covert to list
		def add_rarefied_sample(self, rarefiedSample):
			self.rarefiedSample=rarefiedSample
			
	def add_all_samples(self):
		allSamps = self.table[1:]
		species = self.table[0][1:]
		spArray = np.array(species)
		for s in allSamps:
			currSamp = s[0]
			self.sampleList.append(currSamp)
			currAbund = [float(x) for x in s[1:len(s)]]
			currList = np.repeat(spArray, currAbund, axis=0)
			setattr(self, currSamp, self.sample(currSamp, currList))
			
	def add_new_sample(self, name, speciesList):
		newSamp = self.sample(name, speciesList)
		setattr(self, name, newSamp)
		# Add checking if all species are present
		
	def generate_coverage_range_depths(self, reps, nmin, nmax, by):
		self.coverageTables = np.array([])
		allMats = []
		self.sampling_depths = list(range(nmin, nmax, by))
		pbar = ProgressBar()
		for n in pbar(list(range(nmin, nmax, by))):
			sleep(0.2)
			allMats.append(self.generate_coverage_single_depth(n, reps))
		self.coverageTables = np.stack(allMats, axis=1)
		
	def generate_coverage_single_depth(self, n=1, reps=1, already_rarefied = 'False'):
		cover_mat = np.array([ [0.0]*reps ]*(len(self.sampleList)))
		for sn,s in enumerate(self.sampleList):
			for r in range(reps):
				if already_rarefied == 'False':
					samp = choices(list(getattr(self, s).speciesList),k=n)
				elif already_rarefied == 'True':
					tempsamp = getattr(self,s).rarefiedSample
					samp = []
					for s in list(tempsamp.keys()):
						samp+=[s]*int(tempsamp[s])
				cover_mat[sn][r] = self.func_calculate_coverage(samp)
		return(cover_mat)
	
	def remove_low_coverage_samples(self, coverage_min=0):
		idx_max_sampling_depths = [max(i for i,x in enumerate(self.sampling_depths) if x < D) for D in self.sampleDepths]
		ave_cover_mat = np.average(self.coverageTables, axis=2)
		#toKeep = list(ave_cover_mat[:,maxIndex] < coverage_min)
		toKeep = [j>=coverage_min for j in [ave_cover_mat[i,x] for i,x in enumerate(idx_max_sampling_depths)] ]
		ndiscard = len(toKeep)-sum(toKeep)
		self.coverageTables = self.coverageTables[toKeep]
		self.sampleList = [ x for i,x in enumerate(self.sampleList) if toKeep[i]]
		print("REMOVED LOW COVERAGE SAMPLES:"+str(ndiscard))
		
	def rarefy(self, custom_depth = 'automatic'):
		if custom_depth=='automatic':
			print('AUTOMATIC RAREFACTION DEPTHS')
			depths = self.func_calculate_rarefaction_depth()
		else:
			depths = custom_depth
		for idx,samp in enumerate(self.sampleList):
			templist = sample(list(getattr(self, samp).speciesList), depths[idx])
			rarefiedSample = {x : templist.count(x) for x in set(templist) }
			getattr(self, samp).add_rarefied_sample(rarefiedSample)
		self.depths = depths
		
	def make_rarefied_table(self, sampleList = 'all'):
		if sampleList == 'all':
			sampleList = self.sampleList
		listAllSpecies = []
		for samp in sampleList:
			listAllSpecies += list(getattr(self, samp).rarefiedSample.keys())
		self.listRarefiedSpecies = list(set(listAllSpecies))
		rarefiedTable=np.zeros((len(sampleList), len(self.listRarefiedSpecies)))
		for idx,samp in enumerate(sampleList):
			tempSp = np.zeros((len(self.listRarefiedSpecies)))
			tempSamp = getattr(self, samp).rarefiedSample
			for key in tempSamp:
				pos = [i for i,x in enumerate(self.listRarefiedSpecies) if key == x]
				tempSp[pos] = tempSamp[key]
			rarefiedTable[idx,] = tempSp
		self.rarefiedTable = rarefiedTable.transpose()
			
	def func_calculate_coverage(self, samp):
		counts = [[i, samp.count(i)] for i in set(samp)]
		f1 = float(sum([1 for i in counts if i[1]==1]))
		f2 = float(sum([1 for i in counts if i[1]==2]))
		n = float(len(samp)+1)
		if (( (n-1)*f1 + 2*f2) == 0) :
			coverage = 1
		else:
			coverage = 1 - (f1/n)*( ( (n-1)*f1 ) / ( (n-1)*f1 + 2*f2) )
		return(coverage)
		
	def func_calculate_rarefaction_depth(self):
		maxIndex = len(self.sampling_depths)-1
#		print(self.coverageTables)
		ave_cover_mat = np.average(self.coverageTables, axis=2)
		lastRare = list(ave_cover_mat[:,maxIndex])
		minCoverage = math.floor(min(lastRare)*100)/100 # Get nearest percent coverage
		# Get rarefaction depth for all samples
		diff_from_target_mat = abs(minCoverage-ave_cover_mat)
		index_best_depth = [ list(x).index(min(list(x))) for x in diff_from_target_mat ]
		best_depths = [ self.sampling_depths[i] for i in index_best_depth ]
		beyondDepth = [ i for i,x in enumerate(best_depths) if x>self.sampleDepths[i] ]
		if len(beyondDepth)==0:
			final_rare_depths = best_depths
		else:
			print("WARNING: some samples have less depth than recommended coverage score. Adjusting coverage to value that includes all sampling depths")
			# Get new minimum coverage
			min_depths_temp = [max_depths_allowed[i] for i in beyondDepth ]
			diff_list = [[i-x for x in self.sampling_depths] for i in min_depths_temp ]
			minDiff_temp = [min([x for x in v if x>=0]) for v in diff_list]
			index_new_coverage = [diff_list[i].index(m) for i,m in enumerate(minDiff_temp)]
			# Get new coverages of samples that have depth less than required sampling depth
			new_min_coverage = min(list(ave_cover_mat[beyondDepth,index_new_coverage]))
			###### Now, re-calculate minimum coverage.
			diff_from_target_mat = abs(new_min_coverage-ave_cover_mat)
			index_best_depth = [ list(x).index(min(list(x))) for x in diff_from_target_mat ]
			final_rare_depths = [ self.sampling_depths[i] for i in index_best_depth ]
		return(final_rare_depths)	
	
	def plot_coverage_curves(self, outputdir):
		self.func_mkdir(outputdir)
		ave_cover_mat = np.average(self.coverageTables, axis=2)
		for i in ave_cover_mat:
			plt.plot(self.sampling_depths, i)
		plt.xlabel('Sampling Depth')
		plt.ylabel('Estimate coverage')
		plt.savefig(outputdir+'/coverage_curves.png', transparent=False)
		plt.clf()
	
	def print_coverage_table(self, outputdir):
		self.func_mkdir(outputdir)
		ave_cover_mat = np.round_(np.average(self.coverageTables, axis=2),3)
		ave_cover_mat_wheaders = np.vstack((self.sampling_depths, ave_cover_mat))
		toPrint = 'depth\t'+'\t'.join(self.sampleList) + '\n'
		for i in range(len(ave_cover_mat_wheaders[0])):
			toPrint += '\t'.join([str(x) for x in ave_cover_mat_wheaders[:,i]])+'\n'
		toPrint = toPrint.strip()
		toWrite = open(outputdir+'/rarefied_coverage.txt','w')
		toWrite.write(toPrint)
		toWrite.close()
		
	def print_rarefied_table(self, outputdir):
		self.func_mkdir(outputdir)
		self.def_print_table(self.rarefiedTable, self.listRarefiedSpecies, self.sampleList, outputdir+'/rarefiedTable.txt')
		
	def print_depth_and_coverage(self, outputdir):
		cov = self.generate_coverage_single_depth(already_rarefied='True')
		tab_depth_cov = np.vstack((self.depths,cov[:,0])).transpose()
		self.def_print_table(tab_depth_cov, rown=self.sampleList, coln=['rare_depth','coverage'], output=outputdir+'/depth_and_coverage.txt')
		
	def def_print_table(self, arr, rown, coln, output):
		toPrint = 'X\t' + '\t'.join(coln) + '\n'
		for i,r in enumerate(arr):
			toPrint += rown[i] + '\t' + '\t'.join(str(x) for x in list(r)) + '\n'
		f = open(output, 'w')
		f.write(toPrint)
		f.close()
	
	def func_mkdir(self, newdir):
		if not os.path.exists(newdir):
			os.makedirs(newdir)
		
		 ######## TO BE CONTINUED

    		
# Functions #
def import_table(inputFP):
	dat = open(inputFP, 'r')
	contents = [re.split(d,i.strip()) for i in dat.readlines()]
	dat.close()
	byCol = list(zip(*contents))
	return(byCol)

#
# TESTING AREA
#inputFP='05_run_QTAG/campbell2013/OTUTable_counts.txt'
#d='\t'
#outputFP='coverage_testing'
#r = 3
#m = 100
#M = 300
#b = 10

########## RUN #################
print("Importing table")
TABLE = speciestable(import_table(inputFP))
TABLE.add_all_samples()

print("Generating rarefaction curves for coverage")
TABLE.generate_coverage_range_depths(reps=r, nmin=m,nmax=M,by=b)
TABLE.remove_low_coverage_samples(coverage_min=c)
print("Rarefying")
TABLE.rarefy()
TABLE.make_rarefied_table()
print("Making tables and plots")
TABLE.print_rarefied_table(outputFP)
TABLE.plot_coverage_curves(outputFP)
TABLE.print_coverage_table(outputFP)
TABLE.print_depth_and_coverage(outputFP)

print("\n\n DONE")



