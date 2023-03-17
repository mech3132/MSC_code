#!/bin/bash python3

####### IMPORTS #########

import argparse
from PIL import Image
import numpy as np  #for expanding array
#import re  for splitting via variable
#from random import choices
#from random import sample
import matplotlib.pyplot as plt
import os  # to make directories
import colorsys # for converting to hsv?
#import math  for rounding
#from progressbar import ProgressBar  Progress tracker
#from time import sleep

########## Set up ##########
parser = argparse.ArgumentParser(description = "For conducting coverage-based rarefaction. Finds minimum coverage amongst samples and rarefies all samples to that depth")
parser.add_argument('-i', '--input', type=str
					, help='File path to image to be imported')
#parser.add_argument('-r', '--red_threshold', type=int
#					, help='Threshold of color intensity of red (max 255) [ Default: 100]', default=0.5)
parser.add_argument('-g', '--green_threshold', type=float
					, help='Threshold of color intensity ratio of green to total (max 1) [ Default: 0.36', default=0.36)
#parser.add_argument('-b', '--blue_threshold', type=int
#					, help='Threshold of color intensity of blue (max 255) [ Default: 100]', default=0.5)
parser.add_argument('-s', '--save_subimages', type=bool
					, help='Save sub-images?', default=True)
parser.add_argument('-o', '--output', type=str
					, help='Directory to save final table and images', default="RGB_processing")				

args = parser.parse_args()
inputFP = args.input
#red_threshold = args.red_threshold
green_threshold = args.green_threshold
#blue_threshold = args.blue_threshold
save_subimages = args.save_subimages
outputFP = args.output






class image():
	def __init__(self, inputFP):
		# Load in image
		# Create classes A1 through H12
		# set thresholds
		self.im = Image.open(inputFP)
		self.imarray = np.array(self.im)
	
	def make_custom_mask(self, maskname, maskarray):
		## set image mask
		imarray_zeros = np.zeros(self.imarray.shape, dtype='uint8')
		imarray_zeros[maskarray,:] = [255,255,255]
		# For visualization
		setattr(self, str(maskname), imarray_zeros)
	
	def array_to_image(self, attr):
		temp = getattr(self,attr)
		return(Image.fromarray(temp))
	
	def divide_image(self, imagearray,  colNames, rowNames):
		self.colNames = colNames
		self.rowNames = rowNames
		ncol=len(colNames)
		nrow=len(rowNames)
		rpix = imagearray.shape[0]
		cpix = imagearray.shape[1]
		# n pix per group
		rwidth = rpix/nrow
		cwidth = cpix/ncol
		# Get ranges
		cranges =[[int(round(i*cwidth)),int(round((i+1)*cwidth-1))] for i in range(ncol)]
		rranges =[[int(round(i*rwidth)),int(round((i+1)*rwidth-1))] for i in range(nrow)]
		allSubimageNames = []
		for nr,r in enumerate(rranges):
			for nc,c in enumerate(cranges):
				subimage = imagearray[r[0]:r[1], c[0]:c[1],:]
				subimageName = str(colNames[nc])+str(rowNames[nr])
				setattr(self, subimageName, subimage)
				allSubimageNames.append(subimageName)
		self.subimageNames = allSubimageNames
				
	def extract_pixels_all_subimages(self, function, params):
		# If no attribute with list of masks, set here
		if hasattr(self, 'list_pixel_data'):
			getattr(self, 'list_pixel_data').append(function.__name__)
		else:
			newlist = [function.__name__]
			setattr(self, 'list_pixel_data', newlist)
		temparray = np.zeros([len(self.subimageNames), 2]).tolist()
		for i,s in enumerate(self.subimageNames):
			temp = getattr(self, s)
			filt = function(temp, params)
			npixels = sum(filt.flatten())
			temparray[i] = [s,npixels]
		tempname = 'pixeldata_' + function.__name__
		setattr(self, tempname, temparray)
	
	def print_subimages(self,output):
		try:
  			os.makedirs(outputFP+'/subimages')
		except FileExistsError:
  			# directory already exists
 			 pass
		for s in self.subimageNames:
			tempsubimage = getattr(self,s)
			Image.fromarray(tempsubimage).save(output+'/subimages/'+s+'.png')
	
	def print_pixeldata(self, toprint, output):
		toPrint = 'subimageName\tpixels\n'
		for r in getattr(self, 'pixeldata_'+toprint):
			toPrint += '\t'.join([str(i) for i in r]) + '\n'
		toPrint.strip()
		f = open(output+'/pixeldata_'+toprint+'.txt', 'w')
		f.write(toPrint)
		f.close()


def green_ratio_mask(imagearray, green_threshold):
	## Make mask
	arrayr = np.array(imagearray[:,:,0], dtype=float)
	arrayg = np.array(imagearray[:,:,1], dtype=float)
	arrayb = np.array(imagearray[:,:,2], dtype=float)
	# Calculate ratio of green to all pixels
	# Get single rgb arrays
	allpix = arrayg+arrayr+arrayb
	gratio = (arrayg+1)/(allpix +1)
	# Set mask
	gmask = gratio>green_threshold
	return(gmask)

########## UNDERCONTRUCTIONS
def gray_ratio_amsk(imagearray, tolerance):
	## Make mask
	# Single arrays
	arrayr = np.array(imagearray[:,:,0], dtype=float)
	arrayg = np.array(imagearray[:,:,1], dtype=float)
	arrayb = np.array(imagearray[:,:,2], dtype=float)
	# Calculate pariwise ratio of all colours to get 'grey' values
ratioRG = arrayg/arrayr 
ratioGR = arrayr/arrayg 
ratiogr = np.zeros(imobj.imarray.shape[0:2], dtype='uint8')

[i for i in ratioRG[list(j) for j,x in enumerate(ratioRG)]]
[min() for i,j in mat]

[min(X, X) for i in zip(ratioRG[j for j,x in enumerate(ratiogr),:], ratioGR[j for j,x in enumerate(ratiogr),:])]

ratioRB = arrayr/arrayb
ratioBR = arrayb/arrayr
ratioGB = arrayg/arrayb
ratioBG = arrayb/arrayg



	
	
	
	
	allpix = arrayg+arrayr+arrayb
	gratio = (arrayg+1)/(allpix +1)
	# Set mask
	gmask = gratio>green_threshold
	return(gmask)


########## UNDERCONTRUCTIONS



### Process athal scans
#PATH
#inputFP='~/Documents/PostDoc_UBC/Research/Project_highthroughput/processing_scans/images/test_chmyc001_cropped.jpg'
##red_threshold
#green_threshold = 0.36
#output='~/Documents/PostDoc_UBC/Research/Project_highthroughput/processing_scans/images/process'

inputFP='test_chmyc001_cropped.jpg'
inputFP='scans/polymyxa_plate002_cropped.tif'
#red_threshold
green_threshold = 0.36
outputFP='process'


## VARIABLES

colNames = ['H', 'G', 'F', 'E', 'D', 'C', 'B','A']
rowNames = ['1','2','3','4','5','6','7','8','9','10','11','12']

#blue_threshold

# GHOST CODE
# Create and load image
try:
   os.makedirs(outputFP)
except FileExistsError:
   # directory already exists
   pass
   
imobj = image(inputFP)
imobj.divide_image(imobj.imarray, colNames, rowNames)
imobj.extract_pixels_all_subimages(green_ratio_mask, green_threshold)

# Print subimages
if save_subimages:
	imobj.print_subimages(outputFP)
# Print table
imobj.print_pixeldata('green_ratio_mask',outputFP)


# Print b&w for double check
bwimagearray = np.zeros(imobj.imarray.shape, dtype='uint8')
mask = green_ratio_mask(imobj.imarray, green_threshold)
bwimagearray[mask,:] = [255]*imobj.imarray.shape[2]
bwimage = Image.fromarray(bwimagearray)
bwimage.save('bw_image.png')

# Save list of pixels


#imshape = imobj.imarray.shape
#im_long = imobj.imarray.reshape(imshape[0]*imshape[1],imshape[2])
#tf_im_long = mask.reshape(imshape[0]*imshape[1],1)
#test = [i.tolist() for (i, v) in zip(im_long, tf_im_long) if v]

#a_file = open("im_pixel_long.txt", "w")
#for row in im_long:
#    np.savetxt(a_file, row)
#a_file.close()
### FOR LATER:
# Plans are to be able to load multiple masks
# And then run all of them with a single command
# Make green mask
#green_mask = green_ratio_mask(imobj.imarray, green_threshold)
# load green mask into object
#imobj.make_custom_mask('green_mask', green_mask)




########### Scrastch

#arrayred = np.zeros(imarray.shape, dtype='uint8')
#arrayred[:,:,0] = imarray[:,:,0]
#arraygreen = np.zeros(imarray.shape, dtype='uint8')
#arraygreen[:,:,1] = imarray[:,:,1]
#arrayblue = np.zeros(imarray.shape, dtype='uint8')
#arrayblue[:,:,2] = imarray[:,:,2]

#n=150

#imarray_line = imarray.copy()
#imarray_line[n,:,:] = 0
#Image.fromarray(imarray_line).show()

##imarray_copy = imarray.copy()
#arrayr = imarray[:,:,0]
#arrayg = imarray[:,:,1]
#arrayb = imarray[:,:,2]

#function = green_ratio_mask
#params = green_threshold
#temparray = np.zeros([len(imobj.subimageNames), 2]).tolist()
#for i,s in enumerate(imobj.subimageNames):
#	temp = getattr(imobj, s)
#	filt = function(temp, params)
#	npixels = sum(filt.flatten())
#	temparray[i] = [s,npixels]
#return(temparray)
#		
#		
#		
#		
#ncol = len(colNames)
#nrow = len(rowNames)

#rpix = imagearray.shape[0]
#cpix = imagearray.shape[1]
## n pix per group
#rwidth = rpix/nrow
#cwidth = cpix/ncol
## Get ranges
#cranges =[[int(round(i*cwidth)),int(round((i+1)*cwidth-1))] for i in range(ncol)]
#rranges =[[int(round(i*rwidth)),int(round((i+1)*rwidth-1))] for i in range(nrow)]

#for nr,r in enumerate(rranges):
#	for nc,c in enumerate(cranges):
#		subimage = imagearray[r[0]:r[1], c[0]:c[1],:]
#		subimageName = str(colNames[nc])+str(rowNames[nr])
#		setattr(self, subimageName, subimage)
#		

