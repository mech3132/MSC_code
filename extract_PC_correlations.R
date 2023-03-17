#!bin/bash

######## Code to extract variables that are strongly correlated with a PC in a PCA #####
# Load tidyverse (not technically required, but I use it to manipulate the data)
# library(tidyverse)

# Load table where rows are ASVs or genes and columns are samples
dat <- read.csv("raw_data/alldataDavid.csv", row.names = 1)
# Get rid of rows where ASV/gene sums are zero
dat <- dat[rowSums(dat)!=0,]
# Transpose so that rows are samples instead
datt <- t(dat)

### Make PCA
pca_datt <- prcomp(datt)

#### To calculate variance explained by each PC (principal component),
# divide variance of each axis by the sum of all variance
# Variance is sdev^2
var_explained_by_each_PC <- pca_datt$sdev^2/sum(pca_datt$sdev^2)

#### To calculate how much of each PC is explained by each ASV/gene,
# Square the "rotation" matrix
var_within_PCs <- pca_datt$rotation^2
# Note: if you sum var_within_PCs by column (sum by PC), you should get 1
colSums(var_within_PCs)
# because 100% of the variation in PC1 is explained by all your predictors


### To calculate the total variance in the dataset explained by each ASV/gene
### split up by principal components,
# you can use matrix multiplcation to multiply variance explained by each PC
# by variance explained by each ASV/gene within each PC
total_var_explained_by_predictors <- t(var_within_PCs) * var_explained_by_each_PC

# Summing this entire matrix should yield 1 because all the variance
# should be explained by adding together all components from all ASVs/genes
sum(total_var_explained_by_predictors)

#### Now, there are two ways to ask the data which ASVs/genes are important for each PC
## Option 1: What is the gene that explains the most variation in each PC?
# First, rank ASVs/genes within each PC
predictors_ranked_within_PC <- apply(var_within_PCs,MARGIN=2, function(x) rank(-x) )
# next, choose the PC you want to look at. 
# Let's look at the top 10 explanatory ASV/genes in PC1
top10ASVs_PC1 <- names(sort(predictors_ranked_within_PC[,"PC1"])[1:10])
top10ASVs_PC1
# How much variation do these genes explain of PC1?
var_within_PCs[top10ASVs_PC1,"PC1"]



## Option 2: What PC does each ASV/gene correlate with strongest?
# Rank each PC within each ASV/gene by their "loading"
PC_rankings_by_predictor <- t(apply(var_within_PCs, MARGIN = 1, function(x) rank(-x)))
# Choose a PC to look at; let's look at PC1
names(which(PC_rankings_by_predictor[,"PC1"]==1))
# Notice that there are many ASV/genes with rank "1"-- 
# That is because multiple genes can correlate most strongly with PC1

#### Finally, are there any genes that are the best predictor of a PC 
### AND that PC is also its strongest correlate?
mat1 <- predictors_ranked_within_PC
mat1[mat1!=1] <- 0
mat2 <- PC_rankings_by_predictor
mat2[mat2!=1] <- 0
mat1x2 <- mat1*mat2
mat1x2 <- mat1x2[rowSums(mat1x2)>0,]
matfinal  <- mat1x2[,colSums(mat1x2)>0]
matfinal
