# merge together all the csv files, one for support vectors and the other for coefficients
# this script goes into each of the deconvolved sample directories
# and merges the coefficient files and the support vector files
# the best coefficient will then be chosen locally 

import pandas as pd
import numpy as np
import scanpy as sc

import os

# dataset name
datasetName = "dec_cfRNAmeta" # name of dataset
fend = "_results" # name of file ending...

# path to your data
sampData = "sample_cpm.csv"
basePath = os.getcwd() + "/"

#### do not edit below this line

# read in samples
samps = pd.read_csv(sampData, sep=",", index_col=0, nrows=0)
sampNames = samps.columns.tolist()
os.system("echo DONE")
supVecBase = "_coarsegrain___supportVectors.csv"
coefBase = "_deconvolutionCoefs.csv"
os.system("echo DONE")
allSuppVec = pd.DataFrame()
allCoefs = pd.DataFrame()
os.system("echo DONE")
# concatenate all the samples over all hyperparam combinations
for bioRep in sampNames:
    os.system("echo BIOREP")
    sampSV = pd.read_csv(basePath + bioRep + "/" + bioRep +  supVecBase, sep = ",", index_col = 0)
    sampCoef = pd.read_csv(basePath + bioRep + "/" + bioRep + coefBase, sep = ",", index_col = 0)
    allCoefs = pd.concat([sampCoef, allCoefs], ignore_index = False)
    allSuppVec = pd.concat([sampSV, allSuppVec], axis = "columns", ignore_index = False)
os.system("echo DONE")
# pick the best hyperparam combination per sample
bestCoef = pd.DataFrame()

for samp in sampNames:
    print(samp)
    allHyperThisSamp = [i for i in allCoefs.index if samp in i]
    hyperCombosThisSamp = allCoefs.loc[allHyperThisSamp]
    best = hyperCombosThisSamp[hyperCombosThisSamp['rmse'] == np.min(hyperCombosThisSamp['rmse'])]
    bestCoef = pd.concat([bestCoef, best])

performance = bestCoef.iloc[:,-2:].T
bestCoef = bestCoef.iloc[:,:-2]

# send negative coefs to zero
bestCoef[bestCoef < 0] = 0 # this is a samp x cells matrix
bestCoef = bestCoef.T

# normalize by total to get the fractions of cell type specific RNA
fracs = bestCoef.div(bestCoef.sum(axis = 0), axis = 1) 

fracs = pd.concat([fracs, performance], axis = "index")

# get the support vectors corresponding to the best coefficient pair
suppvecs = allSuppVec[bestCoef.columns.tolist()]

# strip the hyperparameter information so it's just the sample names
bestCoef.columns = [i.split("-NUSVR")[0] for i in bestCoef.columns.tolist()]
suppvecs.columns = [i.split("-NUSVR")[0] for i in bestCoef.columns.tolist()]

# write out the support vectors and the fractions
suppvecs.to_csv(datasetName  + "support_vectors_" + fend + ".csv", sep = ",", header = True, index = False)
fracs.to_csv(datasetName + "_fractions" + fend + ".csv", sep = ",", header = True, index = True)

