#In[0]

from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
import numpy as np
import time
from operator import itemgetter
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform


#In[1]

def GENIE3(expr_data,gene_names,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,nthreads=1):

    '''Terra based analysis
    
    INPUT
    expr_data:  gene expression values of single-cell data
    gene_names: List of features (number of columns in expr_data), containing the names of the genes
    regulators: TF
    tree-method: 'RF'- random forest 
    K: 'sqrt', 'all' or a positive integer- optional category: number of selcted attributes
    ntrees: number of trees grown in an ensemble.
    nthreads: number of threads used for parallel computing
    
    '''

    print('Tree method: {}, K: {}, Number of trees: {} \n'.format(tree_method,K,ntrees),flush= True)
    ngenes = expr_data.shape[1]

    # Get the indices of the candidate TF
    if regulators == 'all':
        input_idx = list(np.arange(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]

    # Learn an ensemble of trees (RF) for each target gene, and compute scores for candidate TF
    VIM = np.zeros((ngenes,ngenes))
    if nthreads > 1:
        
        input_data = list()
        for i in range(ngenes):
            input_data.append( [expr_data,i,input_idx,tree_method,K,ntrees] )

        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_single, input_data)

        for (i,vi) in alloutput:
            VIM[i,:] = vi

    else:
        for i in range(ngenes):
            if (i+1)%20==0:
                print('Gene %d/%d...' % (i+1,ngenes))

            vi = GENIE3_single(expr_data,i,input_idx,tree_method,K,ntrees)
            VIM[i,:] = vi

    return np.transpose(VIM)


#In[2]
def wr_GENIE3_single(args):
    return([args[1], GENIE3_single(args[0], args[1], args[2], args[3], args[4], args[5])])


def GENIE3_single(expr_data,output_idx,input_idx,tree_method,K,ntrees):
    ngenes = expr_data.shape[1]

    # output expression of target gene
    output = expr_data[:,output_idx]

    # Normalize 
    if np.std(output) == 0:
        output = np.zeros(len(output))
    else:
        output = output / np.std(output)

    # Removal of target gene from candidate TF
    input_idx = input_idx[:]
    if output_idx in input_idx:
        input_idx.remove(output_idx)

    expr_data_input = expr_data[:,input_idx]

    if (K == 'all') or (isinstance(K,int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K

    if tree_method == 'RF':
        treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=max_features)
    elif tree_method == 'ET':
        treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=max_features)

    # Learning ensemble of trees
    treeEstimator.fit(expr_data_input,output)

    #importance scores
    importances = [e.tree_.compute_feature_importances(normalize=False) for e in treeEstimator.estimators_]
    importances = np.asarray(importances)
    feature_importances = np.sum(importances,axis=0) / len(treeEstimator)

    vi = np.zeros(ngenes)
    vi[input_idx] = feature_importances

    return vi


#In[3]
def prepare_data_for_GENIE3(adata, cell_type, transcription_factors):
    """
    Prepare the expression data and transcription factor list for GENIE3 from an AnnData object for a specific cell type.
    INPUT
    - adata: AnnData object containing single-cell expression data.
    - cell_type: The cell type for which to infer the regulatory network.
    - transcription_factors: A list of transcription factor names to be considered as potential regulators.

    OUTPUT
    - expr_data: A numpy array of expression data for the specified cell type.
    - gene_names: A list of gene names corresponding to the columns in expr_data.
    - regulators: A list of transcription factor names filtered to include only those present in gene_names.
    """

    # Filtering the AnnData object for the specified cell type
    adata_filtered = adata[adata.obs['celltypes'] == cell_type]

    # Extracting the gene expression matrix and convert to numpy array
    expr_data = adata_filtered.X
    if isinstance(expr_data, np.ndarray) == False:
        expr_data = expr_data.toarray()  

    # Extracting gene names from anndata
    gene_names = list(adata_filtered.var_names)

    # Filtering the list of transcription factors to those present in the gene_names list- Animal TFDB
    regulators = [tf for tf in transcription_factors if tf in gene_names]

    return expr_data, gene_names, regulators


#In[4]
# importing libraries
import numpy as np
import scanpy as sc
import scipy as sci
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import autogenes as ag

from sklearn.svm import NuSVR
import pickle
import h5py
import anndata as ad


#In[5]


# Reading adnndata fike
adata = ad.read_h5ad('/work/ABG/MUSKAN/Deconvolution/AutoGeneS/test_muskan/data/SC-7/sc.h5ad')

# The name of the cell type you're interested in
cell_type = 'CD4+ ab T cells'  

# Path to the file containing the list of transcription factors
transcription_factors_file = '/work/ABG/MUSKAN/Terra-Test/TF.txt'

# Read the file and create a list
with open(transcription_factors_file, 'r') as file:
    transcription_factors = [line.strip() for line in file]


# calling the function
expr_data, gene_names, regulators = prepare_data_for_GENIE3(adata, cell_type, transcription_factors)

#In[6]

# Now you can use the expr_data, gene_names, and regulators with the GENIE3 function
results = GENIE3(expr_data, gene_names=gene_names, regulators=regulators, tree_method='RF', K='sqrt', ntrees=1000, nthreads=1)

#saving the output in a .txt file
output= "/work/ABG/MUSKAN/Terra-Test/grn-results-12Aug'24.txt"
np.savetxt(output, results, delimiter= '\t')

