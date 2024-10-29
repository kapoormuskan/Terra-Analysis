
# Terra-Analysis
Analysis using Terra workspace for gene regulatory network analysis deploying GENIE3 algorithm for the ingested porcine scRNA-seq data. 
=======
 Building a FAIR data ecosystem for incorporating single-cell transcriptomics data into agricultural genome to phenome research


Pyhton version: 3.12.0
Packages: Genie3

Libraries
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
import numpy as np
import time
from operator import itemgetter
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
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

Run
python script.py



(Initial commit with Terra based anlysis files and an interactive output)
