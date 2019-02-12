# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 22:45:42 2018

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
import os

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')

results = pd.read_csv("results_bpass_agn_he2sel2/RESOLVE_param_estimates.csv")
inputfile = 'C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl'
#dat = Table.read(inputfile, format='fits')
#infile = dat.to_pandas()
infile = pd.read_pickle(inputfile)

Z_index = (results['Parameter'] == 'LOGZ')
results_Z = results[Z_index]
results_Z.index = results_Z['Galaxy Name']

agn_index = (results['Parameter'] == 'AGNFRAC')
results_agn = results[agn_index]
results_agn.index = results_agn['Galaxy Name']

results_Z['CI68_low'] = results_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
results_Z['CI68_high'] = results_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)

pd.set_option('display.max_columns', 500)
outliers = results_Z['Estimate'] == min(results_Z['Estimate'])
print results_agn[outliers]['Estimate'], results_Z[outliers]['CI68_low'], results_Z[outliers]['CI68_high']
bpt = pd.read_csv(r'C:\Users\mugdhapolimera\github\BPT\BPT\BPT\resolve_emlineclass_new.csv')
ambig = bpt.agntosf | bpt.ambigagn | bpt.composite | bpt.sftoagn

ambig_ndx = [results['Galaxy Name'][x] in list(bpt['galname'][ambig]) for x in range(len(results))]
