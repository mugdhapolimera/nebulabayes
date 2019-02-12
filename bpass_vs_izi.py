import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('C:/Users/mugdhapolimera/github/nebulabayes/')

izi = pd.read_csv("results_izi_filtered/RESOLVE_param_estimates.csv")
#bpass = pd.read_csv("results_bpass_gooddata_limited/RESOLVE_param_estimates.csv")
bpass = pd.read_csv("results_izi_bpt1filter/RESOLVE_param_estimates.csv")

Z_index = (bpass['Parameter'] == 'LOGZ')
bpass_Z = bpass[Z_index]
bpass_Z.index = range(len(bpass_Z))
Z_bpass = bpass_Z['Estimate']
Z_bpass.index = bpass_Z['Galaxy Name']

#bpass_Z['CI68_low'] = bpass_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
#bpass_Z['CI68_high'] = bpass_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)

#Z_index1 = (bpass1['Parameter'] == 'LOGZ')
#bpass_Z1 = bpass1[Z_index1]
#bpass_Z1.index = range(len(bpass_Z1))
#Z_bpass1 = bpass_Z1['Estimate']
bpass_Z['CI68_low'] = bpass_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
bpass_Z['CI68_high'] = bpass_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)
errup_bpass = bpass_Z['CI68_high']- bpass_Z['Estimate']
errdown_bpass = bpass_Z['Estimate'] - bpass_Z['CI68_low']

#izi = pd.read_csv("results_izi_gooddata/RESOLVE_param_estimates.csv")

Z_index = (izi['Parameter'] == 'LOGZ')
izi_Z = izi[Z_index]
izi_Z.index = range(len(izi_Z))
Z_izi = izi_Z['Estimate'] 
Z_izi.index = izi_Z['Galaxy Name']
Z_izi = Z_izi.loc[Z_bpass.index.values]
izi_Z['CI68_low'] = izi_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
izi_Z['CI68_high'] = izi_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)

errup_izi = izi_Z['CI68_high'] - izi_Z['Estimate']
errup_izi.index = izi_Z['Galaxy Name']
errup_izi = errup_izi.loc[Z_bpass.index.values]
errdown_izi = izi_Z['Estimate'] - izi_Z['CI68_low']
errdown_izi.index = izi_Z['Galaxy Name']
errdown_izi = errdown_izi.loc[Z_bpass.index.values]
plt.figure()
#plt.plot(Z_izi, Z_bpass, 'o')
plt.plot(np.arange(min(Z_izi), max(Z_izi), 0.01), np.arange(min(Z_izi), max(Z_izi), 0.01), 'g')
plt.errorbar(Z_izi, Z_bpass, xerr = [errdown_izi, errup_izi], yerr = [errdown_bpass, errup_bpass], fmt = 'o')
plt.xlim(7.45, 9.0)
plt.xlabel('NebulaBayes + IZI Default Grid (Good Lines): Z estimate')
plt.ylabel('NebulaBayes + Limited BPASS Grid (Good Lines): Z estimate')

results = np.column_stack((Z_bpass.index.values, Z_bpass, errup_bpass, errdown_bpass))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra')
np.savetxt('RESOLVE_NB_bpt1filter.txt', results, fmt = ['%s', '%f', '%f', '%f'])   