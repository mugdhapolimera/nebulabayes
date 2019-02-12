# -*- coding: utf-8 -*-
"""
Created on Wed Oct 03 21:52:30 2018

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
#import matplotlib.mlab as mlab
import os
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
os.chdir('C:/Users/mugdhapolimera/github/nebulabayes/')

results = pd.read_csv("results_bpass_agn_full_interp_SEL/RESOLVE_param_estimates.csv")
#inputfile = 'C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl'
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter.pkl'

infile = pd.read_pickle(inputfile)
#inputfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\RESOLVE_SDSS_dext.fits'
#dat = Table.read(inputfile, format='fits')
#infile = dat.to_pandas()

Z_index = (results['Parameter'] == 'LOGZ')
results_Z = results[Z_index]
results_Z.index = results_Z['Galaxy Name'] #range(len(results_Z))
results_Z["Estimate"] = results_Z["Estimate"]

agn_index = (results['Parameter'] == 'AGNFRAC')
results_agn = results[agn_index]
#results_agn.index = range(len(results_agn))
results_agn.index = results_agn['Galaxy Name']

#results_Z['CI68_low'] = results_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
#results_Z['CI68_high'] = results_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)

Z_results = results_Z['Estimate']
#rrup_results = results_Z['CI68_high'] - results_Z['Estimate']
#errdown_results = results_Z['Estimate'] - results_Z['CI68_low']
#neg_errs = (errup_results < 10**-6) & (errup_results != 0)

#neg_errs_ndx = [infile['NAME'][x] in list(results_Z['Galaxy Name'][neg_errs]) for x in range(len(infile))]
pd.set_option('display.max_columns', 500)
#outliers = (results_Z['Estimate'] == min(results_Z['Estimate'])) | (results_Z['Estimate'] == max(results_Z['Estimate']))
#results_Z = results_Z[~outliers]
#results_agn = results_agn[~outliers]
#infile = infile[~outliers]

bpt = pd.read_csv(r'C:\Users\mugdhapolimera\github\BPT\resolve_emlineclass_full.csv')
bpt.index = bpt.galname #range(len(bpt))
#bpt = bpt[~outliers]
ambig = bpt.agntosf | bpt.ambigagn | bpt.composite | bpt.sftoagn

ambig_ndx = [results['Galaxy Name'][x] in list(bpt['galname'][ambig]) for x in range(len(results))]
mstars = infile['logmstar']

if 'heiisel' in bpt.keys():
    ambig = ambig | bpt.heiisel    
    heii_ndx = bpt.index #[results_Z['Galaxy Name'][x] in list(bpt['galname']) for x in range(len(results_Z))]
    results_Z = results_Z.loc[heii_ndx]
    #results_Z.index = range(len(results_Z))

    results_agn = results_agn.loc[heii_ndx]
    mstars = infile.loc[heii_ndx]['logmstar']
#results_agn.index = range(len(results_agn))

    #odd = bpt.heiisel & (results_agn['Estimate'] == 0)
    #print results_agn[odd]
    #odd_ndx = [x for x in range(len(infile)) if infile['NAME'][x] in list(bpt['galname'][odd]) ]
    #print infile.iloc[odd_ndx]['heii_4685_flux_port_ext_err']
    #print infile.iloc[odd_ndx]['heii_4685_flux_port_ext']


#print results_Z[ambig_ndx]['Estimate'] + 8.76
#print results_agn[ambig_ndx]['Estimate']
'''for i in range(len(ambig_ndx)):
    if ambig_ndx[i]:
        if i%3 == 0:
            print results.iloc[i]['Galaxy Name']
        print results.iloc[i]['Parameter'], results.iloc[i]['Estimate']
'''
plt.figure()
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'g', linewidth = 5,label = 'Tremonti+04')
m = np.linspace(8.5 - 10, max(mstars) - 10, 100)
Z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,Z, color = 'brown', linewidth = 5, label = 'Manucci+10')
plt.ylabel('Z Estimate (12 + Log[O/H])')
plt.xlabel('Stellar Mass')
#for i in np.arange(7.5,11.5 , 0.2):
masses = np.arange(min(mstars), max(mstars), 0.1)
sel = []#(mstars < 0.0)
mass = []
met = []
for i in range(len(masses) - 1):
     
    sel.extend(list(np.where((mstars >= masses[i]) & (mstars < masses[i+1]))[0]))
    #print masses[i], masses[i+1], len(np.where(sel)[0]) #len(np.where((mstars >= masses[i]) & (mstars < masses[i+1]))[0])
    #if len(np.where(sel)[0])> 25:

    if len(sel) > 50:
        Z = [x for x in results_Z.iloc[sel]['Estimate'] if x != 7.45897]
        #plt.figure()
        #Z_dist = np.hist(Z, bins = 'fd')
        med,sig = norm.fit(Z) #np.median(Z)
        mass.append(masses[i])
        met.append(med)
        sel = []#(mstars < 0.0)
    
    #m = mstars[sel]
    
plt.plot(mass, met, 'k', linewidth = 5, label = 'This Work')
plt.legend(loc = 4)

#mstars = infile.loc[results_Z['Galaxy Name']]['logmstar']
#mstars.index = range(len(mstars))
plt.plot(mstars[bpt.defstarform], results_Z[bpt.defstarform]['Estimate'],'k.',alpha = 0.5)
#plt.plot(mstars[bpt.defseyf | bpt.defliner], results_Z[bpt.defseyf | bpt.defliner]['Estimate'],'ro',alpha = 0.5)        
plt.plot(mstars[bpt.agntosf], results_Z[bpt.agntosf]['Estimate'],'g^', markersize = 8)
plt.plot(mstars[(bpt.defseyf | bpt.defliner | bpt.ambigagn)], 
        results_Z[(bpt.defseyf | bpt.defliner | bpt.ambigagn)]['Estimate'],
        'rs', markersize = 8)
plt.plot(mstars[bpt.composite], results_Z[bpt.composite]['Estimate'],'bs', markersize = 8)
if ('heiisel' in bpt.keys()):
    plt.plot(mstars[bpt.heiisel], results_Z[bpt.heiisel]['Estimate'],'ks',mfc ='none', mew = 2, markersize = 8)
plt.plot(mstars[bpt.sftoagn], results_Z[bpt.sftoagn]['Estimate'],'m*', markersize = 8)
        
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'g', linewidth = 5)
m = np.linspace(8.5 - 10, max(mstars) - 10, 100)
Z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
plt.plot(m+10,Z, color = 'brown', linewidth = 5)
plt.plot(mass, met, 'k', linewidth = 5)

plt.figure()
plt.plot(mstars[bpt.defstarform], results_agn[bpt.defstarform]['Estimate'],'k.',alpha = 0.5)
#plt.plot(mstars[bpt.defseyf | bpt.defliner], results_agn[bpt.defseyf | bpt.defliner]['Estimate'],'ro',alpha = 0.5)        
plt.plot(mstars[bpt.agntosf], results_agn[bpt.agntosf]['Estimate'],'g^', markersize = 8)
plt.plot(mstars[(bpt.defseyf | bpt.defliner | bpt.ambigagn)], 
         results_agn[(bpt.defseyf | bpt.defliner | bpt.ambigagn)]['Estimate'],'rs', markersize = 8)
plt.plot(mstars[bpt.composite], results_agn[bpt.composite]['Estimate'],'bs', markersize = 8)
if ('heiisel' in bpt.keys()):
    plt.plot(mstars[bpt.heiisel], results_agn[bpt.heiisel]['Estimate'],'ks',mfc ='none', mew = 2, markersize = 8)
plt.plot(mstars[bpt.sftoagn], results_agn[bpt.sftoagn]['Estimate'],'m*', markersize = 8)
plt.ylim(-0.2, 1.2)
plt.xlabel('Stellar Mass')
plt.ylabel('AGN Fraction')

plt.figure()
#plt.plot(results_Z[~ ambig]['Estimate'], results_agn[~ambig]['Estimate'],'k.',alpha = 0.5)
plt.plot(results_Z[bpt.defstarform]['Estimate'], results_agn[bpt.defstarform]['Estimate'],'k.',alpha = 0.5)
plt.plot(results_Z[bpt.agntosf]['Estimate'], results_agn[bpt.agntosf]['Estimate'],'g^', markersize = 8)
plt.plot(results_Z[(bpt.defseyf | bpt.defliner | bpt.ambigagn)]['Estimate'], 
         results_agn[(bpt.defseyf | bpt.defliner | bpt.ambigagn) ]['Estimate'],'rs', markersize = 8)
plt.plot(results_Z[bpt.composite]['Estimate'], results_agn[bpt.composite]['Estimate'],'bs', markersize = 8)
#plt.plot(results_Z[(bpt.defseyf | bpt.defliner)]['Estimate'], results_agn[(bpt.defseyf | bpt.defliner)]['Estimate'],'ro', alpha = 0.5)

if 'heiisel' in bpt.keys():
    plt.plot(results_Z[bpt.heiisel]['Estimate'], results_agn[bpt.heiisel]['Estimate'],'ks' ,mfc ='none', mew = 2, markersize = 8)
plt.plot(results_Z[bpt.sftoagn]['Estimate'], results_agn[bpt.sftoagn]['Estimate'],'m*', markersize = 8)

plt.ylim(-0.2, 1.2)
plt.xlabel('Z Estimate (12 + Log[O/H])')
plt.ylabel('AGN Fraction')
#print infile[neg_errs_ndx]
#x = np.linspace(min(Z_results)-1, max(Z_results)+1, 100)

#plt.figure()
#plt.plot(x,mlab.normpdf(x, Z_results[0], sigma))

'''
fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(results_Z[bpt.defstarform]['Estimate'], results_agn[bpt.defstarform]['Estimate'], 
        mstars[bpt.defstarform],c = 'k', marker = 'o' ,alpha = 0.5)
ax.scatter(results_Z[bpt.agntosf]['Estimate'], results_agn[bpt.agntosf]['Estimate'],
           mstars[bpt.agntosf],c = 'g', marker = '^', s = 50)
ax.scatter(results_Z[bpt.sftoagn]['Estimate'], results_agn[bpt.sftoagn]['Estimate'],
        mstars[bpt.sftoagn],c = 'm', marker = '*', s = 50)
ax.scatter(results_Z[bpt.ambigagn]['Estimate'], results_agn[bpt.ambigagn]['Estimate'],
        mstars[bpt.ambigagn],c = 'r', marker = 's', s= 50)
ax.scatter(results_Z[bpt.composite]['Estimate'], results_agn[bpt.composite]['Estimate'],
        mstars[bpt.composite],c = 'b', marker = 's', s= 50)
ax.scatter(results_Z[(bpt.defseyf | bpt.defliner)]['Estimate'], 
                  results_agn[(bpt.defseyf | bpt.defliner)]['Estimate'],
                  mstars[bpt.defseyf | bpt.defliner],c = 'r', marker = 'o', alpha = 0.5)
#ax.scatter(results_Z['Estimate'], results_agn['Estimate'], mstars)
ax.set_xlabel('Metallicity (12 + Log[O/H])')
ax.set_ylabel('AGN Fraction')
ax.set_zlabel('Mass (log)')
'''