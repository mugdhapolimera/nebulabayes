# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 23:31:37 2018

@author: mugdhapolimera
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def mz(m):
    if (m < 8.8):
        p = np.poly1d([0.53518733, 3.66817274])
        z = p(m)
    else:
        m = m-10    
        z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
    #p = np.poly1d(np.polyfit(x, z[np.where((m<9.6) & (m>9.0))[0]], 1))    
    return z

def zprior(logmstar):
    Z = []    
    for m in logmstar:
        if (m < 8.8):
            p = np.poly1d([0.53518733, 3.66817274])
            Z.append(p(m))
        else:
            m = m-10    
            Z.append(8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4))
    Z_hist = np.histogram(Z, bins = 'fd')
    logzprior = np.poly1d(np.polyfit(Z_hist[1][:-1],Z_hist[0],5))    
    plt.hist(Z, bins = 'fd')
    zlim = np.linspace(min(Z_hist[1])-0.1, max(Z_hist[1])+0.1, 100)
    plt.plot(zlim,logzprior(zlim)) 
    return logzprior

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.pkl'
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_raw.pkl'
infile = pd.read_pickle(inputfile)

logzprior = zprior(infile.logmstar)
#plt.plot(np.linspace(7.5,11,100), logzprior(np.linspace(7.5,11,100)))
#zlim = np.linspace(min(Z_hist[1])-0.1, max(Z_hist[1])+0.1, 100)
#plt.figure()    
#plt.plot(zlim,logzprior(zlim))
    