# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:07:49 2019

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity # see Ivezic+ pp. 251-255
import os 

rdat = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.csv")
#rdat = pd.read_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter.pkl")
# a simple Te method calib. is the NII/Halpha PP04, for ???2.5 < N2 < ???0.3
# Pettini & Pagel 2004

N2 = np.log10(rdat.nii_6584_flux/rdat.h_alpha_flux)

#selgd<-which(N2>(-2.5) & N2<(-0.3))
selbd = ((N2<-2.5) | (N2>-0.3))
PP04_12logOH = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
#O3N2 = np.log10((rdat.oiii_5007_flux/rdat.h_beta_flux)/(rdat.nii_6584_flux/rdat.h_alpha_flux))
#PP04_12logOH = 8.73-0.32*O3N2
rep = np.isnan(PP04_12logOH)
PP04_12logOH[rep] = (-99.)
PP04_12logOH[selbd] = (-99.)

PP04_12logOH = PP04_12logOH[PP04_12logOH > 0]
(n, bins) = np.histogram(PP04_12logOH, bins = 'fd')

bw = (bins[2]-bins[1]) 
# initially using 0.5*Knuth binsize from above as bandwidth; should test other values
kde = KernelDensity(kernel='epanechnikov',bandwidth=bw).fit(PP04_12logOH[:,np.newaxis])
xx = np.linspace(7.45897,9.061029996,244)[:,np.newaxis]
logdens = kde.score_samples(xx)
plt.figure()
plt.plot(xx,np.exp(logdens),color='green',label='kde', linewidth = 5)
plt.legend(loc="best")
plt.hist(PP04_12logOH, bins = 'fd', normed = True, histtype = 'step')

prior = np.column_stack((xx,np.exp(logdens)))
os.chdir(r'C:\Users\mugdhapolimera\github\nebulabayes')
#np.savetxt('NB_prior_O3N2.txt',prior)