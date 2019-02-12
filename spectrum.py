# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 12:37:41 2018

@author: mugdhapolimera
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def gauss(A,mu,sigma,color):
    x = np.linspace(mu-5*sigma, mu+5*sigma, 100)    
    plt.plot(x,A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.))), color = color)

df = pd.read_csv(r'C:\Anaconda2\Lib\site-packages\NebulaBayes\docs\results_bpass_agn_snr\RESOLVE_best_model.csv',sep=",")
linelist = {'oii3726': 3726, 'oii3729':3729, 
                 'neiii3869':3869, 'oiii4363':4363, 'hgamma':4340, 
                 'hbeta':4861, 'heii4685':4685, 'ariv4711':4711,
                 'oiii4959':4959, 'oiii5007':5007, 'hei5875':5875, 
                 'oi6300':6300, 'nii6548':6548, 'halpha':6563, 
                 'nii6584':6584, 'sii6717':6717, 'sii6731':6731, 
                 'ariii7136':7136}
#HeII selected AGN but with 0% AGN fraction
#names = ['rf0058',
# 'rf0110',
# 'rs0378',
# 'rs0521',
# 'rs0598',
# 'rs0754',
# 'rs0809',
# 'rs0945',
# 'rs1049',
# 'rs1191',
# 'rs1197',
# 'rs1210',
# 'rs1294']
names = ['rs1191']#, 'rs1197']
for name in names:
    gal = df[df['Galaxy Name'] == name]
    gal.index = range(len(gal))
    noise = gal['Obs']/gal['Obs_S/N']
    print gal
    plt.figure()   
    plt.title(gal['Galaxy Name'][0])
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    for i in range(1,len(gal)):
        gauss(gal['Obs'][i], linelist[gal['Line'][i]], noise[i], 'b')
        gauss(gal['Model'][i], linelist[gal['Line'][i]]+10, 0.1,'r')

#Artificially converging to lowest possible Z
names = ['rf0028',
 'rf0102',
 'rf0110',
 'rf0120',
 'rf0151',
 'rf0201',
 'rf0269',
 'rf0365',
 'rf0367',
 'rf0426',
 'rf0443',
 'rf0480',
 'rf0545',
 'rf0670',
 'rf0700',
 'rf0783',
 'rs0047',
 'rs0095',
 'rs0105',
 'rs0124',
 'rs0216',
 'rs0271',
 'rs0325',
 'rs0330',
 'rs0332',
 'rs0402',
 'rs0416',
 'rs0492',
 'rs0545',
 'rs0626',
 'rs0713',
 'rs0752',
 'rs0785',
 'rs0799',
 'rs0805',
 'rs0853',
 'rs0892',
 'rs0973',
 'rs0977',
 'rs0982',
 'rs1002',
 'rs1009',
 'rs1041',
 'rs1073',
 'rs1078',
 'rs1105',
 'rs1135',
 'rs1163',
 'rs1164',
 'rs1166',
 'rs1196',
 'rs1263',
 'rs1283',
 'rs1298',
 'rs1368']
        