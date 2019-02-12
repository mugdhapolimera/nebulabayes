# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 12:09:41 2018

@author: mugdhapolimera
"""

from __future__ import print_function, division
import os
from astropy.table import Table
import pandas as pd
from NebulaBayes import NB_Model
import time
import pickle
import numpy as np

# By default save the output files in the NebulaBayes/docs subdirectory,
# assuming this file is still in that directory.
DOCS_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = DOCS_PATH+"/results_bpass_agn_snr_old/"

time_start = time.time()

##############################################################################
linelist_full = ['oii3726', 'oii3729', 
                 'neiii3869', 'oiii4363', 'hgamma', 
                 'hbeta', 'heii4685', 'ariv4711',
                 'oiii4959', 'oiii5007', 'hei5875', 
                 'oi6300', 'nii6548', 'halpha', 
                 'nii6584', 'sii6717', 'sii6731', 
                 'ariii7136']
fluxnames = ['oii_3726_flux_ext', 'oii_3729_flux_ext', 
             'neiii_3869_flux_ext','oiii_4363_flux_ext', 'h_gamma_flux_ext', 
             'h_beta_flux_ext', 'heii_4685_flux_port_ext', 'ariv_4711_flux_port_ext',
             'oiii_4959_flux_ext', 'oiii_5007_flux_ext', 'hei_5876_flux_ext', 
             'oi_6300_flux_ext', 'nii_6548_flux_ext', 'h_alpha_flux_ext', 
             'nii_6584_flux_ext', 'sii_6717_flux_ext', 'sii_6731_flux_ext', 
             'ariii_7135_flux_ext'] 

errornames = ['oii_3726_flux_port_ext_err', 'oii_3729_flux_port_ext_err', 
              'neiii_3869_flux_ext_err', 'oiii_4363_flux_ext_err', 'h_gamma_flux_ext_err',
              'h_beta_flux_ext_err', 'heii_4685_flux_port_ext_err', 'ariv_4711_flux_port_ext_err',
              'oiii_4959_flux_ext_err', 'oiii_5007_flux_ext_err', 'hei_5876_flux_ext_err', 
              'oi_6300_flux_ext_err', 'nii_6548_flux_ext_err', 'h_alpha_flux_ext_err',
              'nii_6584_flux_ext_err', 'sii_6717_flux_ext_err', 'sii_6731_flux_ext_err', 
              'ariii_7135_flux_ext_err']
  
#os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes')
inputfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
df = dat.to_pandas()

#df = pd.read_pickle('C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl')
#ra=df.radeg
#dec=df.dedeg
#flinsample = df.fl_insample
#grpcz = df.grpcz
#cz = df.cz
#infall = (ra > 22*15.) | (ra < 3*15.)
#inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
#mgas = df.logmgas
#mstars = df.logmstar
#mbary = 10**mgas + 10**mstars
#inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#df = df[inobssample]
df = df[~np.isnan(df.h_alpha_flux_ext)]
df = df[~np.isnan(df.oiii_5007_flux_ext)]
df = df[~np.isnan(df.nii_6584_flux_ext)]
df = df[~np.isnan(df.nii_6548_flux_ext)]
df = df[~np.isnan(df.h_beta_flux_ext)]
df = df[~np.isnan(df.oi_6300_flux_ext)]
df = df[~np.isnan(df.sii_6717_flux_ext)]
df = df[~np.isnan(df.sii_6731_flux_ext)]

nii_sum = (df['nii_6584_flux_ext'] + df['nii_6548_flux_ext'])*3./4
oiii = df['oiii_5007_flux_ext']
h_alpha = df['h_alpha_flux_ext']
h_beta = df['h_beta_flux_ext']
h_beta_err = df['h_beta_flux_ext_err']
oi = df['oi_6300_flux_ext']

sii_sum = df['sii_6717_flux_ext'] + df['sii_6731_flux_ext']
#heii = df['heii_4685_flux_port_ext']
#heii_err = df['heii_4685_flux_port_ext_err']
heii = df['Flux_HeII_4685_ext']
heii_err = df['Flux_HeII_4685__ext_Err']
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) & 
(sii_sum > 0)  & (h_beta/h_beta_err > 3.0))# & (heii/heii_err >=3) & (heii_err > 0))
infile = df[gooddata]
infile.index = range(len(infile))         
 
# Fluxes/errors will be normalised to the flux of the default norm_line, Hbeta
# The emission line names match those in the grid (see the NebulaBayes/grids
# directory)

# Initialise the NB_Model, which loads and interpolates the model flux grids:
#gridfile = r'C:\Users\mugdhapolimera\github\izi\Richardson_bpass_binary_csf_n1e2_40Myr.fits'
gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_binary_csf_n1e2_40.0Myr_lines.fits'
grid0 = Table.read(gridfile, format='fits')
    
# Set outputs:

# Run parameter estimation 
for gal in range(759,len(infile['NAME'])): #679, 759
    t1 = time.time()    
    obs_fluxes = []
    obs_errs = []
    linelist = []
    print (gal, infile['NAME'][gal])
    for i in range(0,len(fluxnames)):
        if (infile[fluxnames[i]][gal] > 0) & (infile[errornames[i]][gal] > 0) & (infile[fluxnames[i]][gal]/infile[errornames[i]][gal] > 2)  : 
            obs_fluxes.append(infile[fluxnames[i]][gal])
            obs_errs.append(infile[errornames[i]][gal])
            linelist.append(linelist_full[i])
    NB_Model_HII = NB_Model(gridfile, grid_params = ["AGNFRAC", "LOGZ", "LOGQ"],line_list=linelist)
    #if not os.path.exists(OUT_DIR+'lineplots/'+infile['NAME'][gal]):
    #    os.mkdir(OUT_DIR+'lineplots/'+infile['NAME'][gal])
    kwargs = {"prior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"_prior_plot.pdf"),
          "likelihood_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"likelihood_plot.pdf"),
          "posterior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"posterior_plot.pdf"),
          "estimate_table": os.path.join(OUT_DIR, "RESOLVE_param_estimates.csv"),
          "best_model_table": os.path.join(OUT_DIR, "RESOLVE_best_model.csv"),
          #"line_plot_dir": os.path.join(OUT_DIR+'lineplots'+infile['NAME'][gal]),
          "gal_name": infile['NAME'][gal]}
    
    Result_HII = NB_Model_HII(obs_fluxes, obs_errs, linelist, **kwargs)
    with open(OUT_DIR+'result/'+infile['NAME'][gal]+'.pkl', 'wb') as output:  # Overwrites any existing file.
        pickle.dump(Result_HII, output, pickle.HIGHEST_PROTOCOL)
    t2 = time.time()
    print ("Time for galaxy {} is {}".format(infile['NAME'][gal], t2-t1))
# NB_Model_HII may be called repeatedly to do Bayesian parameter estimation
# on different sets of observed fluxes with the same grid.

time_end = time.time()
print ("Total time for {} galaxies is {}".format(len(infile['NAME']), time_end - time_start))
