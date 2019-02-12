# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:02:33 2018

@author: mugdhapolimera
"""

from __future__ import print_function, division
import os
from astropy.io import fits
from astropy.table import Table  # Used in converting to pandas DataFrame 
import numpy as np
# import pandas as pd
import NebulaBayes
from NebulaBayes import NB_Model
import time
import pandas as pd
import matplotlib.pyplot as plt
# By default save the output files in the NebulaBayes/docs subdirectory,
# assuming this file is still in that directory.
DOCS_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = DOCS_PATH+"/resolve_nicholls_bpt1filter_SEL/"
time_start = time.time()

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
    
    logzprior = np.poly1d(np.polyfit(Z_hist[1][:-1],Z_hist[0],6))    
    zlim = np.linspace(min(Z_hist[1])-0.1, max(Z_hist[1])+0.1, 100)
    #plt.figure()    
    #plt.plot(zlim,logzprior(zlim))
    #plt.hist(Z, bins = 'fd')    
    return logzprior
def mz(m):
    if (m < 8.8):
        p = np.poly1d([0.53518733, 3.66817274])
        z = p(m)
    else:
        m = m-10    
        z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
    #p = np.poly1d(np.polyfit(x, z[np.where((m<9.6) & (m>9.0))[0]], 1))    
    return z

'''def filter_HII_grid():
    """
    Filter the built-in HII grid to reduce the covered parameter space.
    Returns the grid table as a pandas DataFrame.
    """
    # First load the binary grid table and convert to a pandas DataFrame table
    NB_dir = os.path.dirname(os.path.realpath(NebulaBayes.__file__))
    grid_table_file = os.path.join(NB_dir, "grids", "NB_HII_grid.fits.gz")
    BinTableHDU_0 = fits.getdata(grid_table_file, 0)
    DF_grid = Table(BinTableHDU_0).to_pandas()

    # Fix log P/k to the value log10 P/k = 6.6, reducing the grid from
    # 3 dimensions to 2 (remaining parameters: "log U" and "12 + log O/H")
    print("Original HII grid P/k values:", np.unique(DF_grid["log P/k"]))
    DF_grid = DF_grid[DF_grid["AGNFRAC"] == 0]

    return DF_grid
'''
#############################################################################
interp_shape = [50,50]  # Number of interpolated points along each dimension
# Specify non-dereddened HII region emisson lines and fluxes  
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.pkl'
#inputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_bpt1filter.pkl"
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter.pkl'
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter.pkl'

infile = pd.read_pickle(inputfile)

#Only SELs
linelist_full = ['hbeta', 'oiii4959', 'oiii5007', 'oi6300', 'nii6548', 
                 'halpha','nii6584', 'sii6717', 'sii6731']

#All Lines for Levesque
#linelist_full = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
#        'hbeta', 'oiii4959', 'oiii5007', 'hei5875', 'oi6300', 'nii6548', 
#        'halpha', 'nii6584', 'sii6717', 'sii6731', 'ariii7136']

#All Lines for BPASS/STB99
#linelist_full = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
#        'heii4685', 'ariv4711', 'hbeta', 'oiii4959', 'oiii5007', 'hei5875', 
#   'oi6300', 'nii6548', 'halpha', 'nii6584', 'sii6717', 'sii6731', 'ariii7136']


fluxnames = {'oii3726' : 'Flux_OII_3726',
'oii3729' : 'Flux_OII_3728',
'neiii3869' : 'neiii_3869_flux',
'oiii4363' : 'oiii_4363_flux',
'hgamma' : 'h_gamma_flux',
'heii4685' : 'Flux_HeII_4685',
'ariv4711' : 'Flux_ArIV_4711',
'hbeta' : 'h_beta_flux',
'oiii4959' : 'oiii_4959_flux',
'oiii5007' : 'oiii_5007_flux',
'hei5875' : 'hei_5876_flux',
'oi6300' : 'oi_6300_flux',
'nii6548' : 'nii_6548_flux',
'halpha' : 'h_alpha_flux',
'nii6584' : 'nii_6584_flux',
'sii6717' : 'sii_6717_flux',
'sii6731' : 'sii_6731_flux',
'ariii7136' : 'ariii7135_flux'}

errornames = {'oii3726' : 'Flux_OII_3726_Err',
'oii3729' : 'Flux_OII_3728_Err',
'neiii3869' : 'neiii_3869_flux_err',
'oiii4363' : 'oiii_4363_flux_err',
'hgamma' : 'h_gamma_flux_err',
'heii4685' : 'Flux_HeII_4685_Err',
'ariv4711' : 'Flux_ArIV_4711_Err',
'hbeta' : 'h_beta_flux_err',
'oiii4959' : 'oiii_4959_flux_err',
'oiii5007' : 'oiii_5007_flux_err',
'hei5875' : 'hei_5876_flux_err',
'oi6300' : 'oi_6300_flux_err',
'nii6548' : 'nii_6548_flux_err',
'halpha' : 'h_alpha_flux_err',
'nii6584' : 'nii_6584_flux_err',
'sii6717' : 'sii_6717_flux_err',
'sii6731' : 'sii_6731_flux_err',
'ariii7136' : 'ariii7135_flux_err'}


# Initialise the NB_Model, which loads and interpolates the model flux grids:
#gridfile = r'C:\Users\mugdhapolimera\github\izi\Richardson_bpass_binary_csf_n1e2_40Myr.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_binary_csf_n1e2_40.0Myr_lines.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_binary_csf_n1e2_40.0Myr_newgrid.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr_newgrid.fits'
#gridfile = 'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/l09_high_csf_n1e2_6.0Myr_new.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr-GrovesDep_new.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr_new2.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_BPASSbinary_csf_n1e2_40.0Myr-GrovesDep.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_BPASSbinary_csf_n1e2_40.0Myr-GrovesDepNoN.fits'
gridfile = r'C:/Users/mugdhapolimera/github/izi/Richardson-0-0_1-0agn-BPASS-Binary-CSF-n=1e2-40.0Myr-NichollsCE.fits'
#grid = pd.read_csv(gridfile)
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_bpass_secular_csf_n1e2_10_0Myr.fits'
grid0 = Table.read(gridfile, format='fits')
grid = grid0.to_pandas()
grid['LOGZ'] += 8.76
grid = grid[grid["AGNFRAC"] == 0]
grid = grid[(grid["LOGQ"] > 6.9) & (grid["LOGQ"] < 8.9)]

logzprior = zprior(infile.logmstar)

def calculate_custom_prior(infile_obs, obs_flux_arr_dict, obs_err_arr_dict,
                           grids_dict, grid_spec, grid_rel_err, logzprior=logzprior):
    """
    Example callback function passed to NebulaBayes to calculate a custom prior.

    Calculate a prior Pinfile over the n-D interpolated grid, combining a
    contribution from a line-ratio prior with a contribution from information
    on a particular grid parameter (the ionisation parameter).

    This function shows the required inputs and output.

    Arguments for the observed data:
        infile_obs:     A pandas DataFrame table holding the observed fluxes, with
                    a row for each line.
        obs_flux_arr_dict: A dictionary mapping line names to n-D arrays of
                    observed fluxes over the entire grid.  The fluxes will be
                    the same everywhere unless deredden=True, in which case the
                    fluxes were dereddened to match the predicted Balmer
                    decrement at each point in the interpolated parameter space.
                    The fluxes have been normalised to the input norm_line.
        obs_err_arr_dict:  Same as obs_flux_arr_dict, but for the flux errors.
    Arguments for the model data:
        grids_dict: Dictionary that maps line names to n-D interpolated flux
                    arrays.  These interpolated arrays are based on the input
                    grid and have been normalised to the input norm_line.
        grid_spec:  A NB1_Process_grids.Grid_description instance holding basic
                    information for the interpolated grids, such as the
                    parameter names, parameter values and the grid shape.
        grid_rel_err: The systematic relative error on grid fluxes, as a linear
                    proportion between 0 and 1.  This is an input into
                    NebulaBayes.

    Return a numpy array of the value of the prior over the full interpolated
    parameter space.  The prior is in linear probability space (not log).
    """
    # Next calculate a contribution to the prior that varies only with the
    # ionisation parameter, log U.
    # Find the index of the "log U" grid parameter (each parameter corresponds
    # to a grid dimension; this is the index of that grid dimension):
    param_ind = grid_spec.param_names.index("LOGZ")
    # The list of interpolated parameter values for the "log U" grid parameter:        
    allZ_values = grid_spec.param_values_arrs[param_ind]  # Sorted 1D array
    #zlim = np.linspace(min(allZ_values), max(allZ_values), grid_spec.shape[0])
    #logZ_1D_prior = norm.pinfile(zlim, loc = mz(infile.logmstar[gal]), scale = 10**0.08)
    #print(zlim, logZ_1D_prior)    
    #logZ_1D_prior     = logzprior(allZ_values)
    
    #logZ_1D_prior[np.where(allZ_values < 7.95)[0]] = logzprior(7.95)
    #logZ_1D_prior[np.where(allZ_values > 8.97)[0]] = logzprior(8.97)
    
    #plt.figure()    
    #plt.plot(allZ_values, logZ_1D_prior)
    # This array has only one dimension.  Construct a slice to use numpy
    # "broadcasting" to apply the 1D prior over the whole 2D grid:
    # (The following method works in nD, although here it's a 2D grid)
    logZ_1D_prior = np.loadtxt(r'C:\Users\mugdhapolimera\github\nebulabayes\NB_prior_PP04.txt')[:,1]
    slice_NLR = [np.newaxis for _ in grid_spec.shape]
    slice_NLR[param_ind] = slice(None)  # "[slice(None)]" means "[:]"
    slice_NLR = tuple(slice_NLR)
    contribution_U = np.ones(grid_spec.shape) * logZ_1D_prior[slice_NLR]
    # Combine the two contributions to the prior, weighted equally
    prior = contribution_U
    # The prior Pinfile will be properly normalised by NebulaBayes later
    return prior

infile = infile.rename(index=str, columns={"name": "NAME"})
# Run parameter estimation 
for gal in range(len(infile['NAME'])):
    t1 = time.time()    
    obs_fluxes = []
    obs_errs = []
    linelist = []
    wavelengths = []
    print (gal, infile['NAME'].iloc[gal])
    for i in range(0,len(linelist_full)):
        #if (infile[fluxnames[i]][gal] > 0) & (infile[errornames[i]][gal] > 0) & (infile[fluxnames[i]][gal]/infile[errornames[i]][gal] > 2)  : 
            #if fluxnames[i] == 'sii_6717_flux':
            #    obs_fluxes.append((infile[fluxnames[i]][gal]+infile['sii_6731_flux'][gal])/2)
            #    obs_errs.append((infile[errornames[i]][gal]**2 + infile['sii_6731_flux_err'][gal]**2)**0.5)
            #else:   
            obs_fluxes.append(infile[fluxnames[linelist_full[i]]][gal])
            obs_errs.append(infile[errornames[linelist_full[i]]][gal])
            linelist.append(linelist_full[i])
            
            #wavelengths.append(wavelength_list[i])
    #normline = linelist[np.where(obs_fluxes == max(obs_fluxes))[0][0]]

    maxflux = linelist[np.where(obs_fluxes == max(obs_fluxes))[0][0]]
    #NB_Model_HII = NB_Model(griinfileile, grid_params = ["LOGZ", "LOGQ"],line_list=linelist)
#    NB_Model_HII = NB_Model(gridfile, grid_params = ["LOGZ", "LOGQ"], 
#                    line_list = linelist, interpd_grid_shape=interp_shape, 
#                    grid_error=0.5)
    NB_Model_HII = NB_Model(grid, grid_params = ["LOGZ", "LOGQ"], #, "AGNFRAC"], 
                 line_list = linelist, grid_error=0.5)
    #interpd_grid_shape=interp_shape,
    kwargs = {#"prior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"_prior_plot.pdf"),
          #"likelihood_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"likelihood_plot.pdf"),
          #"posterior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"posterior_plot.pdf"),
          "estimate_table": os.path.join(OUT_DIR, "RESOLVE_param_estimates.csv"),
          "best_model_table": os.path.join(OUT_DIR, "RESOLVE_best_model.csv"),
          #"deredden": True,  # Match model Balmer decrement everywhere in grid
          #"obs_wavelengths": wavelengths,  # Needed for dereddening
          "norm_line": 'hbeta',  # Obs and model fluxes normalised to Hbeta
          #"prior": calculate_custom_prior,  # The callback function
          "gal_name": infile['NAME'][gal],
          "plot_configs": [{}, {}, {}, {}]} # Prior, like, posterior, per-line
         
    # Add "Best model" table as text on the prior plot:
    kwargs["plot_configs"][0]["table_on_plot"] = True
    
    # Do parameter estimation for just the one set of fluxes:
    Result_HII = NB_Model_HII(obs_fluxes, obs_errs, linelist, **kwargs)
    t2 = time.time()
    print ("Time for galaxy {} is {}".format(infile['NAME'][gal], t2-t1))

time_end = time.time()
print ("Total time for {} galaxies is {}".format(len(infile['NAME']), time_end - time_start))

