from __future__ import print_function, division
import os
# import numpy as np
from astropy.table import Table
import pandas as pd
from NebulaBayes import NB_Model
import time
import numpy as np


# By default save the output files in the NebulaBayes/docs subdirectory,
# assuming this file is still in that directory.
DOCS_PATH = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = DOCS_PATH+"/results_izi_prior/"
time_start = time.time()

##############################################################################
print("\nRunning basic HII example...")
# In this example the NebulaBayes built-in HII region grid is used to constrain
# the three parameter of the grid - oxygen abundance, ionisation parameter
# and pressure.

# These HII-region optical emission-line fluxes have already been dereddened
#linelist = ["OII3726_29", "Hbeta", "OIII5007", "OI6300", "Halpha", "NII6583",
#            "SII6716", "SII6731"]
#obs_fluxes = [8.151, 4.634, 1.999, 0.09562, 13.21, 5.116, 1.377, 1.249]
#obs_errs = [0.09008, 0.04013, 0.01888, 0.00222, 0.07635, 0.03159, 0.00999,
#            0.00923]


#linelist_full = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
#        'hbeta', 'heii4685', 'ariv4711', 'oiii4959', 'oiii5007', 'hei5875', 
#        'oi6300', 'nii6548', 'halpha', 'nii6584', 'sii6717', 'sii6731', 'ariii7136']
#fluxnames = ['oii_3726_flux_ext', 'oii_3729_flux_ext', 
#             'neiii_3869_flux_ext','oiii_4363_flux_ext', 'h_gamma_flux_ext', 
#             'h_beta_flux_ext', 'heii_4685_flux_port_ext', 'ariv_4711_flux_port_ext',
#             'oiii_4959_flux_ext', 'oiii_5007_flux_ext', 'hei_5876_flux_ext', 
#             'oi_6300_flux_ext', 'nii_6548_flux_ext', 'h_alpha_flux_ext', 
#             'nii_6584_flux_ext', 'sii_6717_flux_ext', 'sii_6731_flux_ext', 
#             'ariii_7135_flux_ext'] 
#
#errornames = ['oii_3726_flux_ext_err', 'oii_3729_flux_ext_err', 
#              'neiii_3869_flux_ext_err',  
#              'oiii_4363_flux_ext_err', 'h_gamma_flux_ext_err',
#              'h_beta_flux_ext_err', 'heii_4685_flux_port_ext_err', 
#              'ariv_4711_flux_port_ext_err', 'oiii_4959_flux_ext_err', 
#              'oiii_5007_flux_ext_err', 'hei_5876_flux_ext_err', 
#              'oi_6300_flux_ext_err', 'nii_6548_flux_ext_err', 
#              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 
#              'sii_6717_flux_ext_err', 'sii_6731_flux_ext_err', 
#              'ariii_7135_flux_ext_err']

# Excluded 'heii_3203_flux_port_ext'- only 4 galaxies with non-zero value
# Excluded 'heii4685', 'ariv4711',- not in IZI default grid

#os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes')
#inputfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\RESOLVE_SDSS_dext.fits'
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_all_dext.fits'
#dat = Table.read(inputfile, format='fits')
#df = dat.to_pandas()
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_filtered.pkl'
infile = pd.read_pickle(inputfile)

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
#df = df[~np.isnan(df.h_alpha_flux_ext)]
#df = df[~np.isnan(df.oiii_5007_flux_ext)]
#df = df[~np.isnan(df.nii_6584_flux_ext)]
#df = df[~np.isnan(df.nii_6548_flux_ext)]
#df = df[~np.isnan(df.h_beta_flux_ext)]
#df = df[~np.isnan(df.oi_6300_flux_ext)]
#df = df[~np.isnan(df.sii_6717_flux_ext)]
#df = df[~np.isnan(df.sii_6731_flux_ext)]
#
#nii_sum = (df['nii_6584_flux_ext'] + df['nii_6548_flux_ext'])*3./4
#oiii = df['oiii_5007_flux_ext']
#h_alpha = df['h_alpha_flux_ext']
#h_beta = df['h_beta_flux_ext']
#h_beta_err = df['h_beta_flux_ext_err']
#oi = df['oi_6300_flux_ext']
#
#sii_sum = df['sii_6717_flux_ext'] + df['sii_6731_flux_ext']
#heii = df['heii_4685_flux_port_ext']
#heii_err = df['heii_4685_flux_port_ext_err']
#
#gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) & 
#(sii_sum > 0) & (h_beta > 0) & (h_beta > 3*h_beta_err))# & (heii_err > 0) & (heii > 3*heii_err))
#infile = df[gooddata]
#infile.index = range(len(infile))         
'''linelist_full = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
        'hbeta', 'oiii4959', 'oiii5007', 'hei5875', 'oi6300', 'nii6548', 
        'halpha', 'nii6584', 'sii6717', 'sii6731', 'ariii7136']
#linelist_SEL = ['hbeta', 'oiii5007', 'oi6300', 'halpha', 'nii6584', 'sii6717', 'sii6731']
linelist_SEL = ['hbeta', 'nii6584']

fluxnames = ['oii_3726_flux_ext',
 'oii_3729_flux_ext',
 'neiii_3869_flux_ext',
 'oiii_4363_flux_ext',
 'h_gamma_flux_ext',
 'h_beta_flux_ext',
 'oiii_4959_flux_ext',
 'oiii_5007_flux_ext',
 'hei_5876_flux_ext',
 'oi_6300_flux_ext',
 'nii_6548_flux_ext',
 'h_alpha_flux_ext',
 'nii_6584_flux_ext',
 'sii_6717_flux_ext',
 'sii_6731_flux_ext',
 'ariii_7135_flux_ext']
 
 errornames = ['oii_3726_flux_ext_err',
 'oii_3729_flux_ext_err',
 'neiii_3869_flux_ext_err',
 'oiii_4363_flux_ext_err',
 'h_gamma_flux_ext_err',
 'h_beta_flux_ext_err',
 'oiii_4959_flux_ext_err',
 'oiii_5007_flux_ext_err',
 'hei_5876_flux_ext_err',
 'oi_6300_flux_ext_err',
 'nii_6548_flux_ext_err',
 'h_alpha_flux_ext_err',
 'nii_6584_flux_ext_err',
 'sii_6717_flux_ext_err',
 'sii_6731_flux_ext_err',
 'ariii_7135_flux_ext_err']
 

fluxnames_SEL = [ 'h_beta_flux_ext',
 'oiii_5007_flux_ext',
 'oi_6300_flux_ext',
 'h_alpha_flux_ext',
 'nii_6584_flux_ext',
 'sii_6717_flux_ext',
 'sii_6731_flux_ext']
 
errornames_SEL = [ 'h_beta_flux_ext_err',
 'oiii_5007_flux_ext_err',
 'oi_6300_flux_ext_err',
 'h_alpha_flux_ext_err',
 'nii_6584_flux_ext_err',
 'sii_6717_flux_ext_err',
 'sii_6731_flux_ext_err']
'''
#fluxnames_SEL = [ 'h_beta_flux_ext',
# 'nii_6584_flux_ext']
# 
#errornames_SEL = [ 'h_beta_flux_ext_err',
# 'nii_6584_flux_ext_err']

fluxnames = ['oiii_4363_flux_ext', 
             'oiii_4959_flux_ext','h_beta_flux_ext', 'oiii_5007_flux_ext',
             'nii_6548_flux_ext', 'h_alpha_flux_ext', 'nii_6584_flux_ext', 
             'sii_6717_flux_ext'] 

errornames = ['oiii_4363_flux_ext_err', 
              'oiii_4959_flux_ext_err', 'h_beta_flux_ext_err',
              'oiii_5007_flux_ext_err', 'nii_6548_flux_ext_err', 
              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 
              'sii_6717_flux_ext_err']

#Defining IDs for different Flux Lines
linelist_full = ['oiii4363' , 'oiii4959', 'hbeta','oiii5007',  'nii6584', 
                 'halpha', 'nii6548', 'sii6731','sii6731' ]

 
# Fluxes/errors will be normalised to the flux of the default norm_line, Hbeta
# The emission line names match those in the grid (see the NebulaBayes/grids
# directory)

# Initialise the NB_Model, which loads and interpolates the model flux grids:
gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\l09_high_csf_n1e2_6.0Myr_new.fits'
#gridfile = r'C:\Users\mugdhapolimera\github\izi\Richardson_bpass_binary_csf_n1e2_40Myr.fits'
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_binary_csf_n1e2_40.0Myr_lines.fits'
grid0 = Table.read(gridfile, format='fits')
#grid = grid0.to_pandas()
    
#NB_Model_HII = NB_Model("HII",line_list=linelist)

# Set outputs:

# Run parameter estimation 
for gal in range(len(infile['NAME'])):
    t1 = time.time()    
    obs_fluxes = []
    obs_errs = []
    linelist = []
    print (gal, infile['NAME'][gal])
    for i in range(0,len(fluxnames)):
        if (infile[fluxnames[i]][gal] > 0) & (infile[errornames[i]][gal] > 0) & (infile[fluxnames[i]][gal]/infile[errornames[i]][gal] > 2)  : 
            if fluxnames[i] == 'sii_6717_flux_ext':
                obs_fluxes.append((infile[fluxnames[i]][gal]+infile['sii_6731_flux_ext'][gal])/2)
                obs_errs.append((infile[errornames[i]][gal]**2 + infile['sii_6731_flux_ext_err'][gal]**2)**0.5)
            else:   
                obs_fluxes.append(infile[fluxnames[i]][gal])
                obs_errs.append(infile[errornames[i]][gal])
            linelist.append(linelist_full[i])
        
    
    NB_Model_HII = NB_Model(gridfile, grid_params = ["LOGZ", "LOGQ"],line_list=linelist)

    kwargs = {"prior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"_prior_plot.pdf"),
          "likelihood_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"likelihood_plot.pdf"),
          "posterior_plot": os.path.join(OUT_DIR, infile['NAME'][gal]+"posterior_plot.pdf"),
          "estimate_table": os.path.join(OUT_DIR, "RESOLVE_param_estimates.csv"),
          "best_model_table": os.path.join(OUT_DIR, "RESOLVE_best_model.csv"),
          "gal_name": infile['NAME'][gal]}#,
          #"interpd_grid_shape" : [35,35]}
    
    Result_HII = NB_Model_HII(obs_fluxes, obs_errs, linelist, **kwargs)
    t2 = time.time()
    print ("Time for galaxy {} is {}".format(infile['NAME'][gal], t2-t1))
# NB_Model_HII may be called repeatedly to do Bayesian parameter estimation
# on different sets of observed fluxes with the same grid.

time_end = time.time()
print ("Total time for {} galaxies is {}".format(len(infile['NAME']), time_end - time_start))
'''
##############################################################################
print("\nRunning basic NLR example...")
# In this example we use the NebulaBayes built-in AGN narrow-line region (NLR)
# grid to constrain the four parameter of the grid - oxygen abundance,
# ionisation parameter, pressure, and the hardness of the ionising continuum.

# These NLR optical emission-line fluxes have already been dereddened
linelist = ["OII3726_29", "NeIII3869", "Hgamma", "HeII4686", "Hbeta",
            "OIII5007", "HeI5876", "Halpha", "NII6583", "SII6716", "SII6731"]
obs_fluxes = [4.2162, 1.159, 0.7161, 0.3970, 1.292, 12.88, 0.1597, 3.747,
              5.027, 1.105, 1.198]
obs_errs = [0.5330, 0.2073, 0.1172, 0.0630, 0.1864, 1.759, 0.0225, 0.3919,
            0.5226,  0.1156, 0.1248]

NB_Model_NLR = NB_Model("NLR", line_list=linelist)

kwargs = {"prior_plot": os.path.join(OUT_DIR, "1_NLR_prior_plot.pdf"),
          "posterior_plot": os.path.join(OUT_DIR, "1_NLR_posterior_plot.pdf"),
          "estimate_table": os.path.join(OUT_DIR, "1_NLR_param_estimates.csv"),
          "prior": [("SII6716", "SII6731")],
          }

Result_NLR = NB_Model_NLR(obs_fluxes, obs_errs, linelist, **kwargs)


##############################################################################
print("Basic example script complete.")

'''