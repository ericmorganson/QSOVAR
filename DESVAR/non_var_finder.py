#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
#import seaborn as sns
import pandas as pd
import time
import sys
import astropy.io.fits as pyfit

fits_dir = "/home/thrush2/caps_dir/"
csvpath = "/home/thrush2/QSOVAR/DESVAR/"
data_frame = pd.DataFrame()

def get_vals(args, ROW):
    flux_len = {}
    lc_median = {}
    for color in "griz":
        FITS = pyfit.open(args)
        # Obtain flux, err, and time arrays
        try:
            flux = FITS[1].data['LC_FLUX_PSF_'+color][ROW]
            flux_len[color] = len(flux[(flux!=0)])
            coadd_ID = FITS[1].data['COADD_OBJECT_ID'][ROW]
            lc_median[color] = float(FITS[1].data['MEDIAN_PSF_'+color][ROW])
        except IndexError:
            print('Error')
            exit()
        flux_len = {x:y for x,y in flux_len.items() if y!=0}
        lc_median = {x:y for x,y in lc_median.items() if y!=0}
    return coadd_ID, flux_len, lc_median


def get_a_b_chi2(fit, color, a_b_vals):
    #print(a_b_vals)
    field_name = fit.split("/")[-1].split("_")[0]
    for item in a_b_vals:
        if item[0] == field_name and item[1] == color.upper():
            a, b = item[2], item[3]
            #print(item)
    return a, b

if __name__ == "__main__":
  Fits_list = []
  Fits_list.append(str(sys.argv[1])) #["C3_lc.fits"]#, "C2_lc.fits", "C3_lc.fits",
               #"E1_lc.fits", "E2_lc.fits",
               #"S1_lc.fits", "S2_lc.fits",
               #"X1_lc.fits", "X2_lc.fits", "X3_lc.fits"]
  tmp_new = []
  print(Fits_list) 
  for Fits_file in Fits_list:
      Fits_f = fits_dir + Fits_file
      print(Fits_file)
      with fits.open(Fits_f) as hdul:
          max_rows = hdul[1].data.shape[0]
      print(max_rows)
      for ROW in range(0, int(max_rows)):
          if ROW % 100 ==0:
             if ROW ==0: 
                tick = time.time()
             else:
                tock = time.time()
                print(tock-tick)
                tick=tock 
             print(ROW)
              
          coadd_ID, flux_len, mu = get_vals(Fits_f, ROW)
          mu_list = [22.5-2.5*np.log10(mu[i]) for i in mu]
          for color in 'grizA':
              if color == "A":
                  avg_mu = sum(mu_list)/len(mu_list)
                  flux_list = [int(flux_len[i]) for i in flux_len]
                  sum_flux = sum(flux_list)
                  tmp_new.append({"band": ['All'],
                             "fits": [Fits_file[0:2]],
                             "object": [ROW],
                             "COADD_OBJECT_ID": [coadd_ID],
                             "Num_Obs": [sum_flux],
                             "Mu_Bright": [avg_mu]})
              elif color in mu.keys():                  
                  mu_col = 22.5-2.5*np.log10(mu[color])
                  
                  tmp_new.append({"band": [color],
                             "fits": [Fits_file[0:2]],
                             "object": [ROW],
                             "COADD_OBJECT_ID": [coadd_ID],
                             "Num_Obs": [flux_len[color]], 
                             "Mu_Bright": [mu_col]})
              else: 
                  continue
          if ROW % 5000 ==0:
              data_frame = data_frame.append(tmp_new, ignore_index=True)
              tmp_new = []
          elif ROW == int(max_rows)-1:
              data_frame = data_frame.append(tmp_new, ignore_index=True)
              tmp_new = []
  data_frame.to_csv(csvpath+"all_obj_mu_band_objId_obs_{}.csv".format(str(sys.argv[2])))   
  exit()
