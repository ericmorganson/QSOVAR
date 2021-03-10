#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
import seaborn as sns
import astropy.io.fits as pyfit



def get_vals(args, ROW):
    color_dict = {"g": 1, "r": 2, "i": 3, "z": 4}
    lc_median = {}
    time = np.array([])
    flux_norm = np.array([])
    err_norm = np.array([])
    array_org = np.array([])
    fig = plt.figure(figsize=(12, 6))
    ax5 = fig.add_subplot(111)
    color_plt = iter(cm.rainbow(np.linspace(0, 1, 5)))
    #a_b_list = read_a_b_chi2()
    #print(a_b_list)
    for color in "griz":
        FITS = pyfit.open(args)
        # Obtain flux, err, and time arrays
        try:
            lc_flux = FITS[1].data['LC_FLUX_PSF_'+color][ROW]
        except IndexError:
            print('Error')
            exit()

        lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
        lc_time = FITS[1].data['LC_MJD_'+color][ROW]
        limit = len(lc_time[(lc_time != 0)*(lc_flux_err <2)])
        lc_median[color] = FITS[1].data['MEDIAN_PSF_'+color][ROW] #change back to median
        if limit < 3:
            lc_median[color] = 0.0
            #print("Not enough " + color + " observations!")
            continue

        lcur_time = lc_time[(lc_time != 0)*(lc_flux_err <2)]
        lcur_flux = lc_flux[(lc_time != 0)*(lc_flux_err <2)] #:limit
        lcur_flux_err = lc_flux_err[(lc_time != 0)*(lc_flux_err <2)]  # :limit

        col = next(color_plt)

        mean_err = FITS[1].data['MEANERR_PSF_'+color][ROW]

        if (lc_flux_err>2).any():
            mean_err = np.mean(lcur_flux_err)


        a, b = get_a_b_chi2(args, color, a_b_list)
        #print(a, b)
        #print("OG ERR "+color)
        #print(lc_flux_err)

        flux_err_corr = np.sqrt(np.abs((a**2)*lcur_flux*np.sqrt(np.abs(lcur_flux_err**2 - mean_err**2)) + (b**2)*(lcur_flux_err**2 - mean_err**2)))
        normed_err = flux_err_corr/lc_median[color] #lc_flux_err/lc_median[color]
        #print("AB CORR ERR "+color)
        #print(flux_err_corr)
        ax5.scatter(lcur_time,
                    (lcur_flux - lc_median[color])/lc_median[color],
                    label=color, c=np.array([col]))
        ax5.errorbar(lcur_time,
                    (lcur_flux - lc_median[color])/lc_median[color],
                    yerr=normed_err, ecolor=np.array(col),
                    linestyle="None")

        ax5.legend()
        ax5.set_title("Pre-correction light curve: Row "+ str(ROW))
        ax5.set_ylabel("Median-corrected Flux")
        ax5.set_xlabel("time [MJD]")
        normed_flux = (lcur_flux - lc_median[color])/lc_median[color]

        time = np.append(time, lcur_time)  # remove the zeros
        flux_norm = np.append(flux_norm, normed_flux)
        err_norm = np.append(err_norm, normed_err)
        array_org = np.append(array_org, color_dict[color]*np.ones(limit))

    time, flux_norm, err_norm, array_org = map(list,
                                               zip(*sorted(zip(time,
                                                               flux_norm,
                                                               err_norm,
                                                               array_org))))
    return flux_norm, err_norm, time, lc_median, array_org, FITS

def read_a_b_chi2():
    # must have a_b_calc.txt for this to work
    # a_b_calc.txt is made with the following command:
    # python chi_squared_LCs_all_band_data.py > a_b_calc.txt
    a_b_vals = []
    with open("a_b_calc.txt","r") as infile:
        for line in infile:
            #print(line)
            line = line.split('\n')[0]

            field, col, a, b = line.strip('()').split(',')

            a_b_vals.append((field.replace("'", "").strip(), col.replace("'", "").strip(), float(a), float(b)))
    return a_b_vals


def get_a_b_chi2(fit, color, a_b_vals):
    #print(a_b_vals)
    field_name = fit.split("/")[-1].split("_")[0]
    for item in a_b_vals:
        if item[0] == field_name and item[1] == color.upper():
            a, b = item[2], item[3]
            #print(item)
    return a, b

if __name__ == "__main__":
  #VARIABILITY CHI2 & STD SCATTER PLOTS
  with fits.open('../../C3_lc.fits') as hdul:
      max_rows = hdul[1].data.shape[0]

  Fits_file = r'../../C3_lc.fits'
  var_strict = 3
  var_chi = []
  var_std = []
  non_var_chi = []
  non_var_std = []
  goodrownum = 0
  a_b_list = read_a_b_chi2()

  for ROW in range(0, int(max_rows/5)):
      if ROW%1000 == 0:
          print("Running object "+str(ROW))
      #fig = plt.figure(figsize=(10, 10))
      flux, err, time, mu, color_sort, FITS = get_vals(Fits_file, ROW)
      plt.close("all")
      #print(min(np.array(time[1:]) - np.array(time[:-1])))

      # DOESN'T MAKE SENSE TO LOOK AT ROWS WITH NO FLUX MEASUREMENTS
      if len(flux) == 0:
          print("Flux length is zero")
          continue
      if len(flux) < 750:
          #print("Too sparse :"+ str(len(flux)))
          continue

      # ONLY LET POSITIVE FLUXES AND ERRORS THROUGH
      # NOTE: FLUXES HERE ARE NORMALIZED, SO IF FLUX = 0, THEN norm(FLUX) = -1
      zip_tfe = zip(time, flux, err, color_sort)
      filter_tfe = [(t, f, e, col) for t, f, e, col in zip_tfe if f > -1 and e > 0]
      time, flux, err, color_sort = zip(*filter_tfe)
      time = np.array(time)
      flux = np.array(flux)
      err = np.array(err)
      color_sort = np.array(color_sort)
      unique, count = np.unique(color_sort, return_counts=True)
      color_sort_dict = {1: "g", 2: "r", 3: "i", 4: "z"}
      color_counts = dict(zip([color_sort_dict[u] for u in unique], count))

      #print("Counts: " + str(color_counts))

      # ONLY LOOK AT BRIGHT OBJECTS (WITHOUT OVERSATURATION)
      dim = [i for i in mu if 22.5-2.5*np.log10(mu[i])>23]
      dim_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])>23]
      bright = [i for i in mu if 22.5-2.5*np.log10(mu[i])<16]
      bright_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])<16]
      spread = [i for i in "GRIZ" if np.abs(FITS[1].data['SPREAD_MODEL_'+i][ROW]) > .003]
      if dim:
      #    print("Row is too dim in band: "+ str(dim)+str(ROW))
      #    print(dim_val)
      #    plt.close("all")
          continue
      if bright:
      #    print("Row is HELLA bright in band: "+str(bright)+str(ROW))
      #    print(bright_val)
      #    plt.close("all")
          continue
      if spread:
          #print("Too spread out in band: "+str(spread)+str(ROW))
          continue
      goodrownum +=1
      #print(str(ROW)+ " is good")
  print("Possible stars in fifth of C3:     " + str(goodrownum))
  exit()
