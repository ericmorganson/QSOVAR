
# coding: utf-8

# # Searching for "a" and "b" values with $\chi^2$ vs median error #
# 
# So, here's how this is going down:
# 1. I'm trying to find $\chi^2$ for each object's light curve .
# 2. Bin these magnitudes by the median error of the magnitudes (x value)
# 3. Take the median of each y bin and take the 85\% and 15\% of each bin to find the error of the the y value mentioned below.
# 4. Will end with x, y, and y_error; try to fit something to this (like $a^2 + b^2$).
# 
# Note: if $1/\chi^2$ vs the median error's slope is 2, then b will be 2.
# 
# Now, let's define some things to make this easier.
# 
# $\chi^2_{reduced} = \frac{\sum{\frac{flux_i - \mu}{\sigma_i^2}}}{n-1}$
# 
# n = number of flux observations for that observation
# 
# $\sigma_{median, i} = median(\sigma_i)$
# 
# x = binning $\sigma_{median, i}$
# 
# y = median $\chi^2_{reduced}$ in that bin
# 
# $\sigma_y = \frac{percentile(\chi^2_{reduced}, 85) - percentile(\chi^2_{reduced}, 15)}{2\sqrt{n}}$
# 
# Finally, I will then be plotting x vs y, with $\sigma_y$ as the error on y.

#import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import astropy.io.fits as pyfit
from scipy import stats, optimize
np.set_printoptions(precision=5)
import sys
import os


file_default = '/home/sam/Documents/Morganson_research/C1_lc.fits'

if not sys.argv[1]:
    filename = input("Give me a file, dummy.") or file_default
elif sys.argv[1] == "default":
    filename = file_default
else:
    filename = sys.argv[1]


if not sys.argv[2]:
    image_place = input("Where do you want your images in this folder?") or ''
elif sys.argv[2] == "default":
    image_place = ''
else:    
    image_place = sys.argv[2]


while os.path.exists(image_place):
    if image_place == '':
        continue
    else: 
        image_place = input("That image folder exists, try again, silly!")
        sys.exit()
if not os.path.exists(image_place) and image_place != '':
    os.makedirs(image_place)
    image_place = image_place + str('/')
print("You are using folder "+filename)
print("The images are going to "+image_place)


#initialize variable
fits = pyfit.open(filename)[1].data


clip= 10
high_m = 24
low_m = 16
chi_val = []
sigma_med_val = []
sigma_is = []
flux_arr = []
filters = ['G', 'R', 'I', 'Z']
for filt in filters:
    print("STARTING FILTER " + filt+ "\n")
    for fit in fits:
        if (fit['MAG_AUTO_G']- fit['MAG_AUTO_R']) > 0.4:
            if fit['MAG_AUTO_'+filt] < high_m:
                if np.sum(fit["LC_FLUX_PSF_"+filt] > 0) > 10:
                    if np.abs(fit['SPREAD_MODEL_'+filt]) < .003:
                        if fit['MAG_AUTO_'+filt] > low_m:
                            fluxes_arr = fit["LC_FLUX_PSF_"+filt]
                            if any(np.where(fluxes_arr > 1e10)):
                                print(fluxes_arr[np.where(fluxes_arr > 1e10)])
                            ind_errs_arr = fit['LC_FLUXERR_PSF_'+filt]
                            fluxes, ind_errs = zip(*((x, y) for x, y in zip(fluxes_arr, ind_errs_arr) if x!= 0 and y!=0 and x<1e10))
                            fluxes = np.array(fluxes)
                            ind_errs = np.array(ind_errs)                            
                            mean_flux = fit['MEAN_PSF_'+filt]
                            norm_err = (ind_errs/fluxes)
                            chi = np.sum(((fluxes - mean_flux)**2)/(ind_errs**2))/(len(fluxes)-1)
                            if chi > clip:
                                continue
                            for ind_err in ind_errs:
                                sigma_is.append(ind_err)
                            for flux in fluxes:
                                flux_arr.append(flux)
                            chi_val.append(chi)
                            sigma_med_val.append(np.median(norm_err))
                            #err.append(np.median(norm_err)/np.sqrt(len(fluxes)))
                            

    df = pd.DataFrame({'chi2r':chi_val, 'sig_m':sigma_med_val})
    #print(df)

    m, s = stats.norm.fit(chi_val)
    print("Regular values")
    print("Mean="+str(m))
    print("std="+str(s))

    m, s = stats.norm.fit(np.log10(np.array(chi_val))) # get mean and standard deviation 
    print("\nLogged values")
    print("Mean="+str(m))
    print("std="+str(s))

    plt.hist((np.array(chi_val)), bins=30)

    #plt.show()

    bins = np.linspace(df.sig_m.min(), df.sig_m.max(), 30)
    groups = df.groupby(pd.cut(df.sig_m, bins))
    #print(groups.size())

    x = np.array(groups.median().sig_m)
    #print(x)
    y = np.array(groups.median().chi2r)
    #print(y)

    #print(groups.chi2r)
    dy = np.array((groups.chi2r.quantile(.85) - groups.chi2r.quantile(.15))/(2*np.sqrt(groups.size())))
    #print(dy)

    j = x # sigma
    k = y # chi
    dk =dy # error
    jkVar = zip(j, k, dk)

    jkFiltered = [(j, k, dk) for j, k, dk in jkVar if j <= 0.12 and j >= 0.02]

    j, k, dk = zip(*jkFiltered)
    j = np.array(j)
    k = np.array(k)
    dk = np.array(dk)

    plt.figure()
    plt.errorbar(j, k, yerr=dk, label="clipped")

    def test_func(sigma, a, b):
        return a**2/(sigma**1) + b**2

    params, params_covariance = optimize.curve_fit(test_func, j, k, #p0=[1, 3], 
                                                   sigma=dk)

    print(params)

    # calculate new x's and y's
    x_new = np.linspace(j[0], j[-1], 50)
    plt.plot(x_new, test_func(x_new, params[0], params[1]), label='a='+str(round(params[0], 4))+'\n'+'b='+str(round(params[1], 4)))
    plt.legend(loc='best')
    #plt.ylim(-5, 10)
    plt.title("Median binned $\sigma_{i}$ vs Median binned $\chi^2_{reduced}$ for " +str(filt) + " Band")
    plt.xlabel("median binned $\sigma_i$")
    plt.ylabel("median $\chi^2_{reduced}$ in bin")
    plt.savefig(image_place+"filter_"+str(filt)+"_chi2_clip_"+str(clip)+"_sigma_plot_dy_1_sqrt_n_line_no_last_1.pdf")
    #plt.show()

    plt.figure()
    plt.errorbar(j, k - test_func(j, params[0], params[1]), yerr=dk, label='a='+str(round(params[0], 4))+'\n'+'b='+str(round(params[1], 4)))
    plt.legend(loc='best')
    #plt.ylim(-5, 10)
    plt.title("Median binned $\sigma_{i}$ vs Trend-subtracted, Median binned $\chi^2_{reduced}$ for " +str(filt) + " Band")
    plt.xlabel("median binned $\sigma_i$")
    plt.ylabel("median $\chi^2_{reduced}$ in bin")
    plt.savefig(image_place+"filter_"+str(filt)+"_chi2_clip_"+str(clip)+"_sigma_plot_trend_subtract.pdf")
    #plt.show()

    a = params[0]
    b = params[1]


    clip= 10    
    high_m = 24
    low_m = 16
    chi_val = []
    sigma_med_val = []
    sigma_is = []
    flux_arr = []

    print("STARTING FILTER " + filt+ "\n")
    for fit in fits:
        if (fit['MAG_AUTO_G']- fit['MAG_AUTO_R']) > 0.4:
            if fit['MAG_AUTO_'+filt] < high_m:
                if np.sum(fit["LC_FLUX_PSF_"+filt] > 0) > 10:
                    if np.abs(fit['SPREAD_MODEL_'+filt]) < .003:
                        if fit['MAG_AUTO_'+filt] > low_m:
                            fluxes_arr = fit["LC_FLUX_PSF_"+filt]
                            ind_errs_arr = fit['LC_FLUXERR_PSF_'+filt]
                            fluxes, ind_errs = zip(*((x, y) for x, y in zip(fluxes_arr, ind_errs_arr) if x!= 0 and y!=0 and x<1e10))
                            fluxes = np.array(fluxes)
                            ind_errs = np.array(ind_errs)
                            mean_flux = fit['MEAN_PSF_'+filt]
                            norm_err = (ind_errs/fluxes)
                            chi = np.sum(((fluxes - mean_flux)**2)/(((a**2) *ind_errs * fluxes) +((b**2)* ind_errs**2)))/(len(fluxes)-1)
                            if chi > clip:
                                continue
                            for ind_err in ind_errs:
                                sigma_is.append(ind_err)
                            for flux in fluxes:
                                flux_arr.append(flux)
                            chi_val.append(chi)
                            sigma_med_val.append(np.median(norm_err))
                            #err.append(np.median(norm_err)/np.sqrt(len(fluxes)))
                            

    df = pd.DataFrame({'chi2r':chi_val, 'sig_m':sigma_med_val})
    #print(df)

    m, s = stats.norm.fit(chi_val)
    print("Regular values")
    print("Mean="+str(m))
    print("std="+str(s))

    m, s = stats.norm.fit(np.log10(np.array(chi_val))) # get mean and standard deviation 
    print("\nLogged values")
    print("Mean="+str(m))
    print("std="+str(s))

    plt.hist((np.array(chi_val)), bins=30)

    #plt.show()

    bins = np.linspace(df.sig_m.min(), df.sig_m.max(), 30)
    groups = df.groupby(pd.cut(df.sig_m, bins))
    #print(groups.size())

    x = np.array(groups.median().sig_m)
    #print(x)
    y = np.array(groups.median().chi2r)
    #print(y)

    #print(groups.chi2r)
    dy = np.array((groups.chi2r.quantile(.85) - groups.chi2r.quantile(.15))/(2*np.sqrt(groups.size())))
     
    j = x # sigma
    k = y # chi
    dk =dy # error
    jkVar = zip(j, k, dk)

    jkFiltered = [(j, k, dk) for j, k, dk in jkVar if j <= 0.12 and j >= 0.02]

    j, k, dk = zip(*jkFiltered)

    #j = x[:-1] # sigma
    #k = y[:-1] # chi
    #dk =dy[:-1]
    plt.figure()
    plt.errorbar(j, k, yerr=dk, label="a and b corrected")

    plt.legend(loc='best')
    #plt.ylim(-5, 10)
    plt.title("Median binned $\sigma_{i}$ vs Median binned $\chi^2_{reduced}$ for " +str(filt) + " Band")
    plt.xlabel("median binned $\sigma_i$")
    plt.ylabel("median $\chi^2_{reduced}$ in bin")
    plt.savefig(image_place+"filter_"+str(filt)+"_chi2_clip_"+str(clip)+"_sigma_plot_corrected_no_last_1.pdf")
    #plt.show()


