#!/usr/bin/env python
#import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import astropy.io.fits as pyfit
from scipy import stats, optimize
np.set_printoptions(precision=5)

#initialize variable
files = ['/home/sam/Documents/Morganson_research/C1_lc.fits', '/home/sam/Documents/Morganson_research/C2_lc.fits'] #,
        # '/home/sam/Documents/Morganson_research/E1_lc.fits', '/home/sam/Documents/Morganson_research/E2_lc.fits',
        # '/home/sam/Documents/Morganson_research/S1_lc.fits', '/home/sam/Documents/Morganson_research/S2_lc.fits',
        # '/home/sam/Documents/Morganson_research/X1_lc.fits', '/home/sam/Documents/Morganson_research/X2_lc.fits']

clip= 10
high_m = 24
low_m = 16
chi_val = []
sigma_med_val = []
sigma_is = []
flux_arr = []

def find_a_b(fits, LC_file, color):
    for fit in fits:
        if (fit['MAG_AUTO_G']- fit['MAG_AUTO_R']) > 0.4:
            if fit['MAG_AUTO_'+color] < high_m:
                if np.sum(fit["LC_FLUX_PSF_G"] > 0) > 10:
                    if np.abs(fit['SPREAD_MODEL_'+color]) < .003:
                        if fit['MAG_AUTO_'+color] > low_m:
                            fluxes = fit["LC_FLUX_PSF_"+color]
                            fluxes = fluxes[fluxes!=0]
                            mean_flux = fit['MEAN_PSF_'+color]
                            ind_errs = fit['LC_FLUXERR_PSF_'+color]
                            ind_errs = ind_errs[ind_errs!=0]
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

    #m, s = stats.norm.fit(chi_val)
    #print("Regular values " + LC_file + " " +color)
    #print("Mean="+str(m))
    #print("std="+str(s))

    #m, s = stats.norm.fit(np.log10(np.array(chi_val))) # get mean and standard deviation
    #print("\nLogged values")
    #print("Mean="+str(m))
    #print("std="+str(s))

    #plt.hist((np.array(chi_val)), bins=20)

    #plt.show()

    bins = np.linspace(df.sig_m.min(), df.sig_m.max(), 20)
    groups = df.groupby(pd.cut(df.sig_m, bins))
    #print(groups.size())

    x = np.array(groups.median().sig_m)
    #print(x)
    y = np.array(groups.median().chi2r)
    #print(y)

    #print(groups.chi2r)
    dy = np.array((groups.chi2r.quantile(.85) - groups.chi2r.quantile(.15))/(2*np.sqrt(groups.size())))
    #print(dy)

    j = x[:-1] # sigma
    k = y[:-1] # chi
    dk =dy[:-1]
    plt.figure()
    plt.errorbar(j, k, yerr=dk, label="clipped")
    #plt.errorbar(x, y, yerr=dy, label="all")
    #plt.plot(x[:-3], 0.005**2 + 3*x[:-3]**2)

    z = np.polyfit(j, k, 1)
    #print(z)
    f = np.poly1d(z)

    m = np.polyfit(j, k, 1, w=1/dk)
    #print(m)
    g = np.poly1d(m)

    def test_func(sigma, a, b):
        return a**2/(sigma**1) + b**2

    params, params_covariance = optimize.curve_fit(test_func, j, k, p0=[1, 3],
                                                   sigma=dk)

    print(LC_file, color, params[0], params[1])

    # calculate new x's and y's
    x_new = np.linspace(j[0], j[-1], 50)
    y_new = f(x_new)
    y_new2 = g(x_new)
    #plt.plot(x_new, y_new, label="no weight "+str(z))
    #plt.plot(x_new, y_new2, label="weighted "+str(m))
    plt.plot(x_new, test_func(x_new, params[0], params[1]), label='a='+str(round(params[0], 4))+'\n'+'b='+str(round(params[1], 4)))
    plt.legend(loc='best')
    #plt.ylim(-5, 10)
    plt.title("Median binned $\sigma_{i}$ vs Median binned $\chi^2_{reduced}$")
    plt.xlabel("median binned $\sigma_i$")
    plt.ylabel("median $\chi^2_{reduced}$ in bin")
    plt.savefig("chi2_"+LC_file+"_"+color+".pdf")
    #plt.show()



for file in files:
    fits = pyfit.open(file)[1].data
    LC_file = file.split('/')[5][0:2]
    for color in 'GRIZ':
        find_a_b(fits, LC_file, color)
    fits.close()
    
