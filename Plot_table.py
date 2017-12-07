import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfit


count = 0
too_big = []
lis1 = []
lis2 = []
filter_G_R_V = []
filter_G_R_T = []
FITS = pyfit.open(r'C:\Users\Neal\My Documents\LiClipse Workspace\HW7\C1_lc.fits')
data = pd.read_csv(r'C:\Users\Neal\Desktop\ASTR406\table.txt', sep=" ", header=None)
data.columns = ["str1", "object", "str2", "Tau","str3","V"]


blue = (FITS[1].data['MAG_AUTO_G']- FITS[1].data['MAG_AUTO_R']) < 0.4
bright = FITS[1].data['MAG_AUTO_G'] < 23.5
long = np.sum(FITS[1].data["LC_FLUX_PSF_G"] > 0, axis = 1) > 10
point_like = np.abs(FITS[1].data['SPREAD_MODEL_G']) < .003
unsat = FITS[1].data['MAG_AUTO_G'] > 18

good = bright & long & point_like & unsat
goodobjects = np.array(range(good.size))[good]

print(np.sum(blue),np.sum(long),np.sum(bright),np.sum(point_like),np.sum(unsat))

tablegood = np.zeros(data.size)
num = 0
for objectid in data['object']:
    if objectid in goodobjects:
        tablegood[num] = 1
    num+= 1
tablegood = tablegood > 0
print((tablegood), data.size)
filter_G_R_V = (data['V'][tablegood])
filter_G_R_T = (data['Tau'][tablegood])
lis1 = (data['object'][tablegood])

for i in range(len(data['V'])):
    if i >= -1:
        count += 1
        too_big.append(data["object"][i])
plt.scatter(filter_G_R_T,filter_G_R_V)
#plt.figure()
#cs = plt.scatter(data['Tau'],data['V'], c = lis2, vmin= -1, vmax= 2)
#plt.colorbar(cs)

def boundaries(mjds,mags):
  mjdmin = 100*np.floor(np.min(mjds)*.01)
  mjdmax = 100*np.ceil(np.max(mags)*.01)
  [magmin,magmax] = np.percentile(mags,[2,98])
  magmin = 0.5*np.floor(2.*(magmin-0.1))
  magmax = 0.5*np.ceil(2.*(magmax+0.1))
  return [[mjdmin,mjdmax],[magmax,magmin]]

def getdata(fits,num):
  
  [mags_g, errs_g, mjds_g] = goodrow(22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_G']), fits[num]['LC_FLUXERR_PSF_G']/fits[num]['LC_FLUX_PSF_G'], fits[num]['LC_MJD_G'])
  [mags_r, errs_r, mjds_r] = goodrow(22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_R']), fits[num]['LC_FLUXERR_PSF_R']/fits[num]['LC_FLUX_PSF_R'], fits[num]['LC_MJD_R'])
  [mags_i, errs_i, mjds_i] = goodrow(22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_I']), fits[num]['LC_FLUXERR_PSF_I']/fits[num]['LC_FLUX_PSF_I'], fits[num]['LC_MJD_I'])
  [mags_z, errs_z, mjds_z] = goodrow(22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_Z']), fits[num]['LC_FLUXERR_PSF_Z']/fits[num]['LC_FLUX_PSF_Z'], fits[num]['LC_MJD_Z'])
  return [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z]
  
def goodrow(mags,errs,mjds):
  good = (mags > 0) & (mags != 22.5) & (mags != np.inf) & ~np.isnan(mags) & (errs > 0) & (errs != np.inf) & ~np.isnan(errs) & (mjds > 0) & (mjds != np.inf) & ~np.isnan(mjds)  
  return [mags[good], errs[good], mjds[good]]

def countoutliers(fits):
  fluxs = 22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_G']), fits[num]['LC_FLUXERR_PSF_G']/fits[num]['LC_FLUX_PSF_G'], fits[num]['LC_MJD_G'] 
def plotdata(file,VarRows):
  fits = pyfit.open(file)[1].data 
  for num in range(1):
    
    [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits,VarRows)
    rownum = num
    
    plt.rcParams['font.size'] = 18
    plt.figure()
    plt.subplot(111)
    plt.errorbar(mjds_g-57000, mags_g, yerr = errs_g, fmt = 'og')
    plt.errorbar(mjds_r-57000, mags_r, yerr = errs_r, fmt = 'or')
    plt.errorbar(mjds_i-57000, mags_i, yerr = errs_i, fmt = 'ok')
    plt.errorbar(mjds_z-57000, mags_z, yerr = errs_z, fmt = 'ob')
    [xlim,ylim] = boundaries(np.hstack([mjds_g,mjds_r,mjds_i,mjds_z]),np.hstack([mags_g,mags_r,mags_i,mags_z]))

    plt.ylim(ylim)
    plt.xlabel('MJD-57000')
    plt.ylabel('Mags')
    plt.title('C1'+' Object '+ str(VarRows))
    plt.savefig(r'C:\Users\Neal\Documents\LiClipse Workspace\Hw7astr\Fin_curves\\' + 'C1'+'_lc_'+ str(VarRows) +'.png')

#for item in lis1:
    #plotdata(r'C:\Users\Neal\My Documents\LiClipse Workspace\HW7\C1_lc.fits', item)

plt.axes().set_xlabel("$log_{10}(Tau)$")
plt.axes().set_ylabel("$log_{10}(V)$")

plt.show()
