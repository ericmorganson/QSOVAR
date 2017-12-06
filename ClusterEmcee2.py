#!/usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfit
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op

if len(sys.argv) < 5:
  print("lightcurveplot.py INFITS ROWSTART ROWEND BAND")
  print("INFITS is a fits table of DES Data,")
  print("ROWNUM is a row in that file")
  sys.exit()
# These are the guesses that emcee starts with
V =0.3
C = .01
Tau = 365.0
print(sys.argv[2],sys.argv[3],sys.argv[4],type(sys.argv[3]))
def get_vals(args, ROW):
        FITS = pyfit.open(args[1])
        color = sys.argv[4]
        color = color.upper()
            #get the flux, err, and time arrays
        try:
            lc_flux = FITS[1].data['LC_FLUX_PSF_'+color][ROW]
        except IndexError:
            print('Error')
            exit()
        lc_mean = FITS[1].data['MEAN_PSF_'+color][ROW]
        lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
        lc_time = FITS[1].data['LC_MJD_'+color][ROW]
           #remove the zeros
        lc_time = lc_time[lc_time != 0.]

        limit = len(lc_time)
        lc_flux = lc_flux[:limit]
        lc_flux_err = lc_flux_err[:limit]

        return lc_flux, lc_flux_err, lc_time,lc_mean, FITS

def Make_M(V,Tau,C):
        delta_f = (flux - fmean)/fmean
        sigma_sq = (err * err) / (fmean**2)
        time_array = np.fabs(np.array(time,ndmin=2)-np.array(time,ndmin=2).T)
        M = np.diag(sigma_sq+C**2) + (V ** 2) * np.exp(-1.0 * (time_array)/Tau)
        return M,delta_f,sigma_sq

def lnlike(theta, time, delta_f, sigma_sq):
        logV, logTau, logC = theta
        inv_sigma2 = sigma_sq ** -1
        M,delta_f,sigma_sq = Make_M(10**logV,10**logTau,10**logC)
        sign, value = np.linalg.slogdet(M)
        deter = sign * np.exp(value.item())
        lnP = -.5*np.log(deter)-.5*((np.dot(delta_f,(np.dot(delta_f,np.linalg.inv(M))))))
        return lnP

def lnprior(theta):
        logV, logTau, logC = theta
        if -3 < logV < 2 and 0 < logTau < 4 and -3 < logC < 2:
            return 0.0
        return -np.inf

def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        logprobs.append(lp + lnlike(theta, x, y, yerr))
        logvals.append(theta)
        return lp + lnlike(theta, x, y, yerr)

def preform_emcee(time,delta_f,sigma_sq,ROW):
        flux, err, time, fmean, FITS = get_vals(sys.argv, ROW)
        M,delta_f,sigma_sq = Make_M(V,Tau,C)

        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [np.log10(V), np.log10(Tau),np.log10(C)],args=(time,delta_f,sigma_sq))
        ndim, nwalkers = 3, 100
        pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, delta_f, sigma_sq))

        sampler.run_mcmc(pos, 500)
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        max_theta = logvals[logprobs.index(max(logprobs))]
        V_mcmc, Tau_mcmc, C_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
        #plotdata(sys.argv)
        print('V_mcmc:',V_mcmc, 'Tau_mcmc:',Tau_mcmc, max_theta[0], max_theta[1])
        print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]))
        filename = '/scratch/users/nschwei2/'+ str(ROW) + 'object' + '.txt'
        with open(filename, 'w') as fout:
            fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')
for ROW in range(int(sys.argv[2]),int(sys.argv[3])):
    logprobs = []
    logvals = []
    print(ROW)
    flux, err, time, fmean, FITS = get_vals(sys.argv,ROW)

    M,delta_f,sigma_sq = Make_M(V,Tau,C)

    try:
        preform_emcee(time, delta_f, sigma_sq, ROW)
    except np.linalg.linalg.LinAlgError as err:
        continue
