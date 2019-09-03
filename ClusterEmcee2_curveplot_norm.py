#!/usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfit
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

if len(sys.argv) < 5:
  print("lightcurveplot.py INFITS ROWSTART ROWEND BAND")
  print("INFITS is a fits table of DES Data,")
  print("ROWNUM is a row in that file")
  sys.exit()
# These are the guesses that emcee starts with
V =0.3
#C = .01
Tau = 365.0
#Mu = 0.01
print(sys.argv[2],sys.argv[3],sys.argv[4],type(sys.argv[3]))

def lognorm(state, flux, flux_err_sq):
        # This finds the normal distribution for any measurement (x), 
        # mean(mu), and variance squared (var2).
        fmean, V_sq = state
        #return np.log(np.exp(-((flux-fmean)**2)/(2*(V_sq + flux_err_sq)))/np.sqrt(2*np.pi*(V_sq + flux_err_sq))) #good & stable; didn't use normalized flux or errors
        #return -(flux**2/(2*(V_sq + flux_err_sq))) + (flux*fmean/(V_sq + flux_err_sq)) - (fmean**2/(2*(V_sq + flux_err_sq))) -0.5*np.log(2*np.pi*(V_sq + flux_err_sq)) # assume all values are real; stable; didn't use norm. flux or errors
        return np.log(np.exp(-((flux*fmean)**2)/(2*(V_sq + flux_err_sq*fmean**2)))/np.sqrt(2*np.pi*(V_sq + flux_err_sq*fmean**2)))

def evolvestate(state, Tau, dt, mu, V_sq_old):
        fmean, V_sq = state
        fmean_new = np.exp(-dt/Tau)*(fmean - mu) + mu    
        V_sq_new = V_sq_old*(1 - np.exp(-2*dt/Tau)) + np.exp(-2*dt/Tau)*V_sq
        state = (fmean_new, V_sq_new)
        return state

def weightedmean(state, flux, flux_err_sq):
        fmean, V_sq = state
        fmean_new = ((flux*fmean + fmean)*V_sq + (fmean**3)*flux_err_sq)/(V_sq + flux_err_sq*fmean**2)
        V_sq_new = flux_err_sq*(fmean**2)*V_sq/(flux_err_sq*(fmean**2) + V_sq)
        state = (fmean_new, V_sq_new)
        return state
    

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
        lc_mean = FITS[1].data['MEDIAN_PSF_'+color][ROW] #yes, I know this isn't the mean
        lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
        lc_time = FITS[1].data['LC_MJD_'+color][ROW]
           #remove the zeros
        lc_time = lc_time[lc_time != 0.]

        limit = len(lc_time)
        lc_flux = lc_flux[:limit]
        lc_flux_err = lc_flux_err[:limit]

        lc_flux_norm = (lc_flux - lc_mean)/lc_mean
        lc_flux_err_norm = lc_flux_err/lc_mean

        return lc_flux_norm, lc_flux_err_norm, lc_time,lc_mean, FITS

def lnlike(theta, time, flux, flux_err_sq):
        logV, logTau = theta
        V = 10**logV
        Tau = 10**logTau

        
        state=(mu,V**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0]) 
        state = weightedmean(state, flux[0], flux_err_sq[0])
        for n in range(1,len(flux)):
            state=evolvestate(state, Tau, time[n]-time[n-1], mu, V**2) ### NEED TO WORK ON THIS
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        #print(lnp)
        return lnp


def lnprior(theta):
        logV, logTau = theta #,logC
        #print(logV, logTau)
        if -3 < logV < 2 and 0 < logTau < 4: # and -3 < logMu < 2:
            return 0.0
        return -np.inf


def lnprob(theta, x, y, yerr):
        #lp = lnprior(theta)
        #print(lp)
        #if not np.isfinite(lp):
        #    return -np.inf
        #print(lp + lnlike(theta, x, y, yerr))
        logprobs.append(lnlike(theta, x, y, yerr))
        logvals.append(theta)
        return lnlike(theta, x, y, yerr)

def preform_emcee(time,delta_f,sigma_sq,ROW):
        flux, err, time, mu, FITS = get_vals(sys.argv, ROW)
        #print(flux, err, time, mu)
        #M,delta_f,sigma_sq = Make_M(V,Tau,C)

        X = np.arange(0, 4, .1) #tau
 
        Y = np.arange(-3, 2, .1) #variance
        X, Y = np.meshgrid(X, Y)
                #theta = lnV, lnT
        lprob = lnprob((Y, X), time, flux, err)
                #print(l)

        fig = plt.figure()
        #ax = fig.gca(projection='3d')
        plt.contour(X, Y, lprob, cmap=cm.rainbow)#,
                                #linewidth=0, antialiased=False)
        plt.savefig('/home/sam/Documents/Morganson_research/QSOVAR/DESVAR/'+ str(ROW) + 'logprob_contour_norm' + '.pdf')
        plt.xlabel("Tau")
        plt.ylabel("Variance")
        plt.show()
        
        #result = op.minimize(nll, [np.log10(V), np.log10(Tau)],args=(time,flux, err)) #,np.log10(C)
        #ndim, nwalkers = 2, 100
        #pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, flux, err))

        #sampler.run_mcmc(pos, 500)
        #samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        #print('logprobs', logprobs)
        #max_theta = logvals[logprobs.index(max(logprobs))]
        #V_mcmc, Tau_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
        #plotdata(sys.argv)
        #print('V_mcmc:',V_mcmc, 'Tau_mcmc:',Tau_mcmc, max_theta[0], max_theta[1])
        #print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]))
        #filename ='/home/sam/Documents/Morganson_research/QSOVAR/scratch_new/'+ str(ROW) + 'object' + '.txt' 
        #with open(filename, 'w') as fout:
        #    fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')


for ROW in range(int(sys.argv[2]),int(sys.argv[3])):
    logprobs = []
    logvals = []
    print(ROW)
    flux, err, time, mu, FITS = get_vals(sys.argv,ROW)

    #M,delta_f,sigma_sq = Make_M(V,Tau,C)

    #try:
    preform_emcee(time, flux, err, ROW)
    #except np.linalg.linalg.LinAlgError as err:
    #    continue
