import sys
if len(sys.argv) < 3:
  print("lightcurveplot.py INFITS ROWNUM")
  print("INFITS is a fits table of DES Data,") 
  print("ROWNUM is a row in that file") 
  sys.exit()
import numpy as np
import astropy.io.fits as pyfit
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import corner
# These are the guesses that emcee starts with
V =0.3
C = .01
Tau = 365.0
logprobs = []
logvals = []
#Make the lightcurve

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

def plotdata(args):
    fits = pyfit.open(args[1])[1].data 
    VarRows = int(args[2])
    [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits,17999)
    rownum = args[2]
    
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
    field = args[1].split('_')[0]
    plt.title(field+' Object '+rownum)
 

#Here is where the prep from emcee starts
def get_vals(args):
    FITS = pyfit.open(args[1])
    
    #get the flux, err, and time arrays
    row = int(args[2])
    lc_flux = FITS[1].data['LC_FLUX_PSF_G'][row]
    lc_mean = FITS[1].data['MEDIAN_PSF_G'][row]
    lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_G'][row]
    lc_time = FITS[1].data['LC_MJD_G'][row]
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

flux, err, time, fmean, FITS = get_vals(sys.argv)

M,delta_f,sigma_sq = Make_M(V,Tau,C)

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


def sausageplot(Vari,time,delta_f,Tau,dt,sigma_sq):
    mu = fmean
    err_top = []
    err_bot = []
    Logpr = []
    Vari = 10 ** Vari
    Tau = 10 ** Tau
    times = []
   
    for t in np.arange(min(time),max(time),dt):
        sigma_sq_s = np.append(sigma_sq, 0.0)
        dtime = np.append(time,t)
        delta_time = np.fabs(np.array(dtime,ndmin=2)-np.array(dtime,ndmin=2).T)
        #Make the modified matrix
        S = np.diag(sigma_sq_s) + (Vari**2) * np.exp(-1.0 * delta_time / Tau)
        times.append(t-57000)
        #Set up our guess for the flux
        t0 = [i for i in time if i >= t]
        t1 = [i for i in time if i <= t]
        t0 = t0[0]
        time_ind0 = np.where(time == t0)
        t1 = t1[-1]
        time_ind1 = np.where(time == t1)
        tp_0 = t0 / t
        tp_1 = t1 / t
        Fg1 = (tp_0 * delta_f[time_ind0]) + (tp_1 * delta_f[time_ind1])
        Fg1 = (Fg1 - mu)/mu
        Fg2 = Fg1 - 1
        Fg3 = Fg1 + 1
       
        delf1 = np.append(delta_f,Fg1)
        delf2 = np.append(delta_f,Fg2)
        delf3 = np.append(delta_f,Fg3)
       
        sign, value = np.linalg.slogdet(S)
        deter = sign * np.exp(value.item())
        #calculate the Log Probs
        logP = -.5*np.log((deter))-.5*((np.dot(delf1,(np.dot(delf1,np.linalg.inv(S))))))
        logP1 = -.5*np.log((deter))-.5*((np.dot(delf2,(np.dot((delf2),np.linalg.inv(S))))))
        logP2 = -.5*np.log((deter))-.5*((np.dot(delf3,(np.dot((delf3),np.linalg.inv(S)))))) 

        X = [Fg2,Fg1,Fg3]
        Y = [logP1,logP,logP2]
        X = np.array(X)
        X = (X.T)[0]
        Matr = np.array([[X[0]**2,X[0],1],[X[1]**2,X[1],1],[X[2]**2,X[2],1]])
        #Matr = Matr.astype(np.float64)
        result = np.linalg.solve(Matr, Y) 
        sig_sq = -1/(2*result[0])
        f_0 = -(result[1]/(2*result[0]))
        
        parabola = np.polyfit(X,Y,2)
        f = np.poly1d(parabola)
       
        x_new = np.linspace(X[0], X[-1], 5000)
        y_new = f(x_new)
       
        delf = delta_f
        Fm = op.fmin(lambda x: -f(x),0)
        sig = 22.5-2.5*np.log10((sig_sq))
        logPn = result[0]*((Fm - f_0)**2) + result[2] - (result[1]**2)/(4*result[0])
        
        sig1 = mu * (np.sqrt(sig_sq) + f_0) + mu
        sig2 = mu * (-np.sqrt(sig_sq) + f_0) + mu
        center = 22.5-2.5*np.log10((mu * (f_0))+mu)
        err_t = 22.5-2.5*np.log10((mu * (f_0 + np.sqrt(sig_sq)))+mu)
        err_b = 22.5-2.5*np.log10((mu * (f_0 - np.sqrt(sig_sq)))+mu)
        Logpr.append(center)
        err_top.append(err_t)
        err_bot.append(err_b)
        
    plotdata(sys.argv)
    plt.plot(times,err_top, color = 'g')
    plt.plot(times,Logpr, color = 'b')
    plt.plot(times,err_bot, color = 'g')
    plt.title(' Object '+sys.argv[2])
    plt.savefig('sausage' + 'C1_lc' +  sys.argv[2] + '.png')

    plt.show()
def preform_emcee(time,delta_f,sigma_sq):
    flux, err, time, fmean, FITS = get_vals(sys.argv)
    M,delta_f,sigma_sq = Make_M(V,Tau,C)
    
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [np.log10(V), np.log10(Tau),np.log10(C)],args=(time,delta_f,sigma_sq))
    ndim, nwalkers = 3, 100
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, delta_f, sigma_sq))

    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    max_theta = logvals[logprobs.index(max(logprobs))]
    V_mcmc, Tau_mcmc, C_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                   zip(*np.percentile(samples, [16, 50, 84],
                                                      axis=0)))
    fig = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$",r"log$_{10}$C"],
                                                truths=[max_theta[0], max_theta[1], max_theta[2]])
            #make the png
    fig.savefig("triangle.png")
    plotdata(sys.argv)
    print(max_theta)
    print(V_mcmc,Tau_mcmc,C_mcmc)
    sausageplot(max_theta[0], time, delta_f, max_theta[1], 5, sigma_sq)
preform_emcee(time, delta_f, sigma_sq)
    
