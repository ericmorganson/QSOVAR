#!/usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfit
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as op
import corner
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



if len(sys.argv) < 6:
  print("lightcurveplot.py INFITS ROWSTART ROWEND BAND MCMC_TYPE")
  print("INFITS is a fits table of DES Data,")
  print("ROWNUM is a row in that file")
  print("MCMC_TYPE is how you set up the MCMC search array: 'grid' or 'normal'")
  sys.exit()
# These are the guesses that emcee starts with
V =0.3
Tau = 365.0
dMu = 0.0
print(sys.argv[2],sys.argv[3],sys.argv[4],type(sys.argv[3]))

def lognorm(state, flux, flux_err_sq):
        # This finds the normal distribution for any measurement (x),
        # mean(mu), and variance squared (var2).
        fmean, V_sq = state

        return -0.5*(((flux-fmean)**2)/(V_sq + flux_err_sq) + np.log(2*np.pi*(V_sq + flux_err_sq)))

def evolvestate(state, exp_dt_Tau, mu, V_sq_old):
        fmean, V_sq = state
        fmean_new = exp_dt_Tau*(fmean - mu) + mu
        V_sq_new = V_sq_old*(1 - exp_dt_Tau**2) + (exp_dt_Tau**2)*V_sq
        state = (fmean_new, V_sq_new)
        return state

def weightedmean(state, flux, flux_err_sq):
        fmean, V_sq = state
        denom = 1/(V_sq + flux_err_sq)
        fmean_new = ((flux)*V_sq + fmean*flux_err_sq)*denom
        V_sq_new = flux_err_sq*V_sq*denom
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
        lc_median = FITS[1].data['MEDIAN_PSF_'+color][ROW] #yes, I know this isn't the mean
        lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
        lc_time = FITS[1].data['LC_MJD_'+color][ROW]
        lc_time = lc_time[lc_time != 0.] #remove the zeros

        limit = len(lc_time)
        lc_flux = lc_flux[:limit]
        lc_flux_err = lc_flux_err[:limit]

        lc_flux_norm = (lc_flux - lc_median)/lc_median
        lc_flux_err_norm = lc_flux_err/lc_median

        return lc_flux_norm, lc_flux_err_norm, lc_time,lc_median, FITS

def lnlike(theta, time, flux, flux_err_sq):
        logV, logTau, logdMu = theta

        V_ten = 10**logV
        Tau_ten = 10**logTau
        dMu_ten = logdMu

        state=(mu + dMu_ten ,V_ten**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0])
        state = weightedmean(state, flux[0], flux_err_sq[0])
        #print(Tau_ten)
        #exp_dt_Tau = np.exp(-(np.subtract(time[1:], time[:-1])/Tau_ten)) #DIDN'T WORK DUE TO TAU BEING A 2D ARRAY... NEED TO FIGURE OUT HOW TO DIVIDE TIME ARRAY BY TAU WITHOUT FOR LOOP.
        #print(exp_dt_Tau.shape)
        for n in range(1,len(flux)):
            if time[n]-time[n-1] < 0:
                print('AHHHHH, NEGATIVE TIME!!!')
            exp_dt_Tau = np.exp(-(time[n] - time[n-1])/Tau_ten)
            state=evolvestate(state, exp_dt_Tau, mu + dMu_ten, V_ten**2) ### NEED TO WORK ON THIS
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        return lnp


def lnlike_old(theta, time, flux, flux_err_sq):
        logV, logTau = theta

        V_ten = 10**logV
        Tau_ten = 10**logTau

        state=(mu, V_ten**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0])
        state = weightedmean(state, flux[0], flux_err_sq[0])
        #print(Tau_ten)
        #exp_dt_Tau = np.exp(-(np.subtract(time[1:], time[:-1])/Tau_ten)) #DIDN'T WORK DUE TO TAU BEING A 2D ARRAY... NEED TO FIGURE OUT HOW TO DIVIDE TIME ARRAY BY TAU WITHOUT FOR LOOP.
        #print(exp_dt_Tau.shape)
        for n in range(1,len(flux)):
            if time[n]-time[n-1] < 0:
                print('AHHHHH, NEGATIVE TIME!!!')
            exp_dt_Tau = np.exp(-(time[n] - time[n-1])/Tau_ten)
            state=evolvestate(state, exp_dt_Tau, mu, V_ten**2) ### NEED TO WORK ON THIS
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        return lnp


def lnprior(theta):
        logV, logTau, logdMu = theta #,logC
        #print(logV, logTau)
        if -3 < logV < 2 and 0 < logTau < 4 and -3 < logdMu < 2:
            return 0.0
        return -np.inf

def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        logprobs.append(lp + lnlike(theta, x, y, yerr))
        logvals.append(theta)
        return lp + lnlike(theta, x, y, yerr)

def sausageplot(Vari,time,delta_f,Tau,dt,sigma_sq, ROW):
        #mu = fmean
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
            nsig=2
            sig1 = mu * (np.sqrt(sig_sq) + f_0) + mu
            sig2 = mu * (-np.sqrt(sig_sq) + f_0) + mu
            center = 22.5-2.5*np.log10((mu * (f_0))+mu)
            err_t = 22.5-2.5*np.log10((mu * (f_0 + nsig*np.sqrt(sig_sq)))+mu)
            err_b = 22.5-2.5*np.log10((mu * (f_0 - nsig*np.sqrt(sig_sq)))+mu)
            Logpr.append(center)
            err_top.append(err_t)
            err_bot.append(err_b)

        plotdata(sys.argv, ROW)
        plt.plot(times,err_top, color = 'g')
        plt.plot(times,Logpr, color = 'b')
        plt.plot(times,err_bot, color = 'g')
        plt.title(' Object '+sys.argv[2])
        plt.savefig('figure/sausage_mu_' + 'C1_lc' + sys.argv[5] + sys.argv[2] + '.png')

        #plt.show()

def lnprob_dens(theta, x, y, yerr):
        logprobs_dens.append(lnlike_old(theta, x, y, yerr))
        logvals_dens.append(theta)
        return logprobs_dens

def goodrow(mags,errs,mjds):
        good = (mags > 0) & (mags != 22.5) & (mags != np.inf) & ~np.isnan(mags) & (errs > 0) & (errs != np.inf) & ~np.isnan(errs) & (mjds > 0) & (mjds != np.inf) & ~np.isnan(mjds)
        return [mags[good], errs[good], mjds[good]]

def plotdata(args, ROW):
        fits = pyfit.open(args[1])[1].data
        VarRows = ROW #int(args[2])
        [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits, ROW) #17999
        rownum = str(ROW) #args[2]

        plt.rcParams['font.size'] = 18
        plt.figure()
        plt.subplot(111)
        plt.errorbar(mjds_g-57000, mags_g, yerr = errs_g, fmt = 'og')
        #plt.errorbar(mjds_r-57000, mags_r, yerr = errs_r, fmt = 'or')
        #plt.errorbar(mjds_i-57000, mags_i, yerr = errs_i, fmt = 'ok')
        #plt.errorbar(mjds_z-57000, mags_z, yerr = errs_z, fmt = 'ob')
        [xlim,ylim] = boundaries(np.hstack([mjds_g,mjds_r,mjds_i,mjds_z]),np.hstack([mags_g,mags_r,mags_i,mags_z]))

        plt.ylim(ylim)
        plt.xlabel('MJD-57000')
        plt.ylabel('Mags')
        field = args[1].split('_')[0]
        plt.title(field+' Object '+rownum)

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

def preform_emcee(time,flux,sigma_sq,ROW):
        diff_time = [x - time[i - 1] for i, x in enumerate(time)][1:]
        print(min(diff_time))
        plt.figure()
        plt.errorbar(time, flux, err)
        plt.savefig('figure/'+ str(ROW) + 'LC_mu' + '.pdf')
        #plt.show()

        X = np.arange(-1, 5, .1) #tau
        Y = np.arange(-2.5, 1.5, .1) #variance
        X, Y = np.meshgrid(X, Y)
        lprob_dens = lnprob_dens((Y, X), time, flux, err)
        fig = plt.figure()
        #print(np.array(lprob_dens).shape())
        lprob_dens=np.array(lprob_dens)
        print(lprob_dens.shape)
        plt.pcolormesh(X, Y, lprob_dens.reshape(X.shape), shading='gouraud', cmap=cm.rainbow)
        cbar = plt.colorbar()
        cbar.set_label('log(probability)')
        plt.xlabel("Tau")
        plt.ylabel("Variance")
        plt.savefig('figure/'+ str(ROW) + 'logprob_density_norm_mu' + '.pdf')
        plt.show()

        nll = lambda *args: -lnlike(*args)
        ndim, nwalkers = 3, 100

        if sys.argv[5].lower() == 'normal':
            result = [np.log10(V), np.log10(Tau), dMu]
#            pos = [result + (-0.5+np.random.randn(ndim)) for i in range(nwalkers)]
            pos = (np.random.rand(100,3)-0.5)*np.array([1,1,0.2])+result
        elif sys.argv[5].lower() == 'grid':
            v_grid = np.arange(-1, 0, 0.01)
            t_grid = np.arange(1, 2, 0.01)
            dmu_grid = np.arange(-0.5, 0.5, 0.01)
            VG, TG, MG = np.meshgrid(v_grid, t_grid, dmu_grid)
            result = [np.array(thing) for thing in zip(VG.flatten(), TG.flatten(), MG.flatten())] # for python 2.7
            pos = [result[i] + 1e-7*np.random.randn(ndim) for i in range(nwalkers)]
        elif sys.argv[5].lower() == 'optimal':
            result = op.minimize(nll, [np.log10(V), np.log10(Tau), dMu],args=(time,flux, err**2))
            pos = [result['x'] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        else:
            print("What the hell do you want to do?")
            print("'grid', 'optimal', or 'normal' search through MCMC?")
            exit()

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, flux, err**2))
        print(np.array(pos).shape)
        sampler.run_mcmc(pos, 100)
        samples = sampler.chain[:, 20:, :].reshape((-1, ndim))

        plt.figure()
        plt.plot(logprobs)
        plt.savefig('figure/'+ str(ROW) + sys.argv[5]+ 'logprob_mu' + '.pdf')
        #plt.show()

        max_theta = logvals[logprobs.index(max(logprobs))]
        fig = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$",r"log$_{10}d\mu$"],
                                                    truths=[max_theta[0], max_theta[1], max_theta[2]])

        fig.savefig("figure/"+str(ROW)+sys.argv[5]+"triangle_np_mu.pdf")

        V_mcmc, Tau_mcmc, dMu_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
        print('V_mcmc:',V_mcmc, 'Tau_mcmc:',Tau_mcmc, 'dMu_mcmc:', dMu_mcmc, max_theta[0], max_theta[1], max_theta[2])
        print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]), 'dMu:', str(max_theta[2]))
        filename ='scratch_new/'+ str(ROW) + sys.argv[5]+'object_dMu' + '.txt'
        with open(filename, 'w+') as fout:
            fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')
            fout.write('Object: ' + str(ROW)+ ' ' + 'dMu: ' + str(max_theta[2])+ '\n')
        sausageplot(max_theta[0],time,flux,max_theta[1],5,err**2, ROW)


for ROW in range(int(sys.argv[2]),int(sys.argv[3])):
    logprobs = []
    logvals = []

    logprobs_dens = []
    logvals_dens = []


    print(ROW)
    flux, err, time, mu, FITS = get_vals(sys.argv,ROW)
    if len(flux) == 0:
        print("Flux length is zero")
        continue

    #M,delta_f,sigma_sq = Make_M(V,Tau,C)

    try:
        preform_emcee(time, flux, err, ROW)
    except np.linalg.linalg.LinAlgError as err:
        continue
