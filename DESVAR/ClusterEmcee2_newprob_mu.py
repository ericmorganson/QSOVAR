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
dMu = 0.1
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

        #due to normalization, mu is by definition 0
        state=(dMu_ten ,V_ten**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0])
        state = weightedmean(state, flux[0], flux_err_sq[0])
        #print(Tau_ten)
        #exp_dt_Tau = np.exp(-(np.subtract(time[1:], time[:-1])/Tau_ten)) #DIDN'T WORK DUE TO TAU BEING A 2D ARRAY... NEED TO FIGURE OUT HOW TO DIVIDE TIME ARRAY BY TAU WITHOUT FOR LOOP.
        #print(exp_dt_Tau.shape)
        for n in range(1,len(flux)):
            if time[n]-time[n-1] < 0:
                print('AHHHHH, NEGATIVE TIME!!!')
            exp_dt_Tau = np.exp(-(time[n] - time[n-1])/Tau_ten)
            state=evolvestate(state, exp_dt_Tau, dMu_ten, V_ten**2) ### NEED TO WORK ON THIS
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
        logV, logTau, dMu = theta #,logC
        #print(logV, logTau)
        if -3 < logV < 2 and 0 < logTau < 4 and -0.5 < dMu < 0.5 :
            return 0.0
        return -np.inf

def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        logprobs.append(lp + lnlike(theta, x, y, yerr))
        logvals.append(theta)
        return lp + lnlike(theta, x, y, yerr)

def sausageplot(Vari,time,delta_f,Tau,dt,sigma_sq, ROW, fig):
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
        ax4 = fig.add_subplot(2, 2, 4)
        plotdata(sys.argv, ROW, ax4)
        ax4.plot(times,err_top, color = 'g')
        ax4.plot(times,Logpr, color = 'b')
        ax4.plot(times,err_bot, color = 'g')
        ax4.set_title(' Object '+str(ROW)+ " Sausage Plot")
        #plt.savefig('figure/sausage_mu_' + 'C1_lc' + sys.argv[5] + str(ROW) + '.png')
        #plt.show()

def lnprob_dens(theta, x, y, yerr):
        logprobs_dens.append(lnlike(theta, x, y, yerr))
        logvals_dens.append(theta)
        return logprobs_dens

def goodrow(mags,errs,mjds):
        good = (mags > 0) & (mags != 22.5) & (mags != np.inf) & ~np.isnan(mags) & (errs > 0) & (errs != np.inf) & ~np.isnan(errs) & (mjds > 0) & (mjds != np.inf) & ~np.isnan(mjds)
        return [mags[good], errs[good], mjds[good]]

def plotdata(args, ROW, ax4):
        fits = pyfit.open(args[1])[1].data
        VarRows = ROW #int(args[2])
        [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits, ROW) #17999
        rownum = str(ROW) #args[2]

        #ax4.rcParams['font.size'] = 18
        #plt.figure()
        #plt.subplot(111)
        ax4.errorbar(mjds_g-57000, mags_g, yerr = errs_g, fmt = 'og')
        #plt.errorbar(mjds_r-57000, mags_r, yerr = errs_r, fmt = 'or')
        #plt.errorbar(mjds_i-57000, mags_i, yerr = errs_i, fmt = 'ok')
        #plt.errorbar(mjds_z-57000, mags_z, yerr = errs_z, fmt = 'ob')
        [xlim,ylim] = boundaries(np.hstack([mjds_g]),np.hstack([mags_g]))#,mjds_r,mjds_i,mjds_z #,mags_r,mags_i,mags_z

        ax4.set_ylim(ylim)
        ax4.set_xlabel('MJD-57000')
        ax4.set_ylabel('Mags')
        field = args[1].split('_')[0]
        ax4.set_title(field+' Object '+rownum)

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
        fig = plt.figure(figsize=(10,10))

        nll = lambda *args: -lnlike(*args)
        ndim, nwalkers = 3, 100

        if sys.argv[5].lower() == 'normal':
            result = [np.log10(V), np.log10(Tau), dMu]
#            pos = [result + (-0.5+np.random.randn(ndim)) for i in range(nwalkers)]
            pos = (np.random.rand(100,3)-0.5)*np.array([1,1,0.2])+result
        elif sys.argv[5].lower() == 'grid':
            v_grid = np.arange(-1, 0, 0.01)
            t_grid = np.arange(1, 2, 0.01)
            dmu_grid = np.arange(-0.1, 0.1, 0.01)
            VG, TG, MG = np.meshgrid(v_grid, t_grid, dmu_grid)
            result = [np.array(thing) for thing in zip(VG.flatten(), TG.flatten(), MG.flatten())] # for python 2.7
            #pos = [result[i] + 1e-7*np.random.randn(ndim) for i in range(nwalkers)]
            pos = (np.random.rand(len(result),3)-0.5)*np.array([1,1,0.2])+result
        elif sys.argv[5].lower() == 'optimal':
            result = op.minimize(nll, [np.log10(V), np.log10(Tau), dMu],args=(time,flux, err**2))
            pos = [result['x'] + 1e-4*np.random.rand(ndim) for i in range(nwalkers)]
        else:
            print("What the hell do you want to do?")
            print("'grid', 'optimal', or 'normal' search through MCMC?")
            exit()
        #print(pos.shape())
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, flux, err**2))
        print(np.array(pos).shape)
        sampler.run_mcmc(pos, 200)
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

        #plt.figure()
        ax3 = fig.add_subplot(2, 2, 2)
        ax3.plot(logprobs)
        ax3.set_title("MCMC burn-in log(probability)")
        #plt.savefig('figure/'+ str(ROW) + sys.argv[5]+ 'logprob_mu' + '.pdf')

        max_theta = logvals[logprobs.index(max(logprobs))]
        fig1 = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$",r"$d\mu$"],
                            truths=[max_theta[0], max_theta[1], max_theta[2]])

        fig1.savefig("figure/"+str(ROW)+sys.argv[5]+"triangle_np_mu.pdf")

        V_mcmc, Tau_mcmc, dMu_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84],axis=0)))
        print('V_mcmc:',V_mcmc, 'Tau_mcmc:',Tau_mcmc, 'dMu_mcmc:', dMu_mcmc, max_theta[0], max_theta[1], max_theta[2])
        #stack = np.stack(logvals, axis=0)
        stack = np.asarray(logvals)
        V_stack = [x[0] for x in stack]
        T_stack = [x[1] for x in stack]
        print(stack[-1])
        print(V_stack[-1])
        print(T_stack[-1])
        ax1 = fig.add_subplot(2,2,3)
        a1 = ax1.hexbin(T_stack, V_stack, gridsize=(25,25), cmap=cm.viridis )
        cbar = fig.colorbar(a1, ax=ax1)
        ax1.set_ylabel(r"log$_{10}V$")
        ax1.set_xlabel(r"log$_{10}\tau$")
        ax1.set_title("MCMC heatmap of V and Tau")

        #print(logvals[:-1])

        X = np.arange(-1, 5, .1) #tau
        Y = np.arange(-2.5, 1.5, .1) #variance
        X, Y = np.meshgrid(X, Y)
        lprob_dens = lnprob_dens((Y, X, max_theta[2]), time, flux, err)
        lprob_dens=np.array(lprob_dens)
        #lprob_dens[lprob_dens < -100] = -100


        ax2 = fig.add_subplot(2, 2, 1)
        #print(lprob_dens.shape)
        a2 = ax2.pcolormesh(X, Y, lprob_dens.reshape(X.shape), shading='gouraud', cmap=cm.viridis)
        cbar = fig.colorbar(a2, ax=ax2)
        cbar.set_clim(-100,np.max(lprob_dens))
        cbar.set_label('log(probability)')
        ax2.set_ylabel(r"log$_{10}V$")
        ax2.set_xlabel(r"log$_{10}\tau$")
        ax2.set_title("Probability Density")
        #plt.savefig('figure/'+ str(ROW) + 'logprob_density_norm_mu' + '.pdf')
        #plt.show()

        ax2.scatter(max_theta[1], max_theta[0], label="max_theta")
        ax2.legend()
        print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]), 'dMu:', str(max_theta[2]))
        sausageplot(max_theta[0],time,flux,max_theta[1],5,err**2, ROW, fig)
        fig.savefig("figure/"+str(ROW)+sys.argv[5]+"all_plot_mu.pdf")
        filename ='scratch_new/'+ str(ROW) + sys.argv[5]+'object_dMu' + '.txt'
        with open(filename, 'w+') as fout:
            fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')
            fout.write('Object: ' + str(ROW)+ ' ' + 'dMu: ' + str(max_theta[2])+ '\n')



for ROW in range(int(sys.argv[2]),int(sys.argv[3])):
    logprobs = []
    logvals = []

    logprobs_dens = []
    logvals_dens = []


    print(ROW)
    flux, err, time, mu, FITS = get_vals(sys.argv,ROW)


    print(len(flux))
    if len(flux) == 0:
        print("Flux length is zero")
        continue
    zip_tfe = zip(time, flux, err)
    filter_tfe = [(t, f, e) for t, f, e in zip_tfe if f > -1 and e > 0]
    time, flux, err = zip(*filter_tfe)
    time = np.array(time)
    flux = np.array(flux)
    err = np.array(err)

    #if float(22.5-2.5*np.log10(mu)) > 21:
    #    print("Row is too dim")
    #    continue
    #if float(22.5-2.5*np.log10(mu)) < 16:
    #    print("Row is HELLA bright")
    #    continue

    #M,delta_f,sigma_sq = Make_M(V,Tau,C)

    try:
        preform_emcee(time, flux, err, ROW)
    except np.linalg.linalg.LinAlgError as err:
        continue
