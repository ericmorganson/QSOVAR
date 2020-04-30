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
import itertools


if len(sys.argv) < 6:
  print(len(sys.argv))
  print("lightcurveplot.py INFITS ROWSTART ROWEND MCMC_TYPE NAME")
  print("INFITS is a fits table of DES Data,")
  print("ROWNUM is a row in that file")
  print("MCMC_TYPE is how you set up the MCMC search array: 'grid' or 'normal'")
  print("NAME is the additional identifying name for all output files.")
  sys.exit()
# Initial MCMC guesses
V =0.3
Tau = 365.0
dMu = 0.0
scale = 1 #aka, scaling is the same as r band
print(sys.argv[2],sys.argv[3],sys.argv[4])

def lognorm(state, flux, flux_err_sq):
        # This finds the normal distribution for any measurement (x),
        # mean(mu), and variance squared (var2).
        fmean, V_sq = state

        return -0.5*(((flux-fmean)**2)/(V_sq + flux_err_sq) + np.log(2*np.pi*(V_sq + flux_err_sq)))

def evolvestate(state, exp_dt_Tau, mu, V_sq_old):
        fmean, V_sq = state
        fmean_new = exp_dt_Tau*(fmean - mu)+ mu
        V_sq_new = V_sq_old*(1 - exp_dt_Tau*exp_dt_Tau) + (exp_dt_Tau*exp_dt_Tau)*V_sq
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
        color_dict = {"g":1, "r":2, "i":3, "z":4}
        lc_median = {}
        time = np.array([])
        flux_norm = np.array([])
        flux_err_norm = np.array([])
        array_org = np.array([])
        fig = plt.figure(figsize=(10,10))
        ax5 = fig.add_subplot(111)
        color_plt=iter(cm.rainbow(np.linspace(0,1,5)))
        for color in "griz":
            FITS = pyfit.open(args[1])

            #Obtain flux, err, and time arrays
            try:
                lc_flux = FITS[1].data['LC_FLUX_PSF_'+color][ROW]
            except IndexError:
                print('Error')
                exit()

            lc_median[color] = FITS[1].data['MEDIAN_PSF_'+color][ROW]
            lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
            lc_time = FITS[1].data['LC_MJD_'+color][ROW]
            time = np.append(time, lc_time[lc_time != 0.]) #remove the zeros

            limit = len(lc_time[lc_time!=0])
            array_org = np.append(array_org, color_dict[color]*np.ones(limit))
            lc_flux = lc_flux[:limit]
            lc_flux_err = lc_flux_err[:limit]
            col=next(color_plt)

            ax5.scatter(lc_time[lc_time!=0], (lc_flux[0:limit] - lc_median[color])/lc_median[color], label=color, c=np.array([col]))
            ax5.legend()
            ax5.set_title("Pre-correction light curve")


            flux_norm = np.append(flux_norm, (lc_flux - lc_median[color])/lc_median[color])
            flux_err_norm = np.append(flux_err_norm, lc_flux_err/lc_median[color])

        fig.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_"+"all_plot_mu_before.pdf")
        lc_time, lc_flux_norm, lc_flux_err_norm, lc_array_org = map(list, zip(*sorted(zip(time, flux_norm, flux_err_norm, array_org))))
        return lc_flux_norm, lc_flux_err_norm, lc_time, lc_median, lc_array_org, FITS

def lnlike(theta, time, flux_og, flux_err_sq, color_sort_ones):
        logV, logTau, dMu_g, dMu_r, dMu_i, dMu_z, scale_g, scale_i, scale_z = theta

        V_ten = 10**logV
        Tau_ten = 10**logTau
        color_dict = {"g":0, "r":1, "i":2, "z":3}
        dMu_dict = {"g":dMu_g, "r":dMu_r, "i":dMu_i, "z":dMu_z}
        scale_dict = {"g":scale_g, "r":1, "i":scale_i, "z":scale_z}
        flux = np.zeros_like(flux_og)
        for color in 'griz':
            flux += (flux_og*color_sort_ones[color_dict[color]]-dMu_dict[color])*scale_dict[color]
        # For multiband, this line will be  flux = (flux-dMu_griz)*scale_griz
        # dMu_griz will containt dmu_g, dmu_r, etc.
        # scale_griz will by one for r band and the others will be set relative to other bands
        #due to normalization, mu is by definition 0
        state=(0,V_ten**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0])
        state = weightedmean(state, flux[0], flux_err_sq[0])
        exp_dt_Tau_list = np.exp(-(np.subtract(time[1:], time[:-1])/(Tau_ten)))
        for n in range(1,len(flux)):
            if time[n]-time[n-1] < 0:
                print('AHHHHH, NEGATIVE TIME!!!')
            exp_dt_Tau = exp_dt_Tau_list[n-1]
            state=evolvestate(state, exp_dt_Tau, 0, V_ten**2) ### NEED TO WORK ON THIS
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        return lnp


def lnlike_old(theta, time, flux, flux_err_sq):
        logV, logTau = theta

        V_ten = 10**logV
        Tau_ten = 10**logTau

        #due to normalization, mu is by definition 0
        state=(0 ,V_ten**2)
        lnp = lognorm(state, flux[0], flux_err_sq[0])
        state = weightedmean(state, flux[0], flux_err_sq[0])
        for n in range(1,len(flux)):
            if time[n]-time[n-1] < 0:
                print('AHHHHH, NEGATIVE TIME!!!')
            exp_dt_Tau = np.exp(-(time[n] - time[n-1])/Tau_ten)
            state=evolvestate(state, exp_dt_Tau, 0, V_ten**2) ### NEED TO WORK ON THIS
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        return lnp


def lnprior(theta):
        logV, logTau, dMu_g, dMu_r, dMu_i, dMu_z, scale_g, scale_i, scale_z = theta
        if -3 < logV < 2 and 0 < logTau < 4 and list(x for x in [dMu_g, dMu_r, dMu_i, dMu_z] if -1.0 < x < 1.0) and list(y for y in [scale_g, scale_i, scale_z] if 0 < y <5):
            #print('yay!')
            #return 0.0
            return  0.5*logTau + 1*(-dMu_g**2/0.1**2 - dMu_r**2/0.1**2 - dMu_i**2/0.1**2 - dMu_z**2/0.1**2) #- logV
        #print(logV, logTau, dMu_g, dMu_r, dMu_i, dMu_z, scale_g, scale_i, scale_z)
        return -np.inf

def lnprior_step2(theta):
        logV, logTau= theta
        if -3 < logV < 2 and 0 < logTau < 4:
            return  0.5*logTau  #- logV
        #print(logV, logTau, dMu_g, dMu_r, dMu_i, dMu_z, scale_g, scale_i, scale_z)
        return -np.inf


def lnprob(theta, x, y, yerr, color_sort_ones):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        logprobs.append(lp + lnlike(theta, x, y, yerr, color_sort_ones))
        logvals.append(theta)
        return lp + lnlike(theta, x, y, yerr, color_sort_ones)

def lnprob_step2(theta, x, y, yerr):
        lp = lnprior_step2(theta)
        if not np.isfinite(lp):
            return -np.inf
        logprobs.append(lp + lnlike_old(theta, x, y, yerr))
        logvals.append(theta)
        return lp + lnlike_old(theta, x, y, yerr)

def sausageplot(Vari,time,delta_f,Tau,dt,sigma_sq, ROW, fig, dmu_scales, color_sort_ones, mu, err_sig):
        err_top = []
        err_bot = []
        Logpr = []
        Vari = 10 ** Vari
        Tau = 10 ** Tau
        times = []
        ax4 = fig.add_subplot(111)
        color_dict = {"g":0, "r":1, "i":2, "z":3}
        plot_dict = {"g":"og", "r":"or", "i":"ok", "z":"ob"}
        dMu_dict = {"g":dmu_scales[0], "r":dmu_scales[1], "i":dmu_scales[2], "z":dmu_scales[3]}
        scale_dict = {"g":dmu_scales[4], "r":1, "i":dmu_scales[5], "z":dmu_scales[6]}
        flux = np.zeros_like(delta_f)
        color_array = np.zeros_like(delta_f)
        mag_arr=np.array([])
        for color in 'griz':
            flux += (delta_f*color_sort_ones[color_dict[color]]-dMu_dict[color])*scale_dict[color]+ dMu_dict['r']
            color_array += color_dict[color]*color_sort_ones[color_dict[color]]

            mag = 22.5-2.5*np.log10(((delta_f*color_sort_ones[color_dict[color]]-dMu_dict[color])*scale_dict[color]+ dMu_dict['r'])*mu['r'] + mu['r'])
            t = time*color_sort_ones[color_dict[color]]
            e = err*color_sort_ones[color_dict[color]]
            [mag, e, t] =  goodrow(mag,e,t)
            mag_arr = np.append(mag_arr, mag)
            ax4.errorbar(t-57000, mag, yerr = e, fmt = plot_dict[color], label=color.upper()+', dmu='+str(round(dMu_dict[color], 5))+', scale='+str(round(scale_dict[color], 5)))

        [xlim,ylim] = boundaries(time-57000, np.array(mag_arr))
        ax4.set_ylim(ylim)
        ax4.set_xlabel('MJD-57000')
        ax4.set_ylabel('Mags')
        ax4.legend()

        delta_f = flux

        while np.count_nonzero((time[1:] - time[:-1])<0.004) > 0:
            t_df_sig = list(zip(time, delta_f, sigma_sq, color_array))
            t_new=[]
            df_new=[]
            sig_new=[]
            color_iter = itertools.cycle([0, 1, 2, 3])
            i = 0
            while i < len(t_df_sig)-1:
                col = next(color_iter)
                #print(col)
                #if t_df_sig[i+1][0]- t_df_sig[i][0] < 0.004:

                    #print(str(t_df_sig[i+1][0]- t_df_sig[i][0]) + "\t" + str(t_df_sig[i+1][0]) + "\t" + str(t_df_sig[i][0]))
                    #t_new.append(np.mean((t_df_sig[i][0], t_df_sig[i+1][0])))
                    #df_new.append(np.mean((t_df_sig[i][1], t_df_sig[i+1][1])))
                    #sig_new.append(np.mean((t_df_sig[i][2], t_df_sig[i+1][2])))
                    #i= i+2
                if int(t_df_sig[i][0]) in [int(t) for t in t_new]:
                    i=i+1
                elif col != t_df_sig[i][3]:
                    i=i+1
                else:
                    #print(str(col) + "\t"+ str(t_df_sig[i][3]))
                    #print(t_df_sig[i][0])
                    t_new.append(t_df_sig[i][0])
                    df_new.append(t_df_sig[i][1])
                    sig_new.append(t_df_sig[i][2])
                    i = i+1
            time = np.asarray(t_new)

            #for i in range(len(time)-1):
            #    if((time[i+1]-time[i])<0.004):
            #        print(time[i+1] - time[i])

            delta_f = np.asarray(df_new)
            sigma_sq = np.asarray(sig_new)

        for t in np.arange(min(time),max(time),dt):
            sigma_sq_s = np.append(sigma_sq, 0.0)
            dtime = np.append(time,t)
            delta_time = np.fabs(np.array(dtime,ndmin=2)-np.array(dtime,ndmin=2).T)
            #print(delta_time)
            #Make the modified matrix
            S = np.diag(sigma_sq_s) + (Vari**2) * np.exp(-1.0 * delta_time / Tau)

            #print("sigma_sq_s= "+str(sigma_sq_s))
            #print("Vari= " + str(Vari))
            #print("delta_time= "+ str(delta_time))
            #print("Tau= "+str(Tau))
            #print(np.exp(-1.0 * delta_time / Tau))
            times.append(t-57000)
            #Set up our guess for the flux
            t0 = [i for i in time if i >= t]
            t1 = [i for i in time if i <= t]

            t0 = t0[0]
            time_ind0 = np.where(time == t0)
            if len(list(time_ind0[0]))>1:
                time_ind0 = (np.asarray([time_ind0[0][0]]), )

            t1 = t1[-1]
            time_ind1 = np.where(time == t1)
            if len(list(time_ind1[0]))>1:
                time_ind1 = (np.asarray([time_ind1[0][0]]), )

            tp_0 = t0 / t
            tp_1 = t1 / t
            Fg1 = (tp_0 * delta_f[time_ind0]) + (tp_1 * delta_f[time_ind1])
            #print(mu)
            mu_r = mu['r']
            #print(mu_r)
            Fg1 = (Fg1 - mu_r)/mu_r
            Fg2 = Fg1 - 1
            Fg3 = Fg1 + 1

            delf1 = np.append(delta_f,Fg1)
            delf2 = np.append(delta_f,Fg2)
            delf3 = np.append(delta_f,Fg3)

            sign, value = np.linalg.slogdet(S)
            deter = sign * np.exp(value.item())
            #print(value)
            #print(deter)

            logP  = -.5*np.log((deter))-.5*((np.dot(delf1,(np.dot((delf1),np.linalg.inv(S))))))
            logP1 = -.5*np.log((deter))-.5*((np.dot(delf2,(np.dot((delf2),np.linalg.inv(S))))))
            logP2 = -.5*np.log((deter))-.5*((np.dot(delf3,(np.dot((delf3),np.linalg.inv(S))))))

            X = [Fg2,Fg1,Fg3]
            Y = [logP1,logP,logP2]
            X = np.array(X)
            X = (X.T)[0]
            Matr = np.array([[X[0]**2,X[0],1],[X[1]**2,X[1],1],[X[2]**2,X[2],1]])
            result = np.linalg.solve(Matr, Y)
            sig_sq = -1/(2*result[0])
            f_0 = -(result[1]/(2*result[0]))

            parabola = np.polyfit(X,Y,2)
            f = np.poly1d(parabola)

            x_new = np.linspace(X[0], X[-1], 5000)
            y_new = f(x_new)

            delf = delta_f
            Fm = op.fmin(lambda x: -f(x),0, disp=False)
            sig = 22.5-2.5*np.log10((sig_sq))
            logPn = result[0]*((Fm - f_0)**2) + result[2] - (result[1]**2)/(4*result[0])
            nsig=2
            sig1 = mu_r * (np.sqrt(sig_sq) + f_0) + mu_r
            sig2 = mu_r * (-np.sqrt(sig_sq) + f_0) + mu_r
            center = 22.5-2.5*np.log10((mu_r * (f_0))+mu_r)
            err_t = 22.5-2.5*np.log10((mu_r * (f_0 + nsig*np.sqrt(sig_sq)))+mu_r)
            err_b = 22.5-2.5*np.log10((mu_r * (f_0 - nsig*np.sqrt(sig_sq)))+mu_r)
            Logpr.append(center)
            err_top.append(err_t)
            err_bot.append(err_b)

        #plotdata(sys.argv, ROW, ax4, dmu_scales, color_sort_ones)
        #print(err_top)
        #print(Logpr)
        #print(err_bot)
        ax4.plot(times,err_top, color = 'g')
        ax4.plot(times,Logpr, color = 'b')
        ax4.plot(times,err_bot, color = 'g')
        ax4.set_title(' Object '+str(ROW)+ " Sausage Plot"+"\nV="+str(Vari)+" Tau="+str(Tau))

def sausageplot_step2(Vari,time,delta_f,Tau,dt,sigma_sq, dMu_dict, scale_dict, color_sort_ones, ROW, fig):
        err_top = []
        err_bot = []
        Logpr = []
        Vari = 10 ** Vari
        Tau = 10 ** Tau
        times = []

        ax4 = fig.add_subplot(111)
        color_dict = {"g":0, "r":1, "i":2, "z":3}
        plot_dict = {"g":"og", "r":"or", "i":"ok", "z":"ob"}
        #dMu_dict = {"g":dmu_scales[0], "r":dmu_scales[1], "i":dmu_scales[2], "z":dmu_scales[3]}
        #scale_dict = {"g":dmu_scales[4], "r":1, "i":dmu_scales[5], "z":dmu_scales[6]}
        #flux = np.zeros_like(delta_f)
        color_array = np.zeros_like(delta_f)
        mag_arr=np.array([])
        for color in 'griz':
            #flux += (delta_f*color_sort_ones[color_dict[color]]-dMu_dict[color])*scale_dict[color]+ dMu_dict['r']
            color_array += color_dict[color]*color_sort_ones[color_dict[color]]

            mag = 22.5-2.5*np.log10(delta_f*mu['r'] + mu['r'])
            t = time*color_sort_ones[color_dict[color]]
            e = np.sqrt(sigma_sq)*color_sort_ones[color_dict[color]]
            [mag, e, t] =  goodrow(mag,e,t)
            mag_arr = np.append(mag_arr, mag)
            ax4.errorbar(t-57000, mag, yerr = e, fmt = plot_dict[color], label=color.upper()+', dmu='+str(round(dMu_dict[color], 5))+', scale='+str(round(scale_dict[color], 5)))

        [xlim,ylim] = boundaries(time-57000, np.array(mag_arr))
        ax4.set_ylim(ylim)
        ax4.set_xlabel('MJD-57000')
        ax4.set_ylabel('Mags')
        ax4.set_title(' Object '+str(ROW)+ " Sausage Plot"+"\nV="+str(Vari)+" Tau="+str(Tau))
        ax4.legend()

        while np.count_nonzero((time[1:] - time[:-1])<0.004) > 0:
            t_df_sig = list(zip(time, delta_f, sigma_sq, color_array))
            t_new=[]
            df_new=[]
            sig_new=[]
            i = 0
            color_iter = itertools.cycle([0, 1, 2, 3])
            while i < len(t_df_sig)-1:
                #if t_df_sig[i+1][0]- t_df_sig[i][0] < 0.004:
                #    #print(t_df_sig[i+1][0]- t_df_sig[i][0])
                #    t_new.append(np.mean((t_df_sig[i][0], t_df_sig[i+1][0])))
                #    df_new.append(np.mean((t_df_sig[i][1], t_df_sig[i+1][1])))
                #    sig_new.append(np.mean((t_df_sig[i][2], t_df_sig[i+1][2])))
                #    i= i+2
                #else:
                #    t_new.append(t_df_sig[i][0])
                #    df_new.append(t_df_sig[i][1])
                #    sig_new.append(t_df_sig[i][2])
                #    i=i+1
                col = next(color_iter)
                #print(col)
                #if t_df_sig[i+1][0]- t_df_sig[i][0] < 0.004:

                    #print(str(t_df_sig[i+1][0]- t_df_sig[i][0]) + "\t" + str(t_df_sig[i+1][0]) + "\t" + str(t_df_sig[i][0]))
                    #t_new.append(np.mean((t_df_sig[i][0], t_df_sig[i+1][0])))
                    #df_new.append(np.mean((t_df_sig[i][1], t_df_sig[i+1][1])))
                    #sig_new.append(np.mean((t_df_sig[i][2], t_df_sig[i+1][2])))
                    #i= i+2
                if int(t_df_sig[i][0]) in [int(t) for t in t_new]:
                    i=i+1
                elif col != t_df_sig[i][3]:
                    i=i+1
                else:
                    #print(str(col) + "\t"+ str(t_df_sig[i][3]))
                    #print(t_df_sig[i][0])
                    t_new.append(t_df_sig[i][0])
                    df_new.append(t_df_sig[i][1])
                    sig_new.append(t_df_sig[i][2])
                    i = i+1
            time = np.asarray(t_new)
            delta_f = np.asarray(df_new)
            sigma_sq = np.asarray(sig_new)

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
            if len(list(time_ind0[0]))>1:
                time_ind0 = (np.asarray([time_ind0[0][0]]), )

            t1 = t1[-1]
            time_ind1 = np.where(time == t1)
            if len(list(time_ind1[0]))>1:
                time_ind1 = (np.asarray([time_ind1[0][0]]), )

            tp_0 = t0 / t
            tp_1 = t1 / t

            mu_r = mu['r']

            Fg1 = (tp_0 * delta_f[time_ind0]) + (tp_1 * delta_f[time_ind1])
            Fg1 = (Fg1 - mu_r)/mu_r
            Fg2 = Fg1 - 1
            Fg3 = Fg1 + 1

            delf1 = np.append(delta_f,Fg1)
            delf2 = np.append(delta_f,Fg2)
            delf3 = np.append(delta_f,Fg3)

            sign, value = np.linalg.slogdet(S)
            deter = sign * np.exp(value.item())

            logP  = -.5*np.log((deter))-.5*((np.dot(delf1,(np.dot((delf1),np.linalg.inv(S))))))
            logP1 = -.5*np.log((deter))-.5*((np.dot(delf2,(np.dot((delf2),np.linalg.inv(S))))))
            logP2 = -.5*np.log((deter))-.5*((np.dot(delf3,(np.dot((delf3),np.linalg.inv(S))))))

            X = [Fg2,Fg1,Fg3]
            Y = [logP1,logP,logP2]
            X = np.array(X)
            X = (X.T)[0]
            Matr = np.array([[X[0]**2,X[0],1],[X[1]**2,X[1],1],[X[2]**2,X[2],1]])
            result = np.linalg.solve(Matr, Y)
            sig_sq = -1/(2*result[0])
            f_0 = -(result[1]/(2*result[0]))

            parabola = np.polyfit(X,Y,2)
            f = np.poly1d(parabola)

            x_new = np.linspace(X[0], X[-1], 5000)
            y_new = f(x_new)

            delf = delta_f
            Fm = op.fmin(lambda x: -f(x),0, disp=False)
            sig = 22.5-2.5*np.log10((sig_sq))
            logPn = result[0]*((Fm - f_0)**2) + result[2] - (result[1]**2)/(4*result[0])
            nsig=2
            sig1 = mu_r * (np.sqrt(sig_sq) + f_0) + mu_r
            sig2 = mu_r * (-np.sqrt(sig_sq) + f_0) + mu_r
            center = 22.5-2.5*np.log10((mu_r * (f_0))+mu_r)
            err_t = 22.5-2.5*np.log10((mu_r * (f_0 + nsig*np.sqrt(sig_sq)))+mu_r)
            err_b = 22.5-2.5*np.log10((mu_r * (f_0 - nsig*np.sqrt(sig_sq)))+mu_r)
            Logpr.append(center)
            err_top.append(err_t)
            err_bot.append(err_b)

        #ax4 = fig.add_subplot(111)
        #ax4 = fig.add_subplot(2, 2, 4)
        #plotdata_step2(sys.argv, ROW, ax4)
        ax4.plot(times,err_top, color = 'g')
        ax4.plot(times,Logpr, color = 'b')
        ax4.plot(times,err_bot, color = 'g')
        ax4.set_title(' Object '+str(ROW)+ " Sausage Plot"+"\nV="+str(Vari)+" Tau="+str(Tau))
        #fig.savefig("figure/"+str(ROW)+sys.argv[5]+"_"+sys.argv[4]+"_sausageplot_step2.pdf")


def lnprob_dens(theta, x, y, yerr):
        logprobs_dens.append(lnlike_old(theta, x, y, yerr))
        logvals_dens.append(theta)
        return logprobs_dens

def goodrow(mags,errs,mjds):
        good = (mags > 0) & (mags != 22.5) & (mags != np.inf) & ~np.isnan(mags) & (errs > 0) & (errs != np.inf) & ~np.isnan(errs) & (mjds > 0) & (mjds != np.inf) & ~np.isnan(mjds)
        return [mags[good], errs[good], mjds[good]]

def plotdata_step2(args, ROW, ax4):
        fits = pyfit.open(args[1])[1].data
        VarRows = ROW
        [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits, ROW) #17999
        rownum = str(ROW)

        mean_g = np.mean(mags_g)
        mean_r = np.mean(mags_r)
        mean_i = np.mean(mags_i)
        mean_z = np.mean(mags_z)

        diff_r = mean_g - mean_r
        diff_i = mean_g - mean_i
        diff_z = mean_g - mean_z

        #ax4.rcParams['font.size'] = 18

        ax4.errorbar(mjds_g-57000, mags_g, yerr = errs_g, fmt = 'og', label = "g-band")
        ax4.errorbar(mjds_r-57000, mags_r, yerr = errs_r, fmt = 'or', label = "r-band")
        ax4.errorbar(mjds_i-57000, mags_i, yerr = errs_i, fmt = 'ok', label = "i-band")
        ax4.errorbar(mjds_z-57000, mags_z, yerr = errs_z, fmt = 'ob', label = "z-band")
        [xlim,ylim] = boundaries(np.hstack([mjds_g,mjds_r,mjds_i,mjds_z]),np.hstack([mags_g,mags_r,mags_i,mags_z]))# #

        ax4.set_ylim(ylim)
        ax4.set_xlabel('MJD-57000')
        ax4.set_ylabel('Mags')
        ax4.legend()
        field = args[1].split('_')[0]
        ax4.set_title(field+' Object '+rownum)


def plotdata(args, ROW, ax4, dmu_scales, color_sort_ones):
        fits = pyfit.open(args[1])[1].data
        VarRows = ROW
        [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z] = getdata(fits, ROW) #17999
        rownum = str(ROW)

        mean_g = np.median(mags_g)
        mean_r = np.median(mags_r)
        mean_i = np.median(mags_i)
        mean_z = np.median(mags_z)

        mags_g = (mags_g - mean_g)/mean_g
        mags_r = (mags_r - mean_r)/mean_r
        mags_i = (mags_i - mean_i)/mean_i
        mags_z = (mags_z - mean_z)/mean_z
        #dmu_scales[0:4] = 22.5 - 2.5*np.log10(dmu_scales)[0:4]

        dMu_dict = {"g":dmu_scales[0], "r":dmu_scales[1], "i":dmu_scales[2], "z":dmu_scales[3]}
        scale_dict = {"g":dmu_scales[4], "r":1, "i":dmu_scales[5], "z":dmu_scales[6]}

        mags_g = 22.5-2.5*np.log10(((mags_g-dMu_dict['g'])*scale_dict['g'] + dMu_dict['r'])*mean_r + mean_r) #(mags_g-(dMu_dict['g']) - mean_g)*scale_dict['g'] +mean_r + dMu_dict['r']
        mags_r = 22.5-2.5*np.log10(((mags_r-dMu_dict['r'])*scale_dict['r'] + dMu_dict['r'])*mean_r + mean_r) #(mags_r-(dMu_dict['r']) - mean_r)*scale_dict['r'] +mean_r + dMu_dict['r']
        mags_i = 22.5-2.5*np.log10(((mags_i-dMu_dict['i'])*scale_dict['i'] + dMu_dict['r'])*mean_r + mean_r) #(mags_i-(dMu_dict['i']) - mean_i)*scale_dict['i'] +mean_r + dMu_dict['r']
        mags_z = 22.5-2.5*np.log10(((mags_z-dMu_dict['z'])*scale_dict['z'] + dMu_dict['r'])*mean_r + mean_r) #(mags_z-(dMu_dict['z']) - mean_z)*scale_dict['z'] +mean_r + dMu_dict['r']

        #ax4.rcParams['font.size'] = 18

        ax4.errorbar(mjds_g-57000, mags_g, yerr = errs_g, fmt = 'og', label='G, dmu='+str(round(dMu_dict['g'], 5))+', scale='+str(round(scale_dict['g'], 5)))
        ax4.errorbar(mjds_r-57000, mags_r, yerr = errs_r, fmt = 'or', label='R, dmu='+str(round(dMu_dict['r'], 5))+', scale='+str(round(scale_dict['r'], 5)))
        ax4.errorbar(mjds_i-57000, mags_i, yerr = errs_i, fmt = 'ok', label='I, dmu='+str(round(dMu_dict['i'], 5))+', scale='+str(round(scale_dict['i'], 5)))
        ax4.errorbar(mjds_z-57000, mags_z, yerr = errs_z, fmt = 'ob', label='Z, dmu='+str(round(dMu_dict['z'], 5))+', scale='+str(round(scale_dict['z'], 5)))
        [xlim,ylim] = boundaries(np.hstack([mjds_g,mjds_r,mjds_i,mjds_z]),np.hstack([mags_g,mags_r,mags_i,mags_z]))# #

        ax4.set_ylim(ylim)
        ax4.set_xlabel('MJD-57000')
        ax4.set_ylabel('Mags')
        ax4.legend()
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

        [mags_g, errs_g, mjds_g] = goodrow(fits[num]['LC_FLUX_PSF_G'], fits[num]['LC_FLUXERR_PSF_G']/fits[num]['LC_FLUX_PSF_G'], fits[num]['LC_MJD_G']) #22.5-2.5*np.log10(fits[num]['LC_FLUX_PSF_G'])
        [mags_r, errs_r, mjds_r] = goodrow(fits[num]['LC_FLUX_PSF_R'], fits[num]['LC_FLUXERR_PSF_R']/fits[num]['LC_FLUX_PSF_R'], fits[num]['LC_MJD_R'])
        [mags_i, errs_i, mjds_i] = goodrow(fits[num]['LC_FLUX_PSF_I'], fits[num]['LC_FLUXERR_PSF_I']/fits[num]['LC_FLUX_PSF_I'], fits[num]['LC_MJD_I'])
        [mags_z, errs_z, mjds_z] = goodrow(fits[num]['LC_FLUX_PSF_Z'], fits[num]['LC_FLUXERR_PSF_Z']/fits[num]['LC_FLUX_PSF_Z'], fits[num]['LC_MJD_Z'])
        return [mags_g, errs_g, mjds_g, mags_r, errs_r, mjds_r, mags_i, errs_i, mjds_i, mags_z, errs_z, mjds_z]

def perform_emcee(time, flux, sigma_sq, color_sort_ones, ROW, mu):
        diff_time = [x - time[i - 1] for i, x in enumerate(time)][1:]
        fig = plt.figure(figsize=(10,10))

        nll = lambda *args: -lnlike(*args)
        ndim, nwalkers = 9, 100
        #MAKE POSITION ARRAY ARRAY FOR WALKERS
        if sys.argv[4].lower() == 'normal':
            result = [np.log10(V), np.log10(Tau), dMu, dMu, dMu, dMu, scale, scale, scale]
            pos = (np.random.rand(nwalkers,ndim)-0.5)*np.array([1, 1, 0.1, 0.1, 0.1, 0.1, .5, .5, .5])+result
        elif sys.argv[4].lower() == 'optimal':
            result = op.minimize(nll, [np.log10(V), np.log10(Tau), dMu, dMu, dMu, dMu, scale, scale, scale],args=(time,flux, err**2, color_sort_ones)) #not sure how this will work for the multi-band situation??
            pos = [result['x'] + 1e-1*np.random.rand(ndim) for i in range(nwalkers)]
        else:
            print("What the hell do you want to do?")
            print("'optimal', or 'normal' search through MCMC?")
            exit()

        #print(pos[:,1].reshape((-1)))

        #run sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time, flux, err**2, color_sort_ones))
        # run mcmc
        sampler.run_mcmc(pos, 300)
        samples_pick = sampler.chain[:,:,:].reshape((-1, ndim))
        logprobs_pick_2d = sampler.lnprobability[:,:]
        logprobs_pick = sampler.lnprobability[:,:].reshape((-1))

        logprobs_samp_2d = sampler.lnprobability[:,50:]
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        logprobs_samp = sampler.lnprobability[:,50:].reshape((-1))

        #make corner plot

        max_theta_old = logvals[np.argmax(logprobs)] #old way with non-sampler values; for optimal/normal differentiation

        max_theta = samples[np.argmax(logprobs_samp)]

        print(np.unravel_index(logprobs_pick_2d.argmax(), logprobs_pick_2d.shape))
        print(np.argmax(logprobs_pick))
        print(logprobs_pick[np.argmax(logprobs_pick)])
        print(samples_pick[np.argmax(logprobs_pick)])

        print(np.unravel_index(logprobs_samp_2d.argmax(), logprobs_samp_2d.shape))
        print(np.argmax(logprobs_samp))
        print(logprobs_samp[np.argmax(logprobs_samp)])

        print(max_theta)


        fig1 = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$",r"$d\mu_g$", r"$d\mu_r$", r"$d\mu_i$", r"$d\mu_z$", r"scale$_g$", r"scale$_i$", r"scale$_z$"],
                            truths=[max_theta[0], max_theta[1], max_theta[2], max_theta[3], max_theta[4], max_theta[5], max_theta[6], max_theta[7], max_theta[8]])


        fig1.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_"+"triangle_np_mu_"+str(sys.argv[5])+ ".pdf")

        #PRINT MAX THETA VALUES TO THE SCREEN
        print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]), 'dMu:', str(max_theta[2]))
        #dt = 5; currently 2 below :/
        sausageplot(max_theta[0], time, flux, max_theta[1], 5, err**2, ROW, fig, max_theta[2:], color_sort_ones, mu, sigma_sq)
        fig.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_mu_sausage"+ str(sys.argv[5])+ ".pdf")

        fig3=plt.figure(figsize=(10,10))
        ax6 = fig3.add_subplot(111)
        for i in range(nwalkers):
            ax6.plot(sampler.chain[i, 50:, 0].reshape((-1)))
        ax6.set_title("MCMC V values")   #log(probability)")
        fig3.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_V"+ str(sys.argv[5])+ ".pdf")

        fig4=plt.figure(figsize=(10,10))
        ax7 = fig4.add_subplot(111)
        for i in range(nwalkers):
            ax7.plot(sampler.chain[i, 50:, 1].reshape((-1)))
        ax7.set_title("MCMC tau values")   #log(probability)")
        fig4.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_tau"+ str(sys.argv[5])+ ".pdf")
        plt.close("all")

        #WRITE THE FOUND MAX THETA VALUES TO FILE
        filename ='scratch_new/'+ str(ROW) + sys.argv[4]+"_all_band_"+'object_dMu_' + str(sys.argv[5]) +'.txt'
        with open(filename, 'w+') as fout:
            fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')
            fout.write('Object: ' + str(ROW)+ ' ' + 'dMu: ' + str(max_theta[2])+ '\n')

        return max_theta

def perform_emcee_step2(time, flux, sigma_sq, dMu_dict, scale_dict, color_sort_ones, ROW, mu):
        diff_time = [x - time[i - 1] for i, x in enumerate(time)][1:]
        fig = plt.figure(figsize=(10,10))

        nll = lambda *args: -lnlike_old(*args)
        ndim, nwalkers = 2, 100
        #MAKE POSITION ARRAY ARRAY FOR WALKERS
        if sys.argv[4].lower() == 'normal':
            result = [np.log10(V), np.log10(Tau)]
            pos = (np.random.rand(nwalkers,ndim)-0.5)*np.array([1, 1])+result
        elif sys.argv[4].lower() == 'optimal':
            result = op.minimize(nll, [np.log10(V), np.log10(Tau)],args=(time,flux, err**2)) #not sure how this will work for the multi-band situation??
            pos = [result['x'] + 1e-1*np.random.rand(ndim) for i in range(nwalkers)]
        else:
            print("What the hell do you want to do?")
            print("'optimal', or 'normal' search through MCMC?")
            exit()

        #print(pos[:,1].reshape((-1)))

        #run sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_step2, args=(time, flux, err**2))
        # run mcmc
        sampler.run_mcmc(pos, 300)

        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        logprobs_samp = sampler.lnprobability[:,50:].reshape((-1))

        #make corner plot
        max_theta = samples[np.argmax(logprobs_samp)]

        print(max_theta)

        fig1 = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$"],
                            truths=[max_theta[0], max_theta[1]])


        fig1.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_"+"triangle_step2_linear_"+str(sys.argv[5])+ ".pdf")

        #PRINT MAX THETA VALUES TO THE SCREEN
        print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]))
        #dt = 5; currently 2 below :/
        sausageplot_step2(max_theta[0], time, flux, max_theta[1], 5, err**2, dMu_dict, scale_dict, color_sort_ones, ROW, fig)
        fig.savefig("figure/"+str(ROW)+sys.argv[4]+"_all_band_"+"sausage_step2_linear_"+ str(sys.argv[5])+ ".pdf")

        plt.close("all")

        #WRITE THE FOUND MAX THETA VALUES TO FILE
        filename ='scratch_new/'+ str(ROW) + sys.argv[4]+"_all_band_step2_linear_"+'object_dMu_' + str(sys.argv[5]) +'.txt'
        with open(filename, 'w+') as fout:
            fout.write('Object: ' + str(ROW)+ ' ' + 'Tau: ' + str(max_theta[1])+' ' + 'V: '+ str(max_theta[0]) + '\n')


for ROW in range(int(sys.argv[2]),int(sys.argv[3])):
    logprobs = []
    logvals = []

    logprobs_dens = []
    logvals_dens = []

    print("Running object "+str(ROW))
    flux, err, time, mu, color_sort,  FITS = get_vals(sys.argv,ROW)
    print(min(np.array(time[1:]) - np.array(time[:-1])))

    #DOESN'T MAKE SENSE TO LOOK AT ROWS WITH NO FLUX MEASUREMENTS
    if len(flux) == 0:
        print("Flux length is zero")
        continue

    #ONLY LET POSITIVE FLUXES AND ERRORS THROUGH
    #NOTE: FLUXES HERE ARE NORMALIZED, SO IF FLUX = 0, THEN norm(FLUX) = -1
    zip_tfe = zip(time, flux, err, color_sort)
    filter_tfe = [(t, f, e, col) for t, f, e, col in zip_tfe if f > -1 and e > 0]
    time, flux, err, color_sort = zip(*filter_tfe)
    time = np.array(time)
    flux = np.array(flux)
    err = np.array(err)
    color_sort = np.array(color_sort)
    unique, count = np.unique(color_sort, return_counts=True)
    color_sort_dict = {1:"g", 2:"r", 3:"i", 4:"z"}
    color_counts = dict(zip([color_sort_dict[u] for u in unique], count))
    print(color_counts)

    #ONLY LOOK AT BRIGHT OBJECTS (WITHOUT OVERSATURATION)
    #dim = [i for i in mu if 22.5-2.5*np.log10(mu[i])>21]
    #dim_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])>21]
    #bright = [i for i in mu if 22.5-2.5*np.log10(mu[i])<16]
    #bright_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])<16]
    #if dim:
    #    print("Row is too dim in band: "+ str(dim))
    #    print(dim_val)
    #    plt.close("all")
    #    continue
    #if bright:
    #    print("Row is HELLA bright in band: "+str(bright))
    #    print(bright_val)
    #    plt.close("all")
    #    continue

    print(mu)


    color_sort_ones = [(color_sort == 1).astype(int)]
    for num in range(2,5):
        color_sort_ones = np.concatenate((color_sort_ones,[(color_sort == num).astype(int)]), axis=0)
    np.set_printoptions(threshold=np.inf)

    int_time = time.astype(int)
    flux_dict = {"g":[], "r":0, "i":[], "z":[]}
    g_flux = []
    i_flux = []
    z_flux = []
    prev_time = int_time[0]
    for i in range(0,len(int_time)):
        if abs(int_time[i] - prev_time) <= 1:
            if color_sort_dict[color_sort[i]] == 'r':
                flux_dict["r"] = flux[i]
            else:
                flux_dict[color_sort_dict[color_sort[i]]].append(flux[i])
        else:
            #record results
            #print(prev_time, flux_dict)
            if flux_dict["r"] != 0:
                f_r = flux_dict["r"]
                if len(flux_dict["g"]) !=0:
                    for f_g in flux_dict["g"]:
                        g_flux.append([f_r, f_g])
                if len(flux_dict["i"]) !=0:
                    for f_i in flux_dict["i"]:
                        i_flux.append([f_r, f_i])
                if len(flux_dict["z"]) !=0:
                    #print(len(flux_dict["z"]))
                    for f_z in flux_dict["z"]:
                        z_flux.append([f_r, f_z])
            #start new count
            flux_dict = {"g":[], "r":0, "i":[], "z":[]}
            prev_time = int_time[i]
            if color_sort_dict[color_sort[i]] == 'r':
                flux_dict["r"] = flux[i]
            else:
                flux_dict[color_sort_dict[color_sort[i]]].append(flux[i])

    fig_lin, axs = plt.subplots(3)
    fig_lin.suptitle('r vs flux scatter')

    if g_flux:
        g_flux = np.array(g_flux)
        p_g, res_g = np.polyfit(g_flux[:,0], g_flux[:,1] - g_flux[:,0], 1, cov=True)
        axs[0].scatter(g_flux[:,0], g_flux[:,1]- g_flux[:,0])
        axs[0].plot(g_flux[:,0], g_flux[:,0]*p_g[0] +p_g[1], label="{0:.3f}x + {1:.3f}".format(*p_g))
    else:
        p_g = (1, 0)
    axs[0].set_xlabel('r')
    axs[0].set_ylabel('g - r')
    axs[0].legend()
    #axs[0].set_title('g vs. r')

    #print(i_flux)
    if i_flux:
        i_flux = np.array(i_flux)
        p_i, res_i= np.polyfit(i_flux[:,0], i_flux[:,1] - i_flux[:,0], 1, cov=True)
        axs[1].scatter(i_flux[:,0], i_flux[:,1]- i_flux[:,0])
        axs[1].plot(i_flux[:,0], i_flux[:,0]*p_i[0] +p_i[1], label="{0:.3f}x + {1:.3f}".format(*p_i))
    else:
        p_i = (1, 0)
    axs[1].legend()
    axs[1].set_xlabel('r')
    axs[1].set_ylabel('i - r')
    #axs[1].set_title('i vs. r')

    if z_flux:
        z_flux = np.array(z_flux)
        p_z, res_z= np.polyfit(z_flux[:,0], z_flux[:,1] - z_flux[:,0], 1, cov=True)
        axs[2].scatter(z_flux[:,0], z_flux[:,1]- z_flux[:,0])
        axs[2].plot(z_flux[:,0], z_flux[:,0]*p_z[0] +p_z[1], label="{0:.3f}x + {1:.3f}".format(*p_z))
    else:
        p_z = (1, 0)
    axs[2].legend()
    axs[2].set_xlabel('r')
    axs[2].set_ylabel('z - r')
    #axs[2].set_title('z vs. r')

    slope = [p_g[0]+1, 1, p_i[0]+1, p_z[0]+1]
    slope_err = [np.sqrt(res_g[0][0]), 0, np.sqrt(res_i[0][0]), np.sqrt(res_z[0][0])]
    std_col = []
    for col in 'GRIZ':
        sd = FITS[1].data['STD_PSF_'+col][ROW]
        print("initial STD: "+str(sd))
        arr = np.array(FITS[1].data['LC_FLUX_PSF_'+col][ROW])
        final_list = [x for x in arr if (x > np.mean(arr) - 2 * sd)]
        final_list = [x for x in final_list if (x < np.mean(arr) + 2 * sd)]
        print('final STD: '+str(np.std(np.array(final_list))))

        #std_col.append(0.5*(np.percentile(flux_color, 84) - np.percentile(flux_color, 16)))
        #print(np.percentile(flux_color, 84))
        #print(np.percentile(flux_color, 16))
    intercept = [p_g[1], 0, p_i[1], p_z[1]]


    for i in range(len(slope)):
        if slope_err[i]*2>slope[i]:
            print(color_sort_dict[i+1], ROW)
            print(slope[i])
            slope[i] = std_col[i]/std_col[1]
            print(slope[i])


    print("Slope: "+ str(slope))
    print("intercepts: "+str(intercept))
    print("Slope Error: "+str(slope_err))
    print("StDev: "+ str(std_col))

    fig_lin.savefig("figure/"+str(ROW)+sys.argv[4]+sys.argv[5]+"_linear_scatter.pdf")


    try:
        #print("Trying optimal initialization")
        #sys.argv[4] = 'optimal'
        #max_theta = perform_emcee(time, flux, err, color_sort_ones, ROW, mu)

        color_dict = {"g":0, "r":1, "i":2, "z":3}
        dMu_dict = {"g":intercept[0], "r":intercept[1], "i":intercept[2], "z":intercept[3]}
        scale_dict = {"g":slope[0], "r":slope[1], "i":slope[2], "z":slope[3]}
        flux_mod = np.zeros_like(flux)
        for color in 'griz':
            flux_mod += ((flux-dMu_dict[color])*color_sort_ones[color_dict[color]])/scale_dict[color]+ dMu_dict['r']

        perform_emcee_step2(time, flux_mod, err, dMu_dict, scale_dict, color_sort_ones, ROW, mu)

    #except ValueError:
    #    print("Falling back to normal initialization")
    #    sys.argv[4] = 'normal'
    #    perform_emcee(time, flux, err, color_sort, ROW, mu)
    except np.linalg.linalg.LinAlgError as err:
        continue
