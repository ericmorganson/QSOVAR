#!/usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfit
import matplotlib
matplotlib.use('Agg')
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
    print("MCMC_TYPE is how you set up the MCMC search ")
    print("array: 'optimal' or 'normal'")
    print("NAME is the additional identifying name for all output files.")
    sys.exit()
# Initial MCMC guesses
file_path = "fig_ssh/" #"fig_ssh/" #"figure/"
V = 0.3
Tau = 365.0
dMu = 0.0
scale = 1   # aka, scaling is the same as r band
print(sys.argv[2], sys.argv[3], sys.argv[4])


def lognorm(state, flux, flux_err_sq):
    # This finds the normal distribution for any measurement (x),
    # mean(mu), and variance squared (var2).
    fmean, V_sq = state
    var_err = (V_sq + flux_err_sq)
    lognorm = -0.5*(((flux-fmean)**2)/var_err + np.log(2*np.pi*var_err))
    return lognorm


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
    color_dict = {"g": 1, "r": 2, "i": 3, "z": 4}
    lc_median = {}
    time = np.array([])
    flux_norm = np.array([])
    err_norm = np.array([])
    array_org = np.array([])
    fig = plt.figure(figsize=(10, 10))
    ax5 = fig.add_subplot(111)
    color_plt = iter(cm.rainbow(np.linspace(0, 1, 5)))
    for color in "griz":
        FITS = pyfit.open(args[1])
        # Obtain flux, err, and time arrays
        try:
            lc_flux = FITS[1].data['LC_FLUX_PSF_'+color][ROW]
        except IndexError:
            print('Error')
            exit()

        lc_median[color] = FITS[1].data['MEDIAN_PSF_'+color][ROW]
        lc_flux_err = FITS[1].data['LC_FLUXERR_PSF_'+color][ROW]
        lc_time = FITS[1].data['LC_MJD_'+color][ROW]
        time = np.append(time, lc_time[lc_time != 0.])  # remove the zeros

        limit = len(lc_time[lc_time != 0])
        array_org = np.append(array_org, color_dict[color]*np.ones(limit))
        lc_flux = lc_flux[:limit]
        lc_flux_err = lc_flux_err[:limit]
        col = next(color_plt)

        ax5.scatter(lc_time[lc_time != 0],
                    (lc_flux[0:limit] - lc_median[color])/lc_median[color],
                    label=color, c=np.array([col]))
        ax5.legend()
        ax5.set_title("Pre-correction light curve")
        normed_flux = (lc_flux - lc_median[color])/lc_median[color]

        flux_norm = np.append(flux_norm, normed_flux)
        normed_err = lc_flux_err/lc_median[color]
        err_norm = np.append(err_norm, normed_err)
    # TODO: label plot with more descriptive headers
    fig.savefig(file_path+str(ROW)+sys.argv[4]+"_all_band_scatter_before.pdf")
    time, flux_norm, err_norm, array_org = map(list,
                                               zip(*sorted(zip(time,
                                                               flux_norm,
                                                               err_norm,
                                                               array_org))))
    return flux_norm, err_norm, time, lc_median, array_org, FITS


def lnlike_old(theta, time, flux, flux_err_sq):
    logV, logTau = theta

    V_ten = 10**logV
    Tau_ten = 10**logTau

    # due to normalization, mu is by definition 0
    state = (0, V_ten**2)
    lnp = lognorm(state, flux[0], flux_err_sq[0])
    state = weightedmean(state, flux[0], flux_err_sq[0])
    for n in range(1, len(flux)):
        if time[n]-time[n-1] < 0:
            print('AHHHHH, NEGATIVE TIME!!!')
        exp_dt_Tau = np.exp(-(time[n] - time[n-1])/Tau_ten)
        # TODO: SPEED UP EVOLVESTATE
        state = evolvestate(state, exp_dt_Tau, 0, V_ten**2)
        lnp += lognorm(state, flux[n], flux_err_sq[n])
        state = weightedmean(state, flux[n], flux_err_sq[n])
    return lnp


def lnprior_step2(theta):
    logV, logTau = theta
    if -3 < logV < 2 and 0 < logTau < 4:
        return 0.5*logTau
    return -np.inf


def lnprob_step2(theta, x, y, yerr):
    lp = lnprior_step2(theta)
    if not np.isfinite(lp):
        return -np.inf
    logprobs.append(lp + lnlike_old(theta, x, y, yerr))
    logvals.append(theta)
    return lp + lnlike_old(theta, x, y, yerr)


def sausageplot_step2(Vari, time, delta_f, Tau, dt, sigma_sq, dMu_dict,
                      scale_dict, color_sort_ones, ROW, fig):
    err_top = []
    err_bot = []
    Logpr = []
    Vari = 10 ** Vari
    Tau = 10 ** Tau
    times = []

    ax4 = fig.add_subplot(111)
    color_dict = {"g": 0, "r": 1, "i": 2, "z": 3}
    plot_dict = {"g": "og", "r": "or", "i": "ok", "z": "ob"}

    color_array = np.zeros_like(delta_f)
    mag_arr = np.array([])
    for color in 'griz':
        color_array += color_dict[color]*color_sort_ones[color_dict[color]]
        mag = 22.5-2.5*np.log10(delta_f*mu['r'] + mu['r'])
        t = time*color_sort_ones[color_dict[color]]
        e = np.sqrt(sigma_sq)*color_sort_ones[color_dict[color]]
        [mag, e, t] = goodrow(mag, e, t)
        mag_arr = np.append(mag_arr, mag)
        ax4.errorbar(t-57000, mag, yerr=e, fmt=plot_dict[color],
                     label=color.upper()+', dmu='+str(round(dMu_dict[color], 5))+', scale='+str(round(scale_dict[color], 5)))

    [xlim, ylim] = boundaries(time-57000, np.array(mag_arr))
    ax4.set_ylim(ylim)
    ax4.set_xlabel('MJD-57000')
    ax4.set_ylabel('Mags')
    ax4.set_title('Row ' + str(ROW) + "\nV=" + str(Vari) + " Tau=" + str(Tau))
    handles, labels = ax4.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    ax4.legend(handles, labels)

    while np.count_nonzero((time[1:] - time[:-1]) < 0.004) > 0:
        t_df_sig = list(zip(time, delta_f, sigma_sq, color_array))
        t_new = []
        df_new = []
        sig_new = []
        i = 0
        color_iter = itertools.cycle([0, 1, 2, 3])
        while i < len(t_df_sig)-1:
            col = next(color_iter)
            if int(t_df_sig[i][0]) in [int(t) for t in t_new]:
                i = i+1
            elif col != t_df_sig[i][3]:
                i = i+1
            else:
                t_new.append(t_df_sig[i][0])
                df_new.append(t_df_sig[i][1])
                sig_new.append(t_df_sig[i][2])
                i = i+1
        time = np.asarray(t_new)
        delta_f = np.asarray(df_new)
        sigma_sq = np.asarray(sig_new)

    for t in np.arange(min(time), max(time), dt):
        sigma_sq_s = np.append(sigma_sq, 0.0)
        dtime = np.append(time, t)
        delta_time = np.fabs(np.array(dtime, ndmin=2) - np.array(dtime, ndmin=2).T)
        # Make the modified matrix
        S = np.diag(sigma_sq_s) + (Vari**2) * np.exp(-1.0 * delta_time / Tau)
        times.append(t-57000)
        # Set up our guess for the flux
        t0 = [i for i in time if i >= t]
        t1 = [i for i in time if i <= t]

        t0 = t0[0]
        time_ind0 = np.where(time == t0)
        if len(list(time_ind0[0])) > 1:
            time_ind0 = (np.asarray([time_ind0[0][0]]), )

        t1 = t1[-1]
        time_ind1 = np.where(time == t1)
        if len(list(time_ind1[0])) > 1:
            time_ind1 = (np.asarray([time_ind1[0][0]]), )

        tp_0 = t0 / t
        tp_1 = t1 / t

        mu_r = mu['r']

        Fg1 = (tp_0 * delta_f[time_ind0]) + (tp_1 * delta_f[time_ind1])
        Fg1 = (Fg1 - mu_r)/mu_r
        Fg2 = Fg1 - 1
        Fg3 = Fg1 + 1

        delf1 = np.append(delta_f, Fg1)
        delf2 = np.append(delta_f, Fg2)
        delf3 = np.append(delta_f, Fg3)

        sign, value = np.linalg.slogdet(S)
        deter = sign * np.exp(value.item())

        logP = -.5*np.log((deter))-.5*((np.dot(delf1, (np.dot((delf1), np.linalg.inv(S))))))
        logP1 = -.5*np.log((deter))-.5*((np.dot(delf2, (np.dot((delf2), np.linalg.inv(S))))))
        logP2 = -.5*np.log((deter))-.5*((np.dot(delf3, (np.dot((delf3), np.linalg.inv(S))))))

        X = [Fg2, Fg1, Fg3]
        Y = [logP1, logP, logP2]
        X = np.array(X)
        X = (X.T)[0]
        Matr = np.array([[X[0]**2, X[0], 1],
                         [X[1]**2, X[1], 1],
                         [X[2]**2, X[2], 1]])
        result = np.linalg.solve(Matr, Y)
        sig_sq = -1/(2*result[0])
        f_0 = -(result[1]/(2*result[0]))

        parabola = np.polyfit(X, Y, 2)
        f = np.poly1d(parabola)

        x_new = np.linspace(X[0], X[-1], 5000)
        y_new = f(x_new)

        delf = delta_f
        Fm = op.fmin(lambda x: -f(x), 0, disp=False)
        sig = 22.5-2.5*np.log10((sig_sq))
        logPn = result[0]*((Fm - f_0)**2) + result[2] - (result[1]**2)/(4*result[0])
        nsig = 2
        sig1 = mu_r * (np.sqrt(sig_sq) + f_0) + mu_r
        sig2 = mu_r * (-np.sqrt(sig_sq) + f_0) + mu_r
        center = 22.5-2.5*np.log10((mu_r * (f_0))+mu_r)
        err_t = 22.5-2.5*np.log10((mu_r * (f_0 + nsig*np.sqrt(sig_sq)))+mu_r)
        err_b = 22.5-2.5*np.log10((mu_r * (f_0 - nsig*np.sqrt(sig_sq)))+mu_r)
        Logpr.append(center)
        err_top.append(err_t)
        err_bot.append(err_b)

    ax4.plot(times, err_top, color='g')
    ax4.plot(times, Logpr, color='b')
    ax4.plot(times, err_bot, color='g')
    # TODO: TRUNCATE VARI AND TAU VALUES USING SIGFIGS
    title_text = 'Object ' + str(ROW) + " Sausage Plot \n"
    subtitle_text = "V=" + str(Vari) + " Tau=" + str(Tau)
    ax4.set_title(title_text + subtitle_text)


def goodrow(mags, errs, mjds):
    mag_crit = (mags > 0) & (mags != 22.5) & (mags != np.inf) & ~np.isnan(mags)
    err_crit = (errs > 0) & (errs != np.inf) & ~np.isnan(errs)
    time_crit = (mjds > 0) & (mjds != np.inf) & ~np.isnan(mjds)
    good = mag_crit & err_crit & time_crit
    return [mags[good], errs[good], mjds[good]]


def boundaries(mjds, mags):
    mjdmin = 100*np.floor(np.min(mjds)*.01)
    mjdmax = 100*np.ceil(np.max(mags)*.01)
    [magmin, magmax] = np.percentile(mags, [2, 98])
    magmin = 0.5*np.floor(2.*(magmin-0.1))
    magmax = 0.5*np.ceil(2.*(magmax+0.1))
    return [[mjdmin, mjdmax], [magmax, magmin]]


def perform_emcee_step2(time, flux, sigma_sq, dMu_dict, scale_dict,
                        color_sort_ones, ROW, mu):
    diff_time = [x - time[i - 1] for i, x in enumerate(time)][1:]
    fig = plt.figure(figsize=(10, 10))

    nll = lambda *args: -lnlike_old(*args)
    ndim, nwalkers = 2, 100
    # MAKE POSITION ARRAY ARRAY FOR WALKERS
    if sys.argv[4].lower() == 'normal':
        result = [np.log10(V), np.log10(Tau)]
        pos = (np.random.rand(nwalkers, ndim)-0.5)*np.array([1, 1])+result
    elif sys.argv[4].lower() == 'optimal':
        # not sure how this will work for the multi-band situation??
        result = op.minimize(nll, [np.log10(V), np.log10(Tau)],
                             args=(time, flux, err**2))
        pos = [result['x'] + 1e-1*np.random.rand(ndim) for i in range(nwalkers)]
    else:
        print("What the hell do you want to do?")
        print("'optimal', or 'normal' search through MCMC?")
        exit()

    # run sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_step2, args=(time, flux, err**2))
    # run mcmc
    sampler.run_mcmc(pos, 200)

    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    logprobs_samp = sampler.lnprobability[:, 50:].reshape((-1))

    # make corner plot
    max_theta = samples[np.argmax(logprobs_samp)]
    fig1 = corner.corner(samples, labels=[r"log$_{10}V$", r"log$_{10}\tau$"],
                         truths=[max_theta[0], max_theta[1]])
    fig1.savefig(file_path + str(ROW) + sys.argv[4] + "_all_band_" + "triangle_linear_" + str(sys.argv[5]) + ".pdf")

    # PRINT MAX THETA VALUES TO THE SCREEN
    print('ROW:', ROW, 'Tau:', str(max_theta[1]), 'V:', str(max_theta[0]))

    sausageplot_step2(max_theta[0], time, flux, max_theta[1], 5, err**2, dMu_dict, scale_dict, color_sort_ones, ROW, fig)
    fig.savefig(file_path + str(ROW) + sys.argv[4] + "_all_band_" + "sausage_step2_linear_" + str(sys.argv[5]) + ".pdf")

    plt.close("all")

    # WRITE THE FOUND MAX THETA VALUES TO FILE
    filename = 'scratch_new/' + str(ROW) + sys.argv[4] + "_all_band_step2_linear_" + 'object_dMu_' + str(sys.argv[5]) + '.txt'
    with open(filename, 'w+') as fout:
        fout.write('Object: ' + str(ROW) + ' ' + 'Tau: ' + str(max_theta[1]) + ' ' + 'V: ' + str(max_theta[0]) + '\n')


def lin_color_sort(time, color_sort_dict, color_sort):
    int_time = time.astype(int)
    flux_dict = {"g": [], "r": 0, "i": [], "z": []}
    g_flux = []
    i_flux = []
    z_flux = []
    prev_time = int_time[0]
    for i in range(0, len(int_time)):
        if abs(int_time[i] - prev_time) <= 1:
            if color_sort_dict[color_sort[i]] == 'r':
                flux_dict["r"] = flux[i]
            else:
                flux_dict[color_sort_dict[color_sort[i]]].append(flux[i])
        else:
            # record results
            if flux_dict["r"] != 0:
                f_r = flux_dict["r"]
                if len(flux_dict["g"]) != 0:
                    for f_g in flux_dict["g"]:
                        g_flux.append([f_r, f_g])
                if len(flux_dict["i"]) != 0:
                    for f_i in flux_dict["i"]:
                        i_flux.append([f_r, f_i])
                if len(flux_dict["z"]) != 0:
                    for f_z in flux_dict["z"]:
                        z_flux.append([f_r, f_z])
            # start new count
            flux_dict = {"g": [], "r": 0, "i": [], "z": []}
            prev_time = int_time[i]
            if color_sort_dict[color_sort[i]] == 'r':
                flux_dict["r"] = flux[i]
            else:
                flux_dict[color_sort_dict[color_sort[i]]].append(flux[i])
    return g_flux, i_flux, z_flux


def plot_lin(flux, ax, y_label):
    if flux and len(flux) > 2:
        flux = np.array(flux)

        def linear(x, a, b):
            return a * x + b

        p, res = op.curve_fit(linear, flux[:, 0], flux[:, 1] - flux[:, 0])

        ax.scatter(flux[:, 0], flux[:, 1] - flux[:, 0])
        ax.plot(flux[:, 0], flux[:, 0]*p[0] + p[1], label="{0:.3f}x + {1:.3f}".format(*p))
    else:
        p = (0, 0)
        res = np.zeros((2, 2))

    ax.legend()
    ax.set_xlabel('r')
    ax.set_ylabel(y_label)

    return p, res


for ROW in range(int(sys.argv[2]), int(sys.argv[3])):
    logprobs = []
    logvals = []

    print("Running object "+str(ROW))
    flux, err, time, mu, color_sort,  FITS = get_vals(sys.argv, ROW)
    print(min(np.array(time[1:]) - np.array(time[:-1])))

    # DOESN'T MAKE SENSE TO LOOK AT ROWS WITH NO FLUX MEASUREMENTS
    if len(flux) == 0:
        print("Flux length is zero")
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
    print(color_counts)

    # ONLY LOOK AT BRIGHT OBJECTS (WITHOUT OVERSATURATION)
    # dim = [i for i in mu if 22.5-2.5*np.log10(mu[i])>21]
    # dim_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])>21]
    # bright = [i for i in mu if 22.5-2.5*np.log10(mu[i])<16]
    # bright_val = [22.5-2.5*np.log10(mu[i]) for i in mu if 22.5-2.5*np.log10(mu[i])<16]
    # if dim:
    #    print("Row is too dim in band: "+ str(dim))
    #    print(dim_val)
    #    plt.close("all")
    #    continue
    # if bright:
    #    print("Row is HELLA bright in band: "+str(bright))
    #    print(bright_val)
    #    plt.close("all")
    #    continue

    print(mu)

    color_sort_ones = [(color_sort == 1).astype(int)]
    for num in range(2, 5):
        color_sort_ones = np.concatenate((color_sort_ones, [(color_sort == num).astype(int)]), axis=0)
    np.set_printoptions(threshold=np.inf)

    g_flux, i_flux, z_flux = lin_color_sort(time, color_sort_dict, color_sort)

    fig_lin, axs = plt.subplots(3)
    fig_lin.suptitle('r vs color')

    p_g, res_g = plot_lin(g_flux, axs[0], "g-r")
    p_i, res_i = plot_lin(i_flux, axs[1], "i-r")
    p_z, res_z = plot_lin(z_flux, axs[2], "z-r")

    slope = [p_g[0]+1, 1, p_i[0]+1, p_z[0]+1]
    slope_err = [np.sqrt(res_g[0][0]), 0, np.sqrt(res_i[0][0]), np.sqrt(res_z[0][0])]
    int_err = [np.sqrt(res_g[1][1]), 0, np.sqrt(res_i[1][1]), np.sqrt(res_z[1][1])]
    std_col = []
    for col in 'GRIZ':
        std_col.append(FITS[1].data['STD_PSF_'+col][ROW])
    intercept = [p_g[1], 0, p_i[1], p_z[1]]

    for i in range(len(slope)):
        if np.count_nonzero(color_sort_ones[i]) > 1:
            var_crit = np.sum((std_col[i] - err*color_sort_ones[i])/np.sqrt(np.count_nonzero(color_sort_ones[i])))
        else:
            var_crit = 0
        print("Variable??")
        print(var_crit)
        if var_crit > 6:
            print("Variable Source!")
            # if slope_err[i]*2>slope[i]:
            #    slope[i] = std_col[i]/std_col[1]
        else:  # if nonvariable
            print("Non-Variable source in " + str(color_sort_dict[i+1]))
            if slope_err[i]*2 > slope[i]:
                slope[i] = 1

    print("Slope:       " + str(slope))
    print("Slope Error: " + str(slope_err))
    print("intercepts:  " + str(intercept))
    print("Int Error:   " + str(int_err))
    print("StDev:       " + str(std_col))

    fig_lin.savefig(file_path + str(ROW)+sys.argv[4] + "_" + sys.argv[5]+"_linear_scatter.pdf")

    try:
        color_dict = {"g": 0, "r": 1, "i": 2, "z": 3}
        dMu_dict = {"g": intercept[0], "r": intercept[1], "i": intercept[2], "z": intercept[3]}
        dMu_dict_err = {"g": int_err[0], "r": int_err[1], "i": int_err[2], "z": int_err[3]}
        scale_dict = {"g": slope[0], "r": slope[1], "i": slope[2], "z": slope[3]}
        flux_mod = np.zeros_like(flux)
        err_mod = np.zeros_like(flux)
        for color in 'griz':
            flux_mod += ((flux-dMu_dict[color])*color_sort_ones[color_dict[color]])/scale_dict[color] + dMu_dict['r']
            err_mod += np.sqrt(err**2 + dMu_dict_err[color]**2 + (flux-dMu_dict[color])**2/scale_dict[color]**2)/scale_dict[color]*color_sort_ones[color_dict[color]]

        perform_emcee_step2(time, flux_mod, err_mod, dMu_dict, scale_dict, color_sort_ones, ROW, mu)

    except np.linalg.linalg.LinAlgError as err:
        continue
