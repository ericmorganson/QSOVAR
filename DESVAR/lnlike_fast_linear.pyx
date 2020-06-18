""" Example of wrapping cos function from math.h using Cython. """

cdef extern from "math.h":
    double cos(double arg)
    double log(double arg)
    double exp(double arg)

def cos_func(arg):
    return cos(arg)

def sum_func(fluxes):
    output=0
    for num in range(fluxes.shape[0]):
        output+=fluxes[num]
    return output

def lnlike(theta, time, flux, flux_err_sq):
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
            exp_dt_Tau = exp(-(time[n] - time[n-1])/Tau_ten)
            # TODO: SPEED UP EVOLVESTATE
            state = evolvestate(state, exp_dt_Tau, 0, V_ten**2)
            lnp += lognorm(state, flux[n], flux_err_sq[n])
            state = weightedmean(state, flux[n], flux_err_sq[n])
        return lnp


def lognorm(state, flux, flux_err_sq):
        # This finds the normal distribution for any measurement (x),
        # mean(mu), and variance squared (var2).
        fmean, V_sq = state
        var_err = (V_sq + flux_err_sq)
        lognorm = -0.5*(((flux-fmean)*(flux-fmean))/var_err + log(6.28318530718*var_err))
        return lognorm


def evolvestate(state, exp_dt_Tau, mu, V_sq_old):
        fmean, V_sq = state
        fmean_new = exp_dt_Tau*(fmean - mu) + mu

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
