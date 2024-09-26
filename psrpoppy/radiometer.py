#!/usr/bin/python

import math


def calcFlux(snr,
             beta,
             Trec,
             Tsky,
             gain,
             n_p,
             t_obs,
             bw,
             duty):

    """Calculate flux assuming radiometer equation"""

    signal = signalterm(beta,
                        Trec,
                        Tsky,
                        gain,
                        n_p,
                        t_obs,
                        bw,
                        duty)

    return snr * signal


def calcSNR(flux,
            beta,
            Trec,
            Tsky,
            gain,
            n_p,
            t_obs,
            bw,
            duty):

    """Calculate the S/N ratio assuming radiometer equation"""

    signal = signalterm(beta,
                        Trec,
                        Tsky,
                        gain,
                        n_p,
                        t_obs,
                        bw,
                        duty)

    return (1.0/(1.0+0.0473*duty**(-0.627)))*flux / signal


def signalterm(beta,
               Trec,
               Tsky,
               gain,
               n_p,
               t_obs,
               bw,
               duty):

    """Returns the rest of the radiometer equation (aside from
        SNR/Flux"""

    dterm = duty / (1.0 - duty)
    return beta * (Trec+Tsky) * math.sqrt(dterm) \
        / gain / math.sqrt(n_p * t_obs * bw)
