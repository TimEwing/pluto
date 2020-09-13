# std library imports
import glob
import os

# third party imports
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# custom module imports
from utils import S, Vega # import pysynphot and Vega from utils since they needs some setup
import utils



def plot_spectra(fig, ax, spectrum, filter_names):
    # Load new horizon's MVIC bandpass data
    bandpasses = [utils.get_bandpass(x) for x in filter_names]

    # TODO: Need to normalize

    # Plot observations
    for bandpass in bandpasses:
        # Do a synthetic observation
        observation = S.Observation(vega, bandpass)

        ## Convert to vegamag
        observation.convert('vegamag')

        # Plot observation
        ax.plot(
            observation.binwave, 
            observation.binflux, 
            label=bandpass.name,
            c=utils.COLORMAP.get(bandpass.name, None), # Get colors from colormap, default to None
            linewidth=1,
        )

def fill_between(fig, ax, filter_name_sets)


if __name__ == '__main__':
    # Set telescope aperture; otherwise, default for Hubble is used
    S.setref(area=176.7) # cm^2; using 75mm radius aperture from DOI: 10.1117/12.617901

    # Plot the sample observation
    fig, ax = plt.subplots()

    bandpasses = [
        'NH_RED',
        'NH_BLUE',
        # 'JOHNSON_R',
        # 'JOHNSON_B',
        'HST_F435W',
        'HST_F555W',
    ]
    vega = S.Vega

    plot_spectra(fig, ax, bandpasses)

    # Plot the source spectrum
    ax.plot(
        S.Vega.wave, 
        S.Vega.flux, 
        label="Vega", 
        c='black',
        linewidth=1,
    )

    # Setup plots
    ax.legend()
    ax.set_xlim(0,11000)
    ax.invert_yaxis() # To make vegamag units more intuitive
    plt.show()
