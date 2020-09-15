# std library imports
import glob
import os

# third party imports
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# custom module imports
from utils import S # import pysynphot from utils since it needs some setup
import utils



def plot_spectrum(fig, ax, spectrum, filter_names):
    # Load new horizon's MVIC bandpass data
    bandpasses = [utils.get_bandpass(x) for x in filter_names]

    # Plot observations
    for bandpass in bandpasses:
        spectrum.convert('angstrom')
        bandpass.convert('angstrom')

        # Do a synthetic observation
        observation = S.Observation(spectrum, bandpass)

        ## Convert to vegamag
        observation.convert('vegamag')

        # Plot observation
        ax.plot(
            bandpass.wave, 
            bandpass.throughput, 
            label=bandpass.name,
            c=utils.COLORMAP.get(bandpass.name, None), # Get colors from colormap, default to None
            linewidth=1,
        )

# This if statement runs only if the module is called from the command line
# If it gets imported, __name__ will not be set to '__main__'
if __name__ == '__main__':
    # Set telescope aperture; otherwise, default for Hubble is used
    S.setref(area=176.7) # cm^2; using 75mm radius aperture from DOI: 10.1117/12.617901
    spectrum = S.Vega

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

    plot_spectrum(fig, ax, spectrum, bandpasses)

    # Plot the source spectrum
    ax.plot(
        spectrum.wave, 
        spectrum.flux, 
        label="Vega", 
        c='black',
        linewidth=1,
    )

    # Setup plots
    ax.legend()
    ax.set_xlim(0,11000)
    # ax.invert_yaxis() # To make vegamag units more intuitive
    plt.show()
