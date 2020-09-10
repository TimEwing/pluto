# std library imports
import glob
import os

# third party imports
import matplotlib.pyplot as plt
import numpy as np
os.environ["PYSYN_CDBS"] = os.environ.get("PYSYN_CDBS", "data/pysyn/") # set pysyn env if not set
import pysynphot as S # Importing it as S bad style for python, but used by the docs
from astropy.io import fits

# custom module imports
import utils

def plot_bandpasses_vega(fig, ax, filter_names):
    # Load new horizon's MVIC bandpass data
    bandpasses = [utils.get_bandpass(x) for x in filter_names]

    # Load vega spectrum
    vega = S.Vega

    # TODO: Need to normalize

    # Plot observations
    for bandpass in bandpasses:
        # Do a synthetic observation
        observation = S.Observation(vega, bandpass)

        ## Convert to vegamag
        # In Carly's paper, there's a part where flam at Pluto is converted to flam at Earth.
        # We probably don't actually need to do that; difference is going to be something like 
        # (Pluto's orbital radius) / (distance to Vega), which is about 2.5e-5. 
        observation.convert('vegamag')

        # Plot observation
        ax.plot(
            observation.binwave, 
            observation.binflux, 
            label=bandpass.name,
            c=utils.COLORMAP.get(bandpass.name, None), # Get colors from colormap, default to None
            linewidth=1,
        )

        # Print countrates
        print(f"Countrate:{bandpass.name}, {observation.effstim()}")


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
    plot_bandpasses_vega(fig, ax, bandpasses)

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
