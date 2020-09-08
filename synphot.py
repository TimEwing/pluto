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
    new_horizons_bandpasses = [utils.get_nh_bandpass(x) for x in filter_names]

    # Load vega spectrum
    vega = S.Vega

    # Plot observations
    for bandpass in new_horizons_bandpasses:
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
        print(f"{bandpass.name}: {observation.effstim()}")


if __name__ == '__main__':
    # Set telescope aperture; otherwise, default for Hubble is used
    S.setref(area=176.7) # cm^2; using 75mm radius aperture from DOI: 10.1117/12.617901

    # Plot the sample observation
    fig, ax = plt.subplots()
    plot_bandpasses_vega(fig, ax, ["RED", "BLUE", "NIR", "CH4", "PAN_1"])

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
