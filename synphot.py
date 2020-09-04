import glob
import os

import matplotlib.pyplot as plt
import pysynphot as S # Importing it as S bad style, but used by the docs

if __name__ == '__main__':
    # Load new horizon's MVIC bandpass data
    fits_files = glob.glob(os.path.join('data/nh_filters', '*.FITS'))
    new_horizons_bandpass = [S.FileBandpass(f) for f in fits_files]
    # Load vega spectrum
    vega = S.Vega
    # Do a sample observation
    observation = S.Observation(vega, new_horizons_bandpass)
    # Plot the sample observation
    fig, ax = plt.subplots()
    ax.plot(obs.binwave, obs.binflux)
    plt.show()
