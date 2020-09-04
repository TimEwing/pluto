# std library imports
import os
import glob
import csv
# third party imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def plot_fits(fig, ax, files):
    for file in files:
        with fits.open(file) as data:
            # Parse out data and header
            header = data[0].header
            data = data[0].data
            # Scale data
            x, y = data.transpose()
            x = x*1e4

            # Plot data
            basename = os.path.basename(file)
            shortname = basename.split('_')[0]
            ax.plot(
                x,
                y,
                alpha=0.5,
                label=shortname,
            )

def plot_wfc(fig, ax, files):
    for file in files:
        with open(file) as open_file:
            reader = csv.reader(open_file, delimiter=' ')
            # Read off lines, casting everything to floats
            data = np.array([[float(x), float(y)] for x,y in reader])
            # Plot
            basename = os.path.basename(file)
            shortname, ext = os.path.splitext(basename)
            ax.plot(
                *data.transpose(), # Same as passing data.transpose()[0] then [1] seperately
                label=shortname,
            )

def plot_spectrum(fix, ax, files):
    for file in files:
        with open(file) as open_file:
            reader = csv.reader(
                open_file, 
                delimiter=' ',
                skipinitialspace=True, # Strip leading whitespace to avoid getting many columns
            )
            # Read off lines, casting everything to floats
            data = np.array([[float(x), float(y)] for x,y in reader])
            # Scale data
            data = data.transpose()
            data[0] *= 1e4
            # Plot
            basename = os.path.basename(file)
            shortname, ext = os.path.splitext(basename)
            shortname = shortname.split('_')[0]
            ax.plot(
                *data, # Same as passing data[0] then data[1] seperately
                label=shortname,
            )


if __name__ == '__main__':

    fig, axs = plt.subplots(1, 1)
    filter_ax = axs

    ## Plot filters
    # Plot new horizon's MVIC filter data
    fits_files = glob.glob(os.path.join('data/nh_filters', '*.FITS'))
    plot_fits(fig, filter_ax, fits_files)
    # Plot HST filter data
    wfc_files = glob.glob(os.path.join('data/hst_filters', 'wfc_*.dat'))
    plot_wfc(fig, filter_ax, wfc_files)

    ## Plot spectra
    # Save the range so we can set it later since the spectra is broader than the filters
    xlim = filter_ax.get_xlim()
    spectra_files = glob.glob(os.path.join('data/spectra', 'charon_spectrum.dat'))
    plot_spectrum(fig, filter_ax, spectra_files)

    # set x range
    xlim = filter_ax.set_xlim(xlim)
    filter_ax.legend()
    plt.show()