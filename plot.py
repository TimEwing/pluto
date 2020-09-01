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
                # c='black',
                label=shortname,
            )


if __name__ == '__main__':

    fig, ax = plt.subplots()
    path = 'data/'

    fits_files = glob.glob(os.path.join(path, '*.FITS'))
    plot_fits(fig, ax, fits_files)

    # Make sure wfc_files is a list
    wfc_files = glob.glob(os.path.join(path, 'wfc_*.dat'))
    plot_wfc(fig, ax, wfc_files)

    ax.legend()
    plt.show()