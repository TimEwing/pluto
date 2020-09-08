# std library imports
import os
import glob
# third party imports
os.environ["PYSYN_CDBS"] = os.environ.get("PYSYN_CDBS", "data/pysyn/") # set pysyn env if not set
import pysynphot as S # Importing it as S bad style for python but used by the docs
from astropy.io import fits

# Constants
NH_FILTER_DIR = 'data/nh_filters'
COLORMAP = {
    'RED': 'red',
    'BLUE': 'blue',
    'NIR': 'orange',
    'CH4': 'green',
    'PAN_1': 'violet',
    'PAN_2': 'purple',
}


def get_nh_bandpass(filter_name, path=NH_FILTER_DIR):
    # Search for the file; filter_name should be RED, BLUE, NIR, CH4, PAN_1, or PAN_2
    file = glob.glob(os.path.join(path, f'{filter_name}*.FITS'))
    # If we found anything but 1 file, raise an error
    if not file:
        raise FileNotFoundError(f"No .FITS file found for {filter_name}")
    if not file:
        raise ValueError(f"More than one .FITS file found matching filter name")
    file = file[0]

    # Read fits file
    wavelengths, throughputs, header = get_nh_bandpass_data(file)

    # Same as pysynphot.spectrum.ArraySpectralElement
    bandpass = S.ArrayBandpass(
        name=filter_name,
        wave=wavelengths,
        throughput=throughputs,
        waveunits='Angstrom',
    )
    return bandpass


# The new horizons .FITS bandpass files aren't quite in the right format for pysynphot, so we need
# to load them by hand using astropy
def get_nh_bandpass_data(file):
    # Open file
    print(file)
    with fits.open(file) as data:
        # Parse out data and header
        header = data[0].header
        data = data[0].data
        # Remove low values of throughput so pysynphot can check overlaps better
        data = data[data[:,1] > 10e-7]
        # Scale data
        wavelengths, throughputs = data.transpose()
        wavelengths = wavelengths*1e4

        return wavelengths, throughputs, header
