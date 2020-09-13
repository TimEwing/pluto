# std library imports
import os
import glob
import csv
# third party imports
import numpy as np
from astropy.io import fits
# pysynphot
os.environ["PYSYN_CDBS"] = os.environ.get("PYSYN_CDBS", "data/pysyn/") # set pysyn env if not set
import pysynphot as S # Importing it as S bad style for python, but used by the docs
# Normalize Vega to HST: https://arxiv.org/pdf/1403.6861.pdf



# Constants
NH_FILTER_DIR = 'data/nh_filters'
HST_FILTER_DIR = 'data/hst_filters'
CHARON_SPECTRUM = 'data/spectra/charon_spectrum.dat'
VEGAFILE = 'data/vega_spectrum.fits'




COLORMAP = {
    'JOHNSON_U': 'violet',
    'JOHNSON_B': 'darkblue',
    'JOHNSON_V': 'yellow',
    'JOHNSON_R': 'darkred',
    'JOHNSON_G': 'darkgreen',
    'NH_RED': 'red',
    'NH_BLUE': 'blue',
    'NH_NIR': 'orange',
    'NH_CH4': 'green',
    'NH_PAN_1': 'hotpink',
    'NH_PAN_2': 'purple',
    'HST_F555W': 'orange',
    'HST_F435W': 'cyan',
}


def get_bandpass(filter_name):
    # Be case insensitive
    filter_name = filter_name.upper()

    # New Horizons filters
    if filter_name.startswith('NH_'):
        return get_nh_bandpass(filter_name)

    if filter_name.startswith('JOHNSON_'):
        return get_johnson_bandpass(filter_name)

    if filter_name.startswith('HST_'):
        return get_hst_bandpass(filter_name)

    raise ValueError(f"Filter not found: {filter_name}")


def get_johnson_bandpass(filter_name):
    filter_color = filter_name.lstrip('JOHNSON_')
    filter_color = filter_color.lower()

    bandpass = S.ObsBandpass(f'johnson,{filter_color.lower()}')
    bandpass.name = filter_name
    return bandpass


def get_hst_bandpass(filter_name, path=HST_FILTER_DIR):
    # Search for the file; filter_name should be F435W or F555W
    filter_color = filter_name.lstrip('HST_')
    file = glob.glob(os.path.join(path, f'wfc_{filter_color}.dat'))
    # If we found anything but 1 file, raise an error
    if not file:
        raise FileNotFoundError(f"No file found for {file}")
    if not file:
        raise ValueError(f"More than one .FITS file found matching filter name {filter_name}")
    file = file[0]

    # Read fits file
    wavelengths, throughputs = get_hst_bandpass_data(file)

    # Same as pysynphot.spectrum.ArraySpectralElement
    bandpass = S.ArrayBandpass(
        name=filter_name,
        wave=wavelengths,
        throughput=throughputs,
        waveunits='Angstrom',
    )
    return bandpass

def get_hst_bandpass_data(file):
    with open(file) as open_file:
        reader = csv.reader(open_file, delimiter=' ')
        # Read off lines, casting everything to floats
        data = np.array([[float(x), float(y)] for x,y in reader])
        wavelengths, throughputs = data.transpose()

        return wavelengths, throughputs
    

def get_nh_bandpass(filter_name, path=NH_FILTER_DIR):
    # Search for the file; filter_name should be NH_RED, NH_BLUE, NH_NIR, etc.
    filter_color = filter_name.lstrip('NH_')
    file = glob.glob(os.path.join(path, f'{filter_color}*.FITS'))
    # If we found anything but 1 file, raise an error
    if not file:
        raise FileNotFoundError(f"No .FITS file found for {filter_name}")
    if not file:
        raise ValueError(f"More than one .FITS file found matching filter name {filter_name}")
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

def get_charon_spectrum():
    # No idea what the units of the spectrum are. TODO: ask Cathy
    wavelengths, fluxes = get_charon_spectrum_data()
    spectrum = S.ArraySpectrum(
        wave=wavelengths,
        flux=fluxes,
        name='Charon',
        waveunits='Angstrom',
        fluxunits='???'
    )


def get_charon_spectrum_data():
    with open(CHARON_SPECTRUM) as file:
        reader = csv.reader(
            file, 
            delimiter=' ',
            skipinitialspace=True, # Strip leading whitespace to avoid getting many columns
        )
        # Read off lines, casting everything to floats
        data = np.array([[float(x), float(y)] for x,y in reader])
        # Scale data
        wavelengths, fluxes = data.transpose()
        wavelengths = wavelengths*1e4
        return wavelengths, fluxes