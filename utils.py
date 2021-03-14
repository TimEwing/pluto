# std library imports
import os
import glob
import csv
# third party imports
import numpy as np
from astropy.io import fits
from scipy.io.idl import readsav
# pysynphot
os.environ["PYSYN_CDBS"] = os.environ.get("PYSYN_CDBS", "data/pysyn/") # set pysyn env if not set
import pysynphot as S # Importing it as S bad style for python, but used by the docs

# Set telescope aperture; otherwise, default for Hubble is used
S.setref(area=176.7) # cm^2; using 75mm radius aperture from DOI: 10.1117/12.617901

# Constants
REWRITE_DATA_FILE = False # Set to true to make plots write to output.csv

NH_FILTER_DIR = 'data/nh_filters'
HST_FILTER_DIR = 'data/hst_filters'
CHARON_SPECTRUM = 'data/spectra/charon_spectrum.dat'
PLUTO_SPECTRUM = 'data/spectra/pluto_spectrum.dat'
HD_SPECTRUM = 'data/hd205905_spectrum_withheader.txt'
OUTPUT_DATA_FILE = 'data/output.csv'
# vegafile not used since pysynphot comes with a vega spectrum
# VEGAFILE = 'data/vega_spectrum.fits'
NH_OBSERVATION_FILE = 'data/nh_obs.idlsave'

# For plots
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
    'HST_F555W': 'tomato',
    'HST_F435W': 'cyan',
}
LABELMAP = {
    'pluto': 'Pluto',
    'charon': 'Charon',
    'hd': 'HD205905',
    'vega': 'Vega',
    'JOHNSON_U': 'Johnson U',
    'JOHNSON_B': 'Johnson B',
    'JOHNSON_V': 'Johnson V',
    'JOHNSON_R': 'Johnson R',
    'JOHNSON_G': 'Johnson G',
    'NH_RED': 'Red',
    'NH_BLUE': 'Blue',
    'NH_NIR': 'Nir',
    'NH_CH4': 'Ch4',
    'NH_PAN_1': 'Pan 1',
    'NH_PAN_2': 'Pan 2',
    'HST_F555W': 'F555W',
    'HST_F435W': 'F435W',
}

NH_PIVOT_WAVELENGTH = { # in Angstroms
    'NH_BLUE': 4920.00,
    'NH_RED': 6240.00,
    'NH_NIR': 8610.00,
    'NH_CH4': 8830.00,
}

# New Horizons Adjustment Factors
NH_AF = {
    'NH_BLUE': {0: 1.02, 1: 0.99}, 
    'NH_RED': {0: 1.23, 1: 1.21}, 
    'NH_NIR': {0: 1.38, 1: 1.27}, 
    'NH_CH4': {0: 1.46, 1: 1.44}, 
}


def write_to_output(name, data):
    if REWRITE_DATA_FILE:
        with open(OUTPUT_DATA_FILE, 'a') as out:
            out.write(f"{name},{','.join(data)}\n")


def get_observations(target, *filter_names):
    if target == 'hd':
        return get_hd_observations(*filter_names)
    if target == 'charon' or target == 'pluto':
        return get_nh_observations(target, *filter_names)

def get_hd_observations(*filter_names):

    output_dict = {
        'obs_to_target': 1,
        'sun_to_target': 1,
        'lat': -1,
        'lon': -1,
        'met': 0,
    }

    filter_map = {
        'NH_RED': 'red',
        'NH_BLUE': 'blue',
        'NH_CH4': 'ch4',
        'NH_NIR': 'nir',
    }

    gain = 58.5 # e-/count
    exposure_time = 1.10784 # seconds, from mc0_0330933531_0x536_sci_2.fit

    data_map = {
        "blue": 188.6,
        "red": 527.8,
        "nir": 347.5,
        "ch4": 97.0,
    }

    p_map = {
        'psolar_blue': 2.024E+13,
        'psolar_ch4': 2.703E+13,
        'psolar_nir': 1.105E+14,
        'psolar_red': 7.880E+13,
    }

    for filter_name in filter_names:
        filter_str = filter_map[filter_name]
        flux = data_map[filter_str] / (exposure_time * p_map[f'psolar_{filter_str}'])
        # raise NotImplementedError("What side is the HD observation?")
        side = 1

        output_dict.update({
            f'{filter_name}_counts': data_map[filter_str],
            f'{filter_name}_calib_flux': flux * NH_AF[filter_name][side],
            f'{filter_name}_pivot_wavelength': NH_PIVOT_WAVELENGTH[filter_name],
        })
    # Return output as a length-1 list so it works with other code
    return [output_dict]

def get_nh_observations(
        target, *filter_names, 
        file=NH_OBSERVATION_FILE,):
    save = readsav(file)

    if target == 'pluto':
        prefix = ''
    elif target == 'charon':
        prefix = 'charon_'

    filter_map = {
        'NH_RED': 'red',
        'NH_BLUE': 'blue',
        'NH_CH4': 'ch4',
        'NH_NIR': 'nir',
    }

    output_dict = {
        'obs_to_target': save['sc_range'],
        'sun_to_target': save['targ_sun'],
        'lon': save[f'{target}_lon'],
        'lat': save[f'{target}_lat'],
        'phase': save[f'phase_{target}'],
        'met': save['met'].tolist(), # tolist casts int32s to ints for json
    }

    for filter_name in filter_names:
        filter_str = filter_map[filter_name]

        expected_length = len(save['sc_range'])
        adjustment_factors = np.array([NH_AF[filter_name][side] for side in save['side']])

        output_dict.update({
            f'{filter_name}_counts': save[f'{prefix}{filter_str}_counts'] * adjustment_factors,
            f'{filter_name}_calib_flux': save[f'calib_{prefix}{filter_str}_flux'] * adjustment_factors,
            f'{filter_name}_pivot_wavelength': [NH_PIVOT_WAVELENGTH[filter_name]] * expected_length,
        })

    # Return the transpose of the output_dict, i.e. a list of observations
    return transpose_dict(output_dict, output_dict.keys())


def transpose_dict(input_dict, fields):
    ## Given a dict and a list of fields, return a list of dicts. I.e.
    #   inputs: input_dict={'a': [1,2,3], 'b': [4,5,6]}, fields=['a','b']
    #   output: [{'a': 1, 'b', 4}, {'a': 2, 'b', 5}, {'a': 3, 'b', 6}]

    # Zip lists together
    zipped_lists = zip(*[input_dict[field] for field in fields]) # * operator passes lists as args
    # Convert each line of the zipped lists to a dict
    output_list = [dict(zip(fields, line)) for line in zipped_lists]
    return output_list


def get_bandpass(filter_name):
    # Be case insensitive
    filter_name = filter_name.upper()

    # New Horizons filters
    if filter_name.startswith('NH_'):
        bandpass = get_nh_bandpass(filter_name)
    # Johnson filters
    elif filter_name.startswith('JOHNSON_'):
        bandpass = get_johnson_bandpass(filter_name)
    # HST Filters
    elif filter_name.startswith('HST_'):
        bandpass = get_hst_bandpass(filter_name)
    else:
        raise ValueError(f"Filter not found: {filter_name}")

    bandpass.convert('Angstrom')
    return bandpass


def get_johnson_bandpass(filter_name):
    filter_color = filter_name[8:] # Strip JOHNSON_ prefix
    filter_color = filter_color.lower()

    bandpass = S.ObsBandpass(f'johnson,{filter_color.lower()}')
    bandpass.name = filter_name
    return bandpass


def get_hst_bandpass(filter_name, path=HST_FILTER_DIR):
    # Search for the file; filter_name should be F435W or F555W
    filter_color = filter_name[4:] # Remove HST_ prefix
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
    filter_color = filter_name[3:] # Remove NH_ prefix
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
        waveunits='micron',
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
        # data = data[data[:,1] > 10e-7]
        wavelengths, throughputs = data.transpose()

        return wavelengths, throughputs, header

def get_spectrum(target):
    if target == 'charon':
        return get_charon_spectrum()
    elif target == 'pluto':
        return get_pluto_spectrum()
    elif target == 'hd':
        return get_hd_spectrum()
    elif target == 'vega':
        return S.Vega
    raise ValueError(f"Unknown spectrum name: {target}")

def get_charon_spectrum():
    wavelengths, fluxes = get_charon_spectrum_data()
    spectrum = S.ArraySpectrum(
        wave=wavelengths,
        flux=fluxes,
        name='Charon',
        waveunits='micron',
    )
    return spectrum


def get_charon_spectrum_data():
    with open(CHARON_SPECTRUM) as file:
        reader = csv.reader(
            file, 
            delimiter=' ',
            skipinitialspace=True, # Strip leading whitespace to avoid getting many columns
        )
        # Read off lines, casting everything to floats
        data = np.array([[float(x), float(y)] for x,y in reader])
        wavelengths, fluxes = data.transpose()
        return wavelengths, fluxes


def get_pluto_spectrum():
    wavelengths, fluxes = get_pluto_spectrum_data()
    spectrum = S.ArraySpectrum(
        wave=wavelengths,
        flux=fluxes,
        name='Pluto',
        waveunits='micron',
    )
    return spectrum


def get_pluto_spectrum_data():
    with open(PLUTO_SPECTRUM) as file:
        reader = csv.reader(
            file, 
            delimiter=' ',
            skipinitialspace=True, # Strip leading whitespace to avoid getting many columns
        )
        # Read off lines, casting everything to floats
        data = np.array([[float(x), float(y)] for x,y in reader])
        wavelengths, fluxes = data.transpose()
        return wavelengths, fluxes


def get_hd_spectrum():
    wavelengths, fluxes = get_hd_spectrum_data()
    spectrum = S.ArraySpectrum(
        wave=wavelengths,
        flux=fluxes,
        name='h205905',
        waveunits='Angstrom',
        fluxunits='Flam'
    )
    return spectrum


def get_hd_spectrum_data(filename=HD_SPECTRUM):
    with open(filename) as file:
        reader = csv.reader(
            file, 
            delimiter=' ',
            skipinitialspace=True, # Strip leading whitespace to avoid getting many columns
        )
        next(reader) # Skip first line (header)
        
        # Read off lines, casting everything to floats
        data = np.array([[float(x), float(y)] for x,y in reader])
        wavelengths, fluxes = data.transpose()
        return wavelengths, fluxes