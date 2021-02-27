# third party imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# custom module imports
from utils import S
import utils


# 38.5 AU in km, as corrected in the Buie paper
BUIE_OBS_TO_TARGET = 5.760E9 
# 39.5 AU in km, as corrected in the Buie paper
BUIE_TARGET_TO_SUN = 5.909E9



def get_correlation_factors(
        spectrum, 
        input_bandpass_names, 
        output_bandpass_name):

    output_bandpass = utils.get_bandpass(output_bandpass_name)

    # Get the flux from observing as if you held one filter in front of the other
    correlation_factors = {}
    for input_bandpass_name in input_bandpass_names:
        input_bandpass = utils.get_bandpass(input_bandpass_name)
        # Perform synthetic observations
        observation = S.Observation(spectrum, input_bandpass, force='taper')
        observation = S.Observation(observation, output_bandpass, force='taper')

        # Can't be vegamag since magnitudes don't add linearly; needs to be flam
        try:
            correlation_factors[input_bandpass_name] = observation.effstim('flam')
        except ValueError:
            # 0 effstim
            correlation_factors[input_bandpass_name] = 0

    # Normalize factors
    factor_sum = sum(correlation_factors.values())
    correlation_factors = {k: v/factor_sum for k, v in correlation_factors.items()}

    return correlation_factors




# Convert from one single filter
# Distances are in the same units as the new horizons file, usually km
def plot_single_convert(
        fig, ax,
        target_name, # not fully implemented yet
        input_filter_name, output_filter_name, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, 
        output_target_to_sun=BUIE_TARGET_TO_SUN,
        x_axis='lon', y_axis='vegamag', **kwargs):
    # get recorded observations from new horizons
    nh_observations = utils.get_observations(target_name, input_filter_name)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_spectrum(target_name)

        converted_obs = single_convert(
            spectrum=spectrum,
            output_bandpass_name=output_filter_name,
            pivot_wavelength=nh_observation[f'{input_filter_name}_pivot_wavelength'],
            calib_flux=nh_observation[f'{input_filter_name}_calib_flux'],
            input_obs_to_target=nh_observation['obs_to_target'],
            input_target_to_sun=nh_observation['sun_to_target'],
            output_obs_to_target=output_obs_to_target,
            output_target_to_sun=output_target_to_sun,
        )


        plot_x.append(nh_observation[x_axis])
        plot_y.append(float(converted_obs.effstim(y_axis)))

    # Write out as csv to file for tables
    utils.write_to_output(
        name='.'.join([target_name, input_filter_name, output_filter_name, y_axis]),
        data=[f"{x}:{y}" for x,y in zip([x['met'] for x in nh_observations], plot_y)],
    )
    ax.scatter(plot_x, plot_y, **kwargs)


def multi_convert(
        spectrum, 
        input_bandpass_names, output_bandpass_name,
        pivot_wavelengths, calib_fluxes, 
        input_obs_to_target, input_target_to_sun, 
        output_obs_to_target, output_target_to_sun):
    # Get correlation factors to adjust for multi-filter overlap
    correlation_factors = get_correlation_factors(
        spectrum, 
        input_bandpass_names,
        output_bandpass_name,
    )

    # Convert via the single-convert method and add the result to the output
    output_wave = None
    output_flux = None
    for input_bandpass_name in input_bandpass_names:
        converted_spectrum = single_convert(
            spectrum=spectrum,
            output_bandpass_name=output_bandpass_name,
            pivot_wavelength=pivot_wavelengths[input_bandpass_name],
            calib_flux=calib_fluxes[input_bandpass_name],
            input_obs_to_target=input_obs_to_target,
            input_target_to_sun=input_target_to_sun,
            output_obs_to_target=output_obs_to_target,
            output_target_to_sun=output_target_to_sun,
            return_observation=False, # We need spectrum objects to combine them, not observations
        )
        # It should already be in the right units, but if that changes make sure this code
        # doesn't break
        converted_spectrum.convert('Angstrom')
        converted_spectrum.convert('flam')
        # Add the converted observation, weighted by the correlation factor
        if output_flux is None:
            # if this is the first observation, assign instead of adding
            output_wave = converted_spectrum.wave
            output_flux = converted_spectrum.flux * correlation_factors[input_bandpass_name]
        else:
            output_flux += converted_spectrum.flux * correlation_factors[input_bandpass_name]

    # Construct the output spectrum
    output_spectrum = S.ArraySpectrum(
        wave=output_wave,
        flux=output_flux,
        waveunits='Angstrom',
        fluxunits='flam',
        name=f"{spectrum.name}_weighted",
    )


    # Observe output spectrum using output filter (usually HST or Johnson)
    output_bandpass = utils.get_bandpass(output_bandpass_name)

    output_observation = S.Observation(
        output_spectrum, 
        output_bandpass,
        force='taper',
    )
    return output_observation

# Convert from multiple filters
# Distances are in the same units as the new horizons file, usually km
def plot_multi_convert(
        fig, ax,
        target_name, 
        input_filter_names, output_filter_name, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, 
        output_target_to_sun=BUIE_TARGET_TO_SUN,
        x_axis='lon', y_axis='vegamag', **kwargs):
    # get recorded observations from new horizons
    nh_observations = utils.get_observations(target_name, *input_filter_names)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_spectrum(target_name)
        # Parse out pivot wavelengths and fluxes 
        pivot_wavelengths = {
            input_filter_name: nh_observation[f'{input_filter_name}_pivot_wavelength']
            for input_filter_name in input_filter_names
        }
        fluxes = {
            input_filter_name: nh_observation[f'{input_filter_name}_calib_flux']
            for input_filter_name in input_filter_names
        }

        converted_obs = multi_convert(
            spectrum=spectrum,
            input_bandpass_names=input_filter_names,
            output_bandpass_name=output_filter_name,
            pivot_wavelengths=pivot_wavelengths,
            calib_fluxes=fluxes,
            input_obs_to_target=nh_observation['obs_to_target'],
            input_target_to_sun=nh_observation['sun_to_target'],
            output_obs_to_target=output_obs_to_target,
            output_target_to_sun=output_target_to_sun,
        )
        converted_obs.convert('Angstrom')

        plot_x.append(nh_observation[x_axis])
        plot_y.append(converted_obs.effstim(y_axis))

    # Write out as csv to file for tables
    utils.write_to_output(
        name='.'.join([target_name, 'combined', output_filter_name, y_axis]),
        data=[f"{x}:{y}" for x,y in zip([x['met'] for x in nh_observations], plot_y)],
    )
    ax.scatter(plot_x, plot_y, **kwargs)
    

def get_single_spec_factor_table(spectrum_name, input_filter_names, output_filter_names, **kwargs):
    # Using a pandas dataframe so we can use the builtin to_latex function
    df = pd.DataFrame(
        np.zeros((len(input_filter_names), len(output_filter_names))),
        index=input_filter_names,
        columns=output_filter_names,
        dtype=float,
    )

    for output_filter_name in output_filter_names:
        spectrum = utils.get_spectrum(spectrum_name)
        cfs = get_correlation_factors(
            spectrum,
            input_filter_names, # all filters except filter_name
            output_filter_name,
        )
        # fill in dataframe
        for input_filter_name, cf in cfs.items():
            df.loc[input_filter_name][output_filter_name] = cf

    # Rename the index to the pretty names
    df.rename(index=utils.LABELMAP, columns=utils.LABELMAP, inplace=True)
    # Replace 0 with None
    df[df < 1e-5] = None
    return df.T.to_latex(**kwargs)


def get_factor_table(spectrum_names, input_filter_names, output_filter_names, **kwargs):
    # Get the row labels from the filter names
    input_filter_labels = [utils.LABELMAP[x] for x in input_filter_names]

    index = pd.MultiIndex.from_product(
        [spectrum_names, input_filter_names], 
        names=['Target', 'Input Filter']
    )

    # Using a pandas dataframe so we can use the builtin to_latex function
    df = pd.DataFrame(
        np.zeros((len(spectrum_names)*len(input_filter_names), len(output_filter_names))),
        index=index,
        columns=output_filter_names,
        dtype=float,
    )
    for spectrum_name in spectrum_names:
        # Get target spectrum
        spectrum = utils.get_spectrum(spectrum_name)
        for output_filter_name in output_filter_names:
            # Get cfs
            cfs = get_correlation_factors(
                spectrum,
                input_filter_names, 
                output_filter_name
            )
            # Put cfs in dataframe
            for input_filter_name in cfs:
                df.loc[spectrum_name, input_filter_name][output_filter_name] = cfs[input_filter_name]

    # Rename the index to the pretty names
    df.rename(index=utils.LABELMAP, columns=utils.LABELMAP, inplace=True)
    # Replace 0 (ish) with None
    df[df < 1e-5] = None
    return df.T.to_latex(**kwargs)


def get_full_data_table():
    pluto_data = utils.get_observations('pluto', 'NH_RED', 'NH_BLUE')
    charon_data = utils.get_observations('charon', 'NH_RED', 'NH_BLUE')

    columns = [
        ['counts', 'NH_RED'],
        ['counts', 'NH_BLUE'],
        ['flux', 'NH_RED'],
        ['flux', 'NH_BLUE'],
        ['converted', 'HST_F555W'],
        ['converted', 'HST_F435W'],
    ]
    columns = pd.MultiIndex.from_tuples(columns)

    pluto_index = [('pluto', x['met']) for x in pluto_data]
    charon_index = [('charon', x['met']) for x in charon_data]
    index = pd.MultiIndex.from_tuples(pluto_index + charon_index)
    index.set_names(['Target', 'MET'])

    df = pd.DataFrame(index=index, columns=columns)

    data = {
        'pluto': pluto_data,
        'charon': charon_data,
    }

    for target, data in data.items():
        for obs in data:
            for input_filter in ['NH_RED', 'NH_BLUE']:
                df.loc[target, obs['met']]['counts', input_filter] = obs[f'{input_filter}_counts']
                df.loc[target, obs['met']]['flux', input_filter] \
                    = obs[f'{input_filter}_calib_flux'] #gross linewrap is gross

    # Read output.csv
    with open(utils.OUTPUT_DATA_FILE, 'r') as output_csv:
        for line in output_csv.readlines():
            label, *csv_data = line.split(',')
            target, method, output_filter, units = label.split('.')
            csv_data = [x.split(':') for x in csv_data]
            csv_data = [(int(met), float(x)) for met, x in csv_data]
            # F435W only overlaps with blue; F555W should use combined
            if output_filter == 'HST_F435W':
                desired_method = 'NH_BLUE'
            elif output_filter == 'HST_F555W':
                desired_method = 'combined'
            # write data
            if method == desired_method:
                for met, converted_value in csv_data:
                    df.loc[target, met]['converted', output_filter] = converted_value

    # Rename the index to the pretty names
    df.rename(index=utils.LABELMAP, columns=utils.LABELMAP, inplace=True)
    df.rename(inplace=True, columns={
        "counts": "Counts",
        "flux": "Flux",
        "converted": "Transformed HST Magnitude",
    })
    return df



def plot_buie(
        fig, ax,
        target, band, 
        x_field='Elon', **kwargs):

    kwargs.setdefault('label', f'Buie {utils.LABELMAP[target]} {band}')
    kwargs.setdefault('c', 'black')

    data = np.genfromtxt(f'data/buie_{target}_{band.lower()}.csv', delimiter=',', names=True)
    ax.scatter(
        data[x_field], data[band], 
        **kwargs
    )

def get_fourier(data, x_range=(0,360)):
    f = lambda x: x

    for n, an, bn in data:
        f = lambda x: f(x) + np.cos(x_range[0] + 2.0*np.pi/x_range[1] * x * n) * an
        f = lambda x: f(x) + np.sin(x_range[0] + 2.0*np.pi/x_range[1] * x * n) * bn

    return f


# Data should look like a list of (n, an, bn)
def plot_fourier(fig, ax, data, x_range=(0,360), step=1, offset=0, **kwargs):
    kwargs.setdefault('label', f'fourier series')
    kwargs.setdefault('c', 'green')

    x = np.arange(*x_range, step, dtype=float)
    y = np.zeros(x.shape)
    for n, an, bn in data:
        y += np.cos(2.0*np.pi/x_range[1] * x * n) * an
        y += np.sin(2.0*np.pi/x_range[1] * x * n) * bn

    y += offset

    ax.plot(x, y, **kwargs)
