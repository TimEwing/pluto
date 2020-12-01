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


def normalize_spectrum_to_wavelength(spectrum, pivot_wavelength, observed_value):
    expected_value = spectrum.sample(pivot_wavelength)
    # Adjust the spectrum by the ratio of observed/simulated
    correction_ratio = observed_value/expected_value
    corrected_flux = spectrum.flux * correction_ratio
    # Create a new spectrum because spectrum.flux is readonly
    corrected_spectrum = S.ArraySpectrum(
        wave=spectrum.wave,
        flux=corrected_flux,
        waveunits=spectrum.waveunits,
        fluxunits=spectrum.fluxunits,
        name=f"{spectrum.name}_corrected"
    )
    return corrected_spectrum


def normalize_spectrum_to_distance(spectrum, 
        input_obs_to_target_distance, input_sun_to_target_distance,
        output_obs_to_target_distance, output_sun_to_target_distance):
    input_unit = spectrum.fluxunits
    spectrum.convert('flam') # We need to work in a non-log unit

    # Correct for different observation distances
    obs_ratio = (input_obs_to_target_distance / output_obs_to_target_distance)**2
    # Correct for different target->sun differences (lower illumination further away)
    sun_ratio = (input_sun_to_target_distance / output_sun_to_target_distance)**2
    # Correct the flux
    corrected_flux = spectrum.flux * obs_ratio * sun_ratio
    # Create a new spectrum because spectrum.flux is read-only
    corrected_spectrum = S.ArraySpectrum(
        wave=spectrum.wave,
        flux=corrected_flux,
        waveunits=spectrum.waveunits,
        fluxunits=spectrum.fluxunits,
        name=f"{spectrum.name}_corrected-distance",
    )
    corrected_spectrum.convert(input_unit)
    return corrected_spectrum


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


# Convert an observation to an output filter and distance
def single_convert(
        spectrum, output_bandpass_name,
        pivot_wavelength, calib_flux, 
        input_obs_to_target, input_target_to_sun, 
        output_obs_to_target, output_target_to_sun,
        return_observation=True):
    # get lab-calibrated data of filter throughputs
    output_bandpass = utils.get_bandpass(output_bandpass_name)
    output_bandpass.convert('Angstrom')

    spectrum.convert('Angstrom')
    spectrum.convert('flam')

    ## Normalize spectrum to observed value using builtin pysynphot renorm function
    # This sets the value of the spectrum at the pivot wavelength to calib_flux
    spectrum = normalize_spectrum_to_wavelength(
        spectrum,
        pivot_wavelength,
        calib_flux,
    )

    # Normalize spectrum to output distance
    spectrum = normalize_spectrum_to_distance(
        spectrum,
        input_obs_to_target_distance=input_obs_to_target,
        input_sun_to_target_distance=input_target_to_sun,
        output_obs_to_target_distance=output_obs_to_target,
        output_sun_to_target_distance=output_target_to_sun,
    )

    observation = S.Observation(
        spectrum,
        output_bandpass,
        force='taper',
    )

    # Multi-convert needs spectra, not observations
    if return_observation:
        return observation
    else:
        return spectrum



# Convert from one single filter
# Distances are in the same units as the new horizons file, usually km
def plot_single_convert(
        target_name, # not fully implemented yet
        input_filter_name, output_filter_name, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, 
        output_target_to_sun=BUIE_TARGET_TO_SUN):
    # get recorded observations from new horizons
    nh_observations = utils.get_observations(target_name, input_filter_name)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_spectrum(target)

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


        plot_x.append(nh_observation['obs_to_target'])
        plot_y.append(converted_obs.effstim('vegamag'))

    plt.scatter(plot_x, plot_y, label=f"{input_filter_name}", marker='x', 
        c=utils.COLORMAP[input_filter_name])


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
        target_name, # not fully implemented yet
        input_filter_names, output_filter_name, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, 
        output_target_to_sun=BUIE_TARGET_TO_SUN):
    # get recorded observations from new horizons
    nh_observations = utils.get_observations(target_name, *input_filter_names)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_spectrum(target)
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

        plot_x.append(nh_observation['obs_to_target'])
        plot_y.append(converted_obs.effstim('vegamag'))

    plt.scatter(plot_x, plot_y, label=f"{', '.join(input_filter_names)}", marker='x', 
        c='purple')
    

def get_factor_table(spectrum_names, input_filter_names, output_filter_names, **kwargs):
    # Get the row labels from the filter names
    input_filter_labels = [utils.LABELMAP[x] for x in input_filter_names]

    index = pd.MultiIndex.from_product(
        [spectrum_names, input_filter_names], 
        names=['Target', 'Input Filter'])

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
            # If the cf is 0, replace it with ''
            cfs = {k:(v if v > 1e-5 else None) for k,v in cfs.items()}
            # Put cfs in dataframe
            for input_filter_name in cfs:
                df.loc[spectrum_name, input_filter_name][output_filter_name] = cfs[input_filter_name]

    # Rename the index to the pretty names
    df.rename(index=utils.LABELMAP, columns=utils.LABELMAP, inplace=True)
    return df.T.to_latex(**kwargs)


def plot_buie_pluto_V(x_field='Elon', **kwargs):
    data = np.genfromtxt('data/buie_pluto_v.csv', delimiter=',', names=True)
    plt.scatter(data[x_field], data['V'], label='Buie Pluto V', color='Black')


# This if statement runs only if the module is called from the command line
if __name__ == '__main__':
    # filters = ['NH_BLUE', 'NH_RED']
    # output = 'HST_F555W'
    # target = 'hd'
    # for f in filters:
    #     single_convert(
    #         target_name=target, 
    #         input_filter_name=f,
    #         output_filter_name=output, 
    #         output_obs_to_target=1, # 38.5 AU in km, as corrected in the Buie paper
    #         output_target_to_sun=1, # 39.5 AU in km, as corrected in the Buie paper
    #     )
    # multi_convert(
    #     target_name=target, 
    #     input_filter_names=filters,
    #     output_filter_name=output, 
    #     output_obs_to_target=1, # 38.5 AU in km, as corrected in the Buie paper
    #     output_target_to_sun=1, # 39.5 AU in km, as corrected in the Buie paper
    # )


    filters = ['NH_BLUE', 'NH_RED']
    output = 'HST_F555W'
    target = 'hd'
    for f in filters:
        plot_single_convert(
            target_name=target, 
            input_filter_name=f,
            output_filter_name=output, 
            output_obs_to_target=BUIE_OBS_TO_TARGET, # 38.5 AU in km, as corrected in the Buie paper
            output_target_to_sun=BUIE_TARGET_TO_SUN, # 39.5 AU in km, as corrected in the Buie paper
        )
    plot_multi_convert(
        target_name=target, 
        input_filter_names=filters,
        output_filter_name=output, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, # 38.5 AU in km, as corrected in the Buie paper
        output_target_to_sun=BUIE_TARGET_TO_SUN, # 39.5 AU in km, as corrected in the Buie paper
    )

    # plot_buie_pluto_V()

    plt.legend()
    plt.title(f"{output} using {target} spectrum")
    plt.show()
    
