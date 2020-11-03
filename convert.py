# third party imports
import matplotlib.pyplot as plt
import numpy as np

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


# Normalize a spectrum to match an observation and a distance
def normalize_spectrum(
        spectrum, 
        pivot_wavelength, calib_flux, 
        input_obs_to_target, input_target_to_sun, 
        output_obs_to_target, output_target_to_sun):

    spectrum.convert('Angstrom')
    spectrum.convert('flam')

    # Normalize spectrum to observed value using builtin pysynphot renorm function
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
    return spectrum


# Convert from one single filter
# Distances are in the same units as the new horizons file, usually km
def single_convert(
        target_name, # not fully implemented yet
        input_filter_name, output_filter_name, 
        output_obs_to_target=BUIE_OBS_TO_TARGET, 
        output_target_to_sun=BUIE_TARGET_TO_SUN):
    # get lab-calibrated data of filter throughputs
    input_bandpass = utils.get_bandpass(input_filter_name)
    input_bandpass.convert('Angstrom')
    output_bandpass = utils.get_bandpass(output_filter_name)
    output_bandpass.convert('Angstrom')
    # get recorded observations from new horizons
    nh_observations = utils.get_observations(target_name, input_filter_name)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_spectrum(target)

        # Normalize spectrum to the observation and the distances
        spectrum = normalize_spectrum(
            spectrum=spectrum,
            pivot_wavelength=nh_observation[f'{input_filter_name}_pivot_wavelength'],
            calib_flux=nh_observation[f'{input_filter_name}_calib_flux'],
            input_obs_to_target=nh_observation['obs_to_target'],
            input_target_to_sun=nh_observation['sun_to_target'],
            output_obs_to_target=output_obs_to_target,
            output_target_to_sun=output_target_to_sun,
        )
        
        # Observe normalized spectrum using output filter (usually HST or Johnson)
        syn_observation = S.Observation(
            spectrum, 
            output_bandpass,
            force='taper',
        )

        plot_x.append(nh_observation['obs_to_target'])
        plot_y.append(syn_observation.effstim('vegamag'))

    plt.scatter(plot_x, plot_y, label=f"{input_filter_name}", marker='x', c=utils.COLORMAP[input_filter_name])

def multi_convert_plot(
        target_name,
        input_bandpasses, output_bandpass,
        ):
    pass
    

def multi_convert(
        target_name,
        input_filter_names, output_filter_name,
        output_obs_to_target=BUIE_OBS_TO_TARGET,
        output_target_to_sun=BUIE_TARGET_TO_SUN):
    # get lab-calibrated data of filter throughputs and nh observations
    input_bandpasses = {}
    for input_bandpass_name in input_filter_names:
        input_bandpass = utils.get_bandpass(input_bandpass_name)
        input_bandpass.convert('Angstrom')
        input_bandpasses[input_bandpass_name] = input_bandpass
    # get output filter throughputs
    output_bandpass = utils.get_bandpass(output_filter_name)
    output_bandpass.convert('Angstrom')

    # find values for weighted average of input observations
    input_bandpass_fluxes = {}
    for input_bandpass in input_bandpasses:
        spectrum = utils.get_spectrum(target)
        # first, observe the input bandpass
        observation = S.Observation(spectrum, input_bandpasses[input_bandpass], force='taper')
        observation = S.Observation(observation, output_bandpass, force='taper')
        input_bandpass_fluxes[input_bandpass] = observation.effstim('flam')
    # normalize and find weight
    total_flux = sum([flux for flux in input_bandpass_fluxes.values()])
    input_bandpass_weights = {}
    for input_bandpass, flux in input_bandpass_fluxes.items():
        input_bandpass_weights[input_bandpass] = flux / total_flux

    # setup plot
    plot_x = []
    plot_y = []

    # get new horizons observations
    sc_observations = utils.get_observations(target_name, *input_filter_names)
    # convert from input filters to output filter
    for sc_observation in sc_observations:
        spectra = {}
        # We need to reconstruct the target spectrum often because we often modify it
        spectrum = utils.get_spectrum(target)
        spectrum.convert('Angstrom')
        summed_wave = spectrum.wave # Just using the target spectrum to get the length right
        summed_flux = spectrum.flux * 0 # array of 0s with the same length as the target spectrum
        for input_bandpass_name, input_bandpass in input_bandpasses.items(): 
            spectrum = utils.get_spectrum(target)
            spectrum.convert('Angstrom')
            spectrum.convert('flam')

            # Normalize spectrum to the observation and the distances
            spectrum = normalize_spectrum(
                spectrum=spectrum,
                pivot_wavelength=sc_observation[f'{input_bandpass_name}_pivot_wavelength'],
                calib_flux=sc_observation[f'{input_bandpass_name}_calib_flux'],
                input_obs_to_target=sc_observation['obs_to_target'],
                input_target_to_sun=sc_observation['sun_to_target'],
                output_obs_to_target=output_obs_to_target,
                output_target_to_sun=output_target_to_sun,
            )

            # Add the spectrum * the weight to the summed spectra
            summed_flux += spectrum.flux * input_bandpass_weights[input_bandpass_name]

        spectrum = utils.get_spectrum(target)
        summed_spectra = S.ArraySpectrum(
            wave=summed_wave,
            flux=summed_flux,
            waveunits='Angstrom',
            fluxunits='flam',
            name=f"{spectrum.name}_weighted",
        )
        # Observe normalized spectrum using output filter (usually HST or Johnson)
        syn_observation = S.Observation(
            summed_spectra, 
            output_bandpass,
            force='taper',
        )

        plot_x.append(sc_observation['obs_to_target'])
        plot_y.append(syn_observation.effstim('vegamag'))
    print(input_bandpass_weights)
    plt.scatter(plot_x, plot_y, label=f"{input_filter_names}", alpha=0.5, c='green')


# def filter_plots():
#     # matplotlib setup
#     fig, ax = plt.subplots()
#     ax.set_prop_cycle(linestyle=['-', '--', ':', '-.'])

#     bandpasses = ['NH_RED', 'NH_BLUE']
#     nh_observation = utils.get_observations('charon', *bandpasses)[0]
#     for input_bandpass_name in bandpasses:
#         input_bandpass = utils.get_bandpass(input_bandpass_name)
#         input_spectrum = utils.get_spectrum(target)
#         input_observed_effstim = 0.1 # fake number in photlam, about 50% of what the synthetic obs is
#         output_bandpass = utils.get_bandpass('HST_F555W')

#         # Convert everything to angstrom
#         input_spectrum.convert('Angstrom')
#         input_bandpass.convert('Angstrom')
#         output_bandpass.convert('Angstrom')

#         output_obs_to_target=5.760E9, # 38.5 AU in km, as corrected in the Buie paper
#         output_target_to_sun=5.909E9, # 39.5 AU in km, as corrected in the Buie paper
#         input_spectrum = normalize_spectrum(
#             spectrum=input_spectrum,
#             pivot_wavelength=nh_observation[f'{input_bandpass_name}_pivot_wavelength'],
#             calib_flux=nh_observation[f'{input_bandpass_name}_calib_flux'],
#             input_obs_to_target=nh_observation['obs_to_target'],
#             input_target_to_sun=nh_observation['sun_to_target'],
#             output_obs_to_target=output_obs_to_target,
#             output_target_to_sun=output_target_to_sun,
#         )


#         # Simulate the spectrum through both bandpasses, as if you had held the input filter in front of
#         # the output filter
#         observation = S.Observation(input_spectrum, input_bandpass)
#         observation = S.Observation(observation, output_bandpass)

#         # TODO: Adjust for viewing geometry

#         # Plot output spectra
#         observation.convert('vegamag')
#         ax.plot(
#             observation.binwave, 
#             observation.binflux, 
#             label=observation.name,
#             c='blue',
#         )

#     # Plot vega and input spectrum
#     # Input spectra
#     input_spectrum.convert('vegamag')
#     ax.plot(input_spectrum.wave, input_spectrum.flux, label=input_spectrum.name, c='blue')
#     ax.tick_params(axis='y', labelcolor='blue')
#     ax.legend(loc='upper left')
#     ax.invert_yaxis()
#     ax.set_xlim(observation.wave.min(), observation.wave.max())
#     ax.set_ylim(25, 0)

#     ## Plot bandpasses on same plot with different scale
#     ax_bandpass = ax.twinx()  # instantiate a second axes that shares the same x-axis
#     ax_bandpass.tick_params(axis='y', labelcolor='red')
#     ax_bandpass.set_prop_cycle(linestyle=['-', '--', ':', '-.'])
#     for bandpass in bandpasses + ['HST_F555W']:
#         bandpass = utils.get_bandpass(bandpass)
#         bandpass.convert('angstrom')
#         ax_bandpass.plot(bandpass.wave, bandpass.throughput, label=bandpass, color='red')
#     ax_bandpass.legend(loc='upper right')

#     plt.show()

def plot_buie_pluto_V(x_field='Elon'):
    data = np.genfromtxt('data/buie_pluto_v.csv', delimiter=',', names=True)
    plt.scatter(data[x_field], data['V'], label='Buie Pluto V', color='Black')


# This if statement runs only if the module is called from the command line
if __name__ == '__main__':
    filters = ['NH_BLUE', 'NH_RED']
    output = 'HST_F555W'
    target = 'hd'
    for f in filters:
        single_convert(
            target_name=target, 
            input_filter_name=f,
            output_filter_name=output, 
            output_obs_to_target=1, # 38.5 AU in km, as corrected in the Buie paper
            output_target_to_sun=1, # 39.5 AU in km, as corrected in the Buie paper
        )
    multi_convert(
        target_name=target, 
        input_filter_names=filters,
        output_filter_name=output, 
        output_obs_to_target=1, # 38.5 AU in km, as corrected in the Buie paper
        output_target_to_sun=1, # 39.5 AU in km, as corrected in the Buie paper
    )

    # plot_buie_pluto_V()

    plt.legend()
    plt.title(f"{output} using {target} spectrum")
    plt.show()
    
