# third party imports
import matplotlib.pyplot as plt
import numpy as np

# custom module imports
from utils import S
import utils


def normalize_spectrum_to_observation(spectrum, bandpass, observed, observed_units='counts'):
    # Simulate a synthetic observation of the detector should see
    syn_observation = S.Observation(spectrum, bandpass)
    # Adjust the spectrum by the ratio of observed/simulated
    correction_ratio = observed/syn_observation.effstim(observed_units)
    corrected_flux = spectrum.flux * correction_ratio
    # Create a new spectrum because spectrum.flux is readonly
    corrected_spectrum = S.ArraySpectrum(
        wave=spectrum.wave,
        flux=corrected_flux,
        waveunits=spectrum.waveunits,
        fluxunits=spectrum.fluxunits,
        name=f"{spectrum.name}_corrected-{bandpass.name}"
    )
    return corrected_spectrum


def normalize_spectrum_to_distance(spectrum, 
        input_obs_to_target_distance, input_sun_to_target_distance,
        output_obs_to_target_distance, output_sun_to_target_distance):
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
    corrected_spectrum.convert('vegamag')
    print(obs_ratio, corrected_spectrum.sample(5000))
    return corrected_spectrum


# Simulate an observation through multiple filters
def observe_bandpasses(spectrum, bandpasses):
    # Observe the spectrum for the first bandpass
    observation = S.Observation(spectrum, bandpasses[0]) 
    # Reobserve for the rest of the bandpasses
    for bandpass in bandpasses[1:]:
        observation = S.Observation(observation, bandpass)

    return observation


# Convert from one single filter
# Distances are in the same units as the new horizons file, usually km
def single_convert(
        target_name, # not fully implemented yet
        input_filter_name, output_filter_name, 
        output_obs_to_target, output_target_to_sun):
    # Calibrated data of filter throughputs
    input_bandpass = utils.get_bandpass(input_filter_name)
    input_bandpass.convert('Angstrom')
    output_bandpass = utils.get_bandpass(output_filter_name)
    output_bandpass.convert('Angstrom')
    # Recorded observations from spacecraft
    nh_observations = utils.get_nh_observations(target_name, 'charon_red_counts')

    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_charon_spectrum()
        spectrum.convert('Angstrom')
        spectrum.convert('flam')
        # plt.plot(spectrum.wave, spectrum.flux)
        # plt.plot(output_bandpass.wave, output_bandpass.throughput)
        # plt.show()
        # break

        spectrum = normalize_spectrum_to_observation(
            spectrum,
            input_bandpass, 
            nh_observation['charon_red_counts'], 
            observed_units='counts',
        )
        spectrum = normalize_spectrum_to_distance(
            spectrum,
            input_obs_to_target_distance=nh_observation['obs_to_target'],
            input_sun_to_target_distance=nh_observation['sun_to_target'],
            output_obs_to_target_distance=output_obs_to_target,
            output_sun_to_target_distance=output_target_to_sun,
        )
        print(S.Observation(S.Vega, input_bandpass).effstim('vegamag'), nh_observation['charon_red_counts'])
        syn_observation = S.Observation(
            spectrum, 
            output_bandpass,
        )
        break

        # plot_x.append(nh_observation['charon_red_counts'])
        # plot_y.append(syn_observation.effstim('vegamag'))
    # plt.scatter(plot_x, plot_y)
    # plt.show()


# This if statement runs only if the module is called from the command line
if __name__ == '__main__':
    single_convert('charon', 'NH_RED', 'HST_F555W', 4.413E9, 4.563E9)
    # # matplotlib setup
    # fig, ax = plt.subplots()
    # ax.set_prop_cycle(linestyle=['-', '--', ':', '-.'])

    # bandpasses = ['NH_RED', ]
    # for input_bandpass_name in bandpasses:
    #     input_bandpass = utils.get_bandpass(input_bandpass_name)
    #     input_spectrum = utils.get_charon_spectrum()
    #     input_observed_effstim = 0.1 # fake number in photlam, about 50% of what the synthetic obs is
    #     output_bandpass = utils.get_bandpass('HST_F555W')

    #     # Convert everything to angstrom
    #     input_spectrum.convert('Angstrom')
    #     input_bandpass.convert('Angstrom')
    #     output_bandpass.convert('Angstrom')

    #     # Normalize the spectrum to match the observation
    #     spectrum = normalize_spectrum_to_observation(input_spectrum, input_bandpass, input_observed_effstim)

    #     ## Get distances
    #     # Units don't matter as long as they're all the same
    #     # FAKE DATA
    #     input_obs_to_target_distance = 3 # AU
    #     input_sun_to_target_distance = 46 # AU
    #     output_obs_to_target_distance = 39.5 # AU
    #     output_sun_to_target_distance = 38.5 # AU

    #     # Normalize the spectrum with respect to the sc->pluto and pluto->earth distances
    #     spectrum = normalize_spectrum_to_distance(
    #         spectrum=spectrum,
    #         input_obs_to_target_distance=input_obs_to_target_distance,
    #         input_sun_to_target_distance=input_sun_to_target_distance,
    #         output_obs_to_target_distance=output_obs_to_target_distance,
    #         output_sun_to_target_distance=output_sun_to_target_distance,
    #     )

    #     # Simulate the spectrum through both bandpasses, as if you had held the input filter in front of
    #     # the output filter
    #     observation = observe_bandpasses(spectrum, [input_bandpass, output_bandpass])

    #     # TODO: Adjust for viewing geometry

    #     # Plot output spectra
    #     observation.convert('vegamag')
    #     ax.plot(
    #         observation.binwave, 
    #         observation.binflux, 
    #         label=observation.name,
    #         c='blue',
    #     )

    # # Plot vega and input spectrum
    # # Input spectra
    # input_spectrum.convert('vegamag')
    # ax.plot(input_spectrum.wave, input_spectrum.flux, label=input_spectrum.name, c='blue')
    # ax.tick_params(axis='y', labelcolor='blue')
    # ax.legend(loc='upper left')
    # ax.invert_yaxis()
    # ax.set_xlim(observation.wave.min(), observation.wave.max())
    # ax.set_ylim(25, 0)

    # ## Plot bandpasses on same plot with different scale
    # ax_bandpass = ax.twinx()  # instantiate a second axes that shares the same x-axis
    # ax_bandpass.tick_params(axis='y', labelcolor='red')
    # ax_bandpass.set_prop_cycle(linestyle=['-', '--', ':', '-.'])
    # for bandpass in bandpasses + ['HST_F555W']:
    #     bandpass = utils.get_bandpass(bandpass)
    #     bandpass.convert('angstrom')
    #     ax_bandpass.plot(bandpass.wave, bandpass.throughput, label=bandpass, color='red')
    # ax_bandpass.legend(loc='upper right')

    # plt.show()
