# third party imports
import matplotlib.pyplot as plt
import numpy as np

# custom module imports
from utils import S
import utils


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
        output_obs_to_target, output_target_to_sun):
    # get lab-calibrated data of filter throughputs
    input_bandpass = utils.get_bandpass(input_filter_name)
    input_bandpass.convert('Angstrom')
    output_bandpass = utils.get_bandpass(output_filter_name)
    output_bandpass.convert('Angstrom')
    # get recorded observations from new horizons
    nh_observations = utils.get_nh_observations(target_name, input_filter_name)
    # setup plot
    plot_x = []
    plot_y = []

    for nh_observation in nh_observations:
        spectrum = utils.get_charon_spectrum()

        # Normalize spectrum to the observation and the distances
        spectrum = normalize_spectrum(
            spectrum=spectrum,
            pivot_wavelength=nh_observation['pivot_wavelength'],
            calib_flux=nh_observation['calib_flux'],
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

    plt.scatter(plot_x, plot_y)
    plt.title(f"{input_filter_name} converted to {output_filter_name}")
    plt.show()

def multi_convert(
        target_name
        input_filter_names, output_filter_name
        output_obs_to_target, output_target_to_sun):
    # get lab-calibrated data of filter throughputs
    input_bandpasses = []
    for input_filter_name in input_filter_names:
        input_bandpass = utils.get_bandpass(input_filter_name)
        input_bandpass.convert('Angstrom')
        input_bandpasses.append(input_bandpass)
    output_bandpass = utils.get_bandpass(output_filter_name)
    output_bandpass.convert('Angstrom')
    # get recorded observations from new horizons
    nh_observations = utils.get_nh_observations(target_name, input_filter_name)
    # setup plot
    plot_x = []
    plot_y = []


# Simulate an observation through multiple filters
def observe_bandpasses(spectrum, bandpasses):
    # Observe the spectrum for the first bandpass
    observation = S.Observation(spectrum, bandpasses[0]) 
    # Reobserve for the rest of the bandpasses
    for bandpass in bandpasses[1:]:
        observation = S.Observation(observation, bandpass)

    return observation

def filter_plots():
    # matplotlib setup
    fig, ax = plt.subplots()
    ax.set_prop_cycle(linestyle=['-', '--', ':', '-.'])

    bandpasses = ['NH_RED', ]
    for input_bandpass_name in bandpasses:
        input_bandpass = utils.get_bandpass(input_bandpass_name)
        input_spectrum = utils.get_charon_spectrum()
        input_observed_effstim = 0.1 # fake number in photlam, about 50% of what the synthetic obs is
        output_bandpass = utils.get_bandpass('HST_F555W')

        # Convert everything to angstrom
        input_spectrum.convert('Angstrom')
        input_bandpass.convert('Angstrom')
        output_bandpass.convert('Angstrom')

        # Normalize the spectrum to match the observation
        spectrum = normalize_spectrum_to_observation(input_spectrum, input_bandpass, input_observed_effstim)

        ## Get distances
        # Units don't matter as long as they're all the same
        # FAKE DATA
        input_obs_to_target_distance = 3 # AU
        input_sun_to_target_distance = 46 # AU
        output_obs_to_target_distance = 39.5 # AU
        output_sun_to_target_distance = 38.5 # AU

        # Normalize the spectrum with respect to the sc->pluto and pluto->earth distances
        spectrum = normalize_spectrum_to_distance(
            spectrum=spectrum,
            input_obs_to_target_distance=input_obs_to_target_distance,
            input_sun_to_target_distance=input_sun_to_target_distance,
            output_obs_to_target_distance=output_obs_to_target_distance,
            output_sun_to_target_distance=output_sun_to_target_distance,
        )

        # Simulate the spectrum through both bandpasses, as if you had held the input filter in front of
        # the output filter
        observation = observe_bandpasses(spectrum, [input_bandpass, output_bandpass])

        # TODO: Adjust for viewing geometry

        # Plot output spectra
        observation.convert('vegamag')
        ax.plot(
            observation.binwave, 
            observation.binflux, 
            label=observation.name,
            c='blue',
        )

    # Plot vega and input spectrum
    # Input spectra
    input_spectrum.convert('vegamag')
    ax.plot(input_spectrum.wave, input_spectrum.flux, label=input_spectrum.name, c='blue')
    ax.tick_params(axis='y', labelcolor='blue')
    ax.legend(loc='upper left')
    ax.invert_yaxis()
    ax.set_xlim(observation.wave.min(), observation.wave.max())
    ax.set_ylim(25, 0)

    ## Plot bandpasses on same plot with different scale
    ax_bandpass = ax.twinx()  # instantiate a second axes that shares the same x-axis
    ax_bandpass.tick_params(axis='y', labelcolor='red')
    ax_bandpass.set_prop_cycle(linestyle=['-', '--', ':', '-.'])
    for bandpass in bandpasses + ['HST_F555W']:
        bandpass = utils.get_bandpass(bandpass)
        bandpass.convert('angstrom')
        ax_bandpass.plot(bandpass.wave, bandpass.throughput, label=bandpass, color='red')
    ax_bandpass.legend(loc='upper right')

    plt.show()



# This if statement runs only if the module is called from the command line
if __name__ == '__main__':
    single_convert(
        target_name='charon', 
        input_filter_name='NH_BLUE',
        output_filter_name='HST_F435W', 
        output_obs_to_target=5.760E9, # 38.5 AU in km, as corrected in the Buie paper
        output_target_to_sun=5.909E9, # 39.5 AU in km, as corrected in the Buie paper
    )
    
