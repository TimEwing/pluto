# third party imports
import matplotlib.pyplot as plt

# custom module imports
from utils import S
import utils


def normalize_spectrum_to_bandpass(spectrum, bandpass, observed_effstim):
    # Simulate a synthetic observation of the detector should see
    syn_observation = S.Observation(spectrum, bandpass)
    # Adjust the spectrum by the ratio of observed/simulated
    correction_ratio = observed_effstim/syn_observation.effstim('photlam')
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
        output_obs_to_target_distance, output_sun_to_target_distance,):
    # Correct for different observation distances
    obs_ratio = input_obs_to_target_distance / output_obs_to_target_distance
    # Correct for different target->sun differences (lower illumination further away)
    sun_ratio = input_sun_to_target_distance / output_sun_to_target_distance
    # Correct the flux
    corrected_flux = spectrum.flux * obs_ratio * sun_ratio
    # Create a new spectrum because spectrum.flux is readonly
    corrected_spectrum = S.ArraySpectrum(
        wave=spectrum.wave,
        flux=corrected_flux,
        waveunits=spectrum.waveunits,
        fluxunits=spectrum.fluxunits,
        name=f"{spectrum.name}_corrected-distance",
    )
    return corrected_spectrum


# Simulate an observation through multiple filters
def observe_bandpasses(spectrum, bandpasses):
    # Observe the spectrum for the first bandpass
    observation = S.Observation(spectrum, bandpasses[0]) 
    # Reobserve for the rest of the bandpasses
    for bandpass in bandpasses[1:]:
        observation = S.Observation(observation, bandpass)

    return observation

# This if statement runs only if the module is called from the command line
# If it gets imported, __name__ will not be set to '__main__'
if __name__ == '__main__':
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
        spectrum = normalize_spectrum_to_bandpass(input_spectrum, input_bandpass, input_observed_effstim)

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
