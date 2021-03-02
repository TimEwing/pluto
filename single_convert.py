docstr = """
This script takes in observations in the form of a .json file and an output 
bandpass name and outputs a .json file with the converted observations in flam 
and in vegamag.

The input file must have the following keys:
    target: the name of the target, all lowercase (e.g. pluto, charon). This 
        must match the name of one of the targets in get_spectrum in utils.py
    pivot_wavelength: the pivot wavelength of the filter for this data
    data: a list of dicts, each with the following keys:
        flux: the observed flux, in flam
        obs_to_target: the input spacecraft range to the target
        sun_to_target: the distance from the target to the sun for the 
            input observation

The output file will be a copy of the input file, with the following fields 
added:
    output_bandpass: the name of the output filter, i.e. HST_F435W or JOHNSON_B
    converted_obs_to_target: the new spacecraft range to which all data is normalized
    sun_to_target: the new distance from the target to the sun
    data:
        converted_flux: the converted flux, in flam
        converted_vegamag: the converted flux, in the vegamag system
"""

# stdlib imports
import argparse
import json
from datetime import datetime
# custom module imports
from utils import S
import utils

def convert_json(
        input_filename, output_filename, 
        output_bandpass_name,
        output_obs_to_target, sun_to_target):
    with open(input_filename, 'r') as f:
        data = json.load(f)

    data = convert_dict(
        data, 
        output_bandpass_name, 
        output_obs_to_target, sun_to_target,
    )

    if output_filename is None:
        output_filename = get_output_filename(data)
    with open(output_filename, 'w') as f:
        json.dump(data, f, indent=2) # indent for readable json file

def convert_dict(
        data, 
        output_bandpass_name, 
        output_obs_to_target, sun_to_target):
    # get spectrum of target from previous observations
    spectrum = utils.get_spectrum(data['target'])
    spectrum.convert('Angstrom')
    spectrum.convert('flam')

    # get lab-calibrated filter bandpass
    output_bandpass = utils.get_bandpass(output_bandpass_name)
    output_bandpass.convert('Angstrom')

    # Normalize data
    for datapoint in data['data']:
        # Match the value of the spectrum at the pivot wavelength to the flux.
        # Since we only perform multiplication on the spectrum, we can reuse it
        # between loops. We need the full spectrum to calculate vegamag, so we
        # can't just leave the spectrum out of it however.
        spectrum = normalize_spectrum_to_wavelength(
            spectrum,
            data['pivot_wavelength'],
            datapoint['flux'],
        )

        # Normalize spectrum to output distance
        spectrum = normalize_spectrum_to_distance(
            spectrum,
            input_obs_to_target_distance=datapoint['obs_to_target'],
            input_sun_to_target_distance=datapoint['sun_to_target'],
            output_obs_to_target_distance=output_obs_to_target,
            output_sun_to_target_distance=sun_to_target,
        )

        # Perform synthetic observation, as if looking at the spectrum through
        # the output bandpass 
        observation = S.Observation(
            spectrum,
            output_bandpass,
            force='taper',
        )

        # Add results to the data dict
        datapoint['converted_flux'] = observation.effstim('flam')
        datapoint['converted_vegamag'] = observation.effstim('vegamag')

    data['output_bandpass'] = output_bandpass_name
    data['converted_obs_to_target'] = output_obs_to_target
    data['sun_to_target'] = sun_to_target

    # Ensure that the 'data' entry is last in the json file for readability
    data['data'] = data.pop('data')

    return data


def normalize_spectrum_to_wavelength(spectrum, 
        pivot_wavelength, observed_value):
    expected_value = spectrum.sample(pivot_wavelength)
    # Adjust the spectrum by the ratio of observed/simulated
    correction_ratio = observed_value / expected_value
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
    obs_ratio = (input_obs_to_target_distance/output_obs_to_target_distance)**2
    # Correct for different target->sun differences (lower illumination further
    # away)
    sun_ratio = (input_sun_to_target_distance/output_sun_to_target_distance)**2
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


def get_output_filename(data):
    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S.%f')
    target = data['target']
    bandpass_name = data['bandpass_name']
    return f"cnv.{target}.{bandpass_name}.{timestamp}.json"


# Parse args when run from bash
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert observations from one filter to another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "data", 
        help=".json file containing observation data",
    )
    parser.add_argument(
        "filter", 
        help="output filter name, e.g. HST_F435W or JOHNSON_B",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default=None,
    )
    parser.add_argument(
        "--d", 
        help="spacecraft range to which the data should be normalized; default \
              5.760E9 km (38.5 AU)",
        default=5.760E9,
    )
    parser.add_argument(
        "--r", 
        help="distance from the target (e.g. pluto) to the sun to which the \
              data should be normalized, default 5.909E9 km (39.5 AU)",
        default=5.909E9,
    )
    args = parser.parse_args()

    convert_json(
        input_filename=args.data,
        output_filename=args.output,
        output_bandpass_name=args.filter,
        output_obs_to_target=args.d,
        sun_to_target=args.r,
    )
