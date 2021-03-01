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
    output_filter: the name of the output filter, i.e. HST_F435W or JOHNSON_B
    norm_obs_to_target: the new spacecraft range to which all data is normalized
    sun_to_target: the new distance from the target to the sun
    data:
        norm_flux: the converted flux, in flam
        norm_vegamag: the converted flux, in the vegamag system
"""

# stdlib imports
import argparse
import json
import math
# custom module imports
from utils import S
import utils

def multi_covert(
        input_filename,
        output_filename):
    with open(input_filename, 'r') as f:
        data = json.load(f)

    # Pull input bandpass names out
    input_bandpass_names = [
        data[x] for x in data.keys() 
        if 'bandpass_name' in x
    ]
    # Get output bandpass name
    output_bandpass_name = data['output_bandpass']
    # Get spectrum
    spectrum = utils.get_spectrum(data['target'])
    spectrum.convert('Angstrom')
    spectrum.convert('flam')

    correlation_factors = get_correlation_factors(
        spectrum,
        input_bandpass_names,
        output_bandpass_name,
    )

    # Apply weighted average by correlation factors
    for datapoint in data['data']:
        # flux first
        multi_flux = 0
        for bandpass in input_bandpass_names:
            flux_key = f"converted_flux_{bandpass}"
            multi_flux += datapoint[flux_key] * correlation_factors[bandpass]
        datapoint['multi_flux'] = multi_flux
        # Now vegamag. To scale a vegamag by a constant K we actually need to 
        # add -2.5 log10(K)
        scale_bandpass = input_bandpass_names[0]
        scale_flux_key = f"converted_flux_{scale_bandpass}"
        scale_vegamag_key = f"converted_vegamag_{scale_bandpass}"

        scale_factor = multi_flux / datapoint[scale_flux_key]
        multi_vegamag = datapoint[scale_vegamag_key] \
            - 2.5 * math.log10(scale_factor)

        datapoint['multi_vegamag'] = multi_vegamag


    # Save data to output
    with open(output_filename, 'w') as f:
        json.dump(data, f, indent=2) # indent for readable json file


# Calculate the weighted average factors for each observation
def get_correlation_factors(
        spectrum, 
        input_bandpass_names, 
        output_bandpass_name):

    output_bandpass = utils.get_bandpass(output_bandpass_name)

    # Get the flux from observing as if you held one filter in front of the 
    # other
    correlation_factors = {}
    for input_bandpass_name in input_bandpass_names:
        input_bandpass = utils.get_bandpass(input_bandpass_name)
        # Perform synthetic observations
        observation = S.Observation(spectrum, input_bandpass, force='taper')
        observation = S.Observation(observation, output_bandpass, force='taper')

        # Can't be vegamag since magnitudes don't add linearly; needs to be 
        # flam
        try:
            correlation_factors[input_bandpass_name] = \
                observation.effstim('flam')
        except ValueError:
            # 0 effstim
            correlation_factors[input_bandpass_name] = 0

    # Normalize factors
    factor_sum = sum(correlation_factors.values())
    correlation_factors = {
        k: v/factor_sum 
        for k, v in correlation_factors.items()
    }

    return correlation_factors


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
        "--output", 
        help="output file name",
        default="output.json",
    )
    args = parser.parse_args()

    multi_covert(
        input_filename=args.data,
        output_filename=args.output,
    )
