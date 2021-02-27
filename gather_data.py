docstr = """
This script collects new horizons data from various sources (e.g. 
nh_obs.idlsave) and outputs a .json file with the useful bits.
"""

# stdlib imports
import argparse
import json
# custom module imports
import utils

def gather_data(output_filename, target, bandpass_name):
    output_dict = {
        'target': target,
        'bandpass_name': bandpass_name,
        'pivot_wavelength': utils.NH_PIVOT_WAVELENGTH[bandpass_name],
        'data': []
    }

    observations = utils.get_observations(target, bandpass_name)

    for observation in observations:
        # strip "NH_RED" from key names
        observation = {k.strip(bandpass_name): v for k,v in observation.items()}
        # remove pivot wavelength from data (included by utils for old versions 
        # of this code)
        observation.pop('pivot_wavelength')
        # rename calib_flux to just flux
        observation['flux'] = observation.pop('calib_flux')
        
        output_dict['data'].append(observation)

    with open(output_filename, 'w') as f:
        json.dump(output_dict, f, indent=2) # indent for readable json file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert observations from one filter to another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "target", 
        help="pluto, charon, or hd",
    )
    parser.add_argument(
        "filter", 
        help="NH_RED, NH_BLUE, etc",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default="output.json",
    )

    args = parser.parse_args()
    gather_data(
        output_filename=args.output,
        target=args.target,
        bandpass_name=args.filter,
    )

