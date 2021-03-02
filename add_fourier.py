docstr = """
This script adds a datapoint for each row in the input .json: the output of
a fourier fit to that longitude using the fourier series values in the 
fourier_file.

The fourier file should have 5 columns and a header row, like
    n,an,σa,bn,σb
where n is the fourier series n, an is the cosine coefficient, and bn is the 
sine coefficient. 
"""

# stdlib imports
import argparse
import json
import math
# third party imports
import numpy as np
# custom module imports
import utils


def fourier_from_json(output_filename, data_file, fourier_file):
    with open(data_file, 'r') as f:
        data = json.load(f)

    fourier = np.genfromtxt(
        fourier_file, 
        delimiter=',', 
        names=True
    )
    fourier = [[n, an, bn] for n, an, asn, bn, bsn in fourier]
    fourier_fn = get_fourier_fn(fourier)

    for datapoint in data['data']:
        fourier_vegamag = fourier_fn(datapoint['lon'])

        # add deltas for each vegamag found
        vegamag_keys = [k for k in datapoint.keys() if 'vegamag' in k]
        for key in vegamag_keys:
            datapoint[f"{key}_delta"] = fourier_vegamag - datapoint[key]

        datapoint['fourier_vegamag'] = fourier_vegamag


    with open(output_filename, 'w') as f:
        json.dump(data, f, indent=2) # indent for readable json file


# outer function which returns a fourier function
def get_fourier_fn(data, x_range=(0,360)):
    # inner function that gets returned
    def fourier_fn(x):
        y = 0
        for n, an, bn in data:
            y = y + math.cos(x_range[0] + 2.0*math.pi/x_range[1] * x * n) * an
            y = y + math.sin(x_range[0] + 2.0*math.pi/x_range[1] * x * n) * bn
        return y

    return fourier_fn


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert observations from one filter to another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "data", 
        help="observation file",
    )
    parser.add_argument(
        "fourier", 
        help="file of fourier coefficients",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default="output.json",
    )

    args = parser.parse_args()
    fourier_from_json(
        output_filename=args.output,
        data_file=args.data,
        fourier_file=args.fourier,
    )