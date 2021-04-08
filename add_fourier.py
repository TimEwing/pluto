docstr = """
This script adds a datapoint for each row in the input .json: the output of
a fourier fit to that longitude using the fourier series values in the 
fourier_file.

The fourier file should have 5 columns and a header row, like
    n,an,σa,bn,σb
where n is the fourier series n, an is the cosine coefficient, and bn is the 
sin coefficient. 
"""

# stdlib imports
import argparse
import json
import math
import os
# third party imports
import numpy as np
import pandas as pd
# custom module imports
import utils


def add_fourier(output_filename, data_file, fourier_file):
    with open(data_file, 'r') as f:
        data = json.load(f)

    fourier_fn = get_fourier_fn(fourier_file)

    for datapoint in data['data']:
        fourier_vegamag = fourier_fn(datapoint['lon'])

        # add deltas for each vegamag found
        vegamag_keys = [k for k in datapoint.keys() if 'vegamag' in k]
        for key in vegamag_keys:
            datapoint[f"{key}_buie_delta"] = datapoint[key] - fourier_vegamag

        datapoint['fourier_vegamag'] = fourier_vegamag

    # add some metadata
    data['fourier_file'] = os.path.basename(fourier_file)

    # Ensure that the 'data' entry is last in the json file for readability
    data['data'] = data.pop('data')

    with open(output_filename, 'w') as f:
        json.dump(data, f, indent=2) # indent for readable json file


# outer function which returns a fourier function
def get_fourier_fn(fourier_file, x_range=(0,360)):
    with open(fourier_file, 'r') as f:
        data = json.load(f)
    # load df
    df = pd.DataFrame(data['fourier_coefs'], columns=data['fourier_columns'])
    # select columns
    df = df[['n', 'an', 'bn']]
    # convert to list
    fourier_coefs = df.values.tolist()

    # inner function that gets returned
    def fourier_fn(x):
        y = np.zeros(np.shape(x))
        for n, an, bn in fourier_coefs:
            y = y + np.cos(x_range[0] + 2.0*np.pi/x_range[1] * x * n) * an
            y = y + np.sin(x_range[0] + 2.0*np.pi/x_range[1] * x * n) * bn
        y += data['fourier_offset']
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
    add_fourier(
        output_filename=args.output,
        data_file=args.data,
        fourier_file=args.fourier,
    )