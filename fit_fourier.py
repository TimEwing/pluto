docstr = """
This script fits data to a fourier curve using least-squares. 

Given a .json file with a "data" key, the fourier curve is fit with x as the 
datapoint's "lon" key and y as the key arg.

If 
"""

import json
import argparse
import functools

import pandas as pd
import numpy as np
from scipy import optimize


def fit_fourier(input_filename, key, n, output_filename, error=None):
    # Load input data
    with open(input_filename, 'r') as f:
        data = json.load(f)
    df = pd.DataFrame.from_dict(data['data'])
    df.sort_values('lon', inplace=True)

    p0 = [0] * (n*2) # Create a 2n long list of [0, 0, ... 0]
    x = df.lon
    y = df[key]
    yerr = (df[f"{key}_sigma_plus"], df[f"{key}_sigma_minus"])

    # pick error function
    if error == 'average':
        error_f = average_error_f
    elif error == 'select':
        error_f = select_error_f
    else:
        error_f = base_error_f

    ## perform least squares fit
    fit = optimize.leastsq(error_f, p0, args=(x, y, yerr))
    # transform coefs to look like a list of [(n, an, bn), ...]
    coefs = zip(fit[0][::2], fit[0][1::2])
    coefs = [[n, an, bn] for n, (an, bn) in enumerate(coefs)]

    # Add fourier stuff to the output
    data['fourier_yaxis'] = key
    data['fourier_offset'] = 0
    data['fourier_columns'] = ['n', 'an', 'bn']
    data['fourier_coefs'] = coefs

    # Ensure that the 'data' entry is last in the json file for readability
    data['data'] = data.pop('data')

    with open(output_filename, 'w') as f:
        json.dump(data, f, indent=2) # indent for readable json file


## Error functions for the least squares fit
def base_error_f(p, x, y, yerr):
    # ignore the y errors
    return y - fourier(p, x)
    
def average_error_f(p, x, y, yerr):
    # return the difference of x and y divided by the sum of the upper and lower
    # y errors
    yerr_upper, yerr_lower = yerr
    yerr_avg = (yerr_upper+yerr_lower) / 2
    return (y - fourier(p, x)) / yerr_avg

def select_error_f(p, x, y, yerr):
    # return the difference of x and y divided by the relavent error (i.e. if 
    # y > x, divide by yerr_upper)
    yerr_upper, yerr_lower = yerr
    # use np.select to decide which error use
    condition = y > x
    condition_inv = np.invert(condition)
    selected_err = np.select(
        [condition, condition_inv],
        [yerr_upper, yerr_lower]
    )
    return (y - fourier(p, x)) / selected_err


def fourier(p, x):
    # This function just generates a fourier curve 

    # p[0], p[2]... are a_ns
    # p[1], p[3]... are b_ns
    coefs = zip(p[::2], p[1::2]) 
    y = 0 
    for n, (an, bn) in enumerate(coefs): 
        y = y + np.cos(2.0*np.pi/359 * x * n) * an 
        y = y + np.sin(2.0*np.pi/359 * x * n) * bn 
    return y 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Given a .json with data, calculate a phase curve fit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "input", 
        help="filename containing input data",
    )
    parser.add_argument(
        "n", 
        help="number of fourier terms",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--key", 
        default="multi_vegamag",
        help="name of the field to fit data to",
    )
    parser.add_argument(
        "--error",
        choices=["average", "select", "none"],
        default='average',
        help="name of the field to fit data to",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default="output.json",
    )

    args = parser.parse_args()
    fit_fourier(
        input_filename=args.input,
        key=args.key,
        n=args.n,
        error=args.error,
        output_filename=args.output,
    )