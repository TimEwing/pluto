docstr = """
Merge two .json files with data into one. Behavior is not defined if datasets
contain different keys

There are 3 different types of field that will be handled:

1. 'data' entry
    Datapoints will be merged by the 'met' field. Fields which match across all
    datasets will just be added like 'obs_to_target:value'; fields which do not 
    will be added like 'flux_NH_BLUE:value'.

2. 'metadata' (i.e. not the 'data' field) which is unique per dataset
    New fields will be added in the format bandpass_key, 
    e.g. 'NH_BLUE_pivot_wavelength'.

3. 'metadata' which matches for all datasets
    No change; output will be like 'target:value'
"""

# stdlib imports
import argparse
import json
# third party imports
import pandas as pd
# custom module imports
import utils

def merge(*input_filenames, output_filename):
    input_data = []
    for input_filename in input_filenames:
        with open(input_filename, 'r') as f:
            input_data.append(json.load(f))

    # clean up the data and set up the output data dict
    output_data = {}
    for data in input_data:
        # Add metadata (target, pivot wavelength, etc) into output
        for key, value in data.items():
            # ignore 'data' for now, we'll deal with in a moment
            if key == 'data':
                continue
            # We want our metadata as key:set(values) for now
            try:
                output_data[key].add(value)
            except KeyError:
                # this set hasn't been made yet, so make it
                output_data[key] = {value}

        # Convert 'data' to pandas dataframes so we can do some merging 
        data['data'] = pd.DataFrame.from_dict(data['data']).set_index('met')

    # now clean up the metadata. If the length is 1, all values matched and we
    # can just use that value. Otherwise, we need to keep all the values, which
    # we'll do by adding metadata like {pivot_wavelength_NH_RED: value}

    # find keys where the set has more than one value, i.e. the items aren't
    # unique
    unique_keys = [
        k for k,v in output_data.items() 
        if k != 'data' and len(v) == 1
    ]
    nonunique_keys = [
        k for k,v in output_data.items() 
        if k != 'data' and len(v) > 1
    ]
    for key in nonunique_keys:
        for data in input_data:
            new_key = f"{key}_{data['bandpass_name']}"
            output_data[new_key] = data[key]
        output_data.pop(key)
    # act on the rest of the keys (where the set has exactly 1 value in it)
    for key in unique_keys:
        # the comma here turns the set into just the value of the only element
        output_data[key], = output_data[key]

    ## Join DataFrames
    # copy columns from first dataframe in list for use later
    columns = input_data[0]['data'].columns
    # store first df in list in output_dataframe
    output_dataframe = input_data[0]['data']
    # join output_dataframe with each dataframe in list
    for data in input_data:
        suffix = f"_{data['bandpass_name']}"
        output_dataframe = output_dataframe.join(
            data['data'], 
            lsuffix='', 
            rsuffix=suffix,
        )

    ## Remove columns in 'data' with duplicate data
    suffixes = [f"_{d['bandpass_name']}" for d in input_data]
    for column in columns:
        all_match = all([
            output_dataframe[column].equals(output_dataframe[column+suffix])
            for suffix in suffixes
        ])
        # All columns like lon_NH_RED should match, so only keep one copy
        if all_match:
            output_dataframe[column] = output_dataframe[column+suffix]
            output_dataframe.drop(
                [column+s for s in suffixes], 
                'columns', 
                inplace=True
            )
        else:
            output_dataframe.drop(column, 'columns', inplace=True)

    # Set output data
    output_data['data'] = output_dataframe.reset_index().to_dict('records')

    # Save to output file
    with open(output_filename, 'w') as f:
        json.dump(output_data, f, indent=2) # indent for readable json file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert observations from one filter to another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "inputs", 
        nargs='+',
        help="filenames containing input data",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default="output.json",
    )

    args = parser.parse_args()
    merge(
        *args.inputs,
        output_filename=args.output,
    )

