docstr = """

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
            # ignore 'data'
            if key == 'data':
                continue
            # Add filter name as prefix, like "NH_RED_target"
            new_key = f"{data['bandpass_name']}_{key}"
            output_data[new_key] = value

        # Convert 'data' to pandas dataframes so we can do some merging 
        data['data'] = pd.DataFrame.from_dict(data['data']).set_index('met')

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
    ## Remove columns with duplicate data
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

