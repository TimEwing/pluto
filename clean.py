docstr = """
Clean up .json files by removing useless keys and reordering stuff.
"""


# stdlib imports
import argparse
import json
# third party imports
import pandas as pd


def clean(input_filename, output_filename):
    with open(input_filename, 'r') as f:
        input_file = json.load(f)

    # many try-excepts here to avoid crashing on rerun

    try:
        # these don't really have anything useful in them
        del input_file['errorbar_upper']
        del input_file['errorbar_lower']
        del input_file['errorbar_none']
    except KeyError:
        pass

    # rename columns in the data
    df = pd.DataFrame.from_dict(input_file['data'])
    # get target columns
    target_cols = [x for x in df.columns if '_none' in x]
    # make a dict like "flux_none": "flux"
    column_rename_map = {col: col.replace('_none', '') for col in target_cols}
    # use the dict to rename the columns
    df.rename(columns=column_rename_map, inplace=True)

    # add sigma columns with the upper/lower error as differences
    for column in df.columns:
        if '_upper' in column:
            base_column = column.replace('_upper', '')
            new_column = f"{base_column}_sigma_plus"
            df[new_column] = df[column] - df[base_column]
        if '_lower' in column:
            base_column = column.replace('_lower', '')
            new_column = f"{base_column}_sigma_minus"
            df[new_column] = df[base_column] - df[column]

    # reorder the columns
    first_columns = [
        'met',
        'obs_to_target',
        'sun_to_target',
        'lon',
        'lat',
        'phase',
        'side',
    ]
    second_columns = sorted([x for x in df.columns if x not in first_columns])
    new_columns = first_columns + second_columns
    # do the reordering
    df = df[new_columns]
    # convert back to python dict
    input_file['data'] = df.to_dict('records')

    # Save to output file
    with open(output_filename, 'w') as f:
        json.dump(input_file, f, indent=2) # indent for readable json file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert observations from one filter to another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=docstr,
    )

    parser.add_argument(
        "input", 
        help="filename containing input data",
    )
    parser.add_argument(
        "--output", 
        help="output file name",
        default="output.json",
    )

    args = parser.parse_args()
    clean(
        args.input,
        output_filename=args.output
    )
