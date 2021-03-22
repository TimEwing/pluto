docstr = """
This is a management script for running the full sequence of scripts to convert 
an observation using New Horizons (NH) to the filterset on Hubble Space 
Telescope (HST). 

The following scripts are run in this order:

    gather_data.py      # runs for each input filter
    single_convert.py   # runs for each input filter
    merge.py            # runs once
    multi_convert.py    # runs once
    add_fourier.py      # runs once
"""

# stdlib imports
import argparse
from subprocess import Popen, STDOUT
from datetime import datetime
import os


def main(args):
    intermediate_product_filenames = []

    gather_files = run_gather_data(args)
    intermediate_product_filenames += gather_files

    single_files = run_single_convert(gather_files, args)
    intermediate_product_filenames += single_files

    merge_file = run_merge(single_files, args)
    intermediate_product_filenames.append(merge_file)

    multi_file = run_multi_convert(merge_file, args)

    if args.fourier:
        # if running fourier, multi_file is no longer the final file
        intermediate_product_filenames.append(multi_file)
        fourier_file = run_add_fourier(multi_file, args.fourier, args)
        # rename final file to output_filename
        os.rename(fourier_file, args.output_filename)
    else:
        # rename final file to output_filename
        os.rename(multi_file, args.output_filename)

    if not args.intermediate_products:
        for filename in intermediate_product_filenames:
            os.remove(filename)


def run_gather_data(args):
    target = args.target
    input_bandpass_names = args.input_bandpass_names
    errorbar = args.errorbar

    # We're running gather_data concurrently because I need practice and it 
    # isn't really that much harder
    commands = []
    output_filenames = []
    for input_bandpass_name in input_bandpass_names:
        output_filename = f"{args.JOB_STR}.{args.JOB_ID}.raw.json"
        args.JOB_ID += 1
        commands.append(f"python gather_data.py {target} {input_bandpass_name}\
            --output {output_filename} --errorbar {errorbar}")
        output_filenames.append(output_filename)

    # start subprocesses
    procs = [
        Popen(cmd, shell=True, stderr=STDOUT)
        for cmd in commands
    ]
    # wait for processes to finish
    for proc in procs:
        proc.wait()
        if args.verbose:
            print(proc.stdout.read())

    return output_filenames


def run_single_convert(filenames, args):
    target = args.target
    ouput_bandpass = args.output_bandpass_name

    # We're running concurrently because I need practice and it isn't really 
    # that much harder
    commands = []
    output_filenames = []
    for filename in filenames:
        output_filename = f"{args.JOB_STR}.{args.JOB_ID}.sng.json"
        args.JOB_ID += 1
        commands.append(f"python single_convert.py {filename} {ouput_bandpass} \
            --output {output_filename}")
        output_filenames.append(output_filename)

    # start subprocesses
    procs = [
        Popen(cmd, shell=True, stderr=STDOUT)
        for cmd in commands
    ]
    # wait for processes to finish
    for proc in procs:
        proc.wait()
        if args.verbose:
            print(proc.stdout.read())

    return output_filenames


def run_merge(filenames, args):
    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.mrg.json"
    args.JOB_ID += 1
    command = f"python merge.py {' '.join(filenames)} --output \
        {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return output_filename


def run_multi_convert(filename, args):
    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.mlt.json"
    args.JOB_ID += 1
    command = f"python multi_convert.py {filename} --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return output_filename


def run_add_fourier(filename, fourier, args):
    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.for.json"
    args.JOB_ID += 1
    command = f"python add_fourier.py {filename} {fourier} --output \
        {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return output_filename


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
        "inputs", 
        nargs='+',
        help="NH_RED, NH_BLUE, etc",
    )
    parser.add_argument(
        "output", 
        help="HST_F435W, JOHNSON_B, etc",
    )
    parser.add_argument(
        "--errorbar", 
        help="output file name",
        choices=["upper", "lower", "none"],
        default="none",
    )
    parser.add_argument(
        "--fourier",
        help="optional fourier file",
    )
    parser.add_argument(
        "--o", 
        help="output file name",
        default="output.json",
    )
    parser.add_argument(
        "-i", 
        help="include intermediate products",
        action='store_true',
    )
    parser.add_argument(
        "-v", 
        help="verbose",
        action='store_true',
    )
    args = parser.parse_args()

    args.target = args.target
    args.input_bandpass_names = args.inputs
    args.output_bandpass_name = args.output
    args.output_filename = args.o
    args.intermediate_products = args.i
    args.verbose = args.v

    timestamp = datetime.now().strftime('%y%m%d%H%M%S%f')
    args.JOB_STR = f"{timestamp}_{args.target}"
    args.JOB_ID = 0

    main(args)