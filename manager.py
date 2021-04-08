docstr = """
This is a management script for running the full sequence of scripts to convert 
an observation using New Horizons (NH) to the filterset on Hubble Space 
Telescope (HST). 

The following scripts are run in this order:
    for each errorbar (meaning 'upper'/'lower'/'none'):
        gather_data.py      # runs for each input filter
        single_convert.py   # runs for each input filter
        merge.py            # runs once per errorbar, merge by input filter
        multi_convert.py    # runs once per errorbar
    merge.py        # runs once, merge by errorbar
    add_fourier.py  # runs once
    fit_fourier.py  # runs once
"""

# stdlib imports
import argparse
from subprocess import Popen, STDOUT
from datetime import datetime
import os


def manager(args):
    final_file, intermediate_product_filenames = run_all(args)

    os.rename(final_file, args.output_filename)

    if not args.intermediate_products:
        for filename in intermediate_product_filenames:
            os.remove(filename)


def run_all(args):
    intermediate_product_filenames = []

    error_files = []
    for errorbar in ['upper', 'lower', 'none']:
        args.errorbar = errorbar
        current_files = run_gather_data(args)
        intermediate_product_filenames += current_files

        current_files = run_single_convert(current_files, args)
        intermediate_product_filenames += current_files

        current_files = run_merge_bandpass(current_files, args)
        intermediate_product_filenames += current_files

        current_files = run_multi_convert(current_files, args)
        intermediate_product_filenames += current_files
        error_files += current_files

    current_files = run_merge_errorbars(error_files, args)
    intermediate_product_filenames += current_files

    if args.fourier:
        current_files = run_add_fourier(current_files, args.fourier, args)
        intermediate_product_filenames += current_files

    if not args.nofit:
        # if fitting, we need to clean before to get the sigmas
        current_files = run_clean(current_files, args)
        intermediate_product_filenames += current_files

        current_files = run_fit_fourier(current_files, args)
        intermediate_product_filenames += current_files

    current_files = run_clean(current_files, args)
    intermediate_product_filenames += current_files
    
    # Remove final file from the intermediate product list and rename it
    final_file, = current_files
    intermediate_product_filenames.remove(final_file)
    return final_file, intermediate_product_filenames


def run_gather_data(args):
    target = args.target
    input_bandpass_names = args.input_bandpass_names
    errorbar = args.errorbar

    # We're running concurrently because I need practice and it isn't really 
    # that much harder
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


def run_merge_bandpass(filenames, args):
    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.mbp.json"
    args.JOB_ID += 1
    command = f"python merge.py {' '.join(filenames)} --key bandpass_name \
        --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return [output_filename]


def run_merge_errorbars(filenames, args):
    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.mer.json"
    args.JOB_ID += 1
    command = f"python merge.py {' '.join(filenames)} --key errorbar \
        --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return [output_filename]


def run_multi_convert(filenames, args):
    filename, = filenames

    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.mlt.json"
    args.JOB_ID += 1
    command = f"python multi_convert.py {filename} --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return [output_filename]


def run_add_fourier(filenames, fourier, args):
    filename, = filenames

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

    return [output_filename]


def run_fit_fourier(filenames, args):
    filename, = filenames

    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.fit.json"
    args.JOB_ID += 1

    n = 5 if args.target == 'pluto' else 3

    command = f"python fit_fourier.py {filename} {n} --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return [output_filename]


def run_clean(filenames, args):
    filename, = filenames

    output_filename = f"{args.JOB_STR}.{args.JOB_ID}.cln.json"
    args.JOB_ID += 1

    command = f"python clean.py {filename} --output {output_filename}"

    # start subprocess
    proc = Popen(command, shell=True, stderr=STDOUT)

    # wait for process to finish
    proc.wait()
    if args.verbose:
        print(proc.stdout.read())

    return [output_filename]


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
        "--nofit", 
        help="don't add a fourier fit to the end product",
        action='store_true',
    )
    parser.add_argument(
        "--fourier",
        help="fourier file for Buie's fourier fit",
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

    manager(args)