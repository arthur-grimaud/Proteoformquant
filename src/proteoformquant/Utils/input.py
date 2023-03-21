import argparse
from logging import warning

import os

try:
    from proteoformquant.Utils import parameters_default
except ImportError:
    from Utils import parameters_default


def doArgs(argList, name):
    """Parse argument of the comand line prompt"""
    parser = argparse.ArgumentParser(description=name)

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose debugging", default=False
    )

    parser.add_argument(
        "-i", "--indent", action="store", dest="indent_file", type=str, help="Input file name", required=True
    )

    parser.add_argument(
        "-s",
        "--spectra",
        action="store",
        dest="spectra_file",
        type=str,
        help="Input file name",
        required=True,
    )

    parser.add_argument(
        "-d",
        "--output_dir",
        action="store",
        dest="output_dir",
        type=str,
        help="Path to output directory",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--output_file",
        action="store",
        dest="output_file",
        type=str,
        help="Output file name",
        required=False,
    )

    parser.add_argument(
        "-p",
        "--param_file",
        action="store",
        dest="param_file",
        type=str,
        help="Path to json parameter file",
        required=False,
        default="params.jsonc",
    )

    parser.add_argument(
        "-cp",
        "--create_param_file",
        action="store_true",
        dest="create_param_file",
        help="Generate a json parameter file in working directory",
        required=False,
        default=False,
    )

    parser.add_argument(
        "-r",
        "--rankquant",
        action="store",
        dest="rankquant",
        type=int,
        help="Quantification using the top n ranked peptides (skips validation of psm on elution profile)",
        required=False,
        default=0,
    )

    # check id create_param_file is set to true
    if "-cp" in argList or "--create_param_file" in argList:
        print("Creating default parameter file ('params.jsonc') in working directory")
        parameters_default.generate_default_param_file()
        exit()

    args, unknownargs = parser.parse_known_args(argList)

    print_args(args, unknownargs)

    return args, unknownargs


def checkArgs(args):
    """Verifies if inputs files and output folder exists, wreate output folder if necessary"""
    if not os.path.isfile(args.indent_file):
        raise Exception("Input file: " + args.indent_file + ", doesn't exist")

    if not os.path.isfile(args.spectra_file):
        raise Exception("Input file: " + args.spectra_file + ", doesn't exist")

    if args.output_dir != "" and args.output_dir != None:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        args.output_dir = os.path.abspath(args.output_dir)

    else:
        args.output_dir = "."

    return args


def print_args(args, unknownargs):
    print("\n---Running PFQ with the following parameters:---")
    params = ""
    for arg, val in vars(args).items():
        params = params + f"\n{arg} : {val} "
    print(params, "\n")

    params = ""

    if len(unknownargs) > 0:
        dict_unknwownargs = {unknownargs[i][1:]: unknownargs[i + 1] for i in range(0, len(unknownargs), 2)}
        print("Parameters overwrite:")
        for arg, val in dict_unknwownargs.items():
            params = params + f"\n{arg} : {val} "
        print(params)
