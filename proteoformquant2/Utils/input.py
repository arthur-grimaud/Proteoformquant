import argparse
from logging import warning
import os


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
        help="Output directory path",
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
        "--param",
        action="store",
        dest="param_file",
        type=str,
        help="Json parameter file",
        required=False,
        default="params.jsonc",
    )
    args, unknownargs = parser.parse_known_args(argList)

    return args, unknownargs


def checkArgs(args):

    """Verifies if inputs files and output folder exists, wreate output folder if necessary"""
    if not os.path.isfile(args.indent_file):
        warning("Input file: " + args.indent_file + ", doesn't exist")
        return
    if not os.path.isfile(args.spectra_file):
        warning("Input file: " + args.spectra_file + ", doesn't exist")
        return

    # if args.dbse not in ("mascot", "comet"):
    #     warning("DBSE name: " + args.dbse + ", is not valid")
    #     return
    print(args)
    print("arg output: ", args.output_dir)

    if args.output_dir != "" and args.output_dir != None:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        args.output_dir = os.path.abspath(args.output_dir)

    else:
        args.output_dir = "."

    print("arg output: ", args.output_dir)

    return args
