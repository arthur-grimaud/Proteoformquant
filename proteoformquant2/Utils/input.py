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
        "-o", "--output", action="store", dest="output_dir", type=str, help="Output file name", required=False
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

    return parser.parse_args(argList)


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

    # outputBase = os.path.dirname(output_dir)
    # if outputBase != '' and not os.path.exists(outputBase):
    #     print "Output directory doesn't exist, making output dirs: %s" % (outputBase)
    #     os.makedirs(outputBase)

    return args
