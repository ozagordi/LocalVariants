#!/usr/bin/env python3
"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mlocalvariants` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``localvariants.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``minvar.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse
import sys

from pkg_resources import (DistributionNotFound, get_distribution)

try:
    __version__ = get_distribution('localvariants').version

except DistributionNotFound:
    # package is not installed
    pass


# if __name__ == "__main__" and __package__ is None:

# default action is 'store'
parser = argparse.ArgumentParser(description='Local structure',
                                 epilog='Input are mandatory')
parser.add_argument('-s', '--support', dest='support',
                    help='support file')
parser.add_argument('-r', '--region', dest='region', default=None,
                    help='region in chr:start-stop format')
parser.add_argument('-f', '--reference', dest='reference',
                    help='fasta file with the reference used in \
                    running shorah')

# exit so that log file is not written
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()


def main(args=None):
    """What the main does."""
    import logging
    import logging.handlers

    args = parser.parse_args()

    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'
    logging.basicConfig(filename='locstr.log', level=logging.INFO, format=log_format, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from localvariants import LStructure as LS

    sample_ls = LS.LocalStructure(support_file=args.support, ref=args.reference,
                                  region=args.region)

    sample_ls.alignedvariants(threshold=0.95)

    sample_ls.print_mutations(args.reference, out_format='csv',
                              out_file='mutations_DNA.csv')
