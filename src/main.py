#!/usr/bin/env python
"""Deals with the input files, calls the other modules"""


def parse_com_line():
    """ Tries argparse first (python 2.7),
        if it fails uses optparse
        """
    options, args = None, None

    try:
        import argparse

        parser = argparse.ArgumentParser(description='\
                                         Writes a comprehensive report\
                                         with the local variants detected',
                                         epilog='Input are mandatory')
        parser.add_argument('-i', '--input', dest='input',
                            # default action is 'store',
                            help='input file: the extension\
                            [.sff|.fasta|.fastq|-support.fas]\
                            will indicate the procedure to follow')
        parser.add_argument('-r', '--ref', dest='ref',
                            help='reference file with a single sequence\
                            in fasta format\
                            (not necessary for file -support.fas')
        parser.add_argument('-o', '--output', dest='output',
                            help='output sff file')
        args = parser.parse_args()

    except ImportError:
        import optparse

        optparser = optparse.OptionParser()
        optparser.add_option("-c", "--corrected", type="string", default="",
                             dest="corrected", help="file with the corrected\
                             reads (output of diri_sampler)")
        optparser.add_option("-s", "--sff", type="string", default="",
                             dest="sff", help="original sff file")
        optparser.add_option("-o", "--output", help="output sff file",
                             type="string", dest="output")

        (options, args) = optparser.parse_args()

    if args:
        args = vars(args)
    else:
        args = vars(options)

    return options, args


#def msa_reads(reads, reference, out_file):
#    """Uses s2f.py from ShoRAH to build MSA of the reads"""
#
#    import s2f
#    thresh = 0.8
#    s2f.main(reference, reads, out_file, thresh, \
#             pad_insert=False, keep_files=True)

if __name__ == '__main__':

    import sys
    from Bio import SeqIO
    import LValign

    OPTIONS, ARGS = parse_com_line()

    INFILE = ARGS['input']
    ref_file = ARGS['ref']
    fasta_file = '.'.join(INFILE.split('.')[:-1]) + '.fasta'
    fastq_file = '.'.join(INFILE.split('.')[:-1]) + '.fastq'
    far_file = '.'.join(INFILE.split('.')[:-1]) + '.far'

    if INFILE.endswith('.sff'):
        count = SeqIO.convert(INFILE, 'sff-trim', fastq_file, 'fastq')
        print('Converted and trimmed %d reads from sff to fastq' % count)
        LValign.main(fastq_file, ref_file, far_file)
    elif INFILE.endswith('.fasta'):
        print('This is fasta, not converting')
        LValign.main(fasta_file, ref_file, far_file)
    elif INFILE.endswith('.fastq'):
        print('This is fastq, not converting')
        LValign.main(fastq_file, ref_file, far_file)
    elif INFILE.endswith('-support.fas'):
        print('This is a support file')
    else:
        sys.exit('This file is not supported')
