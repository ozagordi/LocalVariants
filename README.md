Local Variants
======

LStructure was designed to analyze the results from `diri_sampler`. Its input:

- the support file (ending in `reads-support.fas`);
- the reference file (fasta reference used to run any of the `shorah` tools).

It does some postprocessing to improve the reliability and readability of
variants:

- parses haplotypes with high posterior support (default: > 0.9) and minimum
  read count (default: 5 reads);
- corrects the single gaps frameshift inducing indels;
- resolves aminoacid indels to the closest possible correction;
- merges the identical haplotypes;
- eliminates haplotypes displaying SNVs removed looking at the strand bias;
- writes fasta files with frequencies, sorted;
- writes CSV files with the lists of single site mutations and their
  frequencies.

Default reference is HIV protease.
    Usage: LStructure.py -s support_file -r reference_file

    Options:
        -h, --help            show this help message and exit
        -s SUPPORT, --support=SUPPORT
                    support file
        -r REFERENCE, --reference=REFERENCE
                    fasta file with reference

### Requirements
3. Python3
1. [Biopython](http://biopython.org/)
2. `needle` from [EMBOSS](http://emboss.sourceforge.net/) software suite

### Example
    Add the `bin/` directory to your path and run the executable `localvar` as

    [user@host amplicon_1]$ localvar -s your_window-support.fas -r HIV-HXB2.fasta
    Reference is HIV-HXB2 from HIV
    Total reads: 1017.959

    [user@host amplicon_1]$ ls -1rht
    ...
    mutations_DNA.csv
    dna_seqs.fasta

The files listed above are the output of LocalVariants.
