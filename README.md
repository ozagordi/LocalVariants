Local Variants
======

LStructure was designed to analyze the results from `diri_sampler`. It takes
as input the support file and a reference file, then it does some
postprocessing to improve the reliability and readability of variants:

- parses haplotypes with high posterior support (default: > 0.9) and minimum
read count (default: 5 reads);
- corrects the single gaps frameshift inducing indels;
- resolves aminoacid indels to the closest possible correction;
- merges the identical haplotypes;
- writes a fasta file with frequencies, sorted;
- writes CSV files with the lists of single site mutations and their
frequencies both for DNA (now superseded by SNV support in ShoRAH 0.6) and
for aminoacids.

Default reference is HIV protease.
    Usage: LStructure.py -s support_file [options]

    Options:
        -h, --help            show this help message and exit
        -s SUPPORT, --support=SUPPORT
                    support file
        -g GENE, --gene=GENE
                    gene name <protease>
        -o ORGANISM, --organism=ORGANISM
                    organism: HIV, HCV <HIV>
        -r REFERENCE, --reference=REFERENCE
                    fasta file with reference

### Requirements
3. Python (**not** Python 3)
1. [Biopython](http://biopython.org/)
2. `needle` from [EMBOSS](http://emboss.sourceforge.net/) software suite

### Example
    [user@host amplicon_1]$ LStructure.py -s your_window-support.fas
    Support is your_window-support.fas
    Reference is protease from HIV
    Total reads: 1017.959

    [user@host amplicon_1]$ ls -1rht
    ...
    prot_seqs.fasta
    mutations_aa.csv
    mutations_DNA.csv
    dna_seqs.fasta
The files listed above are the output of LStructure
