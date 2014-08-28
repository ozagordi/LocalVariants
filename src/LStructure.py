#!/usr/bin/env python
'''This module parses support files, detects and list mutations
    '''
import sys
import os
import warnings

import logging
import logging.handlers

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez

# Make a global logging object.
lslog = logging.getLogger(__name__)

MINIMUM_READ = 5.0

# If HIV file is not there, fetch it
Entrez.email = None  # To tell NCBI who you are
mod_dir = os.path.split(__file__)[0]
hiv_filename = os.path.join(mod_dir, "HIV-HXB2.fasta")

if not os.path.isfile(hiv_filename):
    # Downloading...
    print >> sys.stderr, 'Downloading HIV reference from PubMed'
    handle = Entrez.efetch(db="nucleotide", id="1906382",
                           rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    SeqIO.write(seq_record, hiv_filename, 'fasta')
    print >> sys.stderr, 'Saved to', hiv_filename
HXB2 = list(SeqIO.parse(hiv_filename, 'fasta'))[0]
HXB2.alphabet = IUPAC.ambiguous_dna

hcv_filename = os.path.join(mod_dir, "HCV1.fasta")
HCV = list(SeqIO.parse(hcv_filename, 'fasta'))[0]
HCV.alphabet = IUPAC.ambiguous_dna

# 64 codons + '---'
translation_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*', '---': '-'}


def find_frame(read):
    '''Frame is the one with the smallest number of stop codons
    '''
    from Bio.Seq import translate
    import Bio

    # use this to cut read at multiple of three length
    rem = len(read) % 3
    last_pos = rem if rem else None
    try:
        read = read[:-last_pos]
    except TypeError:
        pass    
    assert len(read) % 3 == 0, read
    read_len = len(read) - 3
    try:
        counts = [(translate(read[f:read_len + f]).count('*'), f + 1)
                  for f in range(3)]
    except Bio.Data.CodonTable.TranslationError:
        counts = [(gap_translation(read[f:]).count('*'), f + 1)
                  for f in range(3)]
    sor_cnt = sorted(counts)
    stop_codons, frame = sor_cnt[0]
    if stop_codons > 0:
        warnings.warn('The sequence %s contains %dstop codons'
                      % (read, stop_codons))
    if sor_cnt[1][0] == 0:
        warnings.warn('Two frames are possible! %d and %d' %
                      (frame, sor_cnt[1][1]))
    return frame


def reframe(read, n_gaps):
    '''It transforms AC---GTGT or ACGT---GT into ACG---TGT if n_gaps=3
        Read must already be in frame
    '''
    import re
    gap_region = '-' * n_gaps
    # check that read is in frame, except for the
    # large gap region
    check_r = read.replace(gap_region, 'X' * n_gaps)
    # replace single gaps with C's, in order to avoid stop codons
    check_r = check_r.replace('-', 'C')
    check_r = check_r.replace('X', '')[:len(read) / 3 * 3]
    if gap_translation(check_r).count('*') > 0:
        warnings.warn('STOP codon found, seq=' + str(check_r))

    # first the case with 2 + 1 gaps
    # CAG--G-AC
    # CAG-G--AC
    # CA-G--GAC
    match_21 = re.search('\w*(--)\w*(-)\w*', read)
    if match_21:
        first_gap = match_21.start(1)
        print first_gap
        print read
        print read[:first_gap]
        sys.exit()

    # move a single letter right or left, in order to have a
    # correct frame (min edit distance to the coding sequence)
    start = read.find(gap_region)

    stop = start + n_gaps
    if start % 3 == 1:
        flank_left = read[:start - 1]
        flank_right = read[start - 1] + read[stop:]
    elif start % 3 == 2:
        flank_left = read[:start] + read[stop]
        flank_right = read[stop + 1:]
    elif start % 3 == 0:  # nothing to do
        return read
    return flank_left + gap_region + flank_right


def gap_translation(seq_in, frame=1):
    '''Biopython translation does not accept gaps
        '''
    aa = []
    try:
        s = str(seq_in).upper()
    except AttributeError:
        s = seq_in.upper()
    for i in range(frame - 1, len(s), 3):
        try:
            aa.append(translation_table[s[i:i + 3]])
        except KeyError:
            break
    return ''.join(aa)


def str_num_compare(x, y):
    '''Used to sort mutations'''
    return int(x) - int(y)


class Mutation:
    '''Simple container for the mutation'''
    def __init__(self, pos, wt, mut, comment=None):
        self.mutated = mut  # mutatis (to)
        self.original = wt  # mutandis (from)
        self.position = pos
        self.comment = comment

    def __str__(self):
        if self.original:
            return '%s%d%s' % (self.original, self.position,
                               self.mutated)
        else:
            return '%d%s' % (self.position, self.mutated)


class LocalVariant(SeqRecord):
    '''Built on SeqRecord, this adds a few useful attributes'''
    def __init__(self, seq, seq_id, **kwargs):
        SeqRecord.__init__(self, seq=seq, id=seq_id)
        try:
            self.frame = find_frame(seq)
        except ValueError:
            self.frame = None
        self.frequency = None
        self.mutations = []
        for k, v in kwargs.items():
            self.__dict__[k] = v

    def get_mutations(self, r_seq):
        '''Given the variant sequence and the reference,
            aligns them and detects the mutations'''
        import tempfile
        import Alignment

        outfile = tempfile.NamedTemporaryFile()
        out_name = outfile.name
        outfile.close()
        r_str = str(r_seq)
        Alignment.needle_align('asis:%s' % r_str,
                               'asis:%s' % str(self.seq),
                               out_name, go=20.0, ge=2.0)
        tal = Alignment.alignfile2dict([out_name],
                                       'ref_mut_align', 20.0, 2.0)
        os.remove(out_name)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and lists mutations
        this.summary()
        m_start, m_stop = this.start, this.stop
        lslog.info('start: %d  stop: %d' % (m_start, m_stop))

        # Check if alignment starts before or after the reference
        if not this.seq_a.startswith('-') and this.seq_b.startswith('-'):
            # after
            it_pair = zip(r_str[m_start - 1:m_stop],
                          str(self.seq))
        elif this.seq_a.startswith('-') and not this.seq_b.startswith('-'):
            # before
            it_pair = zip(r_str[m_start - 1:m_stop],
                          str(self.seq)[m_start - 1:m_stop])
        elif not this.seq_a.startswith('-') and not this.seq_b.startswith('-'):
            # together
            it_pair = zip(r_str,
                          str(self.seq))
            
        seq_dist = sum(p[0] != p[1] for p in it_pair)
        lslog.info('distance between ref and seq: %d' % seq_dist)
        if seq_dist > 20:
            lslog.warning('distance > 20')
            # print ''.join(p[0] for p in it_pair)
            # print ''.join(p[1] for p in it_pair)

        mut_list = []
        for i, p in enumerate(it_pair):
            # print >> sys.stderr, '%s%d%s' % (p[0], i + m_start, p[1])
            if p[0].upper() != p[1].upper():
                mut_list.append(Mutation(i + m_start, p[0], p[1]))
        self.mutations = mut_list


class LocalStructure:
    '''The main class here, takes the file, computes the consensus, the
        frame, the offset and list the variants'''
    def __init__(self, support_file, ref):
        import re
        import glob

        self.sup_file = os.path.abspath(support_file)
        s_head, self.name = os.path.split(self.sup_file)

        descriptions = [s.description
                        for s in SeqIO.parse(self.sup_file, 'fasta')]
        prog = re.compile('posterior=(\d*.*\d*)[-\s\t]ave_reads=(\d*.*\d*)')
        self.posteriors = \
            [float(prog.search(d).group(1)) for d in descriptions]
        self.ave_reads = \
            [float(prog.search(d).group(2)) for d in descriptions]
        lslog.info('parsed %d sequences' % len(descriptions))

        self.seq_obj = SeqIO.parse(self.sup_file, 'fasta')
        self.ref = ref  # HXB2[g_start - 1:g_stop]
        self.cons = Seq(self.get_cons(), IUPAC.unambiguous_dna)
        self.frame = find_frame(self.cons)
        lslog.info('consensus starts with %s; frame is %d' %
                   (self.cons[:12], self.frame))

        # shorah 0.6 introduced strand bias correction to improve precision
        # in  SNVs calling. Results are in file SNVs_*_final.csv.
        snv_files = glob.glob('*_final.csv')
        assert len(snv_files) == 1
        snv_file = snv_files[0]
        csv_reader = open(snv_file)
        csv_reader.next()
        snvs = [row.split(',')[2] + row.split(',')[1] + row.split(',')[3]
                for row in csv_reader]
        self.snvs = snvs
        lslog.info('%d SNVs found in %s' % (len(snvs), snv_file))

        # now parsing the total number of reads
        try:
            sup_far = os.path.join(s_head, '-'.join(self.name.split('-')[:-1])
                                   + '.far')
            ns = SeqIO.parse(open(sup_far), 'fasta')
            self.ds = len([s for s in ns])
            lslog.info('%d reads; found parsing far.file' % self.ds)
        except IOError:
            cov_h = open('coverage.txt')
            for cov_line in cov_h:
                if cov_line.split()[0].split('.')[0] == \
                        self.name.split('.')[0]:
                    self.n_reads = int(cov_line.split()[4])
                    break
            cov_h.close()
            lslog.info('%d reads; found in coverage.txt' % self.n_reads)
        if not self.n_reads:
            lslog.error('coverage not parsed, n of reads unknown')
            sys.exit('coverage not parsed, n of reads unknown')
        self.dna_vars = []

    def alignedvariants(self, threshold=0.9):
        '''Merges sequences identical at DNA level. Returns a list of
        LocalVariant objects.
        '''

        import numpy as np

        var_dict = {}
        for i, s in enumerate(self.seq_obj):
            post, ave_reads = self.posteriors[i], self.ave_reads[i]
            if post < threshold or ave_reads < MINIMUM_READ:
                continue
            if post > 1.0:
                print >> sys.stderr, 'WARNING: posterior=', post
                lslog.warning('posterior=%f' % post)
            ws = str(s.seq)
            var_dict[ws] = var_dict.get(ws, 0) + ave_reads

        tot_freq = sum(var_dict.values())
        print >> sys.stderr, 'Total reads:', tot_freq

        supported_var_dict = {}  # supported by SNVs
        # first pass excludes the variants unsupported in SNVs final file
        lslog.info(self.snvs)
        for k, v in var_dict.items():
            tsr = LocalVariant(Seq(k, IUPAC.unambiguous_dna),
                               seq_id='reconstructed_hap',
                               description='to_be_confirmed',
                               frequency=0.0)

            tsr.get_mutations(self.ref)
            save = True
            lslog.info('Checking mutations on %s' % k)
            for mt in tsr.mutations:
                if str(mt) not in self.snvs:
                    save = False
                    lslog.debug('%s not supported' % str(mt))
            if save:
                supported_var_dict[k] = v

        if len(supported_var_dict) == 0:
            raise Exception("No supported variants!")

        i = 1
        tot_freq = sum(supported_var_dict.values())
        for k, v in supported_var_dict.items():
            freq_here = 100 * v / tot_freq
            tsr = LocalVariant(Seq(k, IUPAC.unambiguous_dna),
                               seq_id='reconstructed_hap_%d' % i,
                               description='Local nt haplo freq=%2.1f'
                               % freq_here,
                               frequency=freq_here)
            self.dna_vars.append(tsr)
            i += 1

        # sort according to freq
        dna_freqs = [s.frequency for s in self.dna_vars]
        dna_argsort = np.argsort(dna_freqs)[::-1]
        self.dna_vars = [self.dna_vars[i] for i in dna_argsort]

        wd = os.getcwd()
        SeqIO.write(self.dna_vars, wd + '/dna_seqs.fasta', 'fasta')
        lslog.info('Sequences written to file')


    def print_mutations(self, rm_seq, out_format='csv', out_file=sys.stdout):
        '''As the name says
        '''

        from operator import itemgetter
        if out_format == 'csv':
            import csv
            fh = open(out_file, 'wb')
            writer = csv.writer(fh, dialect='excel')

        to_print = self.dna_vars
        all_mut = set([])
        for v in to_print:
            v.get_mutations(rm_seq)
            for m in v.mutations:
                all_mut.add(str(m))
        all_mut = sorted(list(all_mut), key=itemgetter(slice(1, -1)),
                         cmp=str_num_compare)
        all_mut.insert(0, 'freq\\mut')

        if out_format == 'csv':
            writer.writerow(all_mut)
        elif out_format == 'human':
            print '\t'.join(all_mut)

        for v in to_print:
            mut_str = [str(m) for m in v.mutations]
            row = ['X' if put in mut_str else ' ' for put in all_mut[1:]]
            row.insert(0, round(v.frequency, 1))
            if out_format == 'csv':
                writer.writerow(row)
            elif out_format == 'human':
                print '\t'.join(map(str, row))

    def get_cons(self):
        '''Consensus by weighting the support with frequencies
        '''
        import tempfile
        from itertools import izip
        from collections import Counter
        import Alignment
        from Bio import AlignIO

        alignment = AlignIO.read(self.sup_file, 'fasta')
        sc = []
        for i in range(len(alignment[0])):
            bases = alignment[:, i]
            c = Counter()
            for i, b in enumerate(bases):
                c[b] += int(round(self.posteriors[i] * self.ave_reads[i]))
            sc.append(c.most_common()[0][0].upper())
        strcons = ''.join(sc)
        # Align the consensus to reference sequence
        outfile = tempfile.NamedTemporaryFile()
        out_name = outfile.name
        outfile.close()
        Alignment.needle_align('asis:%s' % str(self.ref),
                               'asis:%s' % strcons,
                               out_name, go=20.0, ge=0.1)
        tal = Alignment.alignfile2dict([out_name],
                                       'ref_cons_alignment', 20.0, 0.1)
        os.remove(out_name)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and fill gaps
        this.summary()
        start, stop = this.start, this.stop
        # ref_file, consensus
        it_pair = izip(this.seq_a[start - 1:stop],
                       this.seq_b[start - 1:stop])
        this_seq = []
        while True:
            try:
                p = it_pair.next()
            except StopIteration:
                break
            if p is None:
                break
            if p[1] == '-' or p[1] == 'N':
                #assert p[0] != '-', 'gap-gap?'
                this_seq.append(p[0])
            elif p[0] != '-':
                this_seq.append(p[1])
        return ''.join(this_seq)


def parse_com_line():
    '''Standard option parsing'''
    options_l, args_l = None, None

    try:
        import argparse

        # default action is 'store'
        parser = argparse.ArgumentParser(description='Local structure',
                                         epilog='Input are mandatory')
        parser.add_argument('-s', '--support', dest='support',
                            help='support file')
        parser.add_argument('-r', '--reference', dest='reference',
                            help='fasta file with the reference used in \
                            running shorah')
        args_l = parser.parse_args()

    except ImportError:
        import optparse

        usage = "usage: %prog -s support_file -r reference"
        optparser = optparse.OptionParser(usage=usage)

        optparser.add_option("-s", "--support", type="string", default="",
                             help="support file", dest="support")
        optparser.add_option("-r", "--reference", type="string", default="",
                             help="fasta file with the reference used in \
                             running shorah", dest="reference")
        (options_l, args_l) = optparser.parse_args()

        if not options_l.support:
            optparser.error("specifying support file is mandatory")

    return options_l, args_l


if __name__ == '__main__':

    options, args = parse_com_line()
    if args:
        args = vars(args)
    else:
        args = vars(options)

    # set logging level
    lslog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './locstr.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=200000, backupCount=5)
    fo = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    hl.setFormatter(fo)
    lslog.addHandler(hl)
    lslog.info(' '.join(sys.argv))

    sup_file = args['support']
    print 'Support is', sup_file

    if args['reference']:
        ref_rec = list(SeqIO.parse(args['reference'], 'fasta'))[0]
        ref_seq = ref_rec.seq
        #ref_seq_aa = ref_seq.translate()
        print >> sys.stderr, 'Reference is %s from file' % ref_rec.id

    sample_ls = LocalStructure(support_file=sup_file, ref=ref_seq)

    sample_ls.alignedvariants(threshold=0.95)

    sample_ls.print_mutations(ref_seq, out_format='csv',
                              out_file='mutations_DNA.csv')
