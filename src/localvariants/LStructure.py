#!/usr/bin/env python3
"""This module parses support files, detects and list mutations."""
import logging
import logging.handlers
import os
import sys
import warnings

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Make a global logging object.
lslog = logging.getLogger(__name__)

MINIMUM_READ = 5.0

# If HIV file is not there, fetch it
mod_dir = os.path.split(__file__)[0]
hiv_filename = os.path.join(mod_dir, "db/HIV-HXB2.fasta")

HXB2 = list(SeqIO.parse(hiv_filename, 'fasta'))[0]
HXB2.alphabet = IUPAC.ambiguous_dna

hcv_filename = os.path.join(mod_dir, "db/HCV1.fasta")
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


def get_region_limits(ref_seq, cons_seq):
    """Cut supported haplotypes according to the coordinates given in region.

    This function takes just the coonsensus and reference, already
    cut according to region chr:start-stop, and returns the coordinate to cut
    consensus and, consequently, all sequences in support.
    """
    import re
    from tempfile import NamedTemporaryFile
    from . import Alignment

    # align with
    outfile = NamedTemporaryFile()
    out_name = outfile.name
    outfile.close()
    print('>cons', file=sys.stderr)
    print(str(cons_seq), file=sys.stderr)
    print('>ref_seq', file=sys.stderr)
    print(str(ref_seq), file=sys.stderr)

    Alignment.needle_align('asis:%s' % str(cons_seq),
                           'asis:%s' % str(ref_seq),
                           out_name, go=20.0, ge=2.0)
    tal = Alignment.alignfile2dict([out_name],
                                   'ref_cons_align', 20.0, 2.0)
    # os.remove(out_name)
    ka = list(tal.keys())[0]
    this = tal[ka]['asis']
    # Extracts only the matching region
    this.summary()

    # check if region starts before or after consensus
    if not this.seq_a.startswith('-') and this.seq_b.startswith('-'):
        # after
        saved_start = this.start
    else:
        saved_start = 0
    # check if region ends before or after consensus
    if this.seq_a.endswith('-') and not this.seq_b.endswith('-'):
        # after
        saved_stop = -1
    else:
        saved_stop = this.stop

    match_start = re.search(str(cons_seq[saved_start:saved_start + 10]), str(cons_seq))
    match_end = re.search(str(cons_seq[saved_stop - 10:saved_stop]), str(cons_seq))

    return match_start.start(), match_end.end()


def find_frame(read):
    """Frame is the one with the smallest number of stop codons."""
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
    """Transform AC---GTGT or ACGT---GT into ACG---TGT if n_gaps=3.

    Read must already be in frame
    """
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
    match_21 = re.search(r'\w*(--)\w*(-)\w*', read)
    if match_21:
        first_gap = match_21.start(1)
        print(first_gap)
        print(read)
        print(read[:first_gap])
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
    """Biopython translation does not accept gaps."""
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
    """Use this to sort mutations."""
    return int(x[1:-1]) - int(y[1:-1])


class Mutation:
    """Simple container for the mutation."""

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
    """Building on SeqRecord, add a few useful attributes."""

    def __init__(self, seq, seq_id, **kwargs):
        SeqRecord.__init__(self, seq=seq, id=seq_id)
        try:
            self.frame = find_frame(seq)
        except ValueError:
            self.frame = None
        self.frequency = None
        self.mutations = []
        for k, v in list(kwargs.items()):
            self.__dict__[k] = v

    def get_mutations(self, r_seq):
        """Given the variant sequence and the reference, aligns them and detects the mutations."""
        from tempfile import NamedTemporaryFile
        from . import Alignment

        outfile = NamedTemporaryFile()
        out_name = outfile.name
        outfile.close()
        r_str = str(r_seq)
        Alignment.needle_align('asis:%s' % r_str,
                               'asis:%s' % str(self.seq),
                               out_name, go=20.0, ge=2.0)
        tal = Alignment.alignfile2dict([out_name],
                                       'ref_mut_align', 20.0, 2.0)

        os.remove(out_name)
        ka = list(tal.keys())[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and lists mutations
        this.summary()
        m_start, m_stop = this.start, this.stop
        lslog.info('start: %d  stop: %d', m_start, m_stop)

        # Check if alignment starts before or after the reference
        if not this.seq_a.startswith('-') and this.seq_b.startswith('-'):
            # after
            it_pair = list(zip(r_str[m_start - 1:m_stop], str(self.seq)))
        elif this.seq_a.startswith('-') and not this.seq_b.startswith('-'):
            # before
            it_pair = list(zip(r_str[m_start - 1:m_stop], str(self.seq)[m_start - 1:m_stop]))
        elif not this.seq_a.startswith('-') and not this.seq_b.startswith('-'):
            # together
            it_pair = list(zip(r_str, str(self.seq)))

        seq_dist = sum(p[0].upper() != p[1].upper() for p in it_pair)
        lslog.info('distance between ref and seq: %d', seq_dist)
        if seq_dist > 20:
            lslog.warning('distance is %d > 20', seq_dist)
            # print ''.join(p[0] for p in it_pair)
            # print ''.join(p[1] for p in it_pair)

        mut_list = []
        for i, p in enumerate(it_pair):
            # print >> sys.stderr, '%s%d%s' % (p[0], i + m_start, p[1])
            if p[0].upper() != p[1].upper():
                mut_list.append(Mutation(i + m_start, p[0], p[1]))
        self.mutations = mut_list


class LocalStructure:
    """This is the main class.

    Takes the file, computes the consensus, the frame, the offset and list the variants.
    """

    def __init__(self, support_file, ref, region):
        import re
        import glob

        self.sup_file = os.path.abspath(support_file)
        s_head, self.name = os.path.split(self.sup_file)
        self.region = region
        descriptions = [s.description
                        for s in SeqIO.parse(self.sup_file, 'fasta')]
        prog = re.compile(r'posterior=(\d*.*\d*)[-\s\t]ave_reads=(\d*.*\d*)')
        self.posteriors = \
            [float(prog.search(d).group(1)) for d in descriptions]
        self.ave_reads = \
            [float(prog.search(d).group(2)) for d in descriptions]
        lslog.info('parsed %d sequences', len(descriptions))

        self.seq_obj = SeqIO.parse(self.sup_file, 'fasta')
        ref_rec = list(SeqIO.parse(ref, 'fasta'))[0]
        if self.region:
            assert ref_rec.id == self.region.split(':')[0]
            reg_start, reg_stop = [int(a) for a in self.region.split(':')[1].split('-')]
            self.ref = ref_rec.seq[reg_start:reg_stop]
        else:
            self.ref = ref_rec.seq
        # ref_seq_aa = ref_seq.translate()
        print('Reference is %s from file' % ref_rec.id, file=sys.stderr)
        # self.get_cons uses self.ref, so if the region is specified only the
        # corresponding consensus will be extracted
        gc_seq, gc_start, gc_stop = self.get_cons()
        self.cons = Seq(gc_seq, IUPAC.unambiguous_dna)
        if self.region:
            self.region_limits = gc_start, gc_stop  # get_region_limits(self.ref, self.cons)
            lslog.info('region specified, limits are start: %d  stop: %d', gc_start, gc_stop)

        self.frame = find_frame(self.cons)
        lslog.info('consensus starts with %s; frame is %d', self.cons[:12], self.frame)

        # shorah 0.6 introduced strand bias correction to improve precision
        # in  SNVs calling. Results are in file SNVs_*_final.csv.
        self.snvs = None
        snv_files = glob.glob('*_final.csv')
        assert len(snv_files) <= 1
        if len(snv_files) == 1:
            snv_file = snv_files[0]
            csv_reader = open(snv_file)
            next(csv_reader)
            snvs = [row.split(',')[2] + row.split(',')[1] + row.split(',')[3]
                    for row in csv_reader]
            self.snvs = snvs
            lslog.info('%d SNVs found in %s', len(snvs), snv_file)

        # now parsing the total number of reads from dbg file
        self.n_reads = 0
        assert len(glob.glob('*.dbg')) == 1
        dbg_file = glob.glob('*.dbg')[0]
        # dbg_file = os.path.join(s_head,
        #                         '-'.join(self.name.split('-')[:-1]) + '.dbg')
        with open(dbg_file) as f:
            for l in f:
                mobj = re.search(r'Number of reads, n = (\d*)', l)
                if mobj:
                    self.n_reads = int(mobj.group(1))

            lslog.info('%d reads; found parsing dbg file', self.n_reads)
        if not self.n_reads:
            lslog.error('coverage not parsed, n of reads unknown')
            sys.exit('coverage not parsed, n of reads unknown')
        self.dna_vars = []

    def alignedvariants(self, threshold=0.9):
        """Merge sequences identical at DNA level.

        Returns a list of LocalVariant objects.
        """
        import numpy as np

        var_dict = {}
        for i, s in enumerate(self.seq_obj):
            post, ave_reads = self.posteriors[i], self.ave_reads[i]
            if post < threshold or ave_reads < MINIMUM_READ:
                continue
            if post > 1.0:
                print('WARNING: posterior=', post, file=sys.stderr)
                lslog.warning('posterior=%f', post)
            if self.region:
                a, b = self.region_limits
                ws = str(s.seq[a:b])
            else:
                ws = str(s.seq)
            var_dict[ws] = var_dict.get(ws, 0) + ave_reads

        tot_freq = sum(var_dict.values())
        print('Total reads:', tot_freq, file=sys.stderr)

        supported_var_dict = {}  # supported by SNVs
        # first pass excludes the variants unsupported in SNVs final file
        lslog.info(self.snvs)
        for k, v in list(var_dict.items()):
            tsr = LocalVariant(Seq(k, IUPAC.unambiguous_dna), seq_id='reconstructed_hap',
                               description='to_be_confirmed', frequency=0.0)

            tsr.get_mutations(self.ref)
            save = True
            lslog.info('Checking mutations on %s', k)
            for mt in tsr.mutations:
                if self.snvs and str(mt) not in self.snvs:
                    save = False
                    lslog.debug('%s not supported', str(mt))
            if save:
                supported_var_dict[k] = v

        if supported_var_dict == {}:
            raise Exception("No supported variants!")

        i = 1
        tot_freq = sum(supported_var_dict.values())
        for k, v in list(supported_var_dict.items()):
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

    def print_mutations(self, rmut, out_format='csv', out_file=sys.stdout):
        """Do as the name says."""
        from functools import cmp_to_key
        if out_format == 'csv':
            import csv
            fh = open(out_file, 'wt')
            writer = csv.writer(fh, dialect='excel')

        rm_seq = list(SeqIO.parse(rmut, 'fasta'))[0].seq
        to_print = self.dna_vars
        all_mut = set([])
        for v in to_print:
            v.get_mutations(rm_seq)
            for m in v.mutations:
                all_mut.add(str(m))
        all_mut = sorted(list(all_mut), key=cmp_to_key(str_num_compare))
        all_mut.insert(0, 'freq\\mut')

        if out_format == 'csv':
            writer.writerow(all_mut)
        elif out_format == 'human':
            print('\t'.join(all_mut))

        for v in to_print:
            mut_str = [str(m) for m in v.mutations]
            row = ['X' if put in mut_str else ' ' for put in all_mut[1:]]
            row.insert(0, round(v.frequency, 1))
            if out_format == 'csv':
                writer.writerow(row)
            elif out_format == 'human':
                print('\t'.join(map(str, row)))

    def get_cons(self):
        """Compute consensus by weighting the support with frequencies."""
        import re
        import tempfile

        from collections import Counter
        from . import Alignment
        from Bio import AlignIO

        alignment = AlignIO.read(self.sup_file, 'fasta')
        sc = []
        for j in range(len(alignment[0])):
            bases = alignment[:, j]
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
        ka = list(tal.keys())[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and fill gaps
        this.summary()
        start, stop = this.start, this.stop
        # ref_file, consensus
        it_pair = zip(this.seq_a[start - 1:stop],
                      this.seq_b[start - 1:stop])
        this_seq = []
        while True:
            try:
                p = next(it_pair)
            except StopIteration:
                break
            if p is None:
                break
            if p[1] == '-' or p[1] == 'N':
                # assert p[0] != '-', 'gap-gap?'
                this_seq.append(p[0])
            elif p[0] != '-':
                this_seq.append(p[1])

        ts_str = ''.join(this_seq)
        match_start = re.search(ts_str[:10], strcons)
        match_end = re.search(ts_str[-10:], strcons)
        return ts_str, match_start.start(), match_end.end()
