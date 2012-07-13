#!/usr/bin/env python
'''This module parses support files, detects and list mutations
    '''
import sys
import os
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez

from Locals import HIV_gene_coord, HCV_gene_coord

MINIMUM_READ = 5.0

# If HIV file is not there, fetch it
Entrez.email = None  # To tell NCBI who you are
mod_dir = os.path.split(__file__)[0]
hiv_filename = os.path.join(mod_dir, "HIV-HXB2.fasta")
print hiv_filename
if not os.path.isfile(hiv_filename):
    # Downloading...
    handle = Entrez.efetch(db="nucleotide", id="1906382",
                           rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    SeqIO.write(seq_record, hiv_filename, 'fasta')
    print "Saved"
HXB2 = list(SeqIO.parse(hiv_filename, 'fasta'))[0]
HXB2.alphabet = IUPAC.ambiguous_dna

hcv_filename = os.path.join(mod_dir, "HCV1.fasta")
HCV = list(SeqIO.parse(hcv_filename, 'fasta'))[0]
HCV.alphabet = IUPAC.ambiguous_dna


translation_table = {  # 64 codons + '---'
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
    try:
        counts = [(translate(read[f:]).count('*'), f + 1) for f in range(3)]
    except Bio.Data.CodonTable.TranslationError:
        counts = [(gap_translation(read[f:]).count('*'), f + 1) \
                  for f in range(3)]
    sor_cnt = sorted(counts)
    stop_codons, frame = sor_cnt[0]
    if stop_codons > 0:
        warnings.warn('The sequence %s contains %dstop codons' \
                    % (read, stop_codons))
    if sor_cnt[1][0] == 0:
        warnings.warn('Two frames are possible!')
    return frame


def reframe(read, n_gaps):
    '''It transforms AC---GTGT or ACGT---GT into ACG---TGT if n_gaps=3
        Read must already be in frame
    '''
    gap_region = '-' * n_gaps
    # check that read is in frame, except for the
    # large gap region
    check_r = read.replace(gap_region, 'X' * n_gaps)
    # replace single gaps with C's, in order to avoid stop codons
    check_r = check_r.replace('-', 'C')
    check_r = check_r.replace('X', '')[:len(read) / 3 * 3]
    if gap_translation(check_r).count('*') > 0:
        warnings.warn('STOP codon found, seq=' + str(check_r))

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
        s = seq_in.tostring().upper()
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
        r_str = r_seq.tostring()
        Alignment.needle_align('asis:%s' % r_str,
                               'asis:%s' % self.seq.tostring(),
                               out_name, go=20.0, ge=0.1)
        tal = Alignment.alignfile2dict([out_name],
                                       'ref_mut_align', 20.0, 0.1)
        os.remove(out_name)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and lists mutations
        this.summary()
        m_start, m_stop = this.start, this.stop
        it_pair = zip(this.seq_a[m_start - 1:m_stop],
                      this.seq_b[m_start - 1:m_stop])
        mut_list = []
        for i, p in enumerate(it_pair):
            if p[0].upper() != p[1].upper():
                mut_list.append(Mutation(i + m_start, p[0], p[1]))
        self.mutations = mut_list


class LocalStructure:
    '''The main class here, takes the file, computes the consensus, the
        frame, the offset and list the variants'''
    def __init__(self, support_file, ref):
        import re

        self.sup_file = os.path.abspath(support_file)
        s_head, self.name = os.path.split(self.sup_file)

        descriptions = [s.description \
                        for s in SeqIO.parse(self.sup_file, 'fasta')]
        self.posteriors = \
            [float(re.search('posterior=(.*)\s*ave_reads=(.*)', d).group(1))
                for d in descriptions]
        self.ave_reads = \
            [float(re.search('posterior=(.*)\s*ave_reads=(.*)', d).group(2))
                for d in descriptions]
        #   self.ref_obj = HXB2
        #   self.gene = gene
        #   g_start, g_stop = gene_coord[self.gene]
        self.seq_obj = SeqIO.parse(self.sup_file, 'fasta')
        self.ref = ref  # HXB2[g_start - 1:g_stop]
        self.cons = Seq(self.get_cons(), IUPAC.unambiguous_dna)

        # try:
        self.frame = find_frame(self.cons)
        # except:
        # self.frame = 1
        # self.get_offset(HXB2)

        sup_far = os.path.join(s_head, '-'.join(self.name.split('-')[:-1]) \
                               + '.far')
        ns = SeqIO.parse(open(sup_far), 'fasta')
        self.n_reads = len([s for s in ns])
        self.dna_vars = []
        self.prot_vars = []

    def alignedvariants(self, threshold=0.9):
        '''Parses posterior, frequency, reframe large deletions,
        replaces single deletions and merges identical sequences,
        both at DNA and amino acids level. Returns a list of
        LocalVariant objects.
        '''
        import tempfile
        import numpy as np

        var_dict = {}
        aa_var_dict = {}
        for i, s in enumerate(self.seq_obj):
            post, ave_reads = self.posteriors[i], self.ave_reads[i]
            if post < threshold or ave_reads < MINIMUM_READ:
                continue
            if post > 1.0:
                print >> sys.stderr, 'WARNING: posterior=', post

            # reframe large deletions
            # frame is 1-based numbered
            read = s.seq.tostring()[self.frame - 1:]
            max_len = len(read.strip('-')) / 3 * 3
            for i in range(max_len, -1, -3):
                gap_region = '-' * i
                if gap_region in read.strip('-'):
                    read = reframe(read, i)
                    read = read.replace(gap_region, 'X' * i)
                    break

            # replaces single gaps with consensus
            this_seq = (p[1] if p[1] is not '-' else p[0] \
                        for p in zip(self.cons[self.frame - 1:], read))
            ws = ''.join(this_seq)
            ws = ws.replace('X', '-')
            ws_aa = gap_translation(ws)
            var_dict[ws] = var_dict.get(ws, 0) + ave_reads
            aa_var_dict[ws_aa] = aa_var_dict.get(ws_aa, 0) + ave_reads

        i = 1
        tot_freq = sum(var_dict.values())
        print >> sys.stderr, 'Total reads:', tot_freq
        for k, v in var_dict.items():
            freq_here = 100 * v / tot_freq
            tsr = LocalVariant(Seq(k, IUPAC.unambiguous_dna),
                               seq_id='reconstructed_hap_%d' % i,
                               description='Local hap freq=%f' % freq_here,
                               frequency=freq_here)
            self.dna_vars.append(tsr)
            i += 1

        i = 1
        tot_freq = sum(aa_var_dict.values())
        for k, v in aa_var_dict.items():
            freq_here = 100 * v / tot_freq
            trans_read = LocalVariant(Seq(k, IUPAC.protein), \
                                seq_id='translated_reconstructed_hap_%d' % i,
                                description='Local amino hap freq=%f' \
                                % freq_here,
                                frequency=freq_here)
            self.prot_vars.append(trans_read)
            i += 1

        # sort according to freq
        dna_freqs = [s.frequency for s in self.dna_vars]
        aa_freqs = [s.frequency for s in self.prot_vars]
        dna_argsort = np.argsort(dna_freqs)[::-1]
        aa_argsort = np.argsort(aa_freqs)[::-1]
        self.dna_vars = [self.dna_vars[i] for i in dna_argsort]
        self.prot_vars = [self.prot_vars[i] for i in aa_argsort]

        # When run in Xcode, working directory is root
        wd = os.getcwd()
        if wd == '/':
            wd = tempfile.gettempdir()
            print >> sys.stderr, 'Writing to directory', wd

        SeqIO.write(self.dna_vars, wd + '/dna_seqs.fasta', 'fasta')
        SeqIO.write(self.prot_vars, wd + '/prot_seqs.fasta', 'fasta')

    def print_mutations(self, rm_seq, seq_type='DNA',
                        out_format='csv', out_file=sys.stdout):
        '''As the name says
            '''

        from operator import itemgetter
        if out_format == 'csv':
            import csv
            fh = open(out_file, 'wb')
            writer = csv.writer(fh, dialect='excel')

        if seq_type == 'aa':
            to_print = self.prot_vars
        else:
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
        Alignment.needle_align('asis:%s' % self.ref.tostring(),
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
    '''Only tries optparse (deprecated in 2.7, in future add argparse)'''
    import optparse

    usage = "usage: %prog -s support_file [options]"
    optparser = optparse.OptionParser(usage=usage)

    optparser.add_option("-s", "--support", type="string", default="",
                         help="support file", dest="support")
    optparser.add_option("-g", "--gene", type="string", default="protease",
                         help="gene name", dest="gene")
    optparser.add_option("-o", "--organism", type="string",
                         help="organism: HIV, HCV", dest="organism")
    optparser.add_option("-r", "--reference", type="string", default="",
                         help="fasta file with reference", dest="reference")

    (opts, args) = optparser.parse_args()
    if not opts.support:
        optparser.error("specifying support file is mandatory")
    if opts.reference and opts.organism:
        optparser.error("options -r and -o are mutually exclusive")

    optparser.set_defaults(organism="HIV")
    (opts, args) = optparser.parse_args()

    return opts, args


if __name__ == '__main__':

    options, arguments = parse_com_line()
    sup_file = options.support
    print 'Support is', sup_file

    if options.reference:
        ref_rec = list(SeqIO.parse(options.reference, 'fasta'))[0]
        ref_seq = ref_rec.seq
        ref_seq_aa = ref_seq.translate()
        print >> sys.stderr, 'Reference is %s from file' % ref_rec.id
    elif options.organism == 'HIV':
        r_start, r_stop = HIV_gene_coord[options.gene]
        ref_seq = HXB2[r_start - 1:r_stop].seq
        ref_seq_aa = ref_seq.translate()
    elif options.organism == 'HCV':
        r_start, r_stop = HCV_gene_coord[options.gene]
        ref_seq = HCV[r_start - 1:r_stop].seq
        ref_seq_aa = ref_seq.translate()
    print >> sys.stderr, 'Reference is %s from %s' % \
                            (options.gene, options.organism)

    sample_ls = LocalStructure(support_file=sup_file, ref=ref_seq)

    sample_ls.alignedvariants(threshold=0.95)
    sample_ls.print_mutations(ref_seq, seq_type='DNA', out_format='csv',
                              out_file='mutations_DNA.csv')

    sample_ls.print_mutations(ref_seq_aa, seq_type='aa',
                              out_format='csv', out_file='mutations_aa.csv')
