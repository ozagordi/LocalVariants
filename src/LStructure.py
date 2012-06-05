#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

hdir = os.path.expanduser('~/')
DBdir = os.path.expanduser('~/Dropbox/Software')
sys.path.append(DBdir)

from pythonlib.HXB2_data import gene_coord

translation_table = { #64 codons + '---'
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
    from Bio.Seq import translate
    counts = [(translate(read[f:]).count('*'), f+1) for f in range(3)]
    sor_cnt = sorted(counts)
    stop_codons, frame = sor_cnt[0]#, key=operator.itemgetter(0))
    if stop_codons > 0: warnings.warn('The sequence %s contains %d stop codons'\
                                      % (read, stop_codons))
    if sor_cnt[1][0] == 0: warnings.warn('Two frames are possible!')
    return frame

def reframe(read, n_gaps):
    '''It transforms AC---GTGT or ACGT---GT into ACG---TGT if n_gaps=3
        Read must already be in frame
        '''
    gap_region = '-'*n_gaps
    # check that read is in frame, except for the
    # large gap region
    check_r = read.replace(gap_region, 'X'*n_gaps)
    # replace single gaps with C's, in order to avoid stop codons
    check_r = check_r.replace('-','C')
    check_r = check_r.replace('X','')[:len(read)/3*3]
    if gap_translation(check_r).count('*') > 0:
        warnings.warn('STOP codon found, seq=' + str(check_r))
    
    # move a single letter right or left, in order to have a
    # correct frame (min edit distance to the coding sequence)
    start = read.find(gap_region)
    stop = start+n_gaps
    if start%3 == 1:
        flank_left = read[:start-1]
        flank_right = read[start-1] + read[stop:]
    elif start%3 == 2:
        flank_left = read[:start] + read[stop]
        flank_right = read[stop+1:]
    elif start%3 == 0: # nothing to do
        return read
    return flank_left + gap_region + flank_right


def gap_translation(seq_in, frame=1):
    aa = []
    try:
        s = seq_in.tostring().upper()
    except:
        s = seq_in.upper()
    for i in range(frame-1, len(s), 3):
        try:
            aa.append(translation_table[s[i:i+3]])
        except:
            break
    return ''.join(aa)

def str_num_compare(x, y):
    return int(x) - int(y)

class LocalStructure:
    '''This class
        '''
    
    def __init__(self, sup_file, ref_file='References/HIV-HXB2.fasta',
    gene='protease'):
    # def __init__(self, sup_file, ref_file='Projects/V1V2/data/cons_all_TP.fasta',
    #   gene='env'):
        '''
            '''
        import os.path
        
        full_ref = os.path.join(hdir, ref_file)
        full_sup = os.path.abspath(sup_file)
        
        self.sup_file = full_sup
        self.ref_file = full_ref
        self.gene = gene
        
        s_head, self.name = os.path.split(full_sup)
        h = open(full_sup)
        self.seq_obj = SeqIO.parse(h, 'fasta')
        start, stop = gene_coord[self.gene]
        
        try:
            h = open(full_ref)
            sr = list(SeqIO.parse(h, 'fasta'))[0].seq.tostring().upper()
            self.ref = Seq(sr, IUPAC.unambiguous_dna)[start-1:stop]
        except IOError:
            self.ref = None
        
        self.cons = Seq(self.get_cons(), IUPAC.unambiguous_dna)
        try:
            self.frame = find_frame(self.cons)
        except:
            self.frame = 1
        self.get_offset()
        sup_far = os.path.join(s_head, '-'.join(self.name.split('-')[:-1]) + '.far')
        ns = SeqIO.parse(open(sup_far), 'fasta')
        self.n_reads = len([s for s in ns])
        
        self.dna_seqs = []
        self.prot_seqs = []
    
    def alignedvariants(self, threshold=0.9):
        import subprocess
        import re
        import itertools
        import hashlib
        import tempfile
        import numpy as np
        from pythonlib import Alignment
        
        files = []
        var_dict = {}
        aa_var_dict = {}
        for i, s in enumerate(self.seq_obj):
            m_obj = re.search('posterior=(.*)\s*ave_reads=(.*)', s.description)
            post, ave_reads = map(float, (m_obj.group(1), m_obj.group(2)))
            if post < threshold or ave_reads < 1.:
                continue
            if post > 1.0:
                print('WARNING: posterior=', post, file=sys.stderr)
            
            #this.summary()
            #start, stop = this.start, this.stop
            #it_pair = itertools.izip(this.seq_a[start-1:stop], this.seq_b[start-1:stop])        
            read = s.seq.tostring()[self.frame-1:]#this.seq_b
            max_len = len(read.strip('-'))/3*3
            for i in range(max_len, -1, -3):
                gap_region = '-'*i
                if gap_region in read.strip('-'):
                    read = reframe(read, i)
                    read = read.replace(gap_region, 'X'*i)
                    break
            
            it_pair = itertools.izip(self.cons, read)
            # p[0] is the consensus, p[1] the read
            this_seq = []
            while True:
                try:
                    p = it_pair.next()
                except StopIteration:
                    break
                if p is None:
                    break
                if p[1] == '-':
                    #assert p[0] != '-', 'gap-gap?'
                    this_seq.append(p[0])
                elif p[0] != '-':
                    this_seq.append(p[1])
            ws = ''.join(this_seq)
            ws = ws.replace('X', '-')
            ws_aa = gap_translation(ws)
            var_dict[ws] = var_dict.get(ws, 0) + ave_reads
            aa_var_dict[ws_aa] = aa_var_dict.get(ws_aa, 0) + ave_reads
        
        for k, v in var_dict.items():
            ts = Seq(k, IUPAC.unambiguous_dna)
            tsr = SeqRecord(ts, id = hashlib.sha224(k).hexdigest(), \
                            name='Reconstructed local hap')
            freq_here = 100*v/self.n_reads
            tsr.description = 'ave_freq=%f' % freq_here
            self.dna_seqs.append(tsr)
        
        for k, v in aa_var_dict.items():
            ts_translated = Seq(k, IUPAC.protein)
            trans_read = SeqRecord(ts_translated, id = hashlib.sha224(k).hexdigest(), \
                                   name='Translated local hap')
            freq_here = 100*v/self.n_reads
            trans_read.description = 'ave_freq=%f' % freq_here
            self.prot_seqs.append(trans_read)

        # sort according to freq
        dna_freqs = [float(s.description.split('=')[1]) for s in self.dna_seqs]
        aa_freqs = [float(s.description.split('=')[1]) for s in self.prot_seqs]
        dna_argsort = np.argsort(dna_freqs)[::-1]
        aa_argsort = np.argsort(aa_freqs)[::-1]
        tmp = [self.dna_seqs[i] for i in dna_argsort]
        self.dna_seqs = tmp
        tmp = [self.prot_seqs[i] for i in aa_argsort]
        self.prot_seqs = tmp

        # When run in Xcode, working directory is root
        wd = os.getcwd()
        if wd == '/':
            wd = tempfile.gettempdir()
            print('Writing to directory', wd, file=sys.stderr)
        
        SeqIO.write(self.dna_seqs, wd + '/dna_seqs.fasta', 'fasta')
        SeqIO.write(self.prot_seqs, wd + '/prot_seqs.fasta', 'fasta')
        
        return self.dna_seqs
    
    
    def get_cons(self, plurality=0.1, identity=1):
        '''Consensus by running EMBOSS cons and manipulation
            '''
        import subprocess
        import os
        import itertools
        import tempfile
        
        from pythonlib import Alignment
        
        # Compute the consensus with EMBOSS cons
        cline = 'cons -sequence %s -stdout -auto' % self.sup_file
        cline += ' -plurality %f -identity %i -name consensus' % (plurality, identity)
        p = subprocess.Popen(cline, shell=True, bufsize=1024, \
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                             close_fds=True)
        sc =  list(SeqIO.parse(p.stdout, 'fasta'))[0].seq.tostring().upper()
        strcons = sc#.replace('N', '')
        # Align the consensus to reference sequence
        outfile = tempfile.NamedTemporaryFile()
        out_name = outfile.name
        outfile.close()
        Alignment.needle_align(self.ref_file, 'asis:%s' % strcons, \
                               out_name, go=20.0, ge=0.1)
        tal = Alignment.alignfile2dict([out_name], 'ref_cons_alignment', 20.0, 0.1)
        os.remove(out_name)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        # Extracts only the matching region and fill gaps
        this.summary()
        start, stop = this.start, this.stop
        # ref_file, consensus
        it_pair = itertools.izip(this.seq_a[start-1:stop], this.seq_b[start-1:stop])
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
    
    def mytranslate(self):
        '''
            '''
        import re
        import hashlib
        from Bio.Seq import translate
        
        prot_dict = {}
        minstop = 10000
        for fr in range(3):
            s1 = self.dna_seqs[1].seq[fr:].translate()
            if minstop >=  s1.count('*'):
                self.frame = fr
                minstop = s1.count('*')
        
        print('Frame is', self.frame+1, file=sys.stderr)
        
        for s in self.dna_seqs:
            cod = Seq(s.seq[self.frame:].tostring())
            # aas = translate(cod).tostring()
            aas = gap_translation(cod)
            ave_reads = re.search('ave_reads=(.*)', s.description).group(1)
            prot_dict[aas] = prot_dict.get(aas, 0) + float(ave_reads)
        
        for k, v in prot_dict.items():
            ts = Seq(k, IUPAC.protein)
            tsr = SeqRecord(ts, id = hashlib.sha224(k).hexdigest(), \
                            name='Reconstructed local hap: translated')
            tsr.description = 'ave_reads=%f' % v
            self.prot_seqs.append(tsr)
        return self.prot_seqs
    
    
    def mutations(self, wh='DNA', out_format='csv', out_file=sys.stdout):
        '''
            '''
        from operator import itemgetter
        from pythonlib.Alignment import dna_code
        
        if wh == 'DNA': seqs = self.dna_seqs
        elif wh == 'aa': seqs = self.prot_seqs
        
        print('#' * 60, file=sys.stderr)
        print(str(' Now %s variants ' % wh).center(60, '#'), file=sys.stderr)
        print('#' * 60, file=sys.stderr)
        dna_offset = self.offset
        aa_offset = (self.offset+1)/3
        if wh == 'DNA': print('DNA offset is', dna_offset, file=sys.stderr)
        elif wh == 'aa': print('aa offset is', aa_offset, file=sys.stderr)
        all_mut = []
        mut_info = []
        # ref_cons_mut = []
        
        trans_ref = self.ref[0:300].translate()#[aa_offset:]
        print('Translated reference is:', file=sys.stderr)
        print(trans_ref, '\n', file=sys.stderr)
        
        # test position numbering        
        
        # parse the mutation in the haplotypes
        topr = True
        for s in seqs:
            mut =[]
            # frequencies are already in percentage
            # if lower than 0.5%, do not consider them
            freq_here = float(s.description.split('=')[1])
            if freq_here < 0.5: continue
            if wh == 'DNA':
                for i, p in enumerate(zip(self.ref[dna_offset-1:dna_offset+len(s)-1], s)):
                    if dna_code[p[0].upper()] & dna_code[p[1].upper()] == set([]):
                        this_mut = '%s%d%s' % (p[0].upper(), i+dna_offset+1 , p[1].upper())
                        if this_mut not in all_mut:
                            all_mut.append(this_mut)
                        mut.append(this_mut)
            if wh == 'aa':
                if topr:
                    print('This is to chech that reference matches query',
                          file=sys.stderr)
                    print(trans_ref[:10], file=sys.stderr)
                    print(s[:10].seq, file=sys.stderr)
                    print('', file=sys.stderr)
                    topr = False
                for i, p in enumerate(zip(trans_ref, s)):
                    if p[0].upper() != p[1].upper():
                        this_mut = '%s%d%s' % (p[0].upper(), i+1+aa_offset, p[1].upper())
                        if this_mut not in all_mut:
                            all_mut.append(this_mut)
                        mut.append(this_mut)
            mut_info.append([s.name, float(s.description.split('=')[1]), mut])
        
        mut_info = sorted(mut_info, key=itemgetter(1), reverse=True)
        all_mut = sorted(all_mut, key=itemgetter(slice(1,-1)), cmp=str_num_compare)
        tr_reads = sum([m[1] for m in mut_info])
        print('Haps account for %3.1f %% of the reads\n' % tr_reads, file=sys.stderr)
        #        print('That is %f%% of the total\n' % (100 * tr_reads/self.n_reads),
        #   file=sys.stderr)
        
        if out_format == 'human':
            for m in mut_info:
                print(m[0], 'ave_reads = %6.2f%%' % (100 * m[1]/tr_reads),
                      file=sys.stderr)
                for mm in m[2]:
                    print(mm, file=sys.stderr)
                print('\n', file=sys.stderr)
        
        elif out_format == 'csv':
            import csv
            fh = open(out_file, 'wb')
            #writer = csv.writer(fh, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer = csv.writer(fh, dialect='excel')#delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['freq\\muts'] + all_mut)
            for m in mut_info:
                freq_here = m[1]#(100 * m[1]/tr_reads)
                to_write =[('%4.2f' % freq_here)]
                for mm in all_mut:
                    if mm in m[2]:
                        to_write.append('X')
                    else:
                        to_write.append(' ')
                writer.writerow(to_write)
        
        elif out_format == 'tex':
            out_file.write('\\begin{sidewaystable}\n')
            out_file.write('\\centering\n')
            out_file.write('\\rowcolors{1}{Apricot}{cyan}\n')
            out_file.write('\\begin{tabular}{c*{%d}{|c}}\n' % len(all_mut))
            out_file.write('\\%& ')
            out_file.write(' &'.join(all_mut))
            out_file.write('\\\\\n')
            for m in mut_info:
                to_write =[('%4.2f' % (100 * m[1]/tr_reads))]
                for mm in all_mut:
                    if mm in m[2]:
                        to_write.append('X')
                    else:
                        to_write.append(' ')
                out_file.write(' &'.join(to_write))
                out_file.write('\\\\\n')
            out_file.write('\\end{tabular}\n')
            out_file.write(
                           '\\caption{Patient PR: haplotypes account for %f \\%% of the reads.}\n'
                           %  (100 * tr_reads/self.n_reads))
            out_file.write('\\end{sidewaystable}\n')
        
        return tr_reads, 100 * tr_reads/self.n_reads
    
    
    def get_offset(self, ref_file='~/References/HIV-HXB2.fasta'):
        
        from pythonlib import Alignment
        import os
        import tempfile
        
        outfile = tempfile.NamedTemporaryFile()
        out_name = outfile.name
        outfile.close()
        start, stop = gene_coord[self.gene]
        usa_seq = ref_file + '[%d:%d]' % (start, stop)
        Alignment.needle_align(usa_seq, 'asis:%s' % self.cons, out_name,\
                               go=10.0, ge=0.5)
        tal = Alignment.alignfile2dict([out_name], 'get_offset', 10.0, 0.5)
        os.remove(out_name)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        this.summary()
        self.offset = this.start
        print('Offset w.r.t', ref_file, 'is', self.offset, file=sys.stderr)
        return



def parse_com_line():
    import optparse
    
    optparser = optparse.OptionParser()
    
    optparser.add_option("-s","--support", type="string", default="", help="support file",
                         dest="support")
    optparser.add_option("-r","--ref", type="string", default="References/HIV-HXB2.fasta",
                         help="reference file", dest="ref")
    optparser.add_option("-g","--gene", type="string", default="protease", help="gene name",
                         dest="gene")
    
    (options, args) = optparser.parse_args()
    
    return options, args


if __name__ == '__main__':

    options, args = parse_com_line()

    ref_file = options.ref
    sup_file = options.support
    gene_name = options.gene
    
    sample_ls = LocalStructure(sup_file=sup_file, gene=gene_name, ref_file=ref_file)
    vd = sample_ls.alignedvariants(threshold=0.95)
    print('There are ', len(vd), ' DNA variants', file=sys.stderr)
    sample_ls.mutations(wh='DNA', out_format='csv', out_file='mutations_DNA.csv')
    # sample_ls.mutations(wh='aa', out_format='csv', out_file='mutations_aa.csv')
    '''
    n = 1
    while len(qualities) != len(nucleotides):
    nucleotides, nsub = re.subn(
    '[+-]%i[ATCGNatcgn]{%i}' % (n, n), '', nucleotides)
    n += 1
    if n > 100 or len(nucleotides) < len(qualities):
    sys.exit("Something is wrong: Less nucleotides than \
    qualities and/or insertions longer than 100nt!") #TRAP!
    '''
