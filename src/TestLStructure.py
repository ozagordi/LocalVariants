#!/usr/bin/env python
'''Write a few test for LStructure
    '''
import unittest
import os
from Bio import SeqIO
import LStructure as LS


class TestLocalStructure(unittest.TestCase):
    '''The test class'''
    def test_find_frame(self):
        '''As in the name. Frame is 1 based'''
        read = 'GTCGTCGTCCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTT'
        assert LS.find_frame(read) == 1, read
        read = 'AGTCGTCGTCCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTT'
        assert LS.find_frame(read) == 2, read
        read = 'TCGTCGTCCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTT'
        assert LS.find_frame(read) == 3, read

    def test_reframe(self):
        '''Move bases around gaps to improve translation'''
        read1 = 'GAAGCG---ACTACT'  # Already in frame
        assert LS.reframe(read1, 3) == read1
        read2 = 'GAAGC---GACTACT'  # Move one base <--
        assert LS.reframe(read2, 3) == read1
        read3 = 'GAAGCGA---CTACT'  # Move one base -->
        assert LS.reframe(read3, 3) == read1, LS.reframe(read3, 3)

    def test_gap_translation(self):
        '''Translation with gaps, frame dependent'''
        assert LS.gap_translation('AAG') == 'K'
        assert LS.gap_translation('AAG---') == 'K-'
        assert LS.gap_translation('AAGAA') == 'K'
        assert LS.gap_translation('CAAGA', frame=2) == 'K'

    def test_main_class(self):
        '''As in the name, test in general'''
        import tempfile
        seqs = ['GTCACTCTTTGGCAACGACCC',
                'GTC---CTTTGGCAACGACCC',  # Triple del: keep it as variant
                'GTCACTCTTTGGCAACTAC-C']  # Single del: should be corrected
        reads = [50, 30, 20]

        sup_file = tempfile.NamedTemporaryFile(suffix='-support.fas',
                                               delete=False)
        sup_name = sup_file.name
        i = 0
        for s, r in zip(seqs, reads):
            sup_file.write('>hap%d|posterior=1 ave_reads=%d\n' % (i, r))
            sup_file.write(s + '\n')
        sup_file.close()

        # far file is only used to count the total number of reads
        far_name = sup_name.split('-support')[0] + '.far'
        far_file = open(far_name, 'w+b')
        for i in range(sum(reads)):
            far_file.write('>read\nCAAAA\n')
        far_file.close()

        ref_seq = list(SeqIO.parse('HIV-HXB2.fasta', 'fasta'))[0].seq

        test_sample = LS.LocalStructure(support_file=sup_file.name,
                                        ref=ref_seq)

        start_cons = 'GTCACTCTTTGGCAACGACCC'
        assert test_sample.cons.tostring().upper().startswith(start_cons), \
            test_sample.cons
        assert test_sample.n_reads == sum(reads)
        for s in test_sample.dna_vars:
            print s.id, s.seq.tostring()

        os.remove(sup_name)
        os.remove(far_name)

        assert 1 == 0, 'Write more tests'

if __name__ == '__main__':
    unittest.main()
