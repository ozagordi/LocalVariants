#!/usr/bin/env python
"""Builds upon bam2msa from ShoRAH package.
    Takes sff or fastq and aligns with an external
    mapper into SAM format. Transforms into BAM,
    sorts and index it. Then it parses the BAM file
    to build meaningful MSA
    """

import sys
from Locals import aligner

SANGER_SCORE_OFFSET = ord("!")
Q_MAPPING = {}
for letter in range(0, 255):
    Q_MAPPING[chr(letter)] = letter - SANGER_SCORE_OFFSET
QC = 6
STRAND = ['+', '-']


def getinsertions(cigar, rpos):
    """Parses all insertions that will then be used
       to build the MSA (gap propagation)
    """

    assert cigar[0][0] == 4 and cigar[-1][0] == 4, 'These must be soft clipped'
    # Initial and final soft clipping is skipped, as read.pos marks
    # the start of the alignment excluding the soft clipping
    ins_position = {}
    for op_cigar, count in cigar[1:-1]:
        if op_cigar == 0:  # Alignment match
            rpos += count
        elif op_cigar == 1:  # Insertion
            ins_position[rpos] = count
            rpos += count
        elif op_cigar == 2:  # Deletion
            pass
        else:  # Something bad happened
            print 'This is bad'

    return ins_position


def getseq(alignedread, start=0, stop=0, **kwargs):
    """Retrieve the sequence of an alignedread object between start
       and stop of the reference position.
       The output will be padded by N's if the region exceeds the read length.
        """

    if alignedread.is_unmapped:
        return
    qcutoff = kwargs['qcutoff']
    ins_position = kwargs.get('ins_position', None)
    max_ins_per_pos = kwargs.get('max_ins_per_pos', None)

    seq = [s if Q_MAPPING[alignedread.qual[i]] >= qcutoff else 'N'\
            for i, s in enumerate(alignedread.seq)]
    # seq still contains soft_clipped bases

    if alignedread.cigar == None:
        return None

    assert alignedread.cigar[0][0] == 4 and alignedread.cigar[-1][0] == 4
    assert alignedread.qstart == alignedread.cigar[0][1]
    assert alignedread.pos + alignedread.alen == alignedread.aend
    assert alignedread.qstart + alignedread.qlen == alignedread.qend
    # should hold ref[pos:pos+10] aligned to read[qstart:qstart+10]
    # and ref[aend-10:aend] aligned to read[qend-10:qend]

    rpos = 0  # position in the read
    apos = alignedread.pos  # position in the alignment
    #gaps = 0  # number of gaps added
    fasta = []  # will hold the aligned read
    print alignedread.qname, apos, alignedread.seq[11:20]
    codes = []
    for op_cigar, count in alignedread.cigar[:-1]:  # Skip final soft clipping
        codes.extend([op_cigar] * count)

    for op_cigar, count in alignedread.cigar[:-1]:  # Skip final soft clipping
        if op_cigar == 0:  # Alignment match
            fasta.extend(seq[rpos:(rpos + count)])
            rpos += count
        elif op_cigar == 1:  # Insertion
            fasta.extend(seq[rpos:(rpos + count)])
            rpos += count
            print >> sys.stderr, 'Doing insertion'
        elif op_cigar == 2:  # Deletion
            fasta.extend(['-'] * count)
        elif op_cigar == 4:  # Soft clipped base
            rpos += count
        else:  # Something bad happened
            print 'This is bad'

    for k, v in max_ins_per_pos.items():
        already = ins_position.get(k, 0)
        if v > already:
            fasta.insert(k - apos, '-' * (v - already))
    # Compute range to output
    begin = max(0, start - alignedread.pos)
    fasta = fasta[begin:]
    fasta = ''.join(fasta)
    # Pad the ends
    fasta += 'n' * (stop - start - len(fasta))

    return fasta


def bam2fasta(samfile, **kwargs):
    """Extract the reads from an samfile and print as fasta.
        """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna

    chrom = kwargs.get('chrom', None)
    start = kwargs.get('start', 0)
    stop = kwargs.get('stop', 0)
    minlen = kwargs.get('minlen', 1)
    strand = kwargs.get('strand', 2)

    ins_per_read = {read.qname: getinsertions(read.cigar, read.pos) \
                    for read in samfile.fetch(chrom, start, stop)}
    max_ins_per_pos = {}
    for i in range(start, stop):
        mx = max([v[i] if i in v.keys() else 0 for k, v in ins_per_read.items()])
        if mx: max_ins_per_pos[i] = mx

    print ins_per_read, max_ins_per_pos
    records = []
    for read in samfile.fetch(chrom, start, stop):
        ref = read.tid
        print ref
        # Count soft_clipped bases (CIGAR 4) plus the insertions (CIGAR 1)
        soft_clipped = sum([count for op, count in read.cigar if op in (4, 1)])
        if read.rlen - start + read.pos + 1 > minlen + soft_clipped and \
                stop - read.pos + 1 >= minlen + soft_clipped and \
                (strand == 2 or read.is_reverse == strand):
            fasta_seq = getseq(read, start=start, stop=stop, \
                               qcutoff=kwargs.get('qcutoff', QC), \
                               ins_position=ins_per_read[read.qname], \
                               max_ins_per_pos=max_ins_per_pos)
            rec = SeqRecord(Seq(fasta_seq, generic_dna), id=read.qname, \
                            description='%s:%i:%s' % \
                            (chrom, read.pos, STRAND[read.is_reverse]))
            records.append(rec)
        else:
            print read.rlen, read.pos, soft_clipped

    total = SeqIO.write(records, kwargs.get('out', sys.stdout), 'fasta')
    print >> sys.stderr, "[samtools] Fetched %i reads from %s:%i-%i." % \
                        (total, chrom, start, stop)


def pairwisealign(reads_file, ref_file):
    """Wrapper for an external aligner, then calls pysam/samtools
       to convert/index/sort SAM file into sorted BAM
    """
    import os
    import subprocess
    import pysam

    stem = '.'.join(os.path.split(reads_file)[1].split('.')[:-1])
    sam_file = stem + '.sam'

    if 'smalt' in aligner:
        retcode = subprocess.call(\
                '%s index -k 8 -s 2 ref_index %s &> aligner.log' \
                % (aligner, ref_file), shell=True)
        if retcode < 0:
            sys.exit('smalt index terminated by signal', retcode)
        retcode = subprocess.call(\
                '%s map -x -f samsoft -y 0.8 -o %s ref_index %s >> \
                aligner.log 2>&1' \
                % (aligner, sam_file, reads_file), shell=True)
        if retcode < 0:
            sys.exit('smalt map terminated by signal', retcode)

    if 'needle' in aligner:
        sys.exit('needle not yet implemented')

    pysam.faidx(ref_file)
    sorted_bam = stem + '_sorted'
    cml = 'samtools view -u -b -t %s %s | samtools sort - %s' % \
            (ref_file, sam_file, sorted_bam)
    subprocess.call(cml, shell=True)
    pysam.index(sorted_bam + '.bam')

    return sorted_bam + '.bam'


def main(reads_file, ref_file, far_file):
    """Reads can be either in fasta or fastq format
    """
    import pysam
    bam_file = pairwisealign(reads_file, ref_file)
    samfile = pysam.Samfile(bam_file, 'rb')
    if far_file == None:
        output = sys.stdout
    else:
        output = open(far_file, 'w')
    bam2fasta(samfile, chrom=None, start=200, stop=300, \
              minlen=1, out=output, qcutoff=QC, strand=2)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], 'fake.far')
