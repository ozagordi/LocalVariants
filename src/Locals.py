aligner = '/Users/ozagordi/Projects/3rd_party_SW/smalt-0.6.1/smalt_MacOSX_i686_64bit'

gene_coord = {
    # gene(lowercase): [start, stop]
    'pol': [2085, 5096],
    'protease': [2253, 2549],
    'env': [6225, 8795],
    'rt': [2550, 3869]
}

dna_code = {'A': set(['A']),
    'C': set(['C']),
    'G': set(['G']),
    'T': set(['T']),
    
    'R': set(['G', 'A']),
    'Y': set(['T', 'C']),
    'M': set(['A', 'C']),
    'K': set(['G', 'T']),
    'S': set(['G', 'C']),
    'W': set(['A', 'T']),
    
    'H': set(['A', 'C', 'T']),
    'B': set(['C', 'G', 'T']),
    'V': set(['A', 'C', 'G']),
    'D': set(['A', 'G', 'T']),
    'N': set(['A', 'C', 'G', 'T']),
    '-': set(['A', 'C', 'G', 'T'])
}