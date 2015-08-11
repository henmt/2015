'''

A simplified copy of SACGF libraries for Shuly Lim's HEN-1 paper

Created on 08/11/2013

@author: dlawrence
'''

from pyfasta import Fasta #@UnresolvedImport

class Reference(object):
    def __init__(self, genome_fasta):
        # @see: https://pypi.python.org/pypi/pyfasta
        key_fn = lambda key : key.split()[0] # Use first value before whitespace as keys
        self.fasta =  Fasta(genome_fasta, key_fn=key_fn)

    def get_sequence_from_iv(self, iv):
        feature_hash = {'chr' : iv.chrom, 'start' : iv.start, 'stop' : iv.end, 'strand' : iv.strand}
        return self.fasta.sequence(feature_hash, one_based=False)

