#!/usr/bin/env python
'''
Created on 08/11/2013

@author: dlawrence
'''
from argparse import ArgumentParser
from collections import defaultdict
from sacgf_library.reference import Reference
from sacgf_library.sacgf_utils import name_from_file_name, GrowingList, \
    mk_path_for_file
import HTSeq
import abc
import lim_pirna_common
import os
import pandas as pd

PIRNA_SEQUENCE_U1_OR_A10 = 'pirna_sequence_u1_or_a10'
PIRNA_SEQUENCE_U1_XOR_A10 = 'pirna_sequence_u1_xor_a10'
PIRNA_SEQUENCE_U1_OR_A10_LENGTH = 'pirna_sequence_u1_or_a10_length'
PIRNA_SEQUENCE_U1_XOR_A10_LENGTH = 'pirna_sequence_u1_xor_a10_length'

PIRNA_CLASSIFIERS = [PIRNA_SEQUENCE_U1_OR_A10,
                     PIRNA_SEQUENCE_U1_XOR_A10,
                     PIRNA_SEQUENCE_U1_OR_A10_LENGTH,
                     PIRNA_SEQUENCE_U1_XOR_A10_LENGTH,
]

class Processor(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, alignment_status_name, pirna_classifier, pirna_type):
        self.name = name
        self.alignment_status_name = alignment_status_name
        self.pirna_classifier = pirna_classifier
        self.pirna_type = pirna_type

    @abc.abstractmethod
    def process(self, *args):
        pass

    @abc.abstractmethod
    def get_data(self):
        ''' DataFrame/Series or an object with a .to_csv() method '''
        pass
    
    def write_output(self):
        df = self.get_data()
        filename = os.path.join(self.name, self.alignment_status_name, self.pirna_classifier, self.pirna_type, self.__class__.__name__) + '.csv'
        mk_path_for_file(filename)
        df.to_csv(filename)

class CounterProcessor(Processor):
    def __init__(self, *args):
        super(CounterProcessor, self).__init__(*args)
        self.count = 0
        
    def process(self, *args):
        self.count += 1

    def get_data(self):
        return pd.Series([self.count], index=['count'])


class ReadLengthProcessor(Processor):
    def __init__(self, *args):
        super(ReadLengthProcessor, self).__init__(*args)
        self.read_lengths = GrowingList(initial_value=0)

    def process(self, *args):
        aln = args[0]
        self.read_lengths[len(aln.read.seq)] += 1

    def get_data(self):
        return pd.Series([i for i in self.read_lengths])

class UnalignedEndHomopolymerProcessor(Processor):
    NUM_HOMOPOLYMERS = 3
    base = None

    def __init__(self, *args):
        super(UnalignedEndHomopolymerProcessor, self).__init__(*args)
        homopolymers = [self.base*i for i in range(1,self.NUM_HOMOPOLYMERS+1)]
        self.threep_homopolymer_counter = pd.Series(index=homopolymers).fillna(0)

    def process(self, *args):
        aln = args[0]
        for c in self.threep_homopolymer_counter.index:
            if aln.read.seq[-len(c):] == c:
                self.threep_homopolymer_counter[c] += 1

    def get_data(self):
        return self.threep_homopolymer_counter


class UnalignedAdenylationProcessor(UnalignedEndHomopolymerProcessor):
    base = 'A'

class UnalignedUridylationProcessor(UnalignedEndHomopolymerProcessor):
    base = 'T'

class AlignedEndBaseModifiedProcessor(Processor):
    base = None

    def __init__(self, *args):
        super(AlignedEndBaseModifiedProcessor, self).__init__(*args)
        self.modification_counter = pd.Series(index=['unknown', 'modified', 'not modified']).fillna(0)

    def process(self, *args):
        aln = args[0]
        reference_sequence = args[1]
        if reference_sequence[-1] == self.base:
            self.modification_counter['unknown'] += 1
        elif aln.read.seq[-1] == self.base:
            self.modification_counter['modified'] += 1
        else:
            self.modification_counter['not modified'] += 1

    def get_data(self):
        return self.modification_counter


class AdenylationProcessor(AlignedEndBaseModifiedProcessor):
    base = 'A'

class UridylationProcessor(AlignedEndBaseModifiedProcessor):
    base = 'T'

PROCESSORS = {
    "counter" : CounterProcessor,
    "read_length" : ReadLengthProcessor,
    "adenylation" : AdenylationProcessor,
    "uridylation" : UridylationProcessor,
    "unaligned_uridylation" : UnalignedUridylationProcessor,
    "unaligned_adenylation" : UnalignedAdenylationProcessor,
}

def displatch_processors(processors, *args):
    for p in processors:
        p.process(*args)

def classify_pirna_then_dispatch_processors(classification_processors, *args):
    seq = args[0].read.seq

    t_at_1 = seq[0] == 'T'
    a_at_10 = seq[9] == 'A'
    pirna_sized = lim_pirna_common.PIRNA_MIN_SIZE <= len(seq) <= lim_pirna_common.PIRNA_MAX_SIZE

    if t_at_1:
        pirna_classification = lim_pirna_common.PRIMARY_PIRNA
    elif a_at_10:
        pirna_classification = lim_pirna_common.SECONDARY_PIRNA
    else:
        pirna_classification = lim_pirna_common.OTHER

    displatch_processors(classification_processors[PIRNA_SEQUENCE_U1_OR_A10][pirna_classification], *args)
    if pirna_sized:
        displatch_processors(classification_processors[PIRNA_SEQUENCE_U1_OR_A10_LENGTH][pirna_classification], *args)
        
    if t_at_1 ^ a_at_10:
        if t_at_1:
            pirna_classification = lim_pirna_common.PRIMARY_PIRNA
        else:
            pirna_classification = lim_pirna_common.SECONDARY_PIRNA
    else:
        pirna_classification = lim_pirna_common.OTHER

    displatch_processors(classification_processors[PIRNA_SEQUENCE_U1_XOR_A10][pirna_classification], *args)
    if pirna_sized:
        displatch_processors(classification_processors[PIRNA_SEQUENCE_U1_XOR_A10_LENGTH][pirna_classification], *args)

def setup_processors(name, alignment_status_name, processor_names):
    processors = defaultdict(lambda : defaultdict(list))
    processors_list = []
    
    for pirna_classifier in PIRNA_CLASSIFIERS:
        for pirna_classification in lim_pirna_common.CLASSIFICATION_TYPES:
            for p in processor_names:
                clazz = PROCESSORS[p]
                processor = clazz(name, alignment_status_name, pirna_classifier, pirna_classification)
                processors[pirna_classifier][pirna_classification].append(processor)
                processors_list.append(processor)

    return (processors, processors_list)


def process_alignments(name, reference, alignment_iterator):
    unaligned_processor_names = ["counter",
                                 "read_length",
                                 "unaligned_uridylation",
                                 "unaligned_adenylation",
    ]
    aligned_processor_names = unaligned_processor_names + ["adenylation", "uridylation"]

    (unaligned_processors, all_unaligned_processors) = setup_processors(name, "unaligned", unaligned_processor_names)
    (aligned_processors, all_aligned_processors) = setup_processors(name, "aligned", aligned_processor_names)

    for aln in alignment_iterator:
        if aln.aligned:
            reference_sequence = reference.get_sequence_from_iv(aln.iv)
            classify_pirna_then_dispatch_processors(aligned_processors, aln, reference_sequence)
        else:
            classify_pirna_then_dispatch_processors(unaligned_processors, aln)
            
    return all_unaligned_processors + all_aligned_processors

    
def handle_args():
    parser = ArgumentParser()
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("bam")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = handle_args()
    name = name_from_file_name(args.bam)

    reference = Reference(args.genome_fasta)
    processors = process_alignments(name, reference, HTSeq.BAM_Reader(args.bam))
    for p in processors:
        p.write_output()
