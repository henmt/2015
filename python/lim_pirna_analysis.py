#!/usr/bin/env python
'''
ACRF Cancer Genomics Facility

This file performs the analyis based on output from lim_pirna_process_bam
run on short RNA bams

Created on 12/11/2013

@author: dlawrence
'''

from lim_pirna_common import SECONDARY_PIRNA, PRIMARY_PIRNA, PIRNA_TYPES, \
    SAMPLE_NAMES, PIRNA_MIN_SIZE, PIRNA_MAX_SIZE
from lim_pirna_process_bam import PIRNA_SEQUENCE_U1_OR_A10_LENGTH
from lim_pirna_sequence_length_analysis import plot_sequence_lengths
from scipy import stats
import lim_pirna_common
import os
import pandas as pd
import sys

TOTAL_READS = 'total reads'
S_PIRNA_RATIO = "s-piRNA ratio"
TOTAL_PIRNA_PERCENT = "total piRNA %"

def load_dataframe_from_files(classification_dir, file_name):
    classification_series = {}
    for subdir in os.listdir(classification_dir):
        fullpath = os.path.join(classification_dir, subdir)
#        print "fullpath = %s" % fullpath
        if not os.path.isdir(fullpath):
            raise ValueError("Expected '%s' to be a directory" % fullpath) 
        series = pd.Series.from_csv(os.path.join(fullpath, file_name))
#        print "series = %s" % series
        classification_series[subdir] = series
            
    return pd.DataFrame.from_dict(classification_series, orient='index')

def load_classification_counts(alignment_type, classification_method, process_bam_output_dir):
    samples = {}
    for subdir in os.listdir(process_bam_output_dir):
        fullpath = os.path.join(process_bam_output_dir, subdir, alignment_type)
        if os.path.isdir(fullpath):
            classification_dir = os.path.join(fullpath, classification_method)
            counts_df = load_dataframe_from_files(classification_dir, 'CounterProcessor.csv')
            name = subdir.split('.')[0]
            samples[name] = counts_df['count']

    return pd.DataFrame.from_dict(samples, orient='index')

def compare_column_across_conditions(column, df):
    df = df.copy() # Don't modify what was passed in
    lim_pirna_common.add_conditions_and_sample_types(df)

    comparison_df = pd.DataFrame(index=SAMPLE_NAMES.values(), columns=['WT mean', 'pin/pin mean', 'p-value'])    

    for (sample_type, sample_df) in df.groupby("sample_type"):
        condition_gb = sample_df.groupby("condition")[column]
        wt = condition_gb.get_group('WT')
        mut = condition_gb.get_group('pin/pin')
    
        (_, p) = stats.ttest_ind(wt, mut)
        sample_stats = comparison_df.ix[sample_type]
        sample_stats['WT mean'] = wt.mean()
        sample_stats['pin/pin mean'] = mut.mean()
        sample_stats['p-value'] = p

    comparison_df.to_csv("%s_comparison.csv" % column.replace(' ', '_'))

def write_classification_stats(alignment_type, df):
    save_df_with_condition_and_samples('%s_pirna_classification_counts.csv' % alignment_type, df)

    percent_column = lambda col : "%s %%" % col
    columns = [percent_column(i) for i in df]
    summary_df = pd.DataFrame(index=df.index, columns=columns)
    total_reads = df.sum(axis=1)
    for col in df:
        summary_df[percent_column(col)] = 100.0 * df[col] / total_reads
    
    summary_df[S_PIRNA_RATIO] = df[SECONDARY_PIRNA] * 1.0 / df[PRIMARY_PIRNA]
    summary_df[TOTAL_PIRNA_PERCENT] = 100.0 * df[PIRNA_TYPES].sum(axis=1) / total_reads
    
    save_df_with_condition_and_samples("pirna_sequence_count_summary.csv", summary_df)

    compare_column_across_conditions(TOTAL_PIRNA_PERCENT, summary_df)
    compare_column_across_conditions(S_PIRNA_RATIO, summary_df)


def load_combined_pirna_data_from_series_files_in_dirs(alignment_type, classification_method, process_bam_output_dir, file_name):
    samples = {}
    for subdir in os.listdir(process_bam_output_dir):
        fullpath = os.path.join(process_bam_output_dir, subdir, alignment_type)
        if os.path.isdir(fullpath):
            classification_dir = os.path.join(fullpath, classification_method)
            read_lengths_per_classification_df = load_dataframe_from_files(classification_dir, file_name)
            name = subdir.split('.')[0]
            samples[name] = read_lengths_per_classification_df.ix[PIRNA_TYPES].sum() # Combine 2 types

    return pd.DataFrame.from_dict(samples, orient='index')

def save_df_with_condition_and_samples(filename, df):
    df = df.copy() # Don't modify what was passed in
    lim_pirna_common.add_conditions_and_sample_types(df)
    df.to_csv(filename)

if __name__ == '__main__':
    if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]):
        print >>sys.stderr, "Usage: %s lim_pirna_process_bam_output_dir" % os.path.basename(sys.argv[0])
        exit(1)
    
    process_bam_output_dir = sys.argv[1]
    
    classification_method = PIRNA_SEQUENCE_U1_OR_A10_LENGTH
    
    aligned_classification_counts_df = load_classification_counts('aligned', classification_method, process_bam_output_dir)
    write_classification_stats('aligned', aligned_classification_counts_df)

    aligned_pirna_read_lengths_df = load_combined_pirna_data_from_series_files_in_dirs('aligned', classification_method, process_bam_output_dir, 'ReadLengthProcessor.csv')
    aligned_pirna_read_lengths_df.to_csv('aligned_pirna_read_lengths.csv')
    plot_sequence_lengths('piRNAs', aligned_pirna_read_lengths_df, PIRNA_MIN_SIZE, PIRNA_MAX_SIZE)

    for end_modification in ['Uridylation', 'Adenylation']:
        file_name = '%sProcessor.csv' % end_modification
        aligned_pirna_end_mod_df = load_combined_pirna_data_from_series_files_in_dirs('aligned', classification_method, process_bam_output_dir, file_name)
        save_df_with_condition_and_samples('aligned_pirna_%s.csv' % end_modification, aligned_pirna_end_mod_df)

        total_potential_modified = aligned_pirna_end_mod_df['modified'] + aligned_pirna_end_mod_df['not modified']
        aligned_pirna_end_mod_df[end_modification] = 100.0 * aligned_pirna_end_mod_df['modified'] / total_potential_modified
        compare_column_across_conditions(end_modification, aligned_pirna_end_mod_df)

