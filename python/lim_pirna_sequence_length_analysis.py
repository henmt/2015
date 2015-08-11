#!/usr/bin/env python
'''
Created on 08/11/2013

@author: dlawrence
'''

from lim_pirna_common import add_conditions_and_sample_types, SAMPLE_TYPE, \
    CONDITION, GROUP_BY_COLUMNS
from matplotlib import pyplot
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from sacgf_library import sacgf_utils
from sacgf_library import fastqc_parser
import numpy as np
import os
import pandas as pd
import sys

# Cell guidelines
JOURNAL_FONT = 'Helvetica'
JOURNAL_GRAPH_DPIS = [200, 1000] # Do a big one and a small one

colour_sets = {
    "pink_blue" : {"WT" : "#f6756b", "pin/pin" : "#01bec3"}
#    "greyscale" : {"WT" : "darkgrey", "pin/pin" : "black"},
#    "greyscale_highcontrast" : {"WT" : "silver", "pin/pin" : "black"},
#    "black_white" : {"WT" : "white", "pin/pin" : "black"},
#    "green_red" : {"WT" : "green", "pin/pin" : "red"},
#    "blue_red" : {"WT" : "blue", "pin/pin" : "red"},
}

MIN_LENGTH = 18
MAX_LENGTH = 35

def load_sequence_lengths_df(fastqc_dir):
    fastqc_data = fastqc_parser.load_fastqc_data_for_dir(fastqc_dir)
    sample_sequence_length_arrays = {k : v[fastqc_parser.SEQUENCE_LENGTH_DISTRIBUTION] for (k,v) in fastqc_data.iteritems()}
    return pd.DataFrame.from_dict(sample_sequence_length_arrays, 'index')

def plot_greyscale_boxplot(graph_image, description, sample_type, sample_df, colors, y_limit, min_length, max_length, plot_mean):
    pyplot.rc('font',**{'family':'sans-serif','sans-serif':[JOURNAL_FONT]})

    for dpi in JOURNAL_GRAPH_DPIS:
        figure = Figure(dpi=dpi)
        figure.patch.set_facecolor('white')
        ax = figure.add_subplot(1, 1, 1)
    
        for (condition, df) in sample_df.groupby(CONDITION):
            numeric_columns = df.columns - GROUP_BY_COLUMNS
            numeric_df = df[numeric_columns]
    
            # Remove 0 length from boxplot - it will use axis 1... for elements 0... 
            boxplot_df = numeric_df[numeric_df.columns[1:]]
    
            data_lists = []
            for column in boxplot_df:
                data_lists.append(boxplot_df[column])
    
            if plot_mean:
                lines = ax.plot(numeric_df.mean(axis=0), linestyle='dashed', zorder=-1)
                pyplot.setp(lines, color=colors[condition])
                
            r = ax.boxplot(data_lists, patch_artist=True)
            pyplot.setp(r.values(), color=colors[condition], lw=1)
            pyplot.setp(r['boxes'], facecolor='white')
    
            patches = []
            labels = []
            for name in ["WT", "pin/pin"]:
                color = colors[name]
                labels.append(name)
                patches.append(Rectangle((0, 0), 1, 1, fc=color))
            
            ax.legend(patches, labels, loc=2)
    
        # Only show range we're interested in
        # Use .5 to get a bit of spacing, & don't show numbers either side
        ax.set_xlim(min_length - .5, max_length + .5) 
        ax.set_ylim(0, y_limit)
        ax.set_xlabel("Read length", weight='bold', size=14)
        ax.set_ylabel("Percentage of %s" % description, weight='bold', size=14)
        ax.set_title(sample_type, weight='bold', size=18)
    
        canvas = FigureCanvasAgg(figure)
        canvas.print_tif(graph_image + '_%s_dpi.tiff' % dpi)
        canvas.print_png(graph_image + '_%d_dpi.png' % dpi)
        #canvas.print_eps(graph_image + '.eps')
    #    canvas.print_pdf(graph_image + '.pdf')


def plot_sequence_lengths(description, sequence_lengths_df, min_length, max_length):
    graph_dir = "%s_read_length_graphs" % description
    sacgf_utils.mk_path(graph_dir)

    total_reads = sequence_lengths_df.sum(axis=1)
    reads_percent = sequence_lengths_df.divide(total_reads / 100.0, axis=0)

    max_percent = reads_percent.max().max()
    y_limit = np.floor(1 + max_percent / 5) * 5 # use multiples of 5
    add_conditions_and_sample_types(reads_percent)

    for (sample_type, sample_df) in reads_percent.groupby(SAMPLE_TYPE):
        for (color_scheme, colors) in colour_sets.iteritems():
            for plot_mean in [False, True]:
                mean_description = "_mean" if plot_mean else ""
                graph_image = "read_lengths_%s%s_%s_boxplot" % (color_scheme, mean_description, sample_type)
                graph_image = os.path.join(graph_dir, graph_image)
                plot_greyscale_boxplot(graph_image, description, sample_type, sample_df, colors, y_limit, min_length, max_length, plot_mean)

def analyse_unaligned_sequence_lengths(fastqc_dir):
    ''' Works from FastQC data '''
    sequence_lengths_df = load_sequence_lengths_df(fastqc_dir)
    plot_sequence_lengths("total reads", sequence_lengths_df, MIN_LENGTH, MAX_LENGTH)

if __name__ == '__main__':
    if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]):
        print >>sys.stderr, "Usage: %s fastqc_dir" % os.path.basename(sys.argv[0])
        exit(1)
    
    analyse_unaligned_sequence_lengths(sys.argv[1])
