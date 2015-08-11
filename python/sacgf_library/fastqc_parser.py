'''
Created on 08/06/2013

@author: dlawrence
'''

import glob
import numpy as np
import os
import re

def parse_sequence_length(lines):
    last_value = int(lines[-1].split()[0])
    array = np.zeros(last_value+1)

    for line in lines:
        if line.startswith("#"):
            continue
        (length, count) = line.split()
        array[int(length)] = count
    return array

def parse_basic_statistics(lines):
    basic_stats = {}
    
    type_converters = {
        "Total Sequences" : int,
        "Filtered Sequences" : int,
        "%GC" : int,
    }
    
    for line in lines:
        if line.startswith("#"):
            continue
        (key, value) = line.split('\t')
        conversion_func = type_converters.get(key, str)
        basic_stats[key] = conversion_func(value)
    return basic_stats

#Section headers
BASIC_STATISTICS = "Basic Statistics"
SEQUENCE_LENGTH_DISTRIBUTION = "Sequence Length Distribution"

# Add lines here 
SECTION_HANDLERS = {BASIC_STATISTICS : parse_basic_statistics,
                    SEQUENCE_LENGTH_DISTRIBUTION : parse_sequence_length,}

def read_fastqc_data(fastqc_data_filename, custom_handlers=None):
    ''' custom_handlers : dict {"section header" : function passed lines, returns parsed data}
        Returns a dict of section header : parsed data '''

    handlers = SECTION_HANDLERS.copy()
    if custom_handlers:
        handlers.update(custom_handlers)
    
    section_pattern = re.compile(r">>(.*)\t")
    fastqc_data = {}
    
    with open(fastqc_data_filename) as f:
        lines = []
        section = None
        for line in f:
            line = line.rstrip()
            if line.startswith(">>"):
                if line[2:].startswith("END_MODULE"):
                    if section in handlers:
                        fastqc_data[section] = handlers[section](lines)
                else:
                    m = section_pattern.match(line)
                    if m:
                        token = m.group(1)
                        section = token
                        lines = []
            else:
                lines.append(line)

    return fastqc_data


def load_fastqc_data_for_dir(fastqc_dir):
    '''Returns a dict of files:fastqc data'''
    fastqc_data_for_dirs = {}
    
    for subdir in glob.glob1(fastqc_dir, "*_fastqc"):
        fullpath = os.path.join(fastqc_dir, subdir)
        if os.path.isdir(fullpath):
            fastqc_data_filename = os.path.join(fullpath, "fastqc_data.txt")
            fastqc_data_for_dirs[subdir] = read_fastqc_data(fastqc_data_filename)

    return fastqc_data_for_dirs

