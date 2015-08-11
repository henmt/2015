'''
Created on 08/11/2013

@author: dlawrence
'''
import csv
import os
import sys

def mk_path(path):
    if path and not os.path.exists(path):
        os.makedirs(path)

def mk_path_for_file(f):
    mk_path(os.path.dirname(f))

def name_from_file_name(file_name):
    '''Gets file name without extension or directory'''
    return os.path.splitext(os.path.basename(file_name))[0]

def write_csv_dict(csv_file, headers, rows, extrasaction=None, dialect=None):
    '''
    default dialect = 'excel', other dialect option: 'excel-tab'
    These are the same optional arguments as csv.DictWriter
    headers=keys for dicts
    rows=list of dicts
    '''
    if extrasaction is None:
        extrasaction = "raise"
    if dialect is None:
        dialect = 'excel'
    
    with open(csv_file, 'wb') as f:
        writer = csv.DictWriter(f, headers, extrasaction=extrasaction, dialect=dialect)
    #        writer.writeheader()
    #        ersa uses python 2.6
        writer.writerow(dict(zip(headers, headers)))
        writer.writerows(rows)

def get_non_zero_start_end_index(array):
    start = sys.maxint
    end = 0
    for (i, c) in enumerate(array):
        if len(c):
            start = min(start, i)
            end = max(end, i)
    return (start, end)

class GrowingList(list):
    def __init__(self, *args, **kwargs):
        super(GrowingList, self).__init__(*args)
        self.initial_value = kwargs.get("initial_value")

    def __getitem__(self, index):
        if index >= len(self):
            self.extend([self.initial_value]*(index + 1 - len(self)))
        return list.__getitem__(self, index)

    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([self.initial_value]*(index + 1 - len(self)))
        list.__setitem__(self, index, value)