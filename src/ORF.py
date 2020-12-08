import pandas as pd
import numpy as np
import re
from collections import Counter

STOP_CODON = ["TAA","TAG","TGA"]
COMPLEMENTS = {"A":"T","T":"A","C":"G","G":"C"}

# organize input data into a class
class GenomeData:
    def __init__(self):
        self.seq_name = None
        self.sequence = ""
        self.seq_len = 0

def read_fna(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    input_data = []

    for l in lines:
        if l[0] == ">":
            try:
                g.seq_len = len(g.sequence)
                input_data.append(g)
            except:
                pass
            g = GenomeData()
            g.seq_name = l.strip()
        else:
            g.sequence = g.sequence + re.sub(r"((?:(?!A|C|T|G)\S))","T",l.strip().upper())

    g.seq_len = len(g.sequence)
    input_data.append(g)

    return input_data

# find stop codons
def find_stops(seq):
    stop_locations = []

    for i in range(3,len(seq)+1,3):
        window = seq[i-3:i]
        
        if window in STOP_CODON:
            stop_locations.append(i)
    return np.array(stop_locations)

# calculate orf locations
def orf_locations(seq, stop_locations, start_reading):
    start_idxs = stop_locations + start_reading
    stop_idxs = start_idxs - 3
    start_idxs = np.insert(start_idxs, 0, [start_reading])
    stop_idxs = np.append(stop_idxs, [len(seq)])
    return list(zip(start_idxs, stop_idxs))

# orf sequences
def orf_seqs(seq, start_reading):
    stop_locations = find_stops(seq[start_reading:])
    orf_idxs = orf_locations(seq, stop_locations, start_reading)
    return [(i,j) for i, j in orf_idxs if j-i > 0], [seq[i:j] for i, j in orf_idxs if j-i > 0]



def background_seqs(trusted_orfs):

    background = []
    for x in trusted_orfs:
        reversed_orf = x[::-1]
        reverse_complements = ""
        for y in reversed_orf:
            reverse_complements += COMPLEMENTS[y]

        background.append(reverse_complements)
        
    return background

class ORF:
    def __init__(self, seq):
        self.seq = seq
        self.long_len = 1400
        self.short_len = 50


        # for reading frame starting at index = 0 (or index = 1 in biology)
        self.idxs0, self.orf0 = orf_seqs(seq, 0)
        # for reading frame starting at index = 1 
        self.idxs1, self.orf1 = orf_seqs(seq, 1)
        # for reading frame starting at index = 2
        self.idxs2, self.orf2 = orf_seqs(seq, 2)
        # total orfs
        self.total_orfs = self.orf0 + self.orf1 + self.orf2
        self.all_orf_locations = self.idxs0 + self.idxs1 + self.idxs2
        # long orfs
        self.long_orfs = [x for x,y in zip(self.total_orfs,self.all_orf_locations) if len(x)> self.long_len]
        self.long_orfs_location = [y for x,y in zip(self.total_orfs,self.all_orf_locations) if len(x)> self.long_len]
        # short orfs
        self.short_orfs = [x for x,y in zip(self.total_orfs,self.all_orf_locations) if len(x)< self.short_len]
        self.short_orfs_location = [y for x,y in zip(self.total_orfs,self.all_orf_locations) if len(x)< self.short_len]
        # calculate background orfs from long_orfs
        self.background = background_seqs(self.long_orfs)

    def __repr__(self):
        return "total number of orfs found: {} \
            \nnumber of long orfs: {} \
            \nnumber of short orfs: {} ".format(len(self.total_orfs), len(self.long_orfs), len(self.short_orfs))
