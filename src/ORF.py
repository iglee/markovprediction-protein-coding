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

# find ORF
def find_orf(seq):
    stop_locations = []

    for i in range(3,len(seq)+1,3):
        window = seq[i-3:i]
        
        if window in STOP_CODON:
            stop_locations.append(i)

    stop_locations.insert(0,0)

    idxs= list(zip(stop_locations[:-1], stop_locations[1:]))

    return idxs, [seq[i:j-3] for i, j in idxs]

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
        # for reading frame starting at index = 0 (or index = 1 in biology)
        self.stop_idxs0, self.orf0 = find_orf(self.seq)
        # for reading frame starting at index = 1 
        self.stop_idxs1, self.orf1 = find_orf(self.seq[1:])
        # for reading frame starting at index = 2
        self.stop_idxs2, self.orf2 = find_orf(self.seq[2:])
        # total orfs
        self.total_orfs = None
        # long orfs
        self.long_orfs = None
        # short orfs
        self.short_orfs = None

    def __repr__(self):
        return "total number of orfs found: {} \
            \nnumber of long orfs: {} \
            \nnumber of short orfs: {} ".format(len(self.total_orfs), len(self.long_orfs), len(self.short_orfs))
