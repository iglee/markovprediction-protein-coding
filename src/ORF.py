import pandas as pd
import numpy as np
import re
from collections import Counter

STOP_CODON = ["TAA","TAG","TGA"]

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

    return zip(stop_locations[:-1], stop_locations[1:])

# markov model
class MarkovModel:
    def __init__(self, k, seq):
        self.seq = seq
        self.k = k
        self.orfs = find_orf(seq)
        self.trusted_orfs = [x for x in self.orfs if len(x) >=1400]

        # counts
        self.kmer_counts = None
        self.kponemer_counts = None
        self.probs = None

        # initialize with given input
        self.kmer_counts = self.count_kmers(self.k)
        self.kponemer_counts = self.count_kmers(self.k+1)
    
    

    def count_kmers(self, k):
        kmers = [self.seq[i-k:i] for i in range(k,len(self.seq)+1)]
        return Counter(kmers)


    def calculate_probs(self):
        return None
