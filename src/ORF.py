import pandas as pd
import numpy as np
import re

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