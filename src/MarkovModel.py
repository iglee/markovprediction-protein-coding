import pandas as pd
import numpy as np
import re
from collections import Counter
from .ORF import ORF, read_fna


class MarkovModel:
    def __init__(self, k, seq):
        self.seq = seq
        self.k = k


        self.orfs = ORF(seq)

        # counts
        self.kmer_counts = self.count_kmers(self.k)
        self.kponemer_counts = self.count_kmers(self.k+1)
        self.start_counts = self.count_starts()

        self.probs = None


    def generate_kmers(self, seq, k):
        kmers = [seq[i-k:i] for i in range(k,len(seq)+1)]
        return kmers

    def count_kmers(self, k):
        total_kmers = []
        for x in self.orfs.long_orfs:
            total_kmers += self.generate_kmers(x,k)

        return Counter(total_kmers)

    def count_starts(self):
        starts = []
        for x in self.orfs.long_orfs:
            starts.append(x[0:self.k])
        return Counter(starts)

    def calculate_probs(self):
        return None


def main():
    data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
    seq = data[0].sequence
    #orf = ORF(seq)
    mm = MarkovModel(5, seq)
    #print(mm.stop_idxs[:5])
    #print(mm.orfs[:5])
    

if __name__ == "__main__":
    main()