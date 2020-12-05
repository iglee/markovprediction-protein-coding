import pandas as pd
import numpy as np
import re
from collections import Counter
from ORF import find_orf, background_seqs, read_fna


class MarkovModel:
    def __init__(self, k, seq):
        self.seq = seq
        self.k = k

        self.stop_idxs, self.orfs = find_orf(self.seq)
        self.trusted_orfs = [x for x in self.orfs if len(x) >=1400]
        self.backgrounds = background_seqs(self.trusted_orfs)

        # counts
        self.kmer_counts = self.count_kmers(self.k)
        self.kponemer_counts = self.count_kmers(self.k+1)

        self.probs = None


    def generate_kmers(self, seq, k):
        kmers = [seq[i-k:i] for i in range(k,len(seq)+1)]
        return kmers

    def count_kmers(self, k):
        total_kmers = []
        for x in self.trusted_orfs:
            total_kmers += self.generate_kmers(x,k)

        return Counter(total_kmers)


    def calculate_probs(self):
        return None


def main():
    data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
    seq = data[0].sequence
    mm = MarkovModel(5, seq)
    #print(mm.stop_idxs[:5])
    #print(mm.orfs[:5])
    

if __name__ == "__main__":
    main()