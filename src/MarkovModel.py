import pandas as pd
import numpy as np
import re
from collections import Counter
from .ORF import ORF, read_fna
from math import log

nucleotides = list("ACGT")

class MarkovModel:
    def __init__(self, k, seq):
        self.seq = seq
        self.k = k
        self.pseudocount = 0.1

        # orfs
        self.orfs = ORF(seq)

        # long orf counts
        self.kmer_counts = self.count_kmers(self.k, self.orfs.long_orfs)
        self.kponemer_counts = self.count_kmers(self.k+1, self.orfs.long_orfs)
        self.start_counts = self.count_starts()

        # background sequences counts
        self.bg_kmer_counts = self.count_kmers(self.k, self.orfs.background)
        self.bg_kponemer_counts = self.count_kmers(self.k+1, self.orfs.background)

        self.probs = None


    def generate_kmers(self, seq, k):
        kmers = [seq[i-k:i] for i in range(k,len(seq)+1)]
        return kmers

    def count_kmers(self, k, seq):
        total_kmers = []
        for x in seq:
            total_kmers += self.generate_kmers(x,k)

        return Counter(total_kmers)

    def count_starts(self):
        starts = []
        for x in self.orfs.long_orfs:
            starts.append(x[0:self.k])
        return Counter(starts)

    def print_count(self, counts):
        
        df = pd.DataFrame(columns=nucleotides, index=nucleotides)

        for x in nucleotides:
            for y in nucleotides:
                key = "AAG"+x+y+"T"
                df[y][x] = counts[key]
        
        return str(df)

    def __repr__(self):
        return "ORFs Found: " + str(self.orfs) + "\n\n" + \
            "P: count(AAGxyT): " + "\n" + self.print_count(self.kponemer_counts) + "\n" + \
            "Q: count(AAGxyT): " + "\n" + self.print_count(self.bg_kponemer_counts)


    def conditional_proba(self, token, kponemer_counts, kmer_counts):
        return log(kponemer_counts[token] / kmer_counts[token[:-1]])

    def sequence_proba(self, seq, kponemer_counts, kmer_counts):
        logprob = 0 # initialize with start proba
        for i in range(self.k,len(seq)+1):
            logprob += self.conditional_proba(seq[i-self.k:i],kponemer_counts, kmer_counts)
        return logprob
        


def main():
    # test set
    goldens = pd.read_csv("data/plusgenes-subset.gff", delimiter="\t", header=None)

    # input data
    data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
    seq = data[0].sequence
    
    # markov model
    mm = MarkovModel(5, seq)

    print(mm)
    

if __name__ == "__main__":
    main()