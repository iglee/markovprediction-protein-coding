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
        self.pseudocount = 1

        # orfs
        self.orfs = ORF(seq)

        # long orf counts
        self.kmer_counts = self.count_kmers(self.k, self.orfs.long_orfs)
        self.kponemer_counts = self.count_kmers(self.k+1, self.orfs.long_orfs)
        self.V = len(self.kmer_counts.keys()) # normalization for add-1 smoothing, use total vocabulary count for unique kmers
        self.start_counts = self.count_starts()

        # background sequences counts
        self.bg_kmer_counts = self.count_kmers(self.k, self.orfs.background)
        self.bg_kponemer_counts = self.count_kmers(self.k+1, self.orfs.background)

        self.start_probs = self.calculate_start_proba(self.start_counts)


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
            "P count(AAGxyT): " + "\n" + self.print_count(self.kponemer_counts) + "\n" + \
            "Q count(AAGxyT): " + "\n" + self.print_count(self.bg_kponemer_counts) 


    def calculate_start_proba(self, start_counts):
        norm = sum(start_counts.values())
        temp = start_counts.copy()
        for key in temp:
            temp[key] /= norm
        return temp

    def conditional_proba(self, token, kponemer_counts, kmer_counts):
        if token[:-1] not in kmer_counts and token not in kponemer_counts:
            return log(self.pseudocount / self.V)
        elif token[:-1] in kmer_counts and token not in kponemer_counts:
            return log(self.pseudocount / (kmer_counts[token[:-1]] + self.V))
        return log( (kponemer_counts[token] + self.pseudocount) / (kmer_counts[token[:-1]] + self.V) )

    def sequence_proba(self, seq, kponemer_counts, kmer_counts):
        logprob = self.start_probs[seq[0:self.k]] # initialize with start proba
        for i in range(self.k+2,len(seq)+1):
            logprob += self.conditional_proba(seq[i-self.k-1:i],kponemer_counts, kmer_counts)
        return logprob
        
    def score(self, seq):
        P_score = self.sequence_proba(seq, self.kponemer_counts, self.kmer_counts)
        Q_score = self.sequence_proba(seq, self.bg_kponemer_counts, self.bg_kmer_counts)
        return P_score - Q_score


    def results(self):
        results = []
        for seq, loc in zip(self.orfs.total_orfs, self.orfs.all_orf_locations):
            results.append({ "start" : loc[0], "end" : loc[1], "length" : len(seq), "score" : self.score(seq)})

        return results

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