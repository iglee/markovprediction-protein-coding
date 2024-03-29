# implementation of kth order markov model
# specifically for gene sequences (i.e. list of nucleotides)
# probabilities calculated using conditional probabilities and MLE approximations

import pandas as pd
import numpy as np
import re
from collections import Counter
from ORF import ORF, read_fna
from math import log

nucleotides = list("ACGT")

class MarkovModel:
    def __init__(self, k, pseudocount, seq, long_len=1400, short_len=50):
        self.seq = seq
        self.k = k
        self.pseudocount = pseudocount

        # orfs: parse given seq into ORF data structure
        self.orfs = ORF(seq, long_len, short_len)

        # long orf counts - count kmers, k+1mer (similar to ngrams but with nucleotides, i.e. example of 3-mer: ATG) for MLE approximations of Markov probabilities
        self.kmer_counts = self.count_kmers(self.k, self.orfs.long_orfs)
        self.kponemer_counts = self.count_kmers(self.k+1, self.orfs.long_orfs)
        self.start_counts = self.count_starts(self.orfs.long_orfs)

        # background sequences counts - process counts similarly, but for "background" orfs. We define background orfs to be reverse complement of a given sequence.
        self.bg_kmer_counts = self.count_kmers(self.k, self.orfs.background)
        self.bg_kponemer_counts = self.count_kmers(self.k+1, self.orfs.background)
        self.bg_start_counts = self.count_starts(self.orfs.background)



    def generate_kmers(self, seq, k):
        """ generate list of kmers given a sequence """
        kmers = [seq[i-k:i] for i in range(k,len(seq)+1)]
        return kmers

    def count_kmers(self, k, seq):
        """ return a Counter object of kmer counts given a sequence """
        total_kmers = []
        for x in seq:
            total_kmers += self.generate_kmers(x,k)

        return Counter(total_kmers)

    def count_starts(self, seqs):
        """ count the starting kmers, i.e. the beginning of sequence """
        starts = []
        for x in seqs:
            starts.append(x[0:self.k])
        return Counter(starts)

    def print_count(self, counts):
        """ this prints a quick report summarizing nucleotide counts as a sanity check """
        df = pd.DataFrame(columns=nucleotides, index=nucleotides)

        for x in nucleotides:
            for y in nucleotides:
                key = "AAG"+x+y+"T"
                df[y][x] = counts[key]
        
        return str(df)

    def __repr__(self):
        """ repr function for reporting, using print_count method """
        return "ORFs Found: " + str(self.orfs) + "\n\n" + \
            "P count(AAGxyT): " + "\n" + self.print_count(self.kponemer_counts) + "\n" + \
            "Q count(AAGxyT): " + "\n" + self.print_count(self.bg_kponemer_counts) 

    # initialize probability, MLE approximation assumptions described in equation (1), page 2 of pdf
    def calculate_start_proba(self, start_counts, token, V):
        """ MLE approximation of starting kmer probabilities, i.e. probabilities from counts of starting kmers only """
        norm = sum(start_counts.values())
        if token not in start_counts:
            return log(self.pseudocount/(self.pseudocount*V))
        else:
            return log((start_counts[token] + self.pseudocount)/(norm + self.pseudocount*V))

    # MLE approximated conditional probability and markov model probability calculation, page 2 of pdf
    def conditional_proba(self, token, kponemer_counts, kmer_counts, V):
        """ MLE approximation of kmer probabilities, i.e. probabilities from counts of kmers """
        if token[:-1] not in kmer_counts and token not in kponemer_counts:
            return log(self.pseudocount / (self.pseudocount*V))
        elif token[:-1] in kmer_counts and token not in kponemer_counts:
            return log(self.pseudocount / (kmer_counts[token[:-1]] + self.pseudocount*V))
        return log( (kponemer_counts[token] + self.pseudocount) / (kmer_counts[token[:-1]] + self.pseudocount*V) )

    def sequence_proba(self, seq, start_counts, kponemer_counts, kmer_counts):
        """ 
        calculate the probability based on counts using MLE approximation
        --
        input: starting kmer counts -> Counter, k+1mer counts -> Counter, kmer counts -> Counter
        output: log probability of a given sequence -> float
        """
        
        V = len(kmer_counts.keys())
        logprob = self.calculate_start_proba(start_counts , seq[0:self.k], V)
        
        for i in range(self.k+2,len(seq)+1):
            logprob += self.conditional_proba(seq[i-self.k-1:i],kponemer_counts, kmer_counts, V)
        return logprob
        
    # score by log likelihood: 
    #   P score = ORF scored by markov probability
    #   Q score = reverse complement of ORF scored by markov probability
    def score(self, seq):
        """ score the sequence based on markov probability """
        P_score = self.sequence_proba(seq, self.start_counts, self.kponemer_counts, self.kmer_counts)
        Q_score = self.sequence_proba(seq, self.bg_start_counts, self.bg_kponemer_counts, self.bg_kmer_counts)
        return P_score - Q_score # since they're log probabilities, we subtract


    def results(self):
        """ format resulting probabilities into a list; start position of ORF, end position of ORF, length of ORF, markov score of ORF """
        results = []
        for seq, loc in zip(self.orfs.total_orfs, self.orfs.all_orf_locations):
            results.append({ "start" : loc[0], "end" : loc[1], "length" : len(seq), "score" : self.score(seq)})

        return results

def main():
    # golden set (our dev set)
    goldens = pd.read_csv("data/plusgenes-subset.gff", delimiter="\t", header=None)

    # input data
    data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
    seq = data[0].sequence
    
    # markov model
    mm = MarkovModel(5, 1, seq)

    print(mm)
    

if __name__ == "__main__":
    main()
