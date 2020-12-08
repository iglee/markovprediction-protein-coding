import argparse
import numpy as np
from statistics import median
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from MarkovModel import MarkovModel
from ORF import ORF, read_fna
import pandas as pd

# input data
data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
seq = data[0].sequence

# golden testset
goldens = pd.read_csv("data/plusgenes-subset.gff", delimiter="\t", header=None)

def check_match(end_idx):
    stop_idx = end_idx + 3
    matches = goldens.index[goldens[4] == stop_idx]
    if len(matches) == 0:
        return False
    else:
        return True

def roc_len_score(fig, df_results, var_name, var):
    fpr, tpr, _ = roc_curve(df_results["matches"], df_results["score"])
    auc = roc_auc_score(df_results["matches"], df_results["score"])
    plt.plot(fpr, tpr, label="{} = {}, auc = {}".format(var_name, var, auc))
    
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC curves of Length and Markov Scores")
    plt.legend()

def zoomin(fig, xlow, xhigh, ylow, yhigh):
    plt.xlim(xlow,xhigh)
    plt.ylim(ylow,yhigh)

# range of long lengths
longls = [400,800,1000,1200, 1400, 2000]
fig = plt.figure()
for longl in longls:
    mm = MarkovModel(k=5, seq=seq, pseudocount=1, long_len=longl)
    results = mm.results()
    df_results = pd.DataFrame(results)
    df_results["matches"] = df_results["end"].apply(lambda x: check_match(int(x)))
    roc_len_score(fig, df_results, "long_len", longl)

plt.savefig("output/longl_roc.png")
zoomin(fig, -0.02, 0.15, 0.75, 1.03)
plt.savefig("output/longl_roc_zoomed.png")
plt.close()


# range of short lengths
shortls = [20,50,70,80,100]
fig = plt.figure()
for shortl in shortls:
    mm = MarkovModel(k=5, seq=seq, pseudocount=1, short_len=shortl)
    results = mm.results()
    df_results = pd.DataFrame(results)
    df_results["matches"] = df_results["end"].apply(lambda x: check_match(int(x)))
    roc_len_score(fig, df_results, "short_len", shortl)

plt.savefig("output/shortl_roc.png")
zoomin(fig, -0.02, 0.15, 0.75, 1.03)
plt.savefig("output/shortl_roc_zoomed.png")
plt.close()



# range of k to try
ks = [3,4,5,6,7]

fig = plt.figure()
for k in ks:
    mm = MarkovModel(k=k, seq=seq, pseudocount=1)
    results = mm.results()
    df_results = pd.DataFrame(results)
    df_results["matches"] = df_results["end"].apply(lambda x: check_match(int(x)))
    roc_len_score(fig, df_results, "k", k)

plt.savefig("output/k_roc.png")
zoomin(fig, -0.02, 0.15, 0.75, 1.03)
plt.savefig("output/k_roc_zoomed.png")
plt.close()


# range of pseudo counts to try
pseudos = [0.1,0.5,1,3,5]

fig = plt.figure()
for p in pseudos:
    mm = MarkovModel(k=5, seq=seq, pseudocount=p)
    results = mm.results()
    df_results = pd.DataFrame(results)
    df_results["matches"] = df_results["end"].apply(lambda x: check_match(int(x)))
    roc_len_score(fig, df_results, "p", p)

plt.savefig("output/p_roc.png")
zoomin(fig, -0.02, 0.15, 0.75, 1.03)
plt.savefig("output/p_roc_zoomed.png")
plt.close()


