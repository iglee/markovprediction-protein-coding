# train the markov model
# this training script is a binding, overview script to run:
#   1. data loading
#   2. process data into ORFs using ORF class
#   3. run model training using MarkovModel class
#   4. generate relevant plots

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

parser = argparse.ArgumentParser()
parser.add_argument("-longl", type=int, action="store")
parser.add_argument("-shortl", type=int, action="store")
parser.add_argument("-k", type=int, action="store")
parser.add_argument("-pseudo", type=float, action="store")
parser.add_argument("-r", type=float, action="store")
args = parser.parse_args()

if args.longl:
    longl = args.longl
else:
    longl = 1400

if args.shortl:
    shortl = args.shortl
else:
    shortl = 50

if args.k:
    k = args.k
else:
    k = 5

if args.pseudo:
    pseudo = args.pseudo
else:
    pseudo = 1

if args.r:
    r = args.r
else:
    r = 0.2

# input data
data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
seq = data[0].sequence

# golden set (out dev set)
goldens = pd.read_csv("data/plusgenes-subset.gff", delimiter="\t", header=None)

# markov model
mm = MarkovModel(k, pseudo, seq, longl, shortl)
results = mm.results()
df_results = pd.DataFrame(results)

# organize results into df
results = mm.results()

# check matches function
def check_match(end_idx):
    stop_idx = end_idx + 3
    matches = goldens.index[goldens[4] == stop_idx]
    if len(matches) == 0:
        return False
    else:
        return True

df_results["matches"] = df_results["end"].apply(lambda x: check_match(int(x)))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def idx_nearest_match(df_results, col, thresholds):
    accuracies = []

    for t in thresholds:
        acc= accuracy_score(col > t, df_results["matches"])
        accuracies.append(acc)
    
    return find_nearest(accuracies, 0.8)

# report things!
print(mm.orfs)
print("\nreading frame 3\n", "total number: ", len(mm.orfs.idxs0), "\n", \
    " first: ", mm.orfs.idxs0[0], " (length of {})".format(len(mm.orfs.orf0[0])), "\n", \
    " last: ", mm.orfs.idxs0[-1], " (length of {})".format(len(mm.orfs.orf0[-1])), "\n")
print("reading frame 1\n", "total number: ", len(mm.orfs.idxs1), "\n", \
    " first: ", mm.orfs.idxs1[0], " (length of {})".format(len(mm.orfs.orf1[0])), "\n", \
    " last: ", mm.orfs.idxs1[-1], " (length of {})".format(len(mm.orfs.orf1[-1])), "\n")
print("reading frame 2\n", "total number: ", len(mm.orfs.idxs2), "\n", \
    " first: ", mm.orfs.idxs2[0], " (length of {})".format(len(mm.orfs.orf2[0])), "\n", \
    " last: ", mm.orfs.idxs2[-1], " (length of {})".format(len(mm.orfs.orf2[-1])), "\n")
print("total number of CDS strands: ", len(goldens[0]))

long = df_results[df_results["length"] > longl]
short = df_results[df_results["length"] < shortl]

print("shortest orfs: ")
print(short.sort_values("start")[:5])
print("longest orfs: ")
print(long.sort_values("start")[:5])


# plot things! the output of the functions below are shown in the pdf.
def roc_len_score(fig, df_results, combined_results):
    fpr, tpr, thresholds = roc_curve(df_results["matches"], df_results["score"])
    auc = roc_auc_score(df_results["matches"], df_results["score"])
    fpr_len, tpr_len, thresholds_len = roc_curve(df_results["matches"], df_results["length"])
    auc_len = roc_auc_score(df_results["matches"], df_results["length"])
    fpr_combined, tpr_combined, thresholds_combined = roc_curve(df_results["matches"], combined_results)
    auc_combined = roc_auc_score(df_results["matches"], combined_results)


    idx_score = idx_nearest_match(df_results, df_results["score"], thresholds)
    idx_len = idx_nearest_match(df_results, df_results["length"], thresholds_len)
    idx_combined = idx_nearest_match(df_results, combined_results, thresholds_combined)


    plt.plot(fpr, tpr, "g-", label="score, auc = {}".format(auc))
    plt.plot(fpr[idx_score], tpr[idx_score], "go", label="threshold at {}".format(thresholds[idx_score]))
    plt.plot(fpr_len, tpr_len, "r-", label="length, auc = {}".format(auc_len))
    plt.plot(fpr_len[idx_len], tpr_len[idx_len], "r*", label="threshold at {}".format(thresholds_len[idx_len]))
    plt.plot(fpr_combined, tpr_combined, "b-", label="combined/flashbulb, auc = {}".format(auc_combined))
    plt.plot(fpr_combined[idx_combined], tpr_combined[idx_combined], "b*", label="threshold at {}".format(thresholds_combined[idx_combined]))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC curves of Length and Markov Scores")
    plt.legend()

def zoomin(fig, xlow, xhigh, ylow, yhigh):
    plt.xlim(xlow,xhigh)
    plt.ylim(ylow,yhigh)

def scatter_len_score(fig, df_results):
    cmap = ["royalblue" if c else "firebrick" for c in df_results["matches"]]
    plt.scatter(df_results["length"],df_results["score"],c=cmap, s=5)
    plt.xlabel("Length")
    plt.ylabel("Markov Scores")
    plt.title("Length vs Markov scores, blue = gene, red = no gene")

# this "flashbulb method" was a rough approximate method to combine length and markov model score for gene prediction
# essentially identifies decision boundary plane of two features.
def flashbulb(fig, df_results,r):
    long = df_results[df_results["length"] > longl]
    short = df_results[df_results["length"] < shortl]
    (Sx,Sy) = median(short["length"]),median(short["score"])
    (Lx,Ly)=median(long["length"]), median(long["score"])

    #calculate a line through the medians:
    m=(Ly-Sy)/(Lx-Sx)
    b=Ly-m*Lx

    #calculate perpendicular lines:
    xcross = Sx + r*(Lx - Sx)
    #print(xcross)
    #print(Lx,Ly)
    #print(Sx,Sy)
    #print(longl)
    ycross = m*xcross+b
    y_intercept = ycross - (1/m)*xcross

    # plot
    scatter_len_score(fig, long)
    scatter_len_score(fig, short)
    plt.plot(median(short["length"]),median(short["score"]),"r*", 30)
    plt.plot(median(long["length"]), median(long["score"]),"b*", 30)

    ax = plt.axes()
    x = np.linspace(-100,8000,1000)
    ax.plot(x, m*x+b)
    ax.plot(x, (-1/m)*x-y_intercept)

    plt.xlim(-100, 8000)
    plt.ylim(-100, 2500)
    return m, y_intercept


fig = plt.figure()
m, y_intercept = flashbulb(fig, df_results, r)
plt.savefig("output/decision_bdy.png")
plt.close()


def combine_mm_len(df_results, m, y_intercept):
    y_pred = []

    for row in df_results.iterrows():
        x = row[1]["length"]

        temp = (-1/m)*x-y_intercept
        y_pred.append(row[1]["score"]-temp)

    return y_pred

combined_results = combine_mm_len(df_results, m, y_intercept)
fig = plt.figure()
roc_len_score(fig, df_results, combined_results)
plt.savefig("output/roc_curve_k{}_pseudo{}_longl{}_shortl{}.png".format(k, pseudo, longl, shortl))
zoomin(fig, -0.02, 0.15, 0.75, 1.03)
plt.savefig("output/roc_curve_k{}_pseudo{}_longl{}_shortl{}_zoomed.png".format(k, pseudo, longl, shortl))
plt.close()

print(mm)