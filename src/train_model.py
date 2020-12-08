import argparse
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from MarkovModel import MarkovModel
from ORF import ORF, read_fna
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-longl", type=int, action="store")
parser.add_argument("-shortl", type=int, action="store")
parser.add_argument("-k", type=float, action="store")
parser.add_argument("-pseudo", type=int, action="store")
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

# input data
data=read_fna("data/GCF_000091665.1_ASM9166v1_genomic.fna")
seq = data[0].sequence

# golden testset
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

def roc_len_score(fig, df_results):
    fpr, tpr, thresholds = roc_curve(df_results["matches"], df_results["score"])
    fpr_len, tpr_len, thresholds_len = roc_curve(df_results["matches"], df_results["length"])

    plt.plot(fpr, tpr, "g-", label="score")
    plt.plot(fpr[64], tpr[64], "go", label="threshold at {}".format(thresholds[64]))
    plt.plot(fpr_len, tpr_len, "r-", label="length")
    plt.plot(fpr_len[243], tpr_len[243], "r*", label="threshold at {}".format(thresholds_len[243]))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC curves of Length and Markov Scores")
    plt.legend()


def scatter_len_score(fig, df_results):
    cmap = ["royalblue" if c else "firebrick" for c in df_results["matches"]]
    plt.scatter(df_results["length"],df_results["score"],c=cmap, s=5)
    plt.xlabel("Length")
    plt.ylabel("Markov Scores")
    plt.title("Length vs Markov scores, blue = gene, red = no gene")

fig = plt.figure()
roc_len_score(fig, df_results)
plt.savefig("output/roc_curve_k{}_pseudo{}_longl{}_shortl{}.png".format(k, pseudo, longl, shortl))